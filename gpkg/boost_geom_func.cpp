/*
 * Copyright 2013 Luciad (http://www.luciad.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "boost_geom_io.hpp"
#include "boost_geom_wrapper.hpp"

extern "C" {
#include "error.h"
#include "geom_func.h"
#include "sql.h"
}

typedef struct {
  gpkg::GeometryPtr geometry;
  int srid;
} boostgeom_geometry_t;

/* SQLite invokes these callbacks across a C boundary, where an escaping C++ exception is undefined
   behaviour and would take the host process down rather than fail the statement. Boost.Geometry
   raises exceptions for invalid input, and any allocation can raise std::bad_alloc, so every entry
   point is registered through one of these guards and reports the failure as a SQL error instead. */
template <void (*Fn)(sqlite3_context *, int, sqlite3_value **)>
static void boostgeom_guarded(sqlite3_context *context, int nbArgs, sqlite3_value **args) {
  try {
    Fn(context, nbArgs, args);
  } catch (const std::exception &e) {
    sqlite3_result_error(context, e.what(), -1);
  } catch (...) {
    sqlite3_result_error(context, "unhandled error in geometry function", -1);
  }
}

template <void (*Fn)(sqlite3_context *)> static void boostgeom_guarded_final(sqlite3_context *context) {
  try {
    Fn(context);
  } catch (const std::exception &e) {
    sqlite3_result_error(context, e.what(), -1);
  } catch (...) {
    sqlite3_result_error(context, "unhandled error in geometry function", -1);
  }
}

static boostgeom_geometry_t *get_boost_geom(sqlite3_context *context, const spatialdb_t *spatialdb,
                                            sqlite3_value *value, errorstream_t *error) {
  geom_blob_header_t header;

  uint8_t *blob = (uint8_t *)sqlite3_value_blob(value);
  size_t blob_length = (size_t)sqlite3_value_bytes(value);

  if (blob == NULL) {
    return NULL;
  }

  binstream_t stream;
  binstream_init(&stream, blob, blob_length);

  /* A malformed blob used to be processed anyway: the failed header read left header.srid as
     uninitialized stack data and the functions computed results over whatever the geometry read
     produced, instead of reporting the error the reader had already recorded. */
  if (spatialdb->read_blob_header(&stream, &header, error) != SQLITE_OK) {
    if (error_count(error) == 0) {
      error_append(error, "Invalid geometry blob header");
    }
    return NULL;
  }

  boostgeom_writer_t writer;
  boostgeom_writer_init_srid(&writer, header.srid);

  if (spatialdb->read_geometry(&stream, boostgeom_writer_geom_consumer(&writer), error) != SQLITE_OK) {
    if (error_count(error) == 0) {
      error_append(error, "Invalid geometry blob");
    }
    boostgeom_writer_destroy(&writer, true);
    return NULL;
  }

  gpkg::GeometryPtr g = boostgeom_writer_getgeometry(&writer);

  if (g.which() == 0) {
    boostgeom_writer_destroy(&writer, true);
    return nullptr;
  } else {
    boostgeom_writer_destroy(&writer, false);
    boostgeom_geometry_t *geom = new boostgeom_geometry_t();
    geom->geometry = g;
    geom->srid = header.srid;
    return geom;
  }
}

static void free_boost_geom(void *data) {
  if (data == NULL) {
    return;
  }

  boostgeom_geometry_t *geom = (boostgeom_geometry_t *)data;
  if (geom != NULL) {
    gpkg::delete_geometry(geom->geometry);
    delete geom;
  }
}

static int set_boost_geom_result(sqlite3_context *context, const spatialdb_t *spatialdb, boostgeom_geometry_t *geom,
                                 errorstream_t *error) {
  int result = SQLITE_OK;

  if (geom == NULL) {
    sqlite3_result_null(context);
    return result;
  } else {
    geom_blob_writer_t writer;
    spatialdb->writer_init_srid(&writer, geom->srid);

    result = boostgeom_read_geometry(geom->geometry, geom_blob_writer_geom_consumer(&writer), error);

    if (result == SQLITE_OK) {
      /* SQLite takes the buffer and frees it, so the writer must not. */
      sqlite3_result_blob(context, geom_blob_writer_getdata(&writer), geom_blob_writer_length(&writer), sqlite3_free);
      spatialdb->writer_destroy(&writer, 0);
    } else {
      sqlite3_result_error(context, error_message(error), -1);
      /* Nobody took the buffer on this path; it used to be leaked. */
      spatialdb->writer_destroy(&writer, 1);
    }

    return result;
  }
}

#define BOOSTGEOM_START(context)                                                                                       \
  const spatialdb_t *spatialdb = (const spatialdb_t *)sqlite3_user_data(context);                                      \
  char error_buffer[256];                                                                                              \
  errorstream_t error;                                                                                                 \
  error_init_fixed(&error, error_buffer, 256)

#define BOOSTGEOM_GET_GEOM(name, i)                                                                                    \
  if (i >= nbArgs) {                                                                                                   \
    sqlite3_result_error(context, "wrong number of arguments", -1);                                                    \
    return;                                                                                                            \
  }                                                                                                                    \
  boostgeom_geometry_t *name = (boostgeom_geometry_t *)sqlite3_get_auxdata(context, i);                                \
  if (name == NULL) {                                                                                                  \
    name = get_boost_geom(context, spatialdb, args[i], &error);                                                        \
    if (name != NULL) {                                                                                                \
      /* Hand ownership to SQLite at once: registering only after the last argument leaked */                          \
      /* everything acquired so far when a later argument failed to parse and returned early, */                       \
      /* and set_auxdata destroys the value if it cannot store it, so re-read it afterwards. */                        \
      sqlite3_set_auxdata(context, i, (void *)name, free_boost_geom);                                                  \
      name = (boostgeom_geometry_t *)sqlite3_get_auxdata(context, i);                                                  \
    }                                                                                                                  \
  }                                                                                                                    \
  if (name == NULL) {                                                                                                  \
    if (error_count(&error) > 0) {                                                                                     \
      sqlite3_result_error(context, error_message(&error), -1);                                                        \
    } else {                                                                                                           \
      sqlite3_result_null(context);                                                                                    \
    }                                                                                                                  \
    return;                                                                                                            \
  }
/* Ownership passes to SQLite in BOOSTGEOM_GET_GEOM; nothing is left to release here. */
#define BOOSTGEOM_FREE_GEOM(name, i)                                                                                   \
  do {                                                                                                                 \
  } while (0)

#define BOOSTGEOM_GEOM__INT(sql_name, boost_name)                                                                      \
  static void ST_##sql_name(sqlite3_context *context, int nbArgs, sqlite3_value **args) {                              \
    BOOSTGEOM_START(context);                                                                                          \
    BOOSTGEOM_GET_GEOM(g1, 0);                                                                                         \
    int result = gpkg::boost_name(g1->geometry);                                                                       \
    sqlite3_result_int(context, result);                                                                               \
                                                                                                                       \
    BOOSTGEOM_FREE_GEOM(g1, 0);                                                                                        \
  }

#define BOOSTGEOM_GEOM__INT2(sql_name, boost_name)                                                                     \
  static void ST_##sql_name(sqlite3_context *context, int nbArgs, sqlite3_value **args) {                              \
    BOOSTGEOM_START(context);                                                                                          \
    BOOSTGEOM_GET_GEOM(g1, 0);                                                                                         \
    BOOSTGEOM_GET_GEOM(g2, 1);                                                                                         \
    int result = gpkg::boost_name(g1->geometry, g2->geometry);                                                         \
    sqlite3_result_int(context, result);                                                                               \
    BOOSTGEOM_FREE_GEOM(g1, 0);                                                                                        \
    BOOSTGEOM_FREE_GEOM(g2, 1);                                                                                        \
  }

#define BOOSTGEOM_GEOM__INT3(sql_name, boost_name)                                                                     \
  static void ST_##sql_name(sqlite3_context *context, int nbArgs, sqlite3_value **args) {                              \
    BOOSTGEOM_START(context);                                                                                          \
    BOOSTGEOM_GET_GEOM(g1, 0);                                                                                         \
    BOOSTGEOM_GET_GEOM(g2, 1);                                                                                         \
    gpkg::Envelope env1 = gpkg::envelope(g1->geometry);                                                                \
    gpkg::Envelope env2 = gpkg::envelope(g2->geometry);                                                                \
    gpkg::GeometryPtr genv1 = boost::variant<gpkg::Envelope *>(&env1);                                                 \
    gpkg::GeometryPtr genv2 = boost::variant<gpkg::Envelope *>(&env2);                                                 \
    int result = gpkg::boost_name(genv1, genv2);                                                                       \
    if (result == 1) {                                                                                                 \
      result = gpkg::boost_name(g1->geometry, g2->geometry);                                                           \
    }                                                                                                                  \
    sqlite3_result_int(context, result);                                                                               \
    BOOSTGEOM_FREE_GEOM(g1, 0);                                                                                        \
    BOOSTGEOM_FREE_GEOM(g2, 1);                                                                                        \
  }

#define BOOSTGEOM_GEOM__DOUBLE(sql_name, boost_name)                                                                   \
  static void ST_##sql_name(sqlite3_context *context, int nbArgs, sqlite3_value **args) {                              \
    BOOSTGEOM_START(context);                                                                                          \
    BOOSTGEOM_GET_GEOM(g1, 0);                                                                                         \
    double result = gpkg::boost_name(g1->geometry);                                                                    \
    sqlite3_result_double(context, result);                                                                            \
                                                                                                                       \
    BOOSTGEOM_FREE_GEOM(g1, 0);                                                                                        \
  }

#define BOOSTGEOM_GEOM__BLOB2(sql_name, boost_name)                                                                    \
  static void ST_##sql_name(sqlite3_context *context, int nbArgs, sqlite3_value **args) {                              \
    BOOSTGEOM_START(context);                                                                                          \
    BOOSTGEOM_GET_GEOM(g1, 0);                                                                                         \
    BOOSTGEOM_GET_GEOM(g2, 1);                                                                                         \
    boostgeom_geometry_t res;                                                                                          \
    res.srid = (g1->srid == g2->srid) ? g1->srid : 0;                                                                  \
    res.geometry = gpkg::boost_name(g1->geometry, g2->geometry);                                                       \
    set_boost_geom_result(context, spatialdb, &res, &error);                                                           \
    BOOSTGEOM_FREE_GEOM(g1, 0);                                                                                        \
    BOOSTGEOM_FREE_GEOM(g2, 1);                                                                                        \
    gpkg::delete_geometry(res.geometry);                                                                               \
  }

#define BOOSTGEOM_GEOM__BLOB1(sql_name)                                                                                \
  static void ST_##sql_name(sqlite3_context *context, int nbArgs, sqlite3_value **args) {                              \
    BOOSTGEOM_START(context);                                                                                          \
    BOOSTGEOM_GET_GEOM(g1, 0);                                                                                         \
    gpkg::Envelope env1 = gpkg::envelope(g1->geometry);                                                                \
    boostgeom_geometry_t res;                                                                                          \
    res.srid = g1->srid;                                                                                               \
    res.geometry = gpkg::clone(&env1);                                                                                 \
    set_boost_geom_result(context, spatialdb, &res, &error);                                                           \
    BOOSTGEOM_FREE_GEOM(g1, 0);                                                                                        \
    gpkg::delete_geometry(res.geometry);                                                                               \
  }

BOOSTGEOM_GEOM__INT(IsValid, is_valid)
BOOSTGEOM_GEOM__INT(IsSimple, is_simple)

BOOSTGEOM_GEOM__DOUBLE(Length, length)
BOOSTGEOM_GEOM__DOUBLE(Perimeter, perimeter)
BOOSTGEOM_GEOM__DOUBLE(Area, area)

BOOSTGEOM_GEOM__INT(NumPoints, num_points)
BOOSTGEOM_GEOM__INT(NumGeometries, num_geometries)

BOOSTGEOM_GEOM__INT2(Intersects, intersects)

//
BOOSTGEOM_GEOM__INT3(Intersects2, intersects)
BOOSTGEOM_GEOM__BLOB1(Envelop)
//

BOOSTGEOM_GEOM__BLOB2(Union, union_)

static void ST_Scale(sqlite3_context *context, int nbArgs, sqlite3_value **args) {
  BOOSTGEOM_START(context);
  BOOSTGEOM_GET_GEOM(g1, 0);
  double factor = sqlite3_value_double(args[1]);
  boostgeom_geometry_t res;
  res.srid = g1->srid;
  res.geometry = gpkg::scale(g1->geometry, factor);
  set_boost_geom_result(context, spatialdb, &res, &error);
  BOOSTGEOM_FREE_GEOM(g1, 0);
  gpkg::delete_geometry(res.geometry);
}

static void ST_Translate(sqlite3_context *context, int nbArgs, sqlite3_value **args) {
  BOOSTGEOM_START(context);
  BOOSTGEOM_GET_GEOM(g1, 0);
  double scalex = sqlite3_value_double(args[1]);
  double scaley = sqlite3_value_double(args[2]);
  boostgeom_geometry_t res;
  res.srid = g1->srid;
  res.geometry = gpkg::translate(g1->geometry, scalex, scaley);
  set_boost_geom_result(context, spatialdb, &res, &error);
  BOOSTGEOM_FREE_GEOM(g1, 0);
  gpkg::delete_geometry(res.geometry);
}

static void ST_Rotate(sqlite3_context *context, int nbArgs, sqlite3_value **args) {
  BOOSTGEOM_START(context);
  BOOSTGEOM_GET_GEOM(g1, 0);
  double angle = sqlite3_value_double(args[1]);
  boostgeom_geometry_t res;
  res.srid = g1->srid;

  if (nbArgs == 2) {
    res.geometry = gpkg::rotate(g1->geometry, angle);
  } else if (nbArgs == 4) {
    double centerx = sqlite3_value_double(args[2]);
    double centery = sqlite3_value_double(args[3]);
    res.geometry = gpkg::rotate(g1->geometry, angle, centerx, centery);
  }

  set_boost_geom_result(context, spatialdb, &res, &error);
  BOOSTGEOM_FREE_GEOM(g1, 0);
  gpkg::delete_geometry(res.geometry);
}

static void ST_Affine(sqlite3_context *context, int nbArgs, sqlite3_value **args) {
  BOOSTGEOM_START(context);
  BOOSTGEOM_GET_GEOM(g1, 0);
  double a = sqlite3_value_double(args[1]);
  double b = sqlite3_value_double(args[2]);
  double d = sqlite3_value_double(args[3]);
  double e = sqlite3_value_double(args[4]);
  double xoff = sqlite3_value_double(args[5]);
  double yoff = sqlite3_value_double(args[6]);
  boostgeom_geometry_t res;
  res.srid = g1->srid;
  res.geometry = gpkg::affine(g1->geometry, a, b, d, e, xoff, yoff);
  set_boost_geom_result(context, spatialdb, &res, &error);
  BOOSTGEOM_FREE_GEOM(g1, 0);
  gpkg::delete_geometry(res.geometry);
}

static void GPKG_BoostGeometryVersion(sqlite3_context *context, int nbArgs, sqlite3_value **args) {
  int boost_major = BOOST_VERSION / 100000;
  int boost_minor = (BOOST_VERSION / 100) % 1000;
  int boost_subminor = BOOST_VERSION % 100;
  char *version = sqlite3_mprintf("%d.%d.%d", boost_major, boost_minor, boost_subminor);
  if (version) {
    sqlite3_result_text(context, version, -1, SQLITE_TRANSIENT);
    sqlite3_free(version);
  } else {
    sqlite3_result_error(context, "Could not obtain Boost.Geometry version number", -1);
  }
}

struct chain_item {
  boostgeom_geometry_t *geom;
  struct chain_item *next;
}; // struct used to store a chain item

struct geom_chain {
  struct chain_item *first;
  struct chain_item *last;
}; // struct used to store a dynamic chain of boost geometries

static void aggregate_geom_step(sqlite3_context *context, int argc, sqlite3_value **argv) {
  if (argc != 1) {
    sqlite3_result_error(context, "wrong number of arguments", -1);
    return;
  }

  if (sqlite3_value_type(argv[0]) != SQLITE_BLOB) {
    sqlite3_result_null(context);
    return;
  }

  BOOSTGEOM_START(context);
  boostgeom_geometry_t *g1 = (boostgeom_geometry_t *)sqlite3_get_auxdata(context, 0);

  if (g1 == NULL) {
    g1 = get_boost_geom(context, spatialdb, argv[0], &error);
  }

  if (g1 == NULL) {
    if (error_count(&error) > 0) {
      sqlite3_result_error(context, error_message(&error), -1);
    } else {
      sqlite3_result_null(context);
    }

    return;
  }

  struct geom_chain *chain;

  struct geom_chain **p = (geom_chain **)sqlite3_aggregate_context(context, sizeof(struct geom_chain **));

  bool firstRow = !(*p);

  if (firstRow) {
    chain = (geom_chain *)malloc(sizeof(struct geom_chain));
    *p = chain;
  } else {
    chain = *p;
  }

  struct chain_item *item = (chain_item *)malloc(sizeof(struct chain_item));

  item->geom = g1;

  item->next = NULL;

  if (firstRow) {
    chain->first = item;
  } else {
    chain->last->next = item;
  }

  chain->last = item;
}

static void free_geom_chain(struct geom_chain *chain) {
  struct chain_item *item = chain->first;
  while (item) {
    struct chain_item *next = item->next;
    free_boost_geom(item->geom);
    free(item);
    item = next;
  }
  free(chain);
}

static void aggregate_geom_Union_final(sqlite3_context *context) {
  boostgeom_geometry_t res;
  res.geometry = ::boost::blank();
  res.srid = 0;
  struct geom_chain **p = (geom_chain **)sqlite3_aggregate_context(context, 0);

  if (!p || !*p) {
    sqlite3_result_null(context);
    return;
  }

  struct geom_chain *chain = *p;

  /* SQLite calls this finalizer exactly once, so everything the walk has accumulated must be
     released even when Boost.Geometry throws; the guard then reports the exception as a SQL
     error. The old code freed items as it went and leaked the whole tail on a throw. */
  try {
    bool first = true;
    for (struct chain_item *item = chain->first; item != NULL; item = item->next) {
      if (first) {
        // copy the geometry value so free_geom_chain can free item->geom
        res.geometry = gpkg::clone(item->geom->geometry);
        res.srid = item->geom->srid;
        first = false;
      } else {
        gpkg::GeometryPtr tmp = gpkg::union_(res.geometry, item->geom->geometry);
        gpkg::delete_geometry(res.geometry);
        res.geometry = tmp;
        res.srid = (item->geom->srid == res.srid) ? res.srid : 0;
      }
    }
  } catch (...) {
    free_geom_chain(chain);
    gpkg::delete_geometry(res.geometry);
    throw;
  }

  free_geom_chain(chain);

  try {
    if (gpkg::is_empty(res.geometry)) {
      sqlite3_result_null(context);
    } else {
      BOOSTGEOM_START(context);
      set_boost_geom_result(context, spatialdb, &res, &error);
    }
  } catch (...) {
    gpkg::delete_geometry(res.geometry);
    throw;
  }

  gpkg::delete_geometry(res.geometry);
}

//
static void extent_step(sqlite3_context *context, int argc, sqlite3_value **argv) {
  // TODO: Take into account 2 hemisphere

  if (argc != 1) {
    sqlite3_result_error(context, "wrong number of arguments", -1);
    return;
  }

  if (sqlite3_value_type(argv[0]) != SQLITE_BLOB) {
    sqlite3_result_null(context);
    return;
  }

  BOOSTGEOM_START(context);
  boostgeom_geometry_t *g1 = (boostgeom_geometry_t *)sqlite3_get_auxdata(context, 0);

  /* Auxdata is not a useful cache across the steps of an aggregate, so a geometry parsed here is
     owned by this step and must be released before it returns. */
  int g1_owned = 0;
  if (g1 == NULL) {
    g1 = get_boost_geom(context, spatialdb, argv[0], &error);
    g1_owned = 1;
  }

  if (g1 == NULL) {
    if (error_count(&error) > 0) {
      sqlite3_result_error(context, error_message(&error), -1);
    } else {
      sqlite3_result_null(context);
    }

    return;
  }

  // Envelop
  gpkg::Envelope env1 = gpkg::envelope(g1->geometry);
  double minx = env1.min_corner().get<0>();
  double miny = env1.min_corner().get<1>();
  double maxx = env1.max_corner().get<0>();
  double maxy = env1.max_corner().get<1>();

  // Extent
  double *max_min;
  int *srid_check;

  double **p = (double **)sqlite3_aggregate_context(context, sizeof(double **));

  if (!(*p)) {
    // this is the first row
    max_min = (double *)malloc((sizeof(double) * 5));
    *(max_min + 0) = minx;
    *(max_min + 1) = miny;
    *(max_min + 2) = maxx;
    *(max_min + 3) = maxy;
    srid_check = (int *)(max_min + 4);
    *(srid_check + 0) = g1->srid;
    *(srid_check + 1) = g1->srid;
    *p = max_min;
  } else {
    // subsequent rows
    max_min = *p;

    if (minx < *(max_min + 0)) {
      *(max_min + 0) = minx;
    }

    if (miny < *(max_min + 1)) {
      *(max_min + 1) = miny;
    }

    if (maxx > *(max_min + 2)) {
      *(max_min + 2) = maxx;
    }

    if (maxy > *(max_min + 3)) {
      *(max_min + 3) = maxy;
    }

    srid_check = (int *)(max_min + 4);

    /* Latch a mismatch instead of recording the current row's SRID: overwriting the slot every
       row only caught a difference when the very last row disagreed with the first, so a run of
       4326, 3857, 4326 was accepted and mixed coordinate systems into one box. */
    if (*(srid_check + 0) != g1->srid) {
      *(srid_check + 1) = g1->srid;
    }
  }

  if (g1_owned) {
    free_boost_geom(g1);
  }
}

static void extent_final(sqlite3_context *context) {
  double **p = (double **)sqlite3_aggregate_context(context, 0);

  if (!p) {
    sqlite3_result_null(context);
    return;
  }

  double *max_min = *p;

  if (!max_min) {
    sqlite3_result_null(context);
    return;
  }

  int *srid_check = (int *)(max_min + 4);

  if (*(srid_check + 0) != *(srid_check + 1)) {
    /* The accumulator belongs to this call whichever way it ends; returning here used to leak it. */
    free(max_min);
    sqlite3_result_null(context);
    return;
  }

  /* Copy the accumulator and release it before anything that can throw; allocating the result
     polygon first used to leak it when the allocation failed. */
  double minx = *(max_min + 0);
  double miny = *(max_min + 1);
  double maxx = *(max_min + 2);
  double maxy = *(max_min + 3);
  int srid = *(srid_check + 0);

  free(max_min);

  boostgeom_geometry_t res;
  res.geometry = new gpkg::Polygon({{gpkg::Point(minx, miny), gpkg::Point(maxx, miny), gpkg::Point(maxx, maxy),
                                     gpkg::Point(minx, maxy), gpkg::Point(minx, miny)}});
  res.srid = srid;

  try {
    if (gpkg::is_empty(res.geometry)) {
      sqlite3_result_null(context);
    } else {
      BOOSTGEOM_START(context);
      set_boost_geom_result(context, spatialdb, &res, &error);
    }
  } catch (...) {
    gpkg::delete_geometry(res.geometry);
    throw;
  }

  gpkg::delete_geometry(res.geometry);
}

static void extent_coord_step(sqlite3_context *context, int argc, sqlite3_value **argv) {
  // TODO: Take into account 2 hemisphere

  if (argc != 4) {
    sqlite3_result_error(context, "Wrong number of arguments", -1);
    return;
  }

  /* Read the coordinates through the API. The previous code cast sqlite3_value* straight to
     double* and dereferenced it, which yielded the right number only because the double happens
     to sit first in SQLite's internal value struct, and it consulted auxdata, which carries no
     meaning in an aggregate step. Integer arguments were rejected outright. */
  double xy[4];

  for (int i = 0; i < 4; ++i) {
    int type = sqlite3_value_type(argv[i]);
    if (type != SQLITE_FLOAT && type != SQLITE_INTEGER) {
      sqlite3_result_null(context);
      return;
    }

    xy[i] = sqlite3_value_double(argv[i]);
  }

  // Extent
  double *max_min;

  double **p = (double **)sqlite3_aggregate_context(context, sizeof(double **));

  double minx = xy[0];
  double miny = xy[1];
  double maxx = xy[2];
  double maxy = xy[3];

  if (!(*p)) {
    // this is the first row
    max_min = (double *)malloc((sizeof(double) * 4));
    *(max_min + 0) = minx;
    *(max_min + 1) = miny;
    *(max_min + 2) = maxx;
    *(max_min + 3) = maxy;
    *p = max_min;
  } else {
    // subsequent rows
    max_min = *p;

    if (minx < *(max_min + 0)) {
      *(max_min + 0) = minx;
    }

    if (miny < *(max_min + 1)) {
      *(max_min + 1) = miny;
    }

    if (maxx > *(max_min + 2)) {
      *(max_min + 2) = maxx;
    }

    if (maxy > *(max_min + 3)) {
      *(max_min + 3) = maxy;
    }
  }
}

static void extent_coord_final(sqlite3_context *context) {
  double **p = (double **)sqlite3_aggregate_context(context, 0);

  if (!p) {
    sqlite3_result_null(context);
    return;
  }

  double *max_min = *p;

  if (!max_min) {
    sqlite3_result_null(context);
    return;
  }

  /* Copy the accumulator and release it before anything that can throw; allocating the result
     polygon first used to leak it when the allocation failed. */
  double minx = *(max_min + 0);
  double miny = *(max_min + 1);
  double maxx = *(max_min + 2);
  double maxy = *(max_min + 3);

  free(max_min);

  boostgeom_geometry_t res;
  res.geometry = new gpkg::Polygon({{gpkg::Point(minx, miny), gpkg::Point(maxx, miny), gpkg::Point(maxx, maxy),
                                     gpkg::Point(minx, maxy), gpkg::Point(minx, miny)}});
  res.srid = 0;

  try {
    if (gpkg::is_empty(res.geometry)) {
      sqlite3_result_null(context);
    } else {
      BOOSTGEOM_START(context);
      set_boost_geom_result(context, spatialdb, &res, &error);
    }
  } catch (...) {
    gpkg::delete_geometry(res.geometry);
    throw;
  }

  gpkg::delete_geometry(res.geometry);
}

#define STR(x) #x

#define BOOSTGEOM_FUNCTION(db, prefix, name, nbArgs, ctx, error)                                                       \
  do {                                                                                                                 \
    sql_create_function(db, STR(prefix##_##name), boostgeom_guarded<prefix##_##name>, nbArgs, SQL_DETERMINISTIC,       \
                        (void *)ctx, NULL, error);                                                                     \
  } while (0)

extern "C" {
void geom_func_init(sqlite3 *db, const spatialdb_t *spatialdb, errorstream_t *error) {
  BOOSTGEOM_FUNCTION(db, GPKG, BoostGeometryVersion, 0, spatialdb, error);
  BOOSTGEOM_FUNCTION(db, ST, IsValid, 1, spatialdb, error);
  BOOSTGEOM_FUNCTION(db, ST, IsSimple, 1, spatialdb, error);

  BOOSTGEOM_FUNCTION(db, ST, Length, 1, spatialdb, error);
  BOOSTGEOM_FUNCTION(db, ST, Perimeter, 1, spatialdb, error);
  BOOSTGEOM_FUNCTION(db, ST, Area, 1, spatialdb, error);

  BOOSTGEOM_FUNCTION(db, ST, NumPoints, 1, spatialdb, error);
  BOOSTGEOM_FUNCTION(db, ST, NumGeometries, 1, spatialdb, error);

  BOOSTGEOM_FUNCTION(db, ST, Intersects, 2, spatialdb, error);

  BOOSTGEOM_FUNCTION(db, ST, Union, 2, spatialdb, error);
  sqlite3_create_function_v2(db, "ST_Union", 1, SQLITE_UTF8 | SQLITE_DETERMINISTIC, (void *)spatialdb, 0,
                             boostgeom_guarded<aggregate_geom_step>,
                             boostgeom_guarded_final<aggregate_geom_Union_final>, 0);

  BOOSTGEOM_FUNCTION(db, ST, Scale, 2, spatialdb, error);
  BOOSTGEOM_FUNCTION(db, ST, Translate, 3, spatialdb, error);
  BOOSTGEOM_FUNCTION(db, ST, Rotate, 2, spatialdb, error);
  BOOSTGEOM_FUNCTION(db, ST, Rotate, 4, spatialdb, error);
  BOOSTGEOM_FUNCTION(db, ST, Affine, 7, spatialdb, error);

  //
  BOOSTGEOM_FUNCTION(db, ST, Intersects2, 2, spatialdb, error);
  BOOSTGEOM_FUNCTION(db, ST, Envelop, 1, spatialdb, error);
  sqlite3_create_function_v2(db, "ST_Extent", 1, SQLITE_UTF8 | SQLITE_DETERMINISTIC, (void *)spatialdb, 0,
                             boostgeom_guarded<extent_step>, boostgeom_guarded_final<extent_final>, 0);
  sqlite3_create_function_v2(db, "ST_Extent", 4, SQLITE_UTF8 | SQLITE_DETERMINISTIC, (void *)spatialdb, 0,
                             boostgeom_guarded<extent_coord_step>, boostgeom_guarded_final<extent_coord_final>, 0);
  //
}
}
