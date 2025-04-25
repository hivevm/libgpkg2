/* Copyright 2025 HiveVM (http://www.hivevm.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include <boost/algorithm/string.hpp>
#include <boost/geometry/geometry.hpp>
#include <boost/variant/variant.hpp>

#include "boost_geometries.hpp"
#include "boost_geom_wrapper.hpp"

namespace gpkg {

class boost_projection {
  const Transform *_trans;

public:
  inline boost_projection(const Transform *trans)
      : _trans(trans) {}

  inline bool apply(gpkg::Point const &p1, gpkg::Point &p2) const {
    namespace coord = boost::geometry;

    double x = coord::get<0>(p1);
    double y = coord::get<1>(p1);

    _trans->reproject(&x, &y);

    coord::set<0>(p2, x);
    coord::set<1>(p2, y);

    return true;
  }
};

namespace detail {
template <typename Visitor> struct geometry_ptr_visitor : boost::static_visitor<typename Visitor::result_type> {
  typename Visitor::result_type null_result;

  geometry_ptr_visitor(typename Visitor::result_type null_result)
      : null_result(null_result) {}

  // arguments: 1 Geometry

  typename Visitor::result_type operator()(const boost::blank &) const { return null_result; }

  template <typename Geometry> typename Visitor::result_type operator()(const Geometry *value) const {
    return value != nullptr ? Visitor()(*value) : null_result;
  }

  // arguments: 2 Geometries

  typename Visitor::result_type operator()(const boost::blank &, const boost::blank &) const { return null_result; }

  template <typename Geom> typename Visitor::result_type operator()(const boost::blank &, const Geom *) const {
    return null_result;
  }

  template <typename Geom> typename Visitor::result_type operator()(const Geom *, const boost::blank &) const {
    return null_result;
  }

  template <typename Geom1, typename Geom2>
  typename Visitor::result_type operator()(const Geom1 *value1, const Geom2 *value2) const {
    return ((value1 != nullptr) && (value2 != nullptr)) ? Visitor()(*value1, *value2) : null_result;
  }

  // arguments: 1 Geometry, 1 int

  typename Visitor::result_type operator()(const boost::blank &, const int &) const { return null_result; }

  template <typename Geom> typename Visitor::result_type operator()(const Geom *value1, const int &value2) const {
    return (value1 != nullptr) ? Visitor()(*value1, value2) : null_result;
  }

  // arguments: 1 Geometry, 1 double

  typename Visitor::result_type operator()(const boost::blank &, const double &) const { return null_result; }

  template <typename Geom> typename Visitor::result_type operator()(const Geom *value1, const double &value2) const {
    return (value1 != nullptr) ? Visitor()(*value1, value2) : null_result;
  }

  // arguments: 1 Geometry, 2 doubles

  typename Visitor::result_type operator()(const boost::blank &, const double &, const double &) const {
    return null_result;
  }

  template <typename Geom>
  typename Visitor::result_type operator()(const Geom *value1, const double &value2, const double &value3) const {
    return (value1 != nullptr) ? Visitor()(*value1, value2, value3) : null_result;
  }

  // arguments: 1 Geometry, 3 doubles

  typename Visitor::result_type operator()(const boost::blank &, const double &, const double &, const double &) const {
    return null_result;
  }

  template <typename Geom>
  typename Visitor::result_type operator()(const Geom *value1, const double &value2, const double &value3,
                                           const double &value4) const {
    return (value1 != nullptr) ? Visitor()(*value1, value2, value3, value4) : null_result;
  }

  // arguments: 1 Geometry, 6 doubles

  typename Visitor::result_type operator()(const boost::blank &, const double &, const double &, const double &,
                                           const double &, const double &, const double &) const {
    return null_result;
  }

  template <typename Geom>
  typename Visitor::result_type operator()(const Geom *value1, const double &value2, const double &value3,
                                           const double &value4, const double &value5, const double &value6,
                                           const double &value7) const {
    return (value1 != nullptr) ? Visitor()(*value1, value2, value3, value4, value5, value6, value7) : null_result;
  }
};

struct clone : boost::static_visitor<GeometryPtr> {
  template <typename Geometry> GeometryPtr operator()(const Geometry &geom) const {
    Geometry *gc = new Geometry();
    *gc = geom;
    return gc;
  }
};

struct equal : boost::static_visitor<int> {

  int operator()(const gpkg::GeometryCollection &gc1, const gpkg::GeometryCollection &gc2) const {
    if (gc1.size() != gc2.size()) {
      return 0;
    }

    auto it1 = gc1.cbegin();
    auto it2 = gc2.cbegin();

    for (; ((it1 < gc1.cend()) && (it2 < gc2.cend())); ++it1, ++it2) {
      if (!boost::apply_visitor(*this, *it1, *it2)) {
        return 0;
      }
    }

    return 1;
  }

  template <typename Geometry> int operator()(const Geometry &geom1, const Geometry &geom2) const {
    return boost::geometry::equals(geom1, geom2);
  }

  template <typename Geometry1, typename Geometry2> int operator()(const Geometry1 &, const Geometry2 &) const {
    return 0;
  }
};

struct envelope : boost::static_visitor<GeometryPtr> {

  GeometryPtr operator()(const gpkg::GeometryCollection &geom) const {
    GeometryPtr envelope;

    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      GeometryPtr part_env = boost::apply_visitor(*this, *it);

      if (it == geom.cbegin()) {
        envelope = part_env;
        continue;
      }

      GeometryPtr env_union = union_(envelope, part_env);

      delete_geometry(envelope);
      delete_geometry(part_env);

      envelope = env_union;
    }

    return envelope;
  }

  template <typename Geometry> GeometryPtr operator()(const Geometry &geom) const {
    Envelope *envelope = new Envelope();
    boost::geometry::envelope(geom, *envelope);
    return envelope;
  }
};

struct wkt : boost::static_visitor<std::string> {

  std::string operator()(const gpkg::GeometryCollection &geom, const int &precision) const {
    std::string collection = "GEOMETRYCOLLECTION (";
    auto bound_visitor = std::bind(*this, std::placeholders::_1, precision);

    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      if (it != geom.cbegin()) {
        collection.append(",");
      }

      collection.append(boost::apply_visitor(bound_visitor, *it));
    }

    collection.append(")");
    return collection;
  }

  std::string operator()(const gpkg::Envelope &geom, const int &precision) const {
    std::stringstream ss;

    if (precision >= 0) {
      ss << std::fixed << std::setprecision(precision);
    }

    ss << "ENVELOPE " << boost::geometry::dsv(geom);
    std::string boostWkt = ss.str();
    boost::replace_all(boostWkt, ", ", ",");

    return boostWkt;
  }

  template <typename Geometry> std::string operator()(const Geometry &geom, const int &precision) const {
    std::stringstream ss;

    if (precision >= 0) {
      ss << std::fixed << std::setprecision(precision);
    }

    ss << boost::geometry::wkt(geom);
    std::string boostWkt = ss.str();
    boost::replace_first(boostWkt, "(", " (");

    return boostWkt;
  }
};

struct is_valid : boost::static_visitor<int> {
  template <typename Geometry> int operator()(const Geometry &geom) const { return boost::geometry::is_valid(geom); }

  int operator()(const gpkg::GeometryCollection &geom) const {
    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      if (!boost::apply_visitor(*this, *it)) {
        return false;
      }
    }
    return true;
  }
};

struct is_simple : boost::static_visitor<int> {
  template <typename Geometry> int operator()(const Geometry &geom) const { return boost::geometry::is_simple(geom); }

  int operator()(const gpkg::GeometryCollection &geom) const {
    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      if (!boost::apply_visitor(*this, *it)) {
        return false;
      }
    }
    return true;
  }
};

struct is_empty : boost::static_visitor<int> {
  template <typename Geometry> int operator()(const Geometry &geom) const { return boost::geometry::is_empty(geom); }

  int operator()(const gpkg::GeometryCollection &geom) const {
    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      if (!boost::apply_visitor(*this, *it)) {
        return false;
      }
    }
    return true;
  }
};

struct area : boost::static_visitor<double> {
  template <typename Geometry> double operator()(const Geometry &geom) const { return boost::geometry::area(geom); }

  double operator()(const gpkg::GeometryCollection &geom) const {
    double total = 0.0;

    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      total += boost::apply_visitor(*this, *it);
    }

    return total;
  }
};

struct distance : boost::static_visitor<double> {

  double operator()(const gpkg::GeometryCollection &geomColl1, const gpkg::GeometryCollection &geomColl2) const {

    double minDist = DBL_MAX;

    for (auto it1 = geomColl1.cbegin(); it1 < geomColl1.cend(); ++it1) {
      for (auto it2 = geomColl2.cbegin(); it2 < geomColl2.cend(); ++it2) {
        double dist = boost::apply_visitor(*this, *it1, *it2);

        if (dist < minDist) {
          minDist = dist;
        }
      }
    }

    return minDist;
  }

  template <typename Geometry> double operator()(const gpkg::GeometryCollection &geomColl, const Geometry &geom) const {
    return this->operator()(geom, geomColl);
  }

  template <typename Geometry> double operator()(const Geometry &geom, const gpkg::GeometryCollection &geomColl) const {
    double minDist = DBL_MAX;
    const gpkg::Geometry &g = static_cast<gpkg::Geometry>(geom);

    for (auto it = geomColl.cbegin(); it < geomColl.cend(); ++it) {
      double dist = boost::apply_visitor(*this, *it, g);

      if (dist < minDist) {
        minDist = dist;
      }
    }

    return minDist;
  }

  template <typename Geometry1, typename Geometry2>
  double operator()(const Geometry1 &geom1, const Geometry2 &geom2) const {
    return boost::geometry::distance(geom1, geom2);
  }
};

struct length : boost::static_visitor<double> {
  template <typename Geometry> double operator()(const Geometry &geom) const { return boost::geometry::length(geom); }

  double operator()(const gpkg::GeometryCollection &geom) const {
    double total = 0.0;

    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      total += boost::apply_visitor(*this, *it);
    }

    return total;
  }
};

struct perimeter : boost::static_visitor<double> {
  template <typename Geometry> double operator()(const Geometry &geom) const {
    return boost::geometry::perimeter(geom);
  }

  double operator()(const gpkg::GeometryCollection &geom) const {
    double total = 0.0;

    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      total += boost::apply_visitor(*this, *it);
    }

    return total;
  }
};

struct num_points : boost::static_visitor<int> {
  template <typename Geometry> int operator()(const Geometry &geom) const { return boost::geometry::num_points(geom); }

  int operator()(const gpkg::GeometryCollection &geom) const {
    int total = 0;

    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      total += boost::apply_visitor(*this, *it);
    }

    return total;
  }
};

struct num_rings : boost::static_visitor<int> {
  template <typename Geometry> int operator()(const Geometry &) const {
    // point or linestring or not yet implemented
    return 0;
  }

  int operator()(const gpkg::Polygon &geom) const {
    return static_cast<int>(boost::geometry::num_interior_rings(geom));
  }
};

struct num_geometries : boost::static_visitor<int> {
  template <typename Geometry> int operator()(const Geometry &geom) const {
    return boost::geometry::num_geometries(geom);
  }

  int operator()(const gpkg::GeometryCollection &geom) const {
    int total = 0;

    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      total++;
    }

    return total;
  }
};

struct contains : boost::static_visitor<int> {

  int operator()(const gpkg::Polygon &geom1, const gpkg::Polygon &geom2) const {
    return boost::geometry::within(geom2, geom1);
  }

  int operator()(const gpkg::MultiPolygon &geom1, const gpkg::Polygon &geom2) const {
    return boost::geometry::within(geom2, geom1);
  }

  int operator()(const gpkg::Polygon &geom1, const gpkg::MultiPolygon &geom2) const {
    return boost::geometry::within(geom2, geom1);
  }

  int operator()(const gpkg::MultiPolygon &geom1, const gpkg::MultiPolygon &geom2) const {
    return boost::geometry::within(geom2, geom1);
  }

  template <typename Geometry1, typename Geometry2> int operator()(const Geometry1 &, const Geometry2 &) const {
    return 0; // not yet implemented (many combinations are unavailable)
  }
};

struct intersects : boost::static_visitor<int> {

  int operator()(const gpkg::GeometryCollection &geomColl1, const gpkg::GeometryCollection &geomColl2) const {
    for (auto it1 = geomColl1.cbegin(); it1 < geomColl1.cend(); ++it1) {
      for (auto it2 = geomColl2.cbegin(); it2 < geomColl2.cend(); ++it2) {
        if (boost::apply_visitor(*this, *it1, *it2)) {
          return 1;
        }
      }
    }

    return 0;
  }

  template <typename Geometry> int operator()(const gpkg::GeometryCollection &geomColl, const Geometry &geom) const {
    return this->operator()(geom, geomColl);
  }

  template <typename Geometry> int operator()(const Geometry &geom, const gpkg::GeometryCollection &geomColl) const {
    const gpkg::Geometry &g = static_cast<gpkg::Geometry>(geom);

    for (auto it = geomColl.cbegin(); it < geomColl.cend(); ++it) {
      if (boost::apply_visitor(*this, *it, g)) {
        return 1;
      }
    }

    return 0;
  }

  template <typename Geometry1, typename Geometry2>
  int operator()(const Geometry1 &geom1, const Geometry2 &geom2) const {
    return boost::geometry::intersects(geom1, geom2);
  }
};

struct union_ : boost::static_visitor<GeometryPtr> {

  GeometryPtr operator()(const gpkg::GeometryCollection &geomColl1, const gpkg::GeometryCollection &geomColl2) const {
    GeometryCollection *gc = new GeometryCollection();
    *gc = geomColl1;
    gc->insert(gc->end(), geomColl2.begin(), geomColl2.end());
    return gc;
  }

  template <typename Geometry>
  GeometryPtr operator()(const Geometry &geom, const gpkg::GeometryCollection &geomColl) const {
    GeometryCollection *gc = new GeometryCollection();
    gc->push_back(geom);
    gc->insert(gc->end(), geomColl.begin(), geomColl.end());
    return gc;
  }

  template <typename Geometry>
  GeometryPtr operator()(const gpkg::GeometryCollection &geomColl, const Geometry &geom) const {
    GeometryCollection *gc = new GeometryCollection();
    *gc = geomColl;
    gc->push_back(geom);
    return gc;
  }

  GeometryPtr operator()(const gpkg::Envelope &, const gpkg::GeometryCollection &) const {
    return boost::blank(); // this is an error! cannot make union of envelope
                           // and geometry
  }

  GeometryPtr operator()(const gpkg::GeometryCollection &geomColl, const gpkg::Envelope &e) const {
    return this->operator()(e, geomColl);
  }

  template <typename Geometry> GeometryPtr operator()(const Geometry &, const gpkg::Envelope &) const {
    return boost::blank(); // this is an error! cannot make union of envelope
                           // and geometry
  }

  template <typename Geometry> GeometryPtr operator()(const gpkg::Envelope &e, const Geometry &geom) const {
    return this->operator()(geom, e);
  }

  GeometryPtr operator()(const gpkg::Envelope &e1, const gpkg::Envelope &e2) const {
    gpkg::Polygon poly1, poly2;
    boost::geometry::convert(e1, poly1);
    boost::geometry::convert(e2, poly2);

    gpkg::MultiPolygon envUnion;
    boost::geometry::union_(poly1, poly2, envUnion);

    gpkg::Envelope *box = new gpkg::Envelope();
    boost::geometry::envelope(envUnion, *box);

    return box;
  }

#define BOOST_UNION(Type)                                                                                              \
  Multi##Type *geomUnion = new Multi##Type();                                                                          \
  boost::geometry::union_(geom1, geom2, *geomUnion);                                                                   \
  return geomUnion;

#define UNION_OPERATORS(Type)                                                                                          \
  GeometryPtr operator()(const Type &geom1, const Type &geom2)                                                         \
      const {BOOST_UNION(Type)} GeometryPtr operator()(const Type &geom1, const Multi##Type &geom2)                    \
          const {BOOST_UNION(Type)} GeometryPtr operator()(const Multi##Type &geom1, const Type &geom2)                \
              const {BOOST_UNION(Type)} GeometryPtr operator()(const Multi##Type &geom1, const Multi##Type &geom2)     \
                  const {BOOST_UNION(Type)}

  UNION_OPERATORS(Point)
  UNION_OPERATORS(LineString)
  UNION_OPERATORS(Polygon)

  template <typename Geometry1, typename Geometry2>
  GeometryPtr operator()(const Geometry1 &geom1, const Geometry2 &geom2) const {
    GeometryCollection *gc = new GeometryCollection();
    gc->push_back(geom1);
    gc->push_back(geom2);
    return gc;
  }
};

struct coordinates : boost::static_visitor<std::vector<gpkg::Point>> {

  std::vector<gpkg::Point> operator()(const gpkg::Point &geom, int) const {
    std::vector<gpkg::Point> vec;
    vec.push_back(geom);
    return vec;
  }

  std::vector<gpkg::Point> operator()(const gpkg::LineString &geom, int) const { return geom; }

  std::vector<gpkg::Point> operator()(const gpkg::Polygon &geom, int ring) const {
    return (ring == 0) ? geom.outer() : geom.inners()[static_cast<std::size_t>(ring) - 1];
  }

  template <typename Geometry> std::vector<gpkg::Point> operator()(const Geometry &, int) const {
    // not yet implemented
    return std::vector<gpkg::Point>();
  }
};

struct subgeometry : boost::static_visitor<GeometryPtr> {

  GeometryPtr operator()(const gpkg::MultiPoint &geom, int idx) const {
    gpkg::Point *p = new gpkg::Point();
    *p = geom[static_cast<std::size_t>(idx)];
    return p;
  }

  GeometryPtr operator()(const gpkg::MultiLineString &geom, int idx) const {
    gpkg::LineString *l = new gpkg::LineString();
    *l = geom[static_cast<std::size_t>(idx)];
    return l;
  }

  GeometryPtr operator()(const gpkg::MultiPolygon &geom, int idx) const {
    gpkg::Polygon *p = new gpkg::Polygon();
    *p = geom[static_cast<std::size_t>(idx)];
    return p;
  }

  GeometryPtr operator()(const gpkg::GeometryCollection &geom, int idx) const {
    return boost::apply_visitor(detail::clone(), geom[static_cast<std::size_t>(idx)]);
  }

  template <typename Geometry> GeometryPtr operator()(const Geometry &, int) const { return boost::blank(); }
};

struct scale : boost::static_visitor<GeometryPtr> {

  GeometryPtr operator()(const gpkg::GeometryCollection &geom, double factor) const {
    GeometryPtr gc = new GeometryCollection();
    auto bound_visitor = std::bind(*this, std::placeholders::_1, factor);

    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      GeometryPtr scaled = boost::apply_visitor(bound_visitor, *it);
      GeometryPtr _union = gpkg::union_(gc, scaled);
      delete_geometry(gc);
      delete_geometry(scaled);
      gc = _union;
    }

    return gc;
  }

  template <typename Geometry> GeometryPtr operator()(const Geometry &geom, double factor) const {
    Geometry *res = new Geometry();
    boost::geometry::strategy::transform::scale_transformer<double, 2, 2> scale(factor);
    boost::geometry::transform(geom, *res, scale);
    return res;
  }
};

struct translate : boost::static_visitor<GeometryPtr> {

  GeometryPtr operator()(const gpkg::GeometryCollection &geom, double tx, double ty) const {
    GeometryPtr gc = new GeometryCollection();
    auto bound_visitor = std::bind(*this, std::placeholders::_1, tx, ty);

    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      GeometryPtr translated = boost::apply_visitor(bound_visitor, *it);
      GeometryPtr _union = gpkg::union_(gc, translated);
      delete_geometry(gc);
      delete_geometry(translated);
      gc = _union;
    }

    return gc;
  }

  template <typename Geometry> GeometryPtr operator()(const Geometry &geom, double tx, double ty) const {
    Geometry *res = new Geometry();
    boost::geometry::strategy::transform::translate_transformer<double, 2, 2> translate(tx, ty);
    boost::geometry::transform(geom, *res, translate);
    return res;
  }
};

struct rotate : boost::static_visitor<GeometryPtr> {

  GeometryPtr operator()(const gpkg::GeometryCollection &geom, double angle, double cx, double cy) const {
    GeometryPtr gc = new GeometryCollection();
    auto bound_visitor = std::bind(*this, std::placeholders::_1, angle, cx, cy);

    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      GeometryPtr rotated = boost::apply_visitor(bound_visitor, *it);
      GeometryPtr _union = gpkg::union_(gc, rotated);
      delete_geometry(gc);
      delete_geometry(rotated);
      gc = _union;
    }

    return gc;
  }

  template <typename Geometry> GeometryPtr operator()(const Geometry &geom, double angle, double cx, double cy) const {
    Geometry *res_t1 = new Geometry();
    Geometry *res_r = new Geometry();
    Geometry *res_t2 = new Geometry();
    // Rotate is calculated with respect to (0,0) coordinates. To rotate a
    // geometry respect to its center, cx and cy should be the coordinates if
    // the center. If this is the case, we translate the geometry in order to
    // have its center in (0,0) position, rotate and then retranslate it in the
    // original position.
    boost::geometry::strategy::transform::translate_transformer<double, 2, 2> translate1(-cx, -cy);
    boost::geometry::strategy::transform::translate_transformer<double, 2, 2> translate2(+cx, +cy);
    boost::geometry::strategy::transform::rotate_transformer<boost::geometry::degree, double, 2, 2> rotate(angle);
    boost::geometry::transform(geom, *res_t1, translate1);
    boost::geometry::transform(*res_t1, *res_r, rotate);
    boost::geometry::transform(*res_r, *res_t2, translate2);
    delete res_t1;
    delete res_r;
    return res_t2;
  }
};

struct reverse : boost::static_visitor<GeometryPtr> {

  GeometryPtr operator()(const gpkg::GeometryCollection &geom) const {
    GeometryPtr gc = new GeometryCollection();
    auto bound_visitor = std::bind(*this, std::placeholders::_1);

    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      GeometryPtr reversed = boost::apply_visitor(bound_visitor, *it);
      GeometryPtr _union = gpkg::union_(gc, reversed);
      delete_geometry(gc);
      delete_geometry(reversed);
      gc = _union;
    }

    return gc;
  }

  template <typename Geometry> GeometryPtr operator()(const Geometry &geom) const {
    Geometry *res = new Geometry(geom);
    boost::geometry::reverse(*res);
    return res;
  }
};

struct affine : boost::static_visitor<GeometryPtr> {

  GeometryPtr operator()(const gpkg::GeometryCollection &geom, double a, double b, double d, double e, double xoff,
                         double yoff) const {
    GeometryPtr gc = new GeometryCollection();
    auto bound_visitor = std::bind(*this, std::placeholders::_1, a, b, d, e, xoff, yoff);

    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      GeometryPtr transformed = boost::apply_visitor(bound_visitor, *it);
      GeometryPtr _union = gpkg::union_(gc, transformed);
      delete_geometry(gc);
      delete_geometry(transformed);
      gc = _union;
    }

    return gc;
  }

  template <typename Geometry>
  GeometryPtr operator()(const Geometry &geom, double a, double b, double d, double e, double xoff, double yoff) const {
    Geometry *res = new Geometry();

    /*  Affine transformation matrix:

        /  a  b  xoff  \
        |  d  e  yoff  |
        \  0  0     1  /
     */

    boost::geometry::strategy::transform::matrix_transformer<double, 2, 2> affine(a, b, xoff, d, e, yoff, 0, 0, 1);
    boost::geometry::transform(geom, *res, affine);
    return res;
  }
};

template <typename Strategy> struct project : boost::static_visitor<GeometryPtr> {

  GeometryPtr operator()(const gpkg::GeometryCollection &geom, const Strategy &strategy) const {
    GeometryPtr gc = new GeometryCollection();
    auto bound_visitor = std::bind(*this, std::placeholders::_1, strategy);

    for (auto it = geom.cbegin(); it < geom.cend(); ++it) {
      GeometryPtr transformed = boost::apply_visitor(bound_visitor, *it);
      GeometryPtr _union = gpkg::union_(gc, transformed);
      delete_geometry(gc);
      delete_geometry(transformed);
      gc = _union;
    }

    return gc;
  }

  template <typename Geometry> GeometryPtr operator()(const Geometry &geom, const Strategy &strategy) const {
    Geometry *res = new Geometry();
    boost::geometry::transform(geom, *res, strategy);
    return res;
  }
};

struct delete_geometry_visitor : boost::static_visitor<> {

  void operator()(boost::blank &) const {}

  template <typename Geometry> void operator()(Geometry *geom) const { delete geom; }
};

} // namespace detail

void delete_geometry(GeometryPtr &ptr) { boost::apply_visitor(detail::delete_geometry_visitor(), ptr); }

double length(const GeometryPtr &geometry) {
  return boost::apply_visitor(detail::geometry_ptr_visitor<detail::length>(0.0), geometry);
}

double perimeter(const GeometryPtr &geometry) {
  return boost::apply_visitor(detail::geometry_ptr_visitor<detail::perimeter>(0.0), geometry);
}

double area(const GeometryPtr &geometry) {
  return boost::apply_visitor(detail::geometry_ptr_visitor<detail::area>(0.0), geometry);
}

double distance(const GeometryPtr &geometry1, const GeometryPtr &geometry2) {
  return boost::apply_visitor(detail::geometry_ptr_visitor<detail::distance>(-1.0), geometry1, geometry2);
}

int is_valid(const GeometryPtr &geometry) {
  return boost::apply_visitor(detail::geometry_ptr_visitor<detail::is_valid>(true), geometry);
}

int is_simple(const GeometryPtr &geometry) {
  return boost::apply_visitor(detail::geometry_ptr_visitor<detail::is_simple>(true), geometry);
}

int is_empty(const GeometryPtr &geometry) {
  return boost::apply_visitor(detail::geometry_ptr_visitor<detail::is_empty>(true), geometry);
}

int num_points(const GeometryPtr &geometry) {
  return boost::apply_visitor(detail::geometry_ptr_visitor<detail::num_points>(0), geometry);
}

int num_int_rings(const GeometryPtr &geometry) {
  return boost::apply_visitor(detail::geometry_ptr_visitor<detail::num_rings>(0), geometry);
}

int num_geometries(const GeometryPtr &geometry) {
  return boost::apply_visitor(detail::geometry_ptr_visitor<detail::num_geometries>(0), geometry);
}

std::vector<gpkg::Point> getCoordinates(const GeometryPtr &geometry, int ring) {
  auto bound_visitor = std::bind(detail::geometry_ptr_visitor<detail::coordinates>(std::vector<gpkg::Point>()),
                                 std::placeholders::_1, ring);
  return boost::apply_visitor(bound_visitor, geometry);
}

GeometryPtr getPart(const GeometryPtr &geometry, int index) {
  auto bound_visitor =
      std::bind(detail::geometry_ptr_visitor<detail::subgeometry>(boost::blank()), std::placeholders::_1, index);
  return (index >= num_geometries(geometry)) ? boost::blank() : boost::apply_visitor(bound_visitor, geometry);
}

GeometryPtr clone(const GeometryPtr &geometry) {
  return boost::apply_visitor(detail::geometry_ptr_visitor<detail::clone>(boost::blank()), geometry);
}

int equal(const GeometryPtr &geometry1, const GeometryPtr &geometry2) {
  if (is_empty(geometry1) && is_empty(geometry2)) {
    return 1;
  }
  return boost::apply_visitor(detail::geometry_ptr_visitor<detail::equal>(0), geometry1, geometry2);
}

GeometryPtr reverse(const GeometryPtr &geometry) {
  auto bound_visitor = std::bind(detail::geometry_ptr_visitor<detail::reverse>(boost::blank()), std::placeholders::_1);
  return boost::apply_visitor(bound_visitor, geometry);
}

Envelope envelope(const GeometryPtr &geometry) {
  GeometryPtr res = boost::apply_visitor(detail::geometry_ptr_visitor<detail::envelope>(boost::blank()), geometry);
  Envelope e;

  try {
    e = *(boost::get<gpkg::Envelope *>(res));
  } catch (std::exception &) { // blank geometry
    boost::geometry::assign_values(e, DBL_MAX, DBL_MAX, -DBL_MAX, -DBL_MAX);
  }

  delete_geometry(res);
  return e;
}

std::string wkt(const GeometryPtr &geometry, const int &precision) {
  auto bound_visitor = std::bind(detail::geometry_ptr_visitor<detail::wkt>(""), std::placeholders::_1, precision);
  std::string geomWkt = boost::apply_visitor(bound_visitor, geometry);

  boost::replace_all(geomWkt, ",", ", ");
  return geomWkt;
}

int intersects(const GeometryPtr &geometry1, const GeometryPtr &geometry2) {
  return boost::apply_visitor(detail::geometry_ptr_visitor<detail::intersects>(0), geometry1, geometry2);
}

int contains(const GeometryPtr &geometry1, const GeometryPtr &geometry2) {
  return boost::apply_visitor(detail::geometry_ptr_visitor<detail::contains>(0), geometry1, geometry2);
}

GeometryPtr union_(const GeometryPtr &geometry1, const GeometryPtr &geometry2) {
  return boost::apply_visitor(detail::geometry_ptr_visitor<detail::union_>(boost::blank()), geometry1, geometry2);
}

GeometryPtr scale(const GeometryPtr &geometry, double factor) {
  auto bound_visitor =
      std::bind(detail::geometry_ptr_visitor<detail::scale>(boost::blank()), std::placeholders::_1, factor);
  return boost::apply_visitor(bound_visitor, geometry);
}

GeometryPtr translate(const GeometryPtr &geometry, double tx, double ty) {
  auto bound_visitor =
      std::bind(detail::geometry_ptr_visitor<detail::translate>(boost::blank()), std::placeholders::_1, tx, ty);
  return boost::apply_visitor(bound_visitor, geometry);
}

GeometryPtr rotate(const GeometryPtr &geometry, double degrees, double cx, double cy) {
  auto bound_visitor =
      std::bind(detail::geometry_ptr_visitor<detail::rotate>(boost::blank()), std::placeholders::_1, degrees, cx, cy);
  return boost::apply_visitor(bound_visitor, geometry);
}

GeometryPtr affine(const GeometryPtr &geometry, double a, double b, double d, double e, double xoff, double yoff) {
  auto bound_visitor = std::bind(detail::geometry_ptr_visitor<detail::affine>(boost::blank()), std::placeholders::_1, a,
                                 b, d, e, xoff, yoff);
  return boost::apply_visitor(bound_visitor, geometry);
}

GeometryPtr project(const GeometryPtr &geometry, const Transform *transform) {
  gpkg::boost_projection strategy(transform);
  auto bound_visitor = std::bind(detail::geometry_ptr_visitor<detail::project<gpkg::boost_projection>>(boost::blank()),
                                 std::placeholders::_1, &strategy);
  return boost::apply_visitor(bound_visitor, geometry);
}
} // namespace gpkg
