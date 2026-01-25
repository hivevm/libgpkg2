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
#include "boost_geom_proj.hpp"
#include <math.h>

#define GLOBAL_MERCATOR 40075016.68
#define WGS84_MAJOR_AXIS 6378137      // WGS84 semi-major axis
#define WGS84_MINOR_AXIS 6356752.3142 // WGS84 semi-minor axis

/******************************************************
 * Spherical radius.
 * http://en.wikipedia.org/wiki/Earth_radius
 * http://en.wikipedia.org/wiki/World_Geodetic_System
 ******************************************************/
// #define WGS84_MAJOR_AXIS            6378137.0
// #define WGS84_INVERSE_FLATTENING    298.257223563
// #define WGS84_MINOR_AXIS            (WGS84_MAJOR_AXIS - WGS84_MAJOR_AXIS /
// WGS84_INVERSE_FLATTENING) #define WGS84_RADIUS                ((2.0 *
// WGS84_MAJOR_AXIS + WGS84_MINOR_AXIS ) / 3.0)

static double realmod(const double a, const double b) {
  unsigned long long div = static_cast<unsigned long long>(a / b);
  return a - static_cast<double>(div) * b;
}

template <typename T> auto qLn(T v) {
  using std::log;
  return log(v);
}

template <typename T> auto qExp(T v) {
  using std::exp;
  return exp(v);
}

template <typename T> auto qTan(T v) {
  using std::tan;
  return tan(v);
}

template <typename T> auto qAtan(T v) {
  using std::atan;
  return atan(v);
}
template <typename T> constexpr inline const T &qMin(const T &a, const T &b) { return (a < b) ? a : b; }
template <typename T> constexpr inline const T &qMax(const T &a, const T &b) { return (a < b) ? b : a; }
template <typename T> constexpr inline const T &qBound(const T &min, const T &val, const T &max) {
  return qMax(min, qMin(max, val));
}

static void wgs84ToMeractor(double *x, double *y) {
  double lon = *y / 360.0 + 0.5;
  double lat = *x;
  lat = 0.5 - qLn(qTan((M_PI / 4.0) + (M_PI / 2.0) * lat / 180.0)) / M_PI / 2.0;
  lat = qBound(0.0, lat, 1.0); // latitudes > 85° are clamped to 85° (and
                               // mirrored on other emisphere)
  *y = (lon - 0.5) * GLOBAL_MERCATOR;
  *x = (0.5 - lat) * GLOBAL_MERCATOR;
}

static void meractorToWGS84(double *x, double *y) {
  double fx = 0.5 + *x / GLOBAL_MERCATOR;
  double fy = 0.5 - *y / GLOBAL_MERCATOR;

  if (fy < 0.0) {
    fy = 0.0;
  } else if (fy > 1.0) {
    fy = 1.0;
  }

  double lon;

  if (fy == 0.0) {
    lon = 90.0;
  } else if (fy == 1.0) {
    lon = -90.0;
  } else {
    lon = (180.0 / M_PI) * (2.0 * qAtan(qExp(M_PI * (1.0 - 2.0 * fy))) - (M_PI / 2.0));
  }

  double lat;

  if (fx >= 0) {
    lat = realmod(fx, 1.0);
  } else {
    lat = realmod(1.0 - realmod(-1.0 * fx, 1.0), 1.0);
  }

  lat = lat * 360.0 - 180.0;
  *x = lat;
  *y = lon;
}

class Wgs84Projection : public Transform {

public:
  explicit Wgs84Projection() = default;

  void reproject(double *x, double *y, double *) const override { meractorToWGS84(x, y); }
};

class GlobalProjection : public Transform {

public:
  explicit GlobalProjection() = default;

  void reproject(double *x, double *y, double *) const override {
    // Implementation of this projection needs to invert x with y!!!
    wgs84ToMeractor(y, x);
  }
};

static Transform *_WGS84ToMercator = new GlobalProjection;
static Transform *_MercatorToWGS84 = new Wgs84Projection;

void Transform::reproject(double *, double *, double *) const {}

Transform *Transform::wgs84ToMercator() {
  _WGS84ToMercator->setSridTarget(srid::SRID_GLOBAL_MERCATOR);
  _WGS84ToMercator->setSridSource(srid::SRID_WGS84);
  return _WGS84ToMercator;
}

Transform *Transform::mercatorToWGS84() {
  _MercatorToWGS84->setSridSource(srid::SRID_GLOBAL_MERCATOR);
  _MercatorToWGS84->setSridTarget(srid::SRID_WGS84);
  return _MercatorToWGS84;
}

int Transform::sridSource() const { return _sridSource; }

void Transform::setSridSource(int sridSource) { _sridSource = sridSource; }

int Transform::sridTarget() const { return _sridTarget; }

void Transform::setSridTarget(int sridTarget) { _sridTarget = sridTarget; }
