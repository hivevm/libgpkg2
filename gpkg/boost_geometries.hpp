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
#ifndef GPKG_BOOST_GEOMETRIES_HPP
#define GPKG_BOOST_GEOMETRIES_HPP

#include <boost/geometry/geometry.hpp>
#include <vector>

namespace gpkg {
typedef ::boost::geometry::model::point<double, 2, ::boost::geometry::cs::cartesian> Point;
typedef ::boost::geometry::model::ring<Point, false, true> LinearRing;
typedef ::boost::geometry::model::linestring<Point> LineString;
typedef ::boost::geometry::model::polygon<Point, false, true> Polygon;
typedef ::boost::geometry::model::multi_point<Point> MultiPoint;
typedef ::boost::geometry::model::multi_linestring<LineString> MultiLineString;
typedef ::boost::geometry::model::multi_polygon<Polygon> MultiPolygon;
typedef ::boost::geometry::model::box<Point> Envelope;
typedef ::boost::variant<Point, LineString, Polygon, LinearRing, MultiPoint, MultiLineString, MultiPolygon, Envelope>
    Geometry;
typedef std::vector<Geometry> GeometryCollection;

typedef ::boost::variant<::boost::blank, Point *, LineString *, Polygon *, LinearRing *, MultiPoint *,
                         MultiLineString *, MultiPolygon *, GeometryCollection *, Envelope *>
    GeometryPtr;

void delete_geometry(GeometryPtr &ptr);
} // namespace gpkg

#endif
