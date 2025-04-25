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
#ifndef BOOST_GEOM_WRAPPER_HPP
#define BOOST_GEOM_WRAPPER_HPP

#include "boost_geometries.hpp"
#include "boost_geom_proj.hpp"

namespace gpkg {
double length(const GeometryPtr &geometry);
double perimeter(const GeometryPtr &geometry);
double area(const GeometryPtr &geometry);
double distance(const GeometryPtr &geometry1, const GeometryPtr &geometry2);

int is_valid(const GeometryPtr &geometry);
int is_simple(const GeometryPtr &geometry);
int is_empty(const GeometryPtr &geometry);

int num_points(const GeometryPtr &geometry);
int num_int_rings(const GeometryPtr &geometry);
int num_geometries(const GeometryPtr &geometry);

std::vector<gpkg::Point> getCoordinates(const GeometryPtr &geometry, int ring);
GeometryPtr getPart(const GeometryPtr &geometry, int index);

GeometryPtr clone(const GeometryPtr &geometry);
int equal(const GeometryPtr &geometry1, const GeometryPtr &geometry2);
GeometryPtr reverse(const GeometryPtr &geometry);

Envelope envelope(const GeometryPtr &geometry);
std::string wkt(const GeometryPtr &geometry, const int &precision);

int intersects(const GeometryPtr &geometry1, const GeometryPtr &geometry2);
int contains(const GeometryPtr &geometry1, const GeometryPtr &geometry2);
GeometryPtr union_(const GeometryPtr &geometry1, const GeometryPtr &geometry2);

GeometryPtr scale(const GeometryPtr &geometry, double factor);
GeometryPtr translate(const GeometryPtr &geometry, double tx, double ty);
GeometryPtr rotate(const GeometryPtr &geometry, double degrees, double cx = 0.0, double cy = 0.0);
GeometryPtr affine(const GeometryPtr &geometry, double a, double b, double d, double e, double xoff, double yoff);

GeometryPtr project(const GeometryPtr &geometry, const Transform *transform);
} // namespace gpkg

#endif
