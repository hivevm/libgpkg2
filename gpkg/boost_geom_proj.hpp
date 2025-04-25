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
#ifndef BOOST_GEOM_PROJ_HPP
#define BOOST_GEOM_PROJ_HPP

namespace srid {
const int SRID_WGS84 = 4326;
const int SRID_GLOBAL_MERCATOR = 3857;
} // namespace srid

class Transform {
public:
  explicit Transform() = default;

  virtual void reproject(double *x, double *y, double *z = nullptr) const;

  static Transform *wgs84ToMercator();
  static Transform *mercatorToWGS84();

  void setSridSource(int sridSource);
  int sridSource() const;
  void setSridTarget(int sridTarget);
  int sridTarget() const;

private:
  int _sridSource;
  int _sridTarget;
};

#endif
