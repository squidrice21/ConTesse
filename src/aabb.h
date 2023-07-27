// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

/*
    aabb.h -- basic axis-aligned bounding box & ray intersection code

    This file is part of the implementation of

        Instant Field-Aligned Meshes
        Wenzel Jakob, Daniele Panozzo, Marco Tarini, and Olga Sorkine-Hornung
        In ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2015)

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE.txt file.
*/

#ifndef _AABB_H
#define _AABB_H

#include "common.h"
#include <igl/predicates/predicates.h>
#include <vector>

struct Ray {
  Vector3f o, d;
  real_t mint, maxt;

  Ray(const Vector3f &o, const Vector3f &d)
      : o(o), d(d), mint(0), maxt(std::numeric_limits<real_t>::infinity()) {}

  Ray(const Vector3f &o, const Vector3f &d, real_t mint, real_t maxt)
      : o(o), d(d), mint(mint), maxt(maxt) {}

  Vector3f operator()(real_t t) const { return o + t * d; }
};

struct AABB {
  Vector3f min, max;

  AABB() { clear(); }

  AABB(const Vector3f &min, const Vector3f &max) : min(min), max(max) {}

  void clear() {
    const real_t inf = std::numeric_limits<real_t>::infinity();
    min.setConstant(inf);
    max.setConstant(-inf);
  }

  void expandBy(const Vector3f &p) {
    min = min.cwiseMin(p);
    max = max.cwiseMax(p);
  }

  void expandBy(const AABB &aabb) {
    min = min.cwiseMin(aabb.min);
    max = max.cwiseMax(aabb.max);
  }

  bool contains(const Vector3f &p) {
    return (p.array() >= min.array()).all() && (p.array() <= max.array()).all();
  }

  bool rayIntersect(const Ray &ray) const {
    real_t nearT = -std::numeric_limits<real_t>::infinity();
    real_t farT = std::numeric_limits<real_t>::infinity();

    for (int i = 0; i < 3; i++) {
      real_t origin = ray.o[i];
      real_t minVal = min[i], maxVal = max[i];

      if (ray.d[i] == 0) {
        if (origin < minVal || origin > maxVal)
          return false;
      } else {
        real_t t1 = (minVal - origin) / ray.d[i];
        real_t t2 = (maxVal - origin) / ray.d[i];

        if (t1 > t2)
          std::swap(t1, t2);

        nearT = std::max(t1, nearT);
        farT = std::min(t2, farT);

        if (!(nearT <= farT))
          return false;
      }
    }

    return ray.mint <= farT && nearT <= ray.maxt;
  }

  bool planeIntersect(const Vector3f &p, const Vector3f &norm) const {
    // Check if corners are on different sides / on the plane.
    Vector3f p2 = norm.cross(Vector3f(1, 0, 0));
    Vector3f p3 = norm.cross(p2);
    std::vector<Vector3f> corners({this->min, this->max});
    bool has_seen_pos = false, has_seen_neg = false;
    for (uint8_t i = 0; i < 8; i++) {
      int x_idx = (i >> 2) & 1;
      int y_idx = (i >> 1) & 1;
      int z_idx = i & 1;
      Vector3f corner(corners[x_idx].x(), corners[y_idx].y(),
                      corners[z_idx].z());
      auto side = igl::predicates::orient3d(p, p2, p3, corner);
      switch (side) {
      case igl::predicates::Orientation::NEGATIVE:
        has_seen_neg = true;
        break;
      case igl::predicates::Orientation::POSITIVE:
        has_seen_pos = true;
        break;
      default: // Orientation::COPLANAR
        has_seen_neg = true;
        has_seen_pos = true;
        break;
      }

      if (has_seen_pos && has_seen_neg) {
        break;
      }
    }

    return has_seen_pos && has_seen_neg;
  }

  bool aabbIntersect(const AABB &aabb) const {
    return (abs(this->center().x() - aabb.center().x()) * 2 <=
            ((this->max.x() - this->min.x()) +
             (aabb.max.x() - aabb.min.x()))) &&
           (abs(this->center().y() - aabb.center().y()) * 2 <=
            ((this->max.y() - this->min.y()) +
             (aabb.max.y() - aabb.min.y()))) &&
           (abs(this->center().z() - aabb.center().z()) * 2 <=
            ((this->max.z() - this->min.z()) + (aabb.max.z() - aabb.min.z())));
  }

  real_t squaredDistanceTo(const Vector3f &p) const {
    real_t result = 0;
    for (int i = 0; i < 3; ++i) {
      real_t value = 0;
      if (p[i] < min[i])
        value = min[i] - p[i];
      else if (p[i] > max[i])
        value = p[i] - max[i];
      result += value * value;
    }
    return result;
  }

  Vector3f diagonal() const { return max - min; }

  int largestAxis() const {
    Vector3f extents = diagonal();

    if (extents[0] >= extents[1] && extents[0] >= extents[2])
      return 0;
    else if (extents[1] >= extents[0] && extents[1] >= extents[2])
      return 1;
    else
      return 2;
  }

  real_t surfaceArea() const {
    Vector3f d = diagonal();
    return 2.0f * (d[0] * d[1] + d[0] * d[2] + d[1] * d[2]);
  }

  Vector3f center() const { return 0.5f * (min + max); }

  static AABB merge(const AABB &aabb1, const AABB &aabb2) {
    return AABB(aabb1.min.cwiseMin(aabb2.min), aabb1.max.cwiseMax(aabb2.max));
  }
};

#endif
