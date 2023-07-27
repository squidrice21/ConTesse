// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include <Eigen/Core>

struct PatchRotationAngle {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Vector2f w1, v, w2;
  int r;
  int ori;

  PatchRotationAngle()
      : w1(Vector2f(0, 0)), v(Vector2f(0, 0)), w2(Vector2f(0, 0)), r(0),
        ori(1) {}
  PatchRotationAngle(const Vector2f &x, const Vector2f &y, const Vector2f &z,
                     bool to_init = true);
  void init_orientation();
};

PatchRotationAngle add_triangle(PatchRotationAngle const &acc,
                                PatchRotationAngle const &tri);

bool operator==(PatchRotationAngle const &alpha,
                PatchRotationAngle const &beta);
bool operator!=(PatchRotationAngle const &alpha,
                PatchRotationAngle const &beta);
