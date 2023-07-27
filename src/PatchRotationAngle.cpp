// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "PatchRotationAngle.h"

#include "shor.h"

double angle_orientation(PatchRotationAngle const &K,
                         Eigen::Matrix<double, 3, 2> const &tri) {
  double ori = orientation(tri);
  // The Weber-Zorin has CW angle as the reference
  // If we have CCW angle, we need to flip the ori (used for the convex/reflex
  // test)
  if (K.ori < 0)
    ori *= -1;
  return ori;
}

void convex_reflex(PatchRotationAngle const &angle, bool &reflex_or_180,
                   bool &is_zero) {
  Eigen::Matrix<double, 3, 2> tri;
  tri << angle.w1.transpose(), angle.v.transpose(),
      angle.w2.transpose(); // w1,v,w2
  double ori = angle_orientation(angle, tri);
  reflex_or_180 = (ori < 0);
  is_zero = false;

  if (ori == 0) // w1,v,w2 are collinear
  {
    Eigen::Matrix<double, 3, 2> T1;
    Eigen::Matrix<double, 3, 2> T2;
    Eigen::RowVector2d ref;
    ref << 0, 0;
    T1 << ref, angle.v, angle.w1;
    if (orientation(T1) == 0) {
      T1.row(0) << 1, 0;
      ref << 1, 0;
      if (orientation(T1) == 0) {
        T1.row(0) << 0, 1;
        ref << 0, 1;
      }
    }
    T2 << ref, angle.v, angle.w2;
    if ((orientation(T1) < 0) !=
        (orientation(T2) < 0)) // w1 and w2 on different sides of v
      reflex_or_180 = true;
    else
      is_zero = true;
  }
}
/* =============================================== */

PatchRotationAngle::PatchRotationAngle(const Vector2f &x, const Vector2f &y,
                                       const Vector2f &z, bool to_init)
    : w1(x), v(y), w2(z), r(0) {
  if (to_init)
    init_orientation();
}

void PatchRotationAngle::init_orientation() {
  // Determine ori (assume the angle is initialized from a triangle)
  Eigen::Matrix<double, 3, 2> tri;
  tri << w1.transpose(), v.transpose(), w2.transpose();
  double tri_ori = orientation(tri);
  ori = (tri_ori == 0) ? 0 : ((tri_ori > 0) ? 1 : -1); // Only keep the sign

  bool reflex_or_180, is_zero;
  convex_reflex(*this, reflex_or_180, is_zero);

  // Set an arbitrary ori when the triangle angle is 0
  if (is_zero)
    ori = 1;
}

PatchRotationAngle weber_zorin_add(PatchRotationAngle const &acc,
                                   PatchRotationAngle const &K) {
  assert(K.w1 == acc.w2);
  assert((acc.ori <= 0) == (K.ori <= 0));

  auto w3 = K.w2;
  PatchRotationAngle S(acc.w1, acc.v, w3, false);
  S.ori = acc.ori;

  Eigen::Matrix<double, 3, 2> tri1, tri2, tri3;
  tri1 << acc.w1.transpose(), acc.v.transpose(), acc.w2.transpose(); // w1,v,w2
  tri2 << acc.w2.transpose(), acc.v.transpose(), w3.transpose();     // w2,v,w3
  tri3 << acc.w1.transpose(), acc.v.transpose(), w3.transpose();     // w1,v,w3

  double ori1 = angle_orientation(acc, tri1);
  double ori2 = angle_orientation(K, tri2);

  // Reset the ori based on the angle CW/CCW as well
  double ori3 = angle_orientation(acc, tri3);

  bool reflex_or_180_1 = (ori1 < 0);
  bool reflex_or_180_2 = (ori2 < 0);
  bool is_zero_1 = false;
  bool is_zero_2 = false;
  if (ori1 == 0) // w1,v,w2 are collinear
  {
    Eigen::Matrix<double, 3, 2> T1;
    Eigen::Matrix<double, 3, 2> T2;
    Eigen::RowVector2d ref;
    ref << 0, 0;
    T1 << ref, acc.v, acc.w1;
    if (orientation(T1) == 0) {
      T1.row(0) << 1, 0;
      ref << 1, 0;
      if (orientation(T1) == 0) {
        T1.row(0) << 0, 1;
        ref << 0, 1;
      }
    }
    T2 << ref, acc.v, acc.w2;
    if ((orientation(T1) < 0) !=
        (orientation(T2) < 0)) // w1 and w2 on different sides of v
      reflex_or_180_1 = true;
    else
      is_zero_1 = true;
  }
  if (ori2 == 0) // w3,v,w2 are collinear
  {
    Eigen::Matrix<double, 3, 2> T1;
    Eigen::Matrix<double, 3, 2> T2;
    Eigen::RowVector2d ref;
    ref << 0, 0;
    T1 << ref, acc.v, w3;
    if (orientation(T1) == 0) {
      T1.row(0) << 1, 0;
      ref << 1, 0;
      if (orientation(T1) == 0) {
        T1.row(0) << 0, 1;
        ref << 0, 1;
      }
    }
    T2 << ref, acc.v, acc.w2;
    if ((orientation(T1) < 0) !=
        (orientation(T2) < 0)) // w1 and w2 on different sides of v
      reflex_or_180_2 = true;
    else
      is_zero_2 = true;
  }
  // COUNTER-EXAMPLE: PI + 0
  bool add_one = (is_zero_1 || is_zero_2)
                     ? false
                     : (((reflex_or_180_1 != reflex_or_180_2) && ori3 >= 0) ||
                        (reflex_or_180_1 && reflex_or_180_2));

  if (add_one)
    S.r = acc.r + K.r + 1;
  else
    S.r = acc.r + K.r;

  if (is_zero_1)
    S.ori = K.ori;

  return S;
}

PatchRotationAngle add_flipped_orientation(PatchRotationAngle const &acc,
                                           PatchRotationAngle const &K) {
  assert(K.w1 == acc.w2);
  assert((acc.ori <= 0) != (K.ori <= 0));

  assert(K.ori != 0 &&
         "Triangle angle is 180 deg. Impossible to determine the ori.");

  auto w3 = K.w2;
  PatchRotationAngle S(acc.w1, acc.v, w3, false);

  bool reflex_or_180_1, reflex_or_180_2;
  bool is_zero_1, is_zero_2;
  convex_reflex(acc, reflex_or_180_1, is_zero_1);
  convex_reflex(K, reflex_or_180_2, is_zero_2);

  bool subtract_one = false;

  // When the accumulated angle is convex, it is possible to have the resulting
  // rotation index subtracted by 1.
  if (!is_zero_2 && !reflex_or_180_1) {
    Eigen::Matrix<double, 3, 2> tri1, tri2, tri3;
    tri3 << acc.w1.transpose(), acc.v.transpose(), w3.transpose(); // w1,v,w3

    double ori3 = orientation(tri3);
    // If the new w3 is degenerated (ori3 == 0), the new angle is 0 deg or 180
    // deg, the rotation index unchanged; Otherwise, if w3 and acc.w2 is not on
    // the same side (including when acc is 0 deg), the rotation index is
    // subtracted by 1.
    if ((ori3 != 0 && is_zero_1) ||
        (ori3 != 0 && (ori3 > 0) != (acc.ori > 0))) {
      subtract_one = true;
    }
  }

  S.r = acc.r - K.r;
  S.ori = acc.ori;
  if (subtract_one) {
    if (S.r > 0)
      S.r -= 1;
    else // S.r <= 0
    {
      S.r = -S.r;
      S.ori = K.ori;
    }
  }

  return S;
}

PatchRotationAngle add_triangle(PatchRotationAngle const &acc,
                                PatchRotationAngle const &tri) {
  PatchRotationAngle result;

  // Uninit accumulater
  if (acc.w1 == acc.v && acc.w1 == acc.w2)
    return tri;

  // If the ori is the same, then use the same logics as Weber-Zorin
  if ((acc.ori <= 0) == (tri.ori <= 0)) {
    result = weber_zorin_add(acc, tri);
  }
  // Difference logics if we need to 'subtract' the triangle angle (assumed to
  // be convex)
  else {
    result = add_flipped_orientation(acc, tri);
  }

  return result;
}

bool operator==(PatchRotationAngle const &alpha,
                PatchRotationAngle const &beta) {
  return alpha.w1 == beta.w1 && alpha.v == beta.v && alpha.w2 == beta.w2 &&
         alpha.r == beta.r && alpha.ori == beta.ori;
}

bool operator!=(PatchRotationAngle const &alpha,
                PatchRotationAngle const &beta) {
  return !(alpha == beta);
}
