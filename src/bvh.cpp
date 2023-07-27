// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

/*
    bvh.cpp -- bounding volume hierarchy for fast ray-intersection queries

    This file is part of the implementation of

        Instant Field-Aligned Meshes
        Wenzel Jakob, Daniele Panozzo, Marco Tarini, and Olga Sorkine-Hornung
        In ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2015)

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE.txt file.
*/

#include "bvh.h"

#include "logger.h"
#include <spdlog/fmt/ostr.h>

#include <unordered_set>

using namespace std;

struct Bins {
  static const int BIN_COUNT = 8;
  Bins() { memset(counts, 0, sizeof(uintp_t) * BIN_COUNT); }
  uintp_t counts[BIN_COUNT];
  AABB bounds[BIN_COUNT];
};

struct BVHBuildTask : public tbb::task {
  enum { SERIAL_THRESHOLD = 32 };
  BVH &bvh;
  uintp_t node_idx;
  uintp_t *start, *end;
  uintp_t *temp;

  /**
   * Create a new build task
   *
   * \param bvh
   *    Reference to the underlying BVH
   *
   * \param node_idx
   *    Index of the BVH node that should be built
   *
   * \param start
   *    Start pointer into a list of triangle indices to be processed
   *
   * \param end
   *    End pointer into a list of triangle indices to be processed
   *
   *  \param temp
   *    Pointer into a temporary memory region that can be used for
   *    construction purposes. The usable length is <tt>end-start</tt>
   *    unsigned integers.
   */
  BVHBuildTask(BVH &bvh, uintp_t node_idx, uintp_t *start, uintp_t *end,
               uintp_t *temp)
      : bvh(bvh), node_idx(node_idx), start(start), end(end), temp(temp) {}

  task *execute() {
    const MapXu F = MapXu(bvh.mF->data(), 3, bvh.mF->size());
    const std::vector<Vector3f> &V = *bvh.mV;
    uintp_t size = end - start, total_size = F.cols();
    BVHNode &node = bvh.mNodes[node_idx];

    if (size < SERIAL_THRESHOLD) {
      tbb::blocked_range<uintp_t> range(start - bvh.mIndices.data(),
                                        end - bvh.mIndices.data());
      const ProgressCallback &progress = bvh.mProgress;
      SHOW_PROGRESS_RANGE(range, total_size,
                          "Constructing Bounding Volume Hierarchy");
#if defined(SINGLE_PRECISION)
      execute_serially(bvh, node_idx, start, end, (real_t *)temp);
#else
      real_t *temp2 = new real_t[size];
      execute_serially(bvh, node_idx, start, end, temp2);
      delete[] temp2;
#endif
      return nullptr;
    }

    int axis = node.aabb.largestAxis();
    real_t min = node.aabb.min[axis], max = node.aabb.max[axis],
           inv_bin_size = Bins::BIN_COUNT / (max - min);

    Bins bins = tbb::parallel_reduce(
        tbb::blocked_range<uintp_t>(0u, size, GRAIN_SIZE), Bins(),
        [&](const tbb::blocked_range<uintp_t> &range, Bins result) {
          for (uintp_t i = range.begin(); i != range.end(); ++i) {
            uintp_t f = start[i];
            real_t centroid =
                ((1.0 / 3.0) *
                 (V[F(0, f)](axis) + V[F(1, f)](axis) + V[F(2, f)](axis)));

            int index =
                std::min(std::max((int)((centroid - min) * inv_bin_size), 0),
                         (Bins::BIN_COUNT - 1));

            result.counts[index]++;
            AABB &bin_bounds = result.bounds[index];
            bin_bounds.expandBy(V[F(0, f)]);
            bin_bounds.expandBy(V[F(1, f)]);
            bin_bounds.expandBy(V[F(2, f)]);
          }
          return result;
        },
        [](const Bins &b1, const Bins &b2) {
          Bins result;
          for (int i = 0; i < Bins::BIN_COUNT; ++i) {
            result.counts[i] = b1.counts[i] + b2.counts[i];
            result.bounds[i] = AABB::merge(b1.bounds[i], b2.bounds[i]);
          }
          return result;
        });

    AABB bounds_left[Bins::BIN_COUNT];
    bounds_left[0] = bins.bounds[0];
    for (int i = 1; i < Bins::BIN_COUNT; ++i) {
      bins.counts[i] += bins.counts[i - 1];
      bounds_left[i] = AABB::merge(bounds_left[i - 1], bins.bounds[i]);
    }
    AABB bounds_right = bins.bounds[Bins::BIN_COUNT - 1];
    int64_t best_index = -1;
    real_t best_cost = BVH::T_tri * size;
    real_t tri_factor = BVH::T_tri / node.aabb.surfaceArea();
    AABB best_bounds_right;

    for (int i = Bins::BIN_COUNT - 2; i >= 0; --i) {
      uintp_t prims_left = bins.counts[i],
              prims_right = (end - start) - bins.counts[i];
      real_t sah_cost =
          2.0 * BVH::T_aabb +
          tri_factor * (prims_left * bounds_left[i].surfaceArea() +
                        prims_right * bounds_right.surfaceArea());
      if (sah_cost < best_cost) {
        best_cost = sah_cost;
        best_index = i;
        best_bounds_right = bounds_right;
      }
      bounds_right = AABB::merge(bounds_right, bins.bounds[i]);
    }

    if (best_index == -1) {
      /* Could not find a good split plane -- retry with
         more careful serial code just to be sure.. */
#if defined(SINGLE_PRECISION)
      execute_serially(bvh, node_idx, start, end, (real_t *)temp);
#else
      real_t *temp2 = new real_t[size];
      execute_serially(bvh, node_idx, start, end, temp2);
      delete[] temp2;
#endif
      return nullptr;
    }

    uintp_t left_count = bins.counts[best_index];
    int node_idx_left = node_idx + 1;
    int node_idx_right = node_idx + 2 * left_count;

    bvh.mNodes[node_idx_left].aabb = bounds_left[best_index];
    bvh.mNodes[node_idx_right].aabb = best_bounds_right;
    node.inner.rightChild = node_idx_right;
    node.inner.unused = 0;

    std::atomic<uintp_t> offset_left(0), offset_right(bins.counts[best_index]);
    tbb::parallel_for(
        tbb::blocked_range<uintp_t>(0u, size, GRAIN_SIZE),
        [&](const tbb::blocked_range<uintp_t> &range) {
          uintp_t count_left = 0, count_right = 0;
          for (uintp_t i = range.begin(); i != range.end(); ++i) {
            uintp_t f = start[i];
            real_t centroid =
                ((1.0 / 3.0) *
                 (V[F(0, f)](axis) + V[F(1, f)](axis) + V[F(2, f)](axis)));
            int index = (int)((centroid - min) * inv_bin_size);
            (index <= best_index ? count_left : count_right)++;
          }
          uintp_t idx_l = offset_left.fetch_add(count_left);
          uintp_t idx_r = offset_right.fetch_add(count_right);
          for (uintp_t i = range.begin(); i != range.end(); ++i) {
            uintp_t f = start[i];
            real_t centroid =
                ((1.0 / 3.0) *
                 (V[F(0, f)](axis) + V[F(1, f)](axis) + V[F(2, f)](axis)));
            int index = (int)((centroid - min) * inv_bin_size);
            if (index <= best_index)
              temp[idx_l++] = f;
            else
              temp[idx_r++] = f;
          }
        });
    memcpy(start, temp, size * sizeof(uintp_t));
    assert(offset_left == left_count && offset_right == size);

    /* Create an empty parent task */
    tbb::task &c = *new (allocate_continuation()) tbb::empty_task;
    c.set_ref_count(2);

    /* Post right subtree to scheduler */
    BVHBuildTask &b = *new (c.allocate_child()) BVHBuildTask(
        bvh, node_idx_right, start + left_count, end, temp + left_count);
    spawn(b);

    /* Directly start working on left subtree */
    recycle_as_child_of(c);
    node_idx = node_idx_left;
    end = start + left_count;

    return this;
  }

  static void execute_serially(BVH &bvh, uintp_t node_idx, uintp_t *start,
                               uintp_t *end, real_t *temp) {
    uintp_t size = end - start;
    BVHNode &node = bvh.mNodes[node_idx];
    const MapXu F = MapXu(bvh.mF->data(), 3, bvh.mF->size());
    const std::vector<Vector3f> &V = *bvh.mV;
    real_t best_cost = BVH::T_tri * size;
    int64_t best_index = -1, best_axis = -1;
    real_t *left_areas = (real_t *)temp;

    for (int axis = 0; axis < 3; ++axis) {

      std::sort(start, end, [&](uintp_t f1, uintp_t f2) {
        return ((V[F(0, f1)](axis) + V[F(1, f1)](axis) + V[F(2, f1)](axis)) <
                (V[F(0, f2)](axis) + V[F(1, f2)](axis) + V[F(2, f2)](axis)));
      });

      AABB aabb;
      for (uintp_t i = 0; i < size; ++i) {
        uintp_t f = *(start + i);
        aabb.expandBy(V[F(0, f)]);
        aabb.expandBy(V[F(1, f)]);
        aabb.expandBy(V[F(2, f)]);
        left_areas[i] = (real_t)aabb.surfaceArea();
      }
      if (axis == 0)
        node.aabb = aabb;

      aabb.clear();

      real_t tri_factor = BVH::T_tri / node.aabb.surfaceArea();
      for (uintp_t i = size - 1; i >= 1; --i) {
        uintp_t f = *(start + i);

        aabb.expandBy(V[F(0, f)]);
        aabb.expandBy(V[F(1, f)]);
        aabb.expandBy(V[F(2, f)]);

        real_t left_area = left_areas[i - 1];
        real_t right_area = aabb.surfaceArea();
        uintp_t prims_left = i;
        uintp_t prims_right = size - i;

        real_t sah_cost =
            2.0 * BVH::T_aabb +
            tri_factor * (prims_left * left_area + prims_right * right_area);
        if (sah_cost < best_cost) {
          best_cost = sah_cost;
          best_index = i;
          best_axis = axis;
        }
      }
    }

    if (best_index == -1) {
      /* Splitting does not reduce the cost, make a leaf */
      node.leaf.flag = 1;
      node.leaf.start = start - bvh.mIndices.data();
      node.leaf.size = size;
      return;
    }

    std::sort(start, end, [&](uintp_t f1, uintp_t f2) {
      return ((V[F(0, f1)](best_axis) + V[F(1, f1)](best_axis) +
               V[F(2, f1)](best_axis)) <
              (V[F(0, f2)](best_axis) + V[F(1, f2)](best_axis) +
               V[F(2, f2)](best_axis)));
    });

    uintp_t left_count = best_index;
    int node_idx_left = node_idx + 1;
    int node_idx_right = node_idx + 2 * left_count;
    node.inner.rightChild = node_idx_right;
    node.inner.unused = 0;

    execute_serially(bvh, node_idx_left, start, start + left_count, temp);
    execute_serially(bvh, node_idx_right, start + left_count, end,
                     temp + left_count);
  }
};

BVH::BVH(std::vector<uint32_t> *F, std::vector<Vector3f> *V,
         std::vector<Vector3f> *N, const AABB &aabb)
    : mF(F), mV(V), mN(N) {
  mNodes.resize(2 * mF->size() / 3);
  memset(mNodes.data(), 0, sizeof(BVHNode) * mNodes.size());
  mNodes[0].aabb = aabb;
  mIndices.resize(mF->size() / 3);
}

void BVH::build(const ProgressCallback &progress) {
  if (mF->size() == 0)
    return;
  mProgress = progress;

#if defined(SINGLE_PRECISION)
  if (sizeof(BVHNode) != 32)
    throw std::runtime_error(
        "BVH Node is not packed! Investigate compiler settings.");
#endif

  if (progress) {
    cout << "Constructing Bounding Volume Hierarchy .. ";
    cout.flush();
  }

  uintp_t total_size = mF->size() / 3;

  for (uintp_t i = 0; i < total_size; ++i)
    mIndices[i] = i;

  Timer<> timer;
  uintp_t *indices = mIndices.data();
  uintp_t *temp = new uintp_t[total_size];
  BVHBuildTask &task = *new (tbb::task::allocate_root()) BVHBuildTask(
      *this, 0u, indices, indices + total_size, temp);
  tbb::task::spawn_root_and_wait(task);
  delete[] temp;

  std::pair<real_t, uintp_t> stats = statistics();

  if (progress) {
    cout << "done. ("
         << "SAH cost = " << stats.first << ", "
         << "nodes = " << stats.second << ", "
         << "took " << timeString(timer.reset()) << ")" << endl;

    cout.precision(4);
    cout << "Compressing BVH node storage to "
         << 100 * stats.second / (real_t)mNodes.size()
         << "% of its original size .. ";
    cout.flush();
  }

  std::vector<BVHNode> compressed(stats.second);
  std::vector<uintp_t> skipped_accum(mNodes.size());

  for (int64_t i = stats.second - 1, j = mNodes.size(), skipped = 0; i >= 0;
       --i) {
    while (mNodes[--j].isUnused())
      skipped++;
    BVHNode &new_node = compressed[i];
    new_node = mNodes[j];
    skipped_accum[j] = skipped;

    if (new_node.isInner()) {
      new_node.inner.rightChild =
          i + new_node.inner.rightChild - j -
          (skipped - skipped_accum[new_node.inner.rightChild]);
    }
  }

  mNodes = std::move(compressed);

  if (progress) {
    cout << "done. (took " << timeString(timer.value()) << ")" << endl;
  }

  mProgress = nullptr;
}

bool BVH::rayIntersect(Ray ray, uint32_t &idx, real_t &t, Vector2f *uv) const {
  if (mNodes.empty())
    return false;

  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;
  bool hit = false;
  t = std::numeric_limits<real_t>::infinity();

  while (true) {
    const BVHNode &node = mNodes[node_idx];

    if (!node.aabb.rayIntersect(ray)) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      stack[stack_idx++] = node.inner.rightChild;
      node_idx++;
      assert(stack_idx < 64);
    } else {
      real_t _t;
      Vector2f _uv;
      for (uintp_t i = node.start(), end = node.end(); i < end; ++i) {
        if (lineIntersectTri(ray, mIndices[i]) &&
            rayIntersectTri(ray, mIndices[i], _t, _uv)) {
          idx = mIndices[i];
          t = ray.maxt = _t;
          hit = true;
          if (uv)
            *uv = _uv;
        }
      }
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }
  }

  return hit;
}

bool BVH::planeIntersect(const Vector3f &p, const Vector3f &norm,
                         std::vector<uint32_t> &idx) const {
  if (mNodes.empty())
    return false;

  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;
  bool hit = false;

  while (true) {
    const BVHNode &node = mNodes[node_idx];

    if (!node.aabb.planeIntersect(p, norm)) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      stack[stack_idx++] = node.inner.rightChild;
      node_idx++;
      assert(stack_idx < 64);
    } else {
      for (uintp_t i = node.start(), end = node.end(); i < end; ++i) {
        if (planeIntersectTri(p, norm, mIndices[i])) {
          idx.push_back(mIndices[i]);
          hit = true;
        }
      }
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
    }
  }

  return hit;
}

bool BVH::rayIntersect(Ray ray, std::vector<uint32_t> &idx,
                       std::vector<Vector2f> *uv) const {
  if (mNodes.empty())
    return false;

  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;
  bool hit = false;

  while (true) {
    const BVHNode &node = mNodes[node_idx];

    if (!node.aabb.rayIntersect(ray)) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      stack[stack_idx++] = node.inner.rightChild;
      node_idx++;
      assert(stack_idx < 64);
    } else {
      real_t _t;
      Vector2f _uv;
      for (uintp_t i = node.start(), end = node.end(); i < end; ++i) {
        if (lineIntersectTri(ray, mIndices[i]) &&
            rayIntersectTri(ray, mIndices[i], _t, _uv)) {
          idx.push_back(mIndices[i]);
          if (uv)
            uv->push_back(_uv);
          hit = true;
        }
      }
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
    }
  }

  return hit;
}

bool BVH::rayIntersect(Ray ray, std::vector<uint32_t> &idx,
                       std::vector<real_t> &ts,
                       std::vector<Vector2f> *uv) const {
  if (mNodes.empty())
    return false;

  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;
  bool hit = false;

  while (true) {
    const BVHNode &node = mNodes[node_idx];

    if (!node.aabb.rayIntersect(ray)) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      stack[stack_idx++] = node.inner.rightChild;
      node_idx++;
      assert(stack_idx < 64);
    } else {
      real_t _t;
      Vector2f _uv;
      for (uintp_t i = node.start(), end = node.end(); i < end; ++i) {
        if (lineIntersectTri(ray, mIndices[i]) &&
            rayIntersectTri(ray, mIndices[i], _t, _uv)) {
          idx.push_back(mIndices[i]);
          ts.push_back(_t);
          if (uv)
            uv->push_back(_uv);
          hit = true;
        }
      }
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
    }
  }

  return hit;
}

bool BVH::rayIntersect(Ray ray) const {
  if (mNodes.empty())
    return false;

  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;

  while (true) {
    const BVHNode &node = mNodes[node_idx];

    if (!node.aabb.rayIntersect(ray)) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      stack[stack_idx++] = node.inner.rightChild;
      node_idx++;
      assert(stack_idx < 64);
    } else {
      real_t t;
      Vector2f uv;
      for (uintp_t i = node.start(), end = node.end(); i < end; ++i)
        if (lineIntersectTri(ray, mIndices[i]) &&
            rayIntersectTri(ray, mIndices[i], t, uv))
          return true;
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }
  }

  return false;
}

void BVH::findNearestWithRadius(const Vector3f &p, real_t radius,
                                std::vector<uint32_t> &result,
                                bool includeSelf) const {
  result.clear();

  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;
  real_t radius2 = radius * radius;

  while (true) {
    const BVHNode &node = mNodes[node_idx];
    if (node.aabb.squaredDistanceTo(p) > radius2) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      stack[stack_idx++] = node.inner.rightChild;
      node_idx++;
      assert(stack_idx < 64);
    } else {
      uintp_t start = node.leaf.start, end = start + node.leaf.size;
      for (uintp_t i = start; i < end; ++i) {
        uintp_t f = mIndices[i];
        Vector3f pointPos = Vector3f::Zero();
        for (int j = 0; j < 3; ++j)
          pointPos += mV->at((*mF)[j + f * 3]);
        pointPos *= 1.0 / 3.0;
        real_t pointDist2 = (pointPos - p).squaredNorm();
        if (pointDist2 < radius2 && (pointDist2 != 0 || includeSelf))
          result.push_back(f);
      }
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }
  }
}

uint32_t BVH::findNearest(const Vector3f &p, real_t &radius,
                          bool includeSelf) const {
  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;
  real_t radius2 = radius * radius;
  uintp_t result = (uintp_t)-1;

  while (true) {
    const BVHNode &node = mNodes[node_idx];
    if (node.aabb.squaredDistanceTo(p) > radius2) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      uintp_t left = node_idx + 1, right = node.inner.rightChild;
      real_t distLeft = mNodes[left].aabb.squaredDistanceTo(p);
      real_t distRight = mNodes[right].aabb.squaredDistanceTo(p);
      if (distLeft < distRight) {
        node_idx = left;
        if (distRight < radius2)
          stack[stack_idx++] = right;
      } else {
        node_idx = right;
        if (distLeft < radius2)
          stack[stack_idx++] = left;
      }
      assert(stack_idx < 64);
    } else {
      uintp_t start = node.leaf.start, end = start + node.leaf.size;
      for (uintp_t i = start; i < end; ++i) {
        uintp_t f = mIndices[i];
        Vector3f pointPos = Vector3f::Zero();
        for (int j = 0; j < 3; ++j)
          pointPos += mV->at((*mF)[j + f * 3]);
        pointPos *= 1.0 / 3.0;
        real_t pointDist2 = (pointPos - p).squaredNorm();

        if (pointDist2 < radius2 && (pointDist2 != 0 || includeSelf)) {
          radius2 = pointDist2;
          result = f;
        }
      }
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }
  }
  radius = std::sqrt(radius2);
  return result;
}

void BVH::findKNearest(const Vector3f &p, uint32_t k, real_t &radius,
                       std::vector<std::pair<real_t, uint32_t>> &result,
                       bool includeSelf) const {
  result.clear();

  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;
  real_t radius2 = radius * radius;
  bool isHeap = false;
  auto comp = [](const std::pair<real_t, uintp_t> &v1,
                 const std::pair<real_t, uintp_t> &v2) {
    return v1.first < v2.first;
  };

  while (true) {
    const BVHNode &node = mNodes[node_idx];
    if (node.aabb.squaredDistanceTo(p) > radius2) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      uintp_t left = node_idx + 1, right = node.inner.rightChild;
      real_t distLeft = mNodes[left].aabb.squaredDistanceTo(p);
      real_t distRight = mNodes[right].aabb.squaredDistanceTo(p);
      if (distLeft < distRight) {
        node_idx = left;
        if (distRight < radius2)
          stack[stack_idx++] = right;
      } else {
        node_idx = right;
        if (distLeft < radius2)
          stack[stack_idx++] = left;
      }
      assert(stack_idx < 64);
    } else {
      uintp_t start = node.leaf.start, end = start + node.leaf.size;
      for (uintp_t i = start; i < end; ++i) {
        uintp_t f = mIndices[i];
        Vector3f pointPos = Vector3f::Zero();
        for (int j = 0; j < 3; ++j)
          pointPos += mV->at((*mF)[j * 3 + f]);
        pointPos *= 1.0 / 3.0;
        real_t pointDist2 = (pointPos - p).squaredNorm();

        if (pointDist2 < radius2 && (pointDist2 != 0 || includeSelf)) {
          if (result.size() < k) {
            result.push_back(std::make_pair(pointDist2, f));
          } else {
            if (!isHeap) {
              /* Establish the max-heap property */
              std::make_heap(result.begin(), result.end(), comp);
              isHeap = true;
            }

            result.push_back(std::make_pair(pointDist2, f));
            std::push_heap(result.begin(), result.end(), comp);
            std::pop_heap(result.begin(), result.end(), comp);
            result.pop_back();

            /* Reduce the search radius accordingly */
            radius2 = result[0].first;
          }
        }
      }
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }
  }
  radius = std::sqrt(radius2);
}

void BVH::findKNearest(const Vector3f &p, const Vector3f &n, uint32_t k,
                       real_t &radius,
                       std::vector<std::pair<real_t, uint32_t>> &result,
                       real_t angleThresh, bool includeSelf) const {
  result.clear();

  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;
  real_t radius2 = radius * radius;
  bool isHeap = false;
  angleThresh = std::cos(angleThresh * M_PI / 180);
  auto comp = [](const std::pair<real_t, uintp_t> &v1,
                 const std::pair<real_t, uintp_t> &v2) {
    return v1.first < v2.first;
  };

  while (true) {
    const BVHNode &node = mNodes[node_idx];
    if (node.aabb.squaredDistanceTo(p) > radius2) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      uintp_t left = node_idx + 1, right = node.inner.rightChild;
      real_t distLeft = mNodes[left].aabb.squaredDistanceTo(p);
      real_t distRight = mNodes[right].aabb.squaredDistanceTo(p);
      if (distLeft < distRight) {
        node_idx = left;
        if (distRight < radius2)
          stack[stack_idx++] = right;
      } else {
        node_idx = right;
        if (distLeft < radius2)
          stack[stack_idx++] = left;
      }
      assert(stack_idx < 64);
    } else {
      uintp_t start = node.leaf.start, end = start + node.leaf.size;
      for (uintp_t i = start; i < end; ++i) {
        uintp_t f = mIndices[i];
        Vector3f pointPos = Vector3f::Zero();
        for (int j = 0; j < 3; ++j)
          pointPos += mV->at((*mF)[j + 3 * f]);
        pointPos *= 1.0 / 3.0;
        Vector3f pointNormal = Vector3f::Zero();
        for (int j = 0; j < 3; ++j)
          pointNormal += mN->at((*mF)[j + 3 * f]);
        real_t pointDist2 = (pointPos - p).squaredNorm();

        if (pointDist2 < radius2 && (pointDist2 != 0 || includeSelf) &&
            pointNormal.dot(n) > angleThresh) {
          if (result.size() < k) {
            result.push_back(std::make_pair(pointDist2, f));
          } else {
            if (!isHeap) {
              /* Establish the max-heap property */
              std::make_heap(result.begin(), result.end(), comp);
              isHeap = true;
            }

            result.push_back(std::make_pair(pointDist2, f));
            std::push_heap(result.begin(), result.end(), comp);
            std::pop_heap(result.begin(), result.end(), comp);
            result.pop_back();

            /* Reduce the search radius accordingly */
            radius2 = result[0].first;
          }
        }
      }
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }
  }
  radius = std::sqrt(radius2);
}

bool BVH::lineIntersectTri(const Ray &ray, uint32_t i) const {
  uint32_t v0 = (*mF)[i * 3];
  uint32_t v1 = (*mF)[1 + i * 3];
  uint32_t v2 = (*mF)[2 + i * 3];
  if (v0 >= mV->size() || v1 >= mV->size() || v2 >= mV->size())
    return false;
  const Vector3f &p0 = mV->at(v0), &p1 = mV->at(v1), &p2 = mV->at(v2);
  std::vector<Vector3f> vv{p0, p1, p2};

  auto shoot_in_triangle = [](std::vector<Vector3f> const &vv,
                              const Vector3f &v, const Vector3f &dir) -> bool {
    std::unordered_set<igl::predicates::Orientation> seen_orient;
    real_t norm_offset = 10;
    for (size_t i = 0; i < vv.size(); ++i) {
      size_t i_next = (i + 1) % vv.size();
      Eigen::Vector3d v1 = vv[i];
      Eigen::Vector3d v2 = vv[i_next];
      Eigen::Vector3d up = vv[i];
      up = up + norm_offset * dir;

      igl::predicates::Orientation ori =
          igl::predicates::orient3d(v1, v2, up, v);

      seen_orient.emplace(ori);
    }
    size_t side_count = 0;
    for (auto ori : seen_orient) {
      if (ori != igl::predicates::Orientation::COPLANAR)
        side_count++;
    }
    return side_count == 1;
  };

  return shoot_in_triangle(vv, ray.o, ray.d);
}

bool BVH::rayIntersectTri(const Ray &ray, uint32_t i, real_t &t,
                          Vector2f &uv) const {
  uint32_t v0 = (*mF)[i * 3];
  uint32_t v1 = (*mF)[1 + i * 3];
  uint32_t v2 = (*mF)[2 + i * 3];
  if (v0 >= mV->size() || v1 >= mV->size() || v2 >= mV->size())
    return false;
  const Vector3f &p0 = mV->at(v0), &p1 = mV->at(v1), &p2 = mV->at(v2);

  Vector3f edge1 = p1 - p0, edge2 = p2 - p0;
  Vector3f pvec = ray.d.cross(edge2);

  real_t det = edge1.dot(pvec);
  if (det == 0.0)
    return false;
  real_t inv_det = 1.0 / det;

  Vector3f tvec = ray.o - p0;
  real_t u = tvec.dot(pvec) * inv_det;
  if (u < 0.0 || u > 1.0)
    return false;

  Vector3f qvec = tvec.cross(edge1);
  real_t v = ray.d.dot(qvec) * inv_det;

  if (v < 0.0 || u + v > 1.0)
    return false;

  real_t tempT = edge2.dot(qvec) * inv_det;
  if (tempT < ray.mint || tempT > ray.maxt)
    return false;

  t = tempT;
  uv << u, v;
  return true;
}

bool BVH::planeIntersectTri(const Vector3f &pp, const Vector3f &norm,
                            uint32_t i) const {
  const Vector3f &p0 = mV->at((*mF)[i * 3]), &p1 = mV->at((*mF)[1 + 3 * i]),
                 &p2 = mV->at((*mF)[2 + 3 * i]);

  // Check if corners are on different sides / on the plane.
  Vector3f pp2 = norm.cross(Vector3f(1, 0, 0));
  Vector3f pp3 = norm.cross(pp2);
  std::vector<Vector3f> corners({p0, p1, p2});
  bool has_seen_pos = false, has_seen_neg = false;
  for (uint8_t i = 0; i < corners.size(); i++) {
    Vector3f corner = corners[i];
    auto side = igl::predicates::orient3d(pp, pp2, pp3, corner);
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

void BVH::printStatistics() const {
  cout << endl;
  cout << "Bounding Volume Hierarchy statistics:" << endl;
  cout << "    Tree nodes         : "
       << memString(sizeof(BVHNode) * mNodes.size()) << endl;
  cout << "    Index buffer       : " << memString(sizeof(uintp_t) * mF->size())
       << endl;
  cout << "    Total              : "
       << memString(sizeof(BVHNode) * mNodes.size() +
                    sizeof(uintp_t) * mF->size())
       << endl;
}

std::pair<real_t, uintp_t> BVH::statistics(uintp_t node_idx) const {
  const BVHNode &node = mNodes[node_idx];
  if (node.isLeaf()) {
    return std::make_pair(T_tri * node.leaf.size, 1u);
  } else {
    std::pair<real_t, uintp_t> stats_left = statistics(node_idx + 1u);
    std::pair<real_t, uintp_t> stats_right = statistics(node.inner.rightChild);
    real_t saLeft = mNodes[node_idx + 1u].aabb.surfaceArea();
    real_t saRight = mNodes[node.inner.rightChild].aabb.surfaceArea();
    real_t saCur = node.aabb.surfaceArea();
    real_t sahCost =
        2 * BVH::T_aabb +
        (saLeft * stats_left.first + saRight * stats_right.first) / saCur;

    return std::make_pair(sahCost, stats_left.second + stats_right.second + 1u);
  }
}

BVH::~BVH() {}
