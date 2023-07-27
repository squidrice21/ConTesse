// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"
#include <limits>
#include <unordered_set>

void inflat_path_sparse(Mesh const &mesh, Mesh const &flat_patch,
                        Camera const &camera, real_t sampling_delta,
                        real_t min_delta_2d,
                        std::vector<std::pair<Vector3f, int>> &inflat_samples,
                        std::unordered_set<int> const &to_inflat_comp_labels =
                            std::unordered_set<int>());

void inflat_path_iterative(
    Mesh const &mesh, Mesh const &flat_patch, Camera const &camera,
    real_t sampling_delta, real_t min_delta_2d,
    std::vector<std::pair<Vector3f, int>> &inflat_samples,
    std::unordered_set<int> const &to_inflat_comp_labels =
        std::unordered_set<int>(),
    real_t max_distance = std::numeric_limits<real_t>::infinity());
