// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"
#include <limits>
#include <unordered_set>
#include <vector>

void get_wso_failure_removal(
    Mesh &mesh, Camera const &camera,
    std::unordered_set<int> const &wso_failures,
    std::vector<std::vector<int>> &patch_removals,
    real_t max_loop_length = std::numeric_limits<real_t>::infinity(),
    bool ff_only = false);

bool remove_wso_failure(Mesh &mesh, Camera const &camera,
                        std::vector<int> const &patch_removals,
                        std::unordered_set<int> &affected_patches,
                        bool ff_only = false);
