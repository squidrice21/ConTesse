// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"
#include <unordered_set>

void to_2d_mesh(Mesh const &mesh_3d, Eigen::MatrixXd const &V,
                Eigen::MatrixXi const &F,
                std::vector<size_t> const &V_idx_original, Mesh &mesh_out);

void triangulate_wso_patches(
    Mesh const &mesh, Camera const &camera,
    std::map<size_t, Mesh> &patch_triangulations, bool do_refine = false,
    std::unordered_set<int> selected_patches = std::unordered_set<int>(),
    bool ff_only = false);

void triangulate_wso_patches_adaptive(
    Mesh &mesh, Camera const &camera,
    std::map<size_t, Mesh> &patch_triangulations, bool do_refine = false,
    std::unordered_set<int> selected_patches = std::unordered_set<int>());

void triangulate_wso_patches_simplified(
    Mesh &mesh, Camera const &camera,
    std::map<size_t, Mesh> &patch_triangulations, bool do_refine = false,
    std::unordered_set<int> selected_patches = std::unordered_set<int>());
