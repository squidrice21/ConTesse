// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"
#include <vector>

bool is_patch_removeable(
    Mesh const &mesh, Camera const &camera, Halfedge const &in_he,
    real_t max_loop_length = std::numeric_limits<real_t>::infinity());
void remove_patch(Mesh &mesh, Halfedge const &in_he);

void check_cut_feasibility(Mesh &mesh, std::vector<Edge> &invalid_cut_edges);
bool remove_cut_infeasible_disk(Mesh &mesh, Camera const &camera);
