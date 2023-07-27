// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"
#include <vector>

void update_patch(
    Mesh &mesh, Camera const &camera,
    std::vector<Face> const &unassigned_faces = std::vector<Face>());
void subdivide_contour_edges(Mesh &mesh, Camera const &camera,
                             std::vector<Edge> &subdiv_edges);
void subdivide_contour_edges(Mesh &mesh, Camera const &camera);

void subdivide_contour_edges_k_times(Mesh &mesh, Camera const &camera,
                                     std::vector<Edge> &subdiv_edges,
                                     size_t subdiv_times = 1,
                                     bool is_even = false,
                                     real_t min_2d_length = -1);
void subdivide_contour_edges_tagged(Mesh &mesh, Camera const &camera,
                                    size_t subdiv_times = 1,
                                    bool is_even = false,
                                    real_t min_2d_length = -1);
