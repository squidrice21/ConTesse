// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

void determine_t_point_2d(Mesh const &mesh, Camera const &camera,
                          Halfedge const &he, real_t t, real_t &mid_t);

void subdivide_contour_edges_even(Mesh &mesh, Camera const &camera,
                                  std::vector<Edge> &subdiv_edges);

void subdivide_contour_intersection_edges_even(Mesh &mesh, Camera const &camera,
                                               Vertex const &intersection_v,
                                               size_t subdiv_times = 1);
