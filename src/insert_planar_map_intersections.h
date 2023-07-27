// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"
#include <vector>

void edges_intersection(
    Mesh const &mesh, Camera const &camera, std::vector<Edge> const &edges,
    std::map<int, std::vector<std::pair<real_t, int>>> &intersections);

void insert_planar_map_intersections(Mesh &mesh, Camera const &camera,
                                     bool to_match_endpoints = false);
