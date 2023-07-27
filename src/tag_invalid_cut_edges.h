// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

bool tag_invalid_cut_edges(Mesh &mesh, Camera const &camera,
                           size_t boundary_neighbor_size,
                           size_t contour_neighbor_size_1d,
                           bool use_simple_rule);
