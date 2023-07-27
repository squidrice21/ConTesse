// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

void tag_subdivision_contour_edges(Mesh &mesh, Camera const &camera,
                                   real_t loop_distance_threshold,
                                   real_t min_2d_length);
