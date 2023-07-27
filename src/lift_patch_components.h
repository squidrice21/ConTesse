// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

void lift_patch_components(Mesh const &mesh, Camera const &camera,
                           real_t sampling_delta, real_t min_delta_2d,
                           Mesh const &decomp_blift_patch, Mesh &inflat_patch,
                           bool sparse_sampling = false);

void lift_patch_components_iterative(Mesh const &mesh, Camera const &camera,
                                     real_t sampling_delta, real_t min_delta_2d,
                                     real_t max_distance,
                                     Mesh const &decomp_blift_patch,
                                     Mesh &inflat_patch,
                                     size_t max_inflation_itr = 1);
