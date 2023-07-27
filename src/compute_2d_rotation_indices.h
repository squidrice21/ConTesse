// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

real_t compute_2d_rotation_indices(Mesh const &mesh, Camera &camera,
                                   Vertex const &v, FacingType facing);
void compute_2d_rotation_indices(Mesh &mesh, Camera &camera);
