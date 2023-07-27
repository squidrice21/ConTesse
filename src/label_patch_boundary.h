// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

void label_patch_boundary(Mesh const &mesh, Mesh &flat_patch);

void get_interior_samples(
    Mesh const &mesh, Camera const &camera, Mesh const &flat_patch,
    std::vector<std::pair<Vector3f, int>> &interior_points);
