// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once
#pragma once

#include "mesh.h"
#include <vector>

void handle_degenerated_component_facing(Mesh const &mesh,
                                         Mesh const &flat_patch,
                                         Mesh &inflated_patch);

void cdt_non_SO_component(
    Mesh const &mesh, Mesh const &flat_patch, Camera const &camera,
    std::vector<std::pair<Vector3f, int>> const &interior_points,
    Mesh &inflated_patch, bool to_flip = true);
