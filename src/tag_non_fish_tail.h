// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

size_t determine_turn_side(Mesh const &mesh, Camera const &camera,
                           std::shared_ptr<Chain> chain, size_t he_outward);

void tag_non_fish_tail(Mesh &mesh, Camera const &camera, bool ff_only = false);
