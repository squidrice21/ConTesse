// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

bool tag_simplification_edges_flipped_small_loops(Mesh &mesh,
                                                  Camera const &camera,
                                                  real_t loop_epsilon,
                                                  bool ff_only = false);
