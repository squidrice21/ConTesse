// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

bool can_flip_edge(Mesh const &mesh, Edge const &flip_e);
bool flip_edge(Mesh &mesh, Edge const &flip_e);
void fix_flipped_faces(Mesh &mesh, Camera const &camera);
