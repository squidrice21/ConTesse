// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

bool collapse_edges(Mesh &mesh, Camera const &camera, Vertex const &v,
                    Halfedge const &h1, Halfedge const &h2);
bool collapse_flipped_edges(Mesh &mesh, Camera const &camera,
                            real_t parallel_offset = 1e-3);
