// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"
#include <vector>

void insert_subdivision_cusps(Mesh &mesh, Subdiv const &subdiv,
                              Camera const &camera);
bool face_contains_cusp(Mesh &mesh, Subdiv const &subdiv, Camera const &camera,
                        Face const &f, bool fast_check = false);
void insert_subdivision_cusp_edges(Mesh &mesh, Subdiv const &subdiv,
                                   Camera const &camera);
