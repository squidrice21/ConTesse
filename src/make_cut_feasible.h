// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

Vertex split_edge(Mesh &mesh, Camera const &camera, Halfedge const &h,
                  Vector3f const &split_p);
bool is_on_correct_side(Mesh const &mesh, Camera const &camera,
                        Halfedge const &c_h1, Halfedge const &c_h2,
                        Vector3f const &p);
void make_cut_feasible(Mesh &mesh, Camera const &camera, bool ff_only = false);
