// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

real_t tag_cusp_line_case(Mesh &mesh, Camera const &camera,
                          Vertex const &cusp_v, Vertex const &added_v,
                          Halfedge const &he_cusp, Halfedge const &he_line,
                          bool to_write = true);

Vertex add_simplification_vertex(Mesh &mesh, Camera const &camera,
                                 Vertex const &cusp_v, Halfedge const &he_cusp,
                                 Halfedge const &he_line, Vector3f &moved_pos);
