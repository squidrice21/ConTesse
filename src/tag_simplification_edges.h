// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

bool is_same_patch(Mesh const &mesh, Halfedge he, Halfedge pair_he,
                   Halfedge hh);

Halfedge find_nondegenerated_contour_halfedge(Mesh const &mesh, Vertex const &v,
                                              Halfedge const &forward);
bool tag_simplification_edges(Mesh &mesh, Camera const &camera);
bool is_boundary_joint(Mesh const &mesh, Vertex const &v);
