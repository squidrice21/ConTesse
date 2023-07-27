// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

bool is_geometric_concave(Mesh const &mesh, Edge const &e);
bool is_back_facing(Mesh const &mesh, Edge const &edge, Camera const &camera);

void tag_concave_edges_mesh(Mesh &mesh, ContourMode mode);
void tag_concave_edges_side(Mesh &mesh, Camera const &camera, ContourMode mode);
void tag_concave_edges_analytic(Mesh &mesh, Camera const &camera);

// A temporary function used to get curtain fold vertices based on the
// convexity tags on faces (for the interpolated contour mode)
void get_curtain_folds_interpolated(Mesh const &mesh,
                                    std::vector<Vector3f> &curtain_folds);

void tag_cusp_given_concave_edges(Mesh &mesh, Camera const &camera);
void tag_cusp_sharp(Mesh &mesh, Camera const &camera);
void tag_cusp_given_concave_edges_walking(Mesh &mesh, Camera const &camera);
