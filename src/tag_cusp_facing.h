// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

#include <unordered_set>

bool point_to_triangle(Vector3f const &p, Vector3f const &a, Vector3f const &b,
                       Vector3f const &c);
bool point_to_triangle(Vector3f const &p, Vector3f const &a, Vector3f const &b,
                       Vector3f const &c, Vector3f const &norm);
Halfedge
find_nondegenerated_contour_halfedge(Mesh const &mesh, Vertex const &v,
                                     std::unordered_set<int> &visited_v);
Halfedge
find_nondegenerated_contour_halfedge_2d(Mesh const &mesh, Camera const &camera,
                                        Vertex const &v,
                                        std::unordered_set<int> &visited_v);
FacingType get_cusp_facing(Mesh &mesh, Camera const &camera,
                           Vertex const &center);
bool is_potential_cusp(Mesh const &mesh, Camera const &camera, Vertex const &v);
bool compute_cusp_laplacian(Mesh &mesh);
bool tag_cusp_facing(Mesh &mesh, Camera const &camera);
void tag_extraordinary(Mesh &mesh);
void tag_extraordinary_quad(Mesh &mesh);
void tag_extraordinary_quad(Mesh const &mesh,
                            std::vector<Vertex> &extraordinary_v);
