// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "common.h"
#include "mesh.h"

// Helper functions shared with ray_cast_QI_pairing
void first_hit_facing_types(Mesh const &mesh, real_t mid_point_t,
                            std::vector<uint32_t> &indices,
                            std::vector<real_t> &ts, FacingType &vbo,
                            FacingType &vbo_f);
void trim_occluders(Mesh const &mesh, Edge const &e, real_t mid_point_t,
                    std::vector<uint32_t> &indices, std::vector<real_t> &ts);
void get_nearer_adjacent_face(Mesh const &mesh, Camera const &camera,
                              Edge const &edge, Face &nearer_f,
                              Vertex &nearer_v);
bool is_nearer_adjacent_face_consistent(Mesh const &mesh, Camera const &camera,
                                        Edge const &edge);
bool planeIntersectTri(const Vector3f &pp, const Vector3f &norm,
                       std::vector<Vector3f> const &tri_p);

void ray_cast_QI(Mesh &mesh, Camera const &camera, bool is_trivial = false);
