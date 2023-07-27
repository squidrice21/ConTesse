#pragma once

#include "mesh.h"

bool tag_simplification_edges_case1(
    Mesh const &mesh, Camera const &camera, Vertex v, Mesh &mesh_labeled,
    real_t &walk_length1, real_t &chain_walk_length1, real_t &walk_length2,
    real_t &chain_walk_length2, size_t &case1_count, size_t &case2_count);

bool tag_simplification_edges_case1_alt(
    Mesh const &mesh, Camera const &camera, Vertex v, Mesh &mesh_labeled,
    real_t &walk_length1, real_t &chain_walk_length1, real_t &walk_length2,
    real_t &chain_walk_length2, size_t &case1_count, size_t &case2_count);

bool tag_simplification_edges_case3(
    Mesh const &mesh, Camera const &camera, Vertex v, Mesh &mesh_labeled,
    std::unordered_map<int, int> &cusp_projections,
    std::unordered_map<int, Vector3f> &moved_proj_positions,
    real_t &furthest_walk_length_min, real_t &distance);
