#pragma once

#include "common.h"
#include "mesh.h"
#include <unordered_set>
#include <utility>

bool planeIntersectEdge(const Vector3f &pp, const Vector3f &norm,
                        std::vector<Vector3f> const &tri_p,
                        Vector3f &intersection);
bool planeIntersectEdge(const Vector3f &pp, const Vector3f &norm,
                        Vector3f const &ev1, Vector3f const &ev2,
                        Vector3f &intersection);

void get_neighbors(Mesh const &mesh, int id, Face src_orig,
                   size_t neighborhood_size,
                   std::unordered_set<int> &neighbor_set);
void get_neighbors(Mesh const &mesh, int id, Vertex src_orig,
                   size_t neighborhood_size,
                   std::unordered_set<int> &neighbor_set);
void walk_on_original(Mesh const &mesh, int id, Face src_orig,
                      Vector3f const &src_pos, Face dst_orig,
                      Vector3f const &dst_pos, Vector3f const &cut_p,
                      Vector3f const &cut_norm, Camera const &camera,
                      real_t sampling_delta, real_t min_delta_2d,
                      std::vector<std::pair<Vector3f, int>> &samples,
                      std::vector<int> to_relax = std::vector<int>());

void inflat_path(Mesh const &mesh, Mesh const &flat_patch, Camera const &camera,
                 real_t sampling_delta, real_t min_delta_2d,
                 Mesh &inflated_patch,
                 std::vector<std::pair<int, int>> &failed_paths,
                 std::unordered_set<BoundaryType> const &to_inflat_labels =
                     std::unordered_set<BoundaryType>());
void inflat_path(Mesh const &mesh, Mesh const &flat_patch, Camera const &camera,
                 real_t sampling_delta, real_t min_delta_2d,
                 std::vector<std::pair<Vector3f, int>> &inflat_samples,
                 std::unordered_set<int> const &to_inflat_comp_labels =
                     std::unordered_set<int>());
