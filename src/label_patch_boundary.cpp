// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "label_patch_boundary.h"
#include "common.h"

// Eigen printout formatting
#include <limits>
#include <spdlog/fmt/ostr.h>

void label_patch_boundary(Mesh const &mesh, Mesh &flat_patch) {
  auto cut = mesh.get_edge_property<int>("e:disk_cut");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  contess_assert_msg(
      patchID && cut,
      "label_patch_boundary: Patch or cut information is missing.");
  auto orig_idx = flat_patch.get_vertex_property<int>("v:orig_idx");
  contess_assert_msg(orig_idx,
                     "label_patch_boundary: Correspondence is missing.");

  auto boundary_type = flat_patch.edge_property<BoundaryType>(
      "e:boundary_type", BoundaryType::NONE);
  boundary_type.vector().assign(boundary_type.vector().size(),
                                BoundaryType::NONE);

  for (size_t i = 0; i < flat_patch.n_edges(); i++) {
    Edge e_flat(i);
    Vertex v1_flat = flat_patch.vertex(e_flat, 0),
           v2_flat = flat_patch.vertex(e_flat, 1);

    Vertex v1_orig(orig_idx[v1_flat]), v2_orig(orig_idx[v2_flat]);
    contess_assert_msg(v1_orig.is_valid() && v2_orig.is_valid(),
                       "label_patch_boundary: Error in 2D-3D correspondence.");

    // 1. Label interior cuts
    Halfedge he_orig = mesh.find_halfedge(v1_orig, v2_orig);
    if (!he_orig.is_valid())
      continue;

    if (flat_patch.is_boundary(e_flat)) {
      boundary_type[e_flat] = BoundaryType::CUT;
    }

    // 2. Patch boundary
    int pid_1 =
        (mesh.face(he_orig).is_valid()) ? patchID[mesh.face(he_orig)] : -1;
    int pid_2 = (mesh.face(mesh.opposite_halfedge(he_orig)).is_valid())
                    ? patchID[mesh.face(mesh.opposite_halfedge(he_orig))]
                    : -1;
    if (pid_1 != pid_2) {
      boundary_type[e_flat] = BoundaryType::PATCH;
    }
  }
}

void get_interior_samples(
    Mesh const &mesh, Camera const &camera, Mesh const &flat_patch,
    std::vector<std::pair<Vector3f, int>> &interior_points) {
  auto flat_vpositions = flat_patch.get_vertex_property<Vector3f>("v:point");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  contess_assert_msg(patchID,
                     "get_interior_samples: Patch information is missing.");
  auto orig_idx = flat_patch.get_vertex_property<int>("v:orig_idx");
  auto orig_f_idx = flat_patch.get_vertex_property<int>("v:orig_f_idx");
  contess_assert_msg(orig_idx && orig_f_idx,
                     "label_patch_boundary: Correspondence is missing.");
  interior_points.clear();
  interior_points.reserve(flat_patch.n_vertices());

  std::vector<Vector2f> boundary_points;
  std::vector<std::pair<Vector3f, int>> dup_interior_points;
  boundary_points.reserve(flat_patch.n_vertices());
  dup_interior_points.reserve(flat_patch.n_vertices());

  for (size_t i = 0; i < flat_patch.n_vertices(); i++) {
    Vertex v(i);

    bool is_interior = true;
    auto hit = flat_patch.halfedges(v), hit_end = hit;
    do {
      Vertex v1_orig(orig_idx[flat_patch.from_vertex(*hit)]),
          v2_orig(orig_idx[flat_patch.to_vertex(*hit)]);
      Halfedge he_orig = mesh.find_halfedge(v1_orig, v2_orig);
      if (!he_orig.is_valid())
        continue;

      if (patchID[mesh.face(he_orig)] !=
          patchID[mesh.face(mesh.opposite_halfedge(he_orig))]) {
        is_interior = false;
        break;
      }
    } while (++hit != hit_end);

    Vector3f pos = flat_vpositions[v];

    if (is_interior)
      dup_interior_points.emplace_back(pos, orig_f_idx[v]);
  }

  std::vector<Vector2f> dup_interior_points_2d;
  for (auto const &v : dup_interior_points) {
    interior_points.emplace_back(v);
  }
}