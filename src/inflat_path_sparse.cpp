// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "inflat_path_sparse.h"

// Eigen printout formatting
#include <spdlog/fmt/ostr.h>

#include "inflat_path.h"
#include <deque>
#include <limits>
#include <vector>

real_t min_dist(Vector2f const &src_2d, Vector2f const &dst_2d,
                std::vector<Vector2f> const &inflat_samples_2d) {
  real_t result_dist = std::numeric_limits<real_t>::infinity();
  for (auto const &v : inflat_samples_2d) {
    real_t dist;
    Vector2f v0 = src_2d;
    Vector2f v1 = dst_2d;
    real_t end_dist = (src_2d - v).norm();
    real_t end_dist2 = (dst_2d - v).norm();
    if (end_dist2 < end_dist) {
      end_dist = end_dist2;
      v0 = dst_2d;
      v1 = src_2d;
    }
    dist = end_dist;

    // Check if the projection would be on the segment
    if ((v - v0).dot(v1 - v0) >= 0) {
      Vector2f tangent = v1 - v0;
      Vector2f norm(-tangent.y(), tangent.x());
      norm.normalize();
      dist = std::abs((v - v0).dot(norm));
    }

    result_dist = std::min(result_dist, dist);
  }

  return result_dist;
}

void inflat_path_sparse(Mesh const &mesh, Mesh const &flat_patch,
                        Camera const &camera, real_t sampling_delta,
                        real_t min_delta_2d,
                        std::vector<std::pair<Vector3f, int>> &inflat_samples,
                        std::unordered_set<int> const &to_inflat_comp_labels) {
  auto flat_vpositions = flat_patch.get_vertex_property<Vector3f>("v:point");
  auto orig_idx = flat_patch.get_vertex_property<int>("v:orig_idx");
  auto orig_f_idx = flat_patch.get_vertex_property<int>("v:orig_f_idx");
  contess_assert_msg(orig_idx && orig_f_idx,
                     "inflat_path_sparse: Missing index correspondence.");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  auto flat_patchID = flat_patch.get_face_property<int>("f:patchID");
  contess_assert_msg(patchID && flat_patchID && (flat_patch.n_faces() > 0),
                     "inflat_path_sparse: Error in the flat patch.");
  int id = flat_patchID[Face(0)];

  auto boundary_type =
      flat_patch.get_edge_property<BoundaryType>("e:boundary_type");
  contess_assert_msg(boundary_type,
                     "inflat_path_sparse: Boundary types are not labeled.");
  auto patch_component = flat_patch.get_face_property<int>("f:patch_component");
  if (!to_inflat_comp_labels.empty())
    contess_assert_msg(
        patch_component,
        "inflat_path_sparse: Missing decomposed patch component.");

  // Iterate the edges in a flood filling manner
  std::unordered_set<int> seen_edges;
  std::unordered_set<int> seen_faces;
  std::deque<Face> process_faces;

  // 1. Find a boundary face as the seed
  Face face_seed;
  for (int i = 0; i < (int)flat_patch.n_edges(); i++) {
    Edge flat_e(i);

    // Check boundary edges
    if (boundary_type[flat_e] == BoundaryType::NONE)
      continue;

    std::vector<Face> adj_faces;
    if (flat_patch.face(flat_e, 0).is_valid())
      adj_faces.emplace_back(flat_patch.face(flat_e, 0));
    if (flat_patch.face(flat_e, 1).is_valid())
      adj_faces.emplace_back(flat_patch.face(flat_e, 1));

    // Check if it's from the desired component
    if (!to_inflat_comp_labels.empty()) {
      for (auto const &f : adj_faces) {
        if (to_inflat_comp_labels.count(patch_component[f]))
          face_seed = f;
      }
    } else {
      if (!adj_faces.empty())
        face_seed = adj_faces.front();
    }

    if (face_seed.is_valid())
      break;
  }

  contess_assert_msg(face_seed.is_valid(),
                     "inflat_path_sparse: Can't find a seed face.");

  // 2. Grow from the seed face while sampling on edges
  Vector2i viewport({camera.vpWidth(), camera.vpHeight()});
  std::vector<Vector2f> inflat_samples_2d;
  process_faces.push_back(face_seed);
  while (!process_faces.empty()) {
    Face cur_f = process_faces.front();
    process_faces.pop_front();

    seen_faces.emplace(cur_f.idx());

    auto hit = flat_patch.halfedges(cur_f), hit_end = hit;
    do {
      // Expand face frontier
      Face adj_f = flat_patch.face(flat_patch.opposite_halfedge(*hit));
      if (adj_f.is_valid() && !seen_faces.count(adj_f.idx()) &&
          (to_inflat_comp_labels.empty() ||
           to_inflat_comp_labels.count(patch_component[adj_f])))
        process_faces.push_back(adj_f);

      // If we can consider sampling this edge
      if (boundary_type[flat_patch.edge(*hit)] != BoundaryType::NONE ||
          seen_edges.count(flat_patch.edge(*hit).idx()))
        continue;
      seen_edges.emplace(flat_patch.edge(*hit).idx());

      // Check if this edge is sufficiently far from existing samples
      Edge flat_e = flat_patch.edge(*hit);
      Vertex src_flat = flat_patch.vertex(flat_e, 0);
      Vertex dst_flat = flat_patch.vertex(flat_e, 1);

      // Determine the 2D projection
      Vector2f src_2d =
          project(flat_vpositions[src_flat], camera.viewMatrix().matrix(),
                  camera.projectionMatrix(), viewport)
              .head<2>();
      Vector2f dst_2d =
          project(flat_vpositions[dst_flat], camera.viewMatrix().matrix(),
                  camera.projectionMatrix(), viewport)
              .head<2>();

      // Trim from the end points since many edges may share the same end vertex
      // which would inevitably cause small distance near the shared vertex
      Vector2f path_dir_2d = dst_2d - src_2d;
      real_t trim_ratio = 2 * min_delta_2d / path_dir_2d.norm();
      trim_ratio = std::min(trim_ratio, 0.4);
      Vector2f src_2d_trim = src_2d + path_dir_2d * trim_ratio;
      Vector2f dst_2d_trim = src_2d + path_dir_2d * (1 - trim_ratio);

      real_t dist = min_dist(src_2d_trim, dst_2d_trim, inflat_samples_2d);
      if (dist < 0.1 * min_delta_2d)
        continue;

      // Enough distance, sample
      Vector3f path_dir = flat_vpositions[src_flat] - flat_vpositions[dst_flat];
      if (path_dir.norm() < sampling_delta)
        continue;

      Vector3f path_norm = (camera.position() - flat_vpositions[dst_flat])
                               .cross(path_dir)
                               .normalized();

      Face src_orig(orig_f_idx[src_flat]);
      Face dst_orig(orig_f_idx[dst_flat]);

      // If the src is at a vertex
      if (orig_idx[src_flat] >= 0) {
        Vertex orig_v(orig_idx[src_flat]);
        auto hit = mesh.halfedges(orig_v), hit_end = hit;
        do {
          // Check adjacent faces
          auto adj_f_h = mesh.opposite_halfedge(*hit);
          Face adj_f = mesh.face(adj_f_h);
          // In case we are at the boundary of the model
          if (adj_f.is_valid() && patchID[adj_f] == id)
            src_orig = adj_f;
        } while (++hit != hit_end);
      }
      if (orig_idx[dst_flat] >= 0) {
        Vertex orig_v(orig_idx[dst_flat]);
        auto hit = mesh.halfedges(orig_v), hit_end = hit;
        do {
          // Check adjacent faces
          auto adj_f_h = mesh.opposite_halfedge(*hit);
          Face adj_f = mesh.face(adj_f_h);
          // In case we are at the boundary of the model
          if (adj_f.is_valid() && patchID[adj_f] == id)
            dst_orig = adj_f;
        } while (++hit != hit_end);
      }

      // Walk on the original mesh
      std::vector<std::pair<Vector3f, int>> samples;
      Vector3f cut_p = camera.position();
      walk_on_original(mesh, id, src_orig, flat_vpositions[src_flat], dst_orig,
                       flat_vpositions[dst_flat], cut_p, path_norm, camera,
                       sampling_delta, min_delta_2d, samples);

      // Filter samples that is too close
      std::vector<std::pair<Vector3f, int>> filtered_samples;
      for (auto const &v3d : samples) {
        real_t dist_3d = std::numeric_limits<real_t>::infinity();
        for (auto const &v3d2 : inflat_samples) {
          dist_3d = std::min(dist_3d, (v3d.first - v3d2.first).norm());
        }
        if (dist_3d < sampling_delta)
          continue;

        real_t dist_2d = std::numeric_limits<real_t>::infinity();
        Vector2f v2d = project(v3d.first, camera.viewMatrix().matrix(),
                               camera.projectionMatrix(), viewport)
                           .head<2>();
        for (auto const &v2d2 : inflat_samples_2d) {
          dist_2d = std::min(dist_2d, (v2d - v2d2).norm());
        }
        if (dist_2d < min_delta_2d)
          continue;

        filtered_samples.emplace_back(v3d);
      }
      samples = filtered_samples;

      inflat_samples.insert(inflat_samples.end(), samples.begin(),
                            samples.end());
      for (auto const &v3d : samples) {
        Vector2f v2d = project(v3d.first, camera.viewMatrix().matrix(),
                               camera.projectionMatrix(), viewport)
                           .head<2>();
        inflat_samples_2d.emplace_back(v2d);
      }

    } while (++hit != hit_end);
  }
}

void inflat_path_iterative(
    Mesh const &mesh, Mesh const &flat_patch, Camera const &camera,
    real_t sampling_delta, real_t min_delta_2d,
    std::vector<std::pair<Vector3f, int>> &inflat_samples,
    std::unordered_set<int> const &to_inflat_comp_labels, real_t max_distance) {
  auto flat_vpositions = flat_patch.get_vertex_property<Vector3f>("v:point");
  auto orig_idx = flat_patch.get_vertex_property<int>("v:orig_idx");
  auto orig_f_idx = flat_patch.get_vertex_property<int>("v:orig_f_idx");
  contess_assert_msg(orig_idx && orig_f_idx,
                     "lift_min_surface: Missing index correspondence.");
  auto flat_patchID = flat_patch.get_face_property<int>("f:patchID");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  contess_assert_msg(patchID && flat_patchID && (flat_patch.n_faces() > 0),
                     "lift_min_surface: Error in the flat patch.");
  int id = flat_patchID[Face(0)];

  auto boundary_type =
      flat_patch.get_edge_property<BoundaryType>("e:boundary_type");
  contess_assert_msg(boundary_type,
                     "lift_min_surface: Boundary types are not labeled.");
  auto patch_component = flat_patch.get_face_property<int>("f:patch_component");
  if (!to_inflat_comp_labels.empty())
    contess_assert_msg(
        patch_component,
        "lift_patch_components: Missing decomposed patch component.");

  // Iterate on flat edges
  Vector2i viewport({camera.vpWidth(), camera.vpHeight()});
  std::vector<Vector2f> inflat_samples_2d;
  for (int i = 0; i < (int)flat_patch.n_edges(); i++) {
    Edge flat_e(i);

    // Check if we want to inflat it
    // No boundary points
    if (!(boundary_type[flat_e] == BoundaryType::NONE)) {
      continue;
    }

    // Check if it's from the desired component
    if (!to_inflat_comp_labels.empty() &&
        ((!flat_patch.face(flat_e, 0).is_valid() ||
          !to_inflat_comp_labels.count(
              patch_component[flat_patch.face(flat_e, 0)])) ||
         (!flat_patch.face(flat_e, 1).is_valid() ||
          !to_inflat_comp_labels.count(
              patch_component[flat_patch.face(flat_e, 1)])))) {
      continue;
    }

    Vertex src_flat = flat_patch.vertex(flat_e, 0);
    Vertex dst_flat = flat_patch.vertex(flat_e, 1);
    Vector3f path_dir = flat_vpositions[src_flat] - flat_vpositions[dst_flat];

    Face src_orig(orig_f_idx[src_flat]);
    Face dst_orig(orig_f_idx[dst_flat]);

    if (path_dir.norm() < sampling_delta) {
      continue;
    }

    Vector3f path_norm = (camera.position() - flat_vpositions[dst_flat])
                             .cross(path_dir)
                             .normalized();

    if (src_orig == dst_orig) {
      continue;
    }

    // If the src is at a vertex
    if (orig_idx[src_flat] >= 0) {
      Vertex orig_v(orig_idx[src_flat]);
      auto hit = mesh.halfedges(orig_v), hit_end = hit;
      do {
        // Check adjacent faces
        auto adj_f_h = mesh.opposite_halfedge(*hit);
        Face adj_f = mesh.face(adj_f_h);
        // In case we are at the boundary of the model
        if (adj_f.is_valid() && patchID[adj_f] == id)
          src_orig = adj_f;
      } while (++hit != hit_end);
    }
    if (orig_idx[dst_flat] >= 0) {
      Vertex orig_v(orig_idx[dst_flat]);
      auto hit = mesh.halfedges(orig_v), hit_end = hit;
      do {
        // Check adjacent faces
        auto adj_f_h = mesh.opposite_halfedge(*hit);
        Face adj_f = mesh.face(adj_f_h);
        // In case we are at the boundary of the model
        if (adj_f.is_valid() && patchID[adj_f] == id)
          dst_orig = adj_f;
      } while (++hit != hit_end);
    }

    // Walk on the original mesh
    std::vector<std::pair<Vector3f, int>> samples;
    Vector3f cut_p = camera.position();
    walk_on_original(mesh, id, src_orig, flat_vpositions[src_flat], dst_orig,
                     flat_vpositions[dst_flat], cut_p, path_norm, camera,
                     sampling_delta, min_delta_2d, samples);

    // Max distance
    bool is_far_enough = false;
    real_t far_dist = -1;
    for (auto const &v3d : samples) {
      Vector3f tangent =
          (flat_vpositions[dst_flat] - flat_vpositions[src_flat]).normalized();
      Vector3f proj =
          (v3d.first - flat_vpositions[src_flat]).dot(tangent) * tangent;

      real_t d =
          std::abs(((v3d.first - flat_vpositions[src_flat]) - proj).norm());
      far_dist = std::max(far_dist, d);
      if (far_dist > max_distance) {
        is_far_enough = true;
        break;
      }
    }

    // If this path is sufficiently sampled
    if (!is_far_enough) {
      continue;
    }

    // Filter samples that is too close
    std::vector<std::pair<Vector3f, int>> filtered_samples;
    std::vector<Vector2f> cur_inflat_samples_2d;
    for (auto const &v3d : samples) {
      real_t dist_3d = std::numeric_limits<real_t>::infinity();
      for (auto const &v3d2 : inflat_samples) {
        dist_3d = std::min(dist_3d, (v3d.first - v3d2.first).norm());
      }
      if (dist_3d < sampling_delta)
        continue;

      real_t dist_2d = std::numeric_limits<real_t>::infinity();
      Vector2f v2d = project(v3d.first, camera.viewMatrix().matrix(),
                             camera.projectionMatrix(), viewport)
                         .head<2>();
      for (auto const &v2d2 : inflat_samples_2d) {
        dist_2d = std::min(dist_2d, (v2d - v2d2).norm());
      }
      if (dist_2d < min_delta_2d)
        continue;
      else
        cur_inflat_samples_2d.emplace_back(v2d);

      filtered_samples.emplace_back(v3d);
    }
    samples = filtered_samples;

    inflat_samples.insert(inflat_samples.end(), samples.begin(), samples.end());
    for (auto const &v2d : cur_inflat_samples_2d) {
      inflat_samples_2d.emplace_back(v2d);
    }
  }
}
