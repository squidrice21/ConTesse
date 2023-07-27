// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "inflat_path.h"

// Eigen printout formatting
#include <limits>
#include <spdlog/fmt/ostr.h>

#include "common.h"
#include "ray_cast_QI.h"

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

bool planeIntersectEdge(const Vector3f &pp, const Vector3f &norm,
                        std::vector<Vector3f> const &tri_p,
                        Vector3f &intersection) {
  Vector3f v1, v2;
  real_t t = -1;
  // Check if corners are on different sides / on the plane.
  Vector3f pp2 = norm.cross(Vector3f(1, 0, 0));
  Vector3f pp3 = norm.cross(pp2);
  pp2 = pp + pp2;
  pp3 = pp + pp3;

  bool has_seen_pos = false, has_seen_neg = false, is_collinear = false;
  for (uint8_t i = 0; i < tri_p.size(); i++) {
    Vector3f corner = tri_p[i];
    auto side = igl::predicates::orient3d(pp, pp2, pp3, corner);
    switch (side) {
    case igl::predicates::Orientation::NEGATIVE:
      has_seen_neg = true;
      v2 = corner;
      break;
    case igl::predicates::Orientation::POSITIVE:
      has_seen_pos = true;
      v1 = corner;
      break;
    default: // Orientation::COPLANAR
      has_seen_neg = true;
      has_seen_pos = true;
      is_collinear = true;
      v1 = corner;
      break;
    }

    if (has_seen_pos && has_seen_neg) {
      break;
    }
  }

  // Not intersecting
  if (!(has_seen_pos && has_seen_neg))
    return false;

  // Intersecting at a corner
  if (is_collinear) {
    intersection = v1;
    return true;
  }

  // Compute t
  Vector3f r = (v2 - v1).normalized();
  t = (pp - v1).dot(norm) / r.dot(norm);
  t = std::max(0., t);
  intersection = v1 + t * r;

  return true;
}

bool planeIntersectEdge(const Vector3f &pp, const Vector3f &norm,
                        Vector3f const &ev1, Vector3f const &ev2,
                        Vector3f &intersection) {

  std::vector<Vector3f> tri_p({ev1, ev2});

  return planeIntersectEdge(pp, norm, tri_p, intersection);
}

void get_neighbors(Mesh const &mesh, int id, Face src_orig,
                   size_t neighborhood_size,
                   std::unordered_set<int> &neighbor_set) {
  auto patchID = mesh.get_face_property<int>("f:patchID");
  std::unordered_map<int, size_t> face_dist;
  std::queue<Face> queue;
  if (src_orig.is_valid()) {
    queue.push(src_orig);
    face_dist[src_orig.idx()] = 0;
  }
  while (!queue.empty()) {
    Face current_f = queue.front();
    queue.pop();

    if (face_dist[current_f.idx()] > neighborhood_size)
      continue;

    neighbor_set.emplace(current_f.idx());
    // add adjacent faces if no cut boundary or contour is crossed
    // And book keep the corresponding source contour edge
    auto hit = mesh.halfedges(current_f), hend = hit;
    do {
      // adjacent_f is in the same patch
      Face adjacent_f = mesh.face(mesh.opposite_halfedge(*hit));
      if (id >= 0 && patchID[adjacent_f] != id)
        continue;
      if (adjacent_f.is_valid() && !face_dist.count(adjacent_f.idx())) {
        face_dist[adjacent_f.idx()] = face_dist[current_f.idx()] + 1;
        queue.push(adjacent_f);
      }
    } while (++hit != hend);
  }
}

void get_neighbors(Mesh const &mesh, int id, Vertex src_orig,
                   size_t neighborhood_size,
                   std::unordered_set<int> &neighbor_set) {
  auto patchID = mesh.get_face_property<int>("f:patchID");
  std::unordered_map<int, size_t> face_dist;
  std::queue<Face> queue;
  auto hit = mesh.halfedges(src_orig), hit_end = hit;
  do {
    Face f = mesh.face(*hit);
    if (f.is_valid()) {
      queue.push(f);
      face_dist[f.idx()] = 1;
    }
  } while (++hit != hit_end);

  while (!queue.empty()) {
    Face current_f = queue.front();
    queue.pop();

    if (face_dist[current_f.idx()] > neighborhood_size)
      continue;

    neighbor_set.emplace(current_f.idx());
    // add adjacent faces if no cut boundary or contour is crossed
    // And book keep the corresponding source contour edge
    auto hit = mesh.halfedges(current_f), hend = hit;
    do {
      // adjacent_f is in the same patch
      Face adjacent_f = mesh.face(mesh.opposite_halfedge(*hit));
      if (id >= 0 && patchID[adjacent_f] != id)
        continue;
      if (adjacent_f.is_valid() && !face_dist.count(adjacent_f.idx())) {
        face_dist[adjacent_f.idx()] = face_dist[current_f.idx()] + 1;
        queue.push(adjacent_f);
      }
    } while (++hit != hend);
  }
}

void walk_on_original(Mesh const &mesh, int id, Face src_orig,
                      Vector3f const &src_pos, Face dst_orig,
                      Vector3f const &dst_pos, Vector3f const &cut_p,
                      Vector3f const &cut_norm, Camera const &camera,
                      real_t sampling_delta, real_t min_delta_2d,
                      std::vector<std::pair<Vector3f, int>> &samples,
                      std::vector<int> to_relax) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto disk_cut = mesh.get_edge_property<int>("e:disk_cut");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  contess_assert_msg(patchID, "inflat_path: Extracted patches is needed.");

  Vector2i viewport({camera.vpWidth(), camera.vpHeight()});
  std::unordered_set<int> visited_face_indices;
  std::queue<Face> face_frontier;
  std::unordered_map<int, int> prev_idx;

  // Initialize with the src face
  std::unordered_set<int> src_relax;
  std::vector<Face> src_faces({src_orig});
  Vertex src_v;
  {
    Vertex vertices[3];
    mesh.verticesOfFace(src_orig, vertices);
    for (auto const &vv : vertices) {
      auto p = vpositions[vv];
      real_t d = (p - src_pos).norm();
      if (d < 1e-5) {
        src_v = vv;
        break;
      }
    }
  }
  if (src_v.is_valid()) {
    auto hit = mesh.halfedges(src_v), hit_end = hit;
    do {
      // Check adjacent faces
      auto adj_f_h = mesh.opposite_halfedge(*hit);
      Face adj_f = mesh.face(adj_f_h);
      // In case we are at the boundary of the model
      if (adj_f.is_valid() && patchID[adj_f] == id)
        src_faces.emplace_back(adj_f);
    } while (++hit != hit_end);
  }
  if (to_relax.size() == 2 && to_relax[0] > 0) {
    get_neighbors(mesh, id, src_orig, to_relax[0], src_relax);
  }
  for (auto f : src_faces) {
    if (patchID[f] != id)
      continue;

    face_frontier.push(f);
    if (f != src_orig)
      prev_idx[f.idx()] = src_orig.idx();
    if (to_relax.size() == 2 && to_relax[0] > 0) {
      src_relax.emplace(f.idx());
    }
  }

  if (face_frontier.empty())
    return;

  std::unordered_set<int> dst_relax;
  if (to_relax.size() == 2 && to_relax[1] > 0) {
    get_neighbors(mesh, id, dst_orig, to_relax[1], dst_relax);
  }

  auto is_destination = [&](Face f) -> bool {
    if (f == dst_orig || dst_relax.count(f.idx()))
      return true;

    Vertex vertices[3];
    mesh.verticesOfFace(f, vertices);
    std::vector<Vector3f> tri_p({vpositions[vertices[0]],
                                 vpositions[vertices[1]],
                                 vpositions[vertices[2]]});
    for (auto const &p : tri_p) {
      real_t d = (p - dst_pos).norm();
      if (d < 1e-5)
        return true;
    }

    return false;
  };

  auto is_on_same_size = [&](Vector3f const &pp, Vector3f const &v1,
                             Vector3f const &v2) -> bool {
    Vector3f pp2 = camera.position();
    Vector3f pp3 = pp + cut_norm;
    if ((igl::predicates::orient3d(pp, pp2, pp3, v1) ==
         igl::predicates::Orientation::POSITIVE) ==
        (igl::predicates::orient3d(pp, pp2, pp3, v2) ==
         igl::predicates::Orientation::POSITIVE))
      return true;
    return false;
  };

  // Walk
  int final_f = -1;
  std::unordered_map<int, Vector3f> f_intersection;
  do {
    Face cur_f = face_frontier.front();
    face_frontier.pop();

    // Visited?
    if (visited_face_indices.find(cur_f.idx()) != visited_face_indices.end())
      continue;

    // Is this face intersecting with the cut plane?
    Vertex vertices[3];
    mesh.verticesOfFace(cur_f, vertices);
    std::vector<Vector3f> tri_p({vpositions[vertices[0]],
                                 vpositions[vertices[1]],
                                 vpositions[vertices[2]]});
    // March along the intersected triangles
    if (planeIntersectTri(cut_p, cut_norm, tri_p) ||
        src_relax.count(cur_f.idx())) {
      // Extend frontier
      auto fhit = mesh.halfedges(cur_f), fhit_end = fhit;
      do {
        // Sample on the edges of the current face
        // and determine the intersecting edge for the next step
        Vector3f intersection;
        bool edge_intersected = planeIntersectEdge(
            cut_p, cut_norm, vpositions[mesh.from_vertex(*fhit)],
            vpositions[mesh.to_vertex(*fhit)], intersection);

        auto adj_f_h = mesh.opposite_halfedge(*fhit);
        Face adj_f = mesh.face(adj_f_h);

        // Don't grow the non-intersecting edge
        if (!edge_intersected)
          continue;

        f_intersection[cur_f.idx()] = intersection;

        // We can't cross the contour edge
        // And cut edge
        if (is_contour[mesh.edge(*fhit)] >= 0 ||
            disk_cut[mesh.edge(*fhit)] >= 0)
          continue;

        // Check adjacent faces
        if (visited_face_indices.find(adj_f.idx()) !=
            visited_face_indices.end())
          continue;

        // In case we are at the boundary of the model
        if (adj_f.is_valid()) {
          face_frontier.push(adj_f);
          prev_idx[adj_f.idx()] = cur_f.idx();
        }
      } while (++fhit != fhit_end);

      // If this face the destination?
      if (is_destination(cur_f)) {
        final_f = cur_f.idx();
        break;
      }
    }

    visited_face_indices.insert(cur_f.idx());
  } while (!face_frontier.empty());

  // Enough distance and on the correct side
  int trace_f = final_f;
  samples.clear();
  samples.emplace_back(std::make_pair(dst_pos, dst_orig.idx()));
  do {
    if (!f_intersection.count(trace_f))
      continue;
    Vector3f intersection = f_intersection[trace_f];
    Vector2f intersection_2d =
        project(intersection, camera.viewMatrix().matrix(),
                camera.projectionMatrix(), viewport)
            .head<2>();
    Vector2f back_2d =
        project(samples.back().first, camera.viewMatrix().matrix(),
                camera.projectionMatrix(), viewport)
            .head<2>();
    if ((intersection - samples.back().first).norm() > sampling_delta &&
        // !is_on_same_size(intersection, src_pos, dst_pos) &&
        !is_on_same_size(intersection, src_pos, samples.back().first) &&
        (intersection_2d - back_2d).norm() > min_delta_2d) {
      samples.emplace_back(std::make_pair(intersection, trace_f));
    }
    if (prev_idx.count(trace_f))
      trace_f = prev_idx[trace_f];
  } while (prev_idx.count(trace_f));

  std::reverse(samples.begin(), samples.end());

  // Add the destination
  Vector2f back_2d = project(samples.back().first, camera.viewMatrix().matrix(),
                             camera.projectionMatrix(), viewport)
                         .head<2>();
  Vector2f front_2d =
      project(samples.front().first, camera.viewMatrix().matrix(),
              camera.projectionMatrix(), viewport)
          .head<2>();
  Vector2f src_2d = project(src_pos, camera.viewMatrix().matrix(),
                            camera.projectionMatrix(), viewport)
                        .head<2>();
  Vector2f dst_2d = project(dst_pos, camera.viewMatrix().matrix(),
                            camera.projectionMatrix(), viewport)
                        .head<2>();
  if (!samples.empty() &&
      ((dst_pos - samples.back().first).norm() < 0.5 * sampling_delta ||
       (src_2d - front_2d).norm() < min_delta_2d)) {
    samples.pop_back();
  }
  if (!samples.empty() &&
      ((src_pos - samples.front().first).norm() < 0.5 * sampling_delta ||
       (dst_2d - back_2d).norm() < min_delta_2d)) {
    samples.erase(samples.begin());
  }
}

void inflat_path(Mesh const &mesh, Mesh const &flat_patch, Camera const &camera,
                 real_t sampling_delta, real_t min_delta_2d,
                 Mesh &inflated_patch,
                 std::vector<std::pair<int, int>> &failed_paths,
                 std::unordered_set<BoundaryType> const &to_inflat_labels) {
  // Init result mesh
  inflated_patch = flat_patch;
  inflated_patch.m_fedges.clear();
  inflated_patch.m_svertices.clear();

  //
  auto flat_vpositions = flat_patch.get_vertex_property<Vector3f>("v:point");
  auto orig_idx = flat_patch.get_vertex_property<int>("v:orig_idx");
  auto orig_f_idx = flat_patch.get_vertex_property<int>("v:orig_f_idx");
  contess_assert_msg(orig_idx && orig_f_idx,
                     "lift_min_surface: Missing index correspondence.");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  auto flat_patchID = flat_patch.get_face_property<int>("f:patchID");
  contess_assert_msg(patchID && flat_patchID && (flat_patch.n_faces() > 0),
                     "lift_min_surface: Error in the flat patch.");
  int id = flat_patchID[Face(0)];

  // Label the newly added vertices for debugging
  auto is_inflated =
      inflated_patch.vertex_property<bool>("v:is_inflated", false);
  is_inflated.vector().assign(is_inflated.vector().size(), false);
  auto is_inflated_edge =
      inflated_patch.edge_property<bool>("e:is_inflated_edge", false);
  is_inflated_edge.vector().assign(is_inflated_edge.vector().size(), false);

  if (!to_inflat_labels.empty() &&
      !flat_patch.get_edge_property<BoundaryType>("e:boundary_type") &&
      !flat_patch.get_face_property<int>("f:patch_component")) {
    logger().warn("lift_min_surface: Patch is not decomposed.");
  }

  auto boundary_type =
      inflated_patch.edge_property<BoundaryType>("e:boundary_type");

  auto patch_component =
      inflated_patch.get_face_property<int>("f:patch_component");
  auto inflated_orig_f_idx =
      inflated_patch.vertex_property<int>("v:orig_f_idx", -1);
  auto inflated_orig_idx =
      inflated_patch.vertex_property<int>("v:orig_idx", -1);

  // Iterate on flat edges
  for (int i = 0; i < (int)flat_patch.n_edges(); i++) {
    Edge flat_e(i);

    // Check if we want to inflat it
    if (!to_inflat_labels.empty() &&
        !to_inflat_labels.count(boundary_type[flat_e]))
      continue;

    Vertex src_flat = flat_patch.vertex(flat_e, 0);
    Vertex dst_flat = flat_patch.vertex(flat_e, 1);
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
    std::vector<int> to_relax({0, 0});
    std::vector<std::pair<Vector3f, int>> samples;
    size_t lift_trial = 0;
    std::vector<int> relax_params({2, 5, 10});

    // Don't lift if one neighboring component is degenerated
    auto flat_comp_valid = flat_patch.get_face_property<bool>("f:comp_valid");
    bool is_valid = true;
    if ((flat_patch.face(flat_e, 0).is_valid() &&
         !flat_comp_valid[flat_patch.face(flat_e, 0)]) ||
        (flat_patch.face(flat_e, 1).is_valid() &&
         !flat_comp_valid[flat_patch.face(flat_e, 1)])) {
      is_valid = false;
      logger().info("Skip invalid patch {}, comp {}, {}", id,
                    (flat_patch.face(flat_e, 0).is_valid())
                        ? patch_component[flat_patch.face(flat_e, 0)]
                        : -1,
                    (flat_patch.face(flat_e, 1).is_valid())
                        ? patch_component[flat_patch.face(flat_e, 1)]
                        : -1);
    }

    if (is_valid) {
      do {
        Vector3f cut_p = camera.position();
        walk_on_original(mesh, id, src_orig, flat_vpositions[src_flat],
                         dst_orig, flat_vpositions[dst_flat], cut_p, path_norm,
                         camera, sampling_delta, min_delta_2d, samples,
                         to_relax);
        break;

        if (samples.empty()) {
          to_relax = std::vector<int>(
              {relax_params[lift_trial], relax_params[lift_trial]});
          logger().info(
              "Trial lift component boundary: #{}, {}-ring, {} -> {}; "
              "Origin: {} -> {}",
              lift_trial, relax_params[lift_trial], src_flat, dst_flat,
              Vertex(orig_idx[src_flat]), Vertex(orig_idx[dst_flat]));
        }
        lift_trial++;
      } while (lift_trial < relax_params.size() && samples.empty());
    }

    if (samples.empty()) {
      failed_paths.emplace_back(
          std::make_pair(std::min(orig_idx[src_flat], orig_idx[dst_flat]),
                         std::max(orig_idx[src_flat], orig_idx[dst_flat])));
      continue;
    }

    // Insert to the result mesh
    Edge to_insert_e = inflated_patch.find_edge(src_flat, dst_flat);
    BoundaryType flat_e_type = boundary_type[to_insert_e];
    Halfedge next_he = inflated_patch.find_halfedge(src_flat, dst_flat);
    int l_comp = -1, r_comp = -1;

    for (size_t i = 0; i < samples.size(); i++) {
      Vector3f pos = samples[i].first;
      Vertex prev_v = (inflated_patch.vertex(to_insert_e, 0) != dst_flat)
                          ? inflated_patch.vertex(to_insert_e, 0)
                          : inflated_patch.vertex(to_insert_e, 1);

      // Update component labeling
      if (patch_component && i == 0) {
        l_comp = patch_component[inflated_patch.face(next_he)];
        r_comp = patch_component[inflated_patch.face(
            inflated_patch.opposite_halfedge(next_he))];
      }

      Vertex new_v = inflated_patch.split(to_insert_e, pos);
      is_inflated[new_v] = true;
      inflated_orig_f_idx[new_v] = samples[i].second;
      inflated_orig_idx[new_v] = -1;

      auto hit = inflated_patch.halfedges(new_v), hit_end = hit;
      Halfedge prev_he;
      next_he = Halfedge(-1);
      do {
        if (inflated_patch.to_vertex(*hit) == dst_flat) {
          next_he = *hit;
          // break;
        } else if (inflated_patch.to_vertex(*hit) == prev_v) {
          prev_he = *hit;
        }

        // Reset all new edges
        boundary_type[inflated_patch.edge(*hit)] = BoundaryType::NONE;

      } while (++hit != hit_end);

      contess_assert_msg(next_he.is_valid() && prev_he.is_valid(),
                         "lift_min_surface: Error inserting vertex samples.");
      to_insert_e = inflated_patch.edge(next_he);

      // Debug label
      Edge to_label_e = inflated_patch.edge(prev_he);
      is_inflated_edge[to_label_e] = true;
      // Mantain the component boundary label
      boundary_type[to_label_e] = flat_e_type;

      // Fill the last edge
      if (i + 1 >= samples.size())
        boundary_type[inflated_patch.edge(next_he)] = flat_e_type;

      // Update component labeling
      if (patch_component) {
        contess_assert_msg(
            l_comp >= 0 && r_comp >= 0 && l_comp != r_comp,
            "lift_min_surface: Error look up component indices.");
        if (i + 1 >= samples.size()) {
          patch_component[inflated_patch.face(next_he)] = l_comp;
          patch_component[inflated_patch.face(
              inflated_patch.opposite_halfedge(next_he))] = r_comp;
        }
        patch_component[inflated_patch.face(prev_he)] = r_comp;
        patch_component[inflated_patch.face(
            inflated_patch.opposite_halfedge(prev_he))] = l_comp;
      }
    }
  }

  auto inflated_patchID = inflated_patch.face_property<int>("f:patchID");
  inflated_patchID.vector().assign(inflated_patchID.vector().size(), id);
}

void inflat_path(Mesh const &mesh, Mesh const &flat_patch, Camera const &camera,
                 real_t sampling_delta, real_t min_delta_2d,
                 std::vector<std::pair<Vector3f, int>> &inflat_samples,
                 std::unordered_set<int> const &to_inflat_comp_labels) {
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
  for (int i = 0; i < (int)flat_patch.n_edges(); i++) {
    Edge flat_e(i);

    // Check if we want to inflat it
    // No boundary points
    if (!(boundary_type[flat_e] == BoundaryType::NONE))
      continue;

    // Check if it's from the desired component
    if (!to_inflat_comp_labels.empty() &&
        ((!flat_patch.face(flat_e, 0).is_valid() ||
          !to_inflat_comp_labels.count(
              patch_component[flat_patch.face(flat_e, 0)])) ||
         (!flat_patch.face(flat_e, 1).is_valid() ||
          !to_inflat_comp_labels.count(
              patch_component[flat_patch.face(flat_e, 1)]))))
      continue;

    Vertex src_flat = flat_patch.vertex(flat_e, 0);
    Vertex dst_flat = flat_patch.vertex(flat_e, 1);
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

    inflat_samples.insert(inflat_samples.end(), samples.begin(), samples.end());
  }
}
