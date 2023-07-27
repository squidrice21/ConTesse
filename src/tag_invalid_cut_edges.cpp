// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "tag_invalid_cut_edges.h"
#include "common.h"
#include <spdlog/fmt/ostr.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "make_cut_feasible.h"

bool edges_intersect(Mesh const &mesh, Camera const &camera,
                     Vector3f const &pp2, Vector3f const &pp3,
                     std::vector<Vector3f> const &tri_p) {
  // Check if corners are on different sides / on the plane.
  bool has_seen_pos = false, has_seen_neg = false;
  for (uint8_t i = 0; i < tri_p.size(); i++) {
    Vector3f corner = tri_p[i];
    auto side = igl::predicates::orient3d(camera.position(), pp2, pp3, corner);
    switch (side) {
    case igl::predicates::Orientation::NEGATIVE:
      has_seen_neg = true;
      break;
    case igl::predicates::Orientation::POSITIVE:
      has_seen_pos = true;
      break;
    default: // Orientation::COPLANAR
      has_seen_neg = true;
      has_seen_pos = true;
      break;
    }

    if (has_seen_pos && has_seen_neg) {
      break;
    }
  }

  return has_seen_pos && has_seen_neg;
}

bool is_edge_valid_cut(Mesh const &mesh, Camera const &camera, Edge const &e,
                       Edge const &contour_e, size_t contour_neighbor_size_1d,
                       bool &successful) {
  successful = true;
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto VBO_f = mesh.get_face_property<FacingType>("f:VBO_f");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  auto is_cusp = mesh.get_vertex_property<bool>("v:cusp");

  // 0. Avoid all edges belonging to any inconsistent triangle
  // Since we would have a flipped cut if it goes through two edges of the face.
  // It may be fine just going through a single one.
  auto check_face_consistency = [&](Face const &f) -> bool {
    return !f.is_valid() || (VBO[f] == VBO_f[f]);
  };
  bool is_face_consistent = check_face_consistency(mesh.face(e, 0)) &&
                            check_face_consistency(mesh.face(e, 1));
  if (!is_face_consistent)
    return false;

  // Assume the edge is in the patch interior
  int e_pid = patchID[mesh.face(e, 0)];

  // 1. Grow the contour edge
  auto grow_contour = [&](Vertex const &v, std::unordered_set<int> &visited_v,
                          std::vector<Edge> &contour_edges) {
    Vertex cur_v = v;
    size_t itr = 0;

    do {
      auto hit = mesh.halfedges(cur_v);
      auto hit_end = hit;
      do {
        if (is_contour[mesh.edge(*hit)] < 0 &&
            !mesh.is_boundary(mesh.edge(*hit)))
          continue;

        // Make sure we are walking along the right patch
        bool seen_right_patch = false;
        if ((*hit).is_valid() && mesh.face(*hit).is_valid() &&
            patchID[mesh.face(*hit)] == e_pid)
          seen_right_patch = true;
        auto oppo_he = mesh.opposite_halfedge(*hit);
        if ((oppo_he).is_valid() && mesh.face(oppo_he).is_valid() &&
            patchID[mesh.face(oppo_he)] == e_pid)
          seen_right_patch = true;
        if (!seen_right_patch)
          continue;

        if (visited_v.find(mesh.to_vertex(*hit).idx()) == visited_v.end()) {
          contour_edges.emplace_back(mesh.edge(*hit));
          visited_v.emplace(cur_v.idx());
          cur_v = mesh.to_vertex(*hit);

          break;
        }

      } while (++hit != hit_end);
      itr++;
      if (is_cusp[cur_v])
        break;
    } while (itr < contour_neighbor_size_1d);
  };
  std::unordered_set<int> visited_v;
  visited_v.emplace(mesh.vertex(contour_e, 0).idx());
  visited_v.emplace(mesh.vertex(contour_e, 1).idx());
  std::vector<Edge> contour_edges;
  contour_edges.emplace_back(contour_e);
  grow_contour(mesh.vertex(contour_e, 0), visited_v, contour_edges);
  grow_contour(mesh.vertex(contour_e, 1), visited_v, contour_edges);

  std::unordered_set<int> contour_edge_lookup;
  for (auto const &c_e : contour_edges) {
    contour_edge_lookup.emplace(c_e.idx());
  }

  for (auto const &c_e : contour_edges) {
    // 2.1 Find if the two edges intersect at an end vertex
    {
      std::unordered_map<int, size_t> vertex_count;
      std::vector<Edge> edges({c_e, e});
      for (auto const &ee : edges) {
        if (!vertex_count.count(mesh.vertex(ee, 0).idx()))
          vertex_count[mesh.vertex(ee, 0).idx()] = 0;
        vertex_count[mesh.vertex(ee, 0).idx()]++;
        if (!vertex_count.count(mesh.vertex(ee, 1).idx()))
          vertex_count[mesh.vertex(ee, 1).idx()] = 0;
        vertex_count[mesh.vertex(ee, 1).idx()]++;
      }
      Vertex c_v;
      for (auto const &count : vertex_count) {
        if (count.second == 2)
          c_v = Vertex(count.first);
      }
      // Intersect at a vertex
      if (c_v.is_valid()) {
        Vertex oppo_e_v =
            (c_v == mesh.vertex(e, 0)) ? mesh.vertex(e, 1) : mesh.vertex(e, 0);

        std::vector<Halfedge> c_h;
        std::vector<bool> orient;
        auto hit = mesh.halfedges(c_v), hit_end = hit;
        do {
          Edge e = mesh.edge(*hit);
          if (is_contour[e] >= 0 || mesh.is_boundary(e)) {
            if (mesh.face(*hit).is_valid() && patchID[mesh.face(*hit)] == e_pid)
              c_h.emplace_back(*hit);
            else if (mesh.face(mesh.opposite_halfedge(*hit)).is_valid() &&
                     patchID[mesh.face(mesh.opposite_halfedge(*hit))] == e_pid)
              c_h.emplace_back(mesh.opposite_halfedge(*hit));
          }
        } while (++hit != hit_end);

        if (c_h.size() != 2) {
          successful = false;
          return false;
        }

        contess_assert_msg(c_h.size() == 2,
                           "is_edge_valid_cut: Contour is broken at v" +
                               std::to_string(c_v.idx()));

        bool on_correct_side = is_on_correct_side(mesh, camera, c_h[0], c_h[1],
                                                  vpositions[oppo_e_v]);
        if (!on_correct_side)
          return false;
        continue; // So we don't go to 2.1 for the current contour edge
      }
    }

    // 2.1 Otherwise, check if they intersect in the edge interior
    {
      std::vector<Vector3f> e_p(
          {vpositions[mesh.vertex(c_e, 0)], vpositions[mesh.vertex(c_e, 1)]});
      std::vector<Vector3f> e_p2(
          {vpositions[mesh.vertex(e, 0)], vpositions[mesh.vertex(e, 1)]});
      // 2.1.1 If the edge intersects with any contour edge
      if (edges_intersect(mesh, camera, vpositions[mesh.vertex(e, 0)],
                          vpositions[mesh.vertex(e, 1)], e_p) &&
          edges_intersect(mesh, camera, vpositions[mesh.vertex(c_e, 0)],
                          vpositions[mesh.vertex(c_e, 1)], e_p2)) {
        return false;
      }
    }
  }

  return true;
}

bool tag_invalid_cut_edges(Mesh &mesh, Camera const &camera,
                           size_t contour_neighbor_size,
                           size_t contour_neighbor_size_1d,
                           bool use_simple_rule) {
  auto VBO_f = mesh.get_face_property<FacingType>("f:VBO_f");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");
  auto is_cusp = mesh.get_vertex_property<bool>("v:cusp");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto near_contour = mesh.get_face_property<int>("f:near_contour");
  if (!near_contour) {
    near_contour = mesh.add_face_property<int>("f:near_contour", -1);
  }
  near_contour.vector().assign(near_contour.vector().size(), -1);

  auto cut_candidate = mesh.get_edge_property<bool>("e:cut_candidate");
  if (!cut_candidate) {
    cut_candidate = mesh.add_edge_property<bool>("e:cut_candidate", true);
  }
  cut_candidate.vector().assign(cut_candidate.vector().size(), true);

  std::unordered_map<int, std::vector<Edge>> close_to;
  std::unordered_map<int, std::unordered_set<int>> recorded_e;

  auto check_face_consistency = [&](Face const &f) -> bool {
    return !f.is_valid() || (VBO[f] == VBO_f[f]);
  };

  // 1. Flood from contour loops
  for (Edge const &e : mesh.edges()) {
    if (!(is_contour[e] >= 0.f || mesh.is_boundary(e))) {
      // For non contour, check if it belongs to an inconsistent face
      bool is_v = check_face_consistency(mesh.face(e, 0)) &&
                  check_face_consistency(mesh.face(e, 1));
      if (!is_v)
        cut_candidate[e] = false;
      continue;
    }

    std::unordered_map<int, size_t> face_dist;
    std::queue<Face> queue;
    if (mesh.face(e, 0).is_valid()) {
      queue.push(mesh.face(e, 0));
      face_dist[mesh.face(e, 0).idx()] = 0;
    }
    if (mesh.face(e, 1).is_valid()) {
      queue.push(mesh.face(e, 1));
      face_dist[mesh.face(e, 1).idx()] = 0;
    }
    while (!queue.empty()) {
      Face current_f = queue.front();
      queue.pop();
      near_contour[current_f] = e.idx();

      if (face_dist[current_f.idx()] >= contour_neighbor_size)
        continue;

      // add adjacent faces if no cut boundary or contour is crossed
      // And book keep the corresponding source contour edge
      auto hit = mesh.halfedges(current_f), hend = hit;
      do {
        Edge hit_e = mesh.edge(*hit);
        if (is_contour[hit_e] >= 0.f || mesh.is_boundary(hit_e))
          continue;

        if (!close_to.count(hit_e.idx())) {
          close_to[hit_e.idx()] = std::vector<Edge>();
          recorded_e[hit_e.idx()] = std::unordered_set<int>();
        }
        if (!recorded_e[hit_e.idx()].count(e.idx())) {
          close_to[hit_e.idx()].emplace_back(e);
          recorded_e[hit_e.idx()].emplace(e.idx());
        }

        // adjacent_f is in the same patch
        Face adjacent_f = mesh.face(mesh.opposite_halfedge(*hit));
        if (adjacent_f.is_valid() && !face_dist.count(adjacent_f.idx())) {
          face_dist[adjacent_f.idx()] = face_dist[current_f.idx()] + 1;
          queue.push(adjacent_f);
        }
      } while (++hit != hend);
    }
  }

  // Add edges directly adjacent to the contour
  for (Edge const &e : mesh.edges()) {
    if (close_to.count(e.idx()) || is_contour[e] >= 0.f ||
        mesh.is_boundary(e) || !cut_candidate[e])
      continue;
    auto is_adjacent_to_contour = [&](Vertex const &v) -> Edge {
      Edge c_e;
      auto hit = mesh.halfedges(v), hend = hit;
      do {
        Edge hit_e = mesh.edge(*hit);
        if (is_contour[hit_e] >= 0.f || mesh.is_boundary(hit_e)) {
          c_e = hit_e;
          break;
        }
      } while (++hit != hend);
      return c_e;
    };
    Edge c_e = is_adjacent_to_contour(mesh.vertex(e, 0));
    if (!c_e.is_valid())
      c_e = is_adjacent_to_contour(mesh.vertex(e, 1));
    if (c_e.is_valid())
      close_to[e.idx()] = std::vector<Edge>({c_e});
  }

  // 2. Test if collected edges are valid. Each edge is tested against all
  // contour/boundary edges within the given range
  for (auto const &close_record : close_to) {
    Edge e(close_record.first);
    // Already determined
    if (!cut_candidate[e])
      continue;

    bool is_valid = true;
    for (auto const &c_e : close_record.second) {
      bool is_v = true;
      // If edge is adjacent to a cusp
      if (is_cusp[mesh.vertex(e, 0)] || is_cusp[mesh.vertex(e, 1)]) {
        is_v = false;
      } else {
        if (!use_simple_rule) {
          bool successful;
          is_v = is_edge_valid_cut(mesh, camera, e, c_e,
                                   contour_neighbor_size_1d, successful);
          if (!successful)
            return false;
        }
        // Simple rule only checks if the edge is adjacent to any inconsistent
        // face
        else {
          is_v = check_face_consistency(mesh.face(e, 0)) &&
                 check_face_consistency(mesh.face(e, 1));
        }
      }

      if (!is_v) {
        is_valid = false;
        break;
      }
    }

    cut_candidate[e] = is_valid;
  }

  return true;
}
