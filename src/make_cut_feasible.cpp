// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "make_cut_feasible.h"
#include "chain_contour.h"
#include "common.h"
#include "cut_patch_to_disk.h"
#include "subdivide_contour_edges.h"
#include <deque>
#include <igl/predicates/predicates.h>
#include <limits>
#include <spdlog/fmt/ostr.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

bool is_on_correct_side(Mesh const &mesh, Camera const &camera,
                        Halfedge const &c_h1, Halfedge const &c_h2,
                        Vector3f const &p) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto VBO_f = mesh.get_face_property<FacingType>("f:VBO_f");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");

  contess_assert_msg(mesh.from_vertex(c_h1) == mesh.to_vertex(c_h2) ||
                         mesh.from_vertex(c_h2) == mesh.to_vertex(c_h1),
                     "is_on_right_side: Mismatched contour edges.");

  std::vector<Halfedge> c_h;
  c_h.emplace_back(c_h1);
  c_h.emplace_back(c_h2);
  std::vector<bool> orient;
  for (size_t i = 0; i < c_h.size(); i++) {
    // 1. Determine the right side wrt c_h1
    auto h = c_h1;

    auto ori = igl::predicates::orient3d(
        vpositions[mesh.from_vertex(h)], vpositions[mesh.to_vertex(h)],
        vpositions[mesh.to_vertex(mesh.next_halfedge(h))], camera.position());
    // If we have two vertices at the same location
    if (!(ori == igl::predicates::Orientation::NEGATIVE ||
          ori == igl::predicates::Orientation::POSITIVE) &&
        (vpositions[mesh.from_vertex(h)] - vpositions[mesh.to_vertex(h)])
                .norm() < std::numeric_limits<real_t>::epsilon()) {
      logger().warn("is_on_right_side: Contour degenerates to a point: {}, {}",
                    mesh.from_vertex(h), mesh.to_vertex(h));
      return true;
    }
    contess_assert_msg(ori == igl::predicates::Orientation::NEGATIVE ||
                           ori == igl::predicates::Orientation::POSITIVE,
                       "is_on_right_side: Generic camera assumption violated.");
    auto side = (ori == igl::predicates::Orientation::NEGATIVE);
    contess_assert_msg(mesh.face(h).is_valid(),
                       "is_on_right_side: Wrong manifold boundary halfedges.");
    if (VBO[mesh.face(h)] != VBO_f[mesh.face(h)])
      side = !side;
    orient.emplace_back(side);
  }

  // 2. Determine the side of c_h2 wrt c_h1
  bool is_union = true;
  Vertex c_v2 = mesh.to_vertex(c_h2);
  auto ori = igl::predicates::orient3d(vpositions[mesh.from_vertex(c_h1)],
                                       vpositions[mesh.to_vertex(c_h1)],
                                       vpositions[c_v2], camera.position());
  if (ori == igl::predicates::Orientation::NEGATIVE ||
      ori == igl::predicates::Orientation::POSITIVE) {
    bool side2 = (ori == igl::predicates::Orientation::NEGATIVE);
    // Boundary is concave, use the intersection of the condition
    if (side2 == orient[0])
      is_union = false;
  }

  auto is_on_right_side = [&](Halfedge const &h, bool const &is_neg,
                              Vector3f split_candidate) -> bool {
    auto ori = igl::predicates::orient3d(vpositions[mesh.from_vertex(h)],
                                         vpositions[mesh.to_vertex(h)],
                                         split_candidate, camera.position());
    if (ori != igl::predicates::Orientation::NEGATIVE &&
        ori != igl::predicates::Orientation::POSITIVE)
      return false;
    auto side = (ori == igl::predicates::Orientation::NEGATIVE);

    return side == is_neg;
  };
  if (!is_union) {
    if (is_on_right_side(c_h[0], orient[0], p) &&
        is_on_right_side(c_h[1], orient[1], p)) {
      return true;
    }
  } else {
    if (is_on_right_side(c_h[0], orient[0], p) ||
        is_on_right_side(c_h[1], orient[1], p)) {
      return true;
    }
  }
  return false;
}

bool find_split_right_side(Mesh const &mesh, Camera const &camera,
                           Face const &f, Halfedge const &h, int patch_id,
                           std::vector<Vector3f> &split_positions) {
  auto patchID = mesh.get_face_property<int>("f:patchID");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");

  Vertex contour_v = mesh.from_vertex(h);

  // 1. Find the two adjacent boundary edges
  std::vector<Halfedge> c_h;
  std::vector<bool> orient;
  auto hit = mesh.halfedges(contour_v), hit_end = hit;
  do {
    Edge e = mesh.edge(*hit);
    if (is_contour[e] >= 0 || mesh.is_boundary(e)) {
      if (mesh.face(*hit).is_valid() && patchID[mesh.face(*hit)] == patch_id)
        c_h.emplace_back(*hit);
      else if (mesh.face(mesh.opposite_halfedge(*hit)).is_valid() &&
               patchID[mesh.face(mesh.opposite_halfedge(*hit))] == patch_id)
        c_h.emplace_back(mesh.opposite_halfedge(*hit));
    }
  } while (++hit != hit_end);

  if (c_h.size() < 2)
    return false;

  // 2. Grid search on the edge
  Halfedge next_h = mesh.next_halfedge(h);
  double steps = 10;
  for (size_t i = 1; i < steps; i++) {
    Vector3f split_candidate =
        (1 - (i / steps)) * vpositions[mesh.from_vertex(next_h)] +
        (i / steps) * vpositions[mesh.to_vertex(next_h)];
    if (is_on_correct_side(mesh, camera, c_h[0], c_h[1], split_candidate))
      split_positions.emplace_back(split_candidate);
  }

  return !split_positions.empty();
}

Vertex split_edge(Mesh &mesh, Camera const &camera, Halfedge const &h,
                  Vector3f const &split_p) {
  auto facing = mesh.get_vertex_property<FacingType>("v:facing");
  auto vnormals = mesh.get_vertex_property<Vector3f>("v:normal");
  auto ndotv = mesh.get_vertex_property<real_t>("v:ndotv");
  auto VBO_f = mesh.get_face_property<FacingType>("f:VBO_f");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");

  auto vpositions = mesh.vertex_property<Vector3f>("v:point");
  auto param_loc = mesh.halfedge_property<Param_loc>("h:param_loc");

  if (facing)
    facing = mesh.vertex_property<FacingType>("v:facing");
  if (vnormals)
    vnormals = mesh.vertex_property<Vector3f>("v:normal");
  if (ndotv)
    ndotv = mesh.vertex_property<real_t>("v:ndotv");
  if (VBO_f)
    VBO_f = mesh.face_property<FacingType>("f:VBO_f");
  if (VBO)
    VBO = mesh.face_property<FacingType>("f:VBO");

  FacingType neighbor_vbo = FacingType::NA;
  if (VBO)
    neighbor_vbo = VBO[mesh.face(h)];

  Edge e = mesh.edge(h);
  Param_loc new_loc[2];

  real_t c =
      (split_p - vpositions[mesh.from_vertex(h)]).norm() /
      (vpositions[mesh.to_vertex(h)] - vpositions[mesh.from_vertex(h)]).norm();
  new_loc[0].ptexIndex = param_loc[h].ptexIndex;
  new_loc[0].uv =
      (1.f - c) * param_loc[h].uv + c * param_loc[mesh.next_halfedge(h)].uv;
  Halfedge h1 = mesh.opposite_halfedge(h);

  // Deal with edges between two ptex faces
  if (!mesh.is_boundary(e) &&
      param_loc[h].ptexIndex != param_loc[h1].ptexIndex) {
    new_loc[1].ptexIndex = param_loc[h1].ptexIndex;
    new_loc[1].uv =
        c * param_loc[h1].uv + (1.f - c) * param_loc[mesh.next_halfedge(h1)].uv;
  } else {
    new_loc[1] = new_loc[0];
  }
  // save parametric locations of the 4 halfedges (3 for boundaries)
  std::map<Vertex, Param_loc> vertices_param;
  for (int i = 0; i < 2; i++) {
    if (!mesh.is_boundary(mesh.halfedge(e, i))) {
      Halfedge hh = mesh.halfedge(e, i);
      vertices_param[mesh.from_vertex(hh)] = param_loc[hh];
      Halfedge nh = mesh.next_halfedge(mesh.next_halfedge(hh));
      vertices_param[mesh.from_vertex(nh)] = param_loc[nh];
    }
  }

  Vertex v = mesh.split(e, split_p);

  if (facing)
    facing[v] = facing[mesh.from_vertex(h)];
  if (vnormals)
    vnormals[v] = (1.f - c) * vnormals[mesh.from_vertex(h)] +
                  c * vnormals[mesh.to_vertex(h)];
  if (ndotv)
    ndotv[v] = (camera.position() - split_p).normalized().dot(vnormals[v]);

  // update parametric locations
  auto chit = mesh.halfedges(v);
  auto chit_end = chit;
  do {
    // Refill face facing
    if (mesh.face(*chit).is_valid() && VBO_f && VBO) {
      VBO[mesh.face(*chit)] = VBO_f[mesh.face(*chit)] = neighbor_vbo;
    }

    if (!mesh.is_boundary(*chit)) {
      if (param_loc[mesh.next_halfedge(*chit)].ptexIndex ==
          new_loc[0].ptexIndex) {
        param_loc[*chit] = new_loc[0];
      } else {
        param_loc[*chit] = new_loc[1];
      }
    }
    Halfedge oh = mesh.opposite_halfedge(*chit);
    if (!mesh.is_boundary(oh)) {
      param_loc[oh] = vertices_param[mesh.from_vertex(oh)];
    }
    ++chit;
  } while (chit != chit_end);

  return v;
}

void make_cut_feasible(Mesh &mesh, Camera const &camera, bool ff_only) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  contess_assert_msg(patchID, "Missing patches.");
  auto VBO_f = mesh.get_face_property<FacingType>("f:VBO_f");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");
  auto cut_candidate = mesh.get_edge_property<bool>("e:cut_candidate");
  contess_assert_msg(cut_candidate,
                     "make_cut_feasible: Needs to label candidate edges.");

  // 0. Build chains if not existing
  if (mesh.get_oriented_chains().empty())
    chain_contour(mesh, camera, true, ff_only);

  // 1. Find patch loops that is adjacent to only inconsistent edges
  std::vector<std::pair<Chain, int>> infeasible_loops;
  std::vector<bool> infeasible_loops_inserted;
  std::unordered_set<int> patch_ids;
  for (auto const &patch : mesh.get_patch_chains()) {
    patch_ids.emplace(patch.first);
  }

  auto is_face_non_candidate_free = [&](Face const &f) -> bool {
    auto hit = mesh.halfedges(f), hit_end = hit;
    do {
      if (!cut_candidate[mesh.edge(*hit)])
        return false;
    } while (++hit != hit_end);
    return true;
  };

  for (auto const &pid : patch_ids) {
    if (ff_only && mesh.get_patch_facing(pid) != FacingType::FRONT)
      continue;
    logger().info("Patch: {}", pid);
    // Check loops
    auto chains = mesh.get_patch_chains().equal_range(pid);
    int patch_id = pid;
    for (auto chain = chains.first; chain != chains.second; ++chain) {
      bool seen_consistent_edges = false;
      Halfedge h1 = (*chain->second).front();

      for (auto const &b_cut_he : *chain->second) {
        Vertex v = mesh.from_vertex(b_cut_he);
        auto hit = mesh.halfedges(v), hit_end = hit;
        do {
          Edge e = mesh.edge(*hit);

          // Skip if it's boundary or it's in a different patch
          if (!mesh.face(e, 0).is_valid() || !mesh.face(e, 1).is_valid())
            continue;
          if (patchID[mesh.face(e, 0)] != patchID[mesh.face(e, 1)])
            continue;
          if (patchID[mesh.face(e, 0)] != patch_id)
            continue;

          if (cut_candidate[e] && is_face_non_candidate_free(mesh.face(e, 0)) &&
              is_face_non_candidate_free(mesh.face(e, 1))) {
            seen_consistent_edges = true;
            break;
          }
        } while (++hit != hit_end);
      }

      if (!seen_consistent_edges) {
        logger().info("=> No consistent edge outward: {}, {}", pid,
                      mesh.from_vertex(h1));
        Chain c = *chain->second;
        infeasible_loops.emplace_back(c, patch_id);
      }
    }
  }

  // 2. Split adjacent consistent triangles to create adjacent consistent edges.
  // When splitting make sure the split edge is on the right side.
  infeasible_loops_inserted.resize(infeasible_loops.size(), false);
  struct SplitRecord {
    Vertex v1, v2;
    std::vector<Vector3f> split_positions;
    // Loop info
    int infeasible_index;
  };
  std::vector<SplitRecord> split_record;
  for (size_t j = 0; j < infeasible_loops.size(); j++) {
    auto const &infeasible_case = infeasible_loops[j];
    auto loop = infeasible_case.first;
    int patch_id = infeasible_case.second;

    // Gather adjacent consistent faces
    std::unordered_set<int> seen_faces;
    // std::vector<std::pair<Face, Halfedge>> consistent_faces;
    std::unordered_map<int, std::deque<Halfedge>> consistent_faces;

    auto find_candidate_halfedge_in_face = [&](Face const &f,
                                               std::deque<Halfedge> &queue) {
      auto hit = mesh.halfedges(f), hit_end = hit;
      do {
        Edge e = mesh.edge(*hit);
        if (is_contour[e] >= 0 || mesh.is_boundary(e))
          continue;
        // This may free some existing candidate edges
        if (!cut_candidate[e])
          queue.push_front(*hit);
        else
          queue.push_back(*hit);
      } while (++hit != hit_end);
    };

    for (auto const &b_cut_he : loop) {
      Vertex v = mesh.from_vertex(b_cut_he);

      auto hit = mesh.halfedges(v), hit_end = hit;
      do {
        Face f = mesh.face(*hit);

        // Skip if it's boundary or it's in a different patch
        if (!f.is_valid())
          continue;
        if (patchID[f] != patch_id)
          continue;

        if (VBO[f] == VBO_f[f] && seen_faces.count(f.idx()) == 0) {
          seen_faces.emplace(f.idx());
          find_candidate_halfedge_in_face(f, consistent_faces[f.idx()]);
        }
      } while (++hit != hit_end);
    }

    if (consistent_faces.empty())
      logger().error(
          "make_cut_feasible: Cannot add feasible cut edge to loop: v" +
          std::to_string(mesh.to_vertex(loop.front()).idx()) + "; p" +
          std::to_string(patch_id));

    // Split all found faces
    for (auto &f_h : consistent_faces) {
      bool found = false;
      while (!f_h.second.empty() && !found) {
        Halfedge h = f_h.second.front();
        Halfedge next_h = mesh.next_halfedge(h);
        f_h.second.pop_front();

        // Avoid inserting to manifold boundary and contour
        if (mesh.is_boundary(mesh.edge(next_h)) ||
            is_contour[mesh.edge(next_h)] >= 0)
          continue;

        std::vector<Vector3f> split_positions;
        found = find_split_right_side(mesh, camera, Face(f_h.first), h,
                                      patch_id, split_positions);
        if (found) {
          SplitRecord record;
          record.v1 = mesh.from_vertex(next_h);
          record.v2 = mesh.to_vertex(next_h);
          record.split_positions = split_positions;
          record.infeasible_index = j;
          split_record.emplace_back(record);
        }
      }
    }
  }

  {
    // Actual insertion
    for (auto const &record : split_record) {
      Halfedge next_h = mesh.find_halfedge(record.v1, record.v2);
      if (!next_h.is_valid()) {
        logger().warn("make_cut_feasible: {}, {} not connected.", record.v1,
                      record.v2);
        break;
      }
      std::vector<Vector3f> split_positions = record.split_positions;

      size_t max_insertion_number = 3;
      for (size_t j = 0; j < max_insertion_number && j < split_positions.size();
           j++) {
        Vertex end_v = mesh.to_vertex(next_h);
        Vertex new_start_v =
            split_edge(mesh, camera, next_h, split_positions[j]);
        if (new_start_v.is_valid())
          next_h = mesh.find_halfedge(new_start_v, end_v);
        else
          next_h = Halfedge();
        if (!next_h.is_valid()) {
          logger().warn("make_cut_feasible: {}, {} not connected.", new_start_v,
                        end_v);
          break;
        }
        infeasible_loops_inserted[record.infeasible_index] = true;
      }
    }

    for (size_t i = 0; i < infeasible_loops_inserted.size(); i++) {
      if (infeasible_loops_inserted[i])
        continue;

      Halfedge h1 = infeasible_loops[i].first.front();
      int patch_id = infeasible_loops[i].second;
      logger().warn(
          "make_cut_feasible: Cannot insert consistency to Patch {}; Start: {}",
          patch_id, mesh.from_vertex(h1));
    }
  }

  // Update patch ID
  update_patch(mesh, camera);
}
