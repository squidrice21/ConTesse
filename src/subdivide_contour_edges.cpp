// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "subdivide_contour_edges.h"
#include "subdivide_contour_edges_even.h"

#include "common.h"
#include <algorithm>
#include <cmath>
#include <deque>
#include <igl/predicates/predicates.h>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// Eigen printout formatting
#include <spdlog/fmt/ostr.h>

#include "insert_interpolated_contours.h"
#include "subdiv_common.h"
#include "surface_mesh/surface_mesh.h"

void update_patch(Mesh &mesh, Camera const &camera,
                  std::vector<Face> const &unassigned_faces) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto VBO = mesh.face_property<FacingType>("f:VBO");
  auto VBO_f = mesh.face_property<FacingType>("f:VBO_f");
  auto patchID = mesh.face_property<int>("f:patchID");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");

  // Note this unassigned faces fetching is based on the default patchID value
  // of -1 (assigned when adding) In this way, we make sure all patches would be
  // filled with the original patchID regardless of its size
  std::deque<Face> face_queue;
  if (unassigned_faces.empty()) {
    for (size_t i = 0; i < mesh.n_faces(); i++) {
      Face f(i);
      if (patchID[f] == -1)
        face_queue.push_back(f);
    }
  } else {
    for (auto const &f : unassigned_faces)
      face_queue.push_back(f);
  }

  std::pair<int, size_t> inf_loop_check;
  inf_loop_check.first = -1;
  while (!face_queue.empty()) {
    Face f = face_queue.front();
    face_queue.pop_front();

    bool updated = false;
    auto hit = mesh.halfedges(f);
    auto hit_end = hit;
    do {
      Face adj_f = mesh.face(mesh.opposite_halfedge(*hit));
      if (adj_f.is_valid() && is_contour[mesh.edge(*hit)] < 0 &&
          patchID[adj_f] >= 0) {
        patchID[f] = patchID[adj_f];
        VBO[f] = mesh.get_patch_facing(patchID[f]);

        Vertex vs[3];
        mesh.verticesOfFace(f, vs);
        VBO_f[f] =
            (igl::predicates::orient3d(vpositions[vs[0]], vpositions[vs[1]],
                                       vpositions[vs[2]], camera.position()) ==
             igl::predicates::Orientation::NEGATIVE)
                ? FacingType::FRONT
                : FacingType::BACK;
        updated = true;
        break;
      }
    } while (++hit != hit_end);

    if (!updated && face_queue.empty())
      break;

    if (!updated) {
      face_queue.push_back(f);

      if (inf_loop_check.first < 0) {
        inf_loop_check.first = f.idx();
        inf_loop_check.second = face_queue.size();
      } else if (inf_loop_check.first == f.idx()) {
        contess_assert_msg(
            inf_loop_check.second != face_queue.size(),
            "update_patch: Detect face unable to be re-assigned " +
                std::to_string(f.idx()));
      }
    } else {
      if (inf_loop_check.first == f.idx()) {
        inf_loop_check.first = -1;
      }
    }
  }

  // Update patch bounary
  if (unassigned_faces.empty()) {
    auto patchBoundary = mesh.halfedge_property<int>("h:patchBoundary");
    parallel_for(
        tbb::blocked_range<uint32_t>(0u, (uint32_t)mesh.n_edges(), GRAIN_SIZE),
        [&](const tbb::blocked_range<uint32_t> &range) {
          for (uint32_t i = range.begin(); i != range.end(); ++i) {
            Edge e(i);
            Halfedge h0 = mesh.halfedge(e, 0);
            Halfedge h1 = mesh.halfedge(e, 1);
            // test if this edge is on an open boundary
            if (mesh.is_boundary(h0) && mesh.face(h1).is_valid()) {
              patchBoundary[h1] = patchID[mesh.face(h1)];
            } else if (mesh.is_boundary(h1) && mesh.face(h0).is_valid()) {
              patchBoundary[h0] = patchID[mesh.face(h0)];
            } else if (mesh.face(h1).is_valid() && mesh.face(h0).is_valid()) {
              // test if this edge is on a patch boundary
              int id0 = patchID[mesh.face(h0)];
              int id1 = patchID[mesh.face(h1)];
              if (id0 != id1) {
                patchBoundary[h0] = id0;
                patchBoundary[h1] = id1;
              }
            }
          }
        });
  }
}

Vertex get_far_vertex(Mesh &mesh, Face const &f, Edge const &imm_e) {
  Vertex vv[3];
  mesh.verticesOfFace(f, vv);

  for (auto const &v : vv) {
    if (mesh.vertex(imm_e, 0) != v && mesh.vertex(imm_e, 1) != v)
      return v;
  }
  contess_assert_msg(0, "Cannot find a far vertex.");
  return Vertex();
}

bool is_triangle_fine(Mesh const &mesh, std::vector<Vertex> vv,
                      real_t min_angle_radius, real_t &best_angle) {
  best_angle = 1e3;
  bool is_fine = true;
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  contess_assert(vv.size() == 3);
  for (size_t i = 0; i < vv.size(); i++) {
    size_t j = (i + 1) % vv.size();
    size_t k = (i + 2) % vv.size();

    auto vec1 = (vpositions[vv[j]] - vpositions[vv[i]]).normalized();
    auto vec2 = (vpositions[vv[j]] - vpositions[vv[k]]).normalized();
    real_t dot_prod = vec1.dot(vec2);
    real_t radius = std::acos(dot_prod);
    best_angle = std::min(best_angle, radius);
    if (radius < min_angle_radius)
      is_fine = false;
  }
  return is_fine;
}

void improve_adjacent_topology(Mesh &mesh, Edge const &orig_e,
                               real_t min_angle_radius,
                               std::vector<bool> &refined) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto param_loc = mesh.halfedge_property<Param_loc>("h:param_loc");
  std::unordered_set<int> seen_faces;
  std::vector<Face> ref_ff;
  if (mesh.face(orig_e, 0).is_valid())
    ref_ff.emplace_back(mesh.face(orig_e, 0));
  if (mesh.face(orig_e, 1).is_valid())
    ref_ff.emplace_back(mesh.face(orig_e, 1));
  contess_assert(ref_ff.size() >= 1);

  real_t best_angle = 1e3;
  for (auto const &f : ref_ff) {
    std::deque<Face> ff;
    ff.push_back(f);

    // Determine if the triangle has the right orientation
    Vertex oppo_v = get_far_vertex(mesh, f, orig_e);
    Vector3f norm = mesh.compute_face_normal(f);
    Vector3f norm_v = vpositions[mesh.vertex(orig_e, 0)] + norm;
    auto orig_orien = igl::predicates::orient3d(
        vpositions[mesh.vertex(orig_e, 0)], vpositions[mesh.vertex(orig_e, 1)],
        norm_v, vpositions[oppo_v]);

    std::unordered_map<int, int> prev_face;
    prev_face[f.idx()] = -1;
    std::unordered_map<int, int> imm_ee;
    imm_ee[f.idx()] = orig_e.idx();
    Face final_f;
    while (!ff.empty()) {
      std::vector<Vertex> vv({mesh.vertex(orig_e, 0), mesh.vertex(orig_e, 1)});
      Face ref_f = ff.front();
      ff.pop_front();

      Edge imm_e(imm_ee[ref_f.idx()]);
      Vertex v = get_far_vertex(mesh, ref_f, imm_e);

      // Check orientation
      auto ref_orien = igl::predicates::orient3d(
          vpositions[mesh.vertex(orig_e, 0)],
          vpositions[mesh.vertex(orig_e, 1)], norm_v, vpositions[v]);
      if (ref_orien != orig_orien)
        continue;

      vv.emplace_back(v);
      real_t f_best_angle;
      bool to_end = is_triangle_fine(mesh, vv, min_angle_radius, f_best_angle);
      if (to_end) {
        final_f = ref_f;
        break;
      }
      // When we can't achieve the goal, keep the best result
      if (f_best_angle > best_angle) {
        final_f = ref_f;
      }

      // Extend
      auto hit = mesh.halfedges(ref_f), hit_end = hit;
      do {
        Edge e = mesh.edge(*hit);
        // Cann't cross parameterization quad...
        auto oppo_he = mesh.opposite_halfedge(*hit);
        Face adj_f = mesh.face(oppo_he);
        if (is_contour[e] >= 0 || e == imm_e || !adj_f.is_valid() ||
            prev_face.count(adj_f.idx()) ||
            param_loc[oppo_he].ptexIndex != param_loc[*hit].ptexIndex)
          continue;

        imm_ee[adj_f.idx()] = e.idx();
        prev_face[adj_f.idx()] = ref_f.idx();
        ff.push_back(adj_f);
      } while (++hit != hit_end);
    }

    if (!final_f.is_valid()) {
      refined.emplace_back(false);
      continue;
    }

    auto update_param = [&](Halfedge const &he, int ptexIndex) {
      Param_loc existing_param;
      auto hit = mesh.halfedges(mesh.from_vertex(he));
      auto hit_end = hit;
      do {
        if (*hit == he)
          continue;
        if (param_loc[*hit].ptexIndex == ptexIndex) {
          existing_param = param_loc[*hit];
          break;
        }
      } while (++hit != hit_end);
      param_loc[he] = existing_param;
    };

    // Flip faces
    Face flip_f = final_f;
    while (imm_ee[flip_f.idx()] != orig_e.idx()) {
      Edge flip_e(imm_ee[flip_f.idx()]);
      Param_loc before_flip = param_loc[mesh.halfedge(flip_e, 0)];
      Vertex v1 = mesh.to_vertex(mesh.next_halfedge(mesh.halfedge(flip_e, 0)));
      Vertex v2 = mesh.to_vertex(mesh.next_halfedge(mesh.halfedge(flip_e, 1)));

      bool is_okay = mesh.is_flip_ok(flip_e);
      if (!is_okay) {
        refined.emplace_back(false);
        break;
      }
      mesh.flip(flip_e);
      // Update parameter
      update_param(mesh.find_halfedge(v1, v2), before_flip.ptexIndex);
      update_param(mesh.find_halfedge(v2, v1), before_flip.ptexIndex);

      flip_f = Face(prev_face[flip_f.idx()]);
    }
    if (imm_ee[flip_f.idx()] == orig_e.idx())
      refined.emplace_back(true);

    if (!refined.back()) {
      logger().warn("improve_adjacent_topology: Unable to improve the topology "
                    "of face {}.",
                    f);
    }
  }
}

void subdivide_contour_edges(Mesh &mesh, Camera const &camera,
                             std::vector<Edge> &subdiv_edges) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto vpositions = mesh.vertex_property<Vector3f>("v:point");
  auto param_loc = mesh.halfedge_property<Param_loc>("h:param_loc");
  auto facing = mesh.vertex_property<FacingType>("v:facing");
  auto ndotv_v = mesh.vertex_property<real_t>("v:ndotv");
  auto vnormals = mesh.get_vertex_property<Vector3f>("v:normal");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");

  auto patchID = mesh.get_face_property<int>("f:patchID");
  auto cut = mesh.get_edge_property<int>("e:disk_cut");
  surface_mesh::Surface_mesh::Edge_property<int> cut_edit;
  if (cut) {
    cut_edit = mesh.edge_property<int>("e:disk_cut");
  }

  auto is_boundary_edge = mesh.get_edge_property<bool>("e:is_boundary");
  if (is_boundary_edge) {
    is_boundary_edge = mesh.edge_property<bool>("e:is_boundary");
  }

  auto evaluate_g_value = [&](Vertex const &v, bool to_update) -> FacingType {
    Halfedge he = mesh.halfedge(v);
    Param_loc param_res;
    param_res.ptexIndex = param_loc[he].ptexIndex;
    param_res.uv = param_loc[he].uv;
    if (!param_res.is_valid()) {
      auto hit = mesh.halfedges(v);
      auto hit_end = hit;
      do {
        if (!param_loc[*hit].is_valid()) {
          continue;
        }
        param_res = param_loc[*hit];
        break;
      } while (++hit != hit_end);
    }
    Vector3f pos_res, normal_res;

    mesh.subdivision().evaluateLimit(param_res, pos_res, normal_res);
    real_t ndotv = normal_res.dot((camera.position() - pos_res).normalized());

    // Update cache
    if (to_update) {
      vnormals[v] = normal_res;
      ndotv_v[v] = ndotv;
    } else {
      ndotv = ndotv_v[v];
    }

    if (ndotv < -CONTOUR_THRESHOLD)
      return BACK;

    if (ndotv > CONTOUR_THRESHOLD)
      return FRONT;

    return CONTOUR;
  };

  std::map<std::pair<int, int>, int> existing_cut_edges;
  if (cut) {
    // Since this function only affects contour edges which are for sure not cut
    // edges, we don't need to save the values of e:disk_cut
    for (size_t i = 0; i < mesh.n_edges(); i++) {
      Edge e(i);
      if (cut[e] >= 0) {
        auto e_p = std::make_pair(
            std::min(mesh.vertex(e, 0).idx(), mesh.vertex(e, 1).idx()),
            std::max(mesh.vertex(e, 0).idx(), mesh.vertex(e, 1).idx()));
        existing_cut_edges.insert(std::make_pair(e_p, cut[Edge(i)]));
      }
    }
  }

  std::vector<Mesh::EdgeTuple> edges_to_split;
  for (auto const &e : subdiv_edges) {
    contess_assert_msg(is_contour[e] >= 0 || mesh.is_boundary(e),
                       "Cannot subdivide non contour edges.");

    // 1. Unset the contour label
    is_contour[e] = -1;

    // 1.5 Improve the local topology to allow for larger moving range
    real_t min_angle_radius = 20. / 180 * M_PI;
    std::vector<bool> refined;
    if (!mesh.is_boundary(e))
      improve_adjacent_topology(mesh, e, min_angle_radius, refined);

    // 2. Split edge in the middle, determine the side containing the
    // zero-crossing
    Vertex v1 = mesh.vertex(e, 0), v2 = mesh.vertex(e, 1);

    Vector3f mid = 0.5 * (vpositions[v1] + vpositions[v2]);
    Param_loc v1_param = param_loc[mesh.find_halfedge(v1, v2)];
    Param_loc v2_param =
        param_loc[mesh.next_halfedge(mesh.find_halfedge(v1, v2))];
    Param_loc v_param(v1_param.ptexIndex, 0.5 * (v1_param.uv + v2_param.uv));

    std::unordered_map<int, Param_loc> updated_params;
    updated_params[v_param.ptexIndex] = v_param;
    {
      Param_loc v1_param = param_loc[mesh.find_halfedge(v2, v1)];
      Param_loc v2_param =
          param_loc[mesh.next_halfedge(mesh.find_halfedge(v2, v1))];
      Param_loc v_param(v2_param.ptexIndex, 0.5 * (v1_param.uv + v2_param.uv));
      updated_params[v_param.ptexIndex] = v_param;
    }

    std::unordered_map<int, Param_loc> updated_params_v;
    {
      std::vector<Halfedge> hes(
          {mesh.find_halfedge(v1, v2), mesh.find_halfedge(v2, v1)});
      for (auto const &he : hes) {
        updated_params_v[mesh.from_vertex(he).idx()] = param_loc[he];
        Halfedge n_he = mesh.next_halfedge(mesh.next_halfedge(he));
        updated_params_v[mesh.from_vertex(n_he).idx()] = param_loc[n_he];
      }
    }

    bool is_e_boundary = mesh.is_boundary(e);
    Vertex v = mesh.split(e, mid);

    Edge invert_e;
    auto hit = mesh.halfedges(v);
    auto hit_end = hit;
    do {
      // Update the parameter
      int n_ptexIndex = param_loc[mesh.next_halfedge(*hit)].ptexIndex;
      param_loc[*hit] = updated_params[n_ptexIndex];

      param_loc[mesh.opposite_halfedge(*hit)] =
          updated_params_v[mesh.from_vertex(mesh.opposite_halfedge(*hit))
                               .idx()];

      if (mesh.to_vertex(*hit) == v1 || mesh.to_vertex(*hit) == v2) {
        continue;
      }

      FacingType v_facing = evaluate_g_value(mesh.from_vertex(*hit), true);
      FacingType v_n_facing = evaluate_g_value(mesh.to_vertex(*hit), false);
      facing[v] = v_facing;

      if (v_facing != v_n_facing && v_facing != FacingType::CONTOUR &&
          v_n_facing != FacingType::CONTOUR) {
        invert_e = mesh.edge(*hit);
      }
    } while (++hit != hit_end);

    // Handle boundary
    if (is_e_boundary) {
      Param_loc param_res;
      auto hit = mesh.halfedges(v);
      auto hit_end = hit;
      do {
        if (!param_loc[*hit].is_valid()) {
          continue;
        }
        param_res = param_loc[*hit];
        break;
      } while (++hit != hit_end);

      Vector3f pos_res, normal_res;
      mesh.subdivision().evaluateLimit(param_res, pos_res, normal_res);
      vpositions[v] = pos_res;
      continue;
    }

    if (!invert_e.is_valid()) {
      if (abs(ndotv_v[v]) > CONTOUR_THRESHOLD)
        logger().warn("Cannot find zero crossing at {} - {} - {}", v1.idx(),
                      v.idx(), v2.idx());
      contess_assert_msg(abs(ndotv_v[v]) <= CONTOUR_THRESHOLD,
                         "No zero-crossing found in the new edges...");

      continue;
    }

    // 3. Insert or move (in the same way as the initial insertion)
    Vector3f p, n;
    Param_loc param;
    real_t c = mesh.bisect_search(camera.position(), &mesh.const_subdivision(),
                                  invert_e, p, n, param);
    edges_to_split.push_back(Mesh::EdgeTuple(invert_e, p, n, param, c));
  }

  // Actually split FB edges
  bool shift_allowed = true;
  for (auto etuple : edges_to_split) {
    mesh.splitEdge(etuple, shift_allowed);
  }

  mesh.garbage_collection();

  // Update the cut flag if the cut paths exist
  if (cut || is_boundary_edge) {
    // New edges are added to the back of an indexed list
    for (size_t i = 0; i < mesh.n_edges(); i++) {
      Edge e(i);

      if (cut) {
        cut_edit[e] = -1;

        auto e_p = std::make_pair(
            std::min(mesh.vertex(e, 0).idx(), mesh.vertex(e, 1).idx()),
            std::max(mesh.vertex(e, 0).idx(), mesh.vertex(e, 1).idx()));
        if (existing_cut_edges.count(e_p))
          cut_edit[e] = existing_cut_edges[e_p];
      }

      // Update boundary debug flag
      if (is_boundary_edge) {
        is_boundary_edge[e] = mesh.is_boundary(e);
      }
    }
  }

  // Update the contour and facing
  for (size_t i = 0; i < mesh.n_edges(); i++) {
    Edge e(i);
    is_contour[e] = -1;
  }
  // Rebuild the face labeling and contour
  consistently_label_interpolated_contours(mesh, camera.position());
  for (size_t i = 0; i < mesh.n_edges(); i++) {
    Edge e(i);
    if (!mesh.face(e, 0).is_valid() || !mesh.face(e, 1).is_valid())
      continue;
    FacingType f1 = VBO[mesh.face(e, 0)];
    FacingType f2 = VBO[mesh.face(e, 1)];
    if (f1 != f2)
      is_contour[e] = 1;
  }

  // Update the patch
  if (patchID)
    update_patch(mesh, camera);
}

void subdivide_contour_edges(Mesh &mesh, Camera const &camera) {
  std::vector<Edge> subdiv_edges;
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  contess_assert(is_contour);

  for (size_t i = 0; i < mesh.n_edges(); i++) {
    Edge e(i);
    if (is_contour[e] < 0)
      continue;

    subdiv_edges.emplace_back(e);
  }

  subdivide_contour_edges(mesh, camera, subdiv_edges);
}

void subdivide_contour_edges_k_times(Mesh &mesh, Camera const &camera,
                                     std::vector<Edge> &subdiv_edges,
                                     size_t subdiv_times, bool is_even,
                                     real_t min_2d_length) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  contess_assert(is_contour);
  Vector2i viewport = Vector2i(camera.vpWidth(), camera.vpHeight());

  auto is_long_enough = [&](Edge const &e) -> bool {
    if (min_2d_length < 0)
      return true;
    Vector2f v1_2d =
        project(vpositions[mesh.vertex(e, 0)], camera.viewMatrix().matrix(),
                camera.projectionMatrix(), viewport)
            .head<2>();
    Vector2f v2_2d =
        project(vpositions[mesh.vertex(e, 1)], camera.viewMatrix().matrix(),
                camera.projectionMatrix(), viewport)
            .head<2>();
    real_t len_2d = (v1_2d - v2_2d).norm();

    if (len_2d <= min_2d_length)
      logger().info("Drop edge: {}", e);
    return len_2d > min_2d_length;
  };

  std::vector<std::pair<Vertex, Vertex>> subdiv_vertices;
  for (size_t i = 0; i < subdiv_times; i++) {
    subdiv_vertices.clear();
    for (auto const &e : subdiv_edges) {
      logger().info("Subdiv {}: {}, {}", i, mesh.vertex(e, 0),
                    mesh.vertex(e, 1));
      subdiv_vertices.emplace_back(mesh.vertex(e, 0), mesh.vertex(e, 1));
    }

    if (!is_even)
      subdivide_contour_edges(mesh, camera, subdiv_edges);
    else
      subdivide_contour_edges_even(mesh, camera, subdiv_edges);

    subdiv_edges.clear();
    auto is_contour_adjacent = [&](Vertex const &from,
                                   Vertex const &to) -> bool {
      auto hit = mesh.halfedges(from);
      auto hit_end = hit;
      do {
        Edge e = mesh.edge(*hit);
        if ((is_contour[e] >= 0 || mesh.is_boundary(e)) &&
            mesh.to_vertex(*hit) == to) {
          return true;
        }
      } while (++hit != hit_end);
      return false;
    };
    for (auto const &v_pair : subdiv_vertices) {
      auto v1 = v_pair.first;
      auto v2 = v_pair.second;

      bool seen_subdivided_edges = false;
      auto hit = mesh.halfedges(v1);
      auto hit_end = hit;
      do {
        Edge e = mesh.edge(*hit);
        if (!is_long_enough(e))
          continue;
        if (is_contour[e] >= 0 || mesh.is_boundary(e)) {
          if (mesh.to_vertex(*hit) == v2 ||
              is_contour_adjacent(mesh.to_vertex(*hit), v2)) {
            subdiv_edges.emplace_back(e);
            seen_subdivided_edges = true;
            break;
          }
        }
      } while (++hit != hit_end);
      hit = mesh.halfedges(v2);
      hit_end = hit;
      do {
        Edge e = mesh.edge(*hit);
        if (!is_long_enough(e))
          continue;
        if (is_contour[e] >= 0 || mesh.is_boundary(e)) {
          if (is_contour_adjacent(mesh.to_vertex(*hit), v1)) {
            subdiv_edges.emplace_back(e);
            seen_subdivided_edges = true;
            break;
          }
        }
      } while (++hit != hit_end);

      contess_assert(min_2d_length >= 0 || seen_subdivided_edges);
      if (min_2d_length >= 0 && !seen_subdivided_edges)
        return;
    }

    if (i + 1 >= subdiv_times) {
      subdiv_vertices.clear();
      for (auto const &e : subdiv_edges) {
        logger().info("Subdiv {}: {}, {}", i, mesh.vertex(e, 0),
                      mesh.vertex(e, 1));
        subdiv_vertices.emplace_back(mesh.vertex(e, 0), mesh.vertex(e, 1));
      }

      auto vpositions = mesh.vertex_property<Vector3f>("v:point");
      auto patchID = mesh.get_face_property<int>("f:patchID");
      contess_assert(patchID);
      std::unordered_map<int, Vector3f> inflated;
      std::unordered_set<int> to_inflat;
      for (auto const &e : subdiv_edges) {
        if (!mesh.is_boundary(e))
          continue;
        subdiv_vertices.emplace_back(mesh.vertex(e, 0), mesh.vertex(e, 1));
        to_inflat.emplace(mesh.vertex(e, 0).idx());
        to_inflat.emplace(mesh.vertex(e, 1).idx());
      }
      for (auto const v_i : to_inflat) {
        Vertex v(v_i);
        std::vector<Vertex> vv;
        auto hit = mesh.halfedges(v);
        auto hit_end = hit;
        Edge e;
        bool seen_contour = false;
        do {
          if (!mesh.is_boundary(mesh.edge(*hit)))
            continue;

          if (is_contour[mesh.edge(*hit)] >= 0)
            seen_contour = true;

          e = mesh.edge(*hit);
          vv.emplace_back(mesh.to_vertex(*hit));
        } while (++hit != hit_end);

        // We don't want to move the boundary-contour joint point
        if (seen_contour)
          continue;

        Vertex v1 = vv[0];
        Vertex v2 = vv[1];
        real_t w1 = 1 / (vpositions[v1] - vpositions[v]).norm();
        real_t w2 = 1 / (vpositions[v2] - vpositions[v]).norm();
        Vector3f laplacian =
            ((w1 * vpositions[v1] + w2 * vpositions[v2]) / (w1 + w2) -
             vpositions[v]);
        real_t laplacian_norm = laplacian.norm();
        laplacian /= laplacian_norm;

        Face adj_f = mesh.face(e, 0);
        if (!adj_f.is_valid())
          adj_f = mesh.face(e, 1);
        Vector3f f_norm = mesh.compute_face_normal(adj_f);
        if (laplacian_norm < 1e-3)
          laplacian = f_norm;
        if (f_norm.dot(laplacian) > 0)
          laplacian *= -1;

        FacingType facing = mesh.get_patch_facing(patchID[adj_f]);
        if (facing == FacingType::BACK)
          laplacian *= -1;

        real_t laplacian_offset = 1e-3;
        inflated[v.idx()] = vpositions[v] - laplacian_offset * laplacian;
        logger().info("Inflat: {} - {}", v.idx(), laplacian.transpose());
      }
      for (auto const &inflated_v : inflated) {
        vpositions[Vertex(inflated_v.first)] = inflated_v.second;
      }
    }
  }
}

void subdivide_contour_edges_tagged(Mesh &mesh, Camera const &camera,
                                    size_t subdiv_times, bool is_even,
                                    real_t min_2d_length) {
  auto subdiv_contour = mesh.get_edge_property<bool>("e:subdiv_contour");
  if (!subdiv_contour) {
    logger().warn(
        "subdivide_contour_edges_tagged: No tagged contour edges. Return.");
  }
  std::vector<Edge> contour_edges;
  for (uint32_t i = 0; i != mesh.n_edges(); ++i) {
    Edge e(i);
    if (!subdiv_contour[e])
      continue;
    contour_edges.emplace_back(e);
  }
  subdivide_contour_edges_k_times(mesh, camera, contour_edges, subdiv_times,
                                  is_even, min_2d_length);
}
