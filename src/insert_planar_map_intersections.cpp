// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "insert_planar_map_intersections.h"

#include "aabb.h"
#include "common.h"
#include "fix_flipped_faces.h"
#include "insert_interpolated_contours.h"
#include "subdivide_contour_edges.h"
#include "sweepLine.h"
#include <algorithm>
#include <limits>
#include <spdlog/fmt/ostr.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

void edges_intersection(
    Mesh const &mesh, Camera const &camera, std::vector<Edge> const &edges,
    std::map<int, std::vector<std::pair<real_t, int>>> &intersections) {
  intersections.clear();
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");

  struct Edge2d {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    int edge_idx;
    Vector2f v1, v2;
    AABB aabb;

    std::vector<std::pair<real_t, int>> intersection_ts;
  };

  std::vector<Edge2d> edge2d;
  for (auto eit = edges.begin(); eit != edges.end(); ++eit) {
    // Project all to 2D
    Edge2d e2d;
    e2d.edge_idx = (*eit).idx();
    e2d.v1 = project(vpositions[mesh.vertex(*eit, 0)],
                     camera.viewMatrix().matrix(), camera.projectionMatrix(),
                     Vector2i(camera.vpWidth(), camera.vpHeight()))
                 .head(2);
    e2d.v2 = project(vpositions[mesh.vertex(*eit, 1)],
                     camera.viewMatrix().matrix(), camera.projectionMatrix(),
                     Vector2i(camera.vpWidth(), camera.vpHeight()))
                 .head(2);

    // AABB
    e2d.aabb.min = Vector3f(std::min(e2d.v1.x(), e2d.v2.x()),
                            std::min(e2d.v1.y(), e2d.v2.y()), 0);
    e2d.aabb.max = Vector3f(std::max(e2d.v1.x(), e2d.v2.x()),
                            std::max(e2d.v1.y(), e2d.v2.y()), 1);
    edge2d.emplace_back(e2d);
  }

  // Intersection
  for (size_t i = 0; i + 1 < edge2d.size(); i++) {
    Edge2d &ei = edge2d[i];
    for (size_t j = i + 1; j < edge2d.size(); j++) {
      Edge2d &ej = edge2d[j];

      // Skip the intersection check if edges are adjacent
      if (mesh.vertex(Edge(ei.edge_idx), 0) ==
              mesh.vertex(Edge(ej.edge_idx), 0) ||
          mesh.vertex(Edge(ei.edge_idx), 0) ==
              mesh.vertex(Edge(ej.edge_idx), 1) ||
          mesh.vertex(Edge(ei.edge_idx), 1) ==
              mesh.vertex(Edge(ej.edge_idx), 0) ||
          mesh.vertex(Edge(ei.edge_idx), 1) ==
              mesh.vertex(Edge(ej.edge_idx), 1))
        continue;

      // First check the AABB
      if (!ei.aabb.aabbIntersect(ej.aabb))
        continue;

      // Real intersection check
      real_t t, u;
      if (intersect2dSeg2dSegParametric(ei.v1.data(), ei.v2.data(),
                                        ej.v1.data(), ej.v2.data(), t, u,
                                        EPSILON) == 1) {
        // Insert if the intersection is not at the end point
        bool t_in = false, u_in = false;
        if (!(std::abs(t - 1) < EPSILON || std::abs(t) < EPSILON) && t > 0 &&
            t < 1)
          t_in = true;
        if (!(std::abs(u - 1) < EPSILON || std::abs(u) < EPSILON) && u > 0 &&
            u < 1)
          u_in = true;

        if (t_in && u_in) {
          ei.intersection_ts.emplace_back(t, ej.edge_idx);
          ej.intersection_ts.emplace_back(u, ei.edge_idx);
        }
      }
    }
  }

  // Read out (not really necessary)
  for (size_t i = 0; i < edge2d.size(); i++) {
    Edge2d &ei = edge2d[i];

    if (ei.intersection_ts.empty())
      continue;

    // Sort the intersections on edge
    std::sort(ei.intersection_ts.begin(), ei.intersection_ts.end());
    intersections[ei.edge_idx] = ei.intersection_ts;
  }
}

void insert_planar_map_intersections(Mesh &mesh, Camera const &camera,
                                     bool to_match_endpoints) {
  auto ndotv = mesh.vertex_property<real_t>("v:ndotv");
  auto concave = mesh.get_edge_property<bool>("e:concave");
  auto is_contour = mesh.edge_property<real_t>("e:contour");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  if (patchID)
    patchID = mesh.face_property<int>("f:patchID");

  auto intersection_2d = mesh.get_vertex_property<int>("v:intersection_2d");
  if (!intersection_2d) {
    intersection_2d = mesh.add_vertex_property<int>("v:intersection_2d", -1);
  }
  intersection_2d.vector().assign(intersection_2d.vector().size(), -1);

  auto intersection_edges =
      mesh.get_vertex_property<Vector2i>("v:intersection_edges");
  if (!intersection_edges) {
    intersection_edges = mesh.add_vertex_property<Vector2i>(
        "v:intersection_edges", -1 * Vector2i::Ones());
  }
  intersection_edges.vector().assign(intersection_edges.vector().size(),
                                     -1 * Vector2i::Ones());

  contess_assert_msg(is_contour,
                     "insert_planar_map_intersections: Contours are missing. "
                     "Planar map intersections can't be created.");

  auto param_loc = mesh.halfedge_property<Param_loc>("h:param_loc");

  std::vector<Edge> contour_edges;
  for (auto eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit) {
    if (is_contour[*eit] >= 0 || mesh.is_boundary(*eit)) {
      contour_edges.emplace_back(*eit);
    }
  }

  std::map<int, std::vector<std::pair<real_t, int>>> intersections;
  edges_intersection(mesh, camera, contour_edges, intersections);

  std::vector<std::pair<Vertex, Vector2f>> intersection_points;
  Vector2i viewport = Vector2i(camera.vpWidth(), camera.vpHeight());

  bool skip_assert = false;

  // Insert
  auto find_matching_param = [&](Vertex const &from, Vertex const &to,
                                 double param_t, Param_loc &param) {
    std::unordered_map<int, std::vector<Param_loc>> param_matching;
    auto hit = mesh.halfedges(from), hit_end = hit;
    do {
      if (!param_loc[*hit].is_valid() ||
          param_matching.count(param_loc[*hit].ptexIndex))
        continue;
      param_matching[param_loc[*hit].ptexIndex] = std::vector<Param_loc>();
      param_matching[param_loc[*hit].ptexIndex].emplace_back(param_loc[*hit]);
    } while (++hit != hit_end);
    hit = mesh.halfedges(to), hit_end = hit;
    do {
      if (!param_loc[*hit].is_valid())
        continue;
      if (!param_matching.count(param_loc[*hit].ptexIndex))
        param_matching[param_loc[*hit].ptexIndex] = std::vector<Param_loc>();
      param_matching[param_loc[*hit].ptexIndex].emplace_back(param_loc[*hit]);
    } while (++hit != hit_end);

    for (auto const &params : param_matching) {
      if (params.second.size() < 2)
        continue;
      param.ptexIndex = params.first;
      param.uv =
          (1 - param_t) * params.second[0].uv + param_t * params.second[1].uv;
      return;
    }

    skip_assert = true;
  };
  for (auto ei : intersections) {
    if (ei.second.empty())
      continue;

    // Getting all the intersection point positions
    int orig_e_idx = ei.first;
    Edge orig_e = Edge(ei.first);
    Vertex init_v = mesh.vertex(orig_e, 0);
    Vertex final_v = mesh.vertex(orig_e, 1);
    Edge to_split = orig_e;

    bool originally_contour = is_contour[to_split] >= 0;
    for (auto int_record : ei.second) {
      real_t t = int_record.first;
      surface_mesh::Point p =
          (1 - t) * vpositions[init_v] + t * vpositions[final_v];
      std::unordered_set<int> split_vertices;
      split_vertices.insert(mesh.vertex(to_split, 0).idx());
      split_vertices.insert(mesh.vertex(to_split, 1).idx());

      // For parameter updates
      Vertex v1 = mesh.vertex(to_split, 0), v2 = mesh.vertex(to_split, 1);
      real_t param_t = (p - vpositions[v1]).norm() /
                       (vpositions[v2] - vpositions[v1]).norm();
      Param_loc v1_param = param_loc[mesh.find_halfedge(v1, v2)];
      Param_loc v2_param =
          param_loc[mesh.next_halfedge(mesh.find_halfedge(v1, v2))];
      Param_loc v_param(v1_param.ptexIndex,
                        (1 - param_t) * v1_param.uv + param_t * v2_param.uv);
      Param_loc v_param_valid;
      find_matching_param(v1, v2, param_t, v_param_valid);

      std::unordered_map<int, Param_loc> updated_params;
      updated_params[v_param.ptexIndex] = v_param;
      {
        Param_loc v1_param = param_loc[mesh.find_halfedge(v2, v1)];
        Param_loc v2_param =
            param_loc[mesh.next_halfedge(mesh.find_halfedge(v2, v1))];
        Param_loc v_param(v2_param.ptexIndex,
                          (1 - param_t) * v2_param.uv + param_t * v1_param.uv);
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

      int is_concave = -1;
      if (concave) {
        is_concave = concave[to_split];
      }
      Vertex new_v = mesh.split(to_split, p);

      // Intersection matching
      intersection_edges[new_v] = Vector2i(orig_e_idx, int_record.second);

      Vector2f v_2d = project(p, camera.viewMatrix().matrix(),
                              camera.projectionMatrix(), viewport)
                          .head<2>();
      intersection_points.emplace_back(new_v, v_2d);

      // Find the new edge to split
      bool found_split_edge = false;
      auto hit = mesh.halfedges(new_v);
      auto hit_end = hit;
      Param_loc param_res;
      do {
        if (originally_contour &&
            split_vertices.find(mesh.to_vertex(*hit).idx()) !=
                split_vertices.end()) {
          is_contour[mesh.edge(*hit)] = 1;
          if (is_concave >= 0)
            concave[mesh.edge(*hit)] = is_concave;
        }
        if (mesh.to_vertex(*hit) == final_v) {
          to_split = mesh.edge(*hit);
          found_split_edge = true;
        }

        // Update the parameter
        int n_ptexIndex = param_loc[mesh.next_halfedge(*hit)].ptexIndex;
        param_loc[*hit] = updated_params[n_ptexIndex];
        param_res = param_loc[*hit];

        param_loc[mesh.opposite_halfedge(*hit)] =
            updated_params_v[mesh.from_vertex(mesh.opposite_halfedge(*hit))
                                 .idx()];
      } while (++hit != hit_end);

      if (!param_res.is_valid() && v_param_valid.is_valid())
        param_res = v_param_valid;

      contess_assert_msg(found_split_edge, "Missing next edge to split.");

      // Update facing information
      if (!skip_assert) {
        contess_assert_msg(param_res.is_valid(),
                           "Parameter error at the new vertex.");
        Vector3f pos_res, normal_res;
        mesh.subdivision().evaluateLimit(param_res, pos_res, normal_res);
        real_t v_ndotv =
            normal_res.dot((camera.position() - pos_res).normalized());
        ndotv[new_v] = v_ndotv;
      }
    }
  }

  // Fill patch ID for new faces
  if (patchID) {
    std::vector<Face> fill_patch_faces;
    for (size_t i = 0; i < mesh.n_faces(); i++) {
      Face f(i);

      if (patchID[f] == -1)
        fill_patch_faces.emplace_back(f);
    }
    update_patch(mesh, camera, fill_patch_faces);
  }

  // Match the intersections
  std::vector<bool> intersection_points_matched;
  intersection_points_matched.resize(intersection_points.size(), false);
  for (size_t i = 0; i + 1 < intersection_points.size(); i++) {
    if (intersection_points_matched[i])
      continue;

    intersection_points_matched[i] = true;

    // Match as the closest 2D point
    size_t matched_i = i;
    real_t min_dist = std::numeric_limits<real_t>::infinity();
    for (size_t j = i + 1; j < intersection_points.size(); j++) {
      if (intersection_points_matched[j])
        continue;

      // Avoid matching two intersection points close in 3D
      if (intersection_edges[intersection_points[j].first][0] !=
          intersection_edges[intersection_points[i].first][1])
        continue;

      real_t dist =
          (intersection_points[i].second - intersection_points[j].second)
              .norm();
      if (dist < min_dist) {
        min_dist = dist;
        matched_i = j;
        if (min_dist < std::numeric_limits<real_t>::epsilon())
          break;
      }
    }

    if (matched_i == i) {
      continue;
    }

    intersection_2d[intersection_points[i].first] =
        intersection_points[matched_i].first.idx();
    intersection_2d[intersection_points[matched_i].first] =
        intersection_points[i].first.idx();
    intersection_points_matched[matched_i] = true;
  }

  // Match existing vertices
  if (to_match_endpoints) {
    real_t MATCH_EPSILON = 1e-2;
    intersection_points_matched.resize(mesh.n_vertices(), false);
    for (size_t i = 0; i + 1 < mesh.n_vertices(); ++i) {
      auto is_adjacent_to_contour = [&](Vertex const &v) -> bool {
        Edge c_e;
        auto hit = mesh.halfedges(v), hend = hit;
        do {
          if (!(*hit).is_valid())
            return false;
          Edge hit_e = mesh.edge(*hit);
          if (!hit_e.is_valid())
            continue;
          if (is_contour[hit_e] >= 0.f || mesh.is_boundary(hit_e)) {
            c_e = hit_e;
            break;
          }
        } while (++hit != hend);
        return c_e.is_valid();
      };
      Vertex v_i(i);
      if (!is_adjacent_to_contour(v_i) || intersection_2d[v_i] >= 0)
        continue;

      // Match as the closest 2D point
      size_t matched_i = i;
      real_t min_dist = std::numeric_limits<real_t>::infinity();
      for (size_t j = i + 1; j < mesh.n_vertices(); j++) {
        Vertex v_j(j);
        if (intersection_points_matched[j] || !is_adjacent_to_contour(v_j))
          continue;

        Vector2f vi_2d = project(vpositions[v_i], camera.viewMatrix().matrix(),
                                 camera.projectionMatrix(), viewport)
                             .head<2>();
        Vector2f vj_2d = project(vpositions[v_j], camera.viewMatrix().matrix(),
                                 camera.projectionMatrix(), viewport)
                             .head<2>();

        real_t dist = (vi_2d - vj_2d).norm();
        if (dist < min_dist &&
            (vpositions[v_i] - vpositions[v_j]).norm() > MATCH_EPSILON) {
          min_dist = dist;
          matched_i = j;
          if (min_dist < std::numeric_limits<real_t>::epsilon())
            break;
        }
      }
      if (matched_i == i) {
        continue;
      }

      if (min_dist > MATCH_EPSILON)
        continue;

      Vertex v_j(matched_i);
      intersection_2d[v_i] = matched_i;
      intersection_2d[v_j] = i;
      intersection_points_matched[matched_i] = true;
    }
  }

  // Since we have new faces, we need to update VBO and VBO_f
  if (patchID) {
    consistently_label_interpolated_contours(mesh, camera.position(), false);
    mesh.markPatchBoundaryEdges(false);
  }
}