// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "insert_subdivision_cusps.h"
#include "common.h"
#include "evaluate_radial_curvature.h"
#include "make_cut_feasible.h"
#include "sweepLine.h"
#include <algorithm>
#include <limits>
#include <queue>
#include <spdlog/fmt/ostr.h>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

bool has_zero_crossing(Mesh &mesh, Subdiv const &subdiv, Camera const &camera,
                       std::vector<Param_loc> const &triangle,
                       real_t *tri_area = nullptr) {
  std::vector<real_t> g_values;
  std::vector<real_t> k_r_values;
  std::vector<Vector3f> positions;
  bool result = false;
  for (size_t i = 0; i < triangle.size(); i++) {
    Vector3f pos, norm;

    subdiv.evaluateLimit(triangle[i], pos, norm);
    real_t k_r = evaluate_radial_curvature_analytic(mesh, subdiv, camera,
                                                    Vertex(), triangle[i]);

    real_t g = (camera.position() - pos).normalized().dot(norm);
    g_values.emplace_back(g);

    positions.emplace_back(pos);

    k_r_values.emplace_back(k_r);

    // Check zero-crossing
    std::vector<int> g_counts({0, 0}), k_r_counts({0, 0});
    bool g_near_zero = false;
    bool k_r_near_zero = false;
    for (size_t j = 0; j <= i; j++) {
      g_counts[(g_values[j] < 0) ? 0 : 1]++;
      k_r_counts[(k_r_values[j] < 0) ? 0 : 1]++;
      // Using 10 * CONTOUR_THRESHOLD here, since we are using inexact contour
      // points due to the curvature at vertices accuracy issue
      if (std::abs(g_values[j]) < CONTOUR_THRESHOLD)
        g_near_zero = true;
      if (std::abs(k_r_values[j]) < CONTOUR_THRESHOLD)
        k_r_near_zero = true;
    }

    if (((g_counts[0] > 0 && g_counts[1] > 0) || g_near_zero) &&
        ((k_r_counts[0] > 0 && k_r_counts[1] > 0) || k_r_near_zero))
      result = true;
  }

  if (tri_area) {
    *tri_area = (positions[1] - positions[0])
                    .cross(positions[2] - positions[0])
                    .norm() /
                2;
  }

  return result;
}

bool find_root_in_face(Mesh &mesh, Subdiv const &subdiv, Camera const &camera,
                       Face const &f, std::vector<Vector3f> &triangle_pos,
                       std::vector<Param_loc> &root_finding_triangle,
                       Vector3f &cusp_pos) {
  if (!has_zero_crossing(mesh, subdiv, camera, root_finding_triangle))
    return false;

  // Now there must be a zero crossing in the face
  size_t i = 0;
  real_t root_tri_area = (triangle_pos[1] - triangle_pos[0])
                             .cross(triangle_pos[2] - triangle_pos[0])
                             .norm() /
                         2;
  struct RootFindingState {
    std::vector<Vector3f> triangle_pos;
    std::vector<Param_loc> root_finding_triangle;
    real_t root_tri_area;
  };
  std::queue<RootFindingState> root_triangle_queue;
  root_triangle_queue.push(
      RootFindingState({triangle_pos, root_finding_triangle, root_tri_area}));
  bool found_cusp = false;
  size_t CUSP_MAX_ROOT_ITERATIONS = 10 * MAX_ROOT_ITERATIONS;
  real_t CUSP_THRESHOLD = 100 * CONTOUR_THRESHOLD;
  for (; i < CUSP_MAX_ROOT_ITERATIONS && !root_triangle_queue.empty(); i++) {
    RootFindingState state = root_triangle_queue.front();
    root_triangle_queue.pop();
    std::vector<Vector3f> cur_triangle_pos = state.triangle_pos;
    std::vector<Param_loc> cur_root_finding_triangle =
        state.root_finding_triangle;

    // Terminate condition: when the triangle is small enough
    cusp_pos =
        1. / 3 *
        (cur_triangle_pos[0] + cur_triangle_pos[1] + cur_triangle_pos[2]);

    bool terminates = state.root_tri_area < 5e-13;
    if (i > 0.9 * CUSP_MAX_ROOT_ITERATIONS && !terminates) {
      real_t max_k_g = -1;
      for (size_t i = 0; i < cur_root_finding_triangle.size(); i++) {
        Vector3f pos, norm;
        subdiv.evaluateLimit(cur_root_finding_triangle[i], pos, norm);
        real_t k_r = evaluate_radial_curvature_analytic(
            mesh, subdiv, camera, Vertex(), cur_root_finding_triangle[i]);

        real_t g = (camera.position() - pos).normalized().dot(norm);
        max_k_g = std::max(max_k_g, std::max(std::fabs(k_r), std::fabs(g)));
      }
      terminates = (max_k_g < CUSP_THRESHOLD);
    }

    if (terminates) {
      found_cusp = true;
      has_zero_crossing(mesh, subdiv, camera, cur_root_finding_triangle,
                        &root_tri_area);
      // Move the cusp to the limit surface
      Vector3f norm;
      std::vector<Vector3f> sub_pos;
      sub_pos.resize(3);
      subdiv.evaluateLimit(cur_root_finding_triangle[0], sub_pos[0], norm);
      subdiv.evaluateLimit(cur_root_finding_triangle[1], sub_pos[1], norm);
      subdiv.evaluateLimit(cur_root_finding_triangle[2], sub_pos[2], norm);

      root_finding_triangle = cur_root_finding_triangle;
      cusp_pos = 1. / 3 * (sub_pos[0] + sub_pos[1] + sub_pos[2]);
      break;
    }

    // Find the longest edge
    std::vector<real_t> tri_edge_lengths;
    for (size_t j = 0; j < cur_triangle_pos.size(); j++) {
      size_t next = (j + 1) % cur_triangle_pos.size();
      tri_edge_lengths.emplace_back(
          (cur_triangle_pos[next] - cur_triangle_pos[j]).norm());
    }
    auto max_ele =
        std::max_element(tri_edge_lengths.begin(), tri_edge_lengths.end());
    size_t max_e = std::distance(tri_edge_lengths.begin(), max_ele);
    size_t max_e_next = (max_e + 1) % cur_triangle_pos.size();

    // Find the sub-triangle that has the zero-crossing
    Param_loc bi_vertex(cur_root_finding_triangle[max_e].ptexIndex,
                        0.5f * cur_root_finding_triangle[max_e].uv +
                            0.5f * cur_root_finding_triangle[max_e_next].uv);
    std::vector<Param_loc> sub_triangle1({bi_vertex});
    std::vector<Param_loc> sub_triangle2({bi_vertex});
    std::vector<Vector3f> cur_triangle_pos1(
        {0.5f * cur_triangle_pos[max_e] + 0.5f * cur_triangle_pos[max_e_next]});
    std::vector<Vector3f> cur_triangle_pos2(
        {0.5f * cur_triangle_pos[max_e] + 0.5f * cur_triangle_pos[max_e_next]});
    for (size_t j = 0; j < cur_root_finding_triangle.size(); j++) {
      if (j != max_e && j != max_e_next) {
        sub_triangle1.emplace_back(cur_root_finding_triangle[j]);
        sub_triangle2.emplace_back(cur_root_finding_triangle[j]);
        cur_triangle_pos1.emplace_back(cur_triangle_pos[j]);
        cur_triangle_pos2.emplace_back(cur_triangle_pos[j]);
      } else if (j == max_e) {
        sub_triangle1.emplace_back(cur_root_finding_triangle[j]);
        cur_triangle_pos1.emplace_back(cur_triangle_pos[j]);
      } else if (j == max_e_next) {
        sub_triangle2.emplace_back(cur_root_finding_triangle[j]);
        cur_triangle_pos2.emplace_back(cur_triangle_pos[j]);
      }
    }

    if (has_zero_crossing(mesh, subdiv, camera, sub_triangle1,
                          &root_tri_area)) {
      root_triangle_queue.push(
          RootFindingState({cur_triangle_pos1, sub_triangle1, root_tri_area}));
    }
    if (has_zero_crossing(mesh, subdiv, camera, sub_triangle2,
                          &root_tri_area)) {
      root_triangle_queue.push(
          RootFindingState({cur_triangle_pos2, sub_triangle2, root_tri_area}));
    }
  }

  if (!found_cusp) {
    if (i == CUSP_MAX_ROOT_ITERATIONS)
      logger().warn("Cusp root finding reached MAX_ROOT_ITERATIONS {}",
                    CUSP_MAX_ROOT_ITERATIONS);
    return false;
  }

  return true;
}

bool insert_subdivision_cusps(
    Mesh &mesh, Subdiv const &subdiv, Camera const &camera, Face const &f,
    Edge const &e, std::vector<Vector3f> &inserted_vertices,
    std::vector<std::pair<Halfedge, real_t>> &edge_insertions) {
  auto param_loc = mesh.halfedge_property<Param_loc>("h:param_loc");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto is_cusp = mesh.get_vertex_property<bool>("v:cusp");
  if (is_cusp)
    is_cusp = mesh.vertex_property<bool>("v:cusp");

  // Root finding: g == 0 && kappa == 0
  std::map<int, Param_loc> f_vertex_param_map;
  std::vector<Halfedge> he;
  he.reserve(3);
  {
    auto hit = mesh.halfedges(f);
    auto hit_end = hit;
    do {
      he.emplace_back(*hit);
      f_vertex_param_map[mesh.from_vertex(*hit).idx()] = param_loc[*hit];
    } while (++hit != hit_end);
  }

  const Param_loc &v1_param_loc = param_loc[he[0]];
  const Param_loc &v2_param_loc = param_loc[he[1]];
  const Param_loc &v3_param_loc = param_loc[he[2]];
  real_t epsilon_offset = 0; // The current curvature has a bug at the
                             // vertex locations, so we offset them a bit
  std::vector<Vector3f> triangle_pos({vpositions[mesh.from_vertex(he[0])],
                                      vpositions[mesh.from_vertex(he[1])],
                                      vpositions[mesh.from_vertex(he[2])]});
  std::vector<Param_loc> root_finding_triangle({
      Param_loc(v1_param_loc.ptexIndex, (1 - epsilon_offset) * v1_param_loc.uv +
                                            epsilon_offset * v2_param_loc.uv),
      Param_loc(v2_param_loc.ptexIndex, (1 - epsilon_offset) * v2_param_loc.uv +
                                            epsilon_offset * v3_param_loc.uv),
      Param_loc(v3_param_loc.ptexIndex, (1 - epsilon_offset) * v3_param_loc.uv +
                                            epsilon_offset * v1_param_loc.uv),
  });
  std::vector<Param_loc> face_parameters = root_finding_triangle;
  std::vector<Vertex> face_vertices({mesh.from_vertex(he[0]),
                                     mesh.from_vertex(he[1]),
                                     mesh.from_vertex(he[2])});

  // Find root
  Vector3f cusp_pos;
  bool found_root = find_root_in_face(mesh, subdiv, camera, f, triangle_pos,
                                      root_finding_triangle, cusp_pos);
  if (!found_root)
    return false;

  Vector2f cusp_uv =
      1. / 3 *
      (root_finding_triangle[0].uv + root_finding_triangle[1].uv +
       root_finding_triangle[2].uv);

  // If the position is close to an existing vertex/edge, then simply label the
  // vertex and return.
  auto project2d = [&](Vector3f const &v3d, Vector2f &pos2D) {
    Vector2i viewport({camera.vpWidth(), camera.vpHeight()});
    pos2D = project(v3d, camera.viewMatrix().matrix(),
                    camera.projectionMatrix(), viewport)
                .head<2>();
  };
  {
    Vertex vv[3];
    mesh.verticesOfFace(f, vv);
    real_t min_dist = std::numeric_limits<real_t>::infinity();
    real_t min_dist_2d = std::numeric_limits<real_t>::infinity();
    Vertex cusp_v, cusp_v2;
    for (size_t i = 0; i < 3; i++) {
      auto const &v = vv[i];
      real_t d = (vpositions[v] - cusp_pos).norm();
      if (d < min_dist) {
        min_dist = d;
        cusp_v = v;
        cusp_v2 = Vertex();

        Vector2f v2d, cusp2d;
        project2d(vpositions[v], v2d);
        project2d(cusp_pos, cusp2d);
        min_dist_2d = (v2d - cusp2d).norm();
      }

      // Insert regardless of distance to edge
      if (e.is_valid())
        continue;

      // Distance to edge
      size_t j = (i + 1) % 3;
      auto const &v2 = vv[j];
      d = std::abs(
          (cusp_pos - vpositions[v]).cross(cusp_pos - vpositions[v2]).norm() /
          (vpositions[v2] - vpositions[v]).norm());
      Vector2f v2d, cusp2d, v2_2d;
      project2d(vpositions[v], v2d);
      project2d(cusp_pos, cusp2d);
      project2d(vpositions[v2], v2_2d);
      Vector2f norm_ij_2d = (v2_2d - v2d).normalized();
      norm_ij_2d = Vector2f(-norm_ij_2d.y(), norm_ij_2d.x());
      real_t dist_2d = std::abs(norm_ij_2d.dot(cusp2d - v2d));
      if (d < min_dist) {
        min_dist = d;
        cusp_v = v;
        cusp_v2 = v2;
        min_dist_2d = dist_2d;
      }
    }
    if (min_dist < 10 * CUSP_ROOT_THRESHOLD ||
        min_dist_2d < 10 * CUSP_ROOT_THRESHOLD) {
      if (is_cusp && !cusp_v2.is_valid())
        is_cusp[cusp_v] = true;
      return true;
    }
  }

  // If the position is close to an existing vertex/edge in UV space
  int closest_edge_index = -1;
  if (!e.is_valid()) {
    real_t min_dist_2d = std::numeric_limits<real_t>::infinity();
    Vertex cusp_v, cusp_v2;
    for (size_t i = 0; i < 3; i++) {
      auto const &param1 = face_parameters[i];

      // Distance to edge
      size_t j = (i + 1) % 3;
      auto const &param2 = face_parameters[j];
      Vector2f norm_ij_2d = (param2.uv - param1.uv).normalized();
      norm_ij_2d = Vector2f(-norm_ij_2d.y(), norm_ij_2d.x());
      real_t dist_2d = std::abs(norm_ij_2d.dot(cusp_uv - param1.uv));
      if (dist_2d < min_dist_2d) {
        closest_edge_index = i;
        min_dist_2d = dist_2d;
        cusp_v = face_vertices[i];
        cusp_v2 = face_vertices[j];
      }
    }

    if (min_dist_2d < 10 * CUSP_ROOT_THRESHOLD) {
      return true;
    }
  }

  // Check if we've inserted this cusp
  auto itr =
      std::find_if(inserted_vertices.begin(), inserted_vertices.end(),
                   [&cusp_pos](const Vector3f &p) {
                     return (cusp_pos - p).norm() < 10 * CUSP_ROOT_THRESHOLD;
                   });
  if (!inserted_vertices.empty() && itr != inserted_vertices.end())
    return false;
  inserted_vertices.emplace_back(cusp_pos);

  // Insert the cusp vertex and replace the contour edges
  Vertex cusp_v;
  Param_loc cusp_param;
  cusp_param.ptexIndex = root_finding_triangle[0].ptexIndex;
  cusp_param.uv = 1. / 3 *
                  (root_finding_triangle[0].uv + root_finding_triangle[1].uv +
                   root_finding_triangle[2].uv);
  if (e.is_valid()) {
    // if (1) {
    cusp_v = mesh.split(f, cusp_pos);
    if (is_cusp)
      is_cusp[cusp_v] = true;

    logger().info("\tFound a cusp: {}, {}", cusp_v, cusp_uv.transpose());
  } else {
    // Instead of spliting the face, split the edge to reduce the edge
    // increasing. (This is used in the first round of cusp insertion)
    contess_assert(closest_edge_index >= 0);
    size_t i = closest_edge_index;
    size_t j = (i + 1) % 3;

    size_t oppo_i = (i + 2) % 3;
    real_t t, u;
    Vector2f e1_uv = face_parameters[i].uv;
    Vector2f e2_uv = face_parameters[j].uv;
    Vector2f ray1 = face_parameters[oppo_i].uv +
                    10 * (cusp_uv - face_parameters[oppo_i].uv).normalized();
    Vector2f ray2 = face_parameters[oppo_i].uv -
                    10 * (cusp_uv - face_parameters[oppo_i].uv).normalized();
    bool found_intersection =
        (intersect2dSeg2dSegParametric(e1_uv.data(), e2_uv.data(), ray1.data(),
                                       ray2.data(), t, u, EPSILON) == 1);
    contess_assert_msg(found_intersection && t >= 0 && t <= 1,
                       "insert_subdivision_cusps: Cannot insert edge.");

    cusp_param.ptexIndex = face_parameters[i].ptexIndex;
    cusp_param.uv = (1 - t) * face_parameters[i].uv + t * face_parameters[j].uv;
    Edge closest_e = mesh.find_edge(face_vertices[i], face_vertices[j]);
    contess_assert_msg(closest_e.is_valid(),
                       "insert_subdivision_cusps: Edge not connected v" +
                           std::to_string(face_vertices[i].idx()) + ", v" +
                           std::to_string(face_vertices[j].idx()));
    if (!mesh.find_halfedge(face_vertices[i], face_vertices[j]).is_valid())
      return false;

    // Save the halfedge as going from the smaller vertex to the larger
    // one if it has a valid face.
    Vertex he_v1 = face_vertices[i];
    Vertex he_v2 = face_vertices[j];
    if (mesh.face(mesh.opposite_halfedge(mesh.find_halfedge(he_v1, he_v2)))
            .is_valid()) {
      if (face_vertices[i].idx() > face_vertices[j].idx())
        t = 1 - t;
      he_v1 = Vertex(std::min(face_vertices[i].idx(), face_vertices[j].idx()));
      he_v2 = Vertex(std::max(face_vertices[i].idx(), face_vertices[j].idx()));
    }
    edge_insertions.emplace_back(
        std::make_pair(mesh.find_halfedge(he_v1, he_v2), t));
    return true;
  }

  if (is_contour && e.is_valid()) {
    is_contour = mesh.edge_property<real_t>("e:contour");
    auto hit = mesh.halfedges(cusp_v);
    auto hit_end = hit;
    do {
      Edge e = mesh.edge(*hit);
      is_contour[e] = -1;
    } while (++hit != hit_end);
    Edge contour_e1 = mesh.find_edge(cusp_v, mesh.vertex(e, 0));
    Edge contour_e2 = mesh.find_edge(cusp_v, mesh.vertex(e, 1));
    is_contour[e] = -1;
    is_contour[contour_e1] = 1;
    is_contour[contour_e2] = 1;
  }

  // Update the subdivision parameterization
  // 1. From the original face vertices to the new cusp vertex
  for (auto &param_v : f_vertex_param_map) {
    Halfedge he = mesh.find_halfedge(Vertex(param_v.first), cusp_v);
    param_loc[he] = param_v.second;
  }

  // 2. From the new cusp vertex to the original face vertices
  {
    auto hit = mesh.halfedges(cusp_v);
    auto hit_end = hit;
    do {
      param_loc[*hit] = cusp_param;
    } while (++hit != hit_end);
  }

  // Flip the original contour edge to avoid the bad triangle containing two
  // contour edges
  if (e.is_valid()) {
    Vertex oppo_v;
    Param_loc oppo_param;
    {
      Halfedge a0 = mesh.halfedge(e, 0);
      Halfedge a1 = mesh.next_halfedge(a0);
      Halfedge a2 = mesh.next_halfedge(a1);

      if (mesh.to_vertex(a1) != cusp_v) {
        oppo_v = mesh.to_vertex(a1);
        oppo_param = param_loc[a2];
      }

      Halfedge b0 = mesh.halfedge(e, 1);
      Halfedge b1 = mesh.next_halfedge(b0);
      Halfedge b2 = mesh.next_halfedge(b1);
      if (mesh.to_vertex(b1) != cusp_v) {
        oppo_v = mesh.to_vertex(b1);
        oppo_param = param_loc[b2];
      }
    }

    mesh.flip(e);
    Halfedge flip_he1 = mesh.find_halfedge(cusp_v, oppo_v);
    Halfedge flip_he2 = mesh.find_halfedge(oppo_v, cusp_v);
    param_loc[flip_he1] = cusp_param;
    param_loc[flip_he2] = oppo_param;
  }

  return true;
}

void insert_subdivision_cusps(Mesh &mesh, Subdiv const &subdiv,
                              Camera const &camera) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto is_cusp = mesh.vertex_property<bool>("v:cusp");
  if (!is_cusp) {
    is_cusp = mesh.add_vertex_property<bool>("v:cusp");
  }
  is_cusp.vector().assign(is_cusp.vector().size(), false);

  // Get p info (we have one rotation index per p)
  std::vector<Edge> contour_edges;
  for (uint32_t i = 0; i < mesh.n_edges(); ++i) {
    Edge e(i);
    if (is_contour[e] >= 0) {
      contour_edges.emplace_back(e);
    }
  }
  std::unordered_set<int> visited_v;
  std::unordered_set<int> visited_f;
  std::unordered_set<int> potential_cusp_v;
  std::vector<Vector3f> inserted_vertices;
  for (auto const &e : contour_edges) {
    Face f1 = mesh.face(e, 0);
    Face f2 = mesh.face(e, 1);
    std::vector<std::pair<Halfedge, real_t>> edge_insertions;
    bool cusp_found = insert_subdivision_cusps(
        mesh, subdiv, camera, f1, e, inserted_vertices, edge_insertions);
    if (!cusp_found)
      cusp_found = insert_subdivision_cusps(mesh, subdiv, camera, f2, e,
                                            inserted_vertices, edge_insertions);
    // Check the faces adjacent to the vertex
    if (!cusp_found) {
      visited_f.emplace(f1.idx());
      visited_f.emplace(f2.idx());
      std::vector<Vertex> end_v({mesh.vertex(e, 0), mesh.vertex(e, 1)});

      for (auto const &v : end_v) {
        if (visited_v.find(v.idx()) != visited_v.end())
          continue;
        visited_v.emplace(v.idx());
        auto hit = mesh.halfedges(v);
        auto hit_end = hit;
        bool f_cusp = false;

        do {
          Face f = mesh.face(*hit);
          if (visited_f.find(f.idx()) != visited_f.end() || !f.is_valid())
            continue;
          f_cusp = face_contains_cusp(mesh, subdiv, camera, f);
          if (f_cusp)
            break;
        } while (++hit != hit_end);
        if (f_cusp) {
          potential_cusp_v.emplace(v.idx());
        }
      }
    }
  }

  // Check the potential cusps (caused by zero crossing not in the two adjacent
  // triangles, which may happen when the cusp is close to the edge endpoint and
  // it's also pointy)
  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");
  for (int v_idx : potential_cusp_v) {
    Vertex v(v_idx);
    std::vector<int> convexity_contour;
    convexity_contour.resize(2, 0);
    auto hit = mesh.halfedges(v);
    auto hit_end = hit;
    do {
      Edge e = mesh.edge(*hit);

      if (is_contour[e] < 0)
        continue;

      // Evaluate k_r at endpoints
      auto he1 = mesh.halfedge(e, 0);
      auto he2 = mesh.halfedge(e, 1);
      real_t kr_1 = evaluate_radial_curvature_analytic(
          mesh, mesh.const_subdivision(), camera, Vertex(), param_loc[he1]);
      real_t kr_2 = evaluate_radial_curvature_analytic(
          mesh, mesh.const_subdivision(), camera, Vertex(), param_loc[he2]);

      real_t kr = (std::fabs(kr_1) > std::fabs(kr_2)) ? kr_1 : kr_2;
      convexity_contour[int(kr < 0)]++;
    } while (++hit != hit_end);

    // Local minima along the contour edge
    if (convexity_contour[0] > 0 && convexity_contour[1] > 0) {
      is_cusp[v] = true;
    }
  }
}

/*===========================================================*/
bool face_contains_cusp(Mesh &mesh, Subdiv const &subdiv, Camera const &camera,
                        Face const &f, bool fast_check) {
  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_cusp = mesh.face_property<bool>("f:cusp");
  contess_assert(is_cusp);

  // Root finding: g == 0 && kappa == 0
  std::map<int, Param_loc> f_vertex_param_map;
  std::vector<Halfedge> he;
  he.reserve(3);
  {
    auto hit = mesh.halfedges(f);
    auto hit_end = hit;
    do {
      he.emplace_back(*hit);
      f_vertex_param_map[mesh.from_vertex(*hit).idx()] = param_loc[*hit];
    } while (++hit != hit_end);
  }

  const Param_loc &v1_param_loc = param_loc[he[0]];
  const Param_loc &v2_param_loc = param_loc[he[1]];
  const Param_loc &v3_param_loc = param_loc[he[2]];
  real_t epsilon_offset = 0; // The current curvature has a bug at the
                             // vertex locations, so we offset them a bit
  std::vector<Vector3f> triangle_pos({vpositions[mesh.from_vertex(he[0])],
                                      vpositions[mesh.from_vertex(he[1])],
                                      vpositions[mesh.from_vertex(he[2])]});
  std::vector<Param_loc> root_finding_triangle({
      Param_loc(v1_param_loc.ptexIndex, (1 - epsilon_offset) * v1_param_loc.uv +
                                            epsilon_offset * v2_param_loc.uv),
      Param_loc(v2_param_loc.ptexIndex, (1 - epsilon_offset) * v2_param_loc.uv +
                                            epsilon_offset * v3_param_loc.uv),
      Param_loc(v3_param_loc.ptexIndex, (1 - epsilon_offset) * v3_param_loc.uv +
                                            epsilon_offset * v1_param_loc.uv),
  });

  // Find root
  Vector3f cusp_pos;
  bool found_root = false;
  if (!fast_check) {
    found_root = find_root_in_face(mesh, subdiv, camera, f, triangle_pos,
                                   root_finding_triangle, cusp_pos);
  } else {
    found_root = has_zero_crossing(mesh, subdiv, camera, root_finding_triangle);
  }

  if (!found_root)
    return false;

  return true;
}

void insert_subdivision_cusp_edges(Mesh &mesh, Subdiv const &subdiv,
                                   Camera const &camera) {
  auto is_cusp = mesh.face_property<bool>("f:cusp");
  if (!is_cusp) {
    is_cusp = mesh.add_face_property<bool>("f:cusp", false);
  }
  is_cusp.vector().assign(is_cusp.vector().size(), false);
  for (uint32_t i = 0; i < mesh.n_faces(); ++i) {
    Face f(i);
    is_cusp[f] = face_contains_cusp(mesh, subdiv, camera, f);
  }
}
