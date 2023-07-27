// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "subdivide_contour_edges_even.h"
// Eigen printout formatting
#include <igl/predicates/predicates.h>
#include <limits>
#include <spdlog/fmt/ostr.h>
#include <unordered_map>
#include <vector>

#include "collapse_flipped_edges.h"
#include "common.h"
#include "fix_flipped_faces.h"
#include "insert_interpolated_contours.h"
#include "subdiv_common.h"
#include "subdivide_contour_edges.h"
#include "surface_mesh/surface_mesh.h"
#include "sweepLine.h"

inline Vector3f simple_project(const Vector3f &obj, const Matrix4f &model,
                               const Matrix4f &proj) {
  Vector4f tmp;
  tmp << obj, 1;

  tmp = model * tmp;
  tmp = proj * tmp;
  tmp = tmp.array() / tmp(3);

  return tmp.head(3);
}

void determine_t_point_2d(Mesh const &mesh, Camera const &camera,
                          Halfedge const &he, real_t t, real_t &mid_t) {
  if (std::abs(t - 1) < std::numeric_limits<real_t>::epsilon() ||
      std::abs(t) < std::numeric_limits<real_t>::epsilon()) {
    mid_t = t;
    return;
  }
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  Vector2i viewport = Vector2i(camera.vpWidth(), camera.vpHeight());

  Vertex v1 = mesh.from_vertex(he), v2 = mesh.to_vertex(he);
  Vector3f v1_pos = vpositions[v1];
  Vector3f v2_pos = vpositions[v2];
  Vector2f v1_2d =
      simple_project(vpositions[mesh.from_vertex(he)],
                     camera.viewMatrix().matrix(), camera.projectionMatrix())
          .head<2>();
  Vector2f v2_2d =
      simple_project(vpositions[mesh.to_vertex(he)],
                     camera.viewMatrix().matrix(), camera.projectionMatrix())
          .head<2>();

  Vector2f mid_2d = v1_2d + t * (v2_2d - v1_2d);
  Matrix4f P = camera.projectionMatrix() * camera.viewMatrix().matrix();
  real_t B_x = mid_2d[0];
  auto P1 = P.row(0).head(3);
  auto P4 = P.row(3).head(3);
  real_t p1 = P.row(0).tail(1)[0];
  real_t p4 = P.row(3).tail(1)[0];
  real_t lhs = (P1 - B_x * P4) * (v2_pos - v1_pos);
  real_t rhs = B_x * (P4 * v1_pos)[0] + B_x * p4 - P1 * v1_pos - p1;
  mid_t = rhs / lhs;
  if (std::isnan(mid_t))
    mid_t = t;
}

void determine_search_direction(Mesh &mesh, Camera const &camera, Edge const &e,
                                real_t t, Face const &f, Halfedge &end_he,
                                real_t &end_t) {
  end_he = Halfedge();
  auto vpositions = mesh.vertex_property<Vector3f>("v:point");
  Vector2i viewport = Vector2i(camera.vpWidth(), camera.vpHeight());

  // Project the edge to 2D
  Vector2f v1_2d =
      project(vpositions[mesh.vertex(e, 0)], camera.viewMatrix().matrix(),
              camera.projectionMatrix(), viewport)
          .head<2>();
  Vector2f v2_2d =
      project(vpositions[mesh.vertex(e, 1)], camera.viewMatrix().matrix(),
              camera.projectionMatrix(), viewport)
          .head<2>();
  Vector2f e_tan_2d = (v2_2d - v1_2d).normalized();
  Vector2f mid_2d = v1_2d + t * (v2_2d - v1_2d);

  // 2D search direction
  Vector2f e_norm_2d(-e_tan_2d.y(), e_tan_2d.x());
  Vector2f search_v1 = mid_2d + 1e6 * e_norm_2d; // Simply make it very long...
  Vector2f search_v2 = mid_2d - 1e6 * e_norm_2d; // Simply make it very long...

  // Determine the endpoint for the search direction by 2D segment intersection
  Halfedge f_he = mesh.halfedge(e, 0);
  if (mesh.face(f_he) != f)
    f_he = mesh.halfedge(e, 1);
  Vertex v3 = mesh.to_vertex(mesh.next_halfedge(f_he));
  Vector2f v3_2d = project(vpositions[v3], camera.viewMatrix().matrix(),
                           camera.projectionMatrix(), viewport)
                       .head<2>();
  std::vector<Vertex> endpoints({mesh.vertex(e, 0), mesh.vertex(e, 1)});
  std::vector<Vector2f> endpoints_2d({v1_2d, v2_2d});
  for (size_t i = 0; i < endpoints.size(); i++) {
    double tt, uu;
    if (intersect2dSeg2dSegParametric(search_v1.data(), search_v2.data(),
                                      v3_2d.data(), endpoints_2d[i].data(), tt,
                                      uu, EPSILON) == 1) {
      end_he = mesh.find_halfedge(v3, endpoints[i]);
      end_t = uu;
      break;
    }
  }

  contess_assert(end_he.is_valid());
}

FacingType evaluate_g_value_param(Mesh &mesh, Camera const &camera,
                                  Param_loc const &param_res) {
  Vector3f pos_res, normal_res;

  mesh.subdivision().evaluateLimit(param_res, pos_res, normal_res);
  real_t ndotv = normal_res.dot((camera.position() - pos_res).normalized());

  if (ndotv < -CONTOUR_THRESHOLD)
    return BACK;

  if (ndotv > CONTOUR_THRESHOLD)
    return FRONT;

  return CONTOUR;
}

void subdivide_contour_edges_even(Mesh &mesh, Camera const &camera,
                                  std::vector<Edge> &subdiv_edges) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto vpositions = mesh.vertex_property<Vector3f>("v:point");
  auto param_loc = mesh.halfedge_property<Param_loc>("h:param_loc");
  auto facing = mesh.vertex_property<FacingType>("v:facing");
  auto ndotv_v = mesh.vertex_property<real_t>("v:ndotv");
  auto vnormals = mesh.get_vertex_property<Vector3f>("v:normal");

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

  auto get_valid_param =
      [&](Vertex const &v,
          std::unordered_map<int, std::vector<Param_loc>> &param_matching) {
        std::unordered_map<int, Param_loc> local_param_matching;
        auto hit = mesh.halfedges(v);
        auto hit_end = hit;
        do {
          if (!param_loc[*hit].is_valid()) {
            continue;
          }
          if (!local_param_matching.count(param_loc[*hit].ptexIndex)) {
            local_param_matching[param_loc[*hit].ptexIndex] = param_loc[*hit];
          }
        } while (++hit != hit_end);

        for (auto const &local_param : local_param_matching) {
          if (!param_matching.count(local_param.first)) {
            param_matching[local_param.first] = std::vector<Param_loc>();
          }
          param_matching[local_param.first].emplace_back(local_param.second);
        }
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

  struct SearchEnd {
    SearchEnd(Halfedge in_search_end_he, real_t in_end_t)
        : search_end_he(in_search_end_he), end_t(in_end_t) {}
    Halfedge search_end_he;
    real_t end_t;
  };

  for (auto const &e : subdiv_edges) {
    contess_assert_msg(is_contour[e] >= 0 || mesh.is_boundary(e),
                       "Cannot subdivide non contour edges.");

    // 1. Unset the contour label
    is_contour[e] = -1;

    // 2. Split edge in the middle, determine the side containing the
    // zero-crossing
    Vertex v1 = mesh.vertex(e, 0), v2 = mesh.vertex(e, 1);

    // Correct for perspective (though the effect of this correction is quite
    // minor)
    real_t mid_t = 0.5;
    determine_t_point_2d(mesh, camera, mesh.halfedge(e, 0), mid_t, mid_t);

    Vector3f mid = vpositions[v1] + mid_t * (vpositions[v2] - vpositions[v1]);
    Param_loc v1_param = param_loc[mesh.find_halfedge(v1, v2)];
    Param_loc v2_param =
        param_loc[mesh.next_halfedge(mesh.find_halfedge(v1, v2))];
    Param_loc v_param(v1_param.ptexIndex,
                      (1 - mid_t) * v1_param.uv + mid_t * v2_param.uv);

    // Record the parameter values for the new vertex
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

    // 3. Before inserting the upsampling vertex, determine the two
    // search directions
    std::vector<SearchEnd> search_ends;
    FacingType v_facing =
        (!is_e_boundary)
            ? evaluate_g_value_param(mesh, camera, v_param)
            : ((facing[v1] == FacingType::CONTOUR) ? facing[v2] : facing[v1]);
    Param_loc search_param;

    if (!is_e_boundary) {
      int flipping_round = 0;
      Param_loc param_v1;
      Param_loc param_v2;
      real_t end_t;

      do {
        search_ends.clear();

        Halfedge end_he;

        determine_search_direction(mesh, camera, e, mid_t, mesh.face(e, 0),
                                   end_he, end_t);
        search_ends.emplace_back(end_he, end_t);
        determine_search_direction(mesh, camera, e, mid_t, mesh.face(e, 1),
                                   end_he, end_t);
        search_ends.emplace_back(end_he, end_t);

        // 4. Determine between the two search directions: which one contains a
        // zero-crossing

        for (auto const &search_end : search_ends) {
          Halfedge end_he = search_end.search_end_he;
          end_t = search_end.end_t;

          std::unordered_map<int, std::vector<Param_loc>> param_matching;
          get_valid_param(mesh.from_vertex(end_he), param_matching);
          get_valid_param(mesh.to_vertex(end_he), param_matching);
          for (auto const &match_param : param_matching) {
            if (match_param.second.size() == 2 &&
                match_param.first == v_param.ptexIndex) {
              param_v1 = match_param.second.front();
              param_v2 = match_param.second.back();
              break;
            }
          }
          contess_assert(param_v1.is_valid() && param_v2.is_valid());

          Param_loc param_search_end;
          param_search_end = param_v1;
          param_search_end.uv =
              param_v1.uv + end_t * (param_v2.uv - param_v1.uv);

          FacingType v_n_facing =
              evaluate_g_value_param(mesh, camera, param_search_end);

          FacingType v_facing1 = evaluate_g_value_param(mesh, camera, param_v1);

          if (v_facing != v_n_facing && v_facing != FacingType::CONTOUR &&
              v_n_facing != FacingType::CONTOUR) {
            search_param = param_search_end;
            break;
          }
          // There's a zero crossing (caused by multiple zero-crossings in the
          // initial mesh?)
          else if (v_facing1 != v_n_facing &&
                   v_n_facing != FacingType::CONTOUR) {
            if (flipping_round > 0) {
              flipping_round++;
              break;
            }
            // Flipping
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
              if (!existing_param.is_valid()) {
                hit = mesh.halfedges(mesh.from_vertex(he));
                hit_end = hit;
                do {
                  if (*hit == he)
                    continue;
                  if (param_loc[*hit].is_valid()) {
                    existing_param = param_loc[*hit];
                    break;
                  }
                } while (++hit != hit_end);
              }
              param_loc[he] = existing_param;
            };
            Edge flip_e = mesh.edge(end_he);
            Param_loc before_flip = param_loc[mesh.halfedge(flip_e, 0)];
            Vertex vv1 =
                mesh.to_vertex(mesh.next_halfedge(mesh.halfedge(flip_e, 0)));
            Vertex vv2 =
                mesh.to_vertex(mesh.next_halfedge(mesh.halfedge(flip_e, 1)));

            // Flipping check if we are flipping param quad boundary
            bool can_flip = can_flip_edge(mesh, flip_e);
            logger().info("Flipping edge: {} - {}, {}", flip_e,
                          mesh.vertex(flip_e, 0).idx(),
                          mesh.vertex(flip_e, 1).idx());
            if (!can_flip) {
              flipping_round = 2;
              break;
            }
            mesh.flip(flip_e);

            // Update parameter
            update_param(mesh.find_halfedge(vv1, vv2), before_flip.ptexIndex);
            update_param(mesh.find_halfedge(vv2, vv1), before_flip.ptexIndex);

            flipping_round++;
            break;
          } else if (flipping_round > 0) {
            flipping_round++;
            break;
          }

          if (search_param.is_valid() || flipping_round >= 2)
            break;
        }

        // Already flipped OR can't flip at all
        if (flipping_round >= 2 || flipping_round == 0)
          break;
      } while (!search_param.is_valid());

      if (!search_param.is_valid()) {
        int substep_number = 10;
        std::vector<real_t> valid_ts;
        real_t t_dist = 1;
        for (int j = 1; j + 1 < substep_number; j++) {
          real_t sub_t = (real_t)(j) / substep_number;
          Param_loc param_search_sub;
          param_search_sub = param_v1;
          param_search_sub.uv =
              param_v1.uv + sub_t * (param_v2.uv - param_v1.uv);

          FacingType sub_facing =
              evaluate_g_value_param(mesh, camera, param_search_sub);
          if (v_facing != sub_facing && v_facing != FacingType::CONTOUR &&
              sub_facing != FacingType::CONTOUR) {
            real_t dist = std::abs(sub_t - end_t);
            if (dist < t_dist) {
              t_dist = dist;
              search_param = param_search_sub;
            }
          }
        }
      }
    }

    // 5. Insert to the current contour
    Vertex v = mesh.split(e, mid);
    facing[v] = v_facing;

    if (!search_param.is_valid() && !is_e_boundary) {
      if (abs(ndotv_v[v]) > CONTOUR_THRESHOLD)
        logger().warn("Cannot find zero crossing at {} - {} - {}", v1.idx(),
                      v.idx(), v2.idx());
      contess_assert_msg(abs(ndotv_v[v]) <= CONTOUR_THRESHOLD,
                         "No zero-crossing found in the new edges...");

      is_contour[mesh.find_edge(v, v1)] = 1;
      is_contour[mesh.find_edge(v, v2)] = 1;
      continue;
    }

    auto hit = mesh.halfedges(v);
    auto hit_end = hit;
    do {
      // Update the parameter
      int n_ptexIndex = param_loc[mesh.next_halfedge(*hit)].ptexIndex;
      param_loc[*hit] = updated_params[n_ptexIndex];

      param_loc[mesh.opposite_halfedge(*hit)] =
          updated_params_v[mesh.from_vertex(mesh.opposite_halfedge(*hit))
                               .idx()];
    } while (++hit != hit_end);

    // Handle boundary: Directly move the point to the limit position
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
      logger().info("Subdiv boundary v: {}", v);

      continue;
    }

    // 6. Move the inserted vertex
    // Update the parameters going out of this vertex
    Vector3f p, n;
    Param_loc param;
    mesh.bisect_search(camera.position(), &mesh.const_subdivision(), v_param,
                       search_param, p, n, param);
    is_contour[mesh.find_edge(v, v1)] = 1;
    is_contour[mesh.find_edge(v, v2)] = 1;

    vpositions[v] = p;
    vnormals[v] = n;
    ndotv_v[v] = 0;
  }

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

  // Repair parameter
  // Hacky fix to the parameterization...
  // Error can only be introduced on halfedge from an existing vertex to a
  // newly-added vertex (which would have parameter on all outward he and they
  // should have the same quad index). So we shouldn't touch parameter quad
  // boundary and this simple fix should work...
  for (auto hit = mesh.halfedges_begin(); hit != mesh.halfedges_end(); ++hit)
    if (!mesh.is_boundary(*hit) && param_loc[*hit].ptexIndex == -1) {
      std::unordered_map<int, Param_loc> updated_params_v;
      int ref_ptex_idx = param_loc[mesh.opposite_halfedge(*hit)].ptexIndex;
      auto v1 = mesh.from_vertex(*hit);
      auto hit = mesh.halfedges(v1);
      auto hit_end = hit;
      do {
        // Update the parameter
        if (param_loc[*hit].is_valid() &&
            (ref_ptex_idx < 0 || param_loc[*hit].ptexIndex == ref_ptex_idx)) {
          updated_params_v[mesh.from_vertex(*hit).idx()] = param_loc[*hit];
          break;
        }
      } while (++hit != hit_end);
      hit = mesh.halfedges(v1);
      hit_end = hit;
      do {
        // Update the parameter
        if (!param_loc[*hit].is_valid()) {
          param_loc[*hit] = updated_params_v[mesh.from_vertex(*hit).idx()];
        }
      } while (++hit != hit_end);
    }

  fix_flipped_faces(mesh, camera);
  return;
}

void subdivide_contour_intersection_edges_even(Mesh &mesh, Camera const &camera,
                                               Vertex const &intersection_v,
                                               size_t subdiv_times) {
  auto is_contour = mesh.edge_property<real_t>("e:contour");
  auto vpositions = mesh.vertex_property<Vector3f>("v:point");
  auto param_loc = mesh.halfedge_property<Param_loc>("h:param_loc");
  auto vnormals = mesh.get_vertex_property<Vector3f>("v:normal");
  auto ndotv_v = mesh.vertex_property<real_t>("v:ndotv");

  // First reset book keeping for the paired intersection vertex
  // Note the left intersection vertex is currently dangling...
  auto intersection_2d = mesh.get_vertex_property<int>("v:intersection_2d");
  auto is_valid_intersection_2d =
      mesh.get_vertex_property<bool>("v:is_valid_intersection_2d");
  Vertex pair_v(intersection_2d[intersection_v]);
  is_valid_intersection_2d[pair_v] = false;
  is_valid_intersection_2d[intersection_v] = false;

  // Move the intersection v to contour
  Param_loc v_param = param_loc[*mesh.halfedges(intersection_v)];
  auto v_facing = evaluate_g_value_param(mesh, camera, v_param);

  bool found_root = false;
  auto hit = mesh.halfedges(intersection_v), hit_end = hit;
  std::vector<Halfedge> contour_he;
  do {
    // Search the root on non-contour edges
    if (!found_root && is_contour[mesh.edge(*hit)] < 0 &&
        !mesh.is_boundary(mesh.edge(*hit))) {
      Vector3f p, n;
      Param_loc search_param = param_loc[mesh.opposite_halfedge(*hit)];
      auto search_facing = evaluate_g_value_param(mesh, camera, search_param);

      if (v_facing != FacingType::CONTOUR && v_facing != search_facing &&
          search_facing != FacingType::CONTOUR) {
        Param_loc param;
        mesh.bisect_search(camera.position(), &mesh.const_subdivision(),
                           v_param, search_param, p, n, param);

        found_root = true;
        vpositions[intersection_v] = p;
        vnormals[intersection_v] = n;
        ndotv_v[intersection_v] = 0;
      }
    }
    if (is_contour[mesh.edge(*hit)] < 0 && !mesh.is_boundary(mesh.edge(*hit)))
      continue;
    contour_he.emplace_back(*hit);
  } while (++hit != hit_end);

  // Upsample the adjacent contour edges
  for (auto const &he : contour_he) {
    Edge e = mesh.edge(he);
    std::vector<Edge> subdiv_edges({e});
    bool is_even = true;
    subdivide_contour_edges_k_times(mesh, camera, subdiv_edges, subdiv_times,
                                    is_even, -1);
  }
}
