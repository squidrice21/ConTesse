// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "add_simplification_vertex.h"
#include "common.h"
#include "cut_patch_to_disk.h"
#include "insert_interpolated_contours.h"
#include "subdivide_contour_edges.h"
#include "subdivide_contour_edges_even.h"
#include <limits>
#include <spdlog/fmt/ostr.h>
#include <string>
#include <unordered_set>
#include <valarray>

Vertex add_contour_vertex(Mesh &mesh, Camera const &camera, Halfedge const &he,
                          real_t const &t) {
  auto vnormals = mesh.vertex_property<Vector3f>("v:normal");
  auto ndotv = mesh.vertex_property<real_t>("v:ndotv");
  auto concave = mesh.get_edge_property<bool>("e:concave");
  auto is_contour = mesh.edge_property<real_t>("e:contour");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  if (patchID)
    patchID = mesh.face_property<int>("f:patchID");

  contess_assert_msg(is_contour,
                     "insert_planar_map_intersections: Contours are missing. "
                     "Planar map intersections can't be created.");

  auto param_loc = mesh.halfedge_property<Param_loc>("h:param_loc");

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
    contess_assert(0);
  };

  // Getting all the intersection point positions
  Edge orig_e = mesh.edge(he);
  Vertex init_v = mesh.from_vertex(he);
  Vertex final_v = mesh.to_vertex(he);
  Edge to_split = orig_e;

  bool originally_contour = is_contour[to_split] >= 0;

  surface_mesh::Point p =
      (1 - t) * vpositions[init_v] + t * vpositions[final_v];
  std::unordered_set<int> split_vertices;
  split_vertices.insert(mesh.vertex(to_split, 0).idx());
  split_vertices.insert(mesh.vertex(to_split, 1).idx());

  // For parameter updates
  Vertex v1 = mesh.vertex(to_split, 0), v2 = mesh.vertex(to_split, 1);
  real_t param_t =
      (p - vpositions[v1]).norm() / (vpositions[v2] - vpositions[v1]).norm();
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
  //

  int is_concave = -1;
  if (concave) {
    is_concave = concave[to_split];
  }
  Vertex new_v = mesh.split(to_split, p);

  // Update the edge properties (parameter, convex)
  auto hit = mesh.halfedges(new_v);
  auto hit_end = hit;
  Param_loc param_res;
  do {
    if (originally_contour && split_vertices.find(mesh.to_vertex(*hit).idx()) !=
                                  split_vertices.end()) {
      is_contour[mesh.edge(*hit)] = 1;
      if (is_concave >= 0)
        concave[mesh.edge(*hit)] = is_concave;
    }

    int n_ptexIndex = param_loc[mesh.next_halfedge(*hit)].ptexIndex;
    param_loc[*hit] = updated_params[n_ptexIndex];
    param_res = param_loc[*hit];

    param_loc[mesh.opposite_halfedge(*hit)] =
        updated_params_v[mesh.from_vertex(mesh.opposite_halfedge(*hit)).idx()];
  } while (++hit != hit_end);

  if (!param_res.is_valid())
    param_res = v_param_valid;

  contess_assert_msg(param_res.is_valid(),
                     "Parameter error at the new vertex.");

  // Update facing information
  {
    Vector3f pos_res, normal_res;
    mesh.subdivision().evaluateLimit(param_res, pos_res, normal_res);
    real_t v_ndotv = normal_res.dot((camera.position() - pos_res).normalized());
    ndotv[new_v] = v_ndotv;
    vnormals[new_v] = normal_res;
  }

  return new_v;
}

void find_2d_projection(Mesh const &mesh, Camera const &camera, Vertex const &v,
                        Halfedge const &he, real_t &t) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  Vector2i viewport({camera.vpWidth(), camera.vpHeight()});
  auto project2d = [&](Vector3f const &v3d, Vector2f &pos2D) {
    pos2D = project(v3d, camera.viewMatrix().matrix(),
                    camera.projectionMatrix(), viewport)
                .head<2>();
  };

  Vector3f p = vpositions[v];
  Vector2f p_2d;
  project2d(p, p_2d);
  Vector3f p1 = vpositions[mesh.from_vertex(he)];
  Vector2f p1_2d;
  project2d(p1, p1_2d);
  Vector3f p2 = vpositions[mesh.to_vertex(he)];
  Vector2f p2_2d;
  project2d(p2, p2_2d);

  Vector2f tan = (p2_2d - p1_2d).normalized();
  Vector2f proj = (p_2d - p1_2d).dot(tan) * tan + p1_2d;
  real_t temp_t = ((proj - p1_2d).norm() / (p2_2d - p1_2d).norm());
  temp_t = std::min(std::max(temp_t, 0.), 1.);

  // Correct for perspective
  determine_t_point_2d(mesh, camera, he, temp_t, t);
}

Vertex add_simplification_vertex(Mesh &mesh, Camera const &camera,
                                 Vertex const &cusp_v, Halfedge const &he_cusp,
                                 Halfedge const &he_line, Vector3f &moved_pos) {
  auto patchID = mesh.get_face_property<int>("f:patchID");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto patchBoundary = mesh.get_halfedge_property<int>("h:patchBoundary");
  contess_assert_msg(patchBoundary,
                     "tag_2d_intersections: Patch needs to be built.");

  Vector2i viewport({camera.vpWidth(), camera.vpHeight()});
  auto project2d = [&](Vector3f const &v3d, Vector2f &pos2D) {
    pos2D = project(v3d, camera.viewMatrix().matrix(),
                    camera.projectionMatrix(), viewport)
                .head<2>();
  };

  // 1. Find the edge that the cusp projects to (in 2D) and the projection
  // position
  Vector3f cusp_p = vpositions[cusp_v];
  Vector2f cusp_p_2d;
  project2d(cusp_p, cusp_p_2d);
  Vertex closest_end_v;
  Vertex projection_v;
  real_t min_end_dist = std::numeric_limits<real_t>::infinity();
  real_t proj_t = -1;
  Halfedge walk_he = he_line;
  std::unordered_set<size_t> seen_vertices;
  do {
    // Find 2D projection
    find_2d_projection(mesh, camera, cusp_v, walk_he, proj_t);

    // Determine the projection type (end? mid?)
    Vector3f end1_p = vpositions[mesh.from_vertex(walk_he)];
    Vector3f end2_p = vpositions[mesh.to_vertex(walk_he)];
    Vector3f proj_p = (1 - proj_t) * end1_p + proj_t * end2_p;
    Vector2f proj_p_2d;
    project2d(proj_p, proj_p_2d);

    if (seen_vertices.count(mesh.to_vertex(walk_he).idx()))
      return Vertex();
    seen_vertices.emplace(mesh.to_vertex(walk_he).idx());

    bool is_end_proj = true;
    auto end1_ori =
        igl::predicates::orient3d(camera.position(), proj_p, cusp_p, end1_p);
    auto end2_ori =
        igl::predicates::orient3d(camera.position(), proj_p, cusp_p, end2_p);
    if ((end1_ori == igl::predicates::Orientation::NEGATIVE &&
         end2_ori == igl::predicates::Orientation::POSITIVE) ||
        (end1_ori == igl::predicates::Orientation::POSITIVE &&
         end2_ori == igl::predicates::Orientation::NEGATIVE))
      is_end_proj = false;

    if (is_end_proj) {
      real_t end_dist = (proj_p_2d - cusp_p_2d).norm();
      if (end_dist < min_end_dist) {
        min_end_dist = end_dist;
        closest_end_v = (proj_t < 0.5) ? mesh.from_vertex(walk_he)
                                       : mesh.to_vertex(walk_he);
      }
    } else {
      // 2. Insert the projection
      projection_v = add_contour_vertex(mesh, camera, walk_he, proj_t);

      // 3. Move the projection to the cusp within the image plane.
      Vector3f cam_ray = (vpositions[cusp_v] - camera.position()).normalized();
      Vector3f added_proj =
          (vpositions[projection_v] - camera.position()).dot(cam_ray) *
              cam_ray +
          camera.position();
      moved_pos = added_proj;

      break;
    }

    // Walk
    Vertex he_v = mesh.to_vertex(walk_he);
    auto hit = mesh.halfedges(he_v), hit_end = hit;
    Halfedge next_he;
    do {
      if ((mesh.is_boundary(mesh.edge(*hit)) ||
           is_contour[mesh.edge(*hit)] >= 0) &&
          patchBoundary[*hit] == patchBoundary[walk_he]) {
        next_he = *hit;
        break;
      }
    } while (++hit != hit_end);

    contess_assert_msg(next_he.is_valid(),
                       "add_simplification_vertex: Cannot walk from v" +
                           std::to_string(he_v.idx()));
    walk_he = next_he;
  } while (std::abs(proj_t) > std::numeric_limits<real_t>::epsilon());
  // Until the cusp is back wrt the walking direction on the other side

  // Fill patch ID for new faces
  if (patchID) {
    std::vector<Face> fill_patch_faces;
    for (size_t i = 0; i < mesh.n_faces(); i++) {
      Face f(i);

      if (patchID[f] == -1)
        fill_patch_faces.emplace_back(f);
    }
    update_patch(mesh, camera, fill_patch_faces);
    consistently_label_interpolated_contours(mesh, camera.position(), false);
    mesh.markPatchBoundaryEdges(false);
  }

  Vertex final_v = (projection_v.is_valid()) ? projection_v : closest_end_v;
  return final_v;
}

real_t tag_cusp_line_case(Mesh &mesh, Camera const &camera,
                          Vertex const &cusp_v, Vertex const &added_v,
                          Halfedge const &he_cusp, Halfedge const &he_line,
                          bool to_write) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto intersection_2d = mesh.vertex_property<int>("v:intersection_2d");
  contess_assert_msg(
      intersection_2d,
      "tag_unmoveable_vertices: 2D intersections are not created.");
  auto is_valid_intersection_2d =
      mesh.vertex_property<bool>("v:is_valid_intersection_2d");
  contess_assert_msg(
      is_valid_intersection_2d,
      "tag_unmoveable_vertices: Valid 2D intersections are not tagged.");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  contess_assert_msg(is_contour, "Requires contours.");
  auto patchBoundary = mesh.get_halfedge_property<int>("h:patchBoundary");
  contess_assert_msg(patchBoundary,
                     "tag_2d_intersections: Patch needs to be built.");

  auto feasibility_collapsing =
      mesh.edge_property<int>("e:feasibility_collapsing");

  intersection_2d[cusp_v] = added_v.idx();
  intersection_2d[added_v] = cusp_v.idx();
  is_valid_intersection_2d[cusp_v] = is_valid_intersection_2d[added_v] = true;

  real_t walk_length = 0;

  // Walk to tag collapsing
  auto tag_collapsing = [&](Halfedge const &he, Vertex const &end_v) {
    contess_assert(mesh.from_vertex(he) != end_v);

    Vertex he_v = mesh.from_vertex(he);
    Halfedge next_he;
    do {
      auto hit = mesh.halfedges(he_v), hit_end = hit;
      do {
        if ((mesh.is_boundary(mesh.edge(*hit)) ||
             is_contour[mesh.edge(*hit)] >= 0) &&
            patchBoundary[*hit] == patchBoundary[he]) {
          if (to_write)
            feasibility_collapsing[mesh.edge(*hit)] =
                mesh.from_vertex(he).idx();
          auto v1 =
              project(vpositions[mesh.from_vertex(*hit)],
                      camera.viewMatrix().matrix(), camera.projectionMatrix(),
                      Vector2i(camera.vpWidth(), camera.vpHeight()));
          auto v2 =
              project(vpositions[mesh.to_vertex(*hit)],
                      camera.viewMatrix().matrix(), camera.projectionMatrix(),
                      Vector2i(camera.vpWidth(), camera.vpHeight()));
          walk_length += (v1 - v2).norm();
          next_he = *hit;
          break;
        }
      } while (++hit != hit_end);
      he_v = mesh.to_vertex(next_he);
    } while (he_v != end_v);
  };
  tag_collapsing(he_cusp, cusp_v);
  tag_collapsing(he_line, added_v);

  return walk_length;
}
