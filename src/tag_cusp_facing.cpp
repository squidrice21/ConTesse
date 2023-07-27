// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "tag_cusp_facing.h"
#include "common.h"

#include <limits>
#include <spdlog/fmt/ostr.h>
#include <unordered_set>
#include <vector>

#include "inflat_path.h"
#include "ray_cast_QI.h"

namespace predicates = igl::predicates;

FacingType walk_on_original_facing(Mesh const &mesh, Vertex src_orig,
                                   Vector3f const &cut_p,
                                   Vector3f const &cut_norm,
                                   Camera const &camera,
                                   Vector3f const &laplacian) {
  FacingType facing = FacingType::UNDEFINED;
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");
  auto VBO_f = mesh.get_face_property<FacingType>("f:VBO_f");

  std::unordered_set<int> visited_face_indices;
  std::queue<Face> face_frontier;

  // Initialize with all neighboring faces
  {
    auto hit = mesh.halfedges(src_orig), hit_end = hit;
    do {
      // Check adjacent faces
      auto adj_f_h = mesh.opposite_halfedge(*hit);
      Face adj_f = mesh.face(adj_f_h);

      // In case we are at the boundary of the model
      if (adj_f.is_valid())
        face_frontier.push(adj_f);
    } while (++hit != hit_end);
  }

  if (face_frontier.empty())
    return facing;

  auto is_destination = [&](Face f) -> bool {
    Vertex vv[3];
    mesh.verticesOfFace(f, vv);

    bool correct_laplacian_side = false;
    real_t far_lap_abs = -1;
    real_t far_dot_prod = 0;
    std::vector<int> vv_signs;
    vv_signs.resize(2, 0);
    for (auto const &v : vv) {
      // Check if the two overlaps in 2D
      if (v == src_orig)
        continue;
      real_t dot_prod = (vpositions[v] - vpositions[src_orig]).dot(laplacian);
      if (std::abs(dot_prod) > far_lap_abs) {
        far_lap_abs = std::abs(dot_prod);
        far_dot_prod = dot_prod;
      }
      vv_signs[(dot_prod > 0) ? 0 : 1]++;
    }

    if (vv_signs[0] + vv_signs[1] != 3 || (vv_signs[0] > 0 && vv_signs[1] > 0))
      return false;

    if (far_dot_prod > 0) {
      correct_laplacian_side = true;
    }

    if (VBO[f] == VBO_f[f] && correct_laplacian_side)
      return true;

    return false;
  };

  // Walk
  std::unordered_map<int, int> prev_idx;
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
    if (planeIntersectTri(cut_p, cut_norm, tri_p)) {
      // Extend frontier
      auto fhit = mesh.halfedges(cur_f), fhit_end = fhit;
      do {
        // Sample on the edges of the current face
        // and determine the intersecting edge for the next step
        Vector3f intersection;
        bool edge_intersected = planeIntersectEdge(
            cut_p, cut_norm, vpositions[mesh.from_vertex(*fhit)],
            vpositions[mesh.to_vertex(*fhit)], intersection);
        // Don't grow the non-intersecting edge
        if (!edge_intersected)
          continue;

        f_intersection[cur_f.idx()] = intersection;

        // We can't cross the contour edge
        // And cut edge
        if (is_contour[mesh.edge(*fhit)] >= 0)
          continue;

        // Check adjacent faces
        auto adj_f_h = mesh.opposite_halfedge(*fhit);
        Face adj_f = mesh.face(adj_f_h);
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
        facing = VBO[cur_f];

        break;
      }
    }

    visited_face_indices.insert(cur_f.idx());
  } while (!face_frontier.empty());

  return facing;
}

Halfedge
find_nondegenerated_contour_halfedge(Mesh const &mesh, Vertex const &v,
                                     std::unordered_set<int> &visited_v) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  Halfedge far_h(-1);

  Vertex cur_v = v;

  do {
    auto hit = mesh.halfedges(cur_v);
    auto hit_end = hit;
    do {
      if (is_contour[mesh.edge(*hit)] < 0 && !mesh.is_boundary(mesh.edge(*hit)))
        continue;

      if (visited_v.find(mesh.to_vertex(*hit).idx()) == visited_v.end()) {
        visited_v.emplace(cur_v.idx());
        cur_v = mesh.to_vertex(*hit);

        if ((vpositions[v] - vpositions[mesh.to_vertex(*hit)]).norm() >=
            EPSILON) {
          far_h = *hit;
        }
        break;
      }

    } while (++hit != hit_end);
  } while (!far_h.is_valid());

  return far_h;
}

Halfedge
find_nondegenerated_contour_halfedge_2d(Mesh const &mesh, Camera const &camera,
                                        Vertex const &v,
                                        std::unordered_set<int> &visited_v) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  Halfedge far_h(-1);
  Halfedge backup_far_h(-1);

  Vertex cur_v = v;
  do {
    bool seen_contour = false;
    auto hit = mesh.halfedges(cur_v);
    auto hit_end = hit;
    do {
      if (is_contour[mesh.edge(*hit)] < 0 && !mesh.is_boundary(mesh.edge(*hit)))
        continue;

      if (visited_v.find(mesh.to_vertex(*hit).idx()) == visited_v.end()) {
        backup_far_h = *hit;
        seen_contour = true;
        visited_v.emplace(cur_v.idx());
        cur_v = mesh.to_vertex(*hit);

        Vector2f v0_2d = project(vpositions[v], camera.viewMatrix().matrix(),
                                 camera.projectionMatrix(),
                                 Vector2i(camera.vpWidth(), camera.vpHeight()))
                             .head(2);
        Vector2f v1_2d =
            project(vpositions[mesh.to_vertex(*hit)],
                    camera.viewMatrix().matrix(), camera.projectionMatrix(),
                    Vector2i(camera.vpWidth(), camera.vpHeight()))
                .head(2);

        if ((v0_2d - v1_2d).norm() >= EPSILON) {
          far_h = *hit;
        }
        break;
      }

    } while (++hit != hit_end);
    if (!seen_contour)
      break;
  } while (!far_h.is_valid());

  if (!far_h.is_valid())
    far_h = backup_far_h;

  return far_h;
}

bool point_to_triangle(Vector3f const &p, Vector3f const &a, Vector3f const &b,
                       Vector3f const &c, Vector3f const &norm) {
  Vector3f up = a + norm;

  if (((predicates::orient3d(up.cast<double>(), b.cast<double>(),
                             a.cast<double>(), p.cast<double>()) ==
        predicates::Orientation::NEGATIVE) &&
       (predicates::orient3d(up.cast<double>(), a.cast<double>(),
                             c.cast<double>(), p.cast<double>()) ==
        predicates::Orientation::NEGATIVE)) ||
      (predicates::orient3d(up.cast<double>(), b.cast<double>(),
                            a.cast<double>(), p.cast<double>()) ==
       predicates::Orientation::COPLANAR) ||
      (predicates::orient3d(up.cast<double>(), a.cast<double>(),
                            c.cast<double>(), p.cast<double>()) ==
       predicates::Orientation::COPLANAR))
    return true;

  return false;
}

bool point_to_triangle(Vector3f const &p, Vector3f const &a, Vector3f const &b,
                       Vector3f const &c) {
  // Assume a is at the vertex position and abc is CCW
  Vector3f norm = (b - a).cross(c - a);
  return point_to_triangle(p, a, b, c, norm);
}

FacingType get_cusp_facing(Mesh &mesh, Camera const &camera,
                           Vertex const &center) {
  FacingType facing = FacingType::UNDEFINED;
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");
  auto contour_laplacian =
      mesh.get_vertex_property<Vector3f>("v:contour_laplacian");

  contess_assert_msg(contour_laplacian,
                     "get_cusp_facing: Requires contour curve laplacian.");

  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");
  contess_assert(param_loc);

  Vertex v = center;
  Vector3f laplacian = contour_laplacian[v];

  Halfedge he;
  auto hit = mesh.halfedges(v);
  auto hit_end = hit;
  do {
    if (!param_loc[*hit].is_valid())
      continue;
    he = *hit;
    break;
  } while (++hit != hit_end);
  contess_assert(he.is_valid());

  const Param_loc &v_param_loc = param_loc[he];
  Vector3f pos_res, normal_res;
  mesh.subdivision().evaluateLimit(v_param_loc, pos_res, normal_res);

  // Check which face this laplacian is pointing to
  Vector3f p = laplacian + vpositions[v];
  hit = mesh.halfedges(v);
  hit_end = hit;
  do {
    Vertex vv2 = mesh.to_vertex(*hit);
    Vertex vv3 = mesh.to_vertex(mesh.next_halfedge(*hit));
    Vector3f f_norm = mesh.compute_face_normal(mesh.face(*hit));

    double side_dot1 = (vpositions[vv2] - vpositions[v]).dot(laplacian);
    double side_dot2 = (vpositions[vv3] - vpositions[v]).dot(laplacian);

    // Trivial check for flipped surface triangles
    double dot_prod = f_norm.dot(normal_res);
    if (dot_prod < 0 || side_dot1 < 0 || side_dot2 < 0)
      continue;

    if (point_to_triangle(p, vpositions[v], vpositions[vv2], vpositions[vv3])) {
      Face f = mesh.face(*hit);
      facing = VBO[f];

      break;
    }
  } while (++hit != hit_end);

  // Use actual surface walking to determine the facing
  if (facing == FacingType::UNDEFINED) {
    Vector3f cut_p = vpositions[v];
    Vector3f path_norm = (normal_res).cross(laplacian).normalized();
    FacingType walk_facing =
        walk_on_original_facing(mesh, v, cut_p, path_norm, camera, laplacian);
    facing = walk_facing;
  }

  return facing;
}

bool compute_cusp_laplacian(Mesh &mesh) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto is_cusp = mesh.get_vertex_property<bool>("v:cusp");

  contess_assert_msg(is_contour && is_cusp, "Requires contours and cusps.");

  auto contour_laplacian =
      mesh.vertex_property<Vector3f>("v:contour_laplacian");
  if (!contour_laplacian) {
    contour_laplacian = mesh.add_vertex_property<Vector3f>(
        "v:contour_laplacian", Vector3f::Zero());
    contour_laplacian.vector().assign(contour_laplacian.vector().size(),
                                      Vector3f::Zero());
  }

  for (uint32_t i = 0; i != mesh.n_vertices(); ++i) {
    Vertex v(i);

    if (!is_cusp[v])
      continue;

    bool is_on_contour = false;
    auto hit = mesh.halfedges(v);
    auto hit_end = hit;
    do {
      Edge e = mesh.edge(*hit);

      if (is_contour[e] >= 0) {
        is_on_contour = true;
        break;
      }
    } while (++hit != hit_end);

    if (is_on_contour) {
      std::unordered_set<int> visited_v;
      visited_v.emplace(v.idx());
      if (!mesh.is_boundary(v)) {
        Halfedge he1 = find_nondegenerated_contour_halfedge(mesh, v, visited_v);
        visited_v.emplace(mesh.to_vertex(he1).idx());
        Halfedge he2 = find_nondegenerated_contour_halfedge(mesh, v, visited_v);

        Vertex v1 = mesh.to_vertex(he1);
        Vertex v2 = mesh.to_vertex(he2);

        real_t w1 = 1 / (vpositions[v1] - vpositions[v]).norm();
        real_t w2 = 1 / (vpositions[v2] - vpositions[v]).norm();

        Vector3f laplacian =
            ((w1 * vpositions[v1] + w2 * vpositions[v2]) / (w1 + w2) -
             vpositions[v]);
        real_t dot_prod = laplacian.norm();
        if (dot_prod < EPSILON)
          return false;
        laplacian.normalize();

        contour_laplacian[v] = laplacian;
      }
    }
  }

  return true;
}

bool is_potential_cusp(Mesh const &mesh, Camera const &camera,
                       Vertex const &v) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");

  bool seen_contour = false;
  auto hit = mesh.halfedges(v);
  auto hit_end = hit;
  do {
    if (!(*hit).is_valid())
      return false;

    if (is_contour[mesh.edge(*hit)] < 0)
      continue;

    seen_contour = true;
    break;
  } while (++hit != hit_end);

  if (!seen_contour)
    return false;

  std::unordered_set<int> visited_v;
  visited_v.emplace(v.idx());
  Halfedge he1 =
      find_nondegenerated_contour_halfedge_2d(mesh, camera, v, visited_v);
  visited_v.emplace(mesh.to_vertex(he1).idx());
  Halfedge he2 =
      find_nondegenerated_contour_halfedge_2d(mesh, camera, v, visited_v);
  std::vector<Halfedge> hes({he1, he2});

  real_t high_dot_prod_threshold = 0.5; // 60 deg
  std::vector<Vector2f> adj_edge_dirs;
  for (auto const &he : hes) {
    Vector2f v1_2d =
        project(vpositions[mesh.from_vertex(he)], camera.viewMatrix().matrix(),
                camera.projectionMatrix(),
                Vector2i(camera.vpWidth(), camera.vpHeight()))
            .head(2);
    Vector2f v2_2d =
        project(vpositions[mesh.to_vertex(he)], camera.viewMatrix().matrix(),
                camera.projectionMatrix(),
                Vector2i(camera.vpWidth(), camera.vpHeight()))
            .head(2);

    Vector2f v1v2 = (v2_2d - v1_2d);
    adj_edge_dirs.emplace_back(v1v2.normalized());
  }

  real_t dot_prod = adj_edge_dirs[0].dot(adj_edge_dirs[1]);
  return (dot_prod > 0 && std::abs(dot_prod) > high_dot_prod_threshold);
}

bool tag_cusp_facing(Mesh &mesh, Camera const &camera) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto is_cusp = mesh.get_vertex_property<bool>("v:cusp");

  contess_assert_msg(is_contour && is_cusp, "Requires contours and cusps.");

  auto cusp_facing = mesh.get_vertex_property<FacingType>("v:cusp_facing");
  if (!cusp_facing) {
    cusp_facing = mesh.add_vertex_property<FacingType>("v:cusp_facing",
                                                       FacingType::UNDEFINED);
  }
  cusp_facing = mesh.vertex_property<FacingType>("v:cusp_facing");
  cusp_facing.vector().assign(cusp_facing.vector().size(),
                              FacingType::UNDEFINED);

  for (uint32_t i = 0; i != mesh.n_vertices(); ++i) {
    Vertex v(i);

    bool is_potential = is_potential_cusp(mesh, camera, v);
    if (!is_cusp[v] && !is_potential)
      continue;

    if (!is_cusp[v] && is_potential) {
      logger().info("Potential cusp: {}", v);
    }

    bool is_on_contour = false;
    auto hit = mesh.halfedges(v);
    auto hit_end = hit;
    do {
      Edge e = mesh.edge(*hit);

      if (is_contour[e] >= 0) {
        is_on_contour = true;
        break;
      }
    } while (++hit != hit_end);

    if (is_on_contour) {
      std::unordered_set<int> visited_v;
      visited_v.emplace(v.idx());
      if (!mesh.is_boundary(v)) {
        cusp_facing[v] = get_cusp_facing(mesh, camera, v);

        if (cusp_facing[v] == FacingType::UNDEFINED) {
          logger().error(
              "tag_cusp_facing: Unable to detect cusp facing at vertex {}", v);
          return false;
        }
      }
    }
  }
  return true;
}

void tag_extraordinary(Mesh &mesh) {
  auto is_extraordinary = mesh.vertex_property<bool>("v:extraordinary");
  if (!is_extraordinary) {
    is_extraordinary = mesh.add_vertex_property<bool>("v:extraordinary", false);
  }
  is_extraordinary.vector().assign(is_extraordinary.vector().size(), false);
  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");

  for (uint32_t i = 0; i != mesh.n_vertices(); ++i) {
    Vertex v(i);

    std::unordered_set<int> control_patches;
    auto hit = mesh.halfedges(v);
    auto hit_end = hit;
    do {
      control_patches.emplace(param_loc[*hit].ptexIndex);
      control_patches.emplace(
          param_loc[mesh.opposite_halfedge(*hit)].ptexIndex);
    } while (++hit != hit_end);

    if (control_patches.size() != 4 && control_patches.size() > 2) {
      is_extraordinary[v] = true;
    }
  }
}

void tag_extraordinary_quad(Mesh const &mesh,
                            std::vector<Vertex> &extraordinary_v) {
  extraordinary_v.clear();

  for (uint32_t i = 0; i != mesh.n_vertices(); ++i) {
    Vertex v(i);

    int valence = 0;
    auto hit = mesh.halfedges(v);
    auto hit_end = hit;
    do {
      valence++;
    } while (++hit != hit_end);

    if (valence != 4) {
      extraordinary_v.emplace_back(v);
    }
  }
}

void tag_extraordinary_quad(Mesh &mesh) {
  auto is_extraordinary = mesh.vertex_property<bool>("v:extraordinary");
  if (!is_extraordinary) {
    is_extraordinary = mesh.add_vertex_property<bool>("v:extraordinary", false);
  }
  is_extraordinary.vector().assign(is_extraordinary.vector().size(), false);

  contess_assert_msg(mesh.is_quad_mesh(),
                     "tag_extraordinary_quad: Only takes quad mesh.");

  std::vector<Vertex> extraordinary_v;
  tag_extraordinary_quad(mesh, extraordinary_v);
  for (auto const &v : extraordinary_v) {
    is_extraordinary[v] = true;
  }
}
