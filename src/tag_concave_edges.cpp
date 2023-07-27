// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include <cmath>
#include <igl/predicates/predicates.h>

// Eigen printout formatting
#include <spdlog/fmt/ostr.h>
#include <vector>

#include "subdiv_osd.h"
#include "tag_concave_edges.h"

#include "common.h"
#include "evaluate_radial_curvature.h"
#include "inflat_path.h"
#include "logger.h"
#include "ray_cast_QI.h"
#include "tag_cusp_facing.h"

namespace predicates = igl::predicates;

bool is_geometric_concave(Mesh const &mesh, Edge const &e) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");

  surface_mesh::Surface_mesh::Halfedge h = mesh.halfedge(e, 0);

  const Vector3f &a = vpositions[mesh.from_vertex(h)];
  const Vector3f &b = vpositions[mesh.to_vertex(h)];

  // concavity test
  const Vector3f &c = vpositions[mesh.to_vertex(mesh.next_halfedge(h))];
  const Vector3f &d =
      vpositions[mesh.to_vertex(mesh.next_halfedge(mesh.opposite_halfedge(h)))];
  if (predicates::orient3d(a.cast<double>(), b.cast<double>(), c.cast<double>(),
                           d.cast<double>()) !=
      predicates::Orientation::POSITIVE) {
    // concave edge
    return true;
  }

  return false;
}

bool is_back_facing(Mesh const &mesh, Vector3f const &contour_v1,
                    Vector3f const &contour_v2, Vector3f const &contour_up,
                    Vertex const &f_v, Vertex const &b_v, Vector3f const &cam) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  const Vector3f &a = contour_v1;
  const Vector3f &b = contour_v2;
  Vector3f e = contour_up;

  Vector3f p = vpositions[f_v];
  Vector3f p_back_f = vpositions[b_v];

  // The front portion vertex can be at the degenerate position
  // (almost colinear with the contour edge)
  // In this case, test the back front face
  if (vertex_to_edge_distance(a, b, p) < FRONT_PORTION_VERTEX_EPSILON) {
    logger().info("Degenerate front vertex: {}",
                  vertex_to_edge_distance(a, b, p));
    if (!(predicates::orient3d(a.cast<double>(), b.cast<double>(),
                               e.cast<double>(), p_back_f.cast<double>()) ==
          predicates::Orientation::POSITIVE) ==
        (predicates::orient3d(a.cast<double>(), b.cast<double>(),
                              e.cast<double>(), cam.cast<double>()) ==
         predicates::Orientation::POSITIVE))
      return false;

    return true;
  }

  // If the vertex on the front portion is on the same side of the camera
  // then the edge is convex
  if ((predicates::orient3d(a.cast<double>(), b.cast<double>(),
                            e.cast<double>(), p.cast<double>()) ==
       predicates::Orientation::POSITIVE) ==
      (predicates::orient3d(a.cast<double>(), b.cast<double>(),
                            e.cast<double>(), cam.cast<double>()) ==
       predicates::Orientation::POSITIVE))
    return false;

  return true;
}

bool is_back_facing_debug(Mesh &mesh, Edge const &edge, Camera const &camera) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");
  auto contour_front = mesh.vertex_property<bool>("v:contour_front");

  surface_mesh::Surface_mesh::Halfedge h = mesh.halfedge(edge, 0);

  // Build a triangle which has the contour edge as a side and one side
  // parallel to the image plane (necessary?)
  const Vector3f &a = vpositions[mesh.from_vertex(h)];
  const Vector3f &b = vpositions[mesh.to_vertex(h)];

  Vector3f ab = b - a;
  Vector3f ac = camera.position() - a;
  Vector3f on_image_plane = ac.cross(ab);
  on_image_plane.normalize();
  Vector3f e = a + on_image_plane;

  // Find the adjacent vertex that is on the front portion
  Face f1 = mesh.face(edge, 0);

  Vertex f1_v = mesh.to_vertex(mesh.next_halfedge(h));
  Vertex f2_v = mesh.to_vertex(mesh.next_halfedge(mesh.opposite_halfedge(h)));

  // Debug: label the vertex on the front portion
  Vertex p_index = (VBO[f1] == FacingType::FRONT) ? f1_v : f2_v;
  Vertex p_back_index = (VBO[f1] == FacingType::BACK) ? f1_v : f2_v;
  if (mesh.is_contour_edge(edge))
    contour_front[p_index] = true;

  return is_back_facing(mesh, a, b, e, p_index, p_back_index,
                        camera.position());
}

bool is_back_facing(Mesh const &mesh, Edge const &edge, Camera const &camera) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");

  surface_mesh::Surface_mesh::Halfedge h = mesh.halfedge(edge, 0);

  // Build a triangle which has the contour edge as a side and one side
  // parallel to the image plane (necessary?)
  const Vector3f &a = vpositions[mesh.from_vertex(h)];
  const Vector3f &b = vpositions[mesh.to_vertex(h)];

  Vector3f ab = b - a;
  Vector3f ac = camera.position() - a;
  Vector3f on_image_plane = ac.cross(ab);
  on_image_plane.normalize();
  Vector3f e = a + on_image_plane;

  // Find the adjacent vertex that is on the front portion
  Face f1 = mesh.face(edge, 0);

  Vertex f1_v = mesh.to_vertex(mesh.next_halfedge(h));
  Vertex f2_v = mesh.to_vertex(mesh.next_halfedge(mesh.opposite_halfedge(h)));

  Vertex p_index = (VBO[f1] == FacingType::FRONT) ? f1_v : f2_v;
  Vertex p_back_index = (VBO[f1] == FacingType::BACK) ? f1_v : f2_v;

  return is_back_facing(mesh, a, b, e, p_index, p_back_index,
                        camera.position());
}

bool is_back_facing_interpolated(Mesh &mesh, Vector3f const &contour_v1,
                                 Vector3f const &contour_v2,
                                 size_t interpolated_f_index,
                                 Camera const &camera) {
  auto contour_front = mesh.vertex_property<bool>("v:contour_front");

  // Build a triangle which has the contour edge as a side and one side
  // parallel to the image plane (necessary?)
  const Vector3f &a = contour_v1;
  const Vector3f &b = contour_v2;

  Vector3f ab = b - a;
  Vector3f ac = camera.position() - a;
  Vector3f on_image_plane = ac.cross(ab);
  on_image_plane.normalize();
  Vector3f e = a + on_image_plane;

  // 1. Find the front and back vertices on the interpolated contour face
  Face f = Face(interpolated_f_index);

  // Find all vertices of the face
  std::vector<Vertex> f_vertices;
  {
    surface_mesh::Surface_mesh::Halfedge_around_face_circulator hit, hit_end;
    hit = mesh.halfedges(f);
    hit_end = hit;
    do {
      f_vertices.emplace_back(mesh.to_vertex(*hit));
    } while (++hit != hit_end);
  }

  // Determine the front and back vertices
  std::map<FacingType, std::vector<Vertex>> facing_vertices;
  facing_vertices[FacingType::FRONT] = std::vector<Vertex>();
  facing_vertices[FacingType::BACK] = std::vector<Vertex>();
  // This function works under the assumption that a interpolated contour face
  // is always adjacent to at least one non-contour (F/B-facing) face.
  auto get_vertex_facing = [](Mesh const &mesh, Vertex const &v) -> FacingType {
    FacingType v_facing = FacingType::UNDEFINED;
    auto VBO = mesh.get_face_property<FacingType>("f:VBO");

    surface_mesh::Surface_mesh::Halfedge_around_vertex_circulator hit, hit_end;
    hit = mesh.halfedges(v);
    hit_end = hit;
    do {
      Face fit = mesh.face(*hit);

      if (VBO[fit] == FacingType::FRONT || VBO[fit] == FacingType::BACK) {
        v_facing = VBO[fit];
      }
    } while (++hit != hit_end);

    return v_facing;
  };
  for (auto v : f_vertices) {
    FacingType v_f = get_vertex_facing(mesh, v);
    assert(v_f != FacingType::UNDEFINED);

    facing_vertices[v_f].emplace_back(v);
  }

  // 2. Keep the further vertex when there are two vertices on one side
  // of the contour
  auto find_further_vertex = [](Mesh const &mesh, Vector3f const &a,
                                Vector3f const &b,
                                std::vector<Vertex> const &vertices) -> Vertex {
    assert(!vertices.empty());
    auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
    Vertex far_v = vertices.front();

    if (vertices.size() == 1)
      return far_v;

    real_t dist1 = vertex_to_edge_distance(a, b, vpositions[vertices[0]]);
    real_t dist2 = vertex_to_edge_distance(a, b, vpositions[vertices[1]]);
    far_v = (dist1 < dist2) ? vertices[1] : vertices[0];

    return far_v;
  };
  facing_vertices[FacingType::FRONT] = std::vector<Vertex>(
      {find_further_vertex(mesh, a, b, facing_vertices[FacingType::FRONT])});
  facing_vertices[FacingType::BACK] = std::vector<Vertex>(
      {find_further_vertex(mesh, a, b, facing_vertices[FacingType::BACK])});
  // Debug: label the vertex on the front portion
  contour_front[facing_vertices[FacingType::FRONT].front()] = true;

  return is_back_facing(
      mesh, a, b, e, facing_vertices[FacingType::FRONT].front(),
      facing_vertices[FacingType::BACK].front(), camera.position());
}

/* =========================================================== */

void tag_concave_edges_mesh(Mesh &mesh, ContourMode mode) {
  auto concave = mesh.get_edge_property<bool>("e:concave");
  if (!concave) {
    concave = mesh.add_edge_property<bool>("e:concave", false);
  }
  concave = mesh.edge_property<bool>("e:concave");
  concave.vector().assign(concave.vector().size(), false);

  // Label convex/concave on face for the interpolated contour
  auto concave_f = mesh.face_property<bool>("f:concave");
  concave_f.vector().assign(concave_f.vector().size(), false);

  parallel_for(
      tbb::blocked_range<uint32_t>(0u, (uint32_t)mesh.n_edges(), GRAIN_SIZE),
      [&](const tbb::blocked_range<uint32_t> &range) {
        for (uint32_t i = range.begin(); i != range.end(); ++i) {
          Edge e = Edge(i);
          concave[e] = is_geometric_concave(mesh, e);
        }
      });
}

void tag_concave_edges_side(Mesh &mesh, Camera const &camera,
                            ContourMode mode) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto concave = mesh.get_edge_property<bool>("e:concave");
  if (!concave) {
    concave = mesh.add_edge_property<bool>("e:concave", false);
  }
  concave = mesh.edge_property<bool>("e:concave");
  concave.vector().assign(concave.vector().size(), false);

  // Label convex/concave on face for the interpolated contour
  auto concave_f = mesh.face_property<bool>("f:concave");
  concave_f.vector().assign(concave_f.vector().size(), false);

  // Set debug flag for the vertices belong to the front portion
  auto contour_front = mesh.vertex_property<bool>("v:contour_front");
  contour_front.vector().assign(contour_front.vector().size(), false);

  // Run on all edges?
  // Distinguish the contour generator edges?
  if (mode != ContourMode::INTERPOLATED_CONTOUR) {
    parallel_for(
        tbb::blocked_range<uint32_t>(0u, (uint32_t)mesh.n_edges(), GRAIN_SIZE),
        [&](const tbb::blocked_range<uint32_t> &range) {
          for (uint32_t i = range.begin(); i != range.end(); ++i) {
            Edge e = Edge(i);

            if (is_contour[e] >= 0)
              concave[e] = is_back_facing_debug(mesh, e, camera);
          }
        });
  } else {
    parallel_for(
        tbb::blocked_range<uint32_t>(
            0u, (uint32_t)mesh.get_const_interpolated_contour_faces().size(),
            GRAIN_SIZE),
        [&](const tbb::blocked_range<uint32_t> &range) {
          for (uint32_t i = range.begin(); i != range.end(); ++i) {
            Face f(mesh.get_interpolated_contour_faces()[i]);
            concave_f[f] = is_back_facing_interpolated(
                mesh, mesh.get_interpolated_contours().col(2 * i),
                mesh.get_interpolated_contours().col(2 * i + 1), f.idx(),
                camera);
          }
        });
  }
}

void tag_concave_edges_analytic(Mesh &mesh, Camera const &camera) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto concave_kr = mesh.edge_property<bool>("e:concave_kr");
  if (!concave_kr)
    concave_kr = mesh.add_edge_property<bool>("e:concave_kr");
  concave_kr.vector().assign(concave_kr.vector().size(), false);

  bool to_write_concave = true;
  auto concave = mesh.get_edge_property<bool>("e:concave");
  auto concave_consistent = mesh.get_edge_property<bool>("e:concave_consist");
  if (concave && !concave_consistent) {
    concave_consistent = mesh.add_edge_property<bool>("e:concave_consist");
    concave_consistent = mesh.edge_property<bool>("e:concave_consist");
    concave_consistent.vector().assign(concave_consistent.vector().size(),
                                       true);
    to_write_concave = false;
  } else if (!concave) {
    concave = mesh.add_edge_property<bool>("e:concave", false);
  }

  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");

  auto compute_kr = [&](uint32_t i) {
    Edge e = Edge(i);

    if (is_contour[e] < 0)
      return;

    // Evaluate k_r at mid point
    Vertex v1 = mesh.vertex(e, 0);
    Vertex v2 = mesh.vertex(e, 1);

    std::unordered_map<int, std::pair<Halfedge, Halfedge>> matched_param;
    auto hit = mesh.halfedges(v1), hit_end = hit;
    do {
      if (!matched_param.count(param_loc[*hit].ptexIndex))
        matched_param[param_loc[*hit].ptexIndex].first = *hit;
    } while (++hit != hit_end);
    hit = mesh.halfedges(v2), hit_end = hit;
    do {
      if (matched_param.count(param_loc[*hit].ptexIndex) &&
          !matched_param[param_loc[*hit].ptexIndex].second.is_valid())
        matched_param[param_loc[*hit].ptexIndex].second = *hit;
    } while (++hit != hit_end);

    Halfedge he1, he2;

    for (auto const &param_p : matched_param) {
      if (param_p.first >= 0 && param_p.second.second.is_valid()) {
        he1 = param_p.second.first;
        he2 = param_p.second.second;
      }
    }

    contess_assert_msg(
        he1.is_valid() && he2.is_valid() &&
            (param_loc[he1].ptexIndex == param_loc[he2].ptexIndex),
        "tag_concave_edges_analytic: Parameter unmatched on edge.");

    Param_loc mid_param;
    mid_param.ptexIndex = param_loc[he1].ptexIndex;
    mid_param.uv = 0.5 * (param_loc[he1].uv + param_loc[he2].uv);
    real_t kr = evaluate_radial_curvature_analytic(
        mesh, mesh.const_subdivision(), camera, Vertex(), mid_param);

    // To handle the case of having a local maxima within the edge
    // Aka, P-N-P, or N-P-N. This may happen when contour root finding
    // fails to converge
    real_t kr_1 = evaluate_radial_curvature_analytic(
        mesh, mesh.const_subdivision(), camera, Vertex(), param_loc[he1]);
    real_t kr_2 = evaluate_radial_curvature_analytic(
        mesh, mesh.const_subdivision(), camera, Vertex(), param_loc[he2]);
    if ((kr_1 > 0) == (kr_2 > 0) && (kr_1 > 0) != (kr > 0))
      kr = 0.5 * (kr_1 + kr_2);

    concave_kr[e] = (kr < 0);
    if (!to_write_concave && concave_consistent)
      concave_consistent[e] = concave_kr[e] == concave[e];
    else
      concave[e] = concave_kr[e];
  };

  for (uint32_t i = 0; i != mesh.n_edges(); ++i) {
    compute_kr(i);
  }
}

/* ======================================================================================
 */
void get_curtain_folds_interpolated(Mesh const &mesh,
                                    std::vector<Vector3f> &curtain_folds) {
  auto concave_f = mesh.get_face_property<bool>("f:concave");
  assert(concave_f);
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");

  curtain_folds.clear();

  parallel_for(
      tbb::blocked_range<uint32_t>(
          0u, (uint32_t)mesh.get_const_interpolated_contour_faces().size(),
          GRAIN_SIZE),
      [&](const tbb::blocked_range<uint32_t> &range) {
        for (uint32_t i = range.begin(); i != range.end(); ++i) {
          Face f(mesh.get_const_interpolated_contour_faces()[i]);

          // Find convexity changing side
          surface_mesh::Surface_mesh::Halfedge_around_face_circulator hit,
              hit_end;
          hit = mesh.halfedges(f);
          hit_end = hit;

          do {
            Face f_n = mesh.face(mesh.opposite_halfedge(*hit));

            if (VBO[f_n] == FacingType::UNDEFINED &&
                concave_f[f_n] != concave_f[f]) {
              // Find the vertex corresponding to this change
              Vertex v1 = mesh.to_vertex(*hit);
              Vertex v2 = mesh.from_vertex(*hit);

              Vector3f c_v1 = mesh.get_const_interpolated_contours().col(2 * i);
              Vector3f c_v2 =
                  mesh.get_const_interpolated_contours().col(2 * i + 1);

              real_t dist1 =
                  vertex_to_edge_distance(vpositions[v1], vpositions[v2], c_v1);
              real_t dist2 =
                  vertex_to_edge_distance(vpositions[v1], vpositions[v2], c_v2);
              Vector3f curtain_v = (dist1 < dist2) ? c_v1 : c_v2;
              curtain_folds.emplace_back(curtain_v);
              break;
            }
          } while (++hit != hit_end);
        }
      });
}

/* ======================================================================================
 */
void tag_cusp_given_concave_edges(Mesh &mesh, Camera const &camera) {
  auto concave = mesh.get_edge_property<bool>("e:concave");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");

  contess_assert_msg(concave, "Requires convex/concave labels on edges.");

  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto is_cusp = mesh.vertex_property<bool>("v:cusp");
  if (!is_cusp) {
    is_cusp = mesh.add_vertex_property<bool>("v:cusp", false);
  }
  is_cusp.vector().assign(is_cusp.vector().size(), false);

  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");
  real_t high_k_r_threshold = 100;
  real_t high_dot_prod_threshold = 0.8;
  if (mesh.subdivision().backend_type() == Subdiv::Backend::LACEWELL) {
    for (uint32_t i = 0; i != mesh.n_vertices(); ++i) {
      Vertex v(i);

      std::vector<Vector2f> adj_edge_dirs;
      adj_edge_dirs.resize(2, Vector2f::Zero());
      bool seen_convex = false;
      bool seen_concave = false;
      auto hit = mesh.halfedges(v);
      auto hit_end = hit;
      Param_loc param;
      do {
        Edge e = mesh.edge(*hit);
        Vertex v1 = mesh.from_vertex(*hit);
        Vertex v2 = mesh.to_vertex(*hit);

        // Avoid extremely short edges
        if (is_contour[e] >= 0 &&
            (vpositions[v1] - vpositions[v2]).norm() > EPSILON) {
          seen_convex |= !concave[e];
          seen_concave |= concave[e];

          Vector2f v1_2d =
              project(vpositions[v1], camera.viewMatrix().matrix(),
                      camera.projectionMatrix(),
                      Vector2i(camera.vpWidth(), camera.vpHeight()))
                  .head(2);
          Vector2f v2_2d =
              project(vpositions[v2], camera.viewMatrix().matrix(),
                      camera.projectionMatrix(),
                      Vector2i(camera.vpWidth(), camera.vpHeight()))
                  .head(2);

          if (!concave[e] && (v2_2d - v1_2d).norm() > 1e-4) {
            adj_edge_dirs[0] = (v2_2d - v1_2d).normalized();
          } else if (concave[e] && (v2_2d - v1_2d).norm() > 1e-4) {
            adj_edge_dirs[1] = (v2_2d - v1_2d).normalized();
          }
        }

        if (param_loc[*hit].is_valid())
          param = param_loc[*hit];
      } while (++hit != hit_end);

      if (seen_convex && seen_concave) {
        // Kr filtering
        if (param.is_valid()) {
          real_t kr = evaluate_radial_curvature_analytic(
              mesh, mesh.const_subdivision(), camera, Vertex(), param);
          real_t dot_prod = adj_edge_dirs[0].dot(adj_edge_dirs[1]);
          if (std::abs(kr) > high_k_r_threshold ||
              (dot_prod < 0 && std::abs(dot_prod) > high_dot_prod_threshold)) {
            logger().warn("tag_cusp_given_concave_edges: Skipping high "
                          "curvature cusp {} - {}, {}",
                          v, kr, dot_prod);
            continue;
          }
        }

        is_cusp[v] = true;
      }
    }
  } else {
    contess_assert(0);
  }
}

void tag_cusp_sharp(Mesh &mesh, Camera const &camera) {
  auto is_cusp = mesh.vertex_property<bool>("v:cusp");
  for (uint32_t i = 0; i != mesh.n_vertices(); ++i) {
    Vertex v(i);

    if (is_cusp[v])
      continue;

    is_cusp[v] = is_potential_cusp(mesh, camera, v);
  }
}

int walk_on_original_cusp(Mesh const &mesh, int id, Vertex src_orig,
                          Vector3f const &cut_p, Vector3f const &cut_norm) {
  int path_count = 0;
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  contess_assert_msg(patchID,
                     "walk_on_original_cusp: Extracted patches is needed.");

  std::unordered_set<int> visited_face_indices;
  std::queue<Face> face_frontier;

  // Initialize with all neighboring faces
  {
    auto hit = mesh.halfedges(src_orig), hit_end = hit;
    do {
      // Check adjacent faces
      auto adj_f_h = mesh.opposite_halfedge(*hit);
      Face adj_f = mesh.face(adj_f_h);

      if (adj_f.is_valid() && patchID[adj_f] != id)
        continue;

      // In case we are at the boundary of the model
      if (adj_f.is_valid()) {
        face_frontier.push(adj_f);
        visited_face_indices.emplace(adj_f.idx());
      }
    } while (++hit != hit_end);
  }

  if (face_frontier.empty())
    return path_count;

  auto is_destination = [&](Face f) -> bool {
    bool is_contour_f = false;
    auto hit = mesh.halfedges(f), hit_end = hit;
    do {
      Edge e = mesh.edge(*hit);
      if ((mesh.vertex(e, 0) != src_orig && mesh.vertex(e, 1) != src_orig) &&
          (is_contour[e] >= 0 || mesh.is_boundary(e))) {
        Vector3f intersection;
        bool edge_intersected = planeIntersectEdge(
            cut_p, cut_norm, vpositions[mesh.from_vertex(*hit)],
            vpositions[mesh.to_vertex(*hit)], intersection);
        if (edge_intersected) {
          is_contour_f = true;
          break;
        }
      }
    } while (++hit != hit_end);

    if (!is_contour_f)
      return false;

    return true;
  };

  // Walk
  int tolerant_steps = 5;
  std::unordered_map<int, int> prev_idx;
  std::unordered_map<int, Vector3f> f_intersection;
  logger().info("Start from: {}", src_orig);
  do {
    Face cur_f = face_frontier.front();
    face_frontier.pop();

    // Near exception
    auto get_step_count = [&](Face const &end_f) -> int {
      int count = 0;
      int count_f = end_f.idx();
      while (prev_idx.count(count_f) != 0) {
        count_f = prev_idx[count_f];
        count++;
      }
      return count;
    };
    bool is_first_step = (get_step_count(cur_f) < tolerant_steps);

    // Is this face intersecting with the cut plane?
    Vertex vertices[3];
    mesh.verticesOfFace(cur_f, vertices);
    std::vector<Vector3f> tri_p({vpositions[vertices[0]],
                                 vpositions[vertices[1]],
                                 vpositions[vertices[2]]});

    // March along the intersected triangles
    logger().info("\tTest: {}", cur_f);
    if (planeIntersectTri(cut_p, cut_norm, tri_p) || is_first_step) {
      logger().info("\t\tMarch on: {}", cur_f.idx());

      // If this face the destination?
      if (is_destination(cur_f)) {
        path_count++;
        logger().info("\t\tReached {} <- {}", cur_f, prev_idx[cur_f.idx()]);
      } else {
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
          if (!edge_intersected && !is_first_step)
            continue;

          f_intersection[cur_f.idx()] = intersection;

          // We can't cross the contour edge
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
      }
    }

    visited_face_indices.insert(cur_f.idx());
  } while (!face_frontier.empty());

  return path_count;
}

void tag_cusp_given_concave_edges_walking(Mesh &mesh, Camera const &camera) {
  auto concave = mesh.get_edge_property<bool>("e:concave");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");

  contess_assert_msg(concave, "Requires convex/concave labels on edges.");

  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto is_cusp = mesh.vertex_property<bool>("v:cusp");
  if (!is_cusp) {
    is_cusp = mesh.add_vertex_property<bool>("v:cusp");
  }
  is_cusp.vector().assign(is_cusp.vector().size(), false);

  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");

  std::vector<Vertex> candidate_cusps;
  for (uint32_t i = 0; i != mesh.n_vertices(); ++i) {
    Vertex v(i);

    bool seen_convex = false;
    bool seen_concave = false;
    auto hit = mesh.halfedges(v);
    auto hit_end = hit;
    Param_loc param;
    do {
      Edge e = mesh.edge(*hit);
      Vertex v1 = mesh.vertex(e, 0);
      Vertex v2 = mesh.vertex(e, 1);

      // Avoid extremely short edges
      if (is_contour[e] >= 0 &&
          (vpositions[v1] - vpositions[v2]).norm() > EPSILON) {
        seen_convex |= !concave[e];
        seen_concave |= concave[e];
      }

      if (param_loc[*hit].is_valid())
        param = param_loc[*hit];
    } while (++hit != hit_end);

    if (seen_convex && seen_concave) {
      candidate_cusps.emplace_back(v);

      // Temporarily set them to be cusps
      is_cusp[v] = true;
    }
  }

  compute_cusp_laplacian(mesh);
  tag_cusp_facing(mesh, camera);

  // Reset cusps
  is_cusp.vector().assign(is_cusp.vector().size(), false);
  auto patchID = mesh.get_face_property<int>("f:patchID");
  auto cusp_facing = mesh.get_vertex_property<FacingType>("v:cusp_facing");
  contess_assert(cusp_facing);
  contess_assert(patchID);

  // Walk to determine if it's a cusp
  for (auto const &v : candidate_cusps) {
    int pidx = -1;

    // Find patch to walk on
    FacingType facing = cusp_facing[v];
    auto hit = mesh.halfedges(v);
    auto hit_end = hit;
    do {
      if (mesh.face(*hit).is_valid() && VBO[mesh.face(*hit)] != facing) {
        pidx = patchID[mesh.face(*hit)];
        break;
      }
    } while (++hit != hit_end);

    // Find walking direction
    std::unordered_set<int> visited_v;
    visited_v.emplace(v.idx());
    Halfedge he1 = find_nondegenerated_contour_halfedge(mesh, v, visited_v);
    visited_v.emplace(mesh.to_vertex(he1).idx());
    Halfedge he2 = find_nondegenerated_contour_halfedge(mesh, v, visited_v);

    auto normalized_2d = [&](Vector3f const &v1, Vector3f const &v2,
                             Vector3f &dir) {
      dir = v2 - v1;
      Vector2i viewport = Vector2i(camera.vpWidth(), camera.vpHeight());
      Vector2f v2_2d = project(v2, camera.viewMatrix().matrix(),
                               camera.projectionMatrix(), viewport)
                           .head<2>();
      Vector2f v1_2d = project(v1, camera.viewMatrix().matrix(),
                               camera.projectionMatrix(), viewport)
                           .head<2>();
      real_t mag_2d = (v2_2d - v1_2d).norm();
      dir /= mag_2d;
    };
    Vertex v1 = mesh.to_vertex(he1);
    Vertex v2 = mesh.to_vertex(he2);

    Vector3f v_he_dir1, v_he_dir2;
    normalized_2d(vpositions[v], vpositions[v1], v_he_dir1);
    normalized_2d(vpositions[v], vpositions[v2], v_he_dir2);

    Vector3f path_dir = v_he_dir1 + v_he_dir2;
    Vector3f cut_p = camera.position();
    Vector3f path_norm =
        (camera.position() - vpositions[v]).cross(path_dir).normalized();
    int path_count = walk_on_original_cusp(mesh, pidx, v, cut_p, path_norm);
    if (path_count <= 2) {
      logger().info("=> Not a cusp: {}", v);
    }
    if (path_count > 2)
      is_cusp[v] = true;
  }
}
