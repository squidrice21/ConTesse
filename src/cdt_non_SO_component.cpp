// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "cdt_non_SO_component.h"
#include "common.h"
#include <igl/triangle/cdt.h>
#include <igl/triangle/triangulate.h>
#include <limits>
#include <unordered_map>
#include <vector>

// Eigen printout formatting
#include <spdlog/fmt/ostr.h>

void handle_degenerated_component_facing(Mesh const &mesh,
                                         Mesh const &flat_patch,
                                         Mesh &inflated_patch) {
  // Determine the target facing
  FacingType target_facing;
  auto flat_patchID = flat_patch.get_face_property<int>("f:patchID");
  contess_assert_msg(flat_patchID,
                     "cdt_non_SO_component: Missing patch index.");
  int patch_id = flat_patchID[*flat_patch.faces_begin()];
  auto flat_patch_component =
      flat_patch.get_face_property<int>("f:patch_component");
  contess_assert_msg(flat_patch_component,
                     "cdt_non_SO_component: Correspondence is missing.");
  int input_comp_idx = flat_patch_component[Face(0)];
  auto flat_vpositions = flat_patch.get_vertex_property<Vector3f>("v:point");

  target_facing = mesh.get_patch_facing(patch_id);

  // Either copy or flip the input patch depending on the facing
  auto inflated_orig_idx =
      inflated_patch.vertex_property<int>("v:orig_idx", -1);
  inflated_orig_idx.vector().assign(inflated_orig_idx.vector().size(), -1);
  auto inflated_orig_f_idx =
      inflated_patch.vertex_property<int>("v:orig_f_idx", -1);
  inflated_orig_f_idx.vector().assign(inflated_orig_f_idx.vector().size(), -1);

  auto flat_orig_idx = flat_patch.get_vertex_property<int>("v:orig_idx");
  contess_assert_msg(
      flat_orig_idx,
      "handle_degenerated_component_facing: Correspondence is missing.");

  auto flat_comp_idx = flat_patch.get_vertex_property<int>("v:comp_idx");
  auto inflated_comp_idx =
      inflated_patch.vertex_property<int>("v:comp_idx", -1);
  if (!flat_comp_idx)
    inflated_comp_idx.vector().assign(inflated_comp_idx.vector().size(), -1);

  auto inflated_patchID = inflated_patch.face_property<int>("f:patchID", -1);
  auto patch_component =
      inflated_patch.face_property<int>("f:patch_component", -1);

  std::unordered_map<int, int> orig2new;
  for (long i = 0; i < flat_patch.n_vertices(); i++) {
    Vertex orig_v(i);
    Vector3f v3d = flat_vpositions[orig_v];
    Vertex v = inflated_patch.add_vertex(v3d);

    int orig_idx, comp_idx = -1;
    orig_idx = flat_orig_idx[orig_v];
    if (flat_comp_idx)
      comp_idx = flat_comp_idx[orig_v];
    // It's a boundary point
    inflated_orig_idx[v] = orig_idx;
    inflated_comp_idx[v] = comp_idx;
    orig2new[i] = v.idx();
  }
  // The result mesh should have the same vertex order as the input one
  for (long i = 0; i < flat_patch.n_faces(); i++) {
    Face f(i);
    Vertex vv[3];
    flat_patch.verticesOfFace(f, vv);

    std::vector<Vertex> tri({Vertex(orig2new[vv[0].idx()]),
                             Vertex(orig2new[vv[1].idx()]),
                             Vertex(orig2new[vv[2].idx()])});
    // Flip the surface if the target patch is back facing
    if (target_facing == FacingType::BACK)
      tri = std::vector<Vertex>({Vertex(orig2new[vv[0].idx()]),
                                 Vertex(orig2new[vv[2].idx()]),
                                 Vertex(orig2new[vv[1].idx()])});
    Face inflated_f = inflated_patch.add_face(tri);
    inflated_patchID[inflated_f] = patch_id;
    patch_component[inflated_f] = input_comp_idx;
  }

  inflated_patch.init();

  auto out_boundary_type = inflated_patch.edge_property<BoundaryType>(
      "e:boundary_type", BoundaryType::NONE);
  for (size_t i = 0; i < inflated_patch.n_edges(); i++) {
    Edge e(i);
    if (inflated_patch.is_boundary(e))
      out_boundary_type[e] = BoundaryType::PATCH;
    else
      out_boundary_type[e] = BoundaryType::NONE;
  }
}

void cdt_non_SO_component(
    Mesh const &mesh, Mesh const &flat_patch, Camera const &camera,
    std::vector<std::pair<Vector3f, int>> const &interior_points,
    Mesh &inflated_patch, bool to_flip) {
  if (flat_patch.n_faces() == 0)
    return;

  // Determine the target facing
  FacingType target_facing;
  auto flat_patchID = flat_patch.get_face_property<int>("f:patchID");
  contess_assert_msg(flat_patchID,
                     "cdt_non_SO_component: Missing patch index.");
  int patch_id = flat_patchID[Face(0)];
  target_facing = mesh.get_patch_facing(flat_patchID[Face(0)]);

  // Convert to eigen
  MatrixXf V;
  MatrixXi E;

  // First create the PSLG
  // https://www.cs.cmu.edu/~quake/triangle.defs.html
  // https://www.cs.cmu.edu/~quake/triangle.switch.html
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto flat_vpositions = flat_patch.get_vertex_property<Vector3f>("v:point");

  // Note that we assume flat_patch is already lifted along edges if it's a
  // component of a flat path

  auto boundary_type =
      flat_patch.get_edge_property<BoundaryType>("e:boundary_type");
  contess_assert_msg(boundary_type,
                     "cdt_non_SO_component: Boundary types are not labeled.");
  auto flat_orig_idx = flat_patch.get_vertex_property<int>("v:orig_idx");
  contess_assert_msg(flat_orig_idx,
                     "cdt_non_SO_component: Correspondence is missing.");
  auto flat_orig_f_idx = flat_patch.get_vertex_property<int>("v:orig_f_idx");
  contess_assert_msg(flat_orig_f_idx,
                     "cdt_non_SO_component: Correspondence is missing.");
  auto flat_comp_idx = flat_patch.get_vertex_property<int>("v:comp_idx");
  auto flat_patch_component =
      flat_patch.get_face_property<int>("f:patch_component");
  contess_assert_msg(flat_patch_component,
                     "cdt_non_SO_component: Correspondence is missing.");
  int comp_idx = flat_patch_component[Face(0)];
  Vector2i viewport({camera.vpWidth(), camera.vpHeight()});

  // 2D dedup
  auto dedup_2d =
      [](Camera const &camera,
         std::vector<std::pair<Vector3f, int>> const &points,
         std::vector<std::pair<Vector3f, int>> const &extra_check_points,
         std::vector<std::pair<Vector3f, int>> &points_dedup) {
        points_dedup.clear();
        points_dedup.reserve(points.size());
        Vector2i viewport({camera.vpWidth(), camera.vpHeight()});
        std::vector<Vector2f> seen_points_projected;

        // Check extra points (for the interior points, we also check against
        // the boundary points)
        for (auto const &pos3D : extra_check_points) {
          Vector2f pos2D = project(pos3D.first, camera.viewMatrix().matrix(),
                                   camera.projectionMatrix(), viewport)
                               .head<2>();
          seen_points_projected.emplace_back(pos2D);
        }

        for (auto const &pos3D : points) {
          Vector2f pos2D = project(pos3D.first, camera.viewMatrix().matrix(),
                                   camera.projectionMatrix(), viewport)
                               .head<2>();
          if (std::find(seen_points_projected.begin(),
                        seen_points_projected.end(),
                        pos2D) == seen_points_projected.end()) {
            seen_points_projected.emplace_back(pos2D);
            points_dedup.emplace_back(pos3D);
          }
        }
      };

  // Assume the polygon vertices are ordered
  std::vector<std::pair<Vector3f, int>> flat_patch_points;
  flat_patch_points.reserve(flat_patch.n_vertices());
  for (size_t i = 0; i < flat_patch.n_vertices(); i++) {
    Vertex v_orig(flat_orig_idx[Vertex(i)]);
    int f_idx;

    // Project
    Vector3f pos3D;
    if (v_orig.is_valid()) {
      pos3D = vpositions[v_orig];
    } else {
      pos3D = flat_vpositions[Vertex(i)];
    }
    f_idx = flat_orig_f_idx[Vertex(i)];

    flat_patch_points.emplace_back(pos3D, f_idx);
  }
  std::vector<std::pair<Vector3f, int>> flat_patch_points_dedup;
  dedup_2d(camera, flat_patch_points, std::vector<std::pair<Vector3f, int>>(),
           flat_patch_points_dedup);

  std::vector<std::pair<Vector3f, int>> interior_points_dedup;
  dedup_2d(camera, interior_points, flat_patch_points_dedup,
           interior_points_dedup);

  auto inflated_orig_idx =
      inflated_patch.vertex_property<int>("v:orig_idx", -1);
  inflated_orig_idx.vector().assign(inflated_orig_idx.vector().size(), -1);
  auto inflated_orig_f_idx =
      inflated_patch.vertex_property<int>("v:orig_f_idx", -1);
  inflated_orig_f_idx.vector().assign(inflated_orig_f_idx.vector().size(), -1);

  // This is the indices used to stitch different component from the same patch
  auto inflated_comp_idx =
      inflated_patch.vertex_property<int>("v:comp_idx", -1);
  if (!flat_comp_idx)
    inflated_comp_idx.vector().assign(inflated_comp_idx.vector().size(), -1);

  ///////////////////////

  V.resize(flat_patch_points_dedup.size() + interior_points_dedup.size(), 2);
  std::vector<int32_t> E_idx;

  std::vector<std::pair<Vector2f, int>> boundary_points_projected;
  boundary_points_projected.reserve(flat_patch_points_dedup.size());
  for (size_t i = 0; i < flat_patch_points_dedup.size(); i++) {
    // Project
    auto pos3D = flat_patch_points_dedup[i];
    Vector2f pos2D = project(pos3D.first, camera.viewMatrix().matrix(),
                             camera.projectionMatrix(), viewport)
                         .head<2>();

    V.row(i) = pos2D.transpose();
    boundary_points_projected.emplace_back(std::make_pair(pos2D, pos3D.second));
  }

  for (size_t i = 0; i < flat_patch_points_dedup.size(); i++) {
    E_idx.emplace_back(i);
    E_idx.emplace_back((i + 1) % flat_patch_points_dedup.size());
  }
  E = Eigen::Map<MatrixXi>(E_idx.data(), 2, E_idx.size() / 2).transpose();

  // Then fill the interior samples
  std::vector<std::pair<Vector2f, int>> interior_points_projected;
  interior_points_projected.reserve(interior_points_dedup.size());
  for (size_t i = 0; i < interior_points_dedup.size(); i++) {
    // Project to 2D
    Vector2f pos2D =
        project(interior_points_dedup[i].first, camera.viewMatrix().matrix(),
                camera.projectionMatrix(), viewport)
            .head<2>();
    V.row(flat_patch_points_dedup.size() + i) = pos2D.transpose();
    interior_points_projected.emplace_back(
        std::make_pair(pos2D, interior_points_dedup[i].second));
  }

  // Call CDT, p flag for PSLG
  MatrixXf V2, H;
  MatrixXi F2;
  MatrixXi E2;
  VectorXi J;

  // H is empty since the patch component doesn't have hole
  // igl::triangle::triangulate(V, E, H, "pVV", V2, F2); // Verbose
  igl::triangle::triangulate(V, E, H, "pQ", V2, F2);

  // Read out
  auto find_3d_point = [&](Vector2f const &v2d, std::pair<Vector3f, int> &v3d,
                           int &orig_idx, int &comp_idx) {
    orig_idx = -1;
    comp_idx = -1;

    int idx = -1;
    auto itr = std::find_if(
        boundary_points_projected.begin(), boundary_points_projected.end(),
        [&v2d](const std::pair<Vector2f, int> &x) { return x.first == v2d; });
    if (itr != boundary_points_projected.end()) {
      idx = itr - boundary_points_projected.begin();
      orig_idx = flat_orig_idx[Vertex(idx)];
      if (flat_comp_idx) {
        comp_idx = flat_comp_idx[Vertex(idx)];
      }
      v3d.first = flat_vpositions[Vertex(idx)];
      v3d.second = itr->second;
      return;
    }

    auto interior_itr = std::find_if(
        interior_points_projected.begin(), interior_points_projected.end(),
        [&v2d](const std::pair<Vector2f, int> &x) { return x.first == v2d; });
    if (interior_itr != interior_points_projected.end()) {
      idx = interior_itr - interior_points_projected.begin();
      v3d = interior_points_dedup[idx];
    } else {
      logger().warn("cdt_non_SO_component: Cannot find the vertex exactly");
      real_t min_dist = std::numeric_limits<real_t>::infinity();
      for (auto itr = boundary_points_projected.begin();
           itr != boundary_points_projected.end(); ++itr) {
        idx = itr - boundary_points_projected.begin();

        real_t dist = (v2d - itr->first).norm();
        if (dist < min_dist) {
          orig_idx = flat_orig_idx[Vertex(idx)];
          if (flat_orig_idx) {
            comp_idx = flat_orig_idx[Vertex(idx)];
          }
          v3d.first = flat_vpositions[Vertex(idx)];
          v3d.second = itr->second;
        }
      }

      for (auto itr = interior_points_projected.begin();
           itr != interior_points_projected.end(); ++itr) {
        idx = itr - interior_points_projected.begin();

        real_t dist = (v2d - itr->first).norm();
        if (dist < min_dist) {
          v3d = interior_points_dedup[idx];
        }
      }

      logger().warn("cdt_non_SO_component: Cannot find the vertex exactly. The "
                    "distance to the closest point is {}.",
                    min_dist);
    }
  };

  for (long i = 0; i < V2.rows(); i++) {
    std::pair<Vector3f, int> v3d;
    int orig_idx, comp_idx;

    find_3d_point(V2.row(i).transpose(), v3d, orig_idx, comp_idx);
    Vertex v = inflated_patch.add_vertex(v3d.first);

    // It's a boundary point
    inflated_orig_idx[v] = orig_idx;
    inflated_comp_idx[v] = comp_idx;
    inflated_orig_f_idx[v] = v3d.second;
  }

  // Fill other attributes
  auto inflated_patchID = inflated_patch.face_property<int>("f:patchID", -1);
  auto patch_component =
      inflated_patch.face_property<int>("f:patch_component", -1);
  for (long i = 0; i < F2.rows(); i++) {
    std::vector<Vertex> tri({Vertex(F2.coeff(i, 0)), Vertex(F2.coeff(i, 1)),
                             Vertex(F2.coeff(i, 2))});
    // Flip the surface if the target patch is back facing
    if (to_flip && target_facing == FacingType::BACK)
      tri = std::vector<Vertex>({Vertex(F2.coeff(i, 0)), Vertex(F2.coeff(i, 2)),
                                 Vertex(F2.coeff(i, 1))});
    auto f = inflated_patch.add_face(tri);
    inflated_patchID[f] = patch_id;
    patch_component[f] = comp_idx;
  }

  inflated_patch.init();

  auto out_boundary_type = inflated_patch.edge_property<BoundaryType>(
      "e:boundary_type", BoundaryType::NONE);
  for (size_t i = 0; i < inflated_patch.n_edges(); i++) {
    Edge e(i);
    if (inflated_patch.is_boundary(e))
      out_boundary_type[e] = BoundaryType::PATCH;
    else
      out_boundary_type[e] = BoundaryType::NONE;
  }
}
