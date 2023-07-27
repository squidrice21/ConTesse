// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "triangulate_wso_patches.h"

// Eigen printout formatting
#include <spdlog/fmt/ostr.h>
#include <unordered_map>
#include <unordered_set>

#include "common.h"
#include "shor.h"
#include "subdivide_contour_edges.h"
#include "tag_concave_edges.h"
#include "tag_cusp_facing.h"

void to_2d_mesh(Mesh const &mesh_3d, Eigen::MatrixXd const &V,
                Eigen::MatrixXi const &F,
                std::vector<size_t> const &V_idx_original, Mesh &mesh_out) {
  auto build_mesh = [&](Eigen::MatrixXd const &V, Eigen::MatrixXi const &F,
                        std::vector<size_t> const &V_idx_original,
                        Mesh &mesh_out) {
    mesh_out.clear();

    std::unordered_map<int, int> idx_to_original;

    for (long i = 0; i < V.rows(); i++) {
      Vertex v = mesh_out.add_vertex(Vector3f(V.coeff(i, 0), V.coeff(i, 1), 0));

      if (!V_idx_original.empty())
        idx_to_original[v.idx()] = V_idx_original[i];
    }
    for (long i = 0; i < F.rows(); i++) {
      std::vector<Vertex> tri({Vertex(F.coeff(i, 0)), Vertex(F.coeff(i, 1)),
                               Vertex(F.coeff(i, 2))});
      mesh_out.add_face(tri);
    }

    mesh_out.init();

    if (V_idx_original.empty())
      return;

    // Add indexin correspondence
    auto orig_idx = mesh_out.vertex_property<int>("v:orig_idx", -1);
    auto orig_f_idx = mesh_out.vertex_property<int>("v:orig_f_idx", -1);

    for (int i = 0; i < (int)mesh_out.n_vertices(); i++) {
      Vertex v(i);
      orig_idx[v] = idx_to_original[i];

      Vertex orig_v(idx_to_original[i]);
      Face v_f;
      auto hit = mesh_3d.halfedges(orig_v), hit_end = hit;
      Vertex next;
      do {
        if (mesh_3d.face(*hit).is_valid()) {
          v_f = mesh_3d.face(*hit);
          break;
        }
      } while (++hit != hit_end);
      contess_assert(v_f.is_valid());
      orig_f_idx[v] = v_f.idx();
    }
  };
  // Scale and center the mesh wrt the origin
  auto scale_mesh = [](Mesh &mesh, double scale) {
    auto center = mesh.boundingBox().center();
    for (auto v_itr = mesh.vertices_begin(); v_itr != mesh.vertices_end();
         ++v_itr) {
      mesh.position(*v_itr) -= center;
      mesh.position(*v_itr) *= scale;
    }
  };

  build_mesh(V, F, V_idx_original, mesh_out);

  // Scale to match the size of the 3D bbox
  double scale = mesh_3d.boundingBox().diagonal().norm() /
                 mesh_out.boundingBox().diagonal().norm();
  scale_mesh(mesh_out, scale);
}

/*==============================*/

void triangulate_wso_patches(Mesh const &mesh, Camera const &camera,
                             std::map<size_t, Mesh> &patch_triangulations,
                             bool do_refine,
                             std::unordered_set<int> selected_patches,
                             bool ff_only) {
  // Triangulate based on the rotation indices
  auto is_cusp = mesh.get_vertex_property<bool>("v:cusp");
  auto cusp_facing = mesh.get_vertex_property<FacingType>("v:cusp_facing");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  Vector2i viewport({camera.vpWidth(), camera.vpHeight()});

  // Make sure the mesh has patches extracted
  contess_assert_msg(!mesh.get_seams_lengths_const().empty(),
                     "Can't triangulate when no patch is extracted.");
  // Make sure the rotation indices are computed
  contess_assert_msg(is_cusp && cusp_facing,
                     "Can't triangulate when no cusp is extracted.");

  Eigen::MatrixXd P;
  Eigen::MatrixXd V_out;
  Eigen::MatrixXi F_out;

  Matrix3Xf polyline_3d;
  std::vector<Vector2f> polyline;
  std::vector<size_t> polyline_indices;
  std::vector<int> polyline_Rs;

  patch_triangulations.clear();

  // Traverse the patch
  for (auto const &patch_contour : mesh.get_const_patch_chains()) {
    polyline.clear();
    polyline_indices.clear();
    polyline_Rs.clear();

    int patch_id = patch_contour.first;

    if (selected_patches.empty() || selected_patches.count(patch_id))
      logger().info("WSO on {}.", patch_id);
    if (!selected_patches.empty() && !selected_patches.count(patch_id)) {
      continue;
    }

    FacingType facing = mesh.get_patch_facing(patch_id);

    if (ff_only && facing != FacingType::FRONT)
      continue;

    // Traverse the chain
    for (auto const &he : *patch_contour.second) {
      Vertex v = mesh.from_vertex(he);
      polyline_indices.emplace_back(v.idx());

      int r = 0;
      if (is_cusp[v] && facing != cusp_facing[v] &&
          cusp_facing[v] != FacingType::UNDEFINED) {
        r = 1;
      }
      polyline_Rs.emplace_back(r);
    }

    if (polyline_indices.empty())
      continue;

    // Deduplicate the polygon vertices
    {
      std::vector<size_t> dedup_polyline_indices;
      std::vector<int> dedup_polyline_Rs;
      dedup_polyline_indices.reserve(polyline_indices.size());
      dedup_polyline_Rs.reserve(polyline_indices.size());

      dedup_polyline_indices.emplace_back(polyline_indices[0]);
      dedup_polyline_Rs.emplace_back(polyline_Rs[0]);
      for (size_t k = 1; k < polyline_indices.size(); k++) {
        if (dedup_polyline_indices.back() == polyline_indices[k])
          continue;
        dedup_polyline_indices.emplace_back(polyline_indices[k]);
        dedup_polyline_Rs.emplace_back(polyline_Rs[k]);
      }
      // The polygon is closed by default.
      // No need to match the first and the last vertices.
      if (dedup_polyline_indices.size() > 1 &&
          (dedup_polyline_indices.front() == dedup_polyline_indices.back())) {
        dedup_polyline_indices.pop_back();
        dedup_polyline_Rs.pop_back();
      }
      polyline_indices = dedup_polyline_indices;
      polyline_Rs = dedup_polyline_Rs;

      polyline_3d.resize(3, polyline_indices.size());
      size_t polyline_3d_i = 0;
      for (auto v : polyline_indices) {
        polyline_3d.col(polyline_3d_i++) = vpositions[Vertex(v)];
      }
    }

    // Project to 2D
    projectToViewport(polyline_3d, polyline, camera.projectionMatrix(),
                      camera.viewMatrix(), viewport);

    // WSO triangulation
    P.resize(polyline.size(), 2);
    for (size_t k = 0; k < polyline.size(); k++) {
      P.row(k) = polyline[k].cast<double>();
    }
    Eigen::VectorXi R =
        Eigen::Map<Eigen::VectorXi>(polyline_Rs.data(), polyline_Rs.size(), 1);

    contess_assert_msg(P.rows() > 0, "Empty polygon.");

    // Tessellate polygon
    bool is_wso_succeeded = Shor_van_wyck(P, R, "", V_out, F_out, do_refine);
    if (!is_wso_succeeded) {
      logger().error("Weber-Zorin failed on patch {}. The rotation indices may "
                     "be incorrect.",
                     patch_id);
    }

    if (do_refine) {
      logger().warn("Resulting flat patch is refined. The original indexing is "
                    "no longer valid.");
    }

    patch_triangulations[patch_id] = Mesh();
    to_2d_mesh(mesh, V_out, F_out, polyline_indices,
               patch_triangulations[patch_id]);

    auto patch_patchID =
        patch_triangulations[patch_id].face_property<int>("f:patchID");
    patch_patchID.vector().assign(patch_patchID.vector().size(), patch_id);

    V_out.resize(0, 0);
    F_out.resize(0, 0);
  }
}

void triangulate_wso_patches_adaptive(
    Mesh &mesh, Camera const &camera,
    std::map<size_t, Mesh> &patch_triangulations, bool do_refine,
    std::unordered_set<int> selected_patches) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");

  // Call the WSO using the current contour polygons
  triangulate_wso_patches(mesh, camera, patch_triangulations, do_refine,
                          selected_patches);

  // Check success and if any one failed, refine the contour of failed ones
  int max_refine_itr = 2;
  int itr = 0;
  bool all_wso_returned = true;
  do {
    all_wso_returned = true;
    std::unordered_set<int> added_edges; // Avoid double adding
    std::vector<Edge> subdiv_edges;
    for (auto const &triangulation : patch_triangulations) {
      // Failed
      if (triangulation.second.n_faces() == 0) {
        all_wso_returned = false;

        // Add all contour edges to the upsampling list
        auto patch_chains =
            mesh.get_const_patch_chains().find(triangulation.first);
        for (auto const chain_edge : *patch_chains->second) {
          Edge e = mesh.edge(chain_edge);

          // It's possible that the patch boundary is a cut edge
          if (is_contour[e] >= 0 && !added_edges.count(e.idx())) {
            // Avoid subdividing edges that are already short
            if ((vpositions[mesh.from_vertex(chain_edge)] -
                 vpositions[mesh.to_vertex(chain_edge)])
                    .norm() < MIN_EDGE_LENGTH) {
              continue;
            }

            added_edges.emplace(e.idx());
            subdiv_edges.emplace_back(e);
          }
        }
      }
    }

    // Refine and rerun chaining and WSO
    if (!subdiv_edges.empty()) {
      logger().info("Itr {}: Contour refinement and WSO.", itr);

      subdivide_contour_edges(mesh, camera, subdiv_edges);
      tag_concave_edges_side(mesh, camera, ContourMode::VBO_CONTOUR);
      tag_cusp_given_concave_edges(mesh, camera);
      tag_cusp_facing(mesh, camera);

      serialization::SerializationObject serial{
          &mesh, &camera, std::vector<serialization::SerializationRay>()};

      logger().info("Serializating model to {}", "./contour_subdivided.json");
      std::ofstream fstr("./contour_subdivided.json");
      serialization::save(fstr, serial, true, "json");
      fstr.flush();
      fstr.close();

      // Rerun contour orientation
      chain_contour(mesh, camera, true, false);
      chain_cut_graph(mesh, camera, false);

      patch_triangulations.clear();
      triangulate_wso_patches(mesh, camera, patch_triangulations, do_refine,
                              selected_patches);
    }

    itr++;
  } while (!all_wso_returned && itr < max_refine_itr);

  if (!all_wso_returned) {
    logger().warn("WSO fails even with {} rounds of refinement.",
                  max_refine_itr);
  }
}

void triangulate_wso_patches_simplified(
    Mesh &mesh, Camera const &camera,
    std::map<size_t, Mesh> &patch_triangulations, bool do_refine,
    std::unordered_set<int> selected_patches) {
  // Triangulate based on the rotation indices
  auto is_cusp = mesh.get_vertex_property<bool>("v:sim_cusp");
  auto cusp_facing = mesh.get_vertex_property<FacingType>("v:sim_cusp_facing");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  Vector2i viewport({camera.vpWidth(), camera.vpHeight()});

  // Make sure the mesh has patches extracted
  contess_assert_msg(!mesh.get_seams_lengths_const().empty(),
                     "Can't triangulate when no patch is extracted.");
  // Make sure the rotation indices are computed
  contess_assert_msg(is_cusp && cusp_facing,
                     "Can't triangulate when no cusp is extracted.");

  Eigen::MatrixXd P;
  Eigen::MatrixXd V_out;
  Eigen::MatrixXi F_out;

  Matrix3Xf polyline_3d;
  std::vector<Vector2f> polyline;
  std::vector<size_t> polyline_indices;
  std::vector<int> polyline_Rs;

  patch_triangulations.clear();

  // Traverse the patch
  for (auto const &patch_contour : mesh.get_const_patch_simplified_chains()) {
    polyline.clear();
    polyline_indices.clear();
    polyline_Rs.clear();

    int patch_id = patch_contour.first;

    if (selected_patches.empty() || selected_patches.count(patch_id))
      logger().info("WSO on {}.", patch_id);
    if (!selected_patches.empty() && !selected_patches.count(patch_id)) {
      // logger().info("\tSkip WSO on {}.", patch_id);
      continue;
    }

    FacingType facing = mesh.get_patch_facing(patch_id);

    // Traverse the chain
    for (auto const &vhe : *patch_contour.second) {
      Vertex v;
      if (vhe.he.is_valid()) {
        v = mesh.from_vertex(vhe.he);
      } else {
        contess_assert_msg(vhe.from_v.is_valid() && vhe.to_v.is_valid(),
                           "Invalid simplified contour edge");
        v = vhe.from_v;
      }
      contess_assert_msg(v.is_valid(), "Invalid simplified contour edge");

      polyline_indices.emplace_back(v.idx());

      int r = 0;
      if (is_cusp[v] && facing != cusp_facing[v]) {
        r = 1;
      }
      polyline_Rs.emplace_back(r);
    }

    if (polyline_indices.empty())
      continue;

    // Deduplicate the polygon vertices
    {
      std::vector<size_t> dedup_polyline_indices;
      std::vector<int> dedup_polyline_Rs;
      dedup_polyline_indices.reserve(polyline_indices.size());
      dedup_polyline_Rs.reserve(polyline_indices.size());

      dedup_polyline_indices.emplace_back(polyline_indices[0]);
      dedup_polyline_Rs.emplace_back(polyline_Rs[0]);
      for (size_t k = 1; k < polyline_indices.size(); k++) {
        if (dedup_polyline_indices.back() == polyline_indices[k])
          continue;
        dedup_polyline_indices.emplace_back(polyline_indices[k]);
        dedup_polyline_Rs.emplace_back(polyline_Rs[k]);
      }
      // The polygon is closed by default.
      // No need to match the first and the last vertices.
      if (dedup_polyline_indices.size() > 1 &&
          (dedup_polyline_indices.front() == dedup_polyline_indices.back())) {
        dedup_polyline_indices.pop_back();
        dedup_polyline_Rs.pop_back();
      }
      polyline_indices = dedup_polyline_indices;
      polyline_Rs = dedup_polyline_Rs;

      polyline_3d.resize(3, polyline_indices.size());
      size_t polyline_3d_i = 0;
      for (auto v : polyline_indices) {
        polyline_3d.col(polyline_3d_i++) = vpositions[Vertex(v)];
      }
    }

    // Project to 2D
    projectToViewport(polyline_3d, polyline, camera.projectionMatrix(),
                      camera.viewMatrix(), viewport);

    // WSO triangulation
    P.resize(polyline.size(), 2);
    for (size_t k = 0; k < polyline.size(); k++) {
      P.row(k) = polyline[k].cast<double>();
    }
    Eigen::VectorXi R =
        Eigen::Map<Eigen::VectorXi>(polyline_Rs.data(), polyline_Rs.size(), 1);

    // Tessellate polygon
    bool is_wso_succeeded = Shor_van_wyck(P, R, "", V_out, F_out, do_refine);
    if (!is_wso_succeeded) {
      logger().error("Weber-Zorin failed on patch {}. The rotation indices may "
                     "be incorrect.",
                     patch_id);
    }

    if (do_refine) {
      logger().warn("Resulting flat patch is refined. The original indexing is "
                    "no longer valid.");
    }

    patch_triangulations[patch_id] = Mesh();
    to_2d_mesh(mesh, V_out, F_out, polyline_indices,
               patch_triangulations[patch_id]);

    auto patch_patchID =
        patch_triangulations[patch_id].face_property<int>("f:patchID");
    patch_patchID.vector().assign(patch_patchID.vector().size(), patch_id);

    V_out.resize(0, 0);
    F_out.resize(0, 0);
  }
}
