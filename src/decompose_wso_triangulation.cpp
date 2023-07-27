// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "decompose_wso_triangulation.h"
#include "common.h"

// Eigen printout formatting
#include <igl/predicates/predicates.h>
#include <limits>
#include <spdlog/fmt/ostr.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

bool is_simple_polygon(Mesh const &mesh, Camera const &camera,
                       Mesh const &flat_patch,
                       std::vector<Edge> const &left_poly,
                       std::vector<Edge> const &right_poly) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto orig_idx = flat_patch.get_vertex_property<int>("v:orig_idx");
  contess_assert_msg(
      orig_idx, "decompose_wso_triangulation: Missing 2D-3D correspondence.");

  Vector2i viewport({camera.vpWidth(), camera.vpHeight()});

  auto project2d = [&](Vector3f const &v3d, Vector2f &pos2D) {
    pos2D = project(v3d, camera.viewMatrix().matrix(),
                    camera.projectionMatrix(), viewport)
                .head<2>();
  };

  for (auto const &l_e : left_poly) {
    for (auto const &r_e : right_poly) {
      // Check if sharing a vertex
      if (flat_patch.vertex(l_e, 0) == flat_patch.vertex(r_e, 0) ||
          flat_patch.vertex(l_e, 1) == flat_patch.vertex(r_e, 0) ||
          flat_patch.vertex(l_e, 0) == flat_patch.vertex(r_e, 1) ||
          flat_patch.vertex(l_e, 1) == flat_patch.vertex(r_e, 1))
        continue;

      // Check intersection
      Vector2f l0, l1, r0, r1;
      project2d(vpositions[Vertex(orig_idx[flat_patch.vertex(l_e, 0)])], l0);
      project2d(vpositions[Vertex(orig_idx[flat_patch.vertex(l_e, 1)])], l1);
      project2d(vpositions[Vertex(orig_idx[flat_patch.vertex(r_e, 0)])], r0);
      project2d(vpositions[Vertex(orig_idx[flat_patch.vertex(r_e, 1)])], r1);

      auto side1 = igl::predicates::orient2d(l0, l1, r0);
      auto side2 = igl::predicates::orient2d(l0, l1, r1);

      auto side21 = igl::predicates::orient2d(r0, r1, l0);
      auto side22 = igl::predicates::orient2d(r0, r1, l1);

      // Four points colinear
      if (side1 == igl::predicates::Orientation::COLLINEAR &&
          side2 == igl::predicates::Orientation::COLLINEAR &&
          side21 == igl::predicates::Orientation::COLLINEAR &&
          side22 == igl::predicates::Orientation::COLLINEAR) {
        Vector2f l_dir = (l1 - l0).normalized();
        Vector2f l_norm(-l_dir.y(), l_dir.x());
        Vector2f l0_up = l0 + l_norm;
        Vector2f l1_up = l1 + l_norm;
        auto side31 = igl::predicates::orient2d(l0, l0_up, r0);
        auto side32 = igl::predicates::orient2d(l1, l1_up, r1);

        if (side31 == igl::predicates::Orientation::COLLINEAR ||
            side32 == igl::predicates::Orientation::COLLINEAR ||
            (side31 == igl::predicates::Orientation::NEGATIVE) !=
                (side32 == igl::predicates::Orientation::NEGATIVE))
          return false;

        continue;
      }

      if (((side1 == igl::predicates::Orientation::COLLINEAR ||
            side2 == igl::predicates::Orientation::COLLINEAR ||
            (side1 == igl::predicates::Orientation::NEGATIVE) !=
                (side2 == igl::predicates::Orientation::NEGATIVE)) &&
           (side21 == igl::predicates::Orientation::COLLINEAR ||
            side22 == igl::predicates::Orientation::COLLINEAR ||
            (side21 == igl::predicates::Orientation::NEGATIVE) !=
                (side22 == igl::predicates::Orientation::NEGATIVE)))) {
        return false;
      }
    }
  }

  return true;
}

void decompose_wso_triangulation(Mesh const &mesh, Camera const &camera,
                                 Mesh &flat_patch) {
  auto boundary_type =
      flat_patch.get_edge_property<BoundaryType>("e:boundary_type");
  contess_assert_msg(boundary_type,
                     "decompose_wso_triangulation: Missing initial patch edge "
                     "labeling from label_patch_boundary.");
  auto patch_component = flat_patch.face_property<int>("f:patch_component", -1);
  patch_component.vector().assign(patch_component.vector().size(), -1);

  // Iteratively merging component triangles
  std::unordered_map<int, std::vector<int>> component_f;
  for (size_t i = 0; i < flat_patch.n_faces(); i++) {
    patch_component[Face(i)] = i;
    component_f[i] = std::vector<int>({Face(i).idx()});
  }

  // Initialize all interior edges to be component boundary
  std::deque<Edge> processEdges;
  for (size_t i = 0; i < flat_patch.n_edges(); i++) {
    Edge e_flat(i);
    if (boundary_type[e_flat] == BoundaryType::PATCH ||
        boundary_type[e_flat] == BoundaryType::CUT)
      continue;
    boundary_type[e_flat] = BoundaryType::COMPONENT;
    processEdges.push_back(e_flat);
  }

  while (!processEdges.empty()) {
    Edge e_test = processEdges.front();
    processEdges.pop_front();

    // Get left and right components
    std::vector<Edge> left_poly, right_poly;
    int left_idx = patch_component[flat_patch.face(e_test, 0)],
        right_idx = patch_component[flat_patch.face(e_test, 1)];
    contess_assert_msg(
        left_idx != right_idx,
        "decompose_wso_triangulation: Edge is not a component boundary.");

    for (size_t i = 0; i < flat_patch.n_edges(); i++) {
      Edge e_flat(i);
      if (e_flat == e_test || boundary_type[e_flat] == BoundaryType::NONE)
        continue;

      Face l_f = flat_patch.face(e_flat, 0);
      Face r_f = flat_patch.face(e_flat, 1);
      if ((l_f.is_valid() && patch_component[l_f] == left_idx) ||
          (r_f.is_valid() && patch_component[r_f] == left_idx)) {
        left_poly.emplace_back(e_flat);
        contess_assert_msg(
            !((l_f.is_valid() && patch_component[l_f] == right_idx) ||
              (r_f.is_valid() && patch_component[r_f] == right_idx)),
            "decompose_wso_triangulation: Components connected by multiple "
            "edges.");
      }
      if ((l_f.is_valid() && patch_component[l_f] == right_idx) ||
          (r_f.is_valid() && patch_component[r_f] == right_idx)) {
        right_poly.emplace_back(e_flat);
        contess_assert_msg(
            !((l_f.is_valid() && patch_component[l_f] == left_idx) ||
              (r_f.is_valid() && patch_component[r_f] == left_idx)),
            "decompose_wso_triangulation: Components connected by multiple "
            "edges.");
      }
    }

    // Test if the resulting polygon is simple
    bool is_simple =
        is_simple_polygon(mesh, camera, flat_patch, left_poly, right_poly);

    if (is_simple) {
      boundary_type[e_test] = BoundaryType::NONE;
      for (auto const f_i : component_f[left_idx]) {
        patch_component[Face(f_i)] = right_idx;
      }
      component_f[right_idx].insert(component_f[right_idx].end(),
                                    component_f[left_idx].begin(),
                                    component_f[left_idx].end());
      component_f.erase(left_idx);
    }
  }

  // Reassign the component index for debugging
  std::unordered_map<int, int> comp_map;
  std::unordered_map<int, std::vector<Face>> comp2faces;
  for (size_t i = 0; i < flat_patch.n_faces(); i++) {
    Face f(i);
    if (!comp_map.count(patch_component[f])) {
      comp_map[patch_component[f]] = comp_map.size();
    }
    patch_component[f] = comp_map[patch_component[f]];
    if (!comp2faces.count(patch_component[f]))
      comp2faces[patch_component[f]] = std::vector<Face>();
    comp2faces[patch_component[f]].emplace_back(f);
  }

  // Check degeneration as areas
  auto vpositions = flat_patch.get_vertex_property<Vector3f>("v:point");
  auto flat_patchID = flat_patch.get_face_property<int>("f:patchID");
  int pid = flat_patchID[*flat_patch.faces_begin()];
  auto flat_comp_valid = flat_patch.get_face_property<bool>("f:comp_valid");
  if (!flat_comp_valid) {
    flat_patch.add_face_property<bool>("f:comp_valid", true);
  }
  flat_comp_valid = flat_patch.face_property<bool>("f:comp_valid");
  flat_comp_valid.vector().assign(flat_comp_valid.vector().size(), true);

  for (auto const &comp : comp2faces) {
    real_t comp_area = 0;
    for (auto const &f : comp.second) {
      Vertex vv[3];
      flat_patch.verticesOfFace(f, vv);
      comp_area += (vpositions[vv[1]] - vpositions[vv[0]])
                       .cross(vpositions[vv[2]] - vpositions[vv[0]])
                       .norm();
    }
    comp_area /= 2;
    logger().info("Decompose patch {}, comp {}: Area {}", pid, comp.first,
                  comp_area);
    if (comp_area < std::numeric_limits<real_t>::epsilon()) {
      for (auto const &f : comp.second) {
        flat_comp_valid[f] = false;
      }
    }
  }
}
