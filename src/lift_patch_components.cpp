// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "lift_patch_components.h"

// Eigen printout formatting
#include <spdlog/fmt/ostr.h>

#include "cdt_non_SO_component.h"
#include "common.h"
#include "inflat_path.h"
#include "inflat_path_sparse.h"
#include "io/serialization.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

void trace_pslg(Mesh const &boundary_inflated_patch, Vertex start, int comp_idx,
                Mesh &patch_comp) {
  auto patch_component =
      boundary_inflated_patch.get_face_property<int>("f:patch_component");
  contess_assert_msg(patch_component,
                     "trace_pslg: Missing decomposed patch component.");
  auto out_patch_component =
      patch_comp.face_property<int>("f:patch_component", -1);

  auto vpositions =
      boundary_inflated_patch.get_vertex_property<Vector3f>("v:point");
  auto comp_boundary_type =
      patch_comp.edge_property<BoundaryType>("e:boundary_type");
  // Assign unique indices for the component boundary vertices
  // so we can stitch them back later
  auto comp_boundary_idx = patch_comp.vertex_property<int>("v:comp_idx");
  auto orig_idx = patch_comp.vertex_property<int>("v:orig_idx");
  auto orig_f_idx = patch_comp.vertex_property<int>("v:orig_f_idx");
  auto in_orig_idx =
      boundary_inflated_patch.get_vertex_property<int>("v:orig_idx");
  auto in_orig_f_idx =
      boundary_inflated_patch.get_vertex_property<int>("v:orig_f_idx");
  contess_assert_msg(in_orig_idx && in_orig_f_idx,
                     "trace_pslg: Missing index correspondence.");

  auto in_flat_patchID =
      boundary_inflated_patch.get_face_property<int>("f:patchID");
  contess_assert_msg(in_flat_patchID, "trace_pslg: Missing patch index.");

  Vertex cur = start;
  Vertex prev_comp = patch_comp.add_vertex(vpositions[cur]);

  comp_boundary_idx[prev_comp] = cur.idx();
  orig_idx[prev_comp] = in_orig_idx[cur];
  orig_f_idx[prev_comp] = in_orig_f_idx[cur];
  Vertex fixed_comp = prev_comp;
  do {
    auto hit = boundary_inflated_patch.halfedges(cur), hit_end = hit;
    Vertex next;
    do {
      if (boundary_inflated_patch.face(*hit).is_valid() &&
          patch_component[boundary_inflated_patch.face(*hit)] == comp_idx) {
        next = boundary_inflated_patch.to_vertex(*hit);
        break;
      }
    } while (++hit != hit_end);
    contess_assert_msg(next.is_valid(), "trace_pslg: Boundary is not loopy.");

    if (next == start) {
      comp_boundary_type[patch_comp.find_edge(prev_comp, fixed_comp)] =
          BoundaryType::PATCH;
      break;
    }

    Vertex cur_comp = patch_comp.add_vertex(vpositions[next]);

    comp_boundary_idx[cur_comp] = next.idx();
    orig_idx[cur_comp] = in_orig_idx[next];
    orig_f_idx[cur_comp] = in_orig_f_idx[next];

    if (patch_comp.n_vertices() > 2 && next != start) {
      auto f = patch_comp.add_face(
          std::vector<Vertex>({fixed_comp, prev_comp, cur_comp}));
      out_patch_component[f] = comp_idx;

      comp_boundary_type[patch_comp.find_edge(prev_comp, cur_comp)] =
          BoundaryType::PATCH;
      comp_boundary_type[patch_comp.find_edge(fixed_comp, cur_comp)] =
          BoundaryType::NONE;
      comp_boundary_type[patch_comp.find_edge(prev_comp, fixed_comp)] =
          BoundaryType::NONE;
    }

    // To handle the degenerated triangles (first and last)
    if (patch_comp.n_vertices() == 3) {
      comp_boundary_type[patch_comp.find_edge(prev_comp, fixed_comp)] =
          BoundaryType::PATCH;
    }

    cur = next;
    prev_comp = cur_comp;
  } while (cur != start);

  auto out_patchID = patch_comp.face_property<int>("f:patchID");
  if (boundary_inflated_patch.n_faces() > 0)
    out_patchID.vector().assign(out_patchID.vector().size(),
                                in_flat_patchID[Face(0)]);
}

void lift_patch_components(Mesh const &mesh, Camera const &camera,
                           real_t sampling_delta, real_t min_delta_2d,
                           Mesh const &decomp_blift_patch, Mesh &inflat_patch,
                           bool sparse_sampling) {
  auto patch_component =
      decomp_blift_patch.get_face_property<int>("f:patch_component");
  contess_assert_msg(
      patch_component,
      "lift_patch_components: Missing decomposed patch component.");
  auto boundary_type =
      decomp_blift_patch.get_edge_property<BoundaryType>("e:boundary_type");
  contess_assert_msg(boundary_type,
                     "lift_patch_components: Missing initial patch edge "
                     "labeling from label_patch_boundary.");
  // Inflat component boundary
  Mesh boundary_inflated_patch;
  std::unordered_set<BoundaryType> to_inflat_labels;
  to_inflat_labels.emplace(BoundaryType::COMPONENT);
  std::vector<std::pair<int, int>> failed_paths;
  inflat_path(mesh, decomp_blift_patch, camera, sampling_delta, min_delta_2d,
              boundary_inflated_patch, failed_paths, to_inflat_labels);

  if (!inflat_patch.get_vertex_property<int>("v:orig_idx"))
    inflat_patch.add_vertex_property<int>("v:orig_idx", -1);
  auto inflated_orig_idx = inflat_patch.vertex_property<int>("v:orig_idx");
  if (!inflat_patch.get_vertex_property<int>("v:orig_f_idx"))
    inflat_patch.add_vertex_property<int>("v:orig_f_idx", -1);
  auto inflated_orig_f_idx = inflat_patch.vertex_property<int>("v:orig_f_idx");

  std::vector<Mesh> comp_inflated_patches;
  std::unordered_set<int> visited_comp;
  for (size_t i = 0; i < decomp_blift_patch.n_edges(); i++) {
    Edge e(i);
    if (boundary_type[e] == BoundaryType::NONE)
      continue;

    std::vector<int> adj_comp;
    int l_comp = (decomp_blift_patch.face(e, 0).is_valid())
                     ? patch_component[decomp_blift_patch.face(e, 0)]
                     : -1;
    int r_comp = (decomp_blift_patch.face(e, 1).is_valid())
                     ? patch_component[decomp_blift_patch.face(e, 1)]
                     : -1;
    if (l_comp >= 0 && !visited_comp.count(l_comp))
      adj_comp.emplace_back(l_comp);
    if (r_comp >= 0 && !visited_comp.count(r_comp))
      adj_comp.emplace_back(r_comp);

    for (int comp_idx : adj_comp) {
      visited_comp.emplace(comp_idx);
      std::unordered_set<int> to_inflat_comp_labels;
      to_inflat_comp_labels.emplace(comp_idx);

      // Get interior samples per component
      std::vector<std::pair<Vector3f, int>> inflat_samples;
      if (sparse_sampling)
        inflat_path_sparse(mesh, decomp_blift_patch, camera, sampling_delta,
                           min_delta_2d, inflat_samples, to_inflat_comp_labels);
      else
        inflat_path(mesh, decomp_blift_patch, camera, sampling_delta,
                    min_delta_2d, inflat_samples, to_inflat_comp_labels);

      // Get corresponding PSLG by tracing from a boundary vertex
      Mesh patch_comp;
      Vertex start = decomp_blift_patch.vertex(e, 0);
      // Note here boundary_inflated_patch and decomp_blift_patch are indexed in
      // the sample way for the shared vertices
      trace_pslg(boundary_inflated_patch, start, comp_idx, patch_comp);

      serialization::SerializationObject pslg_serial{
          &patch_comp, &camera, std::vector<serialization::SerializationRay>()};

      // CDT each component
      Mesh comp_inflated_patch;
      cdt_non_SO_component(mesh, patch_comp, camera, inflat_samples,
                           comp_inflated_patch);

      serialization::SerializationObject serial{
          &comp_inflated_patch, &camera,
          std::vector<serialization::SerializationRay>()};

      comp_inflated_patches.emplace_back(comp_inflated_patch);
    }
  }
  // Stitch back
  std::unordered_map<int, int> stitch_comp_b_v_to_new_v;
  for (auto const &comp_inflated_patch : comp_inflated_patches) {
    auto comp_boundary_idx =
        comp_inflated_patch.get_vertex_property<int>("v:comp_idx");
    contess_assert_msg(comp_boundary_idx,
                       "lift_patch_components: Index correspondence for "
                       "stitching is missing.");
    auto orig_boundary_idx =
        comp_inflated_patch.get_vertex_property<int>("v:orig_idx");
    auto orig_f_idx =
        comp_inflated_patch.get_vertex_property<int>("v:orig_f_idx");
    contess_assert_msg(orig_f_idx && orig_boundary_idx,
                       "lift_patch_components: Index correspondence for "
                       "stitching is missing.");
    auto vpositions =
        comp_inflated_patch.get_vertex_property<Vector3f>("v:point");

    // This one is used to add faces
    std::unordered_map<int, int> comp_v_to_new_v;

    for (size_t j = 0; j < comp_inflated_patch.n_vertices(); j++) {
      Vertex comp_v(j);

      // This vertex needs to be stitched
      if (comp_boundary_idx[comp_v] >= 0) {
        int comp_b_idx = comp_boundary_idx[comp_v];
        if (!stitch_comp_b_v_to_new_v.count(comp_b_idx)) {
          Vertex v = inflat_patch.add_vertex(vpositions[comp_v]);
          comp_v_to_new_v[j] = v.idx();
          stitch_comp_b_v_to_new_v[comp_b_idx] = v.idx();
          inflated_orig_idx[v] = orig_boundary_idx[comp_v];
          inflated_orig_f_idx[v] = orig_f_idx[comp_v];
        }
        comp_v_to_new_v[j] = stitch_comp_b_v_to_new_v[comp_b_idx];
      } else {
        Vertex v = inflat_patch.add_vertex(vpositions[comp_v]);
        comp_v_to_new_v[j] = v.idx();
        inflated_orig_idx[v] = orig_boundary_idx[comp_v];
        inflated_orig_f_idx[v] = orig_f_idx[comp_v];
      }
    }

    for (size_t j = 0; j < comp_inflated_patch.n_faces(); j++) {
      Vertex vv[3];
      comp_inflated_patch.verticesOfFace(Face(j), vv);

      contess_assert_msg(
          comp_v_to_new_v.count(vv[0].idx()) &&
              comp_v_to_new_v.count(vv[1].idx()) &&
              comp_v_to_new_v.count(vv[2].idx()),
          "lift_patch_components: Face contains unseen vertices.");

      std::vector<Vertex> new_f_vv({Vertex(comp_v_to_new_v[vv[0].idx()]),
                                    Vertex(comp_v_to_new_v[vv[1].idx()]),
                                    Vertex(comp_v_to_new_v[vv[2].idx()])});
      inflat_patch.add_face(new_f_vv);
    }
  }

  // Add boundary information for visualization
  if (!inflat_patch.get_edge_property<bool>("e:is_boundary"))
    inflat_patch.add_edge_property<bool>("e:is_boundary", false);
  auto is_boundary = inflat_patch.edge_property<bool>("e:is_boundary");
  is_boundary.vector().assign(is_boundary.vector().size(), false);
  for (size_t j = 0; j < inflat_patch.n_edges(); j++) {
    Edge e(j);
    is_boundary[e] = inflat_patch.is_boundary(e);
  }
}

void lift_patch_components_iterative(Mesh const &mesh, Camera const &camera,
                                     real_t sampling_delta, real_t min_delta_2d,
                                     real_t max_distance,
                                     Mesh const &decomp_blift_patch,
                                     Mesh &inflat_patch,
                                     size_t max_inflation_itr) {
  auto decomp_patchID = decomp_blift_patch.get_face_property<int>("f:patchID");
  auto patchID = inflat_patch.add_face_property<int>("f:patchID", -1);
  auto patch_component =
      decomp_blift_patch.get_face_property<int>("f:patch_component");
  contess_assert_msg(
      patch_component && decomp_patchID,
      "lift_patch_components: Missing decomposed patch component.");
  auto boundary_type =
      decomp_blift_patch.get_edge_property<BoundaryType>("e:boundary_type");
  contess_assert_msg(boundary_type,
                     "lift_patch_components: Missing initial patch edge "
                     "labeling from label_patch_boundary.");
  int decomp_patch_id = decomp_patchID[*(decomp_blift_patch.faces_begin())];

  // Inflat component boundary
  Mesh boundary_inflated_patch;
  std::unordered_set<BoundaryType> to_inflat_labels;
  to_inflat_labels.emplace(BoundaryType::COMPONENT);
  std::vector<std::pair<int, int>> failed_paths;
  inflat_path(mesh, decomp_blift_patch, camera, sampling_delta, min_delta_2d,
              boundary_inflated_patch, failed_paths, to_inflat_labels);

  if (!inflat_patch.get_vertex_property<int>("v:orig_idx"))
    inflat_patch.add_vertex_property<int>("v:orig_idx", -1);
  auto inflated_orig_idx = inflat_patch.vertex_property<int>("v:orig_idx");
  if (!inflat_patch.get_vertex_property<int>("v:orig_f_idx"))
    inflat_patch.add_vertex_property<int>("v:orig_f_idx", -1);
  auto inflated_orig_f_idx = inflat_patch.vertex_property<int>("v:orig_f_idx");

  std::vector<Mesh> comp_inflated_patches;
  std::unordered_set<int> visited_comp;
  for (size_t i = 0; i < decomp_blift_patch.n_edges(); i++) {
    Edge e(i);
    if (boundary_type[e] == BoundaryType::NONE)
      continue;

    std::vector<int> adj_comp;
    int l_comp = (decomp_blift_patch.face(e, 0).is_valid())
                     ? patch_component[decomp_blift_patch.face(e, 0)]
                     : -1;
    int r_comp = (decomp_blift_patch.face(e, 1).is_valid())
                     ? patch_component[decomp_blift_patch.face(e, 1)]
                     : -1;
    if (l_comp >= 0 && !visited_comp.count(l_comp))
      adj_comp.emplace_back(l_comp);
    if (r_comp >= 0 && !visited_comp.count(r_comp))
      adj_comp.emplace_back(r_comp);

    for (int comp_idx : adj_comp) {
      visited_comp.emplace(comp_idx);
      std::unordered_set<int> to_inflat_comp_labels;
      to_inflat_comp_labels.emplace(comp_idx);

      // Get corresponding PSLG by tracing from a boundary vertex
      Mesh patch_comp;
      Vertex start = decomp_blift_patch.vertex(e, 0);
      // Note here boundary_inflated_patch and decomp_blift_patch are indexed
      // in the sample way for the shared vertices
      trace_pslg(boundary_inflated_patch, start, comp_idx, patch_comp);

      serialization::SerializationObject pslg_serial{
          &patch_comp, &camera, std::vector<serialization::SerializationRay>()};
      std::string pslg_file_name = "./test_pslg_p" +
                                   std::to_string(decomp_patch_id) + "_c" +
                                   std::to_string(comp_idx) + ".json";

      Mesh cur_inflated_patch;
      size_t inflation_itr = 0;
      std::vector<std::pair<Vector3f, int>> prev_inflat_samples;

      // Check if the component is degenerated
      auto comp_valid =
          decomp_blift_patch.get_face_property<bool>("f:comp_valid");
      Face check_f =
          (decomp_blift_patch.face(e, 0).is_valid() &&
           patch_component[decomp_blift_patch.face(e, 0)] == comp_idx)
              ? decomp_blift_patch.face(e, 0)
              : decomp_blift_patch.face(e, 1);
      bool is_comp_valid = comp_valid[check_f];

      Mesh comp_inflated_patch;

      if (is_comp_valid) {
        do {
          // Get interior samples per component
          std::vector<std::pair<Vector3f, int>> inflat_samples;
          if (inflation_itr == 0) {
            inflat_path_iterative(mesh, decomp_blift_patch, camera,
                                  sampling_delta, min_delta_2d, inflat_samples,
                                  to_inflat_comp_labels, max_distance);
          } else {
            inflat_path_iterative(mesh, cur_inflated_patch, camera,
                                  sampling_delta, min_delta_2d, inflat_samples,
                                  to_inflat_comp_labels, max_distance);
            auto inflated_boundary_type =
                cur_inflated_patch.get_edge_property<BoundaryType>(
                    "e:boundary_type");
            contess_assert_msg(
                inflated_boundary_type,
                "lift_patch_components: Missing initial patch edge "
                "labeling from label_patch_boundary.");
            auto is_end_sample_needed = [&](Vertex end_v) -> bool {
              auto hit = cur_inflated_patch.halfedges(end_v), hit_end = hit;
              do {
                // Check if the vertex is in the interior
                if (!(*hit).is_valid() ||
                    !cur_inflated_patch.edge(*hit).is_valid()) {
                  logger().warn("lift_patch_components: {} wrongly connected.",
                                end_v);
                  return false;
                }
                if (inflated_boundary_type[cur_inflated_patch.edge(*hit)] !=
                    BoundaryType::NONE)
                  return false;
              } while (++hit != hit_end);
              return true;
            };
            if (inflat_samples.empty()) {
              break;
            }
            auto inflated_orig_f_idx =
                cur_inflated_patch.get_vertex_property<int>("v:orig_f_idx");
            auto inflated_vpost =
                cur_inflated_patch.get_vertex_property<Vector3f>("v:point");
            for (size_t i = 0; i < cur_inflated_patch.n_vertices(); i++) {
              Vertex vv(i);
              if (!is_end_sample_needed(vv))
                continue;
              inflat_samples.emplace_back(inflated_vpost[vv],
                                          inflated_orig_f_idx[vv]);
            }
          }

          // CDT each component
          Mesh comp_inflated_patch;
          cdt_non_SO_component(mesh, patch_comp, camera, inflat_samples,
                               comp_inflated_patch, false);
          cur_inflated_patch = comp_inflated_patch;
          cur_inflated_patch.m_fedges.clear();
          cur_inflated_patch.m_svertices.clear();
          prev_inflat_samples = inflat_samples;

          inflation_itr++;
        } while (inflation_itr < max_inflation_itr);

        // logger().info("Patch {}, Comp {}: #samples {}", decomp_patch_id,
        //               comp_idx, prev_inflat_samples.size());

        cdt_non_SO_component(mesh, patch_comp, camera, prev_inflat_samples,
                             comp_inflated_patch, true);
      } else {
        logger().info("Invalid patch {}, comp {}, skips sampling.",
                      decomp_patch_id, comp_idx);
        handle_degenerated_component_facing(mesh, patch_comp,
                                            comp_inflated_patch);
      }

      comp_inflated_patches.emplace_back(comp_inflated_patch);
      logger().info("\t=> #F {}", comp_inflated_patch.n_faces());
    }
  }

  // Stitch back
  auto inflat_patch_component =
      inflat_patch.add_face_property<int>("f:patch_component", -1);
  std::unordered_map<int, int> stitch_comp_b_v_to_new_v;
  for (auto const &comp_inflated_patch : comp_inflated_patches) {
    auto comp_boundary_idx =
        comp_inflated_patch.get_vertex_property<int>("v:comp_idx");
    contess_assert_msg(comp_boundary_idx,
                       "lift_patch_components: Index correspondence for "
                       "stitching is missing.");
    auto orig_boundary_idx =
        comp_inflated_patch.get_vertex_property<int>("v:orig_idx");
    auto orig_f_idx =
        comp_inflated_patch.get_vertex_property<int>("v:orig_f_idx");
    contess_assert_msg(orig_f_idx && orig_boundary_idx,
                       "lift_patch_components: Index correspondence for "
                       "stitching is missing.");
    auto vpositions =
        comp_inflated_patch.get_vertex_property<Vector3f>("v:point");

    // This one is used to add faces
    std::unordered_map<int, int> comp_v_to_new_v;

    for (size_t j = 0; j < comp_inflated_patch.n_vertices(); j++) {
      Vertex comp_v(j);

      // This vertex needs to be stitched
      if (comp_boundary_idx[comp_v] >= 0) {
        int comp_b_idx = comp_boundary_idx[comp_v];
        if (!stitch_comp_b_v_to_new_v.count(comp_b_idx)) {
          Vertex v = inflat_patch.add_vertex(vpositions[comp_v]);
          comp_v_to_new_v[j] = v.idx();
          stitch_comp_b_v_to_new_v[comp_b_idx] = v.idx();
          inflated_orig_idx[v] = orig_boundary_idx[comp_v];
          inflated_orig_f_idx[v] = orig_f_idx[comp_v];
        }
        comp_v_to_new_v[j] = stitch_comp_b_v_to_new_v[comp_b_idx];
      } else {
        Vertex v = inflat_patch.add_vertex(vpositions[comp_v]);
        comp_v_to_new_v[j] = v.idx();
        inflated_orig_idx[v] = orig_boundary_idx[comp_v];
        inflated_orig_f_idx[v] = orig_f_idx[comp_v];
      }
    }

    auto comp_inflat_patch_component =
        comp_inflated_patch.get_face_property<int>("f:patch_component");
    int comp_idx =
        comp_inflat_patch_component[*comp_inflated_patch.faces_begin()];
    for (size_t j = 0; j < comp_inflated_patch.n_faces(); j++) {
      Vertex vv[3];
      comp_inflated_patch.verticesOfFace(Face(j), vv);

      contess_assert_msg(
          comp_v_to_new_v.count(vv[0].idx()) &&
              comp_v_to_new_v.count(vv[1].idx()) &&
              comp_v_to_new_v.count(vv[2].idx()),
          "lift_patch_components: Face contains unseen vertices.");

      std::vector<Vertex> new_f_vv({Vertex(comp_v_to_new_v[vv[0].idx()]),
                                    Vertex(comp_v_to_new_v[vv[1].idx()]),
                                    Vertex(comp_v_to_new_v[vv[2].idx()])});
      Face f = inflat_patch.add_face(new_f_vv);
      if (f.is_valid()) {
        patchID[f] = decomp_patch_id;
        inflat_patch_component[f] = comp_idx;
      }
    }
  }

  // Add boundary information for visualization
  if (!inflat_patch.get_edge_property<bool>("e:is_boundary"))
    inflat_patch.add_edge_property<bool>("e:is_boundary", false);
  auto is_boundary = inflat_patch.edge_property<bool>("e:is_boundary");
  is_boundary.vector().assign(is_boundary.vector().size(), false);
  for (size_t j = 0; j < inflat_patch.n_edges(); j++) {
    Edge e(j);
    is_boundary[e] = inflat_patch.is_boundary(e);
  }
}
