// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.


// Catch unit test framework
#include <catch2/catch.hpp>
#include <csignal>

// Eigen printout formatting
#include <cmath>
#include <limits>
#include <spdlog/fmt/ostr.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <spdlog/common.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_sinks.h>

#include "camera.h"
#include "common.h"
#include "foldertools.h"
#include "io/serialization.h"
#include "logger.h"
#include "mesh.h"
#include "svg.h"
#include "testing.h"
#include "testing_cases.h"

#include "chain_contour.h"
#include "collapse_flipped_edges.h"
#include "cut_patch_to_disk.h"
#include "evaluate_radial_curvature.h"
#include "fix_flipped_faces.h"
#include "insert_interpolated_contours.h"
#include "insert_planar_map_intersections.h"
#include "insert_subdivision_cusps.h"
#include "make_cut_feasible.h"
#include "remove_cut_infeasible_disk.h"
#include "remove_wso_failure.h"
#include "shrink_cut.h"
#include "subdivide_contour_edges.h"
#include "tag_2d_intersections.h"
#include "tag_concave_edges.h"
#include "tag_cusp_facing.h"
#include "tag_invalid_cut_edges.h"
#include "tag_simplification_edges.h"
#include "tag_simplification_edges_small_loops.h"

#include "inflat_path.h"
#include "lift_min_surface.h"
#include "shor.h"
#include "triangulate_wso_patches.h"

#include "cdt_non_SO_component.h"
#include "decompose_wso_triangulation.h"
#include "label_patch_boundary.h"
#include "lift_patch_components.h"
#include "sew_components.h"

#include "ray_cast_QI.h"
#include "stitch_mesh.h"

namespace {

struct WSOParameters {
  size_t contour_neighbor_size = 5;
  size_t contour_neighbor_size_1d = 5;
  bool use_simple_condition = false;
};
bool dump_final_result = false;
bool export_blender = false;
size_t max_subdiv_level = 5;

bool ff_only = false;
bool use_heuristics = true;

real_t max_loop_length = 100;
// real_t max_loop_length = 200;
} // namespace

bool run_WSO(std::string const output_dir, const TestCase &t_case,
             WSOParameters const &wso_param, Mesh &mesh, Camera &camera,
             std::map<size_t, Mesh> &patch_triangulations) {
  auto contour_neighbor_size = wso_param.contour_neighbor_size;
  auto contour_neighbor_size_1d = wso_param.contour_neighbor_size_1d;
  auto use_simple_condition = wso_param.use_simple_condition;

  std::string filename;
  if (export_blender) {
    camera.getInterpolator().setCamera(&camera);
    camera.getInterpolator().addKeyFrame(0);
    camera.save(output_dir + "/" + t_case.model_name + ".cam");
    filename = output_dir + "/" + t_case.model_name + "_subdiv.obj";
    mesh.write(filename);
    logger().info("Blender export. Output: " + filename);
    return false;
  }

  bool to_insert_cusp = true;
  logger().info("insert_contours, to_insert_cusp: {}", to_insert_cusp);
  bool successful =
      insert_contours(mesh, camera, ContourMode::VBO_CONTOUR, to_insert_cusp);
  if (!successful)
    return false;

  logger().info("tag_concave_edges_analytic");
  tag_concave_edges_analytic(mesh, camera);
  logger().info("tag_concave_edges_analytic done.");

  logger().info("tag_cusp_given_concave_edges");
  tag_cusp_given_concave_edges(mesh, camera);
  logger().info("tag_cusp_given_concave_edges done.");

  logger().info("tag_cusp_sharp");
  tag_cusp_sharp(mesh, camera);
  logger().info("tag_cusp_sharp done.");

  logger().info("compute_cusp_laplacian");
  bool not_degenerated = compute_cusp_laplacian(mesh);
  if (!not_degenerated)
    return false;
  filename = output_dir + "/" + t_case.model_name + "_laplacian.json";
  logger().info("compute_cusp_laplacian done. Output: " + filename);

  logger().info("tag_cusp_facing");
  successful = tag_cusp_facing(mesh, camera);
  logger().info("tag_cusp_facing done.");
  if (!successful)
    return false;

  logger().info("tag_extraordinary");
  tag_extraordinary(mesh);
  filename = output_dir + "/" + t_case.model_name + "_cusp.json";
  logger().info("tag_extraordinary done. Output: " + filename);

  logger().info("insert_planar_map_intersections");
  insert_planar_map_intersections(mesh, camera);

  // Tag the candidate cut edges for the following operations
  logger().info("tag_invalid_cut_edges");
  successful =
      tag_invalid_cut_edges(mesh, camera, contour_neighbor_size,
                            contour_neighbor_size_1d, use_simple_condition);
  if (!successful)
    return false;

  // Serialize the mesh and camera
  filename = output_dir + "/" + t_case.model_name + "_cutcand.json";
  logger().info("tag_invalid_cut_edges done. Output: " + filename);

  // Cut the patch
  logger().info("cut_patch_to_disk with modified Gu02");
  cut_patch_to_disk(mesh, camera, false, ff_only);
  logger().info("cut_patch_to_disk with modified Gu02 done.");

  logger().info("chain_cut_graph");
  successful = chain_cut_graph(mesh, camera, ff_only);
  logger().info("chain_cut_graph done.");
  if (!successful)
    return false;

  /////////////
  if (use_heuristics) {
    // Find valid 2D intersections
    logger().info("tag_2d_intersections");
    tag_2d_intersections(mesh, camera);
    logger().info("tag_2d_intersections done.");

    // Remove the remaining non fishtails
    logger().info("tag_non_fish_tail");
    tag_non_fish_tail(mesh, camera, ff_only);
    logger().info("tag_non_fish_tail done.");

    // Find edge chains to collapse
    logger().info("tag_simplification_edges");
    bool successful = tag_simplification_edges(mesh, camera);
    logger().info("tag_simplification_edges done.");
    if (!successful)
      return false;

    // Find additional edge chains to collapse (involving non-fishtails)
    {
      real_t loop_epsilon = 1e-3;
      logger().info("tag_simplification_edges_flipped_small_loops");
      tag_simplification_edges_flipped_small_loops(mesh, camera, loop_epsilon,
                                                   ff_only);
      logger().info("tag_simplification_edges_flipped_small_loops done.");
    }

    // Collapse
    real_t small_epsilon = 1e-6;
    logger().info("collapse_flipped_edges");
    successful = collapse_flipped_edges(mesh, camera, small_epsilon);
    logger().info("collapse_flipped_edges done.");
    if (!successful)
      return false;

    //////
    {
      auto feasibility_collapsing =
          mesh.edge_property<int>("e:feasibility_collapsing");
      feasibility_collapsing.vector().assign(
          feasibility_collapsing.vector().size(), -1);
      real_t loop_epsilon = 1e-3;

      logger().info("tag_simplification_edges_flipped_small_loops 2");
      successful = tag_simplification_edges_flipped_small_loops(
          mesh, camera, loop_epsilon, ff_only);
      logger().info("tag_simplification_edges_flipped_small_loops 2 done.");
      if (!successful)
        return false;
      logger().info("collapse_flipped_edges");
      successful = collapse_flipped_edges(mesh, camera, small_epsilon);
      logger().info("collapse_flipped_edges done.");
      if (!successful)
        return false;
    }
    //////
    mesh.get_patch_chains().clear();
    mesh.get_oriented_chains().clear();

    // Tag the candidate cut edges for the following operations
    logger().info("tag_invalid_cut_edges");
    tag_invalid_cut_edges(mesh, camera, contour_neighbor_size,
                          contour_neighbor_size_1d, use_simple_condition);
    logger().info("tag_invalid_cut_edges done.");

    // Redo chaining using the original orientation.
    logger().info("make_cut_feasible");
    make_cut_feasible(mesh, camera, ff_only);
    logger().info("make_cut_feasible done.");

    // Tag the candidate cut edges for the following operations
    logger().info("tag_invalid_cut_edges");
    tag_invalid_cut_edges(mesh, camera, contour_neighbor_size,
                          contour_neighbor_size_1d, use_simple_condition);
    logger().info("tag_invalid_cut_edges done.");
    // Cut the patch
    logger().info("cut_patch_to_disk with modified Gu02: second time");
    cut_patch_to_disk(mesh, camera, false, ff_only);

    bool removed = remove_cut_infeasible_disk(mesh, camera);
    if (removed) {
      // Cut the patch
      logger().info("cut_patch_to_disk with modified Gu02: third time");
      cut_patch_to_disk(mesh, camera, false, ff_only);
    }

    logger().info("shrink_cut");
    shrink_cut(mesh, camera);

    mesh.get_patch_chains().clear();
    mesh.get_oriented_chains().clear();
    logger().info("chain_cut_graph");
    successful = chain_cut_graph(mesh, camera, ff_only);
    logger().info("chain_cut_graph done.");
    if (!successful)
      return false;
  }
  /////////////

  // Run WSO
  bool wso_adaptive = false;
  logger().info("triangulate_wso_patches - wso_adaptive: {}", wso_adaptive);
  bool do_refine = false;
  triangulate_wso_patches(mesh, camera, patch_triangulations, do_refine,
                          std::unordered_set<int>(), ff_only);

  size_t successful_patch_count = 0;
  for (auto &patch : patch_triangulations) {
    // Check if the previous WSO triangulations succeeded
    // Delay throw to later so we can run to the end on successful patches.
    if (patch.second.n_faces() > 0) {
      successful_patch_count++;

      // Copy patch ID
      auto patchID = patch.second.get_face_property<int>("f:patchID");
      if (!patchID)
        patchID = patch.second.add_face_property<int>("f:patchID", -1);
      else
        patchID = patch.second.face_property<int>("f:patchID");

      for (auto f_itr = patch.second.faces_begin();
           f_itr != patch.second.faces_end(); ++f_itr) {
        Face f = *f_itr;
        patchID[f] = patch.first;
      }
    }
  }
  logger().info("triangulate_wso_patches done. First round: {} / {}",
                successful_patch_count, patch_triangulations.size());
  INFO("triangulate_wso_patches done. First round: "
       << successful_patch_count << " / " << patch_triangulations.size());

  // Second round WSO
  if (successful_patch_count < patch_triangulations.size()) {
    std::unordered_set<int> wso_failures;
    for (auto &patch : patch_triangulations) {
      if (patch.second.n_faces() == 0)
        wso_failures.emplace(patch.first);
    }
    std::vector<std::vector<int>> patch_removals;
    logger().info("get_wso_failure_removal: {}", max_loop_length);
    get_wso_failure_removal(mesh, camera, wso_failures, patch_removals,
                            max_loop_length, ff_only);
    for (size_t i = 0; i < patch_removals.size(); i++) {
      if (patch_removals[i].empty())
        continue;

      logger().info("remove_wso_failure");
      std::unordered_set<int> affected_patches;
      bool removed = remove_wso_failure(mesh, camera, patch_removals[i],
                                        affected_patches, ff_only);
      if (!removed) {
        for (auto pid : wso_failures)
          affected_patches.emplace(pid);
        logger().info("remove_wso_failure done.");
        continue;
      }

      // Cut the patch
      logger().info("cut_patch_to_disk with modified Gu02: 3+{}th time", i + 1);
      cut_patch_to_disk(mesh, camera, false, ff_only);

      logger().info("shrink_cut");
      shrink_cut(mesh, camera);

      mesh.get_patch_chains().clear();
      mesh.get_oriented_chains().clear();
      logger().info("chain_cut_graph");
      successful = chain_cut_graph(mesh, camera, ff_only);
      logger().info("chain_cut_graph done.");
      if (!successful)
        return false;

      logger().info("WSO {}th time", i + 2);
      std::map<size_t, Mesh> patch_triangulations2;
      triangulate_wso_patches(mesh, camera, patch_triangulations2, do_refine,
                              affected_patches, ff_only);

      // Combine patches
      std::map<size_t, Mesh> patch_triangulations_comb;
      for (auto const &p : patch_triangulations) {
        if (std::find(patch_removals[i].begin(), patch_removals[i].end(),
                      p.first) == patch_removals[i].end() &&
            !affected_patches.count(p.first)) {
          logger().info("\tCopy: p{}, #f: {}", p.first, p.second.n_faces());
          patch_triangulations_comb[p.first] = p.second;
        }
      }
      for (auto const &p : patch_triangulations2) {
        logger().info("\tWSO {}th time: p{}, #f: {}", i + 2, p.first,
                      p.second.n_faces());
        patch_triangulations_comb[p.first] = p.second;
      }
      patch_triangulations = patch_triangulations_comb;

      // Check failures
      size_t failure_count = 0;
      for (auto const &patch : patch_triangulations) {
        if (patch.second.n_faces() == 0)
          failure_count++;
      }

      logger().info("WSO {}th time done. Failure: {} / {}", i + 2,
                    failure_count, patch_triangulations.size());
      if (!failure_count) {
        if (i == 0) {
          logger().info("# Aggressive");
        } else {
          logger().info("# Conservative");
        }
        logger().info("#Case 1: {}", patch_removals[i].size());
        break;
      }
    }
  }
  {
    successful_patch_count = 0;
    for (auto const &patch : patch_triangulations) {
      // Check if the previous WSO triangulations succeeded
      // Delay throw to later so we can run to the end on successful patches.
      if (patch.second.n_faces() > 0)
        successful_patch_count++;
    }
    logger().info("triangulate_wso_patches done. Success: {} / {}",
                  successful_patch_count, patch_triangulations.size());
    INFO("triangulate_wso_patches done. Success: "
         << successful_patch_count << " / " << patch_triangulations.size());
  }

  return successful_patch_count == patch_triangulations.size();
}

void mesh_statistics(const Mesh &mesh) {
  // 1. Open vs closed
  bool is_closed = true;
  for (size_t i = 0; i < mesh.n_edges(); ++i) {
    if (mesh.is_boundary(Edge(i))) {
      is_closed = false;
      break;
    }
  }

  logger().info("Input mesh: {}", (is_closed) ? "closed" : "open");

  // 2. Input face count
  logger().info("Input faces: {}", mesh.n_faces());
}

bool adjust_subdivision(std::string const output_dir, TestCase const &t_case,
                        WSOParameters const &wso_param, bool single_wso_trial,
                        Mesh &mesh, Camera &camera,
                        std::map<size_t, Mesh> &patch_triangulations) {
  bool wso_terminates = true;
  size_t wso_itr = 0;
  Mesh wso_mesh;
  bool wso_successful;
  TestCase wso_t_case = t_case;
  logger().info("## WSO round {} ##", wso_itr + 1);
  do {
    logger().info("## Subdivision level {} ##", wso_t_case.subdiv_level);
    wso_terminates = true;
    if (wso_t_case.subdiv_level != t_case.subdiv_level)
      wso_mesh = read_input(t_case.model_path + ".obj", false,
                            wso_t_case.subdiv_level, Subdiv::Backend::LACEWELL);

    Camera wso_camera = camera;

    if (wso_t_case.subdiv_level != t_case.subdiv_level) {
      // Run statistics on input mesh
      mesh_statistics(wso_mesh);
      wso_successful = run_WSO(output_dir, wso_t_case, wso_param, wso_mesh,
                               wso_camera, patch_triangulations);
    } else {
      // Run statistics on input mesh
      mesh_statistics(mesh);
      wso_successful = run_WSO(output_dir, wso_t_case, wso_param, mesh,
                               wso_camera, patch_triangulations);
    }

    if (!wso_successful) {
      wso_terminates = false;
      if (wso_itr == 0 && t_case.subdiv_level > 1) {
        wso_t_case.subdiv_level = t_case.subdiv_level - 1;
      } else if (wso_t_case.subdiv_level != (int)max_subdiv_level) {
        wso_t_case.subdiv_level++;
      } else {
        wso_terminates = true;
      }

      if (single_wso_trial) {
        wso_terminates = true;
        wso_t_case.subdiv_level = t_case.subdiv_level;
      }

      if (!wso_terminates) {
        logger().info("## WSO round {}: Subdiv: {} ##", wso_itr + 2,
                      wso_t_case.subdiv_level);
        patch_triangulations.clear();
      } else
        logger().error("WSO gave up.");
    }

    if (wso_terminates) {
      if (wso_t_case.subdiv_level != t_case.subdiv_level) {
        mesh = wso_mesh;
        mesh.m_fedges.clear();
        mesh.m_svertices.clear();
      }
      camera = wso_camera;
    }

    wso_itr++;
  } while (!wso_terminates && wso_itr < max_subdiv_level);

  return wso_successful;
}

bool run_lifting_rendering(std::string const output_dir, TestCase const &t_case,
                           const Camera &camera, const Mesh &mesh,
                           const std::map<size_t, Mesh> &patch_triangulations,
                           Mesh &stitched_mesh) {
  std::string filename;
  std::map<int, std::string> output_description;
  std::map<int, Mesh> output_meshes;
  for (auto const &patch : patch_triangulations) {
    int patch_id = patch.first;
    FacingType facing = mesh.get_patch_facing(patch_id);
    std::string facing_str = (facing == FacingType::FRONT) ? "f" : "b";
    std::string mesh_key = "p" + std::to_string(patch_id) + facing_str;
    output_meshes[patch_id] = patch.second;
    output_description[patch_id] = mesh_key;
  }

  bool early_terminate = false;
  {
    INFO("Verify WSO patch success: ");
    for (auto &patch : patch_triangulations) {
      std::string description = output_description[patch.first];
      INFO("Verifying patch: " + description);

      // Check if the previous WSO triangulations succeeded
      CHECK(patch.second.n_faces() > 0);
      early_terminate |= (patch.second.n_faces() == 0);
    }
  }

  if (early_terminate)
    return false;

  // Write all triangulated results
  std::vector<Mesh> stitched_patches;
  for (auto const &result : output_meshes) {
    if (result.second.faces_size() == 0) {
      continue;
    }
    std::string description = output_description[result.first];

    logger().info("----------------------------------------------");
    INFO("Processing patch: " + description + ".");

    // Lifted as json outputs
    Mesh flat_patch = result.second;
    flat_patch.m_fedges.clear();
    flat_patch.m_svertices.clear();

    // Decompose into non-self-overlapping components
    logger().info("label_patch_boundary");
    label_patch_boundary(mesh, flat_patch);
    logger().info("label_patch_boundary done.");

    logger().info("decompose_wso_triangulation");
    decompose_wso_triangulation(mesh, camera, flat_patch);

    // Lift the patch/manifold boundary to 3D using the correspondence to
    // the original mesh
    logger().info("lift_min_surface");
    Mesh boundary_lifted_patch;
    lift_min_surface(mesh, flat_patch, boundary_lifted_patch);

    // Lift the component boundaries before filling in interior samples
    // (within each component)
    Mesh comp_inflated_patch;
    std::unordered_set<BoundaryType> to_inflat_labels;
    to_inflat_labels.emplace(BoundaryType::COMPONENT);
    real_t sampling_delta =
        std::sqrt(mesh.boundingBox().diagonal().norm()) / 100;
    real_t min_delta_2d = 5;
    logger().info("inflat_path (Component boundary) - sampling_delta: {}; "
                  "min_delta_2d: {}",
                  sampling_delta, min_delta_2d);
    std::vector<std::pair<int, int>> failed_paths;
    inflat_path(mesh, boundary_lifted_patch, camera, sampling_delta,
                min_delta_2d, comp_inflated_patch, failed_paths,
                to_inflat_labels);

    // Fill in interior samples (within each component)
    stitched_patches.emplace_back();
    bool sparse_sampling = true;
    logger().info("inflat_path (Interior samples) - sparse_sampling: {}",
                  sparse_sampling);

    real_t max_distance = 0.001;
    size_t max_inflation_itr = 3;
    lift_patch_components_iterative(mesh, camera, sampling_delta, min_delta_2d,
                                    max_distance, boundary_lifted_patch,
                                    stitched_patches.back(), max_inflation_itr);
    logger().info("inflat_path (Interior samples) done. #faces: {}",
                  stitched_patches.back().n_faces());

    logger().info("sew_components: #paths: {}", failed_paths.size());
    sew_components(stitched_patches.back(), mesh, camera, failed_paths);
    logger().info("sew_components done.");

    logger().info("deduplicate_contour_faces");
    deduplicate_contour_faces(stitched_patches.back(), mesh, camera);

    // Automatic verification
    INFO("verify_inconsistency");
    bool consistent = verify_inconsistency(stitched_patches.back(), camera,
                                           mesh.get_patch_facing(result.first));
    if (!consistent) {
      return false;
    }
  }

  INFO("Verify WSO patch success: ");
  for (auto &patch : patch_triangulations) {
    std::string description = output_description[patch.first];
    INFO("Verifying patch: " + description);

    // Check if the previous WSO triangulations succeeded
    CHECK(patch.second.n_faces() > 0);
  }

  logger().info("==============================================");

  // Render SVGs
  std::vector<Mesh *> patches;
  // First add FF patches
  for (auto &patch : stitched_patches) {
    auto patchID = patch.get_face_property<int>("f:patchID");
    int pid = patchID[*patch.faces_begin()];
    FacingType facing = mesh.get_patch_facing(pid);
    if (facing == FacingType::FRONT)
      patches.emplace_back(&patch);
  }
  // Then BF patches
  for (auto &patch : stitched_patches) {
    auto patchID = patch.get_face_property<int>("f:patchID");
    int pid = patchID[*patch.faces_begin()];
    FacingType facing = mesh.get_patch_facing(pid);
    if (facing == FacingType::BACK)
      patches.emplace_back(&patch);
  }

  logger().info("stitch_mesh");
  stitch_mesh(mesh, camera, patches, stitched_mesh);

  logger().info("Contess remeshing done.");
  logger().info("Output faces: {}", stitched_mesh.n_faces());

  {
    stitched_mesh.updateBoundingBox();
    stitched_mesh.build_bvh();
    stitched_mesh.update_patch_facing(mesh);

    // Insert and match 2D intersections.
    insert_planar_map_intersections(stitched_mesh, camera, true);

    stitched_mesh.updateBoundingBox();
    stitched_mesh.build_bvh();
    logger().info("ray_cast_QI.");
    ray_cast_QI(stitched_mesh, camera, true);
  }

  logger().info("chain_rendered_contour.");
  chain_rendered_contour(camera, stitched_mesh);

  if (dump_final_result) {
    filename = output_dir + "/" + t_case.model_name + "_stitched_chain.json";
    serialize_mesh(stitched_mesh, camera,
                   std::vector<serialization::SerializationRay>(), filename);
    logger().info("chain_rendered_contour done.");
  }

  return true;
}

void run_pipeline(std::string const output_dir, TestCase const &in_t_case) {
  // bool single_wso_trial = true;
  bool single_wso_trial = false;
  // size_t num_pert = 4;
  size_t num_pert = 1;
  TestCase t_case = in_t_case;

  std::string log_name = t_case.get_short_name() + ".log";

  // Set up individual log files
  // logger().sinks().push_back(
  //     std::make_shared<spdlog::sinks::basic_file_sink_mt>(
  //         output_dir + "/" + log_name, true));
  // logger().sinks().back()->set_level(spdlog::level::info);

  std::string description = t_case.get_description();

  // Parameters
  WSOParameters wso_param;
  wso_param.contour_neighbor_size = 5;
  wso_param.contour_neighbor_size_1d = 5;
  wso_param.use_simple_condition = false;

  SECTION(description) {
    logger().info("Logging to {}", output_dir + "/" + log_name);
    std::string model_path =
        CONTOURS_TESSELATION_DATA_DIR "/meshes/" + t_case.model_path + ".obj";
    {
      std::ifstream f(model_path);
      if (!f.is_open())
        model_path =
            CONTOURS_TESSELATION_DATA_DIR "/" + t_case.model_path + ".obj";
    }
    logger().info("Reading {}", model_path);

    if (export_blender) {
      std::string filename =
          output_dir + "/" + t_case.model_name + "_control.obj";
      Mesh control_mesh;
      control_mesh.read(model_path);
      control_mesh.write(filename);
      logger().info("Blender export. Output: " + filename);
    }
    {
      Mesh control_mesh;
      control_mesh.read(model_path);
      logger().info("Control mesh faces: {}", control_mesh.n_faces());
    }
    Mesh mesh = read_input(t_case.model_path + ".obj", false,
                           t_case.subdiv_level, Subdiv::Backend::LACEWELL);

    // Set up the default camera
    Camera camera;
    Vector2i viewport({2048, 1600});
    camera.setViewport(viewport.x(), viewport.y());
    setup_camera(camera, mesh, t_case);

    camera.projectionMatrix();
    camera.viewMatrix();

    std::string filename = output_dir + "/" + t_case.model_name + "_init.json";

    std::map<size_t, Mesh> patch_triangulations;

    size_t pert_itr = 0;
    Camera pert_camera = camera;
    Mesh pert_mesh;
    // Mesh pert_mesh;
    bool wso_successful;
    real_t cam_offset_scale = 1e-2;
    real_t mesh_offset_scale = 1e-4;

    pert_mesh = mesh;
    pert_mesh.m_fedges.clear();
    pert_mesh.m_svertices.clear();
    do {
      logger().info("## Perturbation round {} ##", pert_itr + 1);

      wso_successful =
          adjust_subdivision(output_dir, t_case, wso_param, single_wso_trial,
                             pert_mesh, pert_camera, patch_triangulations);
      if (!wso_successful && pert_itr + 1 < num_pert) {
        pert_camera = perturb_camera(pert_itr, cam_offset_scale, camera);
        pert_mesh = mesh;
        pert_mesh.m_fedges.clear();
        pert_mesh.m_svertices.clear();
        perturb_mesh(pert_mesh, pert_itr, mesh_offset_scale);
      } else {
        camera = pert_camera;
        mesh = pert_mesh;
        pert_mesh.m_fedges.clear();
        pert_mesh.m_svertices.clear();
        if (!wso_successful)
          logger().error("Perturbation fails");
      }

      pert_itr++;

      //////////////////////////////
      // Call the following pipeline
      // Add more readable logger().information
      if (wso_successful) {
        Mesh stitched_mesh;

        wso_successful =
            run_lifting_rendering(output_dir, t_case, camera, mesh,
                                  patch_triangulations, stitched_mesh);
        if (wso_successful) {
          write_chained_contours(
              stitched_mesh, camera,
              output_dir + "/" + t_case.model_name + "_out.svg", 0);

          // Ray cast the trivial solution
          mesh.updateBoundingBox();
          mesh.build_bvh();
          insert_planar_map_intersections(mesh, camera);
          mesh.updateBoundingBox();
          mesh.build_bvh();
          ray_cast_QI(mesh, camera, true);

          write_contours(mesh, camera,
                         output_dir + "/" + t_case.model_name + "_trivial.svg",
                         0);

          logger().info("Verify QI changes.");
          verify_qi_changes(stitched_mesh);

          logger().info("==============================================");
        } else {
          pert_itr--;
          t_case.subdiv_level++;
          pert_mesh =
              read_input(t_case.model_path + ".obj", false, t_case.subdiv_level,
                         Subdiv::Backend::LACEWELL);
          patch_triangulations.clear();
        }
      }
    } while (!wso_successful && pert_itr < num_pert &&
             t_case.subdiv_level < (int)max_subdiv_level);
  }

  logger().sinks().pop_back();
}

/*================= Actual test cases ====================*/
TEST_CASE("basic_pipeline", "[pipeline]") {
  // Create a folder to put the results inside
  if (output_folder.empty()) {
    output_folder = test_output_dir();
    output_folder += "/basic_pipeline";
  }
  const std::string output_dir = output_folder;
  // foldertools::makedir(output_dir.c_str());

  for (auto t_case_path : TestingSingleton::instance().pipeline_test_cases) {
    TestCase t_case;
    if (!t_case_path.is_absolute()) {
      t_case_path = TestingSingleton::instance().spec_root_path / t_case_path;
    }

    serialization::load(t_case, t_case_path);

    run_pipeline(output_dir, t_case);
  }

  SUCCEED();
}

TEST_CASE("multicamera_pipeline", "[pipeline]") {
  // Create a folder to put the results inside
  if (output_folder.empty()) {
    output_folder = test_output_dir();
    output_folder += "/multicamera_pipeline";
  }
  const std::string output_dir = output_folder;
  // foldertools::makedir(output_dir.c_str());

  size_t camera_steps_minus_one = 2;
  double extension_ratio = 2;

  for (auto t_case_path : TestingSingleton::instance().multicamera_test_cases) {
    // Read camera filter if existing
    std::string t_case_path_str = t_case_path.string();
    size_t pos = 0;
    std::vector<std::string> t_case_subs;
    while ((pos = t_case_path_str.find(":")) != std::string::npos) {
      t_case_subs.emplace_back(t_case_path_str.substr(0, pos));
      t_case_path_str.erase(0, pos + 1);
    }
    t_case_subs.emplace_back(t_case_path_str);

    CHECK(!t_case_subs.empty());

    t_case_path = std::filesystem::path(t_case_subs.front());
    t_case_subs.erase(t_case_subs.begin());
    std::unordered_set<size_t> camera_digits;
    for (auto digit : t_case_subs) {
      camera_digits.emplace(std::stoul(digit));
    }

    TestCase t_case;
    if (!t_case_path.is_absolute()) {
      t_case_path = TestingSingleton::instance().spec_root_path / t_case_path;
    }
    serialization::load(t_case, t_case_path);

    // Set up multiple camera locations
    // 0. Check is there exists a front version
    std::string orig_path = t_case.model_path;
    {
      if (t_case.model_path.find("/") != std::string::npos) {
        t_case.model_path.erase(0, t_case.model_path.find("/"));
      }
      t_case.model_path = "front_mesh/" + t_case.model_path;

      std::string input_filename =
          CONTOURS_TESSELATION_DATA_DIR "/meshes/" + t_case.model_path + ".obj";
      {
        std::ifstream f(input_filename);
        if (!f.is_open())
          input_filename =
              CONTOURS_TESSELATION_DATA_DIR "/" + t_case.model_path + ".obj";
      }
      std::ifstream f(input_filename);
      if (!f.is_open())
        t_case.model_path = orig_path;
    }

    // 1. Build and extend the bbox
    Mesh mesh = read_input(t_case.model_path + ".obj", false, -1,
                           Subdiv::Backend::LACEWELL);
    AABB aabb = mesh.boundingBox();
    Vector3f center = aabb.min + 0.5 * (aabb.max - aabb.min);
    Vector3f diag = (aabb.max - aabb.min);

    // Make the AABB into a cube
    Vector3f square_diag = Vector3f(1, 1, 1);
    aabb.min = center - 0.5 * square_diag;
    aabb.max = center + 0.5 * square_diag;
    Vector3f diag_step = (1. / camera_steps_minus_one) * (aabb.max - aabb.min);

    std::map<int, Vector3f> camera_angles;

    for (size_t i = 0; i <= camera_steps_minus_one; i++) {
      for (size_t j = 0; j <= camera_steps_minus_one; j++) {
        for (size_t k = 0; k <= camera_steps_minus_one; k++) {
          if (!camera_digits.empty() &&
              camera_digits.count(i * 100 + j * 10 + k) == 0)
            continue;

          // Skip the model center
          if (std::abs(i - double(camera_steps_minus_one) / 2) <
                  std::numeric_limits<double>::epsilon() &&
              std::abs(j - double(camera_steps_minus_one) / 2) <
                  std::numeric_limits<double>::epsilon() &&
              std::abs(k - double(camera_steps_minus_one) / 2) <
                  std::numeric_limits<double>::epsilon())
            continue;

          // Set camera
          Vector3f cam_pos = aabb.min + i * diag_step.x() * Vector3f(1, 0, 0) +
                             j * diag_step.y() * Vector3f(0, 1, 0) +
                             k * diag_step.z() * Vector3f(0, 0, 1);
          // Move to the sphere
          Vector3f mesh_center = mesh.get_mesh_center();
          cam_pos = (extension_ratio * 0.5) * diag.norm() *
                        (cam_pos - mesh_center).normalized() +
                    mesh_center;
          camera_angles[i * 100 + j * 10 + k] =
              (cam_pos - mesh_center).normalized();

          BasicCameraCase camera_case(cam_pos.x(), cam_pos.y(), cam_pos.z(),
                                      mesh_center);
          TestCase cam_case = t_case;
          cam_case.camera_info = camera_case;
          cam_case.model_name += "_cam" + std::to_string(i) +
                                 std::to_string(j) + std::to_string(k);

          run_pipeline(output_dir, cam_case);
        }
      }
    }
  }

  SUCCEED();
}