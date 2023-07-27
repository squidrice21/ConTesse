// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "tag_simplification_edges.h"
#include <deque>
#include <limits>
#include <spdlog/fmt/ostr.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "common.h"
#include "cut_patch_to_disk.h" // for debug output
#include "subdivide_contour_edges_even.h"
#include "tag_cusp_facing.h"
#include "tag_simplification_edges_cases.h"

// Special case when walking on boundary,
// stops when arriving at a boundary-contour intersection point
bool is_boundary_joint(Mesh const &mesh, Vertex const &v) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");

  auto hit = mesh.halfedges(v), hit_end = hit;
  bool seen_boundary = false, seen_contour = false;
  do {
    if (mesh.is_boundary(mesh.edge(*hit)))
      seen_boundary = true;
    if (is_contour[mesh.edge(*hit)] >= 0)
      seen_contour = true;
    if (seen_boundary && seen_contour)
      return true;
  } while (++hit != hit_end);
  return false;
}

void tag_unmoveable_vertices(Mesh &mesh, Camera const &camera) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto intersection_2d = mesh.get_vertex_property<int>("v:intersection_2d");
  contess_assert_msg(
      intersection_2d,
      "tag_unmoveable_vertices: 2D intersections are not created.");
  auto is_valid_intersection_2d =
      mesh.get_vertex_property<bool>("v:is_valid_intersection_2d");
  contess_assert_msg(
      is_valid_intersection_2d,
      "tag_unmoveable_vertices: Valid 2D intersections are not tagged.");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto is_cusp = mesh.get_vertex_property<bool>("v:cusp");
  auto cusp_facing = mesh.get_vertex_property<FacingType>("v:cusp_facing");
  contess_assert_msg(is_contour && is_cusp && cusp_facing,
                     "Requires contours and cusps.");

  auto feasibility_collapsing =
      mesh.get_edge_property<int>("e:feasibility_collapsing");
  auto stationary = mesh.vertex_property<bool>("v:stationary");

  auto considered_a_cusp = [&](Vertex const &to_v) -> bool {
    return is_cusp[to_v] || cusp_facing[to_v] != FacingType::UNDEFINED;
  };

  // For each valid 2D intersection vertex,
  // walk and tag vertex that can't be moved
  std::unordered_set<int> seen_vertices;
  for (size_t i = 0; i < mesh.n_vertices(); i++) {
    Vertex v(i);

    if (!is_valid_intersection_2d[v] || seen_vertices.count(v.idx()) ||
        stationary[v])
      continue;
    seen_vertices.emplace(v.idx());
    seen_vertices.emplace(intersection_2d[v]);

    // Walk and record cusp and facing
    std::unordered_map<int, FacingType> seen_cusps;
    Halfedge he;
    auto hit = mesh.halfedges(v), hit_end = hit;
    do {
      if (feasibility_collapsing[mesh.edge(*hit)] >= 0) {
        he = *hit;
        break;
      }
    } while (++hit != hit_end);

    // This vertex was determined a false positive
    if (!he.is_valid())
      continue;

    do {
      Vertex to_v = mesh.to_vertex(he);
      if (considered_a_cusp(to_v) || is_boundary_joint(mesh, to_v)) {
        seen_cusps[to_v.idx()] = cusp_facing[to_v];
      }

      hit = mesh.halfedges(to_v), hit_end = hit;
      Halfedge next_he;
      do {
        if (feasibility_collapsing[mesh.edge(*hit)] >= 0 &&
            mesh.to_vertex(*hit) != mesh.from_vertex(he)) {
          next_he = *hit;
          break;
        }
      } while (++hit != hit_end);

      if (!next_he.is_valid())
        break;

      he = next_he;
    } while (he.is_valid() && !is_valid_intersection_2d[mesh.to_vertex(he)]);

    // We have a cusp at one end
    Vertex to_v = mesh.to_vertex(he);
    if (!is_valid_intersection_2d[to_v] &&
        (considered_a_cusp(to_v) || is_boundary_joint(mesh, to_v))) {
      stationary[to_v] = true;
      continue;
    }

    // Determine the cusp to keep
    std::vector<int> facing_count;
    facing_count.resize(2, 0);
    for (auto const &cusp : seen_cusps) {
      if (cusp.second == FacingType::FRONT) {
        facing_count[0]++;
      } else if (cusp.second == FacingType::BACK) {
        facing_count[1]++;
      }
    }

    // Find the furthest point
    Vector2f v_pos = project(vpositions[v], camera.viewMatrix().matrix(),
                             camera.projectionMatrix(),
                             Vector2i(camera.vpWidth(), camera.vpHeight()))
                         .head(2);
    FacingType keep_facing = (facing_count[0] > facing_count[1])
                                 ? FacingType::FRONT
                                 : FacingType::BACK;
    if (facing_count[0] == facing_count[1])
      keep_facing = FacingType::UNDEFINED;
    real_t max_dist = -1;
    int cusp_v = -1;
    for (auto const &cusp : seen_cusps) {
      // Label boundary-contour joints here
      if (cusp.second == FacingType::UNDEFINED) {
        stationary[Vertex(cusp.first)] = true;
        continue;
      }
      if (cusp.second != keep_facing)
        continue;
      Vector2f pos =
          project(vpositions[Vertex(cusp.first)], camera.viewMatrix().matrix(),
                  camera.projectionMatrix(),
                  Vector2i(camera.vpWidth(), camera.vpHeight()))
              .head(2);
      real_t d = (pos - v_pos).norm();
      if (d > max_dist) {
        max_dist = d;
        cusp_v = cusp.first;
      }
    }
    if (cusp_v >= 0)
      stationary[Vertex(cusp_v)] = true;
  }
}

bool is_same_patch(Mesh const &mesh, Halfedge he, Halfedge pair_he,
                   Halfedge hh) {
  auto patchBoundary = mesh.get_halfedge_property<int>("h:patchBoundary");
  contess_assert_msg(patchBoundary, "is_same_patch: Patch needs to be built.");
  bool need_to_flip_side = (patchBoundary[he] < 0);
  auto in_he = he;
  he = (need_to_flip_side) ? mesh.opposite_halfedge(he) : he;

  if (pair_he.is_valid()) {
    std::unordered_set<int> pair_patches;
    if (mesh.face(pair_he).is_valid())
      pair_patches.emplace(patchBoundary[pair_he]);
    if (mesh.face(mesh.opposite_halfedge(pair_he)).is_valid())
      pair_patches.emplace(patchBoundary[mesh.opposite_halfedge(pair_he)]);
    if (!pair_patches.count(patchBoundary[he])) {
      he = mesh.opposite_halfedge(he);
    }
    contess_assert_msg(patchBoundary[he] >= 0 &&
                           pair_patches.count(patchBoundary[he]),
                       "is_same_patch: Cannot find the patch corresponding to "
                       "both paired hes.");
  }

  need_to_flip_side = (he != in_he);

  int patch_id = patchBoundary[he];
  int hh_patch_id = (need_to_flip_side)
                        ? patchBoundary[mesh.opposite_halfedge(hh)]
                        : patchBoundary[hh];
  return patch_id == hh_patch_id;
}

Halfedge find_nondegenerated_contour_halfedge(Mesh const &mesh, Vertex const &v,
                                              Halfedge const &forward) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  std::unordered_set<int> visited_v;
  auto hit = mesh.halfedges(v), hit_end = hit;
  do {
    if ((mesh.is_boundary(mesh.edge(*hit)) ||
         is_contour[mesh.edge(*hit)] >= 0) &&
        *hit != forward) {
      visited_v.emplace(mesh.to_vertex(*hit).idx());
      break;
    }
  } while (++hit != hit_end);

  visited_v.emplace(v.idx());
  Halfedge he_long_forward =
      find_nondegenerated_contour_halfedge(mesh, v, visited_v);
  return he_long_forward;
}

bool tag_simplification_edges(Mesh &mesh, Camera const &camera) {
  auto feasibility_collapsing =
      mesh.get_edge_property<int>("e:feasibility_collapsing");
  if (!feasibility_collapsing) {
    feasibility_collapsing =
        mesh.add_edge_property<int>("e:feasibility_collapsing", -1);
  }
  feasibility_collapsing.vector().assign(feasibility_collapsing.vector().size(),
                                         -1);
  auto stationary = mesh.get_vertex_property<bool>("v:stationary");
  if (!stationary) {
    stationary = mesh.add_vertex_property<bool>("v:stationary", false);
  }
  stationary.vector().assign(stationary.vector().size(), false);

  bool successful = true;

  std::unordered_map<int, int> cusp_projections;
  std::unordered_map<int, Vector3f> moved_proj_positions;

  // For each valid 2D intersection vertex,
  // walk and tag edges that need to be collapsed
  size_t case1_count = 0, case2_count = 0, case3_count = 0;
  for (size_t v_i = 0; v_i < mesh.n_vertices(); v_i++) {
    auto is_cusp = mesh.get_vertex_property<bool>("v:cusp");
    auto intersection_2d = mesh.get_vertex_property<int>("v:intersection_2d");
    contess_assert_msg(
        intersection_2d,
        "tag_simplification_edges: 2D intersections are not created.");
    auto is_valid_intersection_2d =
        mesh.get_vertex_property<bool>("v:is_valid_intersection_2d");
    contess_assert_msg(
        is_valid_intersection_2d,
        "tag_simplification_edges: Valid 2D intersections are not tagged.");

    Vertex v(v_i);
    Vertex inter_v(intersection_2d[v]);

    // The later condition avoid running on created intersection with flipped
    // cusp
    if (!is_valid_intersection_2d[v] ||
        (is_cusp[v] || is_cusp[Vertex(intersection_2d[v])]))
      continue;

    // We run on a pair together
    if (inter_v.idx() < v.idx())
      continue;

    // Default: Case 1&2. Tracing the determined directions.
    Mesh mesh_labeled_case12_def;
    real_t walk_length1_def, chain_walk_length1_def, walk_length2_def,
        chain_walk_length2_def;
    size_t case1_count_def, case2_count_def;
    successful = tag_simplification_edges_case1(
        mesh, camera, v, mesh_labeled_case12_def, walk_length1_def,
        chain_walk_length1_def, walk_length2_def, chain_walk_length2_def,
        case1_count_def, case2_count_def);
    if (!successful)
      return false;

    // Alternative: Case 1&2. Tracing the reversed directions.
    Mesh mesh_labeled_case12_alt;
    real_t walk_length1_alt, chain_walk_length1_alt, walk_length2_alt,
        chain_walk_length2_alt;
    size_t case1_count_alt, case2_count_alt;
    successful = tag_simplification_edges_case1_alt(
        mesh, camera, v, mesh_labeled_case12_alt, walk_length1_alt,
        chain_walk_length1_alt, walk_length2_alt, chain_walk_length2_alt,
        case1_count_alt, case2_count_alt);
    if (!successful)
      return false;

    // Case 3
    Mesh mesh_labeled_case3;
    real_t furthest_walk_length_min, distance;
    bool succeed = tag_simplification_edges_case3(
        mesh, camera, v, mesh_labeled_case3, cusp_projections,
        moved_proj_positions, furthest_walk_length_min, distance);

    // Determine which case it is
    auto is_walk_safe = [&](real_t walk_length,
                            real_t chain_walk_length) -> bool {
      return !(walk_length / chain_walk_length > 0.9 || walk_length > 100);
    };

    // Case 1&2
    if (is_walk_safe(walk_length1_def, chain_walk_length1_def) &&
        is_walk_safe(walk_length2_def, chain_walk_length2_def)) {
      logger().info("walk_length: {} - {}, {} / {} = {}; {} / {} = {}", v,
                    inter_v, walk_length1_def, chain_walk_length1_def,
                    walk_length1_def / chain_walk_length1_def, walk_length2_def,
                    chain_walk_length2_def,
                    walk_length2_def / chain_walk_length2_def);
      mesh = mesh_labeled_case12_def;
      case1_count += case1_count_def;
      case2_count += case2_count_def;
    } else {
      logger().info(
          "walk_length alternative: {} - {}, {} / {} = {}; {} / {} = {} VS {}",
          v, inter_v, walk_length1_alt, chain_walk_length1_alt,
          walk_length1_alt / chain_walk_length1_alt, walk_length2_alt,
          chain_walk_length2_alt, walk_length2_alt / chain_walk_length2_alt,
          furthest_walk_length_min);

      // Case 3
      if ((std::max(walk_length1_alt, walk_length2_alt) /
                   furthest_walk_length_min >
               4 ||
           furthest_walk_length_min < 2) &&
          furthest_walk_length_min < 40 &&
          distance / chain_walk_length1_alt < 0.5) {
        logger().info("\t=> Cusp flipping.");
        if (!succeed) {
          logger().info("\t=> Cusp projection failed.");
          continue;
        } else {
          mesh = mesh_labeled_case3;
          case3_count++;
        }
      }
      // Case 1 & 2, alt
      else if (is_walk_safe(walk_length1_alt, chain_walk_length1_alt) &&
               is_walk_safe(walk_length2_alt, chain_walk_length2_alt)) {
        logger().info("\tDistance {} / {}.", distance, chain_walk_length1_alt);
        mesh = mesh_labeled_case12_alt;
        case1_count += case1_count_alt;
        case2_count += case2_count_alt;
      }
    }
  }

  logger().info("#Case 3: {}", case1_count);
  logger().info("#Case 4: {}", case2_count);
  logger().info("#Case 5: {}", case3_count);

  // Tag vertices that we are not going to move
  tag_unmoveable_vertices(mesh, camera);

  // Move the cusp projections
  auto vpositions = mesh.vertex_property<Vector3f>("v:point");
  for (auto const &proj : moved_proj_positions) {
    Vertex v(proj.first);
    vpositions[v] = proj.second;
  }

  return successful;
}
