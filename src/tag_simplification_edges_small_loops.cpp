// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "tag_simplification_edges_small_loops.h"

#include "tag_non_fish_tail.h"

#include <spdlog/fmt/ostr.h>
#include <unordered_set>
#include <vector>

bool tag_simplification_edges_flipped_small_loops(Mesh &mesh,
                                                  Camera const &camera,
                                                  real_t loop_epsilon,
                                                  bool ff_only) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto vnormals = mesh.get_vertex_property<Vector3f>("v:normal");
  auto is_cusp = mesh.get_vertex_property<bool>("v:cusp");
  auto intersection_2d = mesh.get_vertex_property<int>("v:intersection_2d");
  contess_assert_msg(intersection_2d,
                     "tag_simplification_edges_flipped_small_loops: "
                     "2D intersections are not created.");
  auto is_valid_intersection_2d =
      mesh.get_vertex_property<bool>("v:is_valid_intersection_2d");
  contess_assert_msg(is_valid_intersection_2d,
                     "tag_simplification_edges_flipped_small_loops: Valid 2D "
                     "intersections are not tagged.");

  auto feasibility_collapsing =
      mesh.get_edge_property<int>("e:feasibility_collapsing");
  auto stationary = mesh.get_vertex_property<bool>("v:stationary");
  contess_assert_msg(
      feasibility_collapsing && stationary,
      "tag_simplification_edges_flipped_small_loops: This function must "
      "be called after tag_simplification_edges.");
  feasibility_collapsing = mesh.edge_property<int>("e:feasibility_collapsing");
  stationary = mesh.vertex_property<bool>("v:stationary");

  // Update contour chains
  {
    mesh.get_patch_chains().clear();
    mesh.get_oriented_chains().clear();
    bool successful = chain_contour(mesh, camera, true, ff_only);
    if (!successful)
      return false;
  }

  std::unordered_set<int> seen_vertices;

  size_t case2_count = 0;
  for (auto const &chain : mesh.get_oriented_chains()) {
    // 1. Get all self intersections of the loop
    std::vector<Vertex> chain_vertices;
    chain_vertices.reserve(chain->size());
    std::map<int, int> vid2idx;

    // 2. Precompute the loop length
    std::vector<real_t> chain_length;
    chain_length.reserve(chain->size() + 1);
    chain_length.emplace_back(0);
    for (auto const &he : *chain) {
      vid2idx[mesh.from_vertex(he).idx()] = chain_vertices.size();
      chain_vertices.emplace_back(mesh.from_vertex(he));

      chain_length.emplace_back(
          chain_length.back() +
          (vpositions[mesh.from_vertex(he)] - vpositions[mesh.to_vertex(he)])
              .norm());
    }

    for (auto const &he : *chain) {
      Vertex v = mesh.from_vertex(he);
      // Valid intersections should be handled separatedly
      if (intersection_2d[v] < 0 || is_valid_intersection_2d[v] ||
          seen_vertices.count(v.idx()))
        continue;

      int v1 = v.idx();
      int v2 = intersection_2d[v];

      // Intersecting with a different loop
      if (!vid2idx.count(v2))
        continue;

      seen_vertices.emplace(v1);
      seen_vertices.emplace(v2);

      // Determine the smaller loop
      // From first edge to the second
      real_t v1_len = chain_length[vid2idx[v1]];
      real_t v2_len = chain_length[vid2idx[v2]];

      real_t v1_to_v2_len =
          ((v2_len < v1_len) ? v2_len + chain_length.back() : v2_len) - v1_len;
      real_t v2_to_v1_len =
          ((v1_len < v2_len) ? v1_len + chain_length.back() : v1_len) - v2_len;

      size_t start_i = vid2idx[v2];
      size_t end_i = vid2idx[v1];

      if (v1_to_v2_len < v2_to_v1_len) {
        start_i = vid2idx[v1];
        end_i = vid2idx[v2];
      }

      // Walk on halfedges
      // Count turning points
      Vector3f int_v1 = vpositions[Vertex(v1)];
      Vector3f int_v2 = vpositions[Vertex(v2)];
      real_t h_dist = -1;

      std::vector<std::vector<Vertex>> seen_turns;
      seen_turns.resize(2);

      std::vector<Edge> small_loop;
      for (size_t i = start_i; i != end_i; i = (i + 1) % chain->size()) {
        Vertex vv = mesh.from_vertex(chain->at(i));
        Vector3f p = vpositions[vv];
        Vector3f n = vnormals[vv];
        real_t d = std::fmax(std::abs((int_v1 - p).dot(n)),
                             std::abs((int_v2 - p).dot(n)));
        h_dist = std::fmax(h_dist, d);
        small_loop.emplace_back(mesh.edge(chain->at(i)));

        // Check boundary T-junction
        if (is_cusp[vv]) {
          size_t side = determine_turn_side(mesh, camera, chain, i);
          seen_turns[side].emplace_back(vv);
        }
      }

      if (seen_turns[0].empty() && seen_turns[1].size() == 2 &&
          h_dist < loop_epsilon) {
        case2_count++;
        logger().info(
            "tag_simplification_edges_flipped_small_loops: Found a small "
            "loop at v{} - v{}.",
            v1, v2);
        for (auto const &e : small_loop)
          feasibility_collapsing[e] = v.idx();
      } else if (seen_turns[0].empty() && seen_turns[1].size() == 2) {
        logger().info(
            "tag_simplification_edges_flipped_small_loops: A potential small "
            "loop at v{} - v{}. Distance: {}",
            v1, v2, h_dist);
      }
    }
  }

  logger().info("#Case 2: {}", case2_count);

  mesh.get_patch_chains().clear();
  mesh.get_oriented_chains().clear();

  return true;
}
