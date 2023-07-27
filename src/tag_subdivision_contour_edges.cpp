// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "tag_subdivision_contour_edges.h"

#include <spdlog/fmt/ostr.h>

void tag_subdivision_contour_edges(Mesh &mesh, Camera const &camera,
                                   real_t loop_distance_threshold,
                                   real_t min_2d_length) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto subdiv_contour = mesh.get_edge_property<bool>("e:subdiv_contour");
  if (!subdiv_contour) {
    subdiv_contour = mesh.add_edge_property<bool>("e:subdiv_contour", false);
  }
  subdiv_contour = mesh.edge_property<bool>("e:subdiv_contour");
  subdiv_contour.vector().assign(subdiv_contour.vector().size(), false);

  contess_assert_msg(!mesh.get_oriented_chains().empty(),
                     "Contour chain is missing.");

  auto loop_info = mesh.get_edge_property<Vector4f>("e:loop");
  contess_assert_msg(loop_info, "Contour self-intersection test is missing.");

  Vector2i viewport = Vector2i(camera.vpWidth(), camera.vpHeight());
  for (auto const &patch_contour : mesh.get_const_patch_chains()) {
    auto chain = patch_contour.second;

    // 1. Index the edges wrt the chain
    std::map<int, int> eid2idx;
    for (auto const &he : *chain) {
      eid2idx[mesh.edge(he).idx()] = eid2idx.size();
    }

    // 2. Check the simplification condition at the starting vertices
    std::vector<std::pair<real_t, int>> loop_distance;
    for (auto const &he : *chain) {
      Edge e = mesh.edge(he);
      bool is_start = (loop_info[e][2] == 1);
      if (!is_start)
        continue;

      // Can trim this loop
      real_t dist = loop_info[e][3];
      if (dist < loop_distance_threshold) {
        // The intersection is indeed in this loop (for the loops containing
        // boundaries)
        Edge e_end(loop_info[e][0]);
        if (!eid2idx.count(e_end.idx()))
          continue;

        loop_distance.emplace_back(-dist, e.idx());
      }
    }

    // 3. Remove from the larger loops to avoid duplicate removal of nested
    // loops
    std::sort(loop_distance.begin(), loop_distance.end());
    std::map<size_t, size_t> next_chain_idx;

    auto is_in_removed_range = [&next_chain_idx](size_t idx) -> bool {
      bool is_in = false;
      for (auto const &removed_range : next_chain_idx) {
        if (removed_range.first < removed_range.second)
          is_in = (idx > removed_range.first && idx < removed_range.second);
        else
          is_in = (idx > removed_range.first || idx < removed_range.second);
        if (is_in)
          break;
      }

      return is_in;
    };

    for (auto const &loop : loop_distance) {
      Edge e_start(loop.second);
      size_t start_i = eid2idx[e_start.idx()];

      if (is_in_removed_range(start_i))
        continue;

      Edge e_end(loop_info[e_start][0]);
      size_t end_i = eid2idx[e_end.idx()];

      next_chain_idx[start_i] = (end_i + 1) % chain->size();
    }

    // 4. Label long edges within this range
    for (auto const &loop_range : next_chain_idx) {
      for (size_t i = loop_range.first; i != loop_range.second;
           i = (i + 1) % chain->size()) {
        auto he = chain->at(i);
        Vector2f v1_2d = project(vpositions[mesh.from_vertex(he)],
                                 camera.viewMatrix().matrix(),
                                 camera.projectionMatrix(), viewport)
                             .head<2>();
        Vector2f v2_2d = project(vpositions[mesh.to_vertex(he)],
                                 camera.viewMatrix().matrix(),
                                 camera.projectionMatrix(), viewport)
                             .head<2>();
        real_t len_2d = (v1_2d - v2_2d).norm();
        if (len_2d > min_2d_length) {
          subdiv_contour[mesh.edge(he)] = true;
        }
      }
    }
  }
}
