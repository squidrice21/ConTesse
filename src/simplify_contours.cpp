// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "simplify_contours.h"
#include "common.h"

#include <memory>
#include <spdlog/fmt/ostr.h>
#include <unordered_set>
#include <vector>

#include "tag_cusp_facing.h"

void simplify_contours(Mesh &mesh, Camera const &camera,
                       real_t simplify_epsilon) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_cusp = mesh.get_vertex_property<bool>("v:cusp");
  auto cusp_facing = mesh.get_vertex_property<FacingType>("v:cusp_facing");
  contess_assert_msg(is_cusp && cusp_facing, "Cusp information is missing.");

  auto is_sim_cusp = mesh.vertex_property<bool>("v:sim_cusp");
  if (!is_sim_cusp) {
    is_sim_cusp = mesh.add_vertex_property<bool>("v:sim_cusp");
  }
  is_sim_cusp.vector().assign(is_sim_cusp.vector().size(), false);
  auto sim_cusp_facing = mesh.vertex_property<FacingType>("v:sim_cusp_facing");
  if (!sim_cusp_facing) {
    sim_cusp_facing = mesh.add_vertex_property<FacingType>("v:sim_cusp_facing");
  }
  sim_cusp_facing.vector().assign(sim_cusp_facing.vector().size(),
                                  FacingType::UNDEFINED);

  contess_assert_msg(!mesh.get_oriented_chains().empty(),
                     "Contour chain is missing.");

  auto loop_info = mesh.get_edge_property<Vector4f>("e:loop");
  contess_assert_msg(loop_info, "Contour self-intersection test is missing.");

  auto concave = mesh.get_edge_property<bool>("e:concave");
  contess_assert_msg(concave, "Contour convexity labeling is missing.");

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
      if (dist < simplify_epsilon) {
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

      next_chain_idx[start_i] = end_i;
    }

    // 4. Trace the new chain
    // Find a not deleted starting position whose previous position is also not
    // deleted
    size_t start_i = 0;
    for (; start_i < chain->size() && is_in_removed_range(start_i) &&
           (start_i == 0 || is_in_removed_range(start_i - 1));
         start_i++) {
    }

    // Trace and create new cusps
    auto get_2d_tangent = [&](Vector3f const &v1,
                              Vector3f const &v2) -> Vector2f {
      Vector2f v1_2d =
          project(v1, camera.viewMatrix().matrix(), camera.projectionMatrix(),
                  Vector2i(camera.vpWidth(), camera.vpHeight()))
              .head(2);
      Vector2f v2_2d =
          project(v2, camera.viewMatrix().matrix(), camera.projectionMatrix(),
                  Vector2i(camera.vpWidth(), camera.vpHeight()))
              .head(2);
      return (v2_2d - v1_2d).normalized();
    };
    std::shared_ptr<SimplifiedChain> simp_chain =
        std::make_shared<SimplifiedChain>();
    size_t i = start_i;
    do {
      size_t next_i = (i + 1) % chain->size();
      VirtualHalfedge v_he({Halfedge(-1), Vertex(-1), Vertex(-1)});

      if (next_chain_idx.find(i) == next_chain_idx.end()) {
        v_he.he = chain->at(i);
        simp_chain->emplace_back(v_he);

        // Normal neighborhood, copy the cusp info
        if (next_chain_idx.find(next_i) == next_chain_idx.end()) {
          Vertex to_v = mesh.to_vertex(v_he.he);
          is_sim_cusp[to_v] = is_cusp[to_v];
          sim_cusp_facing[to_v] = cusp_facing[to_v];
        }
      } else {
        // At a simplified connection between vertices
        next_i = next_chain_idx[i];
        v_he.he = Halfedge();
        v_he.from_v = mesh.from_vertex(chain->at(i));
        v_he.to_v = mesh.to_vertex(chain->at(next_i));

        simp_chain->emplace_back(v_he);

        // Determine the cusp
        Vector2f skip_tan =
            get_2d_tangent(vpositions[v_he.from_v], vpositions[v_he.to_v]);
        size_t next_next_i = (next_i + 1) % chain->size();
        size_t prev_i = (i + chain->size() - 1) % chain->size();
        Vector2f prev_tan =
            get_2d_tangent(vpositions[mesh.from_vertex(chain->at(prev_i))],
                           vpositions[v_he.from_v]);
        Vector2f next_tan =
            get_2d_tangent(vpositions[v_he.to_v],
                           vpositions[mesh.to_vertex(chain->at(next_next_i))]);

        // The convexities match, no new cusp
        auto he1 = mesh.find_halfedge(mesh.from_vertex(chain->at(prev_i)),
                                      v_he.from_v);
        auto he2 = mesh.find_halfedge(v_he.to_v,
                                      mesh.to_vertex(chain->at(next_next_i)));
        if (concave[mesh.edge(he1)] != concave[mesh.edge(he2)]) {
          // Convex connected to concave, has a new cusp
          Vertex left, center, right;
          if (skip_tan.dot(prev_tan) < skip_tan.dot(next_tan)) {
            center = v_he.from_v;
            left = mesh.from_vertex(chain->at(prev_i));
            right = v_he.to_v;
          } else {
            center = v_he.to_v;
            right = mesh.to_vertex(chain->at(next_next_i));
            left = v_he.from_v;
          }

          // Determine the facing
          FacingType facing = get_cusp_facing(mesh, camera, center);

          is_sim_cusp[center] = true;
          sim_cusp_facing[center] = facing;
        }

        // Skip one more edge
        if (next_i != start_i)
          next_i = (next_i + 1) % chain->size();
      }

      i = next_i;
    } while (i != start_i);

    mesh.add_simplified_chain(patch_contour.first, simp_chain);
  }
}
