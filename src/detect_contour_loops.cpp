// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "detect_contour_loops.h"

#include "insert_planar_map_intersections.h"
#include <spdlog/fmt/ostr.h>
#include <string>
#include <unordered_set>
#include <vector>

void detect_contour_loops(Mesh &mesh, Camera const &camera) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto vnormals = mesh.get_vertex_property<Vector3f>("v:normal");

  // Vector4f: 0: Intersected edge; 1: Intersected time; 2: Is starting point;
  // 3: Hausdorff distance.
  auto loop_info = mesh.edge_property<Vector4f>("e:loop");
  if (!loop_info) {
    loop_info = mesh.add_edge_property<Vector4f>("e:loop");
  }
  loop_info.vector().assign(loop_info.vector().size(), -1 * Vector4f::Ones());

  contess_assert_msg(!mesh.get_oriented_chains().empty(),
                     "Contour chain is missing.");

  for (auto const &chain : mesh.get_oriented_chains()) {
    // 1. Get all self intersections of the loop
    std::vector<Edge> chain_edges;
    chain_edges.reserve(chain->size());
    std::map<int, int> eid2idx;

    // 2. Precompute the loop length
    std::vector<real_t> chain_length;
    chain_length.reserve(chain->size() + 1);
    chain_length.emplace_back(0);
    for (auto const &he : *chain) {
      eid2idx[mesh.edge(he).idx()] = chain_edges.size();
      chain_edges.emplace_back(mesh.edge(he));

      chain_length.emplace_back(
          chain_length.back() +
          (vpositions[mesh.from_vertex(he)] - vpositions[mesh.to_vertex(he)])
              .norm());
    }
    std::map<int, std::vector<std::pair<real_t, int>>> intersections;
    edges_intersection(mesh, camera, chain_edges, intersections);

    // 3. Compute the Hausdorff distance along the surface normal if the smaller
    // loop is deleted
    std::unordered_set<int> visited_edges;
    for (auto const &intersection : intersections) {
      // For now drop the intersections if there are multiple
      if (intersection.second.size() > 1) {
        logger().warn("Dropping intersections with edge " +
                      std::to_string(intersection.first) + ".");
      }

      if (visited_edges.find(intersection.first) != visited_edges.end() ||
          visited_edges.find(intersection.second.front().second) !=
              visited_edges.end())
        continue;

      int e1 = intersection.first;
      int e2 = intersection.second.front().second;
      visited_edges.insert(e1);
      visited_edges.insert(e2);

      // Book keeping on edges
      loop_info[Edge(e1)][0] = e2;
      loop_info[Edge(e2)][0] = e1;
      loop_info[Edge(e1)][1] = intersection.second.front().first;
      loop_info[Edge(e2)][1] =
          std::find_if(
              intersections[e2].begin(), intersections[e2].end(),
              [&e1](std::pair<real_t, int> const &x) { return x.second == e1; })
              ->first;

      // Determine the smaller loop
      // From first edge to the second
      real_t e1_len = chain_length[eid2idx[e1]];
      real_t e2_len = chain_length[eid2idx[e2]];

      real_t e1_to_e2_len =
          ((e2_len < e1_len) ? e2_len + chain_length.back() : e2_len) - e1_len;
      real_t e2_to_e1_len =
          ((e1_len < e2_len) ? e1_len + chain_length.back() : e1_len) - e2_len;

      size_t start_i = eid2idx[e2];
      size_t end_i = eid2idx[e1];

      loop_info[Edge(e1)][2] = 0;
      loop_info[Edge(e2)][2] = 1;
      if (e1_to_e2_len < e2_to_e1_len) {
        loop_info[Edge(e1)][2] = 1;
        loop_info[Edge(e2)][2] = 0;
        start_i = eid2idx[e1];
        end_i = eid2idx[e2];
      }

      // Walk on halfedges
      Vector3f int_v1 =
          (1 - loop_info[Edge(e1)][1]) * vpositions[mesh.vertex(Edge(e1), 0)] +
          loop_info[Edge(e1)][1] * vpositions[mesh.vertex(Edge(e1), 1)];
      Vector3f int_v2 =
          (1 - loop_info[Edge(e2)][1]) * vpositions[mesh.vertex(Edge(e2), 0)] +
          loop_info[Edge(e2)][1] * vpositions[mesh.vertex(Edge(e2), 1)];
      real_t h_dist = -1;
      for (size_t i = start_i; i != end_i; i = (i + 1) % chain->size()) {
        Vector3f v = vpositions[mesh.to_vertex(chain->at(i))];
        Vector3f n = vnormals[mesh.to_vertex(chain->at(i))];
        real_t d = std::fmax(std::abs((int_v1 - v).dot(n)),
                             std::abs((int_v2 - v).dot(n)));
        h_dist = std::fmax(h_dist, d);
      }
      loop_info[Edge(e1)][3] = h_dist;
      loop_info[Edge(e2)][3] = h_dist;
    }
  }
}
