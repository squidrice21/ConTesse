// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "shrink_cut.h"

#include "common.h"
#include "cut_patch_to_disk.h"
#include <deque>
#include <spdlog/fmt/ostr.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

static real_t PATH_BIG_NUMBER = 1e4;

void trace_cut(Mesh const &mesh, Vertex const &v, int patch_id,
               std::vector<Vertex> &reached_vv) {
  auto trace_cut_record = mesh.get_edge_property<bool>("e:trace_cut");
  contess_assert(trace_cut_record);
  auto cut = mesh.get_edge_property<int>("e:disk_cut");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");

  std::queue<Vertex> vertex_queue;
  std::unordered_set<int> seen_vertices;
  std::unordered_map<int, Vertex> prev_v;
  vertex_queue.push(v);
  prev_v[v.idx()] = Vertex();

  while (!vertex_queue.empty()) {
    Vertex current_v = vertex_queue.front();
    vertex_queue.pop();

    if (seen_vertices.count(current_v.idx()))
      continue;

    seen_vertices.emplace(current_v.idx());

    // Cut branching vertex
    int degree = count_cut_degree(mesh, current_v, patch_id);
    if (degree > 2 && current_v != v) {
      reached_vv.emplace_back(current_v);

      // Record the path
      Vertex vv = current_v;
      while (vv != v) {
        contess_assert(prev_v.count(vv.idx()));
        Vertex p_v = prev_v[vv.idx()];
        Edge e = mesh.find_edge(vv, p_v);
        contess_assert(e.is_valid());

        trace_cut_record[e] = true;
        vv = p_v;
      }
      return;
    }

    // Reached contour/manifold boundary
    auto hit = mesh.halfedges(current_v), hend = hit;
    do {
      Edge hit_e = mesh.edge(*hit);

      if (trace_cut_record[hit_e])
        continue;

      if (seen_vertices.count(mesh.to_vertex(*hit).idx()) ||
          is_contour[mesh.edge(*hit)] >= 0 || mesh.is_boundary(mesh.edge(*hit)))
        continue;

      // Finds a cut edge in the same patch
      if (cut[hit_e] == patch_id) {
        prev_v[mesh.to_vertex(*hit).idx()] = current_v;
        vertex_queue.push(mesh.to_vertex(*hit));
      }
    } while (++hit != hend);
  }
}

// A*
bool shortest_path_image_plane(Mesh &mesh, Camera const &camera,
                               Vertex const &src, Vertex const &dst,
                               int patch_id) {
  auto patchID = mesh.get_face_property<int>("f:patchID");
  auto cut_candidate = mesh.get_edge_property<bool>("e:cut_candidate");
  auto patch_removed = mesh.get_edge_property<bool>("e:patch_removed");
  contess_assert(cut_candidate && mesh.get_edge_property<int>("e:disk_cut2"));
  auto cut2 = mesh.edge_property<int>("e:disk_cut2");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");

  std::unordered_map<int, real_t> vertex_distance;
  std::unordered_set<int> seen_vertices;
  std::unordered_map<int, Vertex> prev_v;

  auto v_cmp = [](std::pair<Vertex, real_t> const &left,
                  std::pair<Vertex, real_t> const &right) -> bool {
    return left.second > right.second;
  };
  std::priority_queue<std::pair<Vertex, real_t>,
                      std::vector<std::pair<Vertex, real_t>>, decltype(v_cmp)>
      vertex_queue(v_cmp);

  auto proj_distance = [&](Vertex v1, Vertex v2) -> real_t {
    auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
    Vector2f vv1 = project(vpositions[v1], camera.viewMatrix().matrix(),
                           camera.projectionMatrix(),
                           Vector2i(camera.vpWidth(), camera.vpHeight()))
                       .head(2);
    Vector2f vv2 = project(vpositions[v2], camera.viewMatrix().matrix(),
                           camera.projectionMatrix(),
                           Vector2i(camera.vpWidth(), camera.vpHeight()))
                       .head(2);
    return (vv2 - vv1).norm();
  };

  // We don't want to pass or double trace an existing cut
  auto is_vertex_valid = [&](Vertex v, bool &cut_valid, bool &boundary_valid) {
    cut_valid = boundary_valid = true;
    if (v == dst) {
      return;
    }
    auto hit = mesh.halfedges(v), hend = hit;
    do {
      Edge hit_e = mesh.edge(*hit);
      if (cut2[hit_e] >= 0) {
        cut_valid = false;
      }
      if (is_contour[hit_e] >= 0 || mesh.is_boundary(hit_e)) {
        boundary_valid = false;
        return;
      }
    } while (++hit != hend);
  };

  prev_v[src.idx()] = Vertex();
  vertex_queue.push(std::make_pair(src, proj_distance(src, dst)));
  vertex_distance[src.idx()] = 0;

  while (!vertex_queue.empty()) {
    Vertex current_v = vertex_queue.top().first;
    vertex_queue.pop();

    if (seen_vertices.count(current_v.idx()))
      continue;

    seen_vertices.emplace(current_v.idx());

    // Reached destination
    if (current_v == dst) {
      break;
    }

    auto hit = mesh.halfedges(current_v), hend = hit;
    do {
      Edge hit_e = mesh.edge(*hit);
      Halfedge oppo = mesh.opposite_halfedge(*hit);

      if (seen_vertices.count(mesh.to_vertex(*hit).idx()) ||
          is_contour[mesh.edge(*hit)] >= 0 || mesh.is_boundary(mesh.edge(*hit)))
        continue;

      bool is_candidate =
          cut_candidate[hit_e] && (!patch_removed || !patch_removed[hit_e]);

      // Grow within the same patch
      if ((mesh.face(*hit).is_valid() &&
           patchID[mesh.face(*hit)] == patch_id) ||
          (oppo.is_valid() && mesh.face(oppo).is_valid() &&
           patchID[mesh.face(oppo)] == patch_id)) {
        Vertex next_v = mesh.to_vertex(*hit);

        // Use a big number to avoid passing this edge
        // (Only pass it if there's no valid solution)
        bool cut_valid, boundary_valid;
        is_vertex_valid(next_v, cut_valid, boundary_valid);

        // We can't pass vertices on boundary
        if (!boundary_valid)
          continue;

        real_t step_length = (is_candidate && cut_valid)
                                 ? proj_distance(next_v, current_v)
                                 : PATH_BIG_NUMBER;

        if (!vertex_distance.count(next_v.idx()) ||
            vertex_distance[current_v.idx()] + step_length <
                vertex_distance[next_v.idx()]) {
          vertex_distance[next_v.idx()] =
              vertex_distance[current_v.idx()] + step_length;
          prev_v[next_v.idx()] = current_v;
        }

        vertex_queue.push(
            std::make_pair(next_v, vertex_distance[next_v.idx()] +
                                       proj_distance(next_v, dst)));
      }
    } while (++hit != hend);
  }

  contess_assert_msg(prev_v.count(dst.idx()),
                     "shrink_cut: Unreachable v" + std::to_string(src.idx()) +
                         " -> v" + std::to_string(dst.idx()));

  // Record the path
  Vertex v = dst;
  logger().info("A* p{}, {} -> {}: {} - {}", patch_id, src, dst,
                vertex_distance[dst.idx()],
                vertex_distance[dst.idx()] < PATH_BIG_NUMBER);
  if (vertex_distance[dst.idx()] >= PATH_BIG_NUMBER)
    return false;
  while (v != src) {
    contess_assert(prev_v.count(v.idx()));
    Vertex p_v = prev_v[v.idx()];
    Edge e = mesh.find_edge(v, p_v);
    contess_assert(e.is_valid());

    cut2[e] = patch_id;
    v = p_v;
  }
  return true;
}

void shrink_cut(Mesh &mesh, Camera const &camera) {
  auto patchID = mesh.get_face_property<int>("f:patchID");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  if (!mesh.get_edge_property<int>("e:disk_cut")) {
    logger().warn("shrink_cut: Must call cut_patch_to_disk first.");
    return;
  }
  auto cut = mesh.edge_property<int>("e:disk_cut");
  auto cut2 = mesh.get_edge_property<int>("e:disk_cut2");
  if (!cut2) {
    cut2 = mesh.add_edge_property<int>("e:disk_cut2", -1);
  } else {
    cut2 = mesh.edge_property<int>("e:disk_cut2");
    cut2.vector().assign(cut2.vector().size(), -1);
  }

  // Debug info
  auto is_cut_branch = mesh.get_vertex_property<bool>("v:cut_branch");
  if (!is_cut_branch) {
    is_cut_branch = mesh.add_vertex_property<bool>("v:cut_branch", false);
  }
  is_cut_branch = mesh.vertex_property<bool>("v:cut_branch", false);
  is_cut_branch.vector().assign(is_cut_branch.vector().size(), false);

  // 1. Find all cut vertices with degree > 2
  std::unordered_set<int> seen_vertices;
  std::unordered_map<int, std::vector<Vertex>> cut_vertices;
  for (size_t i = 0; i < mesh.n_edges(); i++) {
    Edge e(i);
    if (is_contour[e] < 0 && !mesh.is_boundary(e) && cut[e] < 0)
      continue;

    std::unordered_set<int> adj_patch;
    for (size_t j = 0; j < 2; j++) {
      Halfedge he = mesh.halfedge(e, j);
      if (mesh.face(he).is_valid())
        adj_patch.emplace(patchID[mesh.face(he)]);
    }

    std::vector<Vertex> vv({mesh.vertex(e, 0), mesh.vertex(e, 1)});
    for (auto const &v : vv) {
      if (seen_vertices.count(v.idx()))
        continue;
      for (auto patch_id : adj_patch) {
        int degree = count_cut_degree(mesh, v, patch_id);
        if (degree > 2)
          cut_vertices[patch_id].emplace_back(v);
        seen_vertices.emplace(v.idx());
      }
    }
  }

  // 2. Pair vertices by tracing from contour
  // std::unordered_set<std::pair<int, int>> seen_pairs;
  auto trace_cut_record = mesh.add_edge_property<bool>("e:trace_cut", false);
  std::unordered_map<int, std::vector<std::pair<Vertex, Vertex>>> vertex_pairs;
  for (auto patch_vv : cut_vertices) {
    int patch_id = patch_vv.first;
    for (auto v : patch_vv.second) {
      std::vector<Vertex> reached_vv;
      // Trace all the cut paths from the vertex
      do {
        reached_vv.clear();
        trace_cut(mesh, v, patch_id, reached_vv);
        for (auto r_vv : reached_vv) {
          if (!vertex_pairs.count(patch_id))
            vertex_pairs[patch_id] = std::vector<std::pair<Vertex, Vertex>>();
          vertex_pairs[patch_id].emplace_back(std::make_pair(v, r_vv));
        }
      } while (!reached_vv.empty());
    }
  }
  mesh.remove_edge_property(trace_cut_record);

  // 3. Find shortest valid path between pair (image plane Euclidean)
  std::unordered_set<int> successful_patches;
  for (auto const &pair_all : vertex_pairs) {
    int patch_id = pair_all.first;

    // Order the pairs to trace
    std::deque<std::pair<Vertex, Vertex>> vertex_queue;
    for (auto pair : pair_all.second) {
      vertex_queue.push_front(pair);
    }

    bool patch_successful = true;
    size_t itr = 0;
    do {
      logger().info("Shrink p{} round {}", patch_id, itr);
      patch_successful = true;
      std::deque<std::pair<Vertex, Vertex>> vertex_queue2;
      while (!vertex_queue.empty()) {
        auto pair = vertex_queue.front();
        vertex_queue.pop_front();
        Vertex src = pair.first;
        Vertex dst = pair.second;

        if (!patch_successful) {
          vertex_queue2.push_back(pair);
          continue;
        }

        bool path_found =
            shortest_path_image_plane(mesh, camera, src, dst, patch_id);
        if (!path_found) {
          patch_successful = false;
          vertex_queue2.push_front(pair);
        } else {
          vertex_queue2.push_back(pair);
        }
      }

      if (!patch_successful) {
        for (size_t i = 0; i < mesh.n_edges(); i++) {
          Edge e(i);
          if (cut2[e] == patch_id)
            cut2[e] = -1;
        }
        vertex_queue = vertex_queue2;
      }
      itr++;
    } while (!patch_successful && itr < 3);

    if (patch_successful)
      successful_patches.emplace(patch_id);
    else
      logger().warn("shrink_cut: Cannot shrink cuts in patch {}.", patch_id);
  }

  // 4. Swap the stored cut
  for (size_t i = 0; i < mesh.n_edges(); i++) {
    Edge e(i);

    if (cut2[e] == cut[e])
      continue;

    if (is_contour[e] >= 0 || mesh.is_boundary(e) ||
        (cut[e] >= 0 && !successful_patches.count(cut[e]))) {
      cut2[e] = cut[e];
      continue;
    }

    int temp_cut = cut[e];
    cut[e] = cut2[e];
    cut2[e] = temp_cut;
  }
}
