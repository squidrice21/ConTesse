// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "collapse_flipped_edges.h"
#include <deque>
#include <limits>
#include <spdlog/fmt/ostr.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "common.h"
#include "fix_flipped_faces.h"
#include "subdiv_common.h"
#include "subdivide_contour_edges.h"

bool collapse_edges(Mesh &mesh, Camera const &camera, Vertex const &v,
                    Halfedge const &h1, Halfedge const &h2) {
  auto param_loc = mesh.halfedge_property<Param_loc>("h:param_loc");
  auto is_contour = mesh.edge_property<real_t>("e:contour");
  auto disk_cut = mesh.edge_property<int>("e:disk_cut");
  contess_assert_msg(disk_cut,
                     "collapse_fipped_edges: Patch needs to be cut to disk.");

  bool e_contour = is_contour[mesh.edge(h1)];
  bool e_cut = disk_cut[mesh.edge(h1)];

  // Collect 1-ring
  Vertex v1 = mesh.to_vertex(h1), v2 = mesh.to_vertex(h2);

  auto hit = mesh.halfedges(v), hit_end = hit;
  std::vector<Vertex> hole_vertices;
  std::vector<Vertex> hole_vertices_v2;
  do {
    if (mesh.to_vertex(*hit) == v1)
      break;
  } while (++hit != hit_end);
  hit_end = hit;
  do {
    hole_vertices.emplace_back(mesh.to_vertex(*hit));
  } while (++hit != hit_end);

  // Create the backup plan
  hit = mesh.halfedges(v), hit_end = hit;
  do {
    if (mesh.to_vertex(*hit) == v2)
      break;
  } while (++hit != hit_end);
  hit_end = hit;
  do {
    hole_vertices_v2.emplace_back(mesh.to_vertex(*hit));
  } while (++hit != hit_end);
  // collapse_fipped_edges: Cannot order vertices in the 1-ring.
  if (!(hole_vertices.front() == v1 && hole_vertices_v2.front() == v2))
    return false;

  // surface_mesh may encontour problem when adding face connecting a isolated
  // vertex, which would happen when we have boundary edges in the one-ring
  // Hack to solve this is to temporarily add faces to make these boundary edges
  // interior edges, and delete them latter.
  auto is_boundary_edge = [&](int i1, int i2) -> bool {
    Edge e = mesh.find_edge(Vertex(i1), Vertex(i2));
    if (!e.is_valid())
      e = mesh.find_edge(Vertex(i2), Vertex(i1));
    if (e.is_valid() && mesh.is_boundary(e)) {
      return true;
    }
    return false;
  };
  bool seen_boundary = false;
  for (size_t i = 0; i < hole_vertices.size(); i++) {
    if (is_boundary_edge(hole_vertices[i].idx(),
                         hole_vertices[(i + 1) % hole_vertices.size()].idx())) {
      seen_boundary = true;
      break;
    }
  }

  std::unordered_set<int> skip_vertices;
  Vertex orig_v;
  std::vector<Face> pad_faces;
  if (seen_boundary) {
    orig_v = mesh.add_vertex(Vector3f::Ones() * 1e3);
    for (size_t i = 0; i < hole_vertices.size(); i++) {
      if (is_boundary_edge(
              hole_vertices[i].idx(),
              hole_vertices[(i + 1) % hole_vertices.size()].idx())) {
        Face f = mesh.add_face(std::vector<Vertex>(
            {orig_v, hole_vertices[(i + 1) % hole_vertices.size()],
             hole_vertices[i]}));
        pad_faces.emplace_back(f);
      }
    }
  }

  // Parameterization book keeping
  std::unordered_map<int, std::unordered_map<int, Param_loc>> end_param;
  std::unordered_map<int, std::unordered_map<int, Param_loc>> he_param;
  auto add_end_param = [&](Vertex const &vv) {
    auto hit = mesh.halfedges(vv), hit_end = hit;
    do {
      if (!end_param.count(vv.idx())) {
        end_param[vv.idx()] = std::unordered_map<int, Param_loc>();
      }
      Param_loc param = param_loc[*hit];
      end_param[vv.idx()][param.ptexIndex] = param;
    } while (++hit != hit_end);
  };
  for (size_t i = 0; i < hole_vertices.size(); i++) {
    add_end_param(hole_vertices[i]);
  }
  hit = mesh.halfedges(v), hit_end = hit;
  do {
    Face f = mesh.face(*hit);
    if (!f.is_valid())
      continue;
    auto f_hit = mesh.halfedges(f), f_hit_end = f_hit;
    do {
      int from_idx = mesh.from_vertex(*f_hit).idx();
      int to_idx = mesh.to_vertex(*f_hit).idx();
      if (!he_param.count(from_idx)) {
        he_param[from_idx] = std::unordered_map<int, Param_loc>();
      }
      Param_loc param = param_loc[*f_hit];
      he_param[from_idx][to_idx] = param;
    } while (++f_hit != f_hit_end);
  } while (++hit != hit_end);

  // Delete all adjacent faces
  mesh.delete_vertex(v);

  // Fill faces by connecting edges to v1
  for (size_t i = 1; i + 1 < hole_vertices.size(); i++) {
    if (skip_vertices.count(hole_vertices[i].idx()))
      continue;
    Face f = mesh.add_face(
        std::vector<Vertex>({v1, hole_vertices[i], hole_vertices[i + 1]}));
    if (!f.is_valid()) {
      logger().warn("Warning: Cannot add face {}, {}, {}.", v1,
                    hole_vertices[i], hole_vertices[i + 1]);
      continue;
    }

    // Fill the parameterization along the 1-ring boundary
    int ptex_idx = -1;
    auto f_hit = mesh.halfedges(f), f_hit_end = f_hit;
    do {
      int from_idx = mesh.from_vertex(*f_hit).idx();
      int to_idx = mesh.to_vertex(*f_hit).idx();
      if (he_param.count(from_idx) && he_param[from_idx].count(to_idx)) {
        auto param = he_param[from_idx][to_idx];
        param_loc[*f_hit] = param;
        ptex_idx = param.ptexIndex;
      }
    } while (++f_hit != f_hit_end);

    // Note it's fine now since we don't call evaluation after collapsing for
    // now.
    if (ptex_idx < 0)
      continue;

    contess_assert(ptex_idx >= 0);

    // Fill the ones not along the 1-ring boundary
    f_hit = mesh.halfedges(f), f_hit_end = f_hit;
    do {
      int from_idx = mesh.from_vertex(*f_hit).idx();
      if (param_loc[*f_hit].is_valid())
        continue;
      if (end_param.count(from_idx) && end_param[from_idx].count(ptex_idx)) {
        auto param = end_param[from_idx][ptex_idx];
        param_loc[*f_hit] = param;
      }
    } while (++f_hit != f_hit_end);
  }

  // Refill the edge properties
  Edge e12 = mesh.find_edge(v1, v2);
  is_contour[e12] = e_contour;
  disk_cut[e12] = e_cut;

  // Remove padding
  if (orig_v.is_valid())
    mesh.delete_vertex(orig_v);

  return true;
}

void resolve_overlapping_cells(Mesh &mesh, Camera const &camera,
                               real_t parallel_offset) {
  auto is_cusp = mesh.get_vertex_property<bool>("v:cusp");
  auto stationary = mesh.get_vertex_property<bool>("v:stationary");
  auto vpositions = mesh.vertex_property<Vector3f>("v:point");
  auto is_contour = mesh.edge_property<real_t>("e:contour");
  auto disk_cut = mesh.edge_property<int>("e:disk_cut");
  contess_assert_msg(disk_cut,
                     "collapse_fipped_edges: Patch needs to be cut to disk.");
  auto cusp_facing = mesh.get_vertex_property<FacingType>("v:cusp_facing");
  auto intersection_2d = mesh.get_vertex_property<int>("v:intersection_2d");
  contess_assert_msg(
      intersection_2d && cusp_facing && is_cusp,
      "collapse_fipped_edges: 2D intersections are not created.");

  auto feasibility_collapsing =
      mesh.edge_property<int>("e:feasibility_collapsing");
  contess_assert_msg(
      feasibility_collapsing,
      "collapse_fipped_edges: Need to label edges for collapsing first.");

  auto concave = mesh.get_edge_property<bool>("e:concave");
  contess_assert_msg(
      concave,
      "collapse_fipped_edges: Need to label convex/concave edges first.");
  concave = mesh.edge_property<bool>("e:concave");

  std::deque<Edge> to_collapse_edges;
  for (auto eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit) {
    Edge e = *eit;
    Vertex v1 = mesh.vertex(e, 0);
    Vertex v2 = mesh.vertex(e, 1);

    if (!mesh.is_deleted(v1) && !mesh.is_deleted(v2) &&
        feasibility_collapsing[e] >= 0)
      to_collapse_edges.push_back(e);
  }

  std::unordered_set<int> collapsed_edges;
  while (!to_collapse_edges.empty()) {
    Edge e = to_collapse_edges.front();
    to_collapse_edges.pop_front();
    if (collapsed_edges.count(e.idx()))
      continue;

    Vertex v1 = mesh.vertex(e, 0);
    Vertex v2 = mesh.vertex(e, 1);

    // 1. Collapsed small loop
    if (intersection_2d[v1] == v2.idx()) {
      // This case would be reached by loop removal
      feasibility_collapsing[e] = -1;
    }
    // 2. Parallel lines
    else {
      // Find the paired parallel line
      Vertex v3(intersection_2d[v1]);
      Vertex v4(intersection_2d[v2]);
      bool v3_is_valid = v3.is_valid();

      logger().info("Flipping case 2 {}, {}, {}, {}.", v1, v2, v3, v4);

      if (!v3.is_valid() && !v4.is_valid()) {
        logger().warn("False collapsing: {} - {}", v1, v2);
        continue;
      }

      // Collapsed loop with a unmoved vertex
      if (!v3.is_valid()) {
        contess_assert_msg(
            stationary[v1],
            "collapse_fipped_edges: Stationary vertex not labeled: v" +
                std::to_string(v1.idx()));
        v3 = v1;
      }
      if (!v4.is_valid()) {
        contess_assert_msg(
            stationary[v2],
            "collapse_fipped_edges: Stationary vertex not labeled: v" +
                std::to_string(v2.idx()));
        v4 = v2;
      }

      contess_assert_msg(
          v3.is_valid() && v4.is_valid(),
          "collapse_fipped_edges: Intersection data structure is broken.");

      Edge e_parallel = mesh.find_edge(v3, v4);
      if (!e_parallel.is_valid())
        e_parallel = mesh.find_edge(v4, v3);

      if (!e_parallel.is_valid() || feasibility_collapsing[e_parallel] < 0) {
        logger().error("FAILED: collapse_flipped_edges cannot handle "
                       "complicated entangling: {}, {}; {}, {}.",
                       v1, v2, v3, v4);
        continue;
      }

      collapsed_edges.insert(e_parallel.idx());

      // Offset the vertices to resolve the overlapping
      auto move_parallel_vertex = [&](Vertex const &vv, Vertex const &vv_paired,
                                      bool &is_concave) -> Vertex {
        Halfedge he_parallel = mesh.find_halfedge(vv, vv_paired);
        if (!he_parallel.is_valid())
          he_parallel = mesh.find_halfedge(vv_paired, vv);

        Halfedge he_dir;
        Vertex start = vv;
        Vertex prev = vv_paired;
        // Find the direction to move
        // Use adjacent faces to find adjacent vertices
        // (since boundary vertex may have problem finding all adjacent
        // vertices)
        std::unordered_set<int> adj_vv;
        auto hit = mesh.halfedges(start), hit_end = hit;
        do {
          Face f = mesh.face(*hit);
          if (f.is_valid()) {
            Vertex vertices[3];
            mesh.verticesOfFace(f, vertices);
            for (auto const &w : vertices)
              if (w != start && w != prev)
                adj_vv.emplace(w.idx());
          }
        } while (++hit != hit_end);

        for (auto const v_idx : adj_vv) {
          Vertex adj_v(v_idx);
          Edge eit = mesh.find_edge(start, adj_v);
          if (!eit.is_valid())
            eit = mesh.find_edge(adj_v, start);
          if (mesh.is_boundary(eit) || is_contour[eit] >= 0) {
            he_dir = mesh.find_halfedge(start, adj_v);
            if (!he_dir.is_valid())
              he_dir = mesh.find_halfedge(adj_v, start);
            is_concave = concave[eit];
            break;
          }
        }

        if (!(he_dir.is_valid() && he_parallel.is_valid())) {
          logger().error("collapse_fipped_edges: Cannot find the edge to "
                         "offset. Contour is not loopy at {} - {}",
                         vv, vv_paired);
          return Vertex();
        }

        contess_assert_msg(he_dir.is_valid() && he_parallel.is_valid(),
                           "collapse_fipped_edges: Cannot find the edge to "
                           "offset. Contour is not loopy at v" +
                               std::to_string(vv.idx()) + " - v" +
                               std::to_string(vv_paired.idx()));

        return (mesh.to_vertex(he_dir) == start) ? mesh.from_vertex(he_dir)
                                                 : mesh.to_vertex(he_dir);
      };

      bool is_concave_12, is_concave_34;
      Vertex v5 = (stationary[v1])
                      ? move_parallel_vertex(v2, v1, is_concave_12)
                      : move_parallel_vertex(v1, v2, is_concave_12);
      Vertex v6 = (stationary[v3])
                      ? move_parallel_vertex(v4, v3, is_concave_34)
                      : move_parallel_vertex(v3, v4, is_concave_34);

      if (!v5.is_valid() || !v6.is_valid()) {
        logger().error("collapse_fipped_edges: Cannot find the edge to "
                       "offset. Contour is not loopy at {} - {}, {} - {}",
                       v1, v2, v3, v4);
        continue;
      }

      // Point to right of he: v1 -> v2
      Vector3f move_dir = (vpositions[v2] - vpositions[v1])
                              .cross((camera.position() - vpositions[v1]))
                              .normalized();

      std::vector<Vertex> vertices({v1, v5, v6});
      std::vector<Vertex> vertices_3d({v1, v2, v3});
      if (stationary[v1]) {
        vertices = std::vector<Vertex>({v2, v6, v5});
        vertices_3d = std::vector<Vertex>({v2, v3, v1});
      }

      // Determine the orientation without projection to 2D
      auto f_ori = igl::predicates::orient3d(
          vpositions[vertices[0]], vpositions[vertices[1]],
          vpositions[vertices[2]], camera.position());
      contess_assert_msg(f_ori != igl::predicates::Orientation::COPLANAR,
                         "Generic camera assumption violated.");
      if (f_ori == igl::predicates::Orientation::NEGATIVE)
        move_dir *= -1;

      size_t itr = 1;
      bool same_ori = true;
      do {
        Vector3f offset_vec = itr * parallel_offset * move_dir;

        // Not moving cusp to avoid introducing new 2D intersections
        // Unless it's the loop case (which is identified by the stationary)
        if ((!is_cusp[v1] && !is_cusp[v2]) || stationary[v1] ||
            stationary[v2]) {
          vpositions[v1] += offset_vec;
          vpositions[v2] += offset_vec;
        }

        concave[mesh.find_edge(v1, v2)] = is_concave_12;

        offset_vec *= -1;
        if ((!is_cusp[v3] && !is_cusp[v4]) || stationary[v1] ||
            stationary[v2]) {
          vpositions[v3] += offset_vec;
          vpositions[v4] += offset_vec;
        }

        concave[e_parallel] = is_concave_34;

        // Check if they are moved enough
        Vector2f v1_2d = project(vpositions[v1], camera.viewMatrix().matrix(),
                                 camera.projectionMatrix(),
                                 Vector2i(camera.vpWidth(), camera.vpHeight()))
                             .head(2);
        Vector2f v2_2d = project(vpositions[v2], camera.viewMatrix().matrix(),
                                 camera.projectionMatrix(),
                                 Vector2i(camera.vpWidth(), camera.vpHeight()))
                             .head(2);
        Vertex oppo_v = (v3_is_valid) ? v3 : v4;
        Vector2f oppo_v_2d =
            project(vpositions[oppo_v], camera.viewMatrix().matrix(),
                    camera.projectionMatrix(),
                    Vector2i(camera.vpWidth(), camera.vpHeight()))
                .head(2);
        Vector2f norm12 = (v2_2d - v1_2d).normalized();
        norm12 = Vector2f(-norm12.y(), norm12.x());
        real_t moved_dist = norm12.dot(oppo_v_2d - v1_2d);
        auto after_ori = igl::predicates::orient3d(
            vpositions[vertices_3d[0]], vpositions[vertices_3d[1]],
            vpositions[vertices_3d[2]], camera.position());
        same_ori =
            (f_ori == after_ori) ||
            (std::abs(moved_dist) < std::numeric_limits<real_t>::epsilon());
        itr++;
      } while (same_ori && itr < 10);
    }
  }
}

bool collapse_flipped_edges(Mesh &mesh, Camera const &camera,
                            real_t parallel_offset) {
  auto stationary = mesh.get_vertex_property<bool>("v:stationary");
  auto feasibility_collapsing =
      mesh.edge_property<int>("e:feasibility_collapsing");
  contess_assert_msg(
      feasibility_collapsing,
      "collapse_fipped_edges: Need to label edges for collapsing first.");

  struct DeleteInfo {
    Vertex v;
    Halfedge h1, h2;
  };
  auto is_collapsible = [&](Vertex const &v, DeleteInfo &info) -> bool {
    // Cusp or boundary-contour joint
    if (stationary[v])
      return false;
    auto hit = mesh.halfedges(v), hit_end = hit;
    std::vector<Halfedge> to_delete_he;
    do {
      if (feasibility_collapsing[mesh.edge(*hit)] >= 0) {
        to_delete_he.emplace_back(*hit);
      }
    } while (++hit != hit_end);

    if (to_delete_he.size() == 2) {
      info.v = v;
      info.h1 = to_delete_he[0];
      info.h2 = to_delete_he[1];
      return true;
    }
    return false;
  };

  // Collapse adjacent flipped edfges
  std::deque<DeleteInfo> delete_vertices;
  for (size_t i = 0; i < mesh.n_vertices(); i++) {
    Vertex v(i);
    DeleteInfo info;

    if (is_collapsible(v, info)) {
      delete_vertices.push_back(info);
    }
  }

  while (!delete_vertices.empty()) {
    DeleteInfo d_info = delete_vertices.front();
    delete_vertices.pop_front();

    Vertex v = d_info.v;
    Halfedge h1 = d_info.h1, h2 = d_info.h2;
    int feasibility_v = feasibility_collapsing[mesh.edge(h1)];
    Vertex v1 = mesh.to_vertex(h1), v2 = mesh.to_vertex(h2);

    bool successful = collapse_edges(mesh, camera, v, h1, h2);
    if (!successful)
      return false;

    Edge e12 = mesh.find_edge(v1, v2);
    feasibility_collapsing[e12] = feasibility_v;

    // No new vertex is created, we only need to update existing record
    DeleteInfo info;
    if (is_collapsible(v1, info)) {
      for (auto itr = delete_vertices.begin(); itr != delete_vertices.end();
           ++itr) {
        if (itr->v == info.v) {
          *itr = info;
          break;
        }
      }
    }
    if (is_collapsible(v2, info)) {
      for (auto itr = delete_vertices.begin(); itr != delete_vertices.end();
           ++itr) {
        if (itr->v == info.v) {
          *itr = info;
          break;
        }
      }
    }
  }

  // Remove overlapping
  resolve_overlapping_cells(mesh, camera, parallel_offset);
  mesh.garbage_collection(); // Collect here so we keep the same indexing for
                             // resolving above

  // Fill properties: Patch on face
  update_patch(mesh, camera);
  mesh.markPatchBoundaryEdges(false);

  return true;
}
