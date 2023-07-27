// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "remove_cut_infeasible_disk.h"
#include "common.h"
#include "tag_simplification_edges.h"
#include <deque>
#include <limits>
#include <spdlog/fmt/ostr.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

void check_cut_feasibility(Mesh &mesh, std::vector<Edge> &invalid_cut_edges) {
  auto cut_candidate = mesh.get_edge_property<bool>("e:cut_candidate");
  contess_assert_msg(cut_candidate,
                     "check_cut_feasibility: Needs to label candidate edges.");
  auto cut = mesh.get_edge_property<int>("e:disk_cut");
  contess_assert_msg(cut, "check_cut_feasibility: Needs to have cuts.");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  contess_assert_msg(patchID, "check_cut_feasibility: Missing patches.");

  for (size_t i = 0; i < mesh.n_edges(); i++) {
    Edge e(i);

    // Skip boundary (manifold/patch) edges
    if (!mesh.face(e, 1).is_valid() || !mesh.face(e, 0).is_valid() ||
        (mesh.face(e, 0).is_valid() && mesh.face(e, 1).is_valid() &&
         patchID[mesh.face(e, 0)] != patchID[mesh.face(e, 1)]))
      continue;

    if (cut[e] < 0)
      continue;

    // Check if it's adjacent to an inconsistent face
    if (!cut_candidate[e]) {
      logger().error("check_cut_feasibility: edge adjacent to an inconsistent "
                     "face: {}, {}, {} in p{}",
                     e, mesh.vertex(e, 0), mesh.vertex(e, 1),
                     patchID[mesh.face(e, 0)]);
      invalid_cut_edges.emplace_back(e);
    }
  }
}

bool is_patch_removeable(Mesh const &mesh, Camera const &camera,
                         Halfedge const &in_he, real_t max_loop_length) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  auto cut = mesh.get_edge_property<int>("e:disk_cut");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");

  contess_assert(mesh.face(in_he).is_valid());

  // Walk along patch boundary and check if there's any cut
  real_t loop_length = 0;
  std::unordered_set<int> seen_vertices;
  Halfedge he = in_he;
  Vertex current_v = mesh.from_vertex(he);
  auto patchBoundary = mesh.get_halfedge_property<int>("h:patchBoundary");

  auto seen_cut = [&](Vertex const &v) -> bool {
    auto hit = mesh.halfedges(v), hend = hit;
    do {
      Edge hit_e = mesh.edge(*hit);

      // Reaches a contour
      if (is_contour[hit_e] >= 0.f || mesh.is_boundary(hit_e)) {
        continue;
      }

      if (cut[hit_e] == patchID[mesh.face(in_he)])
        return true;
    } while (++hit != hend);
    return false;
  };
  do {
    current_v = mesh.from_vertex(he);
    Vertex next_v = mesh.to_vertex(he);

    {
      Vector2f v1 =
          project(vpositions[mesh.from_vertex(he)],
                  camera.viewMatrix().matrix(), camera.projectionMatrix(),
                  Vector2i(camera.vpWidth(), camera.vpHeight()))
              .head(2);
      Vector2f v2 =
          project(vpositions[mesh.to_vertex(he)], camera.viewMatrix().matrix(),
                  camera.projectionMatrix(),
                  Vector2i(camera.vpWidth(), camera.vpHeight()))
              .head(2);
      loop_length += (v2 - v1).norm();
    }

    // Check
    if (seen_cut(next_v))
      return false;

    seen_vertices.emplace(current_v.idx());
    auto hit = mesh.halfedges(next_v), hend = hit;
    do {
      Edge hit_e = mesh.edge(*hit);

      if ((seen_vertices.count(mesh.to_vertex(*hit).idx()) &&
           mesh.to_vertex(*hit) != mesh.from_vertex(in_he)) ||
          (mesh.from_vertex(*hit) == mesh.to_vertex(in_he) &&
           mesh.to_vertex(*hit) == mesh.from_vertex(in_he)))
        continue;

      if (!is_same_patch(mesh, in_he, Halfedge(), *hit))
        continue;

      // Reaches a contour
      if (is_contour[hit_e] >= 0.f || mesh.is_boundary(hit_e)) {
        he = *hit;
        break;
      }
    } while (++hit != hend);
  } while (mesh.to_vertex(he) != mesh.from_vertex(in_he));

  Vertex next_v = mesh.to_vertex(he);
  current_v = mesh.from_vertex(he);
  {
    Vector2f v1 =
        project(vpositions[mesh.from_vertex(he)], camera.viewMatrix().matrix(),
                camera.projectionMatrix(),
                Vector2i(camera.vpWidth(), camera.vpHeight()))
            .head(2);
    Vector2f v2 =
        project(vpositions[mesh.to_vertex(he)], camera.viewMatrix().matrix(),
                camera.projectionMatrix(),
                Vector2i(camera.vpWidth(), camera.vpHeight()))
            .head(2);
    loop_length += (v2 - v1).norm();
  }

  // Check
  if (seen_cut(next_v))
    return false;

  // Check loop length
  logger().info("is_patch_removeable: p{}, {} vs {}", patchID[mesh.face(in_he)],
                loop_length, max_loop_length);

  if (loop_length > max_loop_length) {
    return false;
  }

  return true;
}

void remove_patch(Mesh &mesh, Halfedge const &in_he) {
  auto patch_removed = mesh.edge_property<bool>("e:patch_removed");
  auto patchID = mesh.face_property<int>("f:patchID");
  auto VBO = mesh.face_property<FacingType>("f:VBO");
  auto is_contour = mesh.edge_property<real_t>("e:contour");
  auto concave = mesh.edge_property<bool>("e:concave");

  auto is_cusp = mesh.vertex_property<bool>("v:cusp");
  auto cusp_facing = mesh.vertex_property<FacingType>("v:cusp_facing");
  auto intersection_2d = mesh.get_vertex_property<int>("v:intersection_2d");
  if (intersection_2d)
    intersection_2d = mesh.vertex_property<int>("v:intersection_2d");
  auto is_valid_intersection_2d =
      mesh.get_vertex_property<bool>("v:is_valid_intersection_2d");
  if (is_valid_intersection_2d)
    is_valid_intersection_2d =
        mesh.vertex_property<bool>("v:is_valid_intersection_2d");

  auto chained_orientation =
      mesh.get_halfedge_property<int>("h:chained_orientation");
  if (chained_orientation)
    chained_orientation = mesh.halfedge_property<int>("h:chained_orientation");
  auto patchBoundary = mesh.get_halfedge_property<int>("h:patchBoundary");
  if (patchBoundary)
    patchBoundary = mesh.halfedge_property<int>("h:patchBoundary");
  auto cut_directed = mesh.get_halfedge_property<bool>("h:cut");
  if (cut_directed)
    cut_directed = mesh.halfedge_property<bool>("h:cut");

  contess_assert(!mesh.is_boundary(mesh.edge(in_he)));
  contess_assert(mesh.face(in_he).is_valid());
  contess_assert(mesh.face(mesh.opposite_halfedge(in_he)).is_valid());
  int patch_id = patchID[mesh.face(in_he)];
  int to_patch_id = patchID[mesh.face(mesh.opposite_halfedge(in_he))];

  // 1. Remove contour (edges, halfedges and vertices)
  std::vector<Halfedge> remove_contour_he;
  std::vector<Halfedge> boundary_he;
  std::unordered_set<int> seen_vertices;
  Halfedge he = in_he;
  Vertex current_v = mesh.from_vertex(he);

  remove_contour_he.emplace_back(in_he);

  do {
    current_v = mesh.from_vertex(he);
    Vertex next_v = mesh.to_vertex(he);

    seen_vertices.emplace(current_v.idx());
    auto hit = mesh.halfedges(next_v), hend = hit;
    do {
      Edge hit_e = mesh.edge(*hit);

      if ((seen_vertices.count(mesh.to_vertex(*hit).idx()) &&
           mesh.to_vertex(*hit) != mesh.from_vertex(in_he)) ||
          (mesh.from_vertex(*hit) == mesh.to_vertex(in_he) &&
           mesh.to_vertex(*hit) == mesh.from_vertex(in_he)))
        continue;

      if (!is_same_patch(mesh, in_he, Halfedge(), *hit))
        continue;

      // Reaches a contour
      if (is_contour[hit_e] >= 0.f || mesh.is_boundary(hit_e)) {
        if (is_contour[hit_e] >= 0 && !mesh.is_boundary(hit_e)) {
          remove_contour_he.emplace_back(*hit);
        } else {
          boundary_he.emplace_back(*hit);
        }
        he = *hit;
        break;
      }
    } while (++hit != hend);
  } while (mesh.to_vertex(he) != mesh.from_vertex(in_he));

  for (auto const &he : remove_contour_he) {
    Vertex v = mesh.from_vertex(he);
    is_cusp[v] = false;
    cusp_facing[v] = FacingType::UNDEFINED;
    if (intersection_2d)
      intersection_2d[v] = -1;
    if (is_valid_intersection_2d)
      is_valid_intersection_2d[v] = false;

    std::vector<Halfedge> hes({he, mesh.opposite_halfedge(he)});
    for (auto const &hh : hes) {
      if (!hh.is_valid())
        continue;
      if (chained_orientation) {
        chained_orientation[hh] = -1;
        chained_orientation[mesh.opposite_halfedge(hh)] = -1;
      }
      if (patchBoundary) {
        patchBoundary[hh] = -1;
        patchBoundary[mesh.opposite_halfedge(hh)] = -1;
      }
      if (cut_directed) {
        cut_directed[hh] = false;
        cut_directed[mesh.opposite_halfedge(hh)] = false;
      }
    }

    is_contour[mesh.edge(he)] = -1;
    concave[mesh.edge(he)] = false;
  }

  for (auto const &he : boundary_he) {
    Halfedge oppo = mesh.opposite_halfedge(he);
    if (patchBoundary) {
      patchBoundary[he] = to_patch_id;
      patchBoundary[oppo] = to_patch_id;

      chained_orientation[he] = -1;
      chained_orientation[oppo] = -1;

      cut_directed[he] = false;
      cut_directed[oppo] = false;
    }
  }

  // 2. Remove faces
  for (auto f_itr = mesh.faces_begin(); f_itr != mesh.faces_end(); ++f_itr) {
    if (patchID[*f_itr] == patch_id) {
      patchID[*f_itr] = to_patch_id;
      VBO[*f_itr] = (VBO[*f_itr] == FacingType::FRONT) ? FacingType::BACK
                                                       : FacingType::FRONT;

      // Label the adjacent edges as removed edges so we don't cut through them
      if (patch_removed) {
        auto hit = mesh.halfedges(*f_itr), hend = hit;
        do {
          Edge hit_e = mesh.edge(*hit);
          patch_removed[hit_e] = true;
        } while (++hit != hend);
      }
    }
  }
}

bool remove_cut_infeasible_disk(Mesh &mesh, Camera const &camera) {
  bool removed = false;

  std::vector<Edge> invalid_cut_edges;
  check_cut_feasibility(mesh, invalid_cut_edges);

  auto patchID = mesh.get_face_property<int>("f:patchID");
  auto cut = mesh.edge_property<int>("e:disk_cut");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");

  std::vector<Halfedge> patch_to_remove;
  std::unordered_set<int> seen_removal;
  for (auto const &e : invalid_cut_edges) {
    // 1. Find closest patch
    std::unordered_map<int, real_t> contour_he_distance;
    {
      std::queue<Edge> edge_queue;
      std::unordered_set<int> seen_vertices;
      std::unordered_map<int, real_t> vertex_distance;
      edge_queue.push(e);
      edge_queue.push(e);
      seen_vertices.emplace(mesh.vertex(e, 0).idx());
      vertex_distance[mesh.vertex(e, 0).idx()] = 0;
      vertex_distance[mesh.vertex(e, 1).idx()] = 0;

      while (!edge_queue.empty()) {
        Edge current_e = edge_queue.front();
        edge_queue.pop();

        Vertex current_v =
            (seen_vertices.count(mesh.vertex(current_e, 0).idx()))
                ? mesh.vertex(current_e, 1)
                : mesh.vertex(current_e, 0);
        Vertex prev_v = (seen_vertices.count(mesh.vertex(current_e, 0).idx()))
                            ? mesh.vertex(current_e, 0)
                            : mesh.vertex(current_e, 1);

        // We need to take the other direction from the beginning edge
        if (seen_vertices.count(mesh.vertex(current_e, 0).idx()) &&
            seen_vertices.count(mesh.vertex(current_e, 1).idx())) {
          current_v = mesh.vertex(current_e, 0);
          prev_v = mesh.vertex(current_e, 1);
        }
        if (!vertex_distance.count(current_v.idx()))
          vertex_distance[current_v.idx()] =
              vertex_distance[prev_v.idx()] +
              (vpositions[current_v] - vpositions[prev_v]).norm();
        seen_vertices.emplace(current_v.idx());
        auto hit = mesh.halfedges(current_v), hend = hit;
        do {
          Edge hit_e = mesh.edge(*hit);

          if (seen_vertices.count(mesh.to_vertex(*hit).idx()))
            continue;

          // Reaches a contour
          if (is_contour[hit_e] >= 0.f || mesh.is_boundary(hit_e)) {
            if (is_contour[hit_e] >= 0.f) {
              Halfedge contour_he = *hit;
              if (!mesh.face(contour_he).is_valid())
                contour_he = mesh.opposite_halfedge(contour_he);
              contess_assert_msg(
                  mesh.face(contour_he).is_valid(),
                  "remove_cut_infeasible_disk: Broken manifold/patch boundary, "
                  "v" +
                      std::to_string(mesh.from_vertex(contour_he).idx()) +
                      ", v" + std::to_string(mesh.to_vertex(contour_he).idx()));
              if (patchID[mesh.face(contour_he)] == cut[e])
                contour_he = mesh.opposite_halfedge(contour_he);
              contour_he_distance[contour_he.idx()] =
                  vertex_distance[current_v.idx()];
            }
            continue;
          }

          // Finds a cut edge in the same patch
          if (cut[hit_e] == cut[e]) {
            edge_queue.push(hit_e);
          }
        } while (++hit != hend);
      }
    }

    Halfedge nearest_patch_he;
    real_t min_patch_dist = std::numeric_limits<real_t>::infinity();
    for (auto const &he_dist : contour_he_distance) {
      Halfedge he(he_dist.first);
      if (he_dist.second < min_patch_dist) {
        min_patch_dist = he_dist.second;
        nearest_patch_he = he;
      }
    }

    contess_assert_msg(
        nearest_patch_he.is_valid(),
        "remove_cut_infeasible_disk: Cut edge not connected to any loop, v" +
            std::to_string(mesh.vertex(e, 0).idx()) + ", v" +
            std::to_string(mesh.vertex(e, 1).idx()));

    // 1.1 Check if we've already recorded this patch
    if (seen_removal.count(patchID[mesh.face(nearest_patch_he)]))
      continue;

    // 2. Check if patch can be removed (if it is a disk)
    if (is_patch_removeable(mesh, camera, nearest_patch_he)) {
      logger().info("remove_cut_infeasible_disk: To remove patch {}",
                    patchID[mesh.face(nearest_patch_he)]);
      patch_to_remove.emplace_back(nearest_patch_he);
      seen_removal.emplace(patchID[mesh.face(nearest_patch_he)]);
    } else
      logger().warn(
          "remove_cut_infeasible_disk: Cannot remove non-disk patch {}.",
          patchID[mesh.face(nearest_patch_he)]);
  }

  // 3. Remove patch
  for (auto const &he_p : patch_to_remove) {
    logger().info("remove_cut_infeasible_disk: Remove patch {}",
                  patchID[mesh.face(he_p)]);
    remove_patch(mesh, he_p);
    removed = true;
  }

  // 4. Remove chains and cuts
  if (removed) {
    mesh.get_patch_chains().clear();
    mesh.get_oriented_chains().clear();
    for (auto e_itr = mesh.edges_begin(); e_itr != mesh.edges_end(); ++e_itr) {
      cut[*e_itr] = -1;
    }
  }

  return removed;
}
