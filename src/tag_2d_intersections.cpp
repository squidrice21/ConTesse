// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "tag_2d_intersections.h"
#include <spdlog/fmt/ostr.h>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common.h"
#include "inflat_path.h"
#include "ray_cast_QI.h"

void walk_on_original_intersection(Mesh const &mesh, int id, Vertex src_orig,
                                   Vector3f const &cut_p,
                                   Vector3f const &cut_norm,
                                   Camera const &camera,
                                   Vertex &intersection_v) {
  intersection_v = Vertex();

  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  contess_assert_msg(
      patchID, "walk_on_original_intersection: Extracted patches is needed.");
  auto intersection_2d = mesh.get_vertex_property<int>("v:intersection_2d");

  contess_assert_msg(
      intersection_2d,
      "walk_on_original_intersection: 2D intersections are not created.");

  std::unordered_set<int> visited_face_indices;
  std::queue<Face> face_frontier;
  Vector2i viewport = Vector2i(camera.vpWidth(), camera.vpHeight());

  // Initialize with all neighboring faces
  {
    auto hit = mesh.halfedges(src_orig), hit_end = hit;
    do {
      // Check adjacent faces
      auto adj_f_h = mesh.opposite_halfedge(*hit);
      Face adj_f = mesh.face(adj_f_h);

      if (adj_f.is_valid() && patchID[adj_f] != id)
        continue;

      // In case we are at the boundary of the model
      if (adj_f.is_valid())
        face_frontier.push(adj_f);
    } while (++hit != hit_end);
  }

  if (face_frontier.empty())
    return;

  auto is_destination = [&](Face f) -> bool {
    Vertex vv[3];
    mesh.verticesOfFace(f, vv);

    for (auto const &v : vv) {
      // Check if the two overlaps in 2D
      if (intersection_2d[v] == src_orig.idx()) {
        intersection_v = v;
        return true;
      }
    }

    return false;
  };

  // Walk
  std::unordered_map<int, int> prev_idx;
  std::unordered_map<int, Vector3f> f_intersection;
  do {
    Face cur_f = face_frontier.front();
    face_frontier.pop();

    // Visited?
    if (visited_face_indices.find(cur_f.idx()) != visited_face_indices.end())
      continue;

    // Is this face intersecting with the cut plane?
    Vertex vertices[3];
    mesh.verticesOfFace(cur_f, vertices);
    std::vector<Vector3f> tri_p({vpositions[vertices[0]],
                                 vpositions[vertices[1]],
                                 vpositions[vertices[2]]});

    // March along the intersected triangles
    if (planeIntersectTri(cut_p, cut_norm, tri_p)) {
      // Extend frontier
      auto fhit = mesh.halfedges(cur_f), fhit_end = fhit;
      do {
        // Sample on the edges of the current face
        // and determine the intersecting edge for the next step
        Vector3f intersection;
        bool edge_intersected = planeIntersectEdge(
            cut_p, cut_norm, vpositions[mesh.from_vertex(*fhit)],
            vpositions[mesh.to_vertex(*fhit)], intersection);

        auto adj_f_h = mesh.opposite_halfedge(*fhit);
        Face adj_f = mesh.face(adj_f_h);
        // Don't grow the non-intersecting edge
        // When the face is not flipped
        if (!edge_intersected) {
          if (adj_f.is_valid()) {
            Vector3f f_adj_norm = mesh.compute_face_normal(adj_f);
            Vector3f f_norm = mesh.compute_face_normal(cur_f);

            // Trivial check for flipped surface triangles
            if (f_adj_norm.dot(f_norm) > 0)
              continue;
          } else {
            continue;
          }
        }

        f_intersection[cur_f.idx()] = intersection;

        // We can't cross the contour edge
        if (is_contour[mesh.edge(*fhit)] >= 0)
          continue;

        // Check adjacent faces
        if (visited_face_indices.find(adj_f.idx()) !=
            visited_face_indices.end())
          continue;

        // In case we are at the boundary of the model
        if (adj_f.is_valid()) {
          face_frontier.push(adj_f);
          prev_idx[adj_f.idx()] = cur_f.idx();
        }
      } while (++fhit != fhit_end);

      // If this face the destination?
      if (is_destination(cur_f)) {
        break;
      }
    }

    visited_face_indices.insert(cur_f.idx());
  } while (!face_frontier.empty());
}

void intersecting_edge_directions(Mesh &mesh, Vertex const &intersecting_v,
                                  int p_idx,
                                  std::vector<Halfedge> &inter_v_hes) {
  auto patchID = mesh.get_face_property<int>("f:patchID");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");

  auto hit = mesh.halfedges(intersecting_v), hit_end = hit;
  do {
    auto oppo_he = mesh.opposite_halfedge(*hit);
    if (!((mesh.face(*hit).is_valid() && patchID[mesh.face(*hit)] == p_idx) ||
          (mesh.face(oppo_he).is_valid() &&
           patchID[mesh.face(oppo_he)] == p_idx)))
      continue;
    if (is_contour[mesh.edge(*hit)] < 0 && !mesh.is_boundary(mesh.edge(*hit)))
      continue;
    inter_v_hes.emplace_back(*hit);
  } while (++hit != hit_end);

  contess_assert_msg(
      inter_v_hes.size() >= 2 || inter_v_hes.empty(),
      "tag_2d_intersections: Cannot find corresponding intersecting edge.");
}

// Assume the input mesh is turned into a 2D planar map
void tag_2d_intersections(Mesh &mesh, Camera const &camera) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto intersection_2d = mesh.get_vertex_property<int>("v:intersection_2d");
  contess_assert_msg(intersection_2d,
                     "tag_2d_intersections: 2D intersections are not created.");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  contess_assert_msg(patchID, "tag_2d_intersections: Patch needs to be built.");
  auto disk_cut = mesh.get_edge_property<int>("e:disk_cut");
  contess_assert_msg(disk_cut,
                     "tag_2d_intersections: Patch needs to be cut to disk.");

  auto loop_info = mesh.halfedge_property<Vector2f>("h:loop");
  if (!loop_info) {
    loop_info = mesh.add_halfedge_property<Vector2f>("h:loop");
  }
  loop_info.vector().assign(loop_info.vector().size(), -1 * Vector2f::Ones());

  // For debug visualization
  auto is_valid_intersection_2d =
      mesh.get_vertex_property<bool>("v:is_valid_intersection_2d");
  if (!is_valid_intersection_2d) {
    is_valid_intersection_2d =
        mesh.add_vertex_property<bool>("v:is_valid_intersection_2d", false);
  }
  is_valid_intersection_2d.vector().assign(
      is_valid_intersection_2d.vector().size(), false);

  std::set<std::pair<int, int>> seen_intersections;

  // Iterate on all 2D intersections and test if we can pair any two
  for (size_t i = 0; i < mesh.n_vertices(); i++) {
    Vertex v(i);

    // Not an intersection
    if (intersection_2d[v] < 0)
      continue;

    // March on surface to check if the intersection is valid (aka, bound
    // surface between them)
    // Check the one/two adjacent patches
    std::unordered_set<int> seen_patch;
    auto hit = mesh.halfedges(v), hit_end = hit;
    do {
      if (!mesh.face(*hit).is_valid())
        continue;
      if (is_contour[mesh.edge(*hit)] < 0 && !disk_cut[mesh.edge(*hit)] &&
          !mesh.is_boundary(mesh.edge(*hit)))
        continue;
      int pidx = patchID[mesh.face(*hit)];
      if (seen_patch.count(pidx))
        continue;
      seen_patch.emplace(pidx);

      // Already paired
      if (seen_intersections.count(std::make_pair(i, pidx)))
        continue;

      // We trace along the two orthogonal mid lines of the two intersecting
      // edges in 2D Assuming no three edges would intersect at the same point
      Vertex intersecting_v(intersection_2d[v]);

      // Find three directional vectors
      Halfedge v_he = *hit;
      std::vector<Halfedge> inter_v_hes;
      intersecting_edge_directions(mesh, intersecting_v, pidx, inter_v_hes);

      // The pair is not intersecting wrt this patch
      if (inter_v_hes.empty()) {
        seen_intersections.emplace(v.idx(), pidx);
        seen_intersections.emplace(intersecting_v.idx(), pidx);
        continue;
      }

      std::vector<Vector3f> path_dirs;
      auto normalized_2d = [&](Vector3f const &v1, Vector3f const &v2,
                               Vector3f &dir) {
        dir = v2 - v1;
        Vector2i viewport = Vector2i(camera.vpWidth(), camera.vpHeight());
        Vector2f v2_2d = project(v2, camera.viewMatrix().matrix(),
                                 camera.projectionMatrix(), viewport)
                             .head<2>();
        Vector2f v1_2d = project(v1, camera.viewMatrix().matrix(),
                                 camera.projectionMatrix(), viewport)
                             .head<2>();
        real_t mag_2d = (v2_2d - v1_2d).norm();
        dir /= mag_2d;
      };
      auto normalized_2d_dir = [&](Vector3f const &v1, Vector3f const &v2,
                                   Vector2f &dir) {
        Vector2i viewport = Vector2i(camera.vpWidth(), camera.vpHeight());
        Vector2f v2_2d = project(v2, camera.viewMatrix().matrix(),
                                 camera.projectionMatrix(), viewport)
                             .head<2>();
        Vector2f v1_2d = project(v1, camera.viewMatrix().matrix(),
                                 camera.projectionMatrix(), viewport)
                             .head<2>();
        real_t mag_2d = (v2_2d - v1_2d).norm();
        dir = (v2_2d - v1_2d);
        dir /= mag_2d;
      };
      Vector3f v_he_dir;
      normalized_2d(vpositions[mesh.from_vertex(v_he)],
                    vpositions[mesh.to_vertex(v_he)], v_he_dir);
      std::vector<real_t> path_dir_angles;
      for (size_t j = 0; j < 2; j++) {
        Halfedge inter_v_he = inter_v_hes[j];
        Vector3f inter_he_dir;
        normalized_2d(vpositions[mesh.from_vertex(inter_v_he)],
                      vpositions[mesh.to_vertex(inter_v_he)], inter_he_dir);
        path_dirs.emplace_back(v_he_dir + inter_he_dir);

        Vector2f v_he_dir2, inter_he_dir2;
        normalized_2d_dir(vpositions[mesh.from_vertex(v_he)],
                          vpositions[mesh.to_vertex(v_he)], v_he_dir2);
        normalized_2d_dir(vpositions[mesh.from_vertex(inter_v_he)],
                          vpositions[mesh.to_vertex(inter_v_he)],
                          inter_he_dir2);
        path_dir_angles.emplace_back(v_he_dir2.dot(inter_he_dir2));
      }

      for (size_t j = 0; j < path_dirs.size(); j++) {
        if (path_dir_angles[j] > 0)
          continue;
        auto path_dir = path_dirs[j];
        Vertex path_intersection_v;
        Vector3f cut_p = camera.position();
        Vector3f path_norm =
            (camera.position() - vpositions[v]).cross(path_dir).normalized();
        walk_on_original_intersection(mesh, pidx, v, cut_p, path_norm, camera,
                                      path_intersection_v);
        // Paired
        if (path_intersection_v.is_valid()) {
          loop_info[*hit][0] = path_intersection_v.idx();
          auto hit_p = mesh.halfedges(path_intersection_v), hit_p_end = hit_p;
          do {
            if (!mesh.face(*hit_p).is_valid())
              continue;
            if (!is_contour[mesh.edge(*hit_p)] && !disk_cut[mesh.edge(*hit_p)])
              continue;
            if (patchID[mesh.face(*hit_p)] == pidx) {
              loop_info[*hit_p][0] = v.idx();
              break;
            }
          } while (++hit_p != hit_p_end);

          seen_intersections.emplace(v.idx(), pidx);
          seen_intersections.emplace(path_intersection_v.idx(), pidx);

          is_valid_intersection_2d[v] = true;
          is_valid_intersection_2d[path_intersection_v] = true;

          break;
        }
      }
    } while (++hit != hit_end);
  }
}
