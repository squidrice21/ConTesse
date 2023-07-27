// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "ray_cast_QI.h"
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <igl/predicates/predicates.h>

#include "common.h"
#include "logger.h"
#include "tag_concave_edges.h"

void first_hit_facing_types(Mesh const &mesh, real_t mid_point_t,
                            std::vector<uint32_t> &indices,
                            std::vector<real_t> &ts, FacingType &vbo,
                            FacingType &vbo_f) {
  contess_assert(!ts.empty() && ts.size() == indices.size());
  auto VBO_f = mesh.get_face_property<FacingType>("f:VBO_f");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");

  vbo = FacingType::NA;
  vbo_f = FacingType::NA;

  // Remove the intersection with the query point
  real_t min_t = 1e3;
  int min_i = -1;
  for (size_t i = 0; i < indices.size(); i++) {
    if (ts[i] < min_t && std::abs(mid_point_t - ts[i]) > EPSILON) {
      min_t = ts[i];
      min_i = i;
    }
  }
  if (min_i >= 0) {
    vbo = VBO[Face(indices[min_i])];
    vbo_f = VBO_f[Face(indices[min_i])];
  }
}

// Trim occluders based by whether it's straightly in front of the contour
// mid point (So the triangle adjacent to the contour is also not included)
void trim_occluders(Mesh const &mesh, Edge const &e, real_t mid_point_t,
                    std::vector<uint32_t> &indices, std::vector<real_t> &ts) {
  // First trim the adjacent triangle
  std::unordered_set<int> adj_faces;
  adj_faces.insert(mesh.face(e, 0).idx());
  adj_faces.insert(mesh.face(e, 1).idx());

  std::vector<size_t> to_trim;
  to_trim.reserve(2);
  for (size_t i = 0; i < indices.size(); i++) {
    if (adj_faces.find(indices[i]) != adj_faces.end()) {
      to_trim.emplace_back(i);
    }
  }
  std::reverse(to_trim.begin(), to_trim.end());
  for (auto trim_i : to_trim) {
    indices.erase(indices.begin() + trim_i);
    ts.erase(ts.begin() + trim_i);
  }

  // Trim based on the time
  // Also trim if the intersection is too close to the midpoint but is not in an
  // adjacent triangle to avoid double counting C based on the adjacent
  // triangles.
  to_trim.clear();
  for (size_t i = 0; i < ts.size(); i++) {
    if (ts[i] >= mid_point_t || std::abs(ts[i] - mid_point_t) < EPSILON) {
      to_trim.emplace_back(i);
    }
  }
  std::reverse(to_trim.begin(), to_trim.end());
  for (auto trim_i : to_trim) {
    indices.erase(indices.begin() + trim_i);
    ts.erase(ts.begin() + trim_i);
  }

  // Trim the faces from back-facing patch
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");
  to_trim.clear();
  for (size_t i = 0; i < indices.size(); i++) {
    if (VBO[Face(indices[i])] == FacingType::BACK) {
      to_trim.emplace_back(i);
    }
  }
  std::reverse(to_trim.begin(), to_trim.end());
  for (auto trim_i : to_trim) {
    indices.erase(indices.begin() + trim_i);
    ts.erase(ts.begin() + trim_i);
  }
}

void get_nearer_adjacent_face(Mesh const &mesh, Camera const &camera,
                              Edge const &edge, Face &nearer_f,
                              Vertex &nearer_v) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");

  surface_mesh::Surface_mesh::Halfedge h = mesh.halfedge(edge, 0);

  // Build a triangle which has the contour edge as a side and one side
  // parallel to the image plane (necessary?)
  const Vector3f &a = vpositions[mesh.from_vertex(h)];
  const Vector3f &b = vpositions[mesh.to_vertex(h)];

  Vector3f ab = b - a;
  Vector3f ac = camera.position() - a;
  Vector3f on_image_plane = ac.cross(ab);
  on_image_plane.normalize();
  Vector3f e = a + on_image_plane;

  // Find the adjacent vertex that is on the front portion
  Face f1 = mesh.face(edge, 0);
  Face f2 = mesh.face(edge, 1);

  Vertex f1_v = mesh.to_vertex(mesh.next_halfedge(h));
  Vertex f2_v = mesh.to_vertex(mesh.next_halfedge(mesh.opposite_halfedge(h)));

  Vector3f p1 = vpositions[f1_v];
  Vector3f p2 = vpositions[f2_v];

  nearer_f = f1;
  nearer_v = f1_v;

  // The front portion vertex can be at the degenerate position
  // (almost colinear with the contour edge)
  // In this case, test the back front face
  if (vertex_to_edge_distance(a, b, p1) < FRONT_PORTION_VERTEX_EPSILON) {
    logger().info("Degenerate contour opposite vertex: {}",
                  vertex_to_edge_distance(a, b, p1));
    if ((igl::predicates::orient3d(a.cast<double>(), b.cast<double>(),
                                   e.cast<double>(), p2.cast<double>()) ==
         igl::predicates::Orientation::POSITIVE) ==
        (igl::predicates::orient3d(a.cast<double>(), b.cast<double>(),
                                   e.cast<double>(),
                                   camera.position().cast<double>()) ==
         igl::predicates::Orientation::POSITIVE)) {
      nearer_f = f2;
      nearer_v = f2_v;
    }
  }

  bool p1_cam_side =
      (igl::predicates::orient3d(a.cast<double>(), b.cast<double>(),
                                 e.cast<double>(), p1.cast<double>()) ==
       igl::predicates::Orientation::POSITIVE) ==
      (igl::predicates::orient3d(a.cast<double>(), b.cast<double>(),
                                 e.cast<double>(),
                                 camera.position().cast<double>()) ==
       igl::predicates::Orientation::POSITIVE);
  bool p2_cam_side =
      (igl::predicates::orient3d(a.cast<double>(), b.cast<double>(),
                                 e.cast<double>(), p2.cast<double>()) ==
       igl::predicates::Orientation::POSITIVE) ==
      (igl::predicates::orient3d(a.cast<double>(), b.cast<double>(),
                                 e.cast<double>(),
                                 camera.position().cast<double>()) ==
       igl::predicates::Orientation::POSITIVE);
  igl::predicates::Orientation p1_side = igl::predicates::orient3d(
      a.cast<double>(), b.cast<double>(), e.cast<double>(), p1.cast<double>());
  igl::predicates::Orientation p2_side = igl::predicates::orient3d(
      a.cast<double>(), b.cast<double>(), e.cast<double>(), p2.cast<double>());

  // If the vertex on the front portion is on the same side of the camera
  // then the edge is convex

  if ((igl::predicates::orient3d(a.cast<double>(), b.cast<double>(),
                                 e.cast<double>(), p1.cast<double>()) ==
       igl::predicates::Orientation::POSITIVE) !=
      (igl::predicates::orient3d(a.cast<double>(), b.cast<double>(),
                                 e.cast<double>(),
                                 camera.position().cast<double>()) ==
       igl::predicates::Orientation::POSITIVE)) {
    nearer_f = f2;
    nearer_v = f2_v;
  }

  auto VBO_f = mesh.get_face_property<FacingType>("f:VBO_f");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");
  auto f_norm = mesh.compute_face_normal(nearer_f);
  auto c_ray = vpositions[nearer_v] - camera.position();
  c_ray.normalized();
  FacingType face_facing =
      (f_norm.dot(c_ray) < 0) ? FacingType::FRONT : FacingType::BACK;

  logger().debug("\tF1: {}, v1: {}; F2: {}, v2: {}; ", f1.idx(), f1_v.idx(),
                 f2.idx(), f2_v.idx());
  logger().debug(
      "\tF1 cam side: {}, f1 side: {}; F2 cam side: {}, f2 side: {}; ",
      p1_cam_side, p1_side, p2_cam_side, p2_side);
  logger().debug("\tNear face: {} out of {}, {}; Facing: VBO: {}; VBO_f: {}; "
                 "Face facing: {}",
                 nearer_f.idx(), f1.idx(), f2.idx(), VBO[nearer_f],
                 VBO_f[nearer_f], face_facing);
}

bool is_nearer_adjacent_face_consistent(Mesh const &mesh, Camera const &camera,
                                        Edge const &edge) {
  Face nearer_f;
  Vertex nearer_v;
  get_nearer_adjacent_face(mesh, camera, edge, nearer_f, nearer_v);

  auto VBO_f = mesh.get_face_property<FacingType>("f:VBO_f");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");

  // If the face belongs to a back facing patch, then we ignore it
  if (VBO[nearer_f] == FacingType::BACK)
    return true;

  return VBO[nearer_f] == VBO_f[nearer_f];
}

bool planeIntersectTri(const Vector3f &pp, const Vector3f &norm,
                       std::vector<Vector3f> const &tri_p) {
  // Check if corners are on different sides / on the plane.
  Vector3f pp2 = norm.cross(Vector3f(1, 0, 0));
  Vector3f pp3 = norm.cross(pp2);
  pp2 = pp + pp2;
  pp3 = pp + pp3;
  bool has_seen_pos = false, has_seen_neg = false;
  for (uint8_t i = 0; i < tri_p.size(); i++) {
    Vector3f corner = tri_p[i];
    auto side = igl::predicates::orient3d(pp, pp2, pp3, corner);
    switch (side) {
    case igl::predicates::Orientation::NEGATIVE:
      has_seen_neg = true;
      break;
    case igl::predicates::Orientation::POSITIVE:
      has_seen_pos = true;
      break;
    default: // Orientation::COPLANAR
      has_seen_neg = true;
      has_seen_pos = true;
      break;
    }

    if (has_seen_pos && has_seen_neg) {
      break;
    }
  }

  return has_seen_pos && has_seen_neg;
}

int ray_cast_contour(Mesh const &mesh, Camera const &camera, Edge const &e,
                     bool is_trivial) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto VBO_f = mesh.get_face_property<FacingType>("f:VBO_f");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");

  Vector3f contour_mid_point =
      (vpositions[mesh.vertex(e, 0)] + vpositions[mesh.vertex(e, 1)]) * 0.5;
  Vector3f dir = contour_mid_point - camera.position();
  real_t mid_point_t = dir.norm();
  Ray ray = Ray(camera.position(), dir.normalized());

  contess_assert_msg(mesh.get_bvh(), "BVH missing.");

  std::unordered_set<uint32_t> adj_faces{(uint32_t)mesh.face(e, 0).idx(),
                                         (uint32_t)mesh.face(e, 1).idx()};

  std::vector<uint32_t> indices;
  std::vector<real_t> ts;
  mesh.get_bvh()->rayIntersect(ray, indices, ts, nullptr);

  if (!is_trivial)
    logger().debug("Edge: {}, {}; E_idx: {}; t: {}", mesh.vertex(e, 0).idx(),
                   mesh.vertex(e, 1).idx(), e.idx(), mid_point_t);

  // Trivial strategy: directly use the number of intersections before hitting
  // the mid-point (included)
  if (is_trivial) {
    int pre_count = 0;

    real_t min_t = mid_point_t;
    for (size_t i = 0; i < ts.size(); ++i) {
      if (adj_faces.count(indices[i]))
        min_t = std::min(min_t, ts[i]);
    }

    std::unordered_set<real_t> times;
    for (size_t i = 0; i < ts.size(); ++i) {
      // Avoid the two adjacent faces
      if (adj_faces.count(indices[i]) || indices[i] >= mesh.n_faces())
        continue;
      times.insert(ts[i]);
    }

    for (auto t : times) {
      if (t < min_t && (min_t - t) > 1e-10)
        pre_count++;
    }

    return pre_count;
  } else { // Advanced strategy: subtract the inconsistent intersections
    auto patchID = mesh.get_face_property<int>("f:patchID");
    contess_assert_msg(
        patchID,
        "Patches are not determined. Can't run the advanced QI computation.");

    logger().debug("\tBefore trimming");
    for (size_t i = 0; i < indices.size(); i++) {
      logger().debug("\t\tF: {}, t: {}", indices[i], ts[i]);
    }

    // Trim occluders based by whether it's straightly in front of the contour
    // mid point (So the triangle adjacent to the contour is also not included)
    trim_occluders(mesh, e, mid_point_t, indices, ts);

    logger().debug("\t--------------");
    logger().debug("\tAfter trimming");
    for (size_t i = 0; i < indices.size(); i++) {
      logger().debug("\t\tF: {}, t: {}", indices[i], ts[i]);
    }

    std::unordered_map<int, std::pair<int, int>> patch_counts;
    for (size_t i = 0; i < indices.size(); i++) {
      // Increment the counters based on whether the intersection is consistent
      Face f(indices[i]);

      // Contour can't be occluded by a back facing patch
      if (VBO[f] == FacingType::BACK)
        continue;

      int p_id = patchID[f];
      if (patch_counts.find(p_id) == patch_counts.end())
        patch_counts[p_id] = std::make_pair(0, 0);
      if (VBO[f] == VBO_f[f]) {
        // When the intersection is consistent,
        patch_counts[p_id].first++;
      } else
        // When the intersection is inconsistent,
        patch_counts[p_id].second++;
    }

    int pre_count = 0;
    for (auto p_count : patch_counts) {
      int p_total = p_count.second.first - p_count.second.second;

      pre_count += p_total;

      logger().debug("\tP: {}; Consis: {}; Inconsis: {}", p_count.first,
                     p_count.second.first, p_count.second.second);
    }

    // Correct the count based on the adjacent face
    size_t adj_inconsist_count = 0;
    bool nearer_consistent =
        is_nearer_adjacent_face_consistent(mesh, camera, e);
    if (!nearer_consistent)
      adj_inconsist_count = 1;

    pre_count -= adj_inconsist_count;
    logger().debug("\tC = {}", adj_inconsist_count);
    logger().debug("\t=> QI: {}", pre_count);
    logger().debug("================================");

    return pre_count;
  }
}

void ray_cast_QI(Mesh &mesh, Camera const &camera, bool is_trivial) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");

  contess_assert_msg(is_contour,
                     "Contours are missing. Can't determine QIs for contours.");

  // Initialize all QIs on edges to be -1 (the default values for non-contour
  // edges)
  auto edge_qi = mesh.edge_property<int>("e:qi");

  if (!edge_qi) {
    edge_qi = mesh.add_edge_property<int>("e:qi");
  }
  edge_qi.vector().assign(edge_qi.vector().size(), -1);

  // Build BVH if non-existing
  if (mesh.get_bvh() == nullptr)
    mesh.build_bvh();

  // Compute QIs on contour edges
  parallel_for(
      tbb::blocked_range<uint32_t>(0u, (uint32_t)mesh.n_edges(), GRAIN_SIZE),
      [&](const tbb::blocked_range<uint32_t> &range) {
        for (uint32_t i = range.begin(); i != range.end(); ++i) {
          Edge e = Edge(i);
          if (is_contour[e] < 0 && !mesh.is_boundary(e))
            continue;
          edge_qi[e] = ray_cast_contour(mesh, camera, e, is_trivial);
        }
      });
}
