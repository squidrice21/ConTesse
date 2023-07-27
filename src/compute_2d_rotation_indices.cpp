// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "compute_2d_rotation_indices.h"

#include <igl/predicates/predicates.h>
#include <math.h>
#include <random>
#include <spdlog/fmt/ostr.h>
#include <vector>

#include "PatchRotationAngle.h"
#include "common.h"
#include "logger.h"
#include "svg.h"

bool is_f1_nearer(Mesh const &mesh, Camera const &camera, Face const &f1,
                  Face const &f2) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");

  auto get_shared_halfedge = [](Mesh const &mesh, Face const &f1,
                                Face const &f2) -> Halfedge {
    surface_mesh::Surface_mesh::Halfedge_around_face_circulator hit, hit_end;
    hit = mesh.halfedges(f1);
    hit_end = hit;

    do {
      if (mesh.face(mesh.opposite_halfedge(*hit)) == f2) {
        return *hit;
      }
    } while (++hit != hit_end);

    // Can't find the adjacent face
    logger().error("Can't find the shared edge: {} - {}.", f1.idx(), f2.idx());
    return *hit;
  };

  surface_mesh::Surface_mesh::Halfedge h = get_shared_halfedge(mesh, f1, f2);

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
  Vertex f1_v = mesh.to_vertex(mesh.next_halfedge(h));
  Vertex f2_v = mesh.to_vertex(mesh.next_halfedge(mesh.opposite_halfedge(h)));

  if (vertex_to_edge_distance(a, b, vpositions[f1_v].cast<double>()) <
      FRONT_PORTION_VERTEX_EPSILON) {
    if ((igl::predicates::orient3d(a.cast<double>(), b.cast<double>(),
                                   e.cast<double>(),
                                   vpositions[f2_v].cast<double>()) ==
         igl::predicates::Orientation::POSITIVE) ==
        (igl::predicates::orient3d(a.cast<double>(), b.cast<double>(),
                                   e.cast<double>(),
                                   camera.position().cast<double>()) ==
         igl::predicates::Orientation::POSITIVE))
      return false;
    return true;
  }

  if ((igl::predicates::orient3d(a.cast<double>(), b.cast<double>(),
                                 e.cast<double>(),
                                 vpositions[f1_v].cast<double>()) ==
       igl::predicates::Orientation::POSITIVE) ==
      (igl::predicates::orient3d(a.cast<double>(), b.cast<double>(),
                                 e.cast<double>(),
                                 camera.position().cast<double>()) ==
       igl::predicates::Orientation::POSITIVE))
    return true;
  return false;
}

void project_one_ring(Mesh const &mesh, Camera &camera, Vertex const &v,
                      FacingType facing, std::vector<PatchRotationAngle> &Ks) {
  auto project_2d = [](Camera &camera, Vector3f const &v) -> Vector2f {
    Vector2f v_2d =
        project(v, camera.viewMatrix().matrix(), camera.projectionMatrix(),
                Vector2i(camera.vpWidth(), camera.vpHeight()))
            .head(2);

    return v_2d;
  };

  auto build_2d_angle =
      [&project_2d](Camera &camera, Vector3f const &w1, Vector3f const &v,
                    Vector3f const &w2) -> PatchRotationAngle {
    Vector2f w1_2d = project_2d(camera, w1);
    Vector2f v_2d = project_2d(camera, v);
    Vector2f w2_2d = project_2d(camera, w2);

    PatchRotationAngle K(w1_2d, v_2d, w2_2d);
    return K;
  };

  struct PatchRotationAngleDepth {
    PatchRotationAngle K;
    Face f;
    bool is_consistent;
  };

  auto VBO = mesh.get_face_property<FacingType>("f:VBO");
  auto VBO_f = mesh.get_face_property<FacingType>("f:VBO_f");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");

  bool is_first_round = true;
  std::vector<PatchRotationAngleDepth> Ks_depth;
  Ks.clear();
  std::mt19937 g_rand(0);

  auto is_invalid = [](Vector2f const &pos) -> bool {
    return pos.array().isNaN().any() || pos.array().isInf().any();
  };

  auto hit = mesh.halfedges(v);
  auto hit_end = hit;
  do {
    Face f = mesh.face(*hit);
    if (VBO[f] != facing)
      continue;

    // CCW move to the next half edge
    auto hit_next = hit;
    ++hit_next;

    // Project to 2D and compute the signed angle (CCW: >0)
    auto K = build_2d_angle(camera, vpositions[mesh.to_vertex(*hit)],
                            vpositions[mesh.from_vertex(*hit)],
                            vpositions[mesh.to_vertex(*hit_next)]);
    PatchRotationAngleDepth k_depth;
    k_depth.K = K;
    k_depth.f = f;
    k_depth.is_consistent = VBO[f] == VBO_f[f];
    Ks_depth.emplace_back(k_depth);

    if (K.ori == 0 || is_invalid(K.w1) || is_invalid(K.w2)) {
      if (is_first_round)
        logger().warn(
            "Degenerated projected triangles caused by non-generic camera");

      auto w1 = vpositions[mesh.to_vertex(*hit)];
      auto v = vpositions[mesh.from_vertex(*hit)];
      auto w2 = vpositions[mesh.to_vertex(*hit_next)];
    }
  } while (++hit != hit_end);

  // Order the Ks to go from one side of the patch to the other side
  {
    int k_start = 0;
    for (; k_start < (int)Ks_depth.size(); k_start++) {
      int k_pre = (k_start - 1 + Ks_depth.size()) % Ks_depth.size();
      auto K_start = Ks_depth[k_start].K;
      auto K_pre = Ks_depth[k_pre].K;
      if (K_pre.w2 != K_start.w1) {
        break;
      }
    }

    std::vector<PatchRotationAngleDepth> Ks_depth_reorder;
    Ks_depth_reorder.reserve(Ks_depth.size());
    int i = k_start;
    do {
      Ks_depth_reorder.emplace_back(Ks_depth[i]);
      i++;
      i = i % Ks_depth.size();
    } while (i != k_start);

    Ks_depth = Ks_depth_reorder;
  }

  // Determine the depth order of the inconsistent range
  std::vector<std::pair<bool, NearFarType>> near_far_record;
  near_far_record.reserve(Ks_depth.size());
  for (size_t i = 0; i < Ks_depth.size(); i++) {
    Face f = Ks_depth[i].f;
    bool is_consistent = (VBO[f] == VBO_f[f]);

    bool is_border_inconsistent = (!is_consistent);
    int prev = i - 1, next = i + 1;
    bool is_consistent_prev =
        (prev <= 0) ? true : (VBO[Ks_depth[prev].f] == VBO_f[Ks_depth[prev].f]);
    bool is_consistent_next =
        (next >= (int)Ks_depth.size())
            ? true
            : (VBO[Ks_depth[next].f] == VBO_f[Ks_depth[next].f]);

    if (is_border_inconsistent) {
      if (!(i + 1 == Ks_depth.size() || i == 0)) {
        if (!is_consistent_prev && !is_consistent_next)
          is_border_inconsistent = false;
      }
    }

    NearFarType near_far = NearFarType::NEAR;
    if (is_border_inconsistent) {
      if (is_consistent_prev && prev >= 0) {
        near_far = (is_f1_nearer(mesh, camera, f, Ks_depth[prev].f))
                       ? NearFarType::NEAR
                       : NearFarType::FAR;
      }
      if (near_far == NearFarType::NEAR && is_consistent_next &&
          next < (int)Ks_depth.size()) {
        near_far = (is_f1_nearer(mesh, camera, f, Ks_depth[next].f))
                       ? NearFarType::NEAR
                       : NearFarType::FAR;
      }
    }
    near_far_record.emplace_back(is_border_inconsistent, near_far);
  }

  // Read out the correct orientations
  auto look_up_near_far = [&](size_t i) -> NearFarType {
    NearFarType nf_type = NearFarType::NEAR;
    for (size_t j = i; j >= 0; j--) {
      if (near_far_record[j].first) {
        nf_type = near_far_record[j].second;
        break;
      }
    }
    if (nf_type == NearFarType::FAR)
      return nf_type;
    for (size_t j = i; j < near_far_record.size(); j++) {
      if (near_far_record[j].first) {
        nf_type = near_far_record[j].second;
        break;
      }
    }
    return nf_type;
  };
  Ks.reserve(Ks_depth.size());
  for (size_t i = 0; i < Ks_depth.size(); i++) {
    Ks.emplace_back(Ks_depth[i].K);
    NearFarType nf_type = NearFarType::NEAR;
    if (!Ks_depth[i].is_consistent) {
      nf_type = look_up_near_far(i);
      if (nf_type == NearFarType::NEAR)
        Ks[i].ori *= -1;
    }
  }
}

real_t compute_2d_rotation_indices(Mesh const &mesh, Camera &camera,
                                   Vertex const &v, FacingType facing) {
  std::vector<PatchRotationAngle> Ks;

  project_one_ring(mesh, camera, v, facing, Ks);

  // Find the correct starting point
  PatchRotationAngle acc;
  for (size_t i = 0; i < Ks.size(); i++) {
    acc = add_triangle(acc, Ks[i]);
  }

  real_t r = acc.r;

  return r;
}

void compute_2d_rotation_indices(Mesh &mesh, Camera &camera) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto rot_idx = mesh.halfedge_property<real_t>("h:rot_idx");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");

  if (!rot_idx) {
    rot_idx = mesh.add_halfedge_property<real_t>("h:rot_idx");
    ;
  } else {
    rot_idx.vector().assign(rot_idx.vector().size(), 0);
  }

  // Get p info (we have one rotation index per p)
  for (uint32_t i = 0; i < mesh.n_vertices(); ++i) {
    Vertex v = Vertex(i);
    auto hit = mesh.halfedges(v);
    auto hit_end = hit;
    do {
      Edge e = mesh.edge(*hit);
      // Only compute the rotation indices at contour vertices
      // Each contour vertex is corresponding to two rotation indices
      // indicated by two half edges.
      if (is_contour[e] >= 0) {
        FacingType facing = VBO[mesh.face(*hit)];
        rot_idx[*hit] = compute_2d_rotation_indices(mesh, camera, v, facing);
      }
    } while (++hit != hit_end);
  }
}
