// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "fix_flipped_faces.h"
#include "insert_interpolated_contours.h"
#include "subdivide_contour_edges.h"

#include <spdlog/fmt/ostr.h>

bool can_flip_edge(Mesh const &mesh, Edge const &flip_e) {
  auto facing = mesh.get_vertex_property<FacingType>("v:facing");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");

  // 0. Check if flipping would generate a face containing three contour
  // vertices
  Halfedge he1 = mesh.halfedge(flip_e, 0);
  Halfedge he2 = mesh.halfedge(flip_e, 1);
  auto is_all_contour = [&](Halfedge const &h1, Halfedge const &h2) -> bool {
    return (facing[mesh.from_vertex(h1)] == FacingType::CONTOUR) &&
           (facing[mesh.to_vertex(mesh.next_halfedge(h1))] ==
            FacingType::CONTOUR) &&
           (facing[mesh.to_vertex(mesh.next_halfedge(h2))] ==
            FacingType::CONTOUR);
  };
  if (mesh.face(he1).is_valid() && mesh.face(he2).is_valid() &&
      (is_all_contour(he1, he2) || is_all_contour(he2, he1)))
    return false;

  // 1. Check parameter quad
  Param_loc before_flip = param_loc[mesh.halfedge(flip_e, 0)];

  // Flipping check if we are flipping param quad boundary
  bool can_flip =
      (param_loc[mesh.next_halfedge(
                     mesh.next_halfedge(mesh.halfedge(flip_e, 0)))]
           .ptexIndex == param_loc[mesh.next_halfedge(mesh.next_halfedge(
                                       mesh.halfedge(flip_e, 1)))]
                             .ptexIndex);
  if (!can_flip)
    return can_flip;

  // 2. Check if surface is flipped
  Vertex vv1 = mesh.to_vertex(mesh.next_halfedge(mesh.halfedge(flip_e, 0)));
  Vertex vv11 = mesh.to_vertex(mesh.halfedge(flip_e, 0));
  Vertex vv12 = mesh.from_vertex(mesh.halfedge(flip_e, 0));

  Vertex vv2 = mesh.to_vertex(mesh.next_halfedge(mesh.halfedge(flip_e, 1)));

  // Reference point
  Vector3f f_norm =
      mesh.compute_face_normal(mesh.face(mesh.halfedge(flip_e, 0)));
  Vector3f v_ref = vpositions[vv1] + f_norm;

  auto ori1 = igl::predicates::orient3d(vpositions[vv1], vpositions[vv2],
                                        vpositions[vv11], v_ref);
  auto ori2 = igl::predicates::orient3d(vpositions[vv2], vpositions[vv1],
                                        vpositions[vv12], v_ref);

  return ori1 == ori2;
}

bool is_param_boundary(Mesh const &mesh, Edge const &e) {
  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");

  // 1. Check parameter quad
  Param_loc before_flip = param_loc[mesh.halfedge(e, 0)];

  // Flipping check if we are flipping param quad boundary
  bool can_flip =
      (param_loc[mesh.next_halfedge(mesh.next_halfedge(mesh.halfedge(e, 0)))]
           .ptexIndex ==
       param_loc[mesh.next_halfedge(mesh.next_halfedge(mesh.halfedge(e, 1)))]
           .ptexIndex);
  return !can_flip;
}

bool is_adj_faces_flipped(Mesh const &mesh, Edge const &e) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");

  // 2. Check if surface is flipped
  Vertex vv1 = mesh.to_vertex(mesh.next_halfedge(mesh.halfedge(e, 0)));
  Vertex vv11 = mesh.to_vertex(mesh.halfedge(e, 0));
  Vertex vv12 = mesh.from_vertex(mesh.halfedge(e, 0));

  Vertex vv2 = mesh.to_vertex(mesh.next_halfedge(mesh.halfedge(e, 1)));
  Vertex vv21 = mesh.to_vertex(mesh.halfedge(e, 1));
  Vertex vv22 = mesh.from_vertex(mesh.halfedge(e, 1));

  // Reference point
  Vector3f f_norm = mesh.compute_face_normal(mesh.face(mesh.halfedge(e, 0)));
  Vector3f v_ref = vpositions[vv11] + 10 * f_norm;

  auto ori1 = igl::predicates::orient3d(vpositions[vv11], vpositions[vv1],
                                        vpositions[vv12], v_ref);
  auto ori2 = igl::predicates::orient3d(vpositions[vv11], vpositions[vv12],
                                        vpositions[vv2], v_ref);

  return !(ori1 == ori2);
}

void fix_flipped_faces(Mesh &mesh, Camera const &camera) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto param_loc = mesh.halfedge_property<Param_loc>("h:param_loc");

  auto contour_adj = [&](Vertex const &v) -> bool {
    auto hit = mesh.halfedges(v);
    auto hit_end = hit;
    do {
      if (is_contour[mesh.edge(*hit)] >= 0)
        return true;
    } while (++hit != hit_end);
    return false;
  };

  for (size_t i = 0; i < mesh.n_edges(); i++) {
    Edge e(i);
    if (is_contour[e] >= 0)
      continue;

    // Adjacent to a contour face
    if (!(is_contour[mesh.edge(mesh.next_halfedge(mesh.halfedge(e, 0)))] >= 0 ||
          is_contour[mesh.edge(mesh.next_halfedge(mesh.halfedge(e, 1)))] >= 0 ||
          contour_adj(
              mesh.to_vertex(mesh.next_halfedge(mesh.halfedge(e, 0)))) ||
          contour_adj(mesh.to_vertex(mesh.next_halfedge(mesh.halfedge(e, 1))))))
      continue;

    if (!mesh.face(e, 0).is_valid() || !mesh.face(e, 1).is_valid() ||
        !is_adj_faces_flipped(mesh, e)) {
      if (!mesh.face(e, 0).is_valid() || !mesh.face(e, 1).is_valid())
        continue;
    }
    // We can't flip the parameter boundary
    if (is_param_boundary(mesh, e) || !can_flip_edge(mesh, e))
      continue;

    // Flipping
    auto update_param = [&](Halfedge const &he, int ptexIndex) {
      Param_loc existing_param;
      auto hit = mesh.halfedges(mesh.from_vertex(he));
      auto hit_end = hit;
      do {
        if (*hit == he)
          continue;
        if (param_loc[*hit].ptexIndex == ptexIndex) {
          existing_param = param_loc[*hit];
          break;
        }
      } while (++hit != hit_end);
      if (!existing_param.is_valid()) {
        hit = mesh.halfedges(mesh.from_vertex(he));
        hit_end = hit;
        do {
          if (*hit == he)
            continue;
          if (param_loc[*hit].is_valid()) {
            existing_param = param_loc[*hit];
            break;
          }
        } while (++hit != hit_end);
      }
      param_loc[he] = existing_param;
    };
    Edge flip_e = e;
    Param_loc before_flip = param_loc[mesh.halfedge(flip_e, 0)];
    Vertex vv1 = mesh.to_vertex(mesh.next_halfedge(mesh.halfedge(flip_e, 0)));
    Vertex vv2 = mesh.to_vertex(mesh.next_halfedge(mesh.halfedge(flip_e, 1)));

    if (!mesh.is_flip_ok(flip_e))
      continue;

    mesh.flip(flip_e);

    // Update parameter
    update_param(mesh.find_halfedge(vv1, vv2), before_flip.ptexIndex);
    update_param(mesh.find_halfedge(vv2, vv1), before_flip.ptexIndex);
  }

  auto patchID = mesh.get_face_property<int>("f:patchID");
  // Update the patch
  if (patchID)
    update_patch(mesh, camera);

  // Rebuild the face labeling
  consistently_label_interpolated_contours(mesh, camera.position(), (!patchID));

  // Check parameterization is correct
  for (auto hit = mesh.halfedges_begin(); hit != mesh.halfedges_end(); ++hit)
    if (!mesh.is_boundary(*hit))
      assert(param_loc[*hit].ptexIndex != -1);
}
