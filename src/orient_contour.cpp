// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "orient_contour.h"

#include <igl/predicates/predicates.h>
#include <limits>
#include <unordered_set>

#include "common.h"
#include "logger.h"
#include "ray_cast_QI.h"

int orient_contour_intrinsic(Mesh &mesh, Camera const &camera,
                             Halfedge const &h) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto VBO_f = mesh.get_face_property<FacingType>("f:VBO_f");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");

  int int_orient = 0;

  // Determine the orientation of the 2D projection of the adjacent face
  Face f = mesh.face(h);
  Vertex vv[3];
  mesh.verticesOfFace(f, vv);
  Vector2f v0_2d = project(vpositions[vv[0]], camera.viewMatrix().matrix(),
                           camera.projectionMatrix(),
                           Vector2i(camera.vpWidth(), camera.vpHeight()))
                       .head(2);
  Vector2f v1_2d = project(vpositions[vv[1]], camera.viewMatrix().matrix(),
                           camera.projectionMatrix(),
                           Vector2i(camera.vpWidth(), camera.vpHeight()))
                       .head(2);
  Vector2f v2_2d = project(vpositions[vv[2]], camera.viewMatrix().matrix(),
                           camera.projectionMatrix(),
                           Vector2i(camera.vpWidth(), camera.vpHeight()))
                       .head(2);

  auto f_ori = igl::predicates::orient2d(v0_2d, v1_2d, v2_2d);
  contess_assert_msg(f_ori != igl::predicates::Orientation::COLLINEAR,
                     "Generic camera assumption violated.");
  // 1:
  int_orient = (f_ori == igl::predicates::Orientation::POSITIVE) ? 1 : -1;

  // Flip the orientation
  if (VBO[f] != VBO_f[f])
    int_orient *= -1;

  return int_orient;
}

int orient_contour_extrinsic(Mesh &mesh, Camera const &camera,
                             Halfedge const &h) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  int ext_orient = 0;

  // Determine which side the patch mesh is on: left or right?
  // By walking on the 3D mesh
  Face f = mesh.face(h);
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");

  Edge end_e = mesh.edge(h);

  // We pick the cut plane that is orthogonal to the contour edge
  Vector3f contour_mid_point =
      (vpositions[mesh.vertex(end_e, 0)] + vpositions[mesh.vertex(end_e, 1)]) *
      0.5;
  Vector3f dir = contour_mid_point - camera.position();
  Vector3f cut_p = camera.position();
  Vector3f edge_dir =
      vpositions[mesh.vertex(end_e, 0)] - vpositions[mesh.vertex(end_e, 1)];
  Vector3f cut_norm = (dir.cross(edge_dir)).cross(dir);
  cut_norm.normalized();

  std::unordered_set<int> visited_face_indices;
  std::queue<Face> face_frontier;
  face_frontier.push(f);
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
    // (cur_f == f is used to avoid extremly small contour edges)
    if (planeIntersectTri(cut_p, cut_norm, tri_p) || cur_f == f) {
      // Extend frontier
      auto fhit = mesh.halfedges(cur_f), fhit_end = fhit;
      do {
        Edge fe = mesh.edge(*fhit);

        // The edge must intersect with the cut plane
        // so we don't dodge the small hole or patch corner
        std::vector<Vector3f> line_p(
            {vpositions[mesh.vertex(fe, 0)], vpositions[mesh.vertex(fe, 1)]});
        if (!planeIntersectTri(cut_p, cut_norm, line_p))
          continue;

        // We can't cross the contour edge
        if (is_contour[fe] >= 0) {
          // Keeps updating the end edge
          // The final results should be the furthest edge
          end_e = fe;
          continue;
        }

        // Check adjacent faces
        auto adj_f_h = mesh.opposite_halfedge(*fhit);
        Face adj_f = mesh.face(adj_f_h);
        if (visited_face_indices.find(adj_f.idx()) !=
            visited_face_indices.end())
          continue;

        // In case we are at the boundary of the model
        if (adj_f.is_valid())
          face_frontier.push(adj_f);
        else {
          end_e = mesh.edge(*fhit);
        }
      } while (++fhit != fhit_end);
    }

    visited_face_indices.insert(cur_f.idx());
  } while (!face_frontier.empty());

  contess_assert_msg(end_e != mesh.edge(h), "Can't walk to the end.");

  // Determine the side of the furthest walkable position
  Vector2f v0_2d =
      project(vpositions[mesh.from_vertex(h)], camera.viewMatrix().matrix(),
              camera.projectionMatrix(),
              Vector2i(camera.vpWidth(), camera.vpHeight()))
          .head(2);
  Vector2f v1_2d =
      project(vpositions[mesh.to_vertex(h)], camera.viewMatrix().matrix(),
              camera.projectionMatrix(),
              Vector2i(camera.vpWidth(), camera.vpHeight()))
          .head(2);
  Vector2f v2_2d =
      project(0.5 * (vpositions[mesh.vertex(end_e, 0)] +
                     vpositions[mesh.vertex(end_e, 1)]),
              camera.viewMatrix().matrix(), camera.projectionMatrix(),
              Vector2i(camera.vpWidth(), camera.vpHeight()))
          .head(2);

  auto f_ori = igl::predicates::orient2d(v0_2d, v1_2d, v2_2d);
  contess_assert_msg(f_ori != igl::predicates::Orientation::COLLINEAR,
                     "Generic camera assumption violated.");
  // 1:
  ext_orient = (f_ori == igl::predicates::Orientation::POSITIVE) ? 1 : -1;
  return ext_orient;
}

void orient_contour(Mesh &mesh, Camera const &camera) {
  auto int_orient = mesh.halfedge_property<int>("h:int_orient");
  if (!int_orient) {
    int_orient = mesh.add_halfedge_property<int>("h:int_orient");
  }
  // 0 means NA, the valid values are 1, -1
  int_orient.vector().assign(int_orient.vector().size(), 0);

  auto ext_orient = mesh.halfedge_property<int>("h:ext_orient");
  if (!ext_orient) {
    ext_orient = mesh.add_halfedge_property<int>("h:ext_orient");
  }
  // 0 means NA, the valid values are 1, -1
  ext_orient.vector().assign(ext_orient.vector().size(), 0);

  auto is_contour = mesh.get_edge_property<real_t>("e:contour");

  for (uint32_t i = 0; i != mesh.n_halfedges(); ++i) {
    Halfedge h = Halfedge(i);
    Edge e = mesh.edge(h);
    if (is_contour[e] < 0)
      continue;

    int_orient[h] = orient_contour_intrinsic(mesh, camera, h);
    ext_orient[h] = orient_contour_extrinsic(mesh, camera, h);
  }
}
