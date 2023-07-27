// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "lift_min_surface.h"
#include "common.h"

#include <Eigen/src/Core/Matrix.h>
#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/slice.h>
#include <igl/slice_into.h>

void lift_min_surface(Mesh const &mesh, Mesh const &flat_patch,
                      Mesh &min_surface_patch) {
  auto orig_idx = flat_patch.get_vertex_property<int>("v:orig_idx");
  contess_assert_msg(orig_idx,
                     "lift_min_surface: Missing index correspondence.");

  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto flat_vpositions = flat_patch.get_vertex_property<Vector3f>("v:point");

  // Init the result surface
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  V.resize(flat_patch.n_vertices(), 3);
  F.resize(flat_patch.n_faces(), 3);
  for (long i = 0; i < (int)flat_patch.n_vertices(); i++) {
    V.row(i) = flat_vpositions[Vertex(i)].transpose();
  }
  for (long i = 0; i < (int)flat_patch.n_faces(); i++) {
    Vertex vv[3];
    flat_patch.verticesOfFace(Face(i), vv);

    F.row(i) = Eigen::RowVector3i(vv[0].idx(), vv[1].idx(), vv[2].idx());
  }
  min_surface_patch = flat_patch;
  min_surface_patch.m_fedges.clear();
  min_surface_patch.m_svertices.clear();

  auto minf_orig_idx = min_surface_patch.get_vertex_property<int>("v:orig_idx");
  auto minf_vpositions = min_surface_patch.vertex_property<Vector3f>("v:point");

  // Move the boundary / cut points to the positions from original mesh
  std::vector<int32_t> fixed_idx;
  for (int i = 0; i < (int)min_surface_patch.n_vertices(); i++) {
    Vertex v(i);
    int v_idx_orig = minf_orig_idx[v];
    if (v_idx_orig < 0)
      continue;

    fixed_idx.emplace_back(i);
    minf_vpositions[v] = vpositions[Vertex(v_idx_orig)];
    V.row(i) = vpositions[Vertex(v_idx_orig)].transpose();
  }
}
