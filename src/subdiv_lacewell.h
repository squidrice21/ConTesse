// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

// https://www.tandfonline.com/doi/abs/10.1080/2151237X.2007.10129243
// Lacewell and Burley

#pragma once

#include "subdiv_common.h"
#include "surface_mesh/surface_mesh.h"

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class Mesh;
class Subd; // form external library.

class SubdivLacewell {
public:
  SubdivLacewell();
  ~SubdivLacewell();

  void load(const std::string &filename, int num_subdiv = 0);
  // after we load the mesh we subdivide it with
  // osubdiv (the Lacewell code does not handle
  // mixed stuff). Then we pass that shape to create
  // this subdiv object.
  void load(Shape, int num_subdiv = 0);

  // Populate a mesh based this subdiv object.
  void create_mesh(Mesh &);

  // Evaluate the values wrt limit surface given the parameter location
  bool evaluateLimit(const Param_loc &p, Vector3f &position,
                     Vector3f &normal) const;
  bool evaluateLimitFrame(const Param_loc &p, Vector3f &position, Vector3f &ds,
                          Vector3f &dt, Vector3f *dsds = nullptr,
                          Vector3f *dsdt = nullptr,
                          Vector3f *dtdt = nullptr) const;

  void determine_extraordinary_vertices(Mesh const &mesh);
  bool is_near_extraordinary(const Param_loc &p) const;
  void round_extraordinary(Param_loc &p) const;

  std::unique_ptr<Subd> m_impl;
  Shape m_shape;
  Mesh *m_mesh;
  std::unordered_map<int, std::vector<Param_loc>> extraordinary_parameters;
  std::unordered_map<int, std::vector<Param_loc>>
      exact_extraordinary_parameters;
  std::unordered_map<int, int> ep_vertices;
  real_t step_size = 1e-7;
};
