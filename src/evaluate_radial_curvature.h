// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

void evaluate_radial_curvature_analytic(Mesh &mesh, Subdiv const &subdiv,
                                        Camera const &camera);
real_t evaluate_radial_curvature_analytic(Mesh &mesh, Subdiv const &subdiv,
                                          Camera const &camera, Vertex const &v,
                                          Param_loc const &param);

void verify_first_partial_derivatives(Mesh &mesh, Subdiv const &subdiv,
                                      real_t step_size = 1e-7);
void verify_second_partial_derivatives(Mesh &mesh, Subdiv const &subdiv,
                                       real_t step_size = 1e-7);

void approach_extraordinary_points(Mesh &mesh, Subdiv const &subdiv,
                                   Camera const &camera);
