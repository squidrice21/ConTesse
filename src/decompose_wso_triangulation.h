// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

void decompose_wso_triangulation(Mesh const &mesh, Camera const &camera,
                                 Mesh &flat_patch);
