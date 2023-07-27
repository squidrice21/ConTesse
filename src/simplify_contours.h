// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

void simplify_contours(Mesh &mesh, Camera const &camera,
                       real_t simplify_epsilon);
