// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

void lift_min_surface(Mesh const &mesh, Mesh const &flat_patch,
                      Mesh &min_surface_patch);
