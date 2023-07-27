// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

bool chain_contour(Mesh &mesh, Camera const &camera, bool to_cache_chain,
                   bool ff_only = false);

bool chain_cut_graph(Mesh &mesh, Camera const &camera, bool ff_only = false);
