// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"
#include <vector>

int count_cut_degree(Mesh const &mesh, Vertex const &v, int patch_id);

void cut_patch_to_disk_gu02(Mesh const &mesh, int patchID,
                            std::vector<int> &cut_halfedge_indices);
void cut_patch_to_disk(Mesh const &mesh, Camera const &camera, int patchID,
                       std::vector<int> &cut_halfedge_indices);
void cut_patch_to_disk(Mesh &mesh, Camera const &camera, bool use_gu02 = false,
                       bool ff_only = false);
