// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"
#include <vector>

void chain_rendered_contour(Camera const &camera, Mesh &stitched_mesh);

void deduplicate_contour_faces(Mesh &mesh, Mesh const &orig_mesh,
                               Camera const &camera);

void stitch_mesh(Mesh const &mesh, Camera const &camera,
                 std::vector<Mesh *> const &patches, Mesh &stitched_mesh);
