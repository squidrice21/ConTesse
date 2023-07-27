// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include "mesh.h"

bool consistently_label_interpolated_contours(Mesh &mesh,
                                              Vector3f const &view_pos,
                                              bool to_update_VBO = true);

void insert_interpolated_contours(Mesh &mesh, Vector3f const &view_pos);
