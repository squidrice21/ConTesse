// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include <string>

namespace surface_mesh {
  class Surface_mesh;
}

/**
Save a mesh with the vtk format. This is not used at the moment.
*/
void save_mesh_vtk(const std::string fname, const surface_mesh::Surface_mesh& mesh);

