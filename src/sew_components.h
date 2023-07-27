#pragma once

#include "mesh.h"

void sew_components(Mesh &mesh, Mesh const &orig_mesh, Camera const &camera,
                    std::vector<std::pair<int, int>> failed_paths);
