// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include <cassert>

#include "save_mesh_vtk.h"
#include "common.h"
#include "logger.h"
#include "surface_mesh/surface_mesh.h"
#include "io/visit_writer.h"

void save_mesh_vtk(const std::string fname, const surface_mesh::Surface_mesh& mesh) {
  using surface_mesh::Surface_mesh;
  using surface_mesh::Point;

  logger().info("Writing vtk mesh to {}", fname);

  Visit_args args;
  args.useBinary = 1;
  args.filename = fname.c_str();

  //vertices
  std::vector<float> points_vector;
  {
    Surface_mesh::Vertex_property<Point> points = mesh.get_vertex_property<Point>("v:point");
    for (Surface_mesh::Vertex_iterator vit=mesh.vertices_begin(); vit!=mesh.vertices_end(); ++vit) {
        const Point& p = points[*vit];
        points_vector.push_back(static_cast<float>(p[0]));
        points_vector.push_back(static_cast<float>(p[1]));
        points_vector.push_back(static_cast<float>(p[2]));
    }
  }
  args.pts = points_vector.data();
  args.npts = static_cast<int>(points_vector.size()/3);

  //faces
  std::vector<int> celltypes;
  std::vector<int> conn;
  for (Surface_mesh::Face_iterator fit=mesh.faces_begin(); fit!=mesh.faces_end(); ++fit) {
    Surface_mesh::Vertex_around_face_circulator fvit=mesh.vertices(*fit), fvend=fvit;
    int num_verts = 0;
    do {
      conn.push_back((*fvit).idx());
      ++num_verts;
    } while (++fvit != fvend);
    if(num_verts==3){
      celltypes.push_back(VISIT_TRIANGLE);
    } else if (num_verts == 4) {
      conn.push_back(VISIT_QUAD);
    } else {
      contess_assert_msg(false, "Face with " << num_verts << " vertices is not supported in vtk");
    }
  }
  args.conn = conn.data();
  args.celltypes = celltypes.data();
  args.ncells = static_cast<int>(celltypes.size());

  args.write();
}

