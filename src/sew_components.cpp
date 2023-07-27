#include "sew_components.h"
#include "common.h"
#include "inflat_path.h"
#include <igl/predicates/predicates.h>

#include <spdlog/fmt/ostr.h>
#include <unordered_set>
#include <utility>
#include <vector>

bool is_delaunay(Mesh const &mesh, Mesh const &orig_mesh, Camera const &camera,
                 Edge const &e) {
  if (!mesh.is_flip_ok(e))
    return true;

  auto f_sew_component = mesh.get_face_property<int>("f:sew_component");
  Face f1 = mesh.face(e, 0);
  Face f2 = mesh.face(e, 1);

  if (!f1.is_valid() || !f2.is_valid())
    return true;

  if (f_sew_component[f1] < 0 || f_sew_component[f2] < 0)
    return true;

  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  Halfedge h1 = mesh.halfedge(e, 0);
  Halfedge h2 = mesh.halfedge(e, 1);
  std::vector<Vertex> incircle_test_vertices(
      {mesh.from_vertex(h1), mesh.to_vertex(h1),
       mesh.to_vertex(mesh.next_halfedge(h1)),
       mesh.to_vertex(mesh.next_halfedge(h2))});

  // Flip the testing triangle if the patch is back facing
  auto patchID = mesh.get_face_property<int>("f:patchID");
  int pid = patchID[*mesh.faces_begin()];
  bool is_back_facing = orig_mesh.get_patch_facing(pid) == FacingType::BACK;
  if (is_back_facing) {
    incircle_test_vertices =
        std::vector<Vertex>({mesh.to_vertex(h1), mesh.from_vertex(h1),
                             mesh.to_vertex(mesh.next_halfedge(h1)),
                             mesh.to_vertex(mesh.next_halfedge(h2))});
  }

  Vector2i viewport({camera.vpWidth(), camera.vpHeight()});
  auto project2d = [&](Vertex const &v, Vector2f &pos2D) {
    pos2D = project(vpositions[v], camera.viewMatrix().matrix(),
                    camera.projectionMatrix(), viewport)
                .head<2>();
  };
  std::vector<Vector2f> proj_v;
  proj_v.resize(incircle_test_vertices.size());
  project2d(incircle_test_vertices[0], proj_v[0]);
  project2d(incircle_test_vertices[1], proj_v[1]);
  project2d(incircle_test_vertices[2], proj_v[2]);
  project2d(incircle_test_vertices[3], proj_v[3]);
  auto orientation =
      igl::predicates::incircle(proj_v[0], proj_v[1], proj_v[2], proj_v[3]);
  return orientation != igl::predicates::Orientation::INSIDE;
}

void sew_components(Mesh &mesh, Mesh const &orig_mesh, Camera const &camera,
                    std::vector<std::pair<int, int>> failed_paths) {
  auto orig_idx = mesh.get_vertex_property<int>("v:orig_idx");
  auto patch_component = mesh.get_face_property<int>("f:patch_component");

  auto f_sew_component = mesh.get_face_property<int>("f:sew_component");
  if (!f_sew_component) {
    mesh.add_face_property<int>("f:sew_component", -1);
  }
  f_sew_component = mesh.face_property<int>("f:sew_component");

  // Get orig->vertex
  std::unordered_map<int, int> orig2vidx;
  for (size_t i = 0; i < mesh.n_vertices(); i++) {
    Vertex v(i);
    if (orig_idx[v] >= 0) {
      orig2vidx[orig_idx[v]] = i;
    }
  }

  // 1. Find the ring to retriangulate
  size_t ring_size = 5;
  for (auto const &path : failed_paths) {
    contess_assert(orig2vidx.count(path.first) && orig2vidx.count(path.second));
    Vertex v1(orig2vidx[path.first]), v2(orig2vidx[path.second]);

    std::vector<Vertex> one_ring_vertices;
    std::vector<Face> one_ring_faces;
    std::unordered_set<int> one_ring_faces_set;
    get_neighbors(mesh, -1, v1, ring_size, one_ring_faces_set);
    get_neighbors(mesh, -1, v2, ring_size, one_ring_faces_set);

    // Label the retriangulation region
    for (auto const &ff : one_ring_faces_set) {
      Face f(ff);
      f_sew_component[f] = patch_component[f];
    }
  }

  // 2. Delaunay flip
  bool all_delaunay = false;
  while (!all_delaunay) {
    all_delaunay = true;
    for (size_t i = 0; i < mesh.n_edges(); i++) {
      Edge e(i);
      if (!is_delaunay(mesh, orig_mesh, camera, e)) {
        Face f1 = mesh.face(e, 0);
        int comp_before = f_sew_component[f1];
        all_delaunay = false;
        if (mesh.is_flip_ok(e)) {
          mesh.flip(e);
          Face f1 = mesh.face(e, 0);
          Face f2 = mesh.face(e, 1);
          f_sew_component[f1] = comp_before;
          f_sew_component[f2] = comp_before;
        }
      }
    }
  }
}
