// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "insert_interpolated_contours.h"

#include "common.h"
#include "logger.h"
#include <spdlog/fmt/ostr.h>
#include <unordered_set>

// Use the face normal to determine the facing type
bool determine_interpolated_facing(Mesh &mesh, Vector3f const &view_pos,
                                   Face const &f, FacingType &facing) {
  bool non_interpolated = true;
  surface_mesh::Surface_mesh::Halfedge hs[3];
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");

  mesh.halfedgesOfFace(f, hs);
  for (auto const &h : hs) {
    Edge e = mesh.edge(h);
    if (is_contour[e] >= 0) {
      non_interpolated = false;
      facing = FacingType::UNDEFINED;
      return non_interpolated;
    }
  }

  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");

  // For stability, only confirm the facing when the no-generic camera
  // assumption holds
  Vector3f normal = mesh.compute_face_normal(f);
  auto vit = mesh.vertices(f);
  Vertex v = *vit;
  Vector3f const &v_pos = vpositions[v];

  real_t ndotv = normal.dot((view_pos - v_pos).normalized());
  if (std::abs(ndotv) < CONTOUR_THRESHOLD) {
    non_interpolated = false;
    facing = FacingType::UNDEFINED;
    return non_interpolated;
  }

  if (ndotv > 0)
    facing = FacingType::FRONT;
  else
    facing = FacingType::BACK;

  return non_interpolated;
}

// Use the g (ndotv) function at vertices to determine the facing type
bool determine_interpolated_facing(Mesh &mesh, Face const &f,
                                   FacingType &facing) {
  auto ndotv = mesh.get_vertex_property<real_t>("v:ndotv");

  // Find the most stable ndotv at vertices (the one with the largest abs value)
  double max_abs_ndotv = 0;

  auto vit = mesh.vertices(f);
  auto vit_end = vit;
  do {
    Vertex v = *vit;
    double ndv = ndotv[v];

    if (std::abs(ndv) > std::abs(max_abs_ndotv))
      max_abs_ndotv = ndv;
  } while (++vit != vit_end);

  if (std::abs(max_abs_ndotv) < CONTOUR_THRESHOLD) {
    facing = FacingType::UNDEFINED;
    return false;
  }

  if (max_abs_ndotv > 0)
    facing = FacingType::FRONT;
  else
    facing = FacingType::BACK;

  return true;
}

bool consistently_label_interpolated_contours(Mesh &mesh,
                                              Vector3f const &view_pos,
                                              bool to_update_VBO) {
  // Store the original face-normal-based facing type
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto VBO_f = mesh.face_property<FacingType>("f:VBO_f");
  if (!VBO_f) {
    VBO_f = mesh.add_face_property<FacingType>("f:VBO_f", FacingType::NA);
  }
  VBO_f.vector().assign(VBO_f.vector().size(), FacingType::NA);

  for (Face f : mesh.faces()) {
    // The facing type based on the face normal
    auto f_norm = mesh.compute_face_normal(f);
    Vertex vs[3];
    mesh.verticesOfFace(f, vs);
    auto c_ray = vpositions[vs[0]] - view_pos;
    c_ray.normalized();
    VBO_f[f] = (f_norm.dot(c_ray) < 0) ? FacingType::FRONT : FacingType::BACK;
    // Use predicate to be mroe accurate
    auto ori = igl::predicates::orient3d(vpositions[vs[0]], vpositions[vs[1]],
                                         vpositions[vs[2]], view_pos);
    VBO_f[f] = (ori == igl::predicates::Orientation::NEGATIVE)
                   ? FacingType::FRONT
                   : FacingType::BACK;
  }

  if (!to_update_VBO)
    return true;

  // Update with the vertex-normal based facing type
  auto VBO = mesh.face_property<FacingType>("f:VBO");
  if (!VBO) {
    VBO = mesh.add_face_property<FacingType>("f:VBO");
  }
  VBO.vector().assign(VBO.vector().size(), FacingType::NA);

  auto patchID = mesh.get_face_property<int>("f:patchID");

  for (Face f : mesh.faces()) {
    if (VBO[f] != FacingType::NA) // already visited
      continue;

    FacingType facing;

    // Determine using the g function at vertices
    bool is_stable = determine_interpolated_facing(mesh, f, facing);

    // Take advantage of patch VBO if exists
    if (!is_stable && patchID) {
      facing = mesh.get_patch_facing(patchID[f]);
      is_stable = true;
    }
    if (!is_stable) {
      logger().warn("ndotv at vertices is not numerically stable at {}.", f);
      // return false;
      facing = FacingType::FRONT;
    }

    VBO[f] = facing;
  }

  return true;
}

void insert_interpolated_contours(Mesh &mesh, Vector3f const &view_pos) {
  mesh.m_contour_edges.clear();

  // 1. Create the map between two vertices adjacent to each other on the
  // contour
  // http://www.alecjacobson.com/weblog/?p=4294
  // https://stackoverflow.com/questions/27048146/using-stdmap-with-eigen-3
  auto cmp = [](VectorXf const &a, VectorXf const &b) -> bool {
    return std::lexicographical_compare(a.data(), a.data() + a.size(), b.data(),
                                        b.data() + b.size());
  };
  std::map<VectorXf, std::vector<VectorXf>, decltype(cmp)>
      interpolated_edges_map(cmp);
  for (size_t i = 0; i < mesh.get_interpolated_contour_faces().size(); i++) {
    Vector3f v1 = mesh.get_interpolated_contours().col(2 * i);
    Vector3f v2 = mesh.get_interpolated_contours().col(2 * i + 1);

    auto insert_end_v = [&interpolated_edges_map](Vector3f const &v,
                                                  VectorXf const &end_v) {
      VectorXf vx = v;
      if (interpolated_edges_map.find(vx) == interpolated_edges_map.end()) {
        interpolated_edges_map[vx] = std::vector<VectorXf>();
      }
      interpolated_edges_map[vx].emplace_back(end_v);
    };

    insert_end_v(v1, v2);
    insert_end_v(v2, v1);
  }

  // 2. Iterate on edges and determine the contour vertex positions.
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  contess_assert_msg(is_contour,
                     "Missing interpolated contours. Cannot insert.");

  std::vector<std::pair<Vector3f, Edge>> interpolated_edges;
  for (auto eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit) {
    if (is_contour[*eit] < 0)
      continue;

    Vertex i0 = mesh.vertex(*eit, 0);
    Vertex i1 = mesh.vertex(*eit, 1);
    const Vector3f &v0 = vpositions[i0];
    const Vector3f &v1 = vpositions[i1];
    real_t w = is_contour[*eit];

    Vector3f v = w * v1 + (1.f - w) * v0;
    interpolated_edges.emplace_back(v, (*eit).idx());
  }

  // 3. Insert all the interpolated contour vertices
  std::unordered_set<int> contour_vertices;
  contour_vertices.reserve(interpolated_edges.size());
  for (size_t i = 0; i < interpolated_edges.size(); i++) {
    auto v_e = interpolated_edges[i];
    Edge e = v_e.second;

    // Split the edge
    auto v = mesh.split(e, v_e.first);
    contour_vertices.insert(v.idx());
  }

  // 4. Label the contour edges
  is_contour.vector().assign(is_contour.vector().size(), -1.f);
  for (int v_idx : contour_vertices) {
    Vertex v(v_idx);

    surface_mesh::Surface_mesh::Halfedge_around_vertex_circulator chit,
        chit_end;
    chit = chit_end = mesh.halfedges(v); // outgoing half edge here
    do {
      Vertex to_v = mesh.to_vertex(*chit);
      // Found a contour edge
      VectorXf v1_pos = vpositions[v];
      VectorXf v2_pos = vpositions[to_v];

      contess_assert_msg(
          interpolated_edges_map[v1_pos].size() == 2 ||
              interpolated_edges_map[v1_pos].size() == 1,
          "Contour vertex is either an endpoint or on a 1-manifold.");

      if (contour_vertices.find(to_v.idx()) != contour_vertices.end() &&
          ((interpolated_edges_map[v1_pos].front() - v2_pos).norm() <
               std::numeric_limits<real_t>::epsilon() ||
           (interpolated_edges_map[v1_pos].back() - v2_pos).norm() <
               std::numeric_limits<real_t>::epsilon())) {
        Edge e = mesh.edge(*chit);
        is_contour[e] = 1.f;
        mesh.m_contour_edges.emplace_back(e);
      }
    } while (++chit != chit_end);
  }

  // 5. Clean up.
  mesh.garbage_collection();

  // 6. Use ndotv at vertices to consistently label
  // (this gives patches with consistent facing types given the mesh
  // has interpolated contours inserted)
  consistently_label_interpolated_contours(mesh, view_pos);
}
