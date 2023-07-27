// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "stitch_mesh.h"
#include "chain_contour.h"
#include "common.h"
#include "insert_interpolated_contours.h"
#include "tag_concave_edges.h"
#include <limits>
#include <random>
#include <spdlog/fmt/ostr.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

void chain_rendered_contour(Camera const &camera, Mesh &stitched_mesh) {
  auto edge_qi = stitched_mesh.get_edge_property<int>("e:qi");
  contess_assert_msg(edge_qi, "chain_rendered_contour: QI is missing.");

  // 1. Chain the contour in the stitched mesh
  chain_contour(stitched_mesh, camera, false);

  // 2. Fill contour vertex types
  // A dead-end vertex is adjacent to a single visible curve (e.g., a visible
  // curtain fold); A connector vertex is adjacent to two visible curves; A
  // junction vertex is adjacent to more than two vertices, i.e., bifurcations
  // and image-space intersection vertices.
  auto contour_type =
      stitched_mesh.vertex_property<ContourVertexType>("v:contour_type");
  if (!contour_type) {
    contour_type = stitched_mesh.add_vertex_property<ContourVertexType>(
        "v:contour_type", ContourVertexType::CONNECTOR);
  }
  contour_type.vector().assign(contour_type.vector().size(),
                               ContourVertexType::CONNECTOR);

  auto stitched_intersection_2d =
      stitched_mesh.get_vertex_property<int>("v:intersection_2d");

  // Label cusps
  auto stitched_is_cusp = stitched_mesh.get_vertex_property<bool>("v:cusp");
  if (!stitched_is_cusp) {
    stitched_is_cusp = stitched_mesh.add_vertex_property<bool>("v:cusp", false);
  }
  stitched_is_cusp.vector().assign(stitched_is_cusp.vector().size(), false);
  auto is_contour = stitched_mesh.get_edge_property<real_t>("e:contour");

  tag_cusp_sharp(stitched_mesh, camera);
  for (size_t i = 0; i < stitched_mesh.n_vertices(); i++) {
    Vertex v(i);
    // Assign manifold boundary T-junction as cusps
    // This is for the automatic QI verification
    auto hit = stitched_mesh.halfedges(v), hit_end = hit;
    bool seen_manifold_boundary = false, seen_contour = false;
    do {
      if (!(*hit).is_valid())
        break;
      if (stitched_mesh.is_boundary(stitched_mesh.edge(*hit))) {
        seen_manifold_boundary = true;
      }
      if (is_contour[stitched_mesh.edge(*hit)] >= 0) {
        seen_contour = true;
      }
    } while (++hit != hit_end);
    if (seen_manifold_boundary && seen_contour)
      stitched_is_cusp[v] = true;
  }

  std::unordered_set<int> visited_patch;
  for (auto const &patch : stitched_mesh.get_patch_chains()) {
    if (visited_patch.count(patch.first))
      continue;

    visited_patch.insert(patch.first);

    // Set boundary cut directions
    auto chains = stitched_mesh.get_patch_chains().equal_range(patch.first);
    for (auto chain = chains.first; chain != chains.second; ++chain) {
      std::vector<Halfedge> chain_hes = *chain->second;
      for (size_t i = 0; i < chain_hes.size(); i++) {
        size_t j = (i + 1) % chain_hes.size();
        Halfedge h1 = chain_hes[i];
        Halfedge h2 = chain_hes[j];
        Vertex v = stitched_mesh.to_vertex(h1);

        // Determine the dead-end by checking QI
        if (edge_qi[stitched_mesh.edge(h1)] !=
            edge_qi[stitched_mesh.edge(h2)]) {
          contour_type[v] = ContourVertexType::DEAD_END;

          // Determine the junction by checking 2D intersections
          if (stitched_intersection_2d[v] >= 0) {
            contour_type[v] = ContourVertexType::JUNCTION;
            Vertex inter_v(stitched_intersection_2d[v]);
            contour_type[inter_v] = ContourVertexType::JUNCTION;
          }
        }
      }
    }
  }
}

void deduplicate_contour_faces(Mesh &mesh, Mesh const &orig_mesh,
                               Camera const &camera) {
  if (mesh.n_faces() == 0)
    return;

  auto patchID = mesh.get_face_property<int>("f:patchID");
  int patch_id = patchID[*mesh.faces_begin()];

  // To avoid inserting duplicated flipped faces
  // This may happend when there are two faces with all vertices on contour
  // from two adjacent patches (FF and BF) surface_mesh doesn't allow
  // inserting the flipped version of an existing face
  auto is_boundary_vertex = [&](Vertex v) -> bool {
    auto hit = mesh.halfedges(v), hit_end = hit;
    do {
      if (mesh.is_boundary(mesh.edge(*hit))) {
        return true;
      }
    } while (++hit != hit_end);
    return false;
  };

  std::unordered_set<int> edge_to_split;
  for (size_t i = 0; i < mesh.n_edges(); i++) {
    Edge e(i);
    if (mesh.is_boundary(e))
      continue;
    if (is_boundary_vertex(mesh.vertex(e, 0)) &&
        is_boundary_vertex(mesh.vertex(e, 1))) {
      edge_to_split.emplace(e.idx());
    }
  }

  // Convert to vertices
  std::vector<std::pair<Vertex, Vertex>> edge_to_split_v;
  for (auto e_i : edge_to_split) {
    Edge e(e_i);
    edge_to_split_v.emplace_back(
        std::make_pair(mesh.vertex(e, 0), mesh.vertex(e, 1)));
  }
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");

  // Use the patch index as seed so we have a different random sequence every
  // patch
  int seed = patch_id;
  std::mt19937 eng(seed);
  // Use the range [0.5, 1] to ensure always have some offset
  std::uniform_real_distribution<real_t> urd(0.5, 1);

  // Perturb the position of the new mid point so we don't get the same
  // insertion from the two patches
  real_t offset_scale = 1e-6;
  FacingType facing = orig_mesh.get_patch_facing(patch_id);

  for (auto e_vv : edge_to_split_v) {
    Vertex v1 = e_vv.first;
    Vertex v2 = e_vv.second;
    Edge split_e = mesh.find_edge(v1, v2);

    if (!split_e.is_valid()) {
      logger().warn("sew_components: Edge no longer exists {}, {}", v1, v2);
    }

    // Offset along the camera direction (FF: toward; BF: away)
    // with a random magnitude
    real_t extra_point_ratio = 0.5;
    Vector3f mid_pos = extra_point_ratio * vpositions[v1] +
                       (1 - extra_point_ratio) * vpositions[v2];
    Vector3f offset = (camera.position() - mid_pos).normalized();
    if (facing == FacingType::BACK)
      offset *= -1;
    real_t offset_rand = urd(eng);
    Vector3f extra_pos = mid_pos + offset_scale * offset * offset_rand;
    Vertex extra_v = mesh.add_vertex(extra_pos);

    if (!extra_v.is_valid()) {
      logger().error("sew_components: Cannot insert to edge {}, {}", v1, v2);
    } else {
      mesh.split(split_e, extra_v);
    }
  }
}

void stitch_mesh(Mesh const &mesh, Camera const &camera,
                 std::vector<Mesh *> const &patches, Mesh &stitched_mesh) {
  std::unordered_map<int, int> stitch_comp_b_v_to_new_v;
  stitched_mesh.add_vertex_property<int>("v:orig_idx", -1);
  auto stitched_orig_idx = stitched_mesh.vertex_property<int>("v:orig_idx");

  auto patchID = stitched_mesh.get_face_property<int>("f:patchID");
  if (!patchID)
    stitched_mesh.add_face_property<int>("f:patchID", -1);
  patchID = stitched_mesh.face_property<int>("f:patchID");

  for (auto const &p : patches) {
    auto p_patchID = p->get_face_property<int>("f:patchID");
    contess_assert_msg(p_patchID, "stitch_mesh: Patch ID is missing.");

    auto orig_boundary_idx = p->get_vertex_property<int>("v:orig_idx");
    contess_assert_msg(orig_boundary_idx,
                       "stitch_mesh: Index correspondence for "
                       "stitching is missing.");
    auto vpositions = p->get_vertex_property<Vector3f>("v:point");

    // This one is used to add faces
    std::unordered_map<int, int> p_v_to_new_v;

    for (size_t j = 0; j < p->n_vertices(); j++) {
      Vertex p_v(j);

      // This vertex needs to be stitched
      if (orig_boundary_idx[p_v] >= 0) {
        int comp_b_idx = orig_boundary_idx[p_v];
        if (!stitch_comp_b_v_to_new_v.count(comp_b_idx)) {
          Vertex v = stitched_mesh.add_vertex(vpositions[p_v]);
          stitch_comp_b_v_to_new_v[comp_b_idx] = v.idx();
          stitched_orig_idx[v] = comp_b_idx;
        }
        p_v_to_new_v[j] = stitch_comp_b_v_to_new_v[comp_b_idx];
      } else {
        Vertex v = stitched_mesh.add_vertex(vpositions[p_v]);
        p_v_to_new_v[j] = v.idx();
      }
    }

    int patch_id = p_patchID[*(p->faces_begin())];
    for (size_t j = 0; j < p->n_faces(); j++) {
      Vertex vv[3];
      p->verticesOfFace(Face(j), vv);

      contess_assert_msg(p_v_to_new_v.count(vv[0].idx()) &&
                             p_v_to_new_v.count(vv[1].idx()) &&
                             p_v_to_new_v.count(vv[2].idx()),
                         "stitch_mesh: Face contains unseen vertices.");

      std::vector<Vertex> new_f_vv({Vertex(p_v_to_new_v[vv[0].idx()]),
                                    Vertex(p_v_to_new_v[vv[1].idx()]),
                                    Vertex(p_v_to_new_v[vv[2].idx()])});
      Face f = stitched_mesh.add_face(new_f_vv);

      if (f.is_valid())
        patchID[f] = patch_id;
      if (!f.is_valid())
        logger().warn(
            "Patch {} face {} / {} has non-manifold edge(s): {}, {}, {} => "
            "{}, {}, {}",
            patch_id, Face(j), p->n_faces(), vv[0].idx(), vv[1].idx(),
            vv[2].idx(), new_f_vv[0].idx(), new_f_vv[1].idx(),
            new_f_vv[2].idx());
    }
  }

  // Create contour
  stitched_mesh.init();
  stitched_mesh.add_edge_property<real_t>("e:contour", -1);
  auto is_stitched_contour = stitched_mesh.edge_property<real_t>("e:contour");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  for (size_t i = 0; i < stitched_mesh.n_edges(); i++) {
    Edge e(i);
    Vertex orig_v1(stitched_orig_idx[stitched_mesh.vertex(e, 0)]);
    Vertex orig_v2(stitched_orig_idx[stitched_mesh.vertex(e, 1)]);

    if (!orig_v1.is_valid() || !orig_v2.is_valid())
      continue;

    Edge orig_e = mesh.find_edge(orig_v1, orig_v2);
    if (orig_e.is_valid())
      is_stitched_contour[e] = is_contour[orig_e];
  }

  // Assign patch id to boundary half edges
  auto is_boundary_edge =
      stitched_mesh.get_edge_property<bool>("e:is_boundary");
  if (!is_boundary_edge) {
    is_boundary_edge =
        stitched_mesh.add_edge_property<bool>("e:is_boundary", false);
  }
  is_boundary_edge = stitched_mesh.edge_property<bool>("e:is_boundary");
  auto patchBoundary =
      stitched_mesh.get_halfedge_property<int>("h:patchBoundary");
  if (!patchBoundary) {
    stitched_mesh.add_halfedge_property<int>("h:patchBoundary", -1);
  }
  patchBoundary = stitched_mesh.halfedge_property<int>("h:patchBoundary");
  patchBoundary.vector().assign(patchBoundary.vector().size(), -1);

  for (uint32_t i = 0; i != stitched_mesh.n_edges(); ++i) {
    Edge e(i);

    if (stitched_mesh.is_boundary(e))
      is_boundary_edge[e] = true;
    else
      is_boundary_edge[e] = false;

    Halfedge h0 = stitched_mesh.halfedge(e, 0);
    Halfedge h1 = stitched_mesh.halfedge(e, 1);

    Vertex orig_v1(stitched_orig_idx[stitched_mesh.vertex(e, 0)]);
    Vertex orig_v2(stitched_orig_idx[stitched_mesh.vertex(e, 1)]);

    // Interior edge
    if (!orig_v1.is_valid() || !orig_v2.is_valid())
      continue;

    // Patch or manifold boundary (would have correspondence to the original
    // mesh)
    bool is_orig_boundary = mesh.find_edge(orig_v1, orig_v2).is_valid() &&
                            mesh.is_boundary(mesh.find_edge(orig_v1, orig_v2));

    // test if this edge is on an open boundary
    if (!stitched_mesh.face(h0).is_valid() && is_orig_boundary) {
      patchBoundary[h1] = patchID[stitched_mesh.face(h1)];
    } else if (!stitched_mesh.face(h1).is_valid() && is_orig_boundary) {
      patchBoundary[h0] = patchID[stitched_mesh.face(h0)];
    } else if (stitched_mesh.face(h0).is_valid() &&
               stitched_mesh.face(h1).is_valid()) {
      // test if this edge is on a patch boundary
      int id0 = patchID[stitched_mesh.face(h0)];
      int id1 = patchID[stitched_mesh.face(h1)];
      if (id0 != id1 && is_stitched_contour[e] >= 0) {
        patchBoundary[h0] = id0;
        patchBoundary[h1] = id1;
      }
    }
  }

  // Assign VBO based on the patch VBO from the original mesh
  // and face orientation
  auto VBO_f = stitched_mesh.face_property<FacingType>("f:VBO_f");
  if (!VBO_f) {
    VBO_f = stitched_mesh.add_face_property<FacingType>("f:VBO_f");
  }
  VBO_f.vector().assign(VBO_f.vector().size(), FacingType::NA);
  auto VBO = stitched_mesh.face_property<FacingType>("f:VBO");
  if (!VBO) {
    VBO = stitched_mesh.add_face_property<FacingType>("f:VBO");
  }
  VBO.vector().assign(VBO.vector().size(), FacingType::NA);
  auto positions = stitched_mesh.get_vertex_property<Vector3f>("v:point");
  for (size_t i = 0; i < stitched_mesh.n_faces(); i++) {
    Face f(i);
    // The facing type based on look up
    int p_id = patchID[f];
    VBO[f] = mesh.get_patch_facing(p_id);

    auto f_norm = stitched_mesh.compute_face_normal(f);
    Vertex vs[3];
    stitched_mesh.verticesOfFace(f, vs);
    auto c_ray = positions[vs[0]] - camera.position();
    c_ray.normalized();
    VBO_f[f] = (f_norm.dot(c_ray) < 0) ? FacingType::FRONT : FacingType::BACK;
  }
}
