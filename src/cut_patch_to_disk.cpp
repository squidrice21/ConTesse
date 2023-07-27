// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "cut_patch_to_disk.h"
#include "common.h"
#include "surface_mesh/surface_mesh.h"

#include <igl/cut_to_disk.h>
#include <spdlog/fmt/ostr.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

int count_cut_degree(Mesh const &mesh, Vertex const &v, int patch_id) {
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto cut = mesh.get_edge_property<int>("e:disk_cut");
  if (!cut)
    return -1;
  size_t cut_degree = 0;
  auto hit = mesh.halfedges(v), hit_end = hit;
  do {
    Edge e = mesh.edge(*hit);
    if (cut[e] == patch_id || mesh.is_boundary(e) || is_contour[e] >= 0)
      cut_degree++;
  } while (++hit != hit_end);
  return cut_degree;
}

void cut_patch_to_disk_gu02(Mesh const &mesh, int patchID,
                            std::vector<int> &cut_halfedge_indices) {
  cut_halfedge_indices.clear();

  // make each patch homeomorphic to a disk by connecting boundary loops
  int start_index = 0;
  int id;
  for (id = 0; id < (int)mesh.get_const_patch_lengths().size() && id != patchID;
       id++) {
    start_index += mesh.get_const_patch_lengths()[id];
  }

  contess_assert_msg(id == patchID, "Nonexisting patch ID.");

  std::vector<std::vector<int32_t>> cuts;
  MatrixXi F;
  F.resize(mesh.get_const_patch_lengths()[id], 3);
  for (int f = 0; f < mesh.get_const_patch_lengths()[id]; f++) {
    F.row(f) =
        Vector3i(mesh.get_const_patch_indices().coeff(0, start_index + f),
                 mesh.get_const_patch_indices().coeff(1, start_index + f),
                 mesh.get_const_patch_indices().coeff(2, start_index + f));
  }
  igl::cut_to_disk(F, cuts);

  // tag edges accordingly
  for (auto &seam : cuts) {
    for (size_t i = 0; i + 1 < seam.size(); i++) {
      Vertex v1(seam[i]);
      Vertex v2(seam[i + 1]);

      Halfedge he = mesh.find_halfedge(v1, v2);
      cut_halfedge_indices.emplace_back(he.idx());
      cut_halfedge_indices.emplace_back(mesh.opposite_halfedge(he).idx());
    }
  }
}

void cut_patch_to_disk(Mesh &mesh, Camera const &camera, int id,
                       std::vector<int> &cut_halfedge_indices) {
  auto patchID = mesh.get_face_property<int>("f:patchID");
  contess_assert_msg(patchID, "Missing patches.");
  auto patch_removed = mesh.get_edge_property<bool>("e:patch_removed");
  auto cut_candidate = mesh.get_edge_property<bool>("e:cut_candidate");
  contess_assert_msg(cut_candidate,
                     "cut_patch_to_disk: Needs to label candidate edges.");

  std::unordered_map<int, std::vector<Face>> disk2faces;

  // 1. Initialize seed face book-keeping
  auto disk = mesh.face_property<int>("f:disk");
  if (!disk) {
    disk = mesh.add_face_property<int>("f:disk");
  }
  disk.vector().assign(disk.vector().size(), -1);

  // 2. Determine the seed faces
  bool seen_seed = false;
  for (int f_idx = 0; f_idx < (int)mesh.n_faces(); f_idx++) {
    Face face(f_idx);
    if (patchID[face] != id)
      continue;

    auto hit = mesh.halfedges(face), hit_end = hit;
    do {
      Edge e = mesh.edge(*hit);
      if (!cut_candidate[e] || (patch_removed && patch_removed[e])) {
        disk[face] = 0;
        seen_seed = true;
        break;
      }
    } while (++hit != hit_end);
  }

  // If there's no inconsistent faces
  if (!seen_seed)
    for (int f_idx = 0; f_idx < (int)mesh.n_faces(); f_idx++) {
      Face face(f_idx);
      if (patchID[face] != id)
        continue;
      disk[face] = 0;
      break;
    }

  // 3. Initialize the multiple seeds
  int disk_seed = 0;
  for (int f_idx = 0; f_idx < (int)mesh.n_faces(); f_idx++) {
    Face face(f_idx);
    if (patchID[face] != id)
      continue;
    if (disk[face] >= 0) {
      disk2faces[disk_seed] = std::vector<Face>({face});
      disk[face] = disk_seed++;
    }
  }
  // Delete the edges within the seeds
  // Add the seed boundary edges to the disk-growing queue
  auto cmp = [](std::pair<Edge, int> const &left,
                std::pair<Edge, int> const &right) -> bool {
    return left.second < right.second;
  };
  std::priority_queue<std::pair<Edge, int>, std::vector<std::pair<Edge, int>>,
                      decltype(cmp)>
      processEdges(cmp);

  std::unordered_set<int> deletededges;
  for (int f_idx = 0; f_idx < (int)mesh.n_faces(); f_idx++) {
    Face face(f_idx);
    if (patchID[face] != id || disk[face] < 0)
      continue;

    auto hit = mesh.halfedges(face), hit_end = hit;
    do {
      Edge e = mesh.edge(*hit);
      if (mesh.is_boundary(e))
        continue;

      int ndeleted = 0;
      if (disk[mesh.face(*hit)] >= 0)
        ndeleted++;
      if (disk[mesh.face(mesh.opposite_halfedge(*hit))] >= 0)
        ndeleted++;
      if (ndeleted == 1)
        processEdges.push(std::make_pair(e, 1));

      // Reached another disk
      if (ndeleted == 2 && (disk[mesh.face(*hit)] !=
                            disk[mesh.face(mesh.opposite_halfedge(*hit))])) {
        // Always handle disk-disk intersections with a non-candidate edge first
        if (cut_candidate[e] && (!patch_removed || !patch_removed[e]))
          processEdges.push(std::make_pair(e, 0));
        // Handle the candidate edge last
        else
          processEdges.push(std::make_pair(e, 2));
      }
    } while (++hit != hit_end);
  }

  // 4. Assuming there's only a single connected component
  // delete all faces adjacent to edges with exactly one adjacent face
  auto propagate_seed = [&](int seed_from, int seed_to) {
    contess_assert_msg(disk2faces.count(seed_from) && disk2faces.count(seed_to),
                       "Error in disk lookup.");
    for (auto f : disk2faces[seed_from]) {
      disk[f] = seed_to;
    }
    // Update lookup map
    disk2faces[seed_to].insert(disk2faces[seed_to].end(),
                               disk2faces[seed_from].begin(),
                               disk2faces[seed_from].end());
    disk2faces.erase(seed_from);
  };

  size_t itr = 0;
  while (!processEdges.empty()) {
    Edge nexte = processEdges.top().first;

    processEdges.pop();

    // All edges in processEdges are valid edges that can be deleted
    Face todelete;
    std::vector<Face> deleted;
    std::vector<Face> adj_faces({mesh.face(nexte, 0), mesh.face(nexte, 1)});
    for (auto f : adj_faces) {
      if (patchID[f] != id)
        continue;
      if (disk[f] < 0)
        todelete = f;
      else
        deleted.emplace_back(f);
    }

    // Merge disks
    if (deleted.size() == 2 &&
        disk[mesh.face(nexte, 0)] != disk[mesh.face(nexte, 1)]) {
      deletededges.insert(nexte.idx());
      int seed_from = disk[deleted.back()], seed_to = disk[deleted.front()];

      // Propagate the seed if both sides are from different disks
      propagate_seed(seed_from, seed_to);
    }

    // Grow the disk
    if (todelete.is_valid() && deleted.size() == 1) {
      deletededges.insert(nexte.idx());

      int disk_idx = disk[deleted.front()];
      disk[todelete] = disk_idx;
      disk2faces[disk_idx].emplace_back(todelete);

      auto hit = mesh.halfedges(todelete), hit_end = hit;
      do {
        Edge e = mesh.edge(*hit);
        if (mesh.is_boundary(e) ||
            patchID[mesh.face(*hit)] !=
                patchID[mesh.face(mesh.opposite_halfedge(*hit))] ||
            deletededges.count(e.idx()))
          continue;

        int ndeleted = 0;
        if (disk[mesh.face(*hit)] >= 0)
          ndeleted++;
        if (disk[mesh.face(mesh.opposite_halfedge(*hit))] >= 0)
          ndeleted++;
        if (ndeleted == 1)
          processEdges.push(std::make_pair(e, 1));

        // Reached another disk
        // Always handle disk-disk intersections first
        if (ndeleted == 2 && (disk[mesh.face(*hit)] !=
                              disk[mesh.face(mesh.opposite_halfedge(*hit))])) {
          // Always handle disk-disk intersections with a non-candidate edge
          // first
          if (cut_candidate[e] && (!patch_removed || !patch_removed[e]))
            processEdges.push(std::make_pair(e, 0));
          // Handle the candidate edge last
          else
            processEdges.push(std::make_pair(e, 2));
        }
      } while (++hit != hit_end);
    }

    itr++;
  }

  // 5. Prune spines
  std::deque<Vertex> processVertices;
  auto is_cut = [&](Edge const &e) -> bool {
    // Boundary is always in cut graph
    if (!mesh.face(e, 0).is_valid() || !mesh.face(e, 1).is_valid())
      return true;
    if (patchID[mesh.face(e, 0)] != patchID[mesh.face(e, 1)])
      return true;
    if (deletededges.count(e.idx()) ||
        (patchID[mesh.face(e, 0)] != id && patchID[mesh.face(e, 1)] != id))
      return false;
    return true;
  };
  auto count_cut_degree = [&](Vertex const &v) -> size_t {
    size_t cut_degree = 0;
    auto hit = mesh.halfedges(v), hit_end = hit;
    do {
      if (is_cut(mesh.edge(*hit)))
        cut_degree++;
    } while (++hit != hit_end);
    return cut_degree;
  };
  // Init the queue with cut-graph-degree-1 vertices
  for (size_t i = 0; i < mesh.n_vertices(); i++) {
    Vertex v(i);

    size_t cut_degree = count_cut_degree(v);
    if (cut_degree == 1)
      processVertices.push_back(v);
  }
  while (!processVertices.empty()) {
    Vertex next_v = processVertices.front();
    processVertices.pop_front();

    auto hit = mesh.halfedges(next_v), hit_end = hit;
    do {
      Edge e = mesh.edge(*hit);

      if (is_cut(e)) {
        contess_assert_msg(patchID[mesh.face(e, 0)] == patchID[mesh.face(e, 1)],
                           "Can't delete boundary cuts.");
        deletededges.emplace(e.idx());

        // If the deletion can continue
        Vertex n_v = mesh.to_vertex(*hit);
        size_t cut_degree = count_cut_degree(n_v);
        if (cut_degree == 1)
          processVertices.push_back(n_v);

        // All vertices in queueh has cut graph degree of 1, so we are done here
        break;
      }
    } while (++hit != hit_end);
  }

  // 6. Label not deleted edges
  for (size_t i = 0; i < mesh.n_edges(); i++) {
    Edge e(i);

    // Skip deleted edges or edges outside the patch
    if (deletededges.count(e.idx()) ||
        (!mesh.face(e, 1).is_valid() && patchID[mesh.face(e, 0)] != id) ||
        (!mesh.face(e, 0).is_valid() && patchID[mesh.face(e, 1)] != id) ||
        (mesh.face(e, 0).is_valid() && mesh.face(e, 1).is_valid() &&
         patchID[mesh.face(e, 0)] != id && patchID[mesh.face(e, 1)] != id))
      continue;

    Halfedge he = mesh.halfedge(e, 0);
    cut_halfedge_indices.emplace_back(he.idx());
    cut_halfedge_indices.emplace_back(mesh.opposite_halfedge(he).idx());
  }
}

void cut_patch_to_disk(Mesh &mesh, Camera const &camera, bool use_gu02,
                       bool ff_only) {
  auto cut_candidate = mesh.get_edge_property<bool>("e:cut_candidate");
  auto patch_removed = mesh.get_edge_property<bool>("e:patch_removed");
  contess_assert_msg(cut_candidate,
                     "cut_patch_to_disk: Needs to label candidate edges.");
  // Store the cut graph in the mesh
  auto cut = mesh.get_edge_property<int>("e:disk_cut");
  if (!cut) {
    cut = mesh.add_edge_property<int>("e:disk_cut", -1);
  } else
    cut = mesh.edge_property<int>("e:disk_cut", -1);

  bool to_init = true;
  for (int id = 0; id < (int)mesh.get_const_patch_lengths().size(); id++) {
    std::vector<int> cut_halfedge_indices;

    auto facing = mesh.get_patch_facing(id);
    if (ff_only && facing != FacingType::FRONT)
      continue;

    if (use_gu02) {
      cut_patch_to_disk_gu02(mesh, id, cut_halfedge_indices);
    } else {
      cut_patch_to_disk(mesh, camera, id, cut_halfedge_indices);
    }

    if (to_init) {
      cut.vector().assign(cut.vector().size(), -1);
      to_init = false;
    }

    for (auto he_idx : cut_halfedge_indices) {
      Edge e = mesh.edge(Halfedge(he_idx));
      cut[e] = id;
    }
  }

  // Verify if all cut edges are validate
  auto patchID = mesh.get_face_property<int>("f:patchID");
  contess_assert_msg(patchID, "Missing patches.");

  for (size_t i = 0; i < mesh.n_edges(); i++) {
    Edge e(i);

    // Skip boundary (manifold/patch) edges
    if (!mesh.face(e, 1).is_valid() || !mesh.face(e, 0).is_valid() ||
        (mesh.face(e, 0).is_valid() && mesh.face(e, 1).is_valid() &&
         patchID[mesh.face(e, 0)] != patchID[mesh.face(e, 1)]))
      continue;

    if (cut[e] < 0)
      continue;
  }
}
