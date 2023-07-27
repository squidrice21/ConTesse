// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "remove_wso_failure.h"

#include "chain_contour.h"
#include "common.h"
#include "remove_cut_infeasible_disk.h"
#include <spdlog/fmt/ostr.h>
#include <unordered_set>
#include <vector>

void get_wso_failure_removal(Mesh &mesh, Camera const &camera,
                             std::unordered_set<int> const &wso_failures,
                             std::vector<std::vector<int>> &patch_removals,
                             real_t max_loop_length, bool ff_only) {
  auto patchID = mesh.get_face_property<int>("f:patchID");
  bool seen_chains = !mesh.get_patch_chains().empty();

  if (mesh.get_patch_chains().empty())
    chain_contour(mesh, camera, true, ff_only);

  // Two removal lists: one aggressive, one conservative
  patch_removals.resize(2);

  // 1. Gather boundary halfedges (indicating the number of holes)
  // Neighboring relationship between patches
  struct PatchInfo {
    int patch_id;
    std::unordered_map<int, Halfedge> neighbors;
    bool is_disk;
    Halfedge self_halfedge;
  };
  std::unordered_map<int, PatchInfo> patch_infos;
  for (auto const &patch : mesh.get_patch_chains()) {
    int patch_id = patch.first;
    if (!wso_failures.count(patch_id))
      continue;

    logger().info("Patch: {}", patch_id);
    PatchInfo info;
    info.patch_id = patch_id;

    auto chains = mesh.get_patch_chains().equal_range(patch_id);
    for (auto chain = chains.first; chain != chains.second; ++chain) {
      for (auto const &b_cut_he : *chain->second) {
        Halfedge he = b_cut_he;
        if (!mesh.face(he).is_valid() || patchID[mesh.face(he)] != patch_id)
          he = mesh.opposite_halfedge(he);

        if (!mesh.face(he).is_valid() || patchID[mesh.face(he)] != patch_id)
          continue;

        // Save a halfedge (this is for the patch with a single loop, which is
        // purely a manifold boundary loop)
        if (!info.self_halfedge.is_valid())
          info.self_halfedge = he;
        Halfedge oppo = mesh.opposite_halfedge(he);
        if (!oppo.is_valid() || !mesh.face(oppo).is_valid())
          continue;
        int oppo_patch_id = patchID[mesh.face(oppo)];

        if (oppo_patch_id < 0 || oppo_patch_id == patch_id)
          continue;
        if (!info.neighbors.count(oppo_patch_id)) {
          logger().info("\tNeighbor: {} - {}", patchID[mesh.face(he)],
                        oppo_patch_id);
          info.neighbors[oppo_patch_id] = oppo;
        }
      }
    }

    patch_infos[patch_id] = info;
  }

  // 2. Check if each patch is a disk by itself
  // If so => Add to both removal lists
  for (auto &p : patch_infos) {
    // If the patch has a single loop, which is
    // purely a manifold boundary loop
    // Or the patch has a single boundary loop
    if (p.second.neighbors.empty() || p.second.neighbors.size() == 1) {
      p.second.is_disk = is_patch_removeable(
          mesh, camera, p.second.self_halfedge, max_loop_length);
      if (p.second.is_disk) {
        logger().info("Patch {}: Disk", p.first);
        patch_removals[0].emplace_back(p.first);
        patch_removals[1].emplace_back(p.first);
      }
    } else {
      p.second.is_disk = false;
    }
  }

  // 3. For the non-disk patch, add all its removeable neighbors to the second
  // list (aggressive one). If none of the neighbors is in the first list, also
  // add them to the first one.
  for (auto const &p : patch_infos) {
    if (!p.second.is_disk) {
      std::unordered_set<int> removeable_neighbors;
      for (auto const &n : p.second.neighbors) {
        bool is_removeable =
            is_patch_removeable(mesh, camera, n.second, max_loop_length);
        if (is_removeable)
          removeable_neighbors.emplace(n.first);
      }

      bool seen_neighbor = false;
      for (auto pid : patch_removals[0])
        if (removeable_neighbors.count(pid)) {
          seen_neighbor = true;
          break;
        }
      if (seen_neighbor) {
        for (auto const &n : removeable_neighbors) {
          logger().info("Patch {} - {}: conservative", p.first, n);
          if (std::find(patch_removals[1].begin(), patch_removals[1].end(),
                        n) == patch_removals[1].end())
            patch_removals[1].emplace_back(n);
        }
      } else {
        for (auto const &n : removeable_neighbors) {
          logger().info("Patch {} - {}: aggressive", p.first, n);
          if (std::find(patch_removals[0].begin(), patch_removals[0].end(),
                        n) == patch_removals[0].end())
            patch_removals[0].emplace_back(n);
          if (std::find(patch_removals[1].begin(), patch_removals[1].end(),
                        n) == patch_removals[1].end())
            patch_removals[1].emplace_back(n);
        }
      }
    }
  }

  if (!seen_chains) {
    mesh.get_patch_chains().clear();
    mesh.get_oriented_chains().clear();
  }
}

bool remove_wso_failure(Mesh &mesh, Camera const &camera,
                        std::vector<int> const &patch_removals,
                        std::unordered_set<int> &affected_patches,
                        bool ff_only) {
  auto patchID = mesh.get_face_property<int>("f:patchID");
  if (patch_removals.empty())
    return false;

  auto patch_removed = mesh.get_edge_property<bool>("e:patch_removed");
  if (!patch_removed) {
    patch_removed = mesh.add_edge_property<bool>("e:patch_removed", false);
  }
  patch_removed.vector().assign(patch_removed.vector().size(), false);

  bool removed = false;

  if (mesh.get_patch_chains().empty())
    chain_contour(mesh, camera, true, ff_only);

  for (auto const &patch : mesh.get_patch_chains()) {
    int patch_id = patch.first;

    if (std::find(patch_removals.begin(), patch_removals.end(), patch_id) ==
        patch_removals.end())
      continue;

    logger().info("Remove Patch: {}", patch_id);

    size_t chain_count = 0;
    auto chains = mesh.get_patch_chains().equal_range(patch_id);
    for (auto chain = chains.first; chain != chains.second; ++chain) {
      // Move to a halfedge with two neighboring patches
      Halfedge remove_he;
      for (auto const &b_cut_he : *chain->second) {
        Halfedge he = b_cut_he;
        if (!mesh.face(he).is_valid() || patchID[mesh.face(he)] != patch_id)
          he = mesh.opposite_halfedge(he);

        if (!mesh.face(he).is_valid() || patchID[mesh.face(he)] != patch_id)
          continue;

        Halfedge oppo = mesh.opposite_halfedge(he);
        if (!oppo.is_valid() || !mesh.face(oppo).is_valid())
          continue;
        int oppo_patch_id = patchID[mesh.face(oppo)];

        if (oppo_patch_id < 0 || oppo_patch_id == patch_id)
          continue;

        remove_he = he;
        break;
      }
      if (remove_he.is_valid()) {
        int to_patch_id = patchID[mesh.face(mesh.opposite_halfedge(remove_he))];
        affected_patches.emplace(to_patch_id);
        remove_patch(mesh, remove_he);
        removed = true;
      } else
        logger().warn("remove_wso_failure: Cannot remove patch {}. It's a "
                      "purely manifold boundary patch.",
                      patch_id);
      chain_count++;
    }

    // Check afterward
    contess_assert(chain_count == 1);
  }

  mesh.get_patch_chains().clear();
  mesh.get_oriented_chains().clear();

  return removed;
}
