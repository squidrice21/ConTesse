// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "chain_contour.h"
#include "common.h"

#include <algorithm>
#include <igl/predicates/predicates.h>
#include <memory>
#include <spdlog/fmt/ostr.h>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

bool chain_contour(Mesh &mesh, Camera const &camera, bool to_cache_chain,
                   bool ff_only) {
  // To avoid misorientation after offseting the flipped loops, we save and
  // reuse the original chaining orientaion in h:chained_orientation.

  // First, we trace from interior edges. Each time it traces two loops
  // bounding the two patches adjacent to the starting edge.
  // Then, we trace from boundary edges (boundary edges may not be traced is a
  // patch is adjacent to more than one patch).
  // Tracing boundary is different from tracing interior edges since we may not
  // have the correct adjacent patch when tracing the target orientation. So
  // we'll have to trace the flipped orientation first then flip the loop.
  std::unordered_set<int> edge_visited;
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");
  auto VBO_f = mesh.get_face_property<FacingType>("f:VBO_f");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto patchBoundary = mesh.get_halfedge_property<int>("h:patchBoundary");
  auto chained_orientation =
      mesh.get_halfedge_property<int>("h:chained_orientation");
  bool existing_chained = (chained_orientation);
  if (!chained_orientation && to_cache_chain) {
    chained_orientation =
        mesh.add_halfedge_property<int>("h:chained_orientation", -1);
  }

  contess_assert_msg(patchBoundary,
                     "Extracted patches is needed for chaining contours.");

  mesh.m_oriented_chains.clear();
  mesh.m_patch_chains.clear();

  auto trace_chain =
      [&edge_visited, &patchBoundary, &is_contour, &chained_orientation](
          Mesh &mesh, Halfedge const &he, std::shared_ptr<Chain> chain,
          int adj_patch_idx, bool patch_on_left,
          std::map<int, Halfedge> &oppo_patch, bool &successful) -> bool {
    successful = true;
    bool seen_manifold_boundary = false;
    oppo_patch.clear();

    std::unordered_set<int>
        seen_vertices; // Safety check (if tracing is trapped)
    std::unordered_map<int, int> temp_orientation;

    // 2. Trace the chain
    auto from_he = he;
    Vertex des = mesh.to_vertex(he);
    Vertex cur = des;
    do {
      edge_visited.emplace(mesh.edge(from_he).idx());

      {
        chain->emplace_back(from_he);
        temp_orientation[from_he.idx()] = 1;
        if (mesh.opposite_halfedge(from_he).is_valid())
          temp_orientation[mesh.opposite_halfedge(from_he).idx()] = 0;
      }
      if (mesh.is_boundary(mesh.edge(from_he)))
        seen_manifold_boundary = true;

      auto hit = mesh.halfedges(cur);
      auto hit_end = hit;
      Halfedge next_he;
      std::map<int, Halfedge> left_patch_to_he;
      do {
        if (patch_on_left && patchBoundary[*hit] != adj_patch_idx)
          continue;
        if (!patch_on_left &&
            patchBoundary[mesh.opposite_halfedge(*hit)] != adj_patch_idx)
          continue;
        if (mesh.opposite_halfedge(*hit) == from_he)
          continue;
        if (is_contour[mesh.edge(*hit)] < 0 &&
            !mesh.is_boundary(mesh.edge(*hit)))
          continue;

        next_he = *hit;
        break;
      } while (++hit != hit_end);

      if (!next_he.is_valid())
        logger().error("chain_contour: Can't find next step at {}.", cur);

      contess_assert_msg(next_he.is_valid(),
                         "chain_contour: Can't find next step.");
      if (seen_vertices.count(mesh.from_vertex(next_he).idx())) {
        successful = false;
        return seen_manifold_boundary;
      }
      contess_assert_msg(!seen_vertices.count(mesh.from_vertex(next_he).idx()),
                         "chain_contour: Tracing is trapped.");
      seen_vertices.emplace(mesh.from_vertex(next_he).idx());

      // Book-keeping for the opposite side patch
      int patchID = patchBoundary[mesh.opposite_halfedge(next_he)];
      if (mesh.is_boundary(mesh.edge(next_he)))
        patchID = -1;
      oppo_patch[patchID] = next_he;

      cur = mesh.to_vertex(next_he);
      from_he = next_he;
    } while (cur != des);

    // We don't record the orientation if it's a fully open manifold boundary
    // loop
    if (!(oppo_patch.count(-1) && oppo_patch.size() == 1) &&
        chained_orientation) {
      for (auto const &ori : temp_orientation) {
        chained_orientation[Halfedge(ori.first)] = ori.second;
      }
    }

    return seen_manifold_boundary;
  };

  // First trace (fully or partially) patch boundary loops
  for (uint32_t i = 0; i != mesh.n_edges(); ++i) {
    Edge e(i);

    Halfedge he = mesh.halfedge(e, 0);

    if (is_contour[e] < 0 || mesh.is_boundary(e) ||
        edge_visited.find(e.idx()) != edge_visited.end())
      continue;

    bool seen_manifold_boundary = false;
    std::shared_ptr<Chain> chain = std::make_shared<Chain>();

    // 1. Find the orientation of the chain
    // CCW for the front portion

    // This edge has been modified. Skip.
    if (existing_chained && chained_orientation &&
        ((he.is_valid() && chained_orientation[he] < 0) ||
         (mesh.opposite_halfedge(he).is_valid() &&
          chained_orientation[mesh.opposite_halfedge(he)] < 0))) {
      continue;
    }

    if (!existing_chained) {
      bool to_flip = false;
      Face f_face = mesh.face(he);
      contess_assert_msg(f_face.is_valid(),
                         "chain_contour: Not a patch boundary edge.");

      if (VBO[f_face] != FacingType::FRONT) {
        he = mesh.halfedge(e, 1);
        f_face = mesh.face(he);
      }

      if (!(VBO[f_face] == FacingType::FRONT &&
            VBO[f_face] != VBO[mesh.face(mesh.opposite_halfedge(he))])) {
        logger().warn(
            "VBO wrongly labeled. Faces: " + std::to_string(f_face.idx()) +
            ": " + std::to_string(VBO[f_face]) + ", " +
            std::to_string(mesh.face(mesh.opposite_halfedge(he)).idx()) + ": " +
            std::to_string(VBO[mesh.face(mesh.opposite_halfedge(he))]) + ".");
        continue;
      }

      // CCW: the opposite vertex should be on the left of the halfedge
      if (igl::predicates::orient3d(
              camera.position(), vpositions[mesh.from_vertex(he)],
              vpositions[mesh.to_vertex(he)],
              vpositions[mesh.to_vertex(mesh.next_halfedge(he))]) ==
          igl::predicates::Orientation::NEGATIVE) {
        to_flip = true;
      }

      // If the adjacent face is inconsistent, then flip
      if (VBO[f_face] != VBO_f[f_face])
        to_flip = !to_flip;

      if (to_flip) {
        he = mesh.opposite_halfedge(he);
      }

      if (chained_orientation) {
        chained_orientation[he] = 1;
        if (mesh.opposite_halfedge(he).is_valid())
          chained_orientation[mesh.opposite_halfedge(he)] = 0;
      }
    }
    // Read the previous orientation to avoid misorientation caused by mesh
    // offseting
    else {
      if (chained_orientation && chained_orientation[he] == 0)
        he = mesh.opposite_halfedge(he);
    }
    ////

    if (ff_only &&
        mesh.get_patch_facing(patchBoundary[he]) != FacingType::FRONT)
      continue;

    std::map<int, Halfedge> oppo_patch;
    bool successful;
    seen_manifold_boundary = trace_chain(mesh, he, chain, patchBoundary[he],
                                         true, oppo_patch, successful);
    if (!successful)
      return false;

    int valid_oppo_patch_count = 0;
    for (auto const &oppo_p : oppo_patch)
      if (oppo_p.first >= 0)
        valid_oppo_patch_count++;

    // This is a pure open boundary loop
    if (valid_oppo_patch_count == 0)
      continue;

    // 3. Book keeping
    mesh.m_oriented_chains.emplace_back(chain);
    mesh.m_patch_chains.insert(std::make_pair(patchBoundary[he], chain));

    // 4. Either keep the same the loop for the patch on the opposite side
    // Or if the current loop contains manifold boundary, retrace based on the
    // right side of the model
    if (!seen_manifold_boundary && valid_oppo_patch_count == 1) {
      // Save the loop with the FF orientation even for BF patch since
      // WSO can only handle FF orientation
      mesh.m_oriented_chains.emplace_back(chain);
      mesh.m_patch_chains.insert(
          std::make_pair(patchBoundary[mesh.opposite_halfedge(he)], chain));
    } else if (seen_manifold_boundary || valid_oppo_patch_count > 1) {
      // Need to rerun the tracing with a different patch index
      for (auto const &oppo_p : oppo_patch)
        if (oppo_p.first >= 0) {
          Halfedge he_oppo = oppo_p.second;

          std::map<int, Halfedge> oppo_patch_oppo;
          chain = std::make_shared<Chain>();
          bool successful;
          trace_chain(mesh, he_oppo, chain,
                      patchBoundary[mesh.opposite_halfedge(he_oppo)], false,
                      oppo_patch_oppo, successful);
          if (!successful)
            return false;

          mesh.m_oriented_chains.emplace_back(chain);
          mesh.m_patch_chains.insert(std::make_pair(
              patchBoundary[mesh.opposite_halfedge(he_oppo)], chain));
        }
    }
  }

  // Then trace the fully open boundary loops and remainig partially open loops
  for (uint32_t i = 0; i != mesh.n_edges(); ++i) {
    Edge e(i);

    if (!mesh.is_boundary(e) ||
        edge_visited.find(e.idx()) != edge_visited.end())
      continue;

    bool seen_manifold_boundary = false;
    std::shared_ptr<Chain> chain = std::make_shared<Chain>();

    // 1. Find the orientation of the chain
    // CCW for the front portion
    bool to_flip = false;
    Halfedge he = mesh.halfedge(e, 0);
    Face f_face = mesh.face(he);
    if (!f_face.is_valid()) {
      he = mesh.halfedge(e, 1);
      f_face = mesh.face(he);
    }

    auto oppo_he = mesh.opposite_halfedge(he);
    if (existing_chained && chained_orientation &&
        ((he.is_valid() && chained_orientation[he] < 0) ||
         (oppo_he.is_valid() && chained_orientation[oppo_he] < 0))) {
      logger().warn("Skip boundary patch: {}", patchBoundary[he]);
      continue;
    }

    if (!existing_chained) {
      bool adj_is_ff = true;

      adj_is_ff = VBO[f_face] == FacingType::FRONT;

      // CCW: the opposite vertex should be on the left of the halfedge
      if (adj_is_ff &&
          igl::predicates::orient3d(
              camera.position(), vpositions[mesh.from_vertex(he)],
              vpositions[mesh.to_vertex(he)],
              vpositions[mesh.to_vertex(mesh.next_halfedge(he))]) ==
              igl::predicates::Orientation::NEGATIVE) {
        to_flip = true;
      } else if (!adj_is_ff &&
                 igl::predicates::orient3d(
                     camera.position(), vpositions[mesh.from_vertex(he)],
                     vpositions[mesh.to_vertex(he)],
                     vpositions[mesh.to_vertex(mesh.next_halfedge(he))]) ==
                     igl::predicates::Orientation::POSITIVE) {
        to_flip = true;
      }
      // If the adjacent face is inconsistent, then flip
      if (VBO[f_face] != VBO_f[f_face])
        to_flip = !to_flip;
      if (!adj_is_ff)
        to_flip = !to_flip;
    } else if (chained_orientation) {
      to_flip = (chained_orientation[he] == 0);
    }

    std::map<int, Halfedge> oppo_patch;
    bool successful = true;
    seen_manifold_boundary = trace_chain(mesh, he, chain, patchBoundary[he],
                                         true, oppo_patch, successful);

    if (!successful)
      return false;

    int valid_oppo_patch_count = 0;
    for (auto const &oppo_p : oppo_patch)
      if (oppo_p.first >= 0)
        valid_oppo_patch_count++;

    // 3. Book keeping
    // By default the open boundary loops need to be flipped
    if (to_flip) {
      std::reverse(chain->begin(), chain->end());
      for (size_t j = 0; j < chain->size(); j++) {
        auto chain_he = (*chain)[j];
        Edge orig_e = mesh.edge(chain_he);
        if (chained_orientation) {
          if (mesh.is_boundary(orig_e)) {
            chained_orientation[chain_he] = 0;
            if (mesh.opposite_halfedge(chain_he).is_valid())
              chained_orientation[mesh.opposite_halfedge(chain_he)] = 1;
          }
        }

        (*chain)[j] = mesh.opposite_halfedge((*chain)[j]);
      }
    } else {
      for (size_t j = 0; j < chain->size(); j++) {
        auto chain_he = (*chain)[j];
        Edge orig_e = mesh.edge(chain_he);
        if (chained_orientation) {
          if (mesh.is_boundary(orig_e)) {
            chained_orientation[chain_he] = 1;
            if (mesh.opposite_halfedge(chain_he).is_valid())
              chained_orientation[mesh.opposite_halfedge(chain_he)] = 0;
          }
        }
      }
    }

    mesh.m_oriented_chains.emplace_back(chain);
    mesh.m_patch_chains.insert(std::make_pair(patchBoundary[he], chain));
  }

  return true;
}

bool chain_cut_graph(Mesh &mesh, Camera const &camera, bool ff_only) {
  // 0. Needs to be cut already
  auto cut = mesh.get_edge_property<int>("e:disk_cut");
  contess_assert_msg(cut, "Cuts are needed for chaining cut graph.");
  auto patchID = mesh.get_face_property<int>("f:patchID");
  contess_assert_msg(patchID,
                     "Extracted patches is needed for chaining contours.");
  auto VBO = mesh.get_face_property<FacingType>("f:VBO");
  auto VBO_f = mesh.get_face_property<FacingType>("f:VBO_f");
  std::unordered_map<int, FacingType> patch_facing;

  // 1. Double the non-bounary cuts for Eulerization
  auto cut_directed = mesh.halfedge_property<bool>("h:cut");
  cut_directed.vector().assign(cut_directed.vector().size(), false);
  std::unordered_map<int, size_t> patch_he_count;

  for (size_t e_idx = 0; e_idx < mesh.n_edges(); e_idx++) {
    Edge e(e_idx);

    // Ignore boundary cuts for now
    if ((!mesh.face(e, 0).is_valid() || !mesh.face(e, 1).is_valid()) ||
        patchID[mesh.face(e, 0)] != patchID[mesh.face(e, 1)]) {
      if (mesh.face(e, 0).is_valid()) {
        if (!patch_he_count.count(patchID[mesh.face(e, 0)]))
          patch_he_count[patchID[mesh.face(e, 0)]] = 0;
        patch_he_count[patchID[mesh.face(e, 0)]]++;
        if (VBO[mesh.face(e, 0)] == VBO_f[mesh.face(e, 0)])
          patch_facing[patchID[mesh.face(e, 0)]] = VBO[mesh.face(e, 0)];
      }
      if (mesh.face(e, 1).is_valid()) {
        if (!patch_he_count.count(patchID[mesh.face(e, 1)]))
          patch_he_count[patchID[mesh.face(e, 1)]] = 0;
        patch_he_count[patchID[mesh.face(e, 1)]]++;

        if (VBO[mesh.face(e, 1)] == VBO_f[mesh.face(e, 1)])
          patch_facing[patchID[mesh.face(e, 1)]] = VBO[mesh.face(e, 1)];
      }

      continue;
    }

    if (cut[e] < 0)
      continue;

    cut_directed[mesh.halfedge(e, 0)] = true;
    cut_directed[mesh.halfedge(e, 1)] = true;

    // The non-boundary cut is assumed to be adjacent to two consistent faces
    patch_facing[patchID[mesh.face(e, 0)]] = VBO[mesh.face(e, 0)];

    if (!patch_he_count.count(patchID[mesh.face(e, 0)]))
      patch_he_count[patchID[mesh.face(e, 0)]] = 0;
    patch_he_count[patchID[mesh.face(e, 0)]] += 2;
  }

  // 2. Chain contour with the correct orientation
  if (mesh.get_oriented_chains().empty()) {
    bool successful = chain_contour(mesh, camera, true, ff_only);
    if (!successful)
      return false;
  }

  // 3. Trace an Eulerian path with Hierholzer’s Algorithm
  std::multimap<int, std::shared_ptr<Chain>> patch_chain;
  std::unordered_set<int> visited_patch;
  for (auto const &patch : mesh.m_patch_chains) {
    if (ff_only && mesh.get_patch_facing(patch.first) != FacingType::FRONT)
      continue;

    if (visited_patch.count(patch.first))
      continue;

    visited_patch.insert(patch.first);

    // Set boundary cut directions
    auto chains = mesh.m_patch_chains.equal_range(patch.first);
    for (auto chain = chains.first; chain != chains.second; ++chain) {
      for (auto const &b_cut_he : *chain->second) {
        cut_directed[b_cut_he] = true;
      }
    }

    auto is_within_patch = [&patchID](Mesh const &mesh, Halfedge const &he,
                                      int patch_idx) -> bool {
      if (!mesh.face(he).is_valid())
        return patchID[mesh.face(mesh.opposite_halfedge(he))] == patch_idx;
      if (!mesh.face(mesh.opposite_halfedge(he)).is_valid())
        return patchID[mesh.face(he)] == patch_idx;

      return patchID[mesh.face(he)] == patch_idx ||
             patchID[mesh.face(mesh.opposite_halfedge(he))] == patch_idx;
    };

    // Trace within the patch
    std::vector<Vertex> final_path;
    // Vertex start_v = mesh.from_vertex(patch.second->front());
    Vertex start_v;

    // Choose the starting vertex to be a vertex with degree of two
    // (aka a patch boundary vertex not adjacent to cut path).
    // So we can trace all vertices with one go
    {
      for (auto const &he : *patch.second) {
        Vertex v = mesh.from_vertex(he);
        auto hit = mesh.halfedges(v), hit_end = hit;

        // Count in and out degree
        size_t in_degree = 0, out_degree = 0;
        do {
          if (cut_directed[*hit] && is_within_patch(mesh, *hit, patch.first))
            out_degree++;

          Halfedge oppo_he = mesh.opposite_halfedge(*hit);
          if (cut_directed[oppo_he] &&
              is_within_patch(mesh, oppo_he, patch.first))
            in_degree++;
        } while (++hit != hit_end);

        // Assume even degree
        contess_assert_msg(
            in_degree == out_degree,
            "Vertex does not have matched degrees: " + std::to_string(v.idx()) +
                " - " + std::to_string(in_degree) + " vs " +
                std::to_string(out_degree) + " in patch " +
                std::to_string(patch.first));

        if (in_degree == 1) {
          start_v = v;
          break;
        }
      }

      contess_assert_msg(start_v.is_valid(),
                         "Cannot find a degree-2 starting vertex in patch " +
                             std::to_string(patch.first));
    }

    int remain_degree = patch_he_count[patch.first];
    std::vector<Vertex> candidate_start_v;
    final_path.reserve(remain_degree + 1);
    final_path.emplace_back(start_v);
    std::vector<Vertex> path;
    path.reserve(remain_degree + 1);
    auto is_step_valid = [&](Vertex const &cur, Vertex const &next) -> bool {
      std::vector<Halfedge> he_remain;
      auto hit = mesh.halfedges(cur), hit_end = hit;
      do {
        if (!cut_directed[*hit] || !is_within_patch(mesh, *hit, patch.first))
          continue;

        // Try to delete
        if (mesh.from_vertex(*hit) == cur && mesh.to_vertex(*hit) == next)
          continue;

        he_remain.emplace_back(*hit);
        if (cut_directed[mesh.opposite_halfedge(*hit)])
          he_remain.emplace_back(mesh.opposite_halfedge(*hit));
      } while (++hit != hit_end);

      // This step results in a back trace in the future.
      // We can't back trace on the doubled edges (not 2-manifold)
      if (he_remain.size() == 2 &&
          mesh.opposite_halfedge(he_remain.front()) == he_remain.back())
        return false;
      return true;
    };
    while (remain_degree > 0) {
      path.clear();

      Vertex cur_v = start_v;
      do {
        std::vector<Vertex> next_v;

        // This circulation is always CCW wrt the surface
        auto hit = mesh.halfedges(cur_v), hit_end = hit;
        // Move to the last step
        if (!path.empty()) {
          do {
            if (mesh.to_vertex(*hit) != path.back())
              continue;
            break;
          } while (++hit != hit_end);
          hit_end = hit;
        }

        do {
          if (!cut_directed[*hit] || !is_within_patch(mesh, *hit, patch.first))
            continue;
          // We can't back trace on the doubled edges (not 2-manifold)
          if (!path.empty() && mesh.to_vertex(*hit) == path.back())
            continue;
          if (!is_step_valid(cur_v, mesh.to_vertex(*hit)))
            continue;
          next_v.emplace_back(mesh.to_vertex(*hit));
        } while (++hit != hit_end);

        if (next_v.empty()) {
          auto hit = mesh.halfedges(cur_v), hit_end = hit;
          // Move to the last step
          if (!path.empty()) {
            do {
              if (mesh.to_vertex(*hit) != path.back())
                continue;
              break;
            } while (++hit != hit_end);
            hit_end = hit;
          }
          logger().info("\tCir center: {}", cur_v);
          do {
            logger().info("\tTo: {}", mesh.to_vertex(*hit));
            if (!cut_directed[*hit] ||
                !is_within_patch(mesh, *hit, patch.first)) {
              logger().info("\t\tNot cut");
              continue;
            }
            // We can't back trace on the doubled edges (not 2-manifold)
            if (!path.empty() && mesh.to_vertex(*hit) == path.back()) {
              logger().info("\t\tBack tracing");
              continue;
            }
            if (!is_step_valid(cur_v, mesh.to_vertex(*hit))) {
              logger().info("\t\tStep not valid");
              continue;
            }
            next_v.emplace_back(mesh.to_vertex(*hit));
          } while (++hit != hit_end);

          logger().flush();
        }

        if (next_v.empty()) {
          logger().info("=====================");
          logger().info("Patch loop: {}", patch.first);
          for (size_t i = 0; i + 1 < path.size(); i++) {
            size_t j = i + 1;
            logger().info("\tHe: {} -> {}", path[i], path[j]);
          }
          return false;
        }

        contess_assert_msg(!next_v.empty(),
                           "Cut graph is not Eulerian. Stuck at v" +
                               std::to_string(cur_v.idx()));

        // Possible to have remaining degrees after this round of tracing
        if (next_v.size() > 1)
          candidate_start_v.emplace_back(cur_v);

        path.emplace_back(cur_v);

        // Always use the out path that is closest to the previous step wrt the
        // patch orientation
        if (patch_facing[patch.first] == FacingType::BACK) {
          cut_directed[mesh.find_halfedge(cur_v, next_v.front())] = false;
          cur_v = next_v.front();
        } else {
          cut_directed[mesh.find_halfedge(cur_v, next_v.back())] = false;
          cur_v = next_v.back();
        }
        remain_degree--;
      } while (cur_v != start_v);
      path.emplace_back(cur_v);

      // Always search for the vertex position along the path
      auto insert_pos =
          std::find(final_path.begin(), final_path.end(), path.front());
      contess_assert_msg(insert_pos != final_path.end(),
                         "Missing vertex from the path.");
      path.erase(path.begin()); // Avoid duplicate
      final_path.insert(++insert_pos, path.begin(), path.end());

      // Move to start from a candidate_start_v
      bool start_v_found = false;
      size_t candid_i = 0;
      for (auto const &v : candidate_start_v) {
        auto hit = mesh.halfedges(v), hit_end = hit;
        do {
          if (!cut_directed[*hit] || !is_within_patch(mesh, *hit, patch.first))
            continue;
          start_v_found = true;
        } while (++hit != hit_end);

        candid_i++;
        if (start_v_found) {
          start_v = v;
          break;
        }
      }
      if (!(start_v_found || remain_degree == 0)) {
        logger().error("=====================");
        logger().error("Unfinished patch loop: {}", patch.first);
        for (size_t i = 0; i + 1 < final_path.size(); i++) {
          size_t j = i + 1;
          logger().error("\tHe: {} -> {}; E: {}", final_path[i], final_path[j],
                         mesh.find_edge(final_path[i], final_path[j]));
        }
        return false;
      }
      contess_assert_msg(start_v_found || remain_degree == 0,
                         "Cut graph is not Eulerian. Remaining degree: " +
                             std::to_string(remain_degree));
      candidate_start_v.erase(candidate_start_v.begin(),
                              candidate_start_v.begin() + candid_i);
    }

    std::shared_ptr<Chain> chain = std::make_shared<Chain>();
    chain->reserve(final_path.size() - 1);
    for (size_t i = 0; i + 1 < final_path.size(); i++) {
      size_t j = i + 1;
      chain->emplace_back(mesh.find_halfedge(final_path[i], final_path[j]));
    }
    patch_chain.insert(std::make_pair(patch.first, chain));
  }

  // 4. Book keeping
  mesh.m_oriented_chains.clear();
  mesh.m_patch_chains.clear();

  for (auto const &patch : patch_chain) {
    mesh.m_patch_chains.insert(patch);
    mesh.m_oriented_chains.emplace_back(patch.second);
  }

  return true;
}
