// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "tag_non_fish_tail.h"
#include "common.h"

#include <igl/predicates/predicates.h>
#include <limits>
#include <spdlog/fmt/ostr.h>
#include <unordered_map>
#include <vector>

real_t max_non_fishtail_length = 2;
real_t max_non_fishtail_length_large = 10;

// When walking CCW. 0: Turn right (correct); 1: Turn left (flipped).
size_t determine_turn_side(Mesh const &mesh, Camera const &camera,
                           std::shared_ptr<Chain> chain, size_t he_outward) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  size_t side = 0;

  // Check the side of the 2D projection
  Vector2f c1_2d =
      project(vpositions[mesh.from_vertex(chain->at(he_outward))],
              camera.viewMatrix().matrix(), camera.projectionMatrix(),
              Vector2i(camera.vpWidth(), camera.vpHeight()))
          .head(2);
  auto prev_start = chain->at((he_outward + chain->size() - 1) % chain->size());
  Vector2f prev_c1_2d =
      project(vpositions[mesh.from_vertex(prev_start)],
              camera.viewMatrix().matrix(), camera.projectionMatrix(),
              Vector2i(camera.vpWidth(), camera.vpHeight()))
          .head(2);
  Vector2f next_c1_2d =
      project(vpositions[mesh.to_vertex(chain->at(he_outward))],
              camera.viewMatrix().matrix(), camera.projectionMatrix(),
              Vector2i(camera.vpWidth(), camera.vpHeight()))
          .head(2);
  // Turns right
  if (igl::predicates::orient2d(c1_2d, next_c1_2d, prev_c1_2d) ==
      igl::predicates::Orientation::NEGATIVE) {
    side = 1;
  }

  return side;
}

// Note: this function needs to be called after removing small flipped loops.
// Otherwise, it may generate false positives.
void tag_non_fish_tail(Mesh &mesh, Camera const &camera, bool ff_only) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto cusp_facing = mesh.vertex_property<FacingType>("v:cusp_facing");
  auto is_cusp = mesh.vertex_property<bool>("v:cusp");
  contess_assert_msg(cusp_facing, "Cusp facing is missing.");

  auto is_valid_intersection_2d =
      mesh.get_vertex_property<bool>("v:is_valid_intersection_2d");
  auto intersection_2d = mesh.get_vertex_property<int>("v:intersection_2d");

  bool to_delete_chain = false;
  if (mesh.get_oriented_chains().empty()) {
    chain_contour(mesh, camera, true, ff_only);
    to_delete_chain = true;
  }
  contess_assert_msg(!mesh.get_oriented_chains().empty(),
                     "Contour chain is missing.");

  auto project2d = [&](Vector3f const &v3d, Vector2f &pos2D) {
    Vector2i viewport({camera.vpWidth(), camera.vpHeight()});
    pos2D = project(v3d, camera.viewMatrix().matrix(),
                    camera.projectionMatrix(), viewport)
                .head<2>();
  };
  for (auto chain : mesh.get_oriented_chains()) {
    // 1. Move to a cusp as the starting point
    size_t start_i = 0;
    for (; start_i < chain->size() &&
           !is_cusp[mesh.from_vertex(chain->at(start_i))];
         start_i++) {
    }

    // No cusp
    if (start_i == chain->size())
      continue;

    // 2. For each convex/concave session, check if it forms a fish tail.
    struct NonFishtailSegment {
      Vertex v1, v2;
      double length;
      size_t valid_2d_intersection_count;
      size_t seen_2d_intersection_count;
    };
    std::vector<NonFishtailSegment> non_fishtail_v;
    size_t moving_start = start_i;
    do {
      // Find the end
      size_t valid_2d_intersection_count = 0;
      size_t seen_2d_intersection_count = 0;
      double segment_length = 0;
      size_t moving_end = moving_start;
      do {
        if (is_cusp[mesh.from_vertex(chain->at(moving_end))] &&
            moving_end != moving_start)
          break;
        Vector2f seg1, seg2;
        project2d(vpositions[mesh.from_vertex(chain->at(moving_end))], seg1);
        project2d(vpositions[mesh.to_vertex(chain->at(moving_end))], seg2);
        segment_length += (seg1 - seg2).norm();

        if (is_valid_intersection_2d &&
            is_valid_intersection_2d[mesh.from_vertex(chain->at(moving_end))])
          valid_2d_intersection_count++;
        if (intersection_2d &&
            intersection_2d[mesh.from_vertex(chain->at(moving_end))] >= 0)
          seen_2d_intersection_count++;

        moving_end = (moving_end + 1) % chain->size();
      } while (moving_end != moving_start);

      if (moving_end == moving_start) {
        logger().warn("tag_non_fish_tail: Saw a single cusp at {}.",
                      mesh.from_vertex(chain->at(moving_start)));
        break;
      }

      contess_assert_msg(moving_end != moving_start,
                         "There must be even number of cusps.");

      // Check the side of the 2D projection
      Vector2f c1_2d =
          project(vpositions[mesh.from_vertex(chain->at(moving_start))],
                  camera.viewMatrix().matrix(), camera.projectionMatrix(),
                  Vector2i(camera.vpWidth(), camera.vpHeight()))
              .head(2);
      Vector2f c2_2d =
          project(vpositions[mesh.from_vertex(chain->at(moving_end))],
                  camera.viewMatrix().matrix(), camera.projectionMatrix(),
                  Vector2i(camera.vpWidth(), camera.vpHeight()))
              .head(2);
      auto prev_start =
          chain->at((moving_start + chain->size() - 1) % chain->size());
      Vector2f prev_c1_2d =
          project(vpositions[mesh.from_vertex(prev_start)],
                  camera.viewMatrix().matrix(), camera.projectionMatrix(),
                  Vector2i(camera.vpWidth(), camera.vpHeight()))
              .head(2);
      Vector2f next_c1_2d =
          project(vpositions[mesh.to_vertex(chain->at(moving_start))],
                  camera.viewMatrix().matrix(), camera.projectionMatrix(),
                  Vector2i(camera.vpWidth(), camera.vpHeight()))
              .head(2);
      auto prev_end =
          chain->at((moving_end + chain->size() - 1) % chain->size());
      Vector2f prev_c2_2d =
          project(vpositions[mesh.from_vertex(prev_end)],
                  camera.viewMatrix().matrix(), camera.projectionMatrix(),
                  Vector2i(camera.vpWidth(), camera.vpHeight()))
              .head(2);
      Vector2f next_c2_2d =
          project(vpositions[mesh.to_vertex(chain->at(moving_end))],
                  camera.viewMatrix().matrix(), camera.projectionMatrix(),
                  Vector2i(camera.vpWidth(), camera.vpHeight()))
              .head(2);
      // Not a fish tail
      if (igl::predicates::orient2d(c1_2d, next_c1_2d, prev_c1_2d) !=
          igl::predicates::orient2d(c2_2d, next_c2_2d, prev_c2_2d)) {
        NonFishtailSegment seg;
        seg.v1 = mesh.from_vertex(chain->at(moving_start));
        seg.v2 = mesh.from_vertex(chain->at(moving_end));
        seg.length = segment_length;
        seg.valid_2d_intersection_count = valid_2d_intersection_count;
        seg.seen_2d_intersection_count = seen_2d_intersection_count;

        non_fishtail_v.emplace_back(seg);
      }

      // Move to start from the end
      moving_start = moving_end;
    } while (moving_start != start_i);

    // 3. Find the segment to delete
    std::unordered_map<int, size_t> cusp_count;
    for (auto const &seg : non_fishtail_v) {
      cusp_count[seg.v1.idx()]++;
      cusp_count[seg.v2.idx()]++;
    }

    std::vector<Vertex> to_delete_vertices;
    for (auto const &v_c : cusp_count) {
      if (v_c.second != 2)
        continue;

      Vertex v(v_c.first);
      // Find the pair connected with the shorter side
      Vertex v_pair;
      double min_length = std::numeric_limits<double>::infinity();
      NonFishtailSegment near_seg;
      for (auto const &seg : non_fishtail_v) {
        if (seg.v1 != v && seg.v2 != v)
          continue;
        if (seg.length < min_length) {
          min_length = seg.length;
          v_pair = (seg.v1 != v) ? seg.v1 : seg.v2;
          near_seg = seg;
        }
      }

      contess_assert(v_pair.is_valid());
      if ((near_seg.length < max_non_fishtail_length &&
           (near_seg.valid_2d_intersection_count % 2 == 0)) ||
          (near_seg.seen_2d_intersection_count == 0 &&
           near_seg.length < max_non_fishtail_length_large)) {
        to_delete_vertices.emplace_back(v);
        to_delete_vertices.emplace_back(v_pair);
      }
    }

    // 4. Book keeping
    for (auto const &v : to_delete_vertices) {
      logger().info("tag_non_fish_tail: Remove non fishtail cusp {}", v);
      cusp_facing[v] = FacingType::UNDEFINED;
      is_cusp[v] = false;
    }
  }

  if (to_delete_chain) {
    mesh.get_patch_chains().clear();
    mesh.get_oriented_chains().clear();
  }
}
