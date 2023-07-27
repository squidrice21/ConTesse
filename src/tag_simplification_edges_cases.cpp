#include "tag_simplification_edges.h"

#include "add_simplification_vertex.h"
#include "tag_simplification_edges.h"

#include <spdlog/fmt/ostr.h>

#include <limits>
#include <unordered_set>

void get_crossing_halfedges(Mesh const &mesh, Vertex const &v,
                            std::vector<Halfedge> &halfedges) {
  auto intersection_2d = mesh.get_vertex_property<int>("v:intersection_2d");
  contess_assert_msg(
      intersection_2d,
      "tag_simplification_edges: 2D intersections are not created.");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto patchBoundary = mesh.get_halfedge_property<int>("h:patchBoundary");
  contess_assert_msg(patchBoundary,
                     "tag_2d_intersections: Patch needs to be built.");
  auto disk_cut = mesh.get_edge_property<int>("e:disk_cut");
  contess_assert_msg(
      disk_cut, "tag_simplification_edges: Patch needs to be cut to disk.");

  auto loop_info = mesh.get_halfedge_property<Vector2f>("h:loop");
  contess_assert_msg(
      loop_info,
      "tag_simplification_edges: Valid 2D intersections are not tagged.");

  auto is_valid_intersection_2d =
      mesh.get_vertex_property<bool>("v:is_valid_intersection_2d");
  contess_assert_msg(
      is_valid_intersection_2d,
      "tag_simplification_edges: Valid 2D intersections are not tagged.");

  Vertex inter_v(intersection_2d[v]);
  contess_assert_msg(
      inter_v.is_valid(),
      "tag_simplification_edges: Intersecting vertices mismatched.");

  auto get_loop_hes = [&](Vertex const &vv, std::vector<Halfedge> &hes) {
    auto hit = mesh.halfedges(vv), hit_end = hit;
    Halfedge out_he;
    do {
      if (loop_info[*hit].x() >= 0) {
        hes.emplace_back(*hit);
        out_he = *hit;
        break;
      }
    } while (++hit != hit_end);
    contess_assert_msg(
        out_he.is_valid(),
        "tag_simplification_edges: Halfedge information is missing.");

    hit = mesh.halfedges(vv), hit_end = hit;
    Halfedge in_he;
    do {
      auto oppo_he = mesh.opposite_halfedge(*hit);
      Edge oppo_e = mesh.edge(oppo_he);
      if ((mesh.is_boundary(oppo_e) || is_contour[oppo_e] >= 0) &&
          patchBoundary[oppo_he] == patchBoundary[out_he]) {
        hes.emplace_back(oppo_he);
        in_he = oppo_he;
        break;
      }
    } while (++hit != hit_end);

    contess_assert_msg(in_he.is_valid(),
                       "tag_simplification_edges: Intersecting edge doesn't "
                       "have matched passing halfedges.");
  };
  get_loop_hes(v, halfedges);
  get_loop_hes(inter_v, halfedges);
}

void is_halfedge_infeasible(Mesh const &mesh, Camera const &camera,
                            Vertex const &v, Vertex const &inter_v,
                            std::vector<Halfedge> const &halfedges,
                            std::vector<bool> &is_infeasible) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto patchBoundary = mesh.get_halfedge_property<int>("h:patchBoundary");

  FacingType facing = mesh.get_patch_facing(patchBoundary[halfedges.front()]);

  // Split based on the original intersecting edges
  std::vector<std::vector<Halfedge>> intersecting_halfedges;
  intersecting_halfedges.resize(2);
  intersecting_halfedges[0].insert(intersecting_halfedges[0].end(),
                                   halfedges.begin(), halfedges.begin() + 2);
  intersecting_halfedges[1].insert(intersecting_halfedges[1].end(),
                                   halfedges.begin() + 2, halfedges.end());

  auto get_nondegenerated_from_to = [&](Halfedge const &h, Vertex const &v,
                                        Vertex &from_v, Vertex &to_v) {
    bool to_flip = (mesh.to_vertex(h) == v);
    Halfedge hh = h;
    if (to_flip)
      hh = mesh.opposite_halfedge(hh);
    Halfedge h1_nondegenerated =
        find_nondegenerated_contour_halfedge(mesh, v, hh);
    if (!to_flip) {
      from_v = v;
      to_v = mesh.to_vertex(h1_nondegenerated);
    } else {
      from_v = mesh.to_vertex(h1_nondegenerated);
      to_v = v;
    }
  };

  // h2 on the correct side wrt h1?
  auto is_one_sided_correct_sided = [&](Halfedge const &h1,
                                        Halfedge const &h2) -> bool {
    Vertex h2_adj_v =
        (mesh.from_vertex(h2) == v || mesh.to_vertex(h2) == v) ? v : inter_v;
    Vertex h2_v_to, h2_v_from, h2_v;
    get_nondegenerated_from_to(h2, h2_adj_v, h2_v_from, h2_v_to);
    h2_v = (mesh.from_vertex(h2) == v || mesh.from_vertex(h2) == inter_v)
               ? h2_v_to
               : h2_v_from;
    Vertex h1_v_to, h1_v_from;
    Vertex h1_adj_v =
        (mesh.from_vertex(h1) == v || mesh.to_vertex(h1) == v) ? v : inter_v;
    get_nondegenerated_from_to(h1, h1_adj_v, h1_v_from, h1_v_to);

    std::vector<Vertex> vertices({h1_v_from, h1_v_to, h2_v});
    Vector2f v0_2d =
        project(vpositions[vertices[0]], camera.viewMatrix().matrix(),
                camera.projectionMatrix(),
                Vector2i(camera.vpWidth(), camera.vpHeight()))
            .head(2);
    Vector2f v1_2d =
        project(vpositions[vertices[1]], camera.viewMatrix().matrix(),
                camera.projectionMatrix(),
                Vector2i(camera.vpWidth(), camera.vpHeight()))
            .head(2);
    Vector2f v2_2d =
        project(vpositions[vertices[2]], camera.viewMatrix().matrix(),
                camera.projectionMatrix(),
                Vector2i(camera.vpWidth(), camera.vpHeight()))
            .head(2);

    auto f_ori = igl::predicates::orient2d(v0_2d, v1_2d, v2_2d);

    if (f_ori == igl::predicates::Orientation::COLLINEAR) {
      logger().error("Generic camera assumption violated: {}, {}, {}; "
                     "Intersecting: {}, {}",
                     vertices[0], vertices[1], vertices[2], v, inter_v);
    }

    contess_assert_msg(f_ori != igl::predicates::Orientation::COLLINEAR,
                       "Generic camera assumption violated.");

    // FF face has CCW orientation; BF has CW
    return (facing == FacingType::FRONT &&
            f_ori == igl::predicates::Orientation::POSITIVE) ||
           (facing == FacingType::BACK &&
            f_ori == igl::predicates::Orientation::NEGATIVE);
  };
  auto is_mutually_correct_sided = [&](Halfedge const &h1,
                                       Halfedge const &h2) -> bool {
    return is_one_sided_correct_sided(h1, h2) &&
           is_one_sided_correct_sided(h2, h1);
  };

  std::vector<std::vector<bool>> intersecting_is_infeasible;
  intersecting_is_infeasible.resize(2);
  intersecting_is_infeasible[0].resize(2, true);
  intersecting_is_infeasible[1].resize(2, true);

  for (size_t i = 0; i < intersecting_halfedges[0].size(); i++) {
    for (size_t j = 0; j < intersecting_halfedges[1].size(); j++) {
      if (is_mutually_correct_sided(intersecting_halfedges[0][i],
                                    intersecting_halfedges[1][j])) {
        intersecting_is_infeasible[0][i] = false;
        intersecting_is_infeasible[1][j] = false;
      }
    }
  }

  is_infeasible.reserve(halfedges.size());
  is_infeasible.insert(is_infeasible.end(),
                       intersecting_is_infeasible[0].begin(),
                       intersecting_is_infeasible[0].end());
  is_infeasible.insert(is_infeasible.end(),
                       intersecting_is_infeasible[1].begin(),
                       intersecting_is_infeasible[1].end());
}

bool test_walk_length(Mesh &mesh, Camera const &camera, Vertex const &v,
                      std::vector<Halfedge> const &halfedges,
                      std::vector<bool> const &is_infeasible,
                      real_t &chain_walk_length, real_t &walk_length) {
  auto intersection_2d = mesh.get_vertex_property<int>("v:intersection_2d");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_valid_intersection_2d =
      mesh.get_vertex_property<bool>("v:is_valid_intersection_2d");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto patchBoundary = mesh.get_halfedge_property<int>("h:patchBoundary");
  contess_assert_msg(patchBoundary,
                     "tag_2d_intersections: Patch needs to be built.");
  auto disk_cut = mesh.get_edge_property<int>("e:disk_cut");
  contess_assert_msg(
      disk_cut, "tag_simplification_edges: Patch needs to be cut to disk.");

  walk_length = 0;
  chain_walk_length = 0;

  // Find the single tracing direction determined by the crossing halfedges and
  // feasibility flags
  Halfedge pair_he;
  Halfedge in_he;
  for (size_t i = 0; i < halfedges.size(); i++) {
    if (mesh.from_vertex(halfedges[i]) != v &&
        mesh.to_vertex(halfedges[i]) != v)
      pair_he = halfedges[i];

    // Avoid walking on halfedges not adjacent to the current vertex
    if (!is_infeasible[i] || (mesh.from_vertex(halfedges[i]) != v &&
                              mesh.to_vertex(halfedges[i]) != v)) {
      continue;
    }
    in_he = halfedges[i];
  }

  if (!pair_he.is_valid() || !in_he.is_valid())
    return false;

  Halfedge he =
      (mesh.from_vertex(in_he) == v) ? in_he : mesh.opposite_halfedge(in_he);
  in_he = he;

  // For safety reason, first compute the length of the collapsing edges
  bool sum_walk = true;
  Vertex ret_v;
  std::deque<int> ret_v_queue;
  std::unordered_set<size_t> seen_vertices;
  while (mesh.to_vertex(he) != v) {
    if (mesh.from_vertex(he) != v &&
        is_valid_intersection_2d[mesh.from_vertex(he)]) {
      sum_walk = false;
      if (is_valid_intersection_2d[mesh.from_vertex(he)] &&
          (ret_v_queue.empty() ||
           mesh.from_vertex(he).idx() != ret_v_queue.front()))
        ret_v_queue.push_front(intersection_2d[mesh.from_vertex(he)]);
    }

    if (!ret_v_queue.empty() &&
        mesh.from_vertex(he).idx() == ret_v_queue.front()) {
      ret_v_queue.pop_front();
      if (ret_v_queue.empty())
        sum_walk = true;
    }

    Vector2f v0_2d =
        project(vpositions[mesh.from_vertex(he)], camera.viewMatrix().matrix(),
                camera.projectionMatrix(),
                Vector2i(camera.vpWidth(), camera.vpHeight()))
            .head(2);
    Vector2f v1_2d =
        project(vpositions[mesh.to_vertex(he)], camera.viewMatrix().matrix(),
                camera.projectionMatrix(),
                Vector2i(camera.vpWidth(), camera.vpHeight()))
            .head(2);
    if (sum_walk) {
      walk_length += (v0_2d - v1_2d).norm();
    }

    // We looped over a traversed vertex
    if (seen_vertices.count(mesh.to_vertex(he).idx()))
      return false;
    seen_vertices.emplace(mesh.to_vertex(he).idx());

    chain_walk_length += (v0_2d - v1_2d).norm();

    auto hit = mesh.halfedges(mesh.to_vertex(he)), hit_end = hit;
    Halfedge next_he;
    do {
      if ((mesh.is_boundary(mesh.edge(*hit)) ||
           is_contour[mesh.edge(*hit)] >= 0) &&
          is_same_patch(mesh, in_he, pair_he, *hit)) {
        next_he = *hit;
        break;
      }
    } while (++hit != hit_end);

    if (!next_he.is_valid()) {
      auto hit = mesh.halfedges(mesh.to_vertex(he)), hit_end = hit;
      Halfedge next_he;
      do {
        if ((mesh.is_boundary(mesh.edge(*hit)) ||
             is_contour[mesh.edge(*hit)] >= 0) &&
            is_same_patch(mesh, in_he, pair_he, *hit)) {
          next_he = *hit;
          break;
        }
      } while (++hit != hit_end);
    }

    contess_assert_msg(
        next_he.is_valid(),
        "tag_simplification_edges: Cannot propagate collapsing label.");
    he = next_he;
  }

  return true;
}

void walk_furthest_cusp(Mesh &mesh, Camera const &camera, Vertex const &v,
                        std::vector<Halfedge> const &halfedges,
                        std::vector<bool> const &is_infeasible, Vertex &cusp_v,
                        real_t &furthest_walk_length,
                        bool write_feasibility_collapsing = true) {
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto is_valid_intersection_2d =
      mesh.get_vertex_property<bool>("v:is_valid_intersection_2d");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto patchBoundary = mesh.get_halfedge_property<int>("h:patchBoundary");
  contess_assert_msg(patchBoundary,
                     "tag_2d_intersections: Patch needs to be built.");
  auto disk_cut = mesh.get_edge_property<int>("e:disk_cut");
  contess_assert_msg(
      disk_cut, "tag_simplification_edges: Patch needs to be cut to disk.");

  auto feasibility_collapsing =
      mesh.edge_property<int>("e:feasibility_collapsing");

  auto dir_2d = [&](Vector3f const &v1, Vector3f const &v2, Vector2f &dir) {
    Vector2f v0_2d =
        project(v1, camera.viewMatrix().matrix(), camera.projectionMatrix(),
                Vector2i(camera.vpWidth(), camera.vpHeight()))
            .head(2);
    Vector2f v1_2d =
        project(v2, camera.viewMatrix().matrix(), camera.projectionMatrix(),
                Vector2i(camera.vpWidth(), camera.vpHeight()))
            .head(2);
    dir = (v1_2d - v0_2d).normalized();
  };

  // Compute chain length
  real_t chain_length = 0;
  for (size_t i = 0; i < halfedges.size(); i++) {
    // Avoid walking on halfedges not adjacent to the current vertex
    if (!is_infeasible[i] || (mesh.from_vertex(halfedges[i]) != v &&
                              mesh.to_vertex(halfedges[i]) != v))
      continue;

    Halfedge he = (mesh.from_vertex(halfedges[i]) == v)
                      ? halfedges[i]
                      : mesh.opposite_halfedge(halfedges[i]);

    // For safety reason, first compute the length of the collapsing edges
    Vector2f prev_dir;
    Vector2f cur_dir;
    while (mesh.to_vertex(he) != v) {
      if (is_valid_intersection_2d[mesh.to_vertex(he)] ||
          is_boundary_joint(mesh, mesh.to_vertex(he)) ||
          (feasibility_collapsing[mesh.edge(he)] >= 0))
        break;

      Vector2f v0_2d =
          project(vpositions[mesh.from_vertex(he)],
                  camera.viewMatrix().matrix(), camera.projectionMatrix(),
                  Vector2i(camera.vpWidth(), camera.vpHeight()))
              .head(2);
      Vector2f v1_2d =
          project(vpositions[mesh.to_vertex(he)], camera.viewMatrix().matrix(),
                  camera.projectionMatrix(),
                  Vector2i(camera.vpWidth(), camera.vpHeight()))
              .head(2);
      if (mesh.from_vertex(he) == v) {
        auto he_forward = find_nondegenerated_contour_halfedge(mesh, v, he);
        dir_2d(vpositions[mesh.from_vertex(he)],
               vpositions[mesh.to_vertex(he_forward)], prev_dir);
      } else {
        prev_dir = cur_dir;
      }
      auto he_forward =
          find_nondegenerated_contour_halfedge(mesh, mesh.from_vertex(he), he);
      dir_2d(vpositions[mesh.from_vertex(he)],
             vpositions[mesh.to_vertex(he_forward)], cur_dir);
      chain_length += (v0_2d - v1_2d).norm();

      auto hit = mesh.halfedges(mesh.to_vertex(he)), hit_end = hit;
      Halfedge next_he;
      do {
        if ((mesh.is_boundary(mesh.edge(*hit)) ||
             is_contour[mesh.edge(*hit)] >= 0) &&
            patchBoundary[*hit] == patchBoundary[he]) {
          next_he = *hit;
          break;
        }
      } while (++hit != hit_end);

      contess_assert_msg(
          next_he.is_valid(),
          "tag_simplification_edges: Cannot propagate collapsing label.");
      he = next_he;
    }

    break;
  }

  furthest_walk_length = -1;
  real_t walk_length = 0;
  Vertex in_cusp_v = cusp_v;

  real_t chain_walk_length = 0;
  bool is_local_maxima = false;
  for (size_t i = 0; i < halfedges.size(); i++) {
    // Avoid walking on halfedges not adjacent to the current vertex
    if (!is_infeasible[i] || (mesh.from_vertex(halfedges[i]) != v &&
                              mesh.to_vertex(halfedges[i]) != v))
      continue;

    Halfedge he = (mesh.from_vertex(halfedges[i]) == v)
                      ? halfedges[i]
                      : mesh.opposite_halfedge(halfedges[i]);

    // For safety reason, first compute the length of the collapsing edges
    Vector2f prev_dir;
    Vector2f cur_dir;
    real_t walk_sign = 1;
    while (mesh.to_vertex(he) != v) {
      if (is_valid_intersection_2d[mesh.to_vertex(he)] ||
          is_boundary_joint(mesh, mesh.to_vertex(he)) ||
          (feasibility_collapsing[mesh.edge(he)] >= 0))
        break;

      if (!in_cusp_v.is_valid()) {
        Vector2f v0_2d =
            project(vpositions[mesh.from_vertex(he)],
                    camera.viewMatrix().matrix(), camera.projectionMatrix(),
                    Vector2i(camera.vpWidth(), camera.vpHeight()))
                .head(2);
        Vector2f v1_2d =
            project(vpositions[mesh.to_vertex(he)],
                    camera.viewMatrix().matrix(), camera.projectionMatrix(),
                    Vector2i(camera.vpWidth(), camera.vpHeight()))
                .head(2);
        if (mesh.from_vertex(he) == v) {
          auto he_forward = find_nondegenerated_contour_halfedge(mesh, v, he);
          dir_2d(vpositions[mesh.from_vertex(he)],
                 vpositions[mesh.to_vertex(he_forward)], prev_dir);
        } else {
          prev_dir = cur_dir;
        }
        auto he_forward = find_nondegenerated_contour_halfedge(
            mesh, mesh.from_vertex(he), he);
        dir_2d(vpositions[mesh.from_vertex(he)],
               vpositions[mesh.to_vertex(he_forward)], cur_dir);
        walk_sign *= (prev_dir.dot(cur_dir) > 0) ? 1 : -1;
        walk_length += walk_sign * (v0_2d - v1_2d).norm();
        chain_walk_length += (v0_2d - v1_2d).norm();
        if (walk_length > furthest_walk_length) {
          furthest_walk_length = walk_length;
          cusp_v = mesh.to_vertex(he);
          is_local_maxima = false;
        } else {
          is_local_maxima = true;
        }

        if (chain_walk_length > 200 || chain_walk_length / chain_length > 0.1) {
          if (furthest_walk_length < 0 || !is_local_maxima)
            furthest_walk_length = std::numeric_limits<real_t>::infinity();
          return;
        }
      } else {
        if (mesh.from_vertex(he) == in_cusp_v) {
          if (furthest_walk_length < 0 || !is_local_maxima)
            furthest_walk_length = std::numeric_limits<real_t>::infinity();
          return;
        }
        if (write_feasibility_collapsing)
          feasibility_collapsing[mesh.edge(he)] = v.idx();
      }

      auto hit = mesh.halfedges(mesh.to_vertex(he)), hit_end = hit;
      Halfedge next_he;
      do {
        if ((mesh.is_boundary(mesh.edge(*hit)) ||
             is_contour[mesh.edge(*hit)] >= 0) &&
            patchBoundary[*hit] == patchBoundary[he]) {
          next_he = *hit;
          break;
        }
      } while (++hit != hit_end);

      contess_assert_msg(
          next_he.is_valid(),
          "tag_simplification_edges: Cannot propagate collapsing label.");
      he = next_he;
    }

    break;
  }

  if (furthest_walk_length < 0 || !is_local_maxima)
    furthest_walk_length = std::numeric_limits<real_t>::infinity();
}

bool label_edges(Mesh &mesh, Vertex const &v,
                 std::vector<Halfedge> const &halfedges,
                 std::vector<bool> const &is_infeasible) {
  auto intersection_2d = mesh.get_vertex_property<int>("v:intersection_2d");
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto is_valid_intersection_2d =
      mesh.get_vertex_property<bool>("v:is_valid_intersection_2d");
  auto patchBoundary = mesh.get_halfedge_property<int>("h:patchBoundary");
  contess_assert_msg(patchBoundary,
                     "tag_2d_intersections: Patch needs to be built.");
  auto feasibility_collapsing =
      mesh.edge_property<int>("e:feasibility_collapsing");

  bool is_back = false;
  for (size_t i = 0; i < halfedges.size(); i++) {
    // Avoid walking on halfedges not adjacent to the current vertex
    if (!is_infeasible[i] || (mesh.from_vertex(halfedges[i]) != v &&
                              mesh.to_vertex(halfedges[i]) != v))
      continue;

    Halfedge he = (mesh.from_vertex(halfedges[i]) == v)
                      ? halfedges[i]
                      : mesh.opposite_halfedge(halfedges[i]);

    // For safety reason, first compute the length of the collapsing edges
    while (feasibility_collapsing[mesh.edge(he)] < 0) {
      feasibility_collapsing[mesh.edge(he)] = v.idx();

      if (is_valid_intersection_2d[mesh.to_vertex(he)] ||
          is_boundary_joint(mesh, mesh.to_vertex(he))) {
        if (intersection_2d[mesh.to_vertex(he)] == v.idx()) {
          is_back = true;
        }
        break;
      }

      auto hit = mesh.halfedges(mesh.to_vertex(he)), hit_end = hit;
      Halfedge next_he;
      do {
        if ((mesh.is_boundary(mesh.edge(*hit)) ||
             is_contour[mesh.edge(*hit)] >= 0) &&
            patchBoundary[*hit] == patchBoundary[he]) {
          next_he = *hit;
          break;
        }
      } while (++hit != hit_end);

      contess_assert_msg(
          next_he.is_valid(),
          "tag_simplification_edges: Cannot propagate collapsing label.");
      he = next_he;
    }
  }

  return is_back;
}

/* ===================================================================== */

bool tag_simplification_edges_case1(
    Mesh const &mesh, Camera const &camera, Vertex v, Mesh &mesh_labeled,
    real_t &walk_length1, real_t &chain_walk_length1, real_t &walk_length2,
    real_t &chain_walk_length2, size_t &case1_count, size_t &case2_count) {
  mesh_labeled = mesh;
  mesh_labeled.m_fedges.clear();
  mesh_labeled.m_svertices.clear();
  case1_count = case2_count = 0;

  auto intersection_2d =
      mesh_labeled.get_vertex_property<int>("v:intersection_2d");
  contess_assert_msg(
      intersection_2d,
      "tag_simplification_edges: 2D intersections are not created.");
  Vertex inter_v(intersection_2d[v]);

  // Determine the direction to walk
  // Get the crossing four halfedges
  std::vector<Halfedge> halfedges;
  get_crossing_halfedges(mesh_labeled, v, halfedges);

  std::vector<bool> is_infeasible;
  is_halfedge_infeasible(mesh_labeled, camera, v, inter_v, halfedges,
                         is_infeasible);

  // Test run
  bool successful =
      test_walk_length(mesh_labeled, camera, v, halfedges, is_infeasible,
                       chain_walk_length1, walk_length1);
  successful &=
      test_walk_length(mesh_labeled, camera, inter_v, halfedges, is_infeasible,
                       chain_walk_length2, walk_length2);
  if (!successful)
    return false;

  bool is_back = label_edges(mesh_labeled, v, halfedges, is_infeasible);
  label_edges(mesh_labeled, inter_v, halfedges, is_infeasible);

  if (is_back)
    case1_count++;
  else
    case2_count++;

  return true;
}

bool tag_simplification_edges_case1_alt(
    Mesh const &mesh, Camera const &camera, Vertex v, Mesh &mesh_labeled,
    real_t &walk_length1, real_t &chain_walk_length1, real_t &walk_length2,
    real_t &chain_walk_length2, size_t &case1_count, size_t &case2_count) {
  mesh_labeled = mesh;
  mesh_labeled.m_fedges.clear();
  mesh_labeled.m_svertices.clear();
  case1_count = case2_count = 0;

  auto intersection_2d =
      mesh_labeled.get_vertex_property<int>("v:intersection_2d");
  contess_assert_msg(
      intersection_2d,
      "tag_simplification_edges: 2D intersections are not created.");
  Vertex inter_v(intersection_2d[v]);

  // Determine the direction to walk
  // Get the crossing four halfedges
  std::vector<Halfedge> halfedges;
  get_crossing_halfedges(mesh_labeled, v, halfedges);

  std::vector<bool> is_infeasible;
  is_halfedge_infeasible(mesh_labeled, camera, v, inter_v, halfedges,
                         is_infeasible);

  std::vector<bool> is_infeasible2;
  is_infeasible2.resize(is_infeasible.size());
  for (size_t j = 0; j < is_infeasible.size(); j++) {
    is_infeasible2[j] = !is_infeasible[j];
  }

  bool successful =
      test_walk_length(mesh_labeled, camera, v, halfedges, is_infeasible2,
                       chain_walk_length1, walk_length1);
  successful &=
      test_walk_length(mesh_labeled, camera, v, halfedges, is_infeasible2,
                       chain_walk_length2, walk_length2);
  if (!successful)
    return false;

  bool is_back = label_edges(mesh_labeled, v, halfedges, is_infeasible2);
  label_edges(mesh_labeled, inter_v, halfedges, is_infeasible2);

  if (is_back)
    case1_count++;
  else
    case2_count++;

  return true;
}

bool tag_simplification_edges_case3(
    Mesh const &mesh, Camera const &camera, Vertex v, Mesh &mesh_labeled,
    std::unordered_map<int, int> &cusp_projections,
    std::unordered_map<int, Vector3f> &moved_proj_positions,
    real_t &furthest_walk_length_min, real_t &distance) {
  mesh_labeled = mesh;
  mesh_labeled.m_fedges.clear();
  mesh_labeled.m_svertices.clear();
  distance = std::numeric_limits<real_t>::infinity();

  auto intersection_2d =
      mesh_labeled.get_vertex_property<int>("v:intersection_2d");
  contess_assert_msg(
      intersection_2d,
      "tag_simplification_edges: 2D intersections are not created.");
  auto is_contour = mesh_labeled.get_edge_property<real_t>("e:contour");
  auto patchBoundary =
      mesh_labeled.get_halfedge_property<int>("h:patchBoundary");
  contess_assert_msg(patchBoundary,
                     "tag_2d_intersections: Patch needs to be built.");
  auto cusp_facing =
      mesh_labeled.get_vertex_property<FacingType>("v:cusp_facing");
  contess_assert_msg(is_contour && cusp_facing, "Requires contours and cusps.");

  Vertex inter_v(intersection_2d[v]);

  // Determine the direction to walk
  // Get the crossing four halfedges
  std::vector<Halfedge> halfedges;
  get_crossing_halfedges(mesh_labeled, v, halfedges);

  std::vector<bool> is_infeasible;
  is_halfedge_infeasible(mesh_labeled, camera, v, inter_v, halfedges,
                         is_infeasible);

  // Try the flipped cusp edge solution
  Vertex cusp_v1, cusp_v2;
  real_t furthest_walk_length1, furthest_walk_length2;
  walk_furthest_cusp(mesh_labeled, camera, v, halfedges, is_infeasible, cusp_v1,
                     furthest_walk_length1);
  walk_furthest_cusp(mesh_labeled, camera, inter_v, halfedges, is_infeasible,
                     cusp_v2, furthest_walk_length2);

  // Compare to see which one is more possible
  // TODO: Better comparison: width
  furthest_walk_length_min =
      std::min(furthest_walk_length1, furthest_walk_length2);

  if (furthest_walk_length_min >= 40)
    return false;

  // Invalid cusp flipping case
  if (!cusp_v1.is_valid() && !cusp_v2.is_valid()) {
    logger().error(
        "tag_simplification_edges: Cannot find cusps from {} - {}; {} - {}", v,
        cusp_v1, inter_v, cusp_v2);
    logger().error("\twalk_length: {}, {}", furthest_walk_length1,
                   furthest_walk_length2);
    return false;
  }

  // Actual flipping labeling
  // Make the unfound cusp side infinitely far
  if (furthest_walk_length1 < 0)
    furthest_walk_length1 = std::numeric_limits<real_t>::infinity();
  if (furthest_walk_length2 < 0)
    furthest_walk_length2 = std::numeric_limits<real_t>::infinity();

  Vertex cusp_v, start_v;
  cusp_v = (furthest_walk_length1 < furthest_walk_length2) ? cusp_v1 : cusp_v2;
  start_v = (furthest_walk_length1 < furthest_walk_length2) ? v : inter_v;

  int pidx = patchBoundary[halfedges[0]];
  FacingType patch_cusp_facing = mesh.get_patch_facing(pidx);
  if (cusp_facing[cusp_v] != FacingType::UNDEFINED &&
      cusp_facing[cusp_v] != patch_cusp_facing) {
    real_t furthest_walk_length;
    walk_furthest_cusp(mesh_labeled, camera, start_v, halfedges, is_infeasible,
                       cusp_v, furthest_walk_length);
  }

  contess_assert_msg(cusp_v.is_valid(),
                     "tag_simplification_edges: Cannot find the cusp to flip.");

  Vertex other_v = (start_v == v) ? inter_v : v;
  Halfedge he_line, he_cusp;
  for (size_t k = 0; k < halfedges.size(); k++) {
    // Avoid walking on halfedges not adjacent to the current vertex
    if (!is_infeasible[k])
      continue;
    if (!(mesh_labeled.from_vertex(halfedges[k]) != other_v &&
          mesh_labeled.to_vertex(halfedges[k]) != other_v)) {
      he_line = halfedges[k];
      if (mesh_labeled.from_vertex(halfedges[k]) != other_v)
        he_line = mesh_labeled.opposite_halfedge(he_line);
    }
    if (!(mesh_labeled.from_vertex(halfedges[k]) != start_v &&
          mesh_labeled.to_vertex(halfedges[k]) != start_v)) {
      he_cusp = halfedges[k];
      if (mesh_labeled.from_vertex(halfedges[k]) != start_v)
        he_cusp = mesh_labeled.opposite_halfedge(he_cusp);
    }
  }

  Vertex added_v;
  if (!cusp_projections.count(cusp_v.idx())) {
    Vector3f moved_pos;
    added_v = add_simplification_vertex(mesh_labeled, camera, cusp_v, he_cusp,
                                        he_line, moved_pos);
    if (!added_v.is_valid())
      return false;
    cusp_projections[cusp_v.idx()] = added_v.idx();
    moved_proj_positions[added_v.idx()] = moved_pos;
  } else {
    added_v = Vertex(cusp_projections[cusp_v.idx()]);
  }

  // Tag here
  if (added_v.is_valid()) {
    // In case the new vertex is added immediately adjacent to the start
    // vertex (the old he invalid)
    auto new_he = mesh_labeled.find_halfedge(other_v, added_v);
    if (new_he.is_valid() &&
        (is_contour[mesh_labeled.edge(new_he)] >= 0 ||
         mesh_labeled.is_boundary(mesh_labeled.edge(new_he))))
      he_line = new_he;

    distance = tag_cusp_line_case(mesh_labeled, camera, cusp_v, added_v,
                                  he_cusp, he_line, false);

    tag_cusp_line_case(mesh_labeled, camera, cusp_v, added_v, he_cusp, he_line);
  } else {
    logger().warn("tag_simplification_edges: Failed to add cusp "
                  "projection at {} on side containing {}.",
                  cusp_v, mesh_labeled.from_vertex(he_line));
    return false;
  }

  return true;
}
