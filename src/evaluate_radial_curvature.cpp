// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "evaluate_radial_curvature.h"

#include "common.h"
#include "subdiv_common.h"
// #include "subdiv_osd.h"
#include "tag_cusp_facing.h"
#include <Eigen/src/Core/Matrix.h>
#include <spdlog/fmt/ostr.h>
#include <unordered_set>
#include <utility>
#include <vector>

real_t evaluate_radial_curvature_analytic(Mesh &mesh, Subdiv const &subdiv,
                                          Camera const &camera, Vertex const &v,
                                          Param_loc const &v_param_loc) {
  auto d2x_dw2_v = mesh.get_vertex_property<Vector3f>("v:d2x_dw2");
  if (v.is_valid()) {
    if (!d2x_dw2_v) {
      d2x_dw2_v =
          mesh.add_vertex_property<Vector3f>("v:d2x_dw2", Vector3f::Zero());
    }
    d2x_dw2_v = mesh.vertex_property<Vector3f>("v:d2x_dw2");
  }

  Vector3f pos, norm;
  // Find the viewDir in UV coordinate
  Vector3f ds, dt, dsds, dsdt, dtdt;
  subdiv.evaluateLimitFrame(v_param_loc, pos, ds, dt, &dsds, &dsdt, &dtdt);
  norm = ds.cross(dt).normalized();
  Vector3f viewDir = (camera.position() - pos).normalized();

  // Get the projected direction in UV by representing it as the linear
  // combination of ds, dt
  MatrixXf M(3, 2);
  M.col(0) = ds;
  M.col(1) = dt;
  Vector2f view_uv = (M.transpose() * M).inverse() * M.transpose() * viewDir;
  view_uv = view_uv.normalized();

  /*=========================================================*/
  // 1. Compute dN/ds
  // Compute A (ds x dt)
  Vector3f A = ds.cross(dt);

  // Compute B (|| A ||)
  real_t B = A.norm();

  // Compute dA/ds
  Vector3f dA_ds = dsds.cross(dt) + ds.cross(dsdt);

  // Compute dB/ds
  real_t dB_ds = 1 / B * dA_ds.dot(A);

  // dN/ds
  Vector3f dN_ds = (B * dA_ds - dB_ds * A) / (B * B);

  // 2. Compute dN/dt
  Vector3f dA_dt = dsdt.cross(dt) + ds.cross(dtdt);
  real_t dB_dt = 1 / B * dA_dt.dot(A);
  Vector3f dN_dt = (B * dA_dt - dB_dt * A) / (B * B);

  // 3. Compose dN/dw
  Vector3f dN_dw = view_uv[0] * dN_ds + view_uv[1] * dN_dt;

  Vector3f df_dw = view_uv[0] * ds + view_uv[1] * dt;

  real_t ndotv = norm.dot(viewDir);
  real_t sin2theta = 1.0 - ndotv * ndotv;

  // View based weighting
  real_t k_r = df_dw.dot(dN_dw) / df_dw.dot(df_dw);
  k_r = k_r / sin2theta;

  return k_r;
}

void evaluate_radial_curvature_analytic(Mesh &mesh, Subdiv const &subdiv,
                                        Camera const &camera) {
  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");
  auto kappa = mesh.vertex_property<real_t>("v:kappa_ana");
  if (!kappa) {
    kappa = mesh.add_vertex_property<real_t>("v:kappa_ana", 0);
  }
  kappa.vector().assign(kappa.vector().size(), 0);

  for (uint32_t i = 0; i < mesh.n_vertices(); ++i) {
    Vertex v(i);

    // Find the valid parameterization (necessary for manifold boundary
    // vertices)
    Param_loc v_param_loc;
    auto hit = mesh.halfedges(v), hit_end = hit;
    do {
      v_param_loc = param_loc[*hit];
      if (v_param_loc.is_valid())
        break;
    } while (++hit != hit_end);

    contess_assert_msg(v_param_loc.is_valid(),
                       "evaluate_radial_curvature_analytic: Can't find valid "
                       "parameterization.");
    real_t k_r = evaluate_radial_curvature_analytic(mesh, subdiv, camera, v,
                                                    v_param_loc);

    // Average when close to EPs
    if (subdiv.backend_type() == Subdiv::Backend::LACEWELL) {
      Param_loc ep_v_param_loc = v_param_loc;
      ep_v_param_loc.uv[0] = (ep_v_param_loc.uv[0] < 0.5) ? 0 : 1;
      ep_v_param_loc.uv[1] = (ep_v_param_loc.uv[1] < 0.5) ? 0 : 1;
      if (subdiv.is_near_extraordinary(ep_v_param_loc) &&
          (v_param_loc.uv[0] < EP_GUARDING_ZONE ||
           v_param_loc.uv[1] < EP_GUARDING_ZONE ||
           1 - v_param_loc.uv[0] < EP_GUARDING_ZONE ||
           1 - v_param_loc.uv[1] < EP_GUARDING_ZONE)) {
        Param_loc param;
        int ep_adj_count = 0;
        real_t k_r_adj = 0;
        if (subdiv.is_near_extraordinary(v_param_loc)) {
          std::unordered_map<int, Param_loc> quads;
          auto hit = mesh.halfedges(v), hit_end = hit;
          do {
            param = param_loc[*hit];
            if (!param.is_valid() || quads.count(param.ptexIndex))
              continue;
            quads[param.ptexIndex] = param;
          } while (++hit != hit_end);
          for (auto const &q : quads) {
            auto param = q.second;
            subdiv.round_extraordinary(param);
            k_r_adj += evaluate_radial_curvature_analytic(mesh, subdiv, camera,
                                                          v, param);
            ep_adj_count++;
          }
        } else {
          Param_loc param = v_param_loc;
          subdiv.round_extraordinary(param);
          k_r_adj += evaluate_radial_curvature_analytic(mesh, subdiv, camera, v,
                                                        param);
          ep_adj_count++;
        }

        k_r_adj /= ep_adj_count;
        k_r = k_r_adj;
      }
    }

    kappa[v] = k_r;
  }
}

void verify_first_partial_derivatives(Mesh &mesh, Subdiv const &subdiv,
                                      real_t step_size) {
  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");
  auto ds_fd = mesh.vertex_property<Vector3f>("v:ds_fd");
  if (!ds_fd) {
    ds_fd = mesh.add_vertex_property<Vector3f>("v:ds_fd");
  }
  auto dt_fd = mesh.vertex_property<Vector3f>("v:dt_fd");
  if (!dt_fd) {
    dt_fd = mesh.add_vertex_property<Vector3f>("v:dt_fd");
  }
  dt_fd.vector().assign(dt_fd.vector().size(), Vector3f::Zero());
  auto ds_osd = mesh.vertex_property<Vector3f>("v:ds_osd");
  if (!ds_osd) {
    ds_osd = mesh.add_vertex_property<Vector3f>("v:ds_osd");
  }
  auto dt_osd = mesh.vertex_property<Vector3f>("v:dt_osd");
  if (!dt_osd) {
    dt_osd = mesh.add_vertex_property<Vector3f>("v:dt_osd");
  }
  dt_osd.vector().assign(dt_osd.vector().size(), Vector3f::Zero());

  for (uint32_t i = 0; i < mesh.n_vertices(); ++i) {
    Vertex v(i);
    Halfedge he = mesh.halfedge(v);
    const Param_loc &v_param_loc = param_loc[he];

    Vector3f pos, norm;
    subdiv.evaluateLimit(v_param_loc, pos, norm);

    // Find the viewDir in UV coordinate
    Vector3f ds, dt;
    subdiv.evaluateLimitFrame(v_param_loc, pos, ds, dt);

    ds_osd[v] = ds;
    dt_osd[v] = dt;

    // Take steps (assume the steps are small so we don't go outside the patch)
    std::vector<Vector2f> uv_steps({Vector2f(1, 0), Vector2f(0, 1)});
    for (size_t j = 0; j < uv_steps.size(); j++) {
      auto dir_uv = uv_steps[j];
      Vector3f pos_diff_pos, pos_diff_neg, norm_diff_pos, norm_diff_neg;
      Param_loc pos_p_pos, pos_p_neg;
      pos_p_pos.ptexIndex = pos_p_neg.ptexIndex = v_param_loc.ptexIndex;
      pos_p_pos.uv = step_size * dir_uv + v_param_loc.uv;
      pos_p_neg.uv = -step_size * dir_uv + v_param_loc.uv;

      real_t actual_stepsize = 2 * step_size;
      if (!pos_p_pos.is_valid()) {
        pos_p_pos = v_param_loc;
        actual_stepsize = step_size;
      }
      if (!pos_p_neg.is_valid()) {
        pos_p_neg = v_param_loc;
        actual_stepsize = step_size;
      }

      subdiv.evaluateLimit(pos_p_pos, pos_diff_pos, norm);
      subdiv.evaluateLimit(pos_p_neg, pos_diff_neg, norm);

      if (j == 0) {
        ds_fd[v] = (pos_diff_pos - pos_diff_neg) / actual_stepsize;
      } else if (j == 1) {
        dt_fd[v] = (pos_diff_pos - pos_diff_neg) / actual_stepsize;
      }
    }
  }
}

void verify_second_partial_derivatives(Mesh &mesh, Subdiv const &subdiv,
                                       real_t step_size) {
  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");
  auto dsds_fd = mesh.vertex_property<Vector3f>("v:dsds_fd");
  if (!dsds_fd) {
    dsds_fd = mesh.add_vertex_property<Vector3f>("v:dsds_fd");
  }
  dsds_fd.vector().assign(dsds_fd.vector().size(), Vector3f::Zero());
  auto dsdt_fd = mesh.vertex_property<Vector3f>("v:dsdt_fd");
  if (!dsdt_fd) {
    dsdt_fd = mesh.add_vertex_property<Vector3f>("v:dsdt_fd");
  }
  dsdt_fd.vector().assign(dsdt_fd.vector().size(), Vector3f::Zero());
  auto dtdt_fd = mesh.vertex_property<Vector3f>("v:dtdt_fd");
  if (!dtdt_fd) {
    dtdt_fd = mesh.add_vertex_property<Vector3f>("v:dtdt_fd");
  }
  dtdt_fd.vector().assign(dtdt_fd.vector().size(), Vector3f::Zero());
  auto dsds_osd = mesh.vertex_property<Vector3f>("v:dsds_osd");
  if (!dsds_osd) {
    dsds_osd = mesh.add_vertex_property<Vector3f>("v:dsds_osd");
  }
  dsds_osd.vector().assign(dsds_osd.vector().size(), Vector3f::Zero());
  auto dsdt_osd = mesh.vertex_property<Vector3f>("v:dsdt_osd");
  if (!dsdt_osd) {
    dsdt_osd = mesh.add_vertex_property<Vector3f>("v:dsdt_osd");
  }
  dsdt_osd.vector().assign(dsdt_osd.vector().size(), Vector3f::Zero());
  auto dtdt_osd = mesh.vertex_property<Vector3f>("v:dtdt_osd");
  if (!dtdt_osd) {
    dtdt_osd = mesh.add_vertex_property<Vector3f>("v:dtdt_osd");
  }
  dtdt_osd.vector().assign(dtdt_osd.vector().size(), Vector3f::Zero());

  for (uint32_t i = 0; i < mesh.n_vertices(); ++i) {
    Vertex v(i);
    Halfedge he = mesh.halfedge(v);
    const Param_loc &v_param_loc = param_loc[he];

    Vector3f pos, norm;
    subdiv.evaluateLimit(v_param_loc, pos, norm);

    // Find the viewDir in UV coordinate
    Vector3f ds, dt, dsds, dsdt, dtdt;
    subdiv.evaluateLimitFrame(v_param_loc, pos, ds, dt, &dsds, &dsdt, &dtdt);

    dsds_osd[v] = dsds;
    dsdt_osd[v] = dsdt;
    dtdt_osd[v] = dtdt;

    // Take steps (assume the steps are small so we don't go outside the patch)
    std::vector<Vector2f> uv_steps({Vector2f(1, 0), Vector2f(0, 1)});
    for (size_t j = 0; j < uv_steps.size(); j++) {
      auto dir_uv = uv_steps[j];
      Vector3f pos_diff_pos, pos_diff_neg, norm_diff_pos, norm_diff_neg;
      Param_loc pos_p_pos, pos_p_neg;
      pos_p_pos.ptexIndex = pos_p_neg.ptexIndex = v_param_loc.ptexIndex;
      pos_p_pos.uv = step_size * dir_uv + v_param_loc.uv;
      pos_p_neg.uv = -step_size * dir_uv + v_param_loc.uv;

      real_t actual_stepsize = 2 * step_size;
      if (!pos_p_pos.is_valid()) {
        pos_p_pos = v_param_loc;
        actual_stepsize = step_size;
      }
      if (!pos_p_neg.is_valid()) {
        pos_p_neg = v_param_loc;
        actual_stepsize = step_size;
      }

      Vector3f ds_diff_pos, ds_diff_neg, dt_diff_pos, dt_diff_neg;
      subdiv.evaluateLimitFrame(pos_p_pos, pos_diff_pos, ds_diff_pos,
                                dt_diff_pos);
      subdiv.evaluateLimitFrame(pos_p_neg, pos_diff_neg, ds_diff_neg,
                                dt_diff_neg);

      if (j == 0) {
        dsds_fd[v] = (ds_diff_pos - ds_diff_neg) / actual_stepsize;
      } else if (j == 1) {
        dsdt_fd[v] = (ds_diff_pos - ds_diff_neg) / actual_stepsize;
        dtdt_fd[v] = (dt_diff_pos - dt_diff_neg) / actual_stepsize;
      }
    }
  }
}

void approach_extraordinary_points(Mesh &mesh, Subdiv const &subdiv,
                                   Camera const &camera) {
  auto is_extraordinary = mesh.get_vertex_property<bool>("v:extraordinary");
  contess_assert(is_extraordinary);
  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");
  contess_assert(param_loc);
  auto kappa = mesh.vertex_property<real_t>("v:kappa_ana");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");

  std::vector<Halfedge> from_eps;
  for (uint32_t i = 0; i < mesh.n_vertices(); ++i) {
    Vertex v(i);
    if (!is_extraordinary[v])
      continue;
    auto hit = mesh.halfedges(v), hit_end = hit;
    do {
      from_eps.emplace_back(*hit);
    } while (++hit != hit_end);
  }

  real_t step_size = 1e-5;
  real_t small_zone = 1e-4;
  for (auto const he : from_eps) {
    Vertex ep_v = mesh.from_vertex(he);
    Vertex v = mesh.to_vertex(he);
    Param_loc ep_param = param_loc[he];
    Param_loc v_param = param_loc[mesh.opposite_halfedge(he)];
    Vertex start_v = ep_v;
    bool dist_printed = false;
    for (real_t step = step_size; step < 1; step += step_size) {
      if (step > small_zone)
        break;
      Param_loc param = ep_param;
      param.uv = ep_param.uv + step * (v_param.uv - ep_param.uv).normalized();
      if (!param.is_valid())
        break;

      Vector3f pos, norm;
      subdiv.evaluateLimit(param, pos, norm);
      real_t k_r = evaluate_radial_curvature_analytic(mesh, subdiv, camera,
                                                      Vertex(), param);
      if (start_v == ep_v && !dist_printed) {
        logger().info("Approaching {} by step: {}, {} => {}", ep_v, step_size,
                      (vpositions[ep_v] - pos).norm(), k_r);
        dist_printed = true;
      }
      start_v = mesh.split(mesh.find_edge(start_v, v), pos);
      kappa[start_v] = k_r;
    }
  }
}
