// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.


#include "subdiv_lacewell.h"
#include "common.h"
#include "mesh.h"

#include <fstream>
#include <limits>
#include <memory>
#include <spdlog/fmt/ostr.h>
#include <subdeval/lib/subdeval/Subd.h>
#include <unordered_map>
#include <vector>

#include "subdiv_common.h"
#include "tag_cusp_facing.h"

SubdivLacewell::SubdivLacewell() = default;
SubdivLacewell::~SubdivLacewell() = default;

void SubdivLacewell::load(const std::string &filename, int num_subdiv) {
  std::string extension;
  if (filename.size() > 4)
    extension = str_tolower(filename.substr(filename.size() - 4));

  std::string str;
  if (extension == ".obj") {
    std::ifstream ifs(filename.c_str());
    if (ifs) {
      std::stringstream ss;
      ss << ifs.rdbuf();
      ifs.close();
      cout << "Reading " << filename << endl;
      str = ss.str();
    } else {
      cout << "Error in reading " << filename << endl;
      return;
    }
  } else {
    str = filename;
  }

  auto shape = std::unique_ptr<Shape>(Shape::parseObj(str.c_str(), kCatmark));
  load(std::move(*shape), num_subdiv);
}

void SubdivLacewell::load(Shape shape, int num_subdiv) {
  m_shape = std::move(shape);
  std::vector<int> faceverts;
  faceverts.reserve(m_shape.faceverts.size());
  for (const auto &v : m_shape.faceverts) {
    faceverts.push_back(static_cast<int>(v));
  }
  // Here use osubdiv to subdiv num_subdiv times (get rid of triangles) and only
  // get the position and connectiviy out. Then use to create this. For now
  // let's ignore the initial subdiv.
  m_impl = std::make_unique<Subd>(
      static_cast<int>(m_shape.verts.size() / 3), m_shape.verts.data(),
      static_cast<int>(m_shape.nvertsPerFace.size()),
      m_shape.nvertsPerFace.data(), faceverts.data());
  m_impl->subdivide(num_subdiv);
}

void SubdivLacewell::create_mesh(Mesh &mesh) {
  mesh.clear();

  std::vector<surface_mesh::Surface_mesh::Vertex> map;
  {
    map.resize(m_impl->nverts());
    const double *v = m_impl->limitverts();
    for (int i = 0; i < m_impl->nverts(); i++, v += 3) {
      map[i] = mesh.add_vertex(Vector3f(v[0], v[1], v[2]));
    }
  }

  std::unordered_map<int, Face> impl2mesh;
  {
    const int *nv = m_impl->nvertsPerFace();
    const int *fv = m_impl->faceverts();
    for (int i = 0; i < m_impl->nfaces(); i++, nv++) {
      assert(*nv == 4);
      Face f = mesh.add_quad(map[fv[0]], map[fv[1]], map[fv[2]], map[fv[3]]);
      impl2mesh[i] = f;
      fv += *nv;
    }
  }

  {
    auto param_loc = mesh.add_halfedge_property<Param_loc>("h:param_loc");

    const int *nv = m_impl->nvertsPerFace();
    const int *fv = m_impl->faceverts();
    double uvs_quad[][2] = {
        {0, 0},
        {1, 0},
        {1, 1},
        {0, 1},
    };
    double uvs_tri[][2] = {
        {0, 0},
        {1, 0},
        {0, 1},
    };

    for (int i = 0; i < m_impl->nfaces(); i++, nv++) {
      Face f = impl2mesh[i];
      Vertex v = map[fv[0]];
      Halfedge h;
      auto hit = mesh.halfedges(f), hit_end = hit;
      do {
        if (v == mesh.from_vertex(*hit)) {
          h = *hit;
          break;
        }
      } while (++hit != hit_end);

      contess_assert(h.is_valid());

      for (int vi = 0; vi < *nv; ++vi) {
        real_t u = uvs_quad[vi][0];
        real_t v = uvs_quad[vi][1];
        if (*nv == 3) {
          u = uvs_tri[vi][0];
          v = uvs_tri[vi][1];
        }
        contess_assert(*nv == 4 || *nv == 3);

        param_loc[h] = Param_loc(i, u, v);
        h = mesh.next_halfedge(h);
      }
      fv += *nv;
    }

    // Check parametric locations
    surface_mesh::Surface_mesh::Halfedge_iterator hit;
    for (hit = mesh.halfedges_begin(); hit != mesh.halfedges_end(); ++hit)
      if (!mesh.is_boundary(*hit))
        assert(param_loc[*hit].ptexIndex != -1);
  }

  // Determine and store EPs
  determine_extraordinary_vertices(mesh);

  m_mesh = &mesh;
}

bool SubdivLacewell::evaluateLimit(const Param_loc &p, Vector3f &position,
                                   Vector3f &normal) const {
  real_t pos[3], dPdU[3], dPdV[3];
  m_impl->eval(p.ptexIndex, p.uv[0], p.uv[1], pos, dPdU, dPdV);

  position << pos[0], pos[1], pos[2];
  Vector3f du, dv;
  du << dPdU[0], dPdU[1], dPdU[2];
  dv << dPdV[0], dPdV[1], dPdV[2];
  normal = du.cross(dv).normalized();

  // Find a near point to avoid zero normal
  auto p_check = p;
  p_check.uv[0] = (p_check.uv[0] < 0.5) ? 0 : 1;
  p_check.uv[1] = (p_check.uv[1] < 0.5) ? 0 : 1;
  if (is_near_extraordinary(p_check) &&
      (p.uv[0] < EP_GUARDING_ZONE || p.uv[1] < EP_GUARDING_ZONE ||
      1 - p.uv[0] < EP_GUARDING_ZONE || 1 - p.uv[1] < EP_GUARDING_ZONE)) {
    auto param_loc = m_mesh->get_halfedge_property<Param_loc>("h:param_loc");
    contess_assert(
        ep_vertices.count(p_check.ptexIndex * 100 + p_check.uv[0] * 10 + p_check.uv[1]));
    Vertex ep_v(ep_vertices.at(p_check.ptexIndex * 100 + p_check.uv[0] * 10 + p_check.uv[1]));

    position << 0, 0, 0;
    normal << 0, 0, 0;
    int ep_adj_count = 0;
    auto hit = m_mesh->halfedges(ep_v), hit_end = hit;
    Vector3f pos_adj;
    Vector3f norm_adj;
    if (is_near_extraordinary(p)) {
      std::unordered_map<int, Param_loc> quads;
      do {
        Param_loc p_near = param_loc[*hit];
        if (!p_near.is_valid() || quads.count(p_near.ptexIndex))
          continue;
        quads[p_near.ptexIndex] = p_near;
      } while (++hit != hit_end);

      for (auto const &q : quads) {
        Param_loc p_near = q.second;
        round_extraordinary(p_near);
        m_impl->eval(p_near.ptexIndex, p_near.uv[0], p_near.uv[1], pos, dPdU,
                     dPdV);

        pos_adj << pos[0], pos[1], pos[2];
        du << dPdU[0], dPdU[1], dPdU[2];
        dv << dPdV[0], dPdV[1], dPdV[2];
        norm_adj = du.cross(dv).normalized();

        position += pos_adj;
        normal += norm_adj;

        ep_adj_count++;
      }
    } else {
      Param_loc p_near = p;
      round_extraordinary(p_near);
      m_impl->eval(p_near.ptexIndex, p_near.uv[0], p_near.uv[1], pos, dPdU,
                   dPdV);

      pos_adj << pos[0], pos[1], pos[2];
      du << dPdU[0], dPdU[1], dPdU[2];
      dv << dPdV[0], dPdV[1], dPdV[2];
      norm_adj = du.cross(dv).normalized();

      position += pos_adj;
      normal += norm_adj;

      ep_adj_count++;
    }
    position /= ep_adj_count;
    normal /= ep_adj_count;
    normal.normalize();
  }

  // save EV locations and return true if at EV exactly
  return is_near_extraordinary(p);
}

bool SubdivLacewell::evaluateLimitFrame(const Param_loc &in_p,
                                        Vector3f &position, Vector3f &ds,
                                        Vector3f &dt, Vector3f *dsds,
                                        Vector3f *dsdt, Vector3f *dtdt) const {
  Param_loc p = in_p;
  real_t local_step_size = step_size;

  if (is_near_extraordinary(p) && (dsds || dsdt || dtdt)) {
    local_step_size = std::max(step_size, 1e-6);
  }
  real_t pos[3], dPdU[3], dPdV[3];
  m_impl->eval(p.ptexIndex, p.uv[0], p.uv[1], pos, dPdU, dPdV);

  position << pos[0], pos[1], pos[2];
  Vector3f du, dv;
  du << dPdU[0], dPdU[1], dPdU[2];
  dv << dPdV[0], dPdV[1], dPdV[2];
  ds = du;
  dt = dv;

  auto fd_1st_order = [&](const Param_loc &p, Vector2f const &step,
                          Vector3f &dPds) {
    Param_loc pos_p_pos, pos_p_neg;
    pos_p_pos.ptexIndex = pos_p_neg.ptexIndex = p.ptexIndex;
    pos_p_pos.uv = local_step_size * step + p.uv;
    pos_p_neg.uv = -local_step_size * step + p.uv;

    // If we are in the intersection point of multiple patches and the
    // current patch doesn't contain both steps
    contess_assert(pos_p_pos.is_valid() || pos_p_neg.is_valid());

    real_t actual_stepsize = 2 * local_step_size;
    if (!pos_p_pos.is_valid()) {
      pos_p_pos = p;
      actual_stepsize = local_step_size;
    }
    if (!pos_p_neg.is_valid()) {
      pos_p_neg = p;
      actual_stepsize = local_step_size;
    }

    Vector3f dP_pos, dP_neg;
    real_t pos[3], dPdU[3], dPdV[3];
    m_impl->eval(pos_p_pos.ptexIndex, pos_p_pos.uv[0], pos_p_pos.uv[1], pos,
                 dPdU, dPdV);
    dP_pos << pos[0], pos[1], pos[2];

    m_impl->eval(pos_p_neg.ptexIndex, pos_p_neg.uv[0], pos_p_neg.uv[1], pos,
                 dPdU, dPdV);
    dP_neg << pos[0], pos[1], pos[2];

    dPds = (dP_pos - dP_neg) / actual_stepsize;
  };

  // Finite difference to get the 2nd order derivatives
  auto fd_in_dir = [&](const Param_loc &p, Vector2f const &step, bool is_ds,
                       Vector3f &dPds) {
    Param_loc pos_p_pos, pos_p_neg;
    pos_p_pos.ptexIndex = pos_p_neg.ptexIndex = p.ptexIndex;
    pos_p_pos.uv = local_step_size * step + p.uv;
    pos_p_neg.uv = -local_step_size * step + p.uv;
    int step_side = 0;

    // If we are in the intersection point of multiple patches and the
    // current patch doesn't contain both steps
    contess_assert(pos_p_pos.is_valid() || pos_p_neg.is_valid());

    real_t actual_stepsize = 2 * local_step_size;
    if (!pos_p_pos.is_valid()) {
      pos_p_pos = p;
      actual_stepsize = local_step_size;
      step_side = -1;
    }
    if (!pos_p_neg.is_valid()) {
      pos_p_neg = p;
      actual_stepsize = local_step_size;
      step_side = 1;
    }

    Vector3f dP_pos, dP_neg;
    real_t pos[3], dPdU[3], dPdV[3];

    // At EPs, finite difference from the position since the 1st derivative is
    // also off
    // if (is_near_extraordinary(p)) {
    if (false) {
      Param_loc pos_p_next;
      pos_p_next.ptexIndex = p.ptexIndex;
      if (step_side == 1) {
        pos_p_next.uv = local_step_size * step + pos_p_pos.uv;
      } else if (step_side == -1) {
        pos_p_next.uv = -local_step_size * step + pos_p_neg.uv;
      } else {
        contess_assert_msg(0, "EPs must be at quad UV corner.");
      }

      // Non cross term
      if ((std::abs(step[0] - 1) < std::numeric_limits<real_t>::epsilon()) ==
          is_ds) {
        Vector3f P_pos, P_neg, P_3;
        m_impl->eval(pos_p_pos.ptexIndex, pos_p_pos.uv[0], pos_p_pos.uv[1], pos,
                     dPdU, dPdV);
        P_pos << pos[0], pos[1], pos[2];

        m_impl->eval(pos_p_neg.ptexIndex, pos_p_neg.uv[0], pos_p_neg.uv[1], pos,
                     dPdU, dPdV);
        P_neg << pos[0], pos[1], pos[2];

        m_impl->eval(pos_p_next.ptexIndex, pos_p_next.uv[0], pos_p_next.uv[1],
                     pos, dPdU, dPdV);
        P_3 << pos[0], pos[1], pos[2];

        // First order accuracy one sided FD
        Vector3f P_1 = (step_side == 1) ? P_neg : P_pos;
        Vector3f P_2 = (step_side == 1) ? P_pos : P_neg;
        dPds = (P_1 - 2 * P_2 + P_3) / (local_step_size * local_step_size);
      }
      // Cross term
      else {
        Vector2f order1_step = Vector2f(1, 1) - step;
        Vector3f P_pos, P_neg;
        fd_1st_order(pos_p_pos, order1_step, P_pos);
        fd_1st_order(pos_p_neg, order1_step, P_neg);
        dPds = (P_pos - P_neg) / actual_stepsize;
      }
    } else {
      m_impl->eval(pos_p_pos.ptexIndex, pos_p_pos.uv[0], pos_p_pos.uv[1], pos,
                   dPdU, dPdV);

      if (is_ds)
        dP_pos << dPdU[0], dPdU[1], dPdU[2];
      else
        dP_pos << dPdV[0], dPdV[1], dPdV[2];

      m_impl->eval(pos_p_neg.ptexIndex, pos_p_neg.uv[0], pos_p_neg.uv[1], pos,
                   dPdU, dPdV);
      if (is_ds)
        dP_neg << dPdU[0], dPdU[1], dPdU[2];
      else
        dP_neg << dPdV[0], dPdV[1], dPdV[2];

      dPds = (dP_pos - dP_neg) / actual_stepsize;
    }

    local_step_size = step_size;
  };
  if (dsds)
    fd_in_dir(p, Vector2f(1, 0), true, *dsds);
  if (dsdt)
    fd_in_dir(p, Vector2f(0, 1), true, *dsdt);
  if (dtdt)
    fd_in_dir(p, Vector2f(0, 1), false, *dtdt);

  return is_near_extraordinary(p);
}

void SubdivLacewell::determine_extraordinary_vertices(Mesh const &mesh) {
  std::vector<Vertex> extraordinary_vertices;
  tag_extraordinary_quad(mesh, extraordinary_vertices);

  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");
  contess_assert(param_loc);

  for (auto const &v : extraordinary_vertices) {
    auto hit = mesh.halfedges(v), hit_end = hit;
    do {
      if (param_loc[*hit].is_valid()) {
        if (!extraordinary_parameters.count(param_loc[*hit].ptexIndex))
          extraordinary_parameters[param_loc[*hit].ptexIndex] =
              std::vector<Param_loc>();
        extraordinary_parameters[param_loc[*hit].ptexIndex].emplace_back(
            param_loc[*hit]);
        auto he2 = mesh.next_halfedge(*hit);
        auto he3 = mesh.next_halfedge(he2);
        auto he4 = mesh.next_halfedge(he3);
        extraordinary_parameters[param_loc[*hit].ptexIndex].emplace_back(
            param_loc[he2]);
        extraordinary_parameters[param_loc[*hit].ptexIndex].emplace_back(
            param_loc[he3]);

        // Second triangle face of the quad
        extraordinary_parameters[param_loc[*hit].ptexIndex].emplace_back(
            param_loc[he3]);
        extraordinary_parameters[param_loc[*hit].ptexIndex].emplace_back(
            param_loc[he4]);
        extraordinary_parameters[param_loc[*hit].ptexIndex].emplace_back(
            param_loc[*hit]);

        // The exact at EP check
        if (!exact_extraordinary_parameters.count(param_loc[*hit].ptexIndex))
          exact_extraordinary_parameters[param_loc[*hit].ptexIndex] =
              std::vector<Param_loc>();
        exact_extraordinary_parameters[param_loc[*hit].ptexIndex].emplace_back(
            param_loc[*hit]);

        // This is a ugly workaround to get vertex in the evaluate frame function...
        // to avoid changing the existing APIs
        // Indexed by face idx, u, v
        ep_vertices[param_loc[*hit].ptexIndex * 100 +
                    param_loc[*hit].uv[0] * 10 + param_loc[*hit].uv[1]] =
            v.idx();
      }
    } while (++hit != hit_end);
  }
}

bool SubdivLacewell::is_near_extraordinary(const Param_loc &p) const {
  // Debug
  return false;
  //

  // TODO: exactly at the EP
  if (exact_extraordinary_parameters.find(p.ptexIndex) ==
      exact_extraordinary_parameters.end())
    return false;

  // Lacewell code rounds points to EPs with distance of 1e-6
  auto rounded_uv = p.uv;
  if (rounded_uv[0] < 1e-6)
    rounded_uv[0] = 0;
  if (rounded_uv[1] < 1e-6)
    rounded_uv[1] = 0;
  for (auto const &param : exact_extraordinary_parameters.at(p.ptexIndex))
    if ((rounded_uv - param.uv).norm() < std::numeric_limits<real_t>::epsilon())
      return true;
  return false;

  // 1-ring check
  if (extraordinary_parameters.find(p.ptexIndex) ==
      extraordinary_parameters.end())
    return false;

  // Check every three UVs (from the same triangle)
  for (size_t i = 0; i < extraordinary_parameters.at(p.ptexIndex).size();
       i += 3) {
    std::vector<Vector2f> uvs;
    for (size_t j = 0; j < 3; j++) {
      uvs.emplace_back(extraordinary_parameters.at(p.ptexIndex).at(i + j).uv);
    }

    // Check if the input uv is within this triangle
    if ((igl::predicates::orient2d(uvs[0], uvs[1], p.uv) ==
             igl::predicates::Orientation::POSITIVE ||
         igl::predicates::orient2d(uvs[0], uvs[1], p.uv) ==
             igl::predicates::Orientation::COLLINEAR) &&
        (igl::predicates::orient2d(uvs[1], uvs[2], p.uv) ==
             igl::predicates::Orientation::POSITIVE ||
         igl::predicates::orient2d(uvs[1], uvs[2], p.uv) ==
             igl::predicates::Orientation::COLLINEAR) &&
        (igl::predicates::orient2d(uvs[2], uvs[0], p.uv) ==
             igl::predicates::Orientation::POSITIVE ||
         igl::predicates::orient2d(uvs[2], uvs[0], p.uv) ==
             igl::predicates::Orientation::COLLINEAR))
      return true;
  }

  return false;
}

void SubdivLacewell::round_extraordinary(Param_loc &p) const {
  Param_loc ep = p;
  ep.uv.setZero();
  // Debug
  return;
  //
  if (!(is_near_extraordinary(ep) &&
        (p.uv[0] < EP_GUARDING_ZONE || p.uv[1] < EP_GUARDING_ZONE)))
    return;
  for (size_t i = 0; i < 2; i++)
    if (p.uv[i] < EP_GUARDING_ZONE)
      p.uv[i] += EP_GUARDING_ZONE;
}