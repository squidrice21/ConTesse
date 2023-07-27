// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#ifndef _SUBDIV_H
#define _SUBDIV_H

#include "subdiv_common.h"

#include <algorithm>
#include <opensubdiv/far/patchMap.h>
#include <opensubdiv/far/patchTableFactory.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/ptexIndices.h>
#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/topologyLevel.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class Mesh;

class SubdivOsd {

  //------------------------------------------------------------------------------
  // Vertex container implementation.
  //
  struct Vertex {

    // Minimal required interface ----------------------
    Vertex() {}

    Vertex(Vertex const &src) {
      position[0] = src.position[0];
      position[1] = src.position[1];
      position[2] = src.position[2];
    }

    void Clear(void * = 0) { position[0] = position[1] = position[2] = 0.0f; }

    void AddWithWeight(Vertex const &src, real_t weight) {
      position[0] += weight * src.position[0];
      position[1] += weight * src.position[1];
      position[2] += weight * src.position[2];
    }

    // Public interface ------------------------------------
    void SetPosition(real_t x, real_t y, real_t z) {
      position[0] = x;
      position[1] = y;
      position[2] = z;
    }

    real_t position[3];
  };

  //------------------------------------------------------------------------------
  // Limit frame container implementation -- this interface is not strictly
  // required but follows a similar pattern to Vertex.
  //
  struct LimitFrame {

    void Clear(void * = 0) {
      position[0] = position[1] = position[2] = 0.0f;
      deriv_s[0] = deriv_s[1] = deriv_s[2] = 0.0f;
      deriv_t[0] = deriv_t[1] = deriv_t[2] = 0.0f;
      deriv2_t[0] = deriv2_t[1] = deriv2_t[2] = 0.0f;
      deriv2_s[0] = deriv2_s[1] = deriv2_s[2] = 0.0f;
      deriv_st[0] = deriv_st[1] = deriv_st[2] = 0.0f;
    }

    void AddWithWeight(Vertex const &src, real_t p_weight, real_t ds_Weight,
                       real_t dt_Weight, real_t d2s_Weight, real_t dst_Weight,
                       real_t d2t_Weight) {

      for (int i = 0; i < 3; i++) {
        position[i] += p_weight * src.position[i];
        deriv_s[i] += ds_Weight * src.position[i];
        deriv_t[i] += dt_Weight * src.position[i];
        deriv2_s[i] += d2s_Weight * src.position[i];
        deriv_st[i] += dst_Weight * src.position[i];
        deriv2_t[i] += d2t_Weight * src.position[i];
      }
    }

    real_t position[3], deriv_s[3], deriv_t[3], deriv2_s[3], deriv2_t[3],
        deriv_st[3];
  };

public:
  SubdivOsd();

  ~SubdivOsd();

  bool load(const std::string &filename);
  bool load(Shape shape);

  Mesh *refineUniform(int maxlevel);
  void refineUniform(int maxlevel, Mesh &mesh);

  size_t get_num_vertices() const { return m_shape->GetNumVertices(); }
  size_t get_num_faces() const { return m_shape->GetNumFaces(); }

  MapXu get_indices() const {
    return MapXu(m_shape->faceverts.data(), 1, m_shape->faceverts.size());
  }
  const std::vector<int> &get_nvertsPerFace() const {
    return m_shape->nvertsPerFace;
  }
  MapXf get_positions() const {
    return MapXf(m_shape->verts.data(), 3, m_shape->GetNumVertices());
  }

  bool evaluateLimit(const Param_loc &p, Vector3f &position,
                     Vector3f &normal) const;
  bool evaluateLimitFrame(const Param_loc &p, Vector3f &position, Vector3f &ds,
                          Vector3f &dt, Vector3f *dsds = nullptr,
                          Vector3f *dsdt = nullptr,
                          Vector3f *dtdt = nullptr) const;

  void determine_extraordinary_vertices(Mesh const &mesh);
  bool is_near_extraordinary(const Param_loc &p) const;

private:
  void adaptiveRefinement();

  std::unique_ptr<Shape> m_shape;

  OpenSubdiv::Sdc::SchemeType m_type = OpenSubdiv::Sdc::SCHEME_CATMARK;
  typedef OpenSubdiv::Sdc::Options SdcOptions;
  SdcOptions m_options;
  OpenSubdiv::Far::TopologyRefiner *m_refiner;
  OpenSubdiv::Far::PatchTable *patchTable;
  OpenSubdiv::Far::PatchMap *patchmap;
  OpenSubdiv::Far::PtexIndices *ptexIndices;
  std::vector<Vertex> verts;
  std::unordered_map<int, std::vector<Param_loc>> extraordinary_parameters;
};

#endif //_SUBDIV_H
