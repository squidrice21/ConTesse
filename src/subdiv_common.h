// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once
#include <memory>

#include "common.h"
#include "shape_utils.h"

static real_t EP_GUARDING_ZONE = 1e-5;

struct Param_loc {
  Param_loc(int idx = -1, real_t u = 0.f, real_t v = 0.f)
      : ptexIndex(idx), uv(Vector2f(u, v)) {}

  Param_loc(int idx, const Vector2f &st) : ptexIndex(idx), uv(st) {}

  int ptexIndex; ///< ptex face index
  Vector2f uv;   ///< parametric location on face

  bool is_valid() const {
    return (uv[0] >= 0 && uv[0] <= 1) && (uv[1] >= 0 && uv[1] <= 1) &&
           ptexIndex >= 0;
  }
  void clamp() {
    uv[0] = std::clamp(uv[0], 0., 1.);
    uv[1] = std::clamp(uv[1], 0., 1.);
  };
};

class SubdivOsd;
class SubdivLacewell;
class Mesh;

// So that we can just put this in a shared_ptr always
class Subdiv {
public:
  enum class Backend {
    LACEWELL,
  };

  Subdiv();
  ~Subdiv();

  bool load(const std::string &filename, Backend backend, int num_subdiv);
  bool load(Shape shape, Backend backend, int num_subdiv);

  bool evaluateLimit(const Param_loc &p, Vector3f &position,
                     Vector3f &normal) const;
  bool evaluateLimitFrame(const Param_loc &p, Vector3f &position, Vector3f &ds,
                          Vector3f &dt, Vector3f *dsds = nullptr,
                          Vector3f *dsdt = nullptr,
                          Vector3f *dtdt = nullptr) const;

  void determine_extraordinary_vertices(Mesh const &mesh);
  bool is_near_extraordinary(const Param_loc &p) const;

  void round_extraordinary(Param_loc &p) const;

  // This should not be used really. Just so that the viewer code compiles.
  SubdivOsd *as_subdiv_osd_workaround();

  Backend backend_type() const { return m_backend; }

private:
  friend void create_subdivided_mesh(std::shared_ptr<Subdiv>, Mesh &);

  Backend m_backend;
  std::unique_ptr<SubdivOsd> m_osd;
  std::unique_ptr<SubdivLacewell> m_lacewell;
};

void create_subdivided_mesh(std::shared_ptr<Subdiv>, Mesh &);
