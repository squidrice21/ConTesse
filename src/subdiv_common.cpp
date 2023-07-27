// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "subdiv_common.h"
#include "mesh.h"
#include "subdiv_lacewell.h"
#include "subdiv_osd.h"

#include <sstream>

Subdiv::Subdiv() = default;
Subdiv::~Subdiv() = default;

bool Subdiv::load(const std::string &filename, Backend backend,
                  int num_subdiv) {
  std::string extension;
  if (filename.size() > 4) {
    extension = str_tolower(filename.substr(filename.size() - 4));
  }

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
      return false;
    }
  } else {
    str = filename;
  }

  auto shape = std::unique_ptr<Shape>(Shape::parseObj(str.c_str(), kCatmark));
  return load(std::move(*shape), backend, num_subdiv);
}

bool Subdiv::load(Shape shape, Backend backend, int num_subdiv) {
  m_backend = backend;

  // Check if the mesh has triangles
  bool is_mixed = false;
  for (auto nverts : shape.nvertsPerFace) {
    if (nverts == 3) {
      is_mixed = true;
    } else if (nverts != 4) {
      is_mixed = true;
    }
  }

  if (m_backend == Backend::LACEWELL) {
    m_lacewell = std::make_unique<SubdivLacewell>();
    contess_assert(num_subdiv >= 1);

    // Lacewell cannot handle mixed meshes. So we need to do a pre subdiv for
    // it.
    if (is_mixed) {
      // since osd does one subdiv, reduce total num subdivs.
      num_subdiv -= 1;

      // Do one subdiv with osd to create quads.
      SubdivOsd subdiv;
      subdiv.load(shape);
      Mesh mesh;
      subdiv.refineUniform(1, mesh);

      // Now use this mesh to create a shape
      std::stringstream obj_as_str;
      surface_mesh::write_obj(mesh, obj_as_str);
      obj_as_str.seekg(0);
      shape = std::move(*std::unique_ptr<Shape>(
          Shape::parseObj(obj_as_str.str().c_str(), kCatmark)));
    }

    m_lacewell->load(std::move(shape), num_subdiv);
    return true;
  } else {
    contess_assert(0);
  }
}

bool Subdiv::evaluateLimit(const Param_loc &p, Vector3f &position,
                           Vector3f &normal) const {
  if (m_backend == Backend::LACEWELL) {
    contess_assert(m_lacewell);
    return m_lacewell->evaluateLimit(p, position, normal);
  } else {
    contess_assert(0);
  }
}

bool Subdiv::evaluateLimitFrame(const Param_loc &p, Vector3f &position,
                                Vector3f &ds, Vector3f &dt, Vector3f *dsds,
                                Vector3f *dsdt, Vector3f *dtdt) const {
  if (m_backend == Backend::LACEWELL) {
    contess_assert(m_lacewell);
    return m_lacewell->evaluateLimitFrame(p, position, ds, dt, dsds, dsdt,
                                          dtdt);
  } else {
    contess_assert(0);
  }
}

void Subdiv::determine_extraordinary_vertices(Mesh const &mesh) {
  if (m_backend == Backend::LACEWELL) {
    // Should already be determined at this point
    // (when the mesh is created)
    // contess_assert(m_lacewell);
    // m_lacewell->determine_extraordinary_vertices(mesh);
    // contess_assert(0);
  } else {
    contess_assert(0);
  }
}

bool Subdiv::is_near_extraordinary(const Param_loc &p) const {
  if (m_backend == Backend::LACEWELL) {
    contess_assert(m_lacewell);
    return m_lacewell->is_near_extraordinary(p);
  } else {
    contess_assert(0);
  }
}

void Subdiv::round_extraordinary(Param_loc &p) const {
  if (m_backend == Backend::LACEWELL) {
    contess_assert(m_lacewell);
    m_lacewell->round_extraordinary(p);
  } else {
    contess_assert(0);
  }
}

SubdivOsd *Subdiv::as_subdiv_osd_workaround() {
  logger().warn("This is only a workaroudn so that the viewer compiles. This "
                "should not be used");
  contess_assert(m_osd);
  return m_osd.get();
}

void create_subdivided_mesh(std::shared_ptr<Subdiv> s, Mesh &mesh) {
  contess_assert(s);
  if (s->m_backend == Subdiv::Backend::LACEWELL) {
    s->m_lacewell->create_mesh(mesh);
  } else {
    contess_assert(0);
  }
  mesh.set_subdivision(s);
}
