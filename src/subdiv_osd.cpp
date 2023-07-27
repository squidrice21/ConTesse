// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "subdiv_osd.h"
#include "common.h"
#include "far_utils.h"
#include "mesh.h"
#include "surface_mesh/surface_mesh.h"

#include <Eigen/Eigenvalues>
#include <igl/predicates/predicates.h>
#include <spdlog/fmt/ostr.h>

#include <assert.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace OpenSubdiv;

SubdivOsd::SubdivOsd()
    : m_refiner(nullptr), patchTable(nullptr), patchmap(nullptr),
      ptexIndices(nullptr) {}

SubdivOsd::~SubdivOsd() {
  if (m_refiner)
    delete m_refiner;
  if (patchTable)
    delete patchTable;
  if (patchmap)
    delete patchmap;
  if (ptexIndices)
    delete ptexIndices;
}

bool SubdivOsd::load(const std::string &filename) {
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
      return false;
    }
  } else {
    str = filename;
  }

  m_shape = std::unique_ptr<Shape>(Shape::parseObj(str.c_str(), kCatmark));
  return this->load(*m_shape);
}

bool SubdivOsd::load(Shape shape) {
  m_shape = std::make_unique<Shape>(std::move(shape));

  // create Far mesh (topology)
  m_type = GetSdcType(*m_shape);
  m_options = GetSdcOptions(*m_shape);

  m_refiner = Far::TopologyRefinerFactory<Shape>::Create(
      *m_shape, Far::TopologyRefinerFactory<Shape>::Options(m_type, m_options));
  assert(m_refiner);

  adaptiveRefinement();
  return true;
}

void RefinedToCoarseUV(Far::PatchParam const &param, real_t &u, real_t &v) {
  // See: Far::PatchParam::Normalize()
  real_t frac = param.GetParamFraction();

  real_t pu = (real_t)param.GetU() * frac;
  real_t pv = (real_t)param.GetV() * frac;

  u = u * frac + pu, v = v * frac + pv;
}

Mesh *SubdivOsd::refineUniform(int maxlevel) {
  // Instantiate a FarTopologyRefiner from the descriptor
  Far::TopologyRefiner *refiner = Far::TopologyRefinerFactory<Shape>::Create(
      *m_shape, Far::TopologyRefinerFactory<Shape>::Options(m_type, m_options));
  delete m_refiner;
  m_refiner = refiner;

  // Uniformly refine the topolgy up to 'maxlevel'
  // note: fullTopologyInLastLevel must be true to work with face-varying data
  Far::TopologyRefiner::UniformOptions refineOptions(maxlevel);
  refineOptions.fullTopologyInLastLevel = true;
  refiner->RefineUniform(refineOptions);

  // Allocate a buffer for vertex primvar data. The buffer length is set to
  // be the sum of all children vertices up to the highest level of refinement.
  std::vector<Vertex> vbuffer(refiner->GetNumVerticesTotal());
  Vertex *verts = &vbuffer[0];

  // Initialize coarse mesh positions
  int nCoarseVerts = m_shape->GetNumVertices();
  for (int i = 0; i < nCoarseVerts; ++i) {
    verts[i].SetPosition(m_shape->verts[i * 3], m_shape->verts[i * 3 + 1],
                         m_shape->verts[i * 3 + 2]);
  }

  // Interpolate vertex primvar data
  Far::PrimvarRefinerReal<real_t> primvarRefiner(*refiner);

  Vertex *src = verts;
  for (int level = 1; level <= maxlevel; ++level) {
    Vertex *dst = src + refiner->GetLevel(level - 1).GetNumVertices();
    primvarRefiner.Interpolate(level, src, dst);
    src = dst;
  }

  // Compute limit normals
  Far::TopologyLevel const &refLastLevel = refiner->GetLevel(maxlevel);
  int nverts = refLastLevel.GetNumVertices();
  // int nfaces = refLastLevel.GetNumFaces();
  int firstOfLastVerts = refiner->GetNumVerticesTotal() - nverts;
  std::vector<Vertex> fineLimitPos(nverts);
  std::vector<Vertex> fineDu(nverts);
  std::vector<Vertex> fineDv(nverts);

  primvarRefiner.Limit(&verts[firstOfLastVerts], fineLimitPos, fineDu, fineDv);

  Mesh *mesh = new Mesh();
  auto vnormals = mesh->add_vertex_property<Vector3f>("v:normal");

  std::vector<surface_mesh::Surface_mesh::Vertex> map;
  map.resize(nverts);

  for (int vert = 0; vert < nverts; ++vert) {

    map[vert] = mesh->add_vertex(Vector3f(fineLimitPos[vert].position[0],
                                          fineLimitPos[vert].position[1],
                                          fineLimitPos[vert].position[2]));

    Vector3f du = Vector3f(fineDu[vert].position[0], fineDu[vert].position[1],
                           fineDu[vert].position[2]);
    Vector3f dv = Vector3f(fineDv[vert].position[0], fineDv[vert].position[1],
                           fineDv[vert].position[2]);

    vnormals[map[vert]] = du.cross(dv).normalized();
  }

  Far::PatchTableFactory::Options ptOptions;
  ptOptions.SetPatchPrecision<real_t>();
  Far::PatchTable const *patchTable =
      Far::PatchTableFactory::Create(*refiner, ptOptions);

  // uv indices for a refined quad face
  real_t uv[4][2] = {
      {0.0, 0.0},
      {1.0, 0.0},
      {1.0, 1.0},
      {0.0, 1.0},
  };

  auto param_loc = mesh->add_halfedge_property<Param_loc>("h:param_loc");

  for (int array = 0; array < patchTable->GetNumPatchArrays(); ++array) {
    for (int patch = 0; patch < patchTable->GetNumPatches(array); ++patch) {

      Far::ConstIndexArray faceVerts =
          patchTable->GetPatchVertices(array, patch);

      auto f = mesh->add_quad(
          map[faceVerts[0] - nCoarseVerts], map[faceVerts[1] - nCoarseVerts],
          map[faceVerts[2] - nCoarseVerts], map[faceVerts[3] - nCoarseVerts]);

      Far::PatchParam param = patchTable->GetPatchParam(array, patch);
      auto h = mesh->halfedge(f);
      for (int quadVert = 0; quadVert < 4; ++quadVert) {
        if (map[faceVerts[0] - nCoarseVerts] == mesh->from_vertex(h))
          break;
        h = mesh->next_halfedge(h);
      }
      assert(map[faceVerts[0] - nCoarseVerts] == mesh->from_vertex(h));
      for (int quadVert = 0; quadVert < 4; ++quadVert) {
        real_t u = uv[quadVert][0];
        real_t v = uv[quadVert][1];
        RefinedToCoarseUV(param, u, v);
        param_loc[h] = Param_loc(param.GetFaceId(), u, v);
        h = mesh->next_halfedge(h);
      }
    }
  }

  // Check parametric locations
  surface_mesh::Surface_mesh::Halfedge_iterator hit;
  for (hit = mesh->halfedges_begin(); hit != mesh->halfedges_end(); ++hit)
    if (!mesh->is_boundary(*hit))
      assert(param_loc[*hit].ptexIndex != -1);

  return mesh;
}

void SubdivOsd::refineUniform(int maxlevel, Mesh &mesh) {
  // Instantiate a FarTopologyRefiner from the descriptor
  Far::TopologyRefiner *refiner = Far::TopologyRefinerFactory<Shape>::Create(
      *m_shape, Far::TopologyRefinerFactory<Shape>::Options(m_type, m_options));

  delete m_refiner;
  m_refiner = refiner;

  // Uniformly refine the topolgy up to 'maxlevel'
  // note: fullTopologyInLastLevel must be true to work with face-varying data
  Far::TopologyRefiner::UniformOptions refineOptions(maxlevel);
  refineOptions.fullTopologyInLastLevel = true;
  refiner->RefineUniform(refineOptions);

  // Allocate a buffer for vertex primvar data. The buffer length is set to
  // be the sum of all children vertices up to the highest level of refinement.
  std::vector<Vertex> vbuffer(refiner->GetNumVerticesTotal());
  Vertex *verts = &vbuffer[0];

  // Initialize coarse mesh positions
  int nCoarseVerts = m_shape->GetNumVertices();
  for (int i = 0; i < nCoarseVerts; ++i) {
    verts[i].SetPosition(m_shape->verts[i * 3], m_shape->verts[i * 3 + 1],
                         m_shape->verts[i * 3 + 2]);
  }

  // Interpolate vertex primvar data
  Far::PrimvarRefinerReal<real_t> primvarRefiner(*refiner);

  Vertex *src = verts;
  for (int level = 1; level <= maxlevel; ++level) {
    Vertex *dst = src + refiner->GetLevel(level - 1).GetNumVertices();
    primvarRefiner.Interpolate(level, src, dst);
    src = dst;
  }

  // Compute limit normals
  Far::TopologyLevel const &refLastLevel = refiner->GetLevel(maxlevel);
  int nverts = refLastLevel.GetNumVertices();
  // int nfaces = refLastLevel.GetNumFaces();
  int firstOfLastVerts = refiner->GetNumVerticesTotal() - nverts;
  std::vector<Vertex> fineLimitPos(nverts);
  std::vector<Vertex> fineDu(nverts);
  std::vector<Vertex> fineDv(nverts);

  primvarRefiner.Limit(&verts[firstOfLastVerts], fineLimitPos, fineDu, fineDv);

  auto vnormals = mesh.add_vertex_property<Vector3f>("v:normal");

  std::vector<surface_mesh::Surface_mesh::Vertex> map;
  map.resize(nverts);

  for (int vert = 0; vert < nverts; ++vert) {

    map[vert] = mesh.add_vertex(Vector3f(fineLimitPos[vert].position[0],
                                         fineLimitPos[vert].position[1],
                                         fineLimitPos[vert].position[2]));

    Vector3f du = Vector3f(fineDu[vert].position[0], fineDu[vert].position[1],
                           fineDu[vert].position[2]);
    Vector3f dv = Vector3f(fineDv[vert].position[0], fineDv[vert].position[1],
                           fineDv[vert].position[2]);

    vnormals[map[vert]] = du.cross(dv).normalized();
  }

  Far::PatchTableFactory::Options ptOptions;
  ptOptions.SetPatchPrecision<real_t>();
  Far::PatchTable const *patchTable =
      Far::PatchTableFactory::Create(*refiner, ptOptions);

  // uv indices for a refined quad face
  real_t uv[4][2] = {
      {0.0, 0.0},
      {1.0, 0.0},
      {1.0, 1.0},
      {0.0, 1.0},
  };

  auto param_loc = mesh.add_halfedge_property<Param_loc>("h:param_loc");

  for (int array = 0; array < patchTable->GetNumPatchArrays(); ++array) {
    for (int patch = 0; patch < patchTable->GetNumPatches(array); ++patch) {

      Far::ConstIndexArray faceVerts =
          patchTable->GetPatchVertices(array, patch);

      auto f = mesh.add_quad(
          map[faceVerts[0] - nCoarseVerts], map[faceVerts[1] - nCoarseVerts],
          map[faceVerts[2] - nCoarseVerts], map[faceVerts[3] - nCoarseVerts]);

      Far::PatchParam param = patchTable->GetPatchParam(array, patch);
      auto h = mesh.halfedge(f);
      for (int quadVert = 0; quadVert < 4; ++quadVert) {
        if (map[faceVerts[0] - nCoarseVerts] == mesh.from_vertex(h))
          break;
        h = mesh.next_halfedge(h);
      }
      assert(map[faceVerts[0] - nCoarseVerts] == mesh.from_vertex(h));
      for (int quadVert = 0; quadVert < 4; ++quadVert) {
        real_t u = uv[quadVert][0];
        real_t v = uv[quadVert][1];
        RefinedToCoarseUV(param, u, v);
        param_loc[h] = Param_loc(param.GetFaceId(), u, v);
        h = mesh.next_halfedge(h);
      }
    }
  }

  // Check parametric locations
  surface_mesh::Surface_mesh::Halfedge_iterator hit;
  for (hit = mesh.halfedges_begin(); hit != mesh.halfedges_end(); ++hit)
    if (!mesh.is_boundary(*hit))
      assert(param_loc[*hit].ptexIndex != -1);
}

void SubdivOsd::adaptiveRefinement() {
  // Adaptively refine the topology
  float max_sharpness = 0.f;
  for (Far::Index i = 0; i < m_refiner->GetLevel(0).GetNumEdges(); i++) {
    float shaprness = m_refiner->GetLevel(0).GetEdgeSharpness(i);
    if (shaprness > max_sharpness)
      max_sharpness = shaprness;
  }
  for (Far::Index i = 0; i < m_refiner->GetLevel(0).GetNumVertices(); i++) {
    float shaprness = m_refiner->GetLevel(0).GetVertexSharpness(i);
    if (shaprness > max_sharpness)
      max_sharpness = shaprness;
  }
  int maxIsolation = int(std::ceil(max_sharpness));

  Far::TopologyRefiner::AdaptiveOptions options(maxIsolation);
  options.considerFVarChannels = !m_shape->uvs.empty();
  options.useInfSharpPatch = false;
  m_refiner->RefineAdaptive(options);

  // Generate a set of Far::PatchTable that we will use to evaluate the surface
  // limit
  Far::PatchTableFactory::Options patchOptions(maxIsolation);
  patchOptions.SetPatchPrecision<real_t>();
  patchOptions.SetFVarPatchPrecision<real_t>();
  patchOptions.SetEndCapType(
      Far::PatchTableFactory::Options::ENDCAP_GREGORY_BASIS);

  patchOptions.useInfSharpPatch = false;
  patchOptions.generateFVarTables = false;
  patchOptions.generateFVarLegacyLinearPatches = false;
  patchOptions.includeFVarBaseLevelIndices = true;

  patchTable = Far::PatchTableFactory::Create(*m_refiner, patchOptions);

  // Compute the total number of points we need to evaluate patchtable.
  // we use local points around extraordinary features.
  int nRefinerVertices = m_refiner->GetNumVerticesTotal();
  int nLocalPoints = patchTable->GetNumLocalPoints();

  // Create a buffer to hold the position of the refined verts and
  // local points, then copy the coarse positions at the beginning.
  verts = std::vector<Vertex>(nRefinerVertices + nLocalPoints);
  for (Far::Index i = 0; i < m_refiner->GetLevel(0).GetNumVertices(); i++) {
    verts[i].SetPosition(m_shape->verts[i * 3], m_shape->verts[i * 3 + 1],
                         m_shape->verts[i * 3 + 2]);
  }

  // Adaptive refinement may result in fewer levels than maxIsolation.
  int nRefinedLevels = m_refiner->GetNumLevels();

  // Interpolate vertex primvar data : they are the control vertices of the
  // limit patches
  Far::PrimvarRefinerReal<real_t> primvarRefiner(*m_refiner);

  Vertex *src = &verts[0];
  for (int level = 1; level < nRefinedLevels; ++level) {
    Vertex *dst = src + m_refiner->GetLevel(level - 1).GetNumVertices();
    primvarRefiner.Interpolate(level, src, dst);
    src = dst;
  }

  // Evaluate local points from interpolated vertex primvars.
  if (nLocalPoints) {
    patchTable->GetLocalPointStencilTable<real_t>()->UpdateValues(
        &verts[0], &verts[nRefinerVertices]);
  }

  // Create a Far::PatchMap to help locating patches in the table
  patchmap = new Far::PatchMap(*patchTable);

  // Create a Far::PtexIndices to help find indices of ptex faces.
  ptexIndices = new Far::PtexIndices(*m_refiner);
}

bool SubdivOsd::evaluateLimit(const Param_loc &p, Vector3f &position,
                              Vector3f &normal) const {
  Vector3f ds, dt;
  this->evaluateLimitFrame(p, position, ds, dt);
  normal = ds.cross(dt).normalized();
  return is_near_extraordinary(p);
}

bool SubdivOsd::evaluateLimitFrame(const Param_loc &p, Vector3f &position,
                                   Vector3f &ds, Vector3f &dt, Vector3f *dsds,
                                   Vector3f *dsdt, Vector3f *dtdt) const {
  // Locate the patch corresponding to the face ptex idx and (s,t)
  Far::PatchTable::PatchHandle const *handle =
      patchmap->FindPatch(p.ptexIndex, p.uv(0), p.uv(1));
  assert(handle);

  real_t pWeights[20], dsWeights[20], dtWeights[20], d2sWeights[20],
      dstWeights[20], d2tWeights[20];
  // Evaluate the patch weights, identify the CVs and compute the limit frame:
  patchTable->EvaluateBasis(*handle, p.uv(0), p.uv(1), pWeights, dsWeights,
                            dtWeights, d2sWeights, dstWeights, d2tWeights);

  Far::ConstIndexArray cvs = patchTable->GetPatchVertices(*handle);

  LimitFrame dst;
  dst.Clear();
  for (int cv = 0; cv < cvs.size(); ++cv) {
    dst.AddWithWeight(verts[cvs[cv]], pWeights[cv], dsWeights[cv],
                      dtWeights[cv], d2sWeights[cv], dstWeights[cv],
                      d2tWeights[cv]);
  }
  position = Vector3f(dst.position[0], dst.position[1], dst.position[2]);

  ds << dst.deriv_s[0], dst.deriv_s[1], dst.deriv_s[2];
  dt << dst.deriv_t[0], dst.deriv_t[1], dst.deriv_t[2];

  if (dsds)
    *dsds << dst.deriv2_s[0], dst.deriv2_s[1], dst.deriv2_s[2];
  if (dsdt)
    *dsdt << dst.deriv_st[0], dst.deriv_st[1], dst.deriv_st[2];
  if (dtdt)
    *dtdt << dst.deriv2_t[0], dst.deriv2_t[1], dst.deriv2_t[2];

  return is_near_extraordinary(p);
}

bool SubdivOsd::is_near_extraordinary(const Param_loc &p) const {
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

void SubdivOsd::determine_extraordinary_vertices(Mesh const &mesh) {
  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");

  if (!param_loc) {
    logger().warn("No parameterization in determine_extraordinary_vertices.");
    return;
  }

  std::map<int, size_t> patch_count;
  for (uint32_t i = 0; i != mesh.n_vertices(); ++i) {
    surface_mesh::Surface_mesh::Vertex v(i);

    std::unordered_map<int, Param_loc> control_patches;
    auto hit = mesh.halfedges(v);
    auto hit_end = hit;
    do {
      control_patches[param_loc[*hit].ptexIndex] = param_loc[*hit];
      control_patches[param_loc[mesh.opposite_halfedge(*hit)].ptexIndex] =
          param_loc[mesh.opposite_halfedge(*hit)];
    } while (++hit != hit_end);

    patch_count[i] = control_patches.size();
  }

  std::unordered_set<int> extraordinary_face_indices;
  std::unordered_set<int> seam_face_indices;

  // Add the between seam neighbors
  for (uint32_t i = 0; i != mesh.n_vertices(); ++i) {
    surface_mesh::Surface_mesh::Vertex v(i);

    // Not extraordinary
    if (!(patch_count[i] != 4 && patch_count[i] > 2))
      continue;

    extraordinary_face_indices.emplace(i);
    auto hit = mesh.halfedges(v);
    auto hit_end = hit;
    do {
      int n_v = mesh.to_vertex(*hit).idx();
      if (patch_count[n_v] == 2) {
        extraordinary_face_indices.emplace(n_v);
        seam_face_indices.emplace(n_v);
      }
    } while (++hit != hit_end);
  }

  // Add the within-patch vertices
  for (uint32_t i = 0; i != mesh.n_vertices(); ++i) {
    surface_mesh::Surface_mesh::Vertex v(i);

    // Not within a patch
    if ((patch_count[i] != 1))
      continue;

    size_t between_seam_neighbors = 0;
    auto hit = mesh.halfedges(v);
    auto hit_end = hit;
    do {
      int n_v = mesh.to_vertex(*hit).idx();
      if (seam_face_indices.find(n_v) != seam_face_indices.end())
        between_seam_neighbors++;
    } while (++hit != hit_end);

    if (between_seam_neighbors == 2) {
      extraordinary_face_indices.emplace(i);
    }
  }

  // Record the extraordinary faces
  for (uint32_t i = 0; i != mesh.n_faces(); ++i) {
    Face f(i);
    std::vector<surface_mesh::Surface_mesh::Vertex> vv;
    auto hit = mesh.halfedges(f);
    auto hit_end = hit;
    do {
      vv.emplace_back(mesh.from_vertex(*hit));
    } while (++hit != hit_end);

    contess_assert_msg(vv.size() == 3, "Only runs on triangule mesh.");

    if (extraordinary_face_indices.find(vv[0].idx()) !=
            extraordinary_face_indices.end() &&
        extraordinary_face_indices.find(vv[1].idx()) !=
            extraordinary_face_indices.end() &&
        extraordinary_face_indices.find(vv[2].idx()) !=
            extraordinary_face_indices.end()) {
      auto hit = mesh.halfedges(f);
      auto hit_end = hit;
      do {
        if (extraordinary_parameters.find(param_loc[*hit].ptexIndex) ==
            extraordinary_parameters.end())
          extraordinary_parameters[param_loc[*hit].ptexIndex] =
              std::vector<Param_loc>();
        extraordinary_parameters[param_loc[*hit].ptexIndex].emplace_back(
            param_loc[*hit]);
      } while (++hit != hit_end);
    }
  }
}
