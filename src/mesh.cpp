// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "mesh.h"

#include <igl/cut_to_disk.h>
#include <igl/triangle/triangulate.h>

#include <limits>
#include <spdlog/fmt/ostr.h>

#include "common.h"
#include "insert_interpolated_contours.h"
#include "insert_planar_map_intersections.h"
#include "logger.h"
#include "shor.h"
#include "svg.h"
#include "sweepLine.h"
#include "tag_concave_edges.h"

using namespace std;

Mesh::Mesh() : m_bvh(nullptr), m_visible(-1), m_subdiv() {}

Mesh::Mesh(const std::string &filename)
    : m_bvh(nullptr), m_visible(-1), m_subdiv() {
  load(filename);
}

Mesh::~Mesh() {
  for (size_t i = 0; i < m_fedges.size(); i++) {
    if (!m_fedges[i] || !m_fedges[i]->m_segment)
      continue;
    delete_pointed_to<FEdge>(m_fedges[i]);
  }
  m_fedges.clear();
  for (size_t i = 0; i < m_svertices.size(); i++) {
    if (!m_svertices[i] || m_svertices[i]->m_fedges.empty() ||
        !m_svertices[i]->m_fedges.front())
      continue;
    delete_pointed_to<SVertex>(m_svertices[i]);
  }
  m_svertices.clear();
}

void Mesh::init() {
  if (!is_triangle_mesh()) {
    triangulate();
    auto param_loc = get_halfedge_property<Param_loc>("h:param_loc");
    if (param_loc) {
      Halfedge_iterator hit;
      for (hit = halfedges_begin(); hit != halfedges_end(); ++hit) {
        if (!is_boundary(*hit) && param_loc[*hit].ptexIndex == -1) {
          // fix parametric location of the new halfedge inserted by the
          // triangulation
          Halfedge oh = opposite_halfedge(*hit);
          param_loc[*hit] = param_loc[next_halfedge(oh)];
          param_loc[oh] = param_loc[next_halfedge(*hit)];
        }
      }
    }
  }

  // Once it's triangulated, we can run the extraordinary determination
  // Store the parameter of the extraordinary vertices
  if (m_subdiv)
    m_subdiv->determine_extraordinary_vertices(*this);

  updateBoundingBox();

  update_face_normals();

  auto vpositions = get_vertex_property<Vector3f>("v:point");
  assert(vpositions);

  auto vnormals = get_vertex_property<Vector3f>("v:normal");
  if (!vnormals) {
    update_vertex_normals();
    vnormals = get_vertex_property<Vector3f>("v:normal");
  }
  assert(vnormals);
  m_normal_segments = Matrix3Xf(3, 2 * n_vertices());
  Vertex_iterator vit;
  int i = 0;
  for (vit = vertices_begin(); vit != vertices_end(); ++vit) {
    Vector3f p = vpositions[*vit];
    Vector3f n = vnormals[*vit];
    m_normal_segments.col(i++) = p;
    m_normal_segments.col(i++) = p + 0.05 * get_dist_max() * n;
  }

  auto colors = get_vertex_property<Vector3f>("v:color");
  if (!colors) {
    colors = add_vertex_property<Vector3f>("v:color");
    this->get_colors().setOnes();
  }

  m_indices.reserve(3 * n_faces());
  // face iterator
  Face_iterator fit, fend = faces_end();
  // vertex circulator
  Vertex vertices[3];
  for (fit = faces_begin(); fit != fend; ++fit) {
    verticesOfFace(*fit, vertices);
    for (Vertex v : vertices)
      m_indices.push_back(v.idx());
  }
  m_sorted_indices = m_indices;
}

bool Mesh::load(const std::string &filename) {
  if (!read(filename))
    return false;

  init();

  return true;
}

void Mesh::updateBoundingBox() {
  m_bbox.clear();
  auto vertices = get_vertex_property<Vector3f>("v:point");
  for (const auto &p : vertices.vector())
    m_bbox.expandBy(p);
}

void Mesh::updateIndicesFromVBO(const Vector3f &view_pos) {
  std::vector<uint32_t> front_indices;
  std::vector<uint32_t> back_indices;
  std::vector<uint32_t> undefined_indices;
  m_indices.clear();
  m_indices.reserve(3 * n_faces());
  m_inconsistent_indices.clear();

  auto VBO = get_face_property<FacingType>("f:VBO");
  auto positions = get_vertex_property<Vector3f>("v:point");
  auto fnormals = get_face_property<Vector3f>("f:normal");

  Vertex v[3];
  for (Face_iterator fit = faces_begin(); fit != faces_end(); ++fit) {
    verticesOfFace(*fit, v);

    switch (VBO[*fit]) {
    case FRONT:
      for (int i = 0; i < 3; i++) {
        front_indices.push_back(v[i].idx());
        m_indices.push_back(v[i].idx());
      }
      break;
    case BACK:
      for (int i = 0; i < 3; i++) {
        back_indices.push_back(v[i].idx());
        m_indices.push_back(v[i].idx());
      }
      break;
    case UNDEFINED:
      for (int i = 0; i < 3; i++) {
        undefined_indices.push_back(v[i].idx());
        m_indices.push_back(v[i].idx());
      }
      break;
    default:
      break;
    }

    // compute consistency
    real_t ndotv =
        fnormals[*fit].dot((view_pos - positions[v[0]]).normalized());
    if ((ndotv > 0 && VBO[*fit] != FRONT) || (ndotv < 0 && VBO[*fit] != BACK)) {
      for (int i = 0; i < 3; i++)
        m_inconsistent_indices.push_back(v[i].idx());
    }
  }
  m_sorted_indices.clear();
  m_sorted_indices.insert(m_sorted_indices.end(), front_indices.begin(),
                          front_indices.end());
  m_sorted_indices.insert(m_sorted_indices.end(), back_indices.begin(),
                          back_indices.end());
  m_sorted_indices.insert(m_sorted_indices.end(), undefined_indices.begin(),
                          undefined_indices.end());
  m_front = front_indices.size();
  m_back = back_indices.size();
}

void Mesh::tagConcaveEdges(ContourMode cont_mode, ConvexityMode mode,
                           Camera const &camera) {
  cont_mode = ContourMode::VBO_CONTOUR;

  switch (mode) {
  case ConvexityMode::MESH_CONVEXITY:
    tag_concave_edges_mesh(*this, cont_mode);
    break;
  case ConvexityMode::CAMERA_CONVEXITY:
    tag_concave_edges_side(*this, camera, cont_mode);
    break;
  default:
    break;
  }
}

inline FacingType getFacing(const Vector3f &viewVec, const Vector3f &normal,
                            real_t *n_dot_v = NULL) {
  // view vec is (position - cameraCenter). always normalized.
  // only produces a normalized result if the input normal is normalized
  assert(viewVec.norm() > 0);

  Vector3f vv = viewVec.normalized();

  real_t ndotv = normal.dot(vv);

  if (n_dot_v != NULL)
    *n_dot_v = ndotv;

  if (ndotv < -CONTOUR_THRESHOLD)
    return BACK;

  if (ndotv > CONTOUR_THRESHOLD)
    return FRONT;

  return CONTOUR;
}

void Mesh::update_patch_facing(const Mesh &mesh) {
  m_patch_facing = mesh.m_patch_facing;
}

FacingType Mesh::get_patch_facing(size_t patch_index) const {
  assert(patch_index < m_patch_facing.size() && "Patch index out of range.");
  return m_patch_facing[patch_index];
}

void Mesh::computeConsistency(const Vector3f &view_pos) {
  // Update ndotv & vertex facing
  auto ndotv = vertex_property<real_t>("v:ndotv");
  auto facing = vertex_property<FacingType>("v:facing");
  auto positions = get_vertex_property<Vector3f>("v:point");
  assert(positions);
  auto vnormals = get_vertex_property<Vector3f>("v:normal");
  assert(vnormals);

  parallel_for(
      tbb::blocked_range<uint32_t>(0u, (uint32_t)n_vertices(), GRAIN_SIZE),
      [&](const tbb::blocked_range<uint32_t> &range) {
        for (uint32_t i = range.begin(); i != range.end(); ++i) {
          Vertex v = Vertex(i);
          facing[v] =
              getFacing(view_pos - positions[v], vnormals[v], &ndotv[v]);
        }
      });

  // Update face Vertex-Based Orientation and sort indices accordingly
  auto VBO = get_face_property<FacingType>("f:VBO");
  if (!VBO)
    VBO = add_face_property<FacingType>("f:VBO");

  parallel_for(
      tbb::blocked_range<uint32_t>(0u, (uint32_t)n_faces(), GRAIN_SIZE),
      [&](const tbb::blocked_range<uint32_t> &range) {
        Vertex v[3];
        for (uint32_t i = range.begin(); i != range.end(); ++i) {
          Face f = Face(i);
          verticesOfFace(f, v);

          FacingType result = UNDEFINED;
          for (int i = 0; i < 3; i++) {
            if (facing[v[i]] == UNDEFINED)
              continue;

            if (result != UNDEFINED && facing[v[i]] != result) {
              result = UNDEFINED;
              break;
            }

            result = facing[v[i]];
          }
          VBO[f] = result;
        }
      });

  // Update the VBO_f (the facing based on face normal)
  auto VBO_f = face_property<FacingType>("f:VBO_f");
  if (!VBO_f) {
    VBO_f = add_face_property<FacingType>("f:VBO_f", FacingType::NA);
  }
  VBO_f.vector().assign(VBO_f.vector().size(), FacingType::NA);

  parallel_for(
      tbb::blocked_range<uint32_t>(0u, (uint32_t)n_faces(), GRAIN_SIZE),
      [&](const tbb::blocked_range<uint32_t> &range) {
        Vertex vs[3];
        for (uint32_t i = range.begin(); i != range.end(); ++i) {
          Face f = Face(i);
          // The facing type based on the face normal
          auto f_norm = compute_face_normal(f);

          verticesOfFace(f, vs);
          auto c_ray = positions[vs[0]] - view_pos;
          c_ray.normalized();
          VBO_f[f] =
              (f_norm.dot(c_ray) < 0) ? FacingType::FRONT : FacingType::BACK;
          // Use predicate to be mroe accurate
          VBO_f[f] = (igl::predicates::orient3d(
                          positions[vs[0]], positions[vs[1]], positions[vs[2]],
                          view_pos) == igl::predicates::Orientation::NEGATIVE)
                         ? FacingType::FRONT
                         : FacingType::BACK;
        }
      });
  updateIndicesFromVBO(view_pos);
}

void Mesh::computeConsistencySubdivided(const Vector3f &view_pos) {
  // Update ndotv & vertex facing
  auto ndotv = vertex_property<real_t>("v:ndotv");
  auto facing = vertex_property<FacingType>("v:facing");
  auto positions = vertex_property<Vector3f>("v:point");
  contess_assert(positions);
  auto vnormals = vertex_property<Vector3f>("v:normal");
  contess_assert(vnormals);
  auto param_loc = get_halfedge_property<Param_loc>("h:param_loc");
  contess_assert(param_loc);

  {
    for (uint32_t i = 0; i != n_vertices(); ++i) {
      Vertex v = Vertex(i);

      // Compute the position and normal from the subdivision
      Halfedge he;
      auto hit = halfedges(v);
      auto hit_end = hit;
      do {
        if (is_boundary(*hit))
          continue;
        he = *hit;
        break;
      } while (++hit != hit_end);
      contess_assert(he.is_valid());

      const Param_loc &v_param_loc = param_loc[he];
      Vector3f pos_res, normal_res;
      m_subdiv->evaluateLimit(v_param_loc, pos_res, normal_res);
      vnormals[v] = normal_res;
      positions[v] = pos_res;

      facing[v] = getFacing(view_pos - pos_res, vnormals[v], &ndotv[v]);
    }
  }

  // Update face Vertex-Based Orientation and sort indices accordingly
  auto VBO = get_face_property<FacingType>("f:VBO");
  if (!VBO)
    VBO = add_face_property<FacingType>("f:VBO");

  parallel_for(
      tbb::blocked_range<uint32_t>(0u, (uint32_t)n_faces(), GRAIN_SIZE),
      [&](const tbb::blocked_range<uint32_t> &range) {
        Vertex v[3];
        for (uint32_t i = range.begin(); i != range.end(); ++i) {
          Face f = Face(i);
          verticesOfFace(f, v);

          FacingType result = UNDEFINED;
          for (int i = 0; i < 3; i++) {
            if (facing[v[i]] == UNDEFINED)
              continue;

            if (result != UNDEFINED && facing[v[i]] != result) {
              result = UNDEFINED;
              break;
            }

            result = facing[v[i]];
          }
          VBO[f] = result;
        }
      });

  // Update the VBO_f (the facing based on face normal)
  auto VBO_f = face_property<FacingType>("f:VBO_f");
  if (!VBO_f) {
    VBO_f = add_face_property<FacingType>("f:VBO_f", FacingType::NA);
  }
  VBO_f.vector().assign(VBO_f.vector().size(), FacingType::NA);

  parallel_for(
      tbb::blocked_range<uint32_t>(0u, (uint32_t)n_faces(), GRAIN_SIZE),
      [&](const tbb::blocked_range<uint32_t> &range) {
        Vertex vs[3];
        for (uint32_t i = range.begin(); i != range.end(); ++i) {
          Face f = Face(i);
          // The facing type based on the face normal
          auto f_norm = compute_face_normal(f);

          verticesOfFace(f, vs);
          auto c_ray = positions[vs[0]] - view_pos;
          c_ray.normalized();
          VBO_f[f] =
              (f_norm.dot(c_ray) < 0) ? FacingType::FRONT : FacingType::BACK;
          // Use predicate to be mroe accurate
          VBO_f[f] = (igl::predicates::orient3d(
                          positions[vs[0]], positions[vs[1]], positions[vs[2]],
                          view_pos) == igl::predicates::Orientation::NEGATIVE)
                         ? FacingType::FRONT
                         : FacingType::BACK;
        }
      });

  updateIndicesFromVBO(view_pos);
}

inline real_t find_zero_linear(real_t val0, real_t val1) {
  return val0 / (val0 - val1);
}

bool Mesh::is_contour_edge(Edge const &e) const {
  auto is_contour = get_edge_property<real_t>("e:contour");

  if (!is_contour) {
    // Warning here
    std::cerr
        << "Warning: Calling contour edge check before extracting contours"
        << std::endl;
    return false;
  }

  return is_contour[e] >= 0;
}

template <typename T> inline int sign(T val) {
  return (T(0) < val) - (val < T(0));
}

void Mesh::extractContours(ContourMode mode, Camera const &camera) {
  Vector3f view_pos = camera.position();
  auto vpositions = get_vertex_property<Vector3f>("v:point");
  auto is_contour = get_edge_property<real_t>("e:contour");
  if (!is_contour) {
    is_contour = add_edge_property<real_t>("e:contour", -1.f);
  } else {
    for (size_t i = 0; i < is_contour.vector().size(); ++i)
      is_contour.vector()[i] = -1.f;
  }

  m_contour_edges.clear();
  m_contour_edges.reserve(std::sqrt(n_faces()));

  auto fnormals = get_face_property<Vector3f>("f:normal");
  auto ndotv = get_vertex_property<real_t>("v:ndotv");
  auto VBO = get_face_property<FacingType>("f:VBO");

  parallel_for(
      tbb::blocked_range<uint32_t>(0u, (uint32_t)n_edges(), GRAIN_SIZE),
      [&](const tbb::blocked_range<uint32_t> &range) {
        for (uint32_t i = range.begin(); i != range.end(); ++i) {
          Edge e = Edge(i);

          if (mode == INTERPOLATED_CONTOUR) {
            Vertex i0 = vertex(e, 0);
            Vertex i1 = vertex(e, 1);
            real_t d0 = ndotv[i0];
            real_t d1 = ndotv[i1];
            if ((d0 > 0 && d1 <= 0) || (d0 < 0 && d1 >= 0)) {
              real_t w = find_zero_linear(d0, d1);
              m_contour_edges.push_back(e);
              is_contour[e] = w;
            }
          } else {
            if (is_boundary(e))
              continue;
            Face f1 = face(e, 0);
            Face f2 = face(e, 1);
            const Vector3f &normal1 = fnormals[f1];
            const Vector3f &normal2 = fnormals[f2];
            const Vector3f &v1 = vpositions[vertex(e, 0)];
            const Vector3f &v2 = vpositions[vertex(e, 1)];

            if ((mode == MESH_CONTOUR &&
                 sign((view_pos - v1).dot(normal1)) !=
                     sign((view_pos - v2).dot(normal2))) ||
                (mode == VBO_CONTOUR && VBO[f1] != VBO[f2])) {
              // silhouette edge found
              m_contour_edges.push_back(e);
              is_contour[e] = 1.f;
            }
          }
        }
      });

  if (mode == MESH_CONTOUR || mode == VBO_CONTOUR) {
    m_mesh_contours_indices = MatrixXu(2, m_contour_edges.size());
    for (size_t i = 0; i < m_contour_edges.size(); ++i) {
      m_mesh_contours_indices.col(i) =
          Vector2u(vertex(m_contour_edges.at(i), 0).idx(),
                   vertex(m_contour_edges.at(i), 1).idx());
    }
  } else {
    tbb::concurrent_vector<std::tuple<Vector3f, Vector3f, int>> contour_pairs;

    parallel_for(
        tbb::blocked_range<uint32_t>(0u, (uint32_t)n_faces(), GRAIN_SIZE),
        [&](const tbb::blocked_range<uint32_t> &range) {
          Halfedge_around_face_circulator hit, hit_end;
          for (uint32_t i = range.begin(); i != range.end(); ++i) {
            Face f = Face(i);
            hit = hit_end = halfedges(f);
            Edge e[2];
            real_t w[2];
            int j = 0;
            do {
              if (is_contour[edge(*hit)] >= 0) {
                e[j] = edge(*hit);
                w[j] = is_contour[edge(*hit)];
                j++;
              }
            } while (++hit != hit_end && j < 2);
            assert(e[0].is_valid() == e[1].is_valid());
            if (!e[0].is_valid() || !e[1].is_valid())
              continue;
            Vector3f vm[2];
            for (j = 0; j < 2; j++) {
              Vertex i0 = vertex(e[j], 0);
              Vertex i1 = vertex(e[j], 1);
              const Vector3f &v0 = vpositions[i0];
              const Vector3f &v1 = vpositions[i1];
              vm[j] = w[j] * v1 + (1.f - w[j]) * v0;
            }
            contour_pairs.push_back(std::make_tuple(vm[0], vm[1], i));
          }
        });

    m_interpolated_contours = Matrix3Xf(3, contour_pairs.size() * 2);
    m_interpolated_contour_faces.resize(contour_pairs.size());
    for (size_t i = 0; i < contour_pairs.size(); ++i) {
      m_interpolated_contours.col(i * 2) = std::get<0>(contour_pairs.at(i));
      m_interpolated_contours.col(i * 2 + 1) = std::get<1>(contour_pairs.at(i));
      m_interpolated_contour_faces[i] = std::get<2>(contour_pairs.at(i));
    }
    // Since m_interpolated_contours changes size and will be no longer
    // corresponding to m_interpolated_contour_faces later...
    // For now, just simply keep a copy of the original m_interpolated_contours
    m_interpolated_contours_orig = m_interpolated_contours;
  }
}

void Mesh::extractBoundaries() {
  m_boundary_indices.clear();
  auto visited = get_edge_property<real_t>("e:visited");
  if (!visited) {
    visited = add_edge_property<real_t>("e:visited", false);
  } else {
    for (size_t i = 0; i < visited.vector().size(); ++i)
      visited.vector()[i] = false;
  }
  auto boundary_fedges = edge_property<FEdge *>("e:fedge", nullptr);
  auto svertices = vertex_property<SVertex *>("v:svertex", nullptr);
  auto vpositions = get_vertex_property<Vector3f>("v:point");

  Halfedge_around_vertex_circulator hit, hit_end;
  int size = 0;
  for (Edge e : edges()) {
    if (is_boundary(e) && !visited[e]) {
      Halfedge sh = halfedge(e, 0);
      if (is_boundary(sh))
        sh = halfedge(e, 1);
      m_boundary_indices.push_back(from_vertex(sh).idx());
      m_boundary_indices.push_back(from_vertex(sh).idx());
      const Vector3f &p = vpositions[from_vertex(sh)];
      SVertex *sv_first, *sv0;
      sv_first = new SVertex(p);
      sv0 = sv_first;
      assert(svertices[from_vertex(sh)] == nullptr);
      svertices[from_vertex(sh)] = sv0;
      m_svertices.push_back(sv0);
      visited[e] = true;
      Halfedge h = sh;
      do {
        hit = hit_end = halfedges(to_vertex(h));
        do {
          if (edge(h) != edge(*hit) && is_boundary(edge(*hit))) {
            h = *hit;
            break;
          }
        } while (++hit != hit_end);
        assert(hit != hit_end);
        visited[edge(h)] = true;
        SVertex *sv1;
        if (svertices[from_vertex(h)]) {
          sv1 = svertices[from_vertex(h)];
        } else {
          const Vector3f &p = vpositions[from_vertex(h)];
          sv1 = new SVertex(p);
          svertices[from_vertex(h)] = sv1;
          m_svertices.push_back(sv1);
        }
        m_fedges.push_back(
            new FEdgeSharp(sv0, sv1, edge(h), EdgeNature::BOUNDARY));
        boundary_fedges[edge(h)] = m_fedges.back();
        sv0->m_fedges.push_back(m_fedges.back());
        sv1->m_fedges.push_back(m_fedges.back());
        sv0 = sv1;
        m_boundary_indices.push_back(from_vertex(h).idx());
      } while (edge(h) != edge(sh));
      m_boundary_indices.push_back(from_vertex(h).idx());
      m_boundaries_lengths.push_back(m_boundary_indices.size() - size);
      size = m_boundary_indices.size();
    }
  }
}

void Mesh::extractBoundaryCurtainFolds(const Vector3f &c) {
  auto vpositions = get_vertex_property<Vector3f>("v:point");
  auto svertices = get_vertex_property<SVertex *>("v:svertex");
  assert(svertices);

  parallel_for(
      tbb::blocked_range<uint32_t>(0u, (uint32_t)m_fedges.size(), GRAIN_SIZE),
      [&](const tbb::blocked_range<uint32_t> &range) {
        for (uint32_t i = range.begin(); i != range.end(); ++i) {
          FEdge *fe = m_fedges[i];
          if (!fe || fe->m_nature != EdgeNature::BOUNDARY)
            continue;
          FEdgeSharp *fs = dynamic_cast<FEdgeSharp *>(fe);
          assert(fs);
          Edge edge = fs->m_edge;
          for (int i = 0; i < 2; ++i) {
            Vertex ie = vertex(edge, i);
            const Vector3f &e = vpositions[ie];
            Vertex ip = vertex(edge, (i + 1) % 2);
            const Vector3f &p = vpositions[ip];
            Halfedge_around_vertex_circulator hit, hit_end;
            hit = hit_end = halfedges(ip);
            do {
              if (is_boundary(*hit))
                continue;
              Vertex iq = to_vertex(*hit);
              Vertex ir = to_vertex(next_halfedge(*hit));
              if (iq == ie || ir == ie || iq == ip || ir == ip)
                continue;
              const Vector3f &q = vpositions[iq];
              const Vector3f &r = vpositions[ir];
              if (!sameSide(p.cast<double>(), q.cast<double>(),
                            r.cast<double>(), c.cast<double>(),
                            e.cast<double>()) &&
                  sameSide(c.cast<double>(), p.cast<double>(), q.cast<double>(),
                           e.cast<double>(), r.cast<double>()) &&
                  sameSide(c.cast<double>(), p.cast<double>(), r.cast<double>(),
                           e.cast<double>(), q.cast<double>())) {
                // boundary curtain fold found
                SVertex *sv = svertices[ip];
                assert(sv);
                sv->m_nature =
                    VertexNature::BOUNDARY_CURTAIN_FOLD | sv->m_nature;
              }
            } while (++hit != hit_end);
          }
        }
      });
}

extern int tri_tri_intersection_test_3d(double p1[3], double q1[3],
                                        double r1[3], double p2[3],
                                        double q2[3], double r2[3],
                                        int *coplanar, double source[3],
                                        double target[3]);

void Mesh::verticesOfFace(Face f, Vertex v[3]) const {
  Vertex_around_face_circulator fvit, fvend;
  fvit = fvend = vertices(f);
  v[0] = *fvit;
  ++fvit;
  v[2] = *fvit;
  do {
    v[1] = v[2];
    ++fvit;
    v[2] = *fvit;
  } while (++fvit != fvend);
}

void Mesh::halfedgesOfFace(Face f, Halfedge h[3]) const {
  Halfedge_around_face_circulator fvit, fvend;
  fvit = fvend = halfedges(f);
  h[0] = *fvit;
  ++fvit;
  h[2] = *fvit;
  do {
    h[1] = h[2];
    ++fvit;
    h[2] = *fvit;
  } while (++fvit != fvend);
}

bool Mesh::are_adjacent(Face f1, Face f2) const {
  Vertex_around_face_circulator vit, vit_end;
  vit = vit_end = vertices(f2);
  do {
    Face_around_vertex_circulator fit, fit_end;
    fit = fit_end = faces(*vit);
    do {
      if (*fit == f1)
        return true;
    } while (++fit != fit_end);
  } while (++vit != vit_end);
  return false;
};

void Mesh::extractSurfaceIntersections() {

  if (m_bvh)
    delete m_bvh;
  m_bvh = new BVH(
      &m_indices, &get_vertex_property<Vector3f>("v:point").vector(),
      &get_vertex_property<Vector3f>("v:normal").vector(), boundingBox());
  m_bvh->build(nullptr);

  auto vpositions = get_vertex_property<Vector3f>("v:point");
  struct SurfaceIntersection {
    SurfaceIntersection(real_t _t, Edge _e, const Vector3f &_p)
        : t(_t), e(_e), visited(false) {
      neighboors.reserve(2);
      sv = new SVertex(_p);
    }
    real_t t;
    Edge e;
    bool visited;
    std::vector<SurfaceIntersection *> neighboors;
    SVertex *sv;
  };
  auto surf_intersect =
      add_edge_property<std::vector<SurfaceIntersection *>>("e:intersection");
  auto boundary_fedges = edge_property<FEdge *>("e:fedge", nullptr);
  auto fedges_f = face_property<std::vector<FEdgeIntersection *>>("f:fedge");
  auto svertices = edge_property<std::vector<SVertex *>>("e:svertex");

  std::set<Edge> intersected_edges;

  auto findSurfaceIntersection = [&](Edge e,
                                     real_t t) -> SurfaceIntersection * {
    SurfaceIntersection *s = nullptr;
    for (size_t j = 0; j < surf_intersect[e].size(); j++) {
      if (std::abs(surf_intersect[e][j]->t - t) < EPSILON) {
        // same intersection point
        s = surf_intersect[e][j];
        break;
      }
    }
    if (!s) {
      const Vector3f &v0 = vpositions[vertex(e, 0)];
      const Vector3f &v1 = vpositions[vertex(e, 1)];
      Vector3f p = v0 + t * (v1 - v0).normalized();
      s = new SurfaceIntersection(t, e, p);
      surf_intersect[e].push_back(s);
      svertices[e].push_back(s->sv);
      m_svertices.push_back(s->sv);

      if (boundary_fedges[e]) { // 3D intersection with a boundary edge
        s->sv->m_nature = VertexNature::INTERSECTION_3D | s->sv->m_nature;
      }
    }
    return s;
  };

  for (Face f1 : faces()) {
    Vertex v1[3];
    verticesOfFace(f1, v1);
    Eigen::Vector3d p1 = vpositions[v1[0]].cast<double>();
    Eigen::Vector3d q1 = vpositions[v1[1]].cast<double>();
    Eigen::Vector3d r1 = vpositions[v1[2]].cast<double>();

    auto intersect = [&](const Eigen::Vector3d &o, const Eigen::Vector3d &d) {
      Ray r = Ray(o.cast<real_t>(), d.normalized().cast<real_t>(), 0, d.norm());
      std::vector<uint32_t> indices;
      if (m_bvh->rayIntersect(r, indices)) {
        for (int idx : indices) {
          Face f2 = Face(idx);

          if (f1 == f2 || are_adjacent(f1, f2))
            continue;
          bool found = false;
          for (FEdgeIntersection *fe : fedges_f[f1]) {
            if (fe->m_f1 == f2 || fe->m_f2 == f2) {
              found = true;
              break;
            }
          }
          if (found) // Face already processed
            continue;

          Vertex v2[3];
          verticesOfFace(f2, v2);
          Eigen::Vector3d p2 = vpositions[v2[0]].cast<double>();
          Eigen::Vector3d q2 = vpositions[v2[1]].cast<double>();
          Eigen::Vector3d r2 = vpositions[v2[2]].cast<double>();
          Eigen::Vector3d i1, i2;
          int colplanar;
          int res = tri_tri_intersection_test_3d(
              p1.data(), q1.data(), r1.data(), p2.data(), q2.data(), r2.data(),
              &colplanar, i1.data(), i2.data());

          if (res == 0) // no intersection
            continue;

          // figure out which intersection goes with which face edge
          Face faces[2] = {f1, f2};
          real_t min_distA, min_distB;
          min_distA = min_distB = std::numeric_limits<real_t>::infinity();
          Edge eA, eB;
          real_t tA, tB;
          for (int i = 0; i < 2; i++) {
            Halfedge_around_face_circulator hit, hit_end;
            hit = hit_end = halfedges(faces[i]);
            do {
              Edge e = edge(*hit);
              const Vector3f &v0 = vpositions[vertex(e, 0)];
              const Vector3f &v1 = vpositions[vertex(e, 1)];
              Vector3f dir = (v1 - v0).normalized();
              Vector3f e1 = v0 - i1.cast<real_t>();
              real_t dist = dir.cross(e1).squaredNorm();
              if (dist <= min_distA) {
                min_distA = dist;
                tA = -dir.dot(e1);
                eA = e;
              }
              Vector3f e2 = v0 - i2.cast<real_t>();
              dist = dir.cross(e2).squaredNorm();
              if (dist <= min_distB) {
                min_distB = dist;
                tB = -dir.dot(e2);
                eB = e;
              }
            } while (++hit != hit_end);
          }

          if (!eA.is_valid() || !eB.is_valid())
            throw std::runtime_error(
                "Cannot find matching surface-surface intersections.");

          SurfaceIntersection *sA = findSurfaceIntersection(eA, tA);
          SurfaceIntersection *sB = findSurfaceIntersection(eB, tB);
          sA->neighboors.push_back(sB);
          sB->neighboors.push_back(sA);

          intersected_edges.insert(eA);
          intersected_edges.insert(eB);

          FEdgeIntersection *fe = new FEdgeIntersection(sA->sv, sB->sv, f1, f2);
          m_fedges.push_back(fe);
          sA->sv->m_fedges.push_back(fe);
          sB->sv->m_fedges.push_back(fe);
          fedges_f[f1].push_back(fe);
          fedges_f[f2].push_back(fe);
        }
      }
    };

    intersect(p1, q1 - p1);
    intersect(p1, r1 - p1);
    intersect(q1, r1 - q1);
  }

  // Chaining
  std::vector<std::list<SVertex *>> intersections_chains;
  int sum = 0;
  for (auto it = intersected_edges.begin(); it != intersected_edges.end();
       it++) {
    Edge e = *it;
    for (SurfaceIntersection *s : surf_intersect[e]) {
      if (s->visited)
        continue;
      SurfaceIntersection *is = s;
      std::list<SVertex *> chain;
      bool forward = true;
      while (s) {
        s->visited = true;

        if (forward)
          chain.push_back(s->sv);
        else
          chain.push_front(s->sv);

        bool found = false;
        bool loop_found = false;
        for (auto ns : s->neighboors) {
          if (!ns->visited) {
            s = ns;
            loop_found = false;
            found = true;
            break;
          }
          if (ns == is)
            loop_found = true;
        }
        if (loop_found) {
          chain.push_front(s->sv);
          break;
        }
        if (!found) {
          if (forward) {
            for (auto ns : is->neighboors) {
              if (!ns->visited) {
                s = ns;
                forward = false;
                found = true;
                break;
              }
            }
          }
          if (!found)
            s = nullptr;
        }
      }
      chain.push_back(chain.back());
      chain.push_front(chain.front());
      intersections_chains.push_back(chain);
      m_surface_intersections_lengths.push_back(chain.size());
      sum += chain.size();
    }
  }

  m_surface_intersections = Matrix3Xf(3, sum);
  int i = 0;
  for (auto chain : intersections_chains) {
    auto it = chain.begin();
    SVertex *sv0 = *it;
    m_surface_intersections.col(i++) = sv0->m_pos3D;
    it++;
    while (it != chain.end()) {
      SVertex *sv1 = *it;
      m_surface_intersections.col(i++) = sv1->m_pos3D;
      sv0 = sv1;
      it++;
    }
  }

  for (auto vec : surf_intersect.vector())
    for (auto s : vec)
      delete s;
  remove_edge_property(surf_intersect);
}

void Mesh::projectSVertices(const Matrix4f &model, const Matrix4f &proj,
                            const Vector2i &viewportSize) {
  for (auto v : m_svertices) {
    v->m_pos2D = project(v->m_pos3D, model, proj, viewportSize);
  }
}

int Mesh::splitEdgesWithZeroCrossing(const Vector3f &view_pos,
                                     const Subdiv *subdiv) {
  auto facing = get_vertex_property<FacingType>("v:facing");
  auto param_loc = get_halfedge_property<Param_loc>("h:param_loc");

  tbb::concurrent_vector<EdgeTuple> edges_to_split;

  for (uint32_t i = 0; i != n_edges(); ++i) {
    Edge e = Edge(i);
    // 1. Check if the zero-crossing exists on the current edge
    // Only consider FF / BB edges
    if ((facing[vertex(e, 1)] != facing[vertex(e, 0)] &&
         facing[vertex(e, 1)] != CONTOUR && facing[vertex(e, 0)] != CONTOUR) ||
        (facing[vertex(e, 1)] == facing[vertex(e, 0)] &&
         facing[vertex(e, 0)] == CONTOUR))
      continue;

    auto he = halfedge(e, 0);
    if (is_boundary(he))
      he = halfedge(e, 1);
    Param_loc param_res;
    const Param_loc &v1_param_loc = param_loc[he];
    const Param_loc &v2_param_loc = param_loc[next_halfedge(he)];
    assert(v1_param_loc.ptexIndex == v2_param_loc.ptexIndex);
    param_res.ptexIndex = v1_param_loc.ptexIndex;

    // If one is extraordinary, the other is not, switching between
    // evaluation schemes may cause bad-shaped contour.
    if (subdiv->is_near_extraordinary(v1_param_loc) !=
        subdiv->is_near_extraordinary(v2_param_loc))
      continue;

    if (subdiv->backend_type() == Subdiv::Backend::LACEWELL &&
        (subdiv->is_near_extraordinary(v1_param_loc) ||
         subdiv->is_near_extraordinary(v2_param_loc)))
      continue;

    // 2. Move from v0 to v1 in NUM_SAMPLES even steps;
    // Check using the limit surface, if the facing of the sample
    // goes from ==v0 to ==v1.
    Vector3f pos_res, normal_res;
    const int NUM_SAMPLES = 10;
    EdgeTuple to_split_tuple(Edge(-1), pos_res, normal_res, param_res, 0);
    float inconsistent_count = 0;

    for (int j = 0; j < NUM_SAMPLES; j++) {
      real_t t = real_t(j + 1) / real_t(NUM_SAMPLES + 1);
      param_res.uv = (1.f - t) * v1_param_loc.uv + t * v2_param_loc.uv;
      subdiv->evaluateLimit(param_res, pos_res, normal_res);

      FacingType ft = getFacing(view_pos - pos_res, normal_res);
      if (ft != facing[vertex(e, 1)] && ft != FacingType::CONTOUR) {
        if (!std::get<0>(to_split_tuple).is_valid())
          to_split_tuple = EdgeTuple(e, pos_res, normal_res, param_res, t);

        inconsistent_count++;
      }
    }

    if (std::get<0>(to_split_tuple).is_valid() &&
        (inconsistent_count / NUM_SAMPLES > MAX_SHIFT_PERCENTAGE)) {
      edges_to_split.push_back(to_split_tuple);
    }
  }

  // Actually split FF or BB edges with zero-crossings
  for (auto etuple : edges_to_split) {
    Edge e = std::get<0>(etuple);
    FacingType n_facing = (facing[vertex(e, 0)] != FacingType::CONTOUR)
                              ? facing[vertex(e, 0)]
                              : facing[vertex(e, 1)];
    FacingType v_facing =
        (n_facing == FacingType::FRONT) ? FacingType::BACK : FacingType::FRONT;
    Vertex new_v = splitEdge(etuple, false);
    facing[new_v] = v_facing;
  }

  garbage_collection();

#ifndef NDEBUG
  // Assert validity of parametric locations
  Halfedge_iterator hit;
  for (hit = halfedges_begin(); hit != halfedges_end(); ++hit)
    if (!is_boundary(*hit))
      assert(param_loc[*hit].ptexIndex != -1);
#endif

  return edges_to_split.size();
}

real_t Mesh::bisect_search(const Vector3f &view_pos, const Subdiv *subdiv,
                           Param_loc const &v1_param_loc,
                           Param_loc const &v2_param_loc, Vector3f &pos_res,
                           Vector3f &normal_res, Param_loc &param_res) {
  // Search along the parameter direction for zero crossing
  assert(v1_param_loc.ptexIndex == v2_param_loc.ptexIndex);
  param_res.ptexIndex = v1_param_loc.ptexIndex;

  Vector3f pos1_res, normal1_res, pos2_res, normal2_res;
  subdiv->evaluateLimit(v1_param_loc, pos1_res, normal1_res);
  subdiv->evaluateLimit(v2_param_loc, pos2_res, normal2_res);

  real_t v1_facing = (view_pos - pos1_res).normalized().dot(normal1_res);
  real_t v2_facing = (view_pos - pos2_res).normalized().dot(normal2_res);
  assert((v1_facing < 0) != (v2_facing < 0));
  real_t res_facing;
  real_t a = 0.f;
  real_t b = 1.f;
  real_t c = 0.f;

  for (int i = 0; i < MAX_ROOT_ITERATIONS; i++) {
    c = b - (b - a) / (v2_facing - v1_facing) * v2_facing;
    param_res.uv = (1.f - c) * v1_param_loc.uv + c * v2_param_loc.uv;
    subdiv->evaluateLimit(param_res, pos_res, normal_res);

    res_facing = (view_pos - pos_res).normalized().dot(normal_res);
    if (std::abs(res_facing) < CONTOUR_THRESHOLD) {
      // We got lucky and found the zero crossing
      return c;
    }

    if ((res_facing < 0) == (v1_facing < 0)) {
      a = c;
      v1_facing = res_facing;
    } else {
      b = c;
      v2_facing = res_facing;
    }
  }

  // Doesn't converge
  logger().warn("Root finding doesn't converge {}: {}, {} - {}, {}", res_facing,
                v1_param_loc.ptexIndex, v1_param_loc.uv, v2_param_loc.ptexIndex,
                v2_param_loc.uv);

  return c;
}

real_t Mesh::bisect_search(const Vector3f &view_pos, const Subdiv *subdiv,
                           Edge e, Vector3f &pos_res, Vector3f &normal_res,
                           Param_loc &param_res) {
  // Search edge for zero crossing
  auto ndotv = get_vertex_property<real_t>("v:ndotv");
  auto param_loc = get_halfedge_property<Param_loc>("h:param_loc");

  auto he = halfedge(e, 0);
  auto v2_he = he;
  if (is_boundary(he)) {
    v2_he = opposite_halfedge(he);
    he = next_halfedge(opposite_halfedge(he));
  } else {
    v2_he = next_halfedge(he);
  }
  const Param_loc &v1_param_loc = param_loc[he];
  const Param_loc &v2_param_loc = param_loc[v2_he];
  assert(v1_param_loc.ptexIndex == v2_param_loc.ptexIndex);
  param_res.ptexIndex = v1_param_loc.ptexIndex;
  real_t v1_facing = ndotv[vertex(e, 1)];
  real_t v2_facing = ndotv[vertex(e, 0)];
  assert((v1_facing < 0) != (v2_facing < 0));
  real_t res_facing;
  real_t a = 0.f;
  real_t b = 1.f;
  real_t c = 0.f;

  for (int i = 0; i < MAX_ROOT_ITERATIONS; i++) {
    c = b - (b - a) / (v2_facing - v1_facing) * v2_facing;
    param_res.uv = (1.f - c) * v1_param_loc.uv + c * v2_param_loc.uv;
    subdiv->evaluateLimit(param_res, pos_res, normal_res);

    res_facing = (view_pos - pos_res).normalized().dot(normal_res);

    if (std::abs(res_facing) < CONTOUR_THRESHOLD) {
      // We got lucky and found the zero crossing
      return c;
    }

    if ((res_facing < 0) == (v1_facing < 0)) {
      a = c;
      v1_facing = res_facing;
    } else {
      b = c;
      v2_facing = res_facing;
    }
  }

  // Doesn't converge
  logger().warn("Root finding doesn't converge {}: {} - {}", res_facing,
                std::min(vertex(e, 0), vertex(e, 1)),
                std::max(vertex(e, 0), vertex(e, 1)));

  return c;
}

bool Mesh::canShift(Vertex vsource, const EdgeTuple &etuple) {
  auto param_loc = get_halfedge_property<Param_loc>("h:param_loc");
  auto facing = get_vertex_property<FacingType>("v:facing");
  Halfedge_around_vertex_circulator hit, hend;
  hit = hend = halfedges(vsource);
  Eigen::Vector2d new_loc = std::get<3>(etuple).uv.cast<double>();
  int numZCs = 0; // number of zero-crossings in the one-ring
  do {
    Halfedge h0 = *hit;
    Halfedge h1 = next_halfedge(h0);
    Halfedge h2 = next_halfedge(h1);

    Vertex v1 = from_vertex(h1);
    Vertex v2 = from_vertex(h2);
    assert(v1 != vsource && v2 != vsource);
    FacingType f1 = facing[v1];
    FacingType f2 = facing[v2];

    if (f1 == CONTOUR && f2 == CONTOUR)
      return false;
    if (f1 == CONTOUR)
      numZCs++;
    if (f1 != CONTOUR && f2 != CONTOUR && f1 != f2)
      numZCs++;
    if (numZCs > 2) {
      // std::cout<< "CANNOT SHIFT" << std::endl;
      return false;
    }

    if (std::get<3>(etuple).ptexIndex != param_loc[h0].ptexIndex) {
      return false;
    }

    if (!sameSide(param_loc[h1].uv.cast<double>(),
                  param_loc[h2].uv.cast<double>(),
                  param_loc[h0].uv.cast<double>(), new_loc))
      return false;

  } while (++hit != hend);

  return true;
}

Vertex Mesh::splitEdge(const EdgeTuple &etuple, bool shift_allowed) {
  Edge e = std::get<0>(etuple);
  Param_loc e_param = std::get<3>(etuple);

  auto facing = vertex_property<FacingType>("v:facing");
  auto vpositions = vertex_property<Vector3f>("v:point");
  auto vnormals = vertex_property<Vector3f>("v:normal");
  auto param_loc = halfedge_property<Param_loc>("h:param_loc");
  auto colors = get_vertex_property<Vector3f>("v:color");
  auto ndotv = vertex_property<real_t>("v:ndotv");

  // If the zero crossing is too close to an endpoint,
  // we enforce shift
  bool endpoint_shift = false;
  double end_dist =
      fmin((std::get<1>(etuple) - vpositions[vertex(e, 0)]).norm(),
           (std::get<1>(etuple) - vpositions[vertex(e, 1)]).norm());
  double edge_ratio =
      end_dist / (vpositions[vertex(e, 1)] - vpositions[vertex(e, 0)]).norm();
  // if (end_dist < 1e-3 ||
  if (end_dist < 1e-3 || edge_ratio < 0.01) {
    shift_allowed = true;
    endpoint_shift = true;
  }

  Vertex v;

  bool bsplit = true;
  // Try to shift instead
  if (shift_allowed) {
    real_t edge_length =
        (vpositions[vertex(e, 1)] - vpositions[vertex(e, 0)]).norm();
    real_t sub_length = (std::get<1>(etuple) - vpositions[vertex(e, 0)]).norm();
    if (sub_length < MAX_SHIFT_PERCENTAGE * edge_length) {
      // shift first vertex
      v = vertex(e, 0);
    } else if (sub_length > (1.f - MAX_SHIFT_PERCENTAGE) * edge_length) {
      // shift second vertex
      v = vertex(e, 1);
    }

    // An additional case of having t ~0/1
    bool overwrite_vertex = false;
    real_t c = std::get<4>(etuple);
    if (std::fabs(std::min(c, 1 - c)) <
        std::numeric_limits<real_t>::epsilon()) {
      v = (sub_length / edge_length < 0.5) ? vertex(e, 0) : vertex(e, 1);
      overwrite_vertex = true;
    }

    if (v.is_valid() && !is_boundary(v)) {
      bool to_shift_geometry = false;
      if (canShift(v, etuple)) {
        to_shift_geometry = true;
      }

      // If the vertex almost co-locates with an existing vertex, even if we
      // can't shift the geometry, always assign the existing vertex as a
      // contour vertex to avoid degeneration downstream.
      if (overwrite_vertex || to_shift_geometry) {
        facing[v] = CONTOUR;
        ndotv[v] = 0;
        bsplit = false;
      }

      if (to_shift_geometry) {
        vpositions[v] = std::get<1>(etuple);
        vnormals[v] = std::get<2>(etuple);

        auto hit = halfedges(v);
        auto hit_end = hit;
        do {
          param_loc[*hit] = e_param;
        } while (++hit != hit_end);
      }
    }
  }

  Halfedge_around_vertex_circulator chit, chit_end;
  if (bsplit) {
    Param_loc new_loc[2];
    new_loc[0] = e_param;
    // Deal with edges between two ptex faces
    if (!is_boundary(e) && param_loc[halfedge(e, 0)].ptexIndex !=
                               param_loc[halfedge(e, 1)].ptexIndex) {
      real_t c = std::get<4>(etuple);
      new_loc[1].ptexIndex = param_loc[halfedge(e, 1)].ptexIndex;
      new_loc[1].uv = c * param_loc[halfedge(e, 1)].uv +
                      (1.f - c) * param_loc[next_halfedge(halfedge(e, 1))].uv;
    } else {
      new_loc[1] = new_loc[0];
    }
    // save parametric locations of the 4 halfedges (3 for boundaries)
    std::map<Vertex, Param_loc> vertices_param;
    for (int i = 0; i < 2; i++) {
      if (!is_boundary(halfedge(e, i))) {
        Halfedge h = halfedge(e, i);
        vertices_param[from_vertex(h)] = param_loc[h];
        Halfedge nh = next_halfedge(next_halfedge(h));
        vertices_param[from_vertex(nh)] = param_loc[nh];
      }
    }

    v = split(e, std::get<1>(etuple));
    facing[v] = CONTOUR;
    vnormals[v] = std::get<2>(etuple);
    colors[v].setOnes();

    // update parametric locations
    chit = chit_end = halfedges(v);
    do {
      if (!is_boundary(*chit)) {
        if (param_loc[next_halfedge(*chit)].ptexIndex == new_loc[0].ptexIndex) {
          param_loc[*chit] = new_loc[0];
        } else {
          param_loc[*chit] = new_loc[1];
        }
      }
      Halfedge oh = opposite_halfedge(*chit);
      if (!is_boundary(oh)) {
        param_loc[oh] = vertices_param[from_vertex(oh)];
      }
      ++chit;
    } while (chit != chit_end);
  }

  return v;
}

void Mesh::insertContours(const Vector3f &view_pos, const Subdiv *subdiv,
                          bool shift_allowed) {
  auto facing = get_vertex_property<FacingType>("v:facing");
  if (!facing)
    computeConsistencySubdivided(view_pos);

  // Detect FB edges
  tbb::concurrent_vector<EdgeTuple> edges_to_split;

  edges_to_split.reserve(std::sqrt(n_faces()));
  for (uint32_t i = 0; i != n_edges(); ++i) {
    Edge e = Edge(i);
    // 1. Check if the zero-crossing exists on the current edge
    // CONTOUR: the face is roughly parallel to the camera ray.
    // (which shouldn't happen under the generic camera assumption.)
    if (facing[vertex(e, 0)] == CONTOUR || facing[vertex(e, 1)] == CONTOUR)
      continue;

    // 2. Bisect search for the contour zero point on the limit surface
    // along the edge. (different from the fixed grid search in the
    // splitEdgesWithZeroCrossing)
    if (facing[vertex(e, 0)] != facing[vertex(e, 1)]) {
      Vector3f p, n;
      Param_loc param;
      real_t c = bisect_search(view_pos, subdiv, e, p, n, param);
      edges_to_split.push_back(EdgeTuple(e, p, n, param, c));
    }
  }

  // Actually split FB edges
  for (auto etuple : edges_to_split) {
    FacingType f1 = facing[vertex(std::get<0>(etuple), 0)];
    FacingType f2 = facing[vertex(std::get<0>(etuple), 1)];

    if (f1 == CONTOUR || f2 == CONTOUR || f1 == f2) {
      // skip if not an FB edge anymore
      continue;
    }
    splitEdge(etuple, shift_allowed);
  }

  garbage_collection();

#ifndef NDEBUG
  auto param_loc = get_halfedge_property<Param_loc>("h:param_loc");
  // Assert validity of parametric locations
  Halfedge_iterator hit;
  for (hit = halfedges_begin(); hit != halfedges_end(); ++hit)
    if (!is_boundary(*hit))
      assert(param_loc[*hit].ptexIndex != -1);
#endif

  updateVBO(this, view_pos);
}

void Mesh::updateVBO(Mesh *mesh, const Vector3f &view_pos) {
  // Update VBO
  auto VBO = mesh->face_property<FacingType>("f:VBO");
  auto facing = mesh->vertex_property<FacingType>("v:facing");
  parallel_for(
      tbb::blocked_range<uint32_t>(0u, (uint32_t)mesh->n_faces(), GRAIN_SIZE),
      [&](const tbb::blocked_range<uint32_t> &range) {
        Vertex v[3];
        for (uint32_t i = range.begin(); i != range.end(); ++i) {
          Face f = Face(i);

          mesh->compute_face_normal(f);
          mesh->verticesOfFace(f, v);

          int nbFacing[4] = {0, 0, 0, 0};

          for (int i = 0; i < 3; i++)
            nbFacing[facing[v[i]]]++;

          if ((nbFacing[FRONT] == 3) ||
              (nbFacing[FRONT] == 2 && nbFacing[CONTOUR] == 1) ||
              (nbFacing[FRONT] == 1 && nbFacing[CONTOUR] == 2))
            VBO[f] = FRONT;
          else if ((nbFacing[BACK] == 3) ||
                   (nbFacing[BACK] == 2 && nbFacing[CONTOUR] == 1) ||
                   (nbFacing[BACK] == 1 && nbFacing[CONTOUR] == 2))
            VBO[f] = BACK;
          else {
            VBO[f] = UNDEFINED;
            logger().error("UNDEFINED VBO AFTER CONTOUR INSERTION: " +
                           std::to_string(nbFacing[FRONT]) +
                           std::to_string(nbFacing[BACK]) +
                           std::to_string(nbFacing[CONTOUR]));
            VBO[f] = FRONT;
          }
        }
      });

  mesh->updateIndicesFromVBO(view_pos);
}

void Mesh::computeContourChains(ContourMode mode) {
  auto is_contour = get_edge_property<real_t>("e:contour");
  auto visited = get_edge_property<real_t>("e:visited");
  if (!visited) {
    visited = add_edge_property<real_t>("e:visited", false);
  } else {
    for (size_t i = 0; i < visited.vector().size(); ++i)
      visited.vector()[i] = false;
  }
  m_chains.clear();
  m_chain_lengths.clear();

  size_t sum = 0;

  for (size_t i = 0; i < m_contour_edges.size(); ++i) {
    Edge e = m_contour_edges[i];
    if (visited[e])
      continue;
    std::list<Edge> chain;
    chain.push_back(e);
    visited[e] = true;
    bool forward = true;
    { // MESH & VBO CONTOURS
      Halfedge nh = next_halfedge(halfedge(e, 0));
      Edge ne = edge(nh);
      Edge ie = ne;
      while (ne != e) {
        if (is_contour[ne] >= 0 && !visited[ne]) {
          if (forward)
            chain.push_back(ne);
          else
            chain.push_front(ne);
          nh = opposite_halfedge(nh);
          ie = edge(nh);
        } else {
          nh = next_halfedge(opposite_halfedge(nh));
          if (edge(nh) == ie) { // full loop around the one-ring
            if (forward) {
              // Chain backward
              forward = false;
              nh = next_halfedge(halfedge(e, 1));
              ie = edge(nh);
            } else {
              break;
            }
          }
        }
        visited[ne] = true;
        ne = edge(nh);
      }
    }
    m_chains.push_back(chain);
    m_chain_lengths.push_back(chain.size());
    sum += chain.size();
  }

  auto vpositions = get_vertex_property<Vector3f>("v:point");
  auto fedges_e = edge_property<FEdge *>("e:fedge", nullptr);
  auto svertices_e = edge_property<std::vector<SVertex *>>("e:svertex");

  int i = 0;
  { // MESH CONTOURS
    auto svertices_v = vertex_property<SVertex *>("v:svertex", nullptr);
    auto concave = get_edge_property<bool>("e:concave");
    // create a fedge between two vertices of the mesh
    auto add_fedge = [&](Edge edge, int i0, SVertex *sv0) {
      if (!sv0) {
        Vertex v0 = vertex(edge, i0);
        if (svertices_v[v0]) {
          sv0 = svertices_v[v0];
          // 3D intersection
          sv0->m_nature = VertexNature::INTERSECTION_3D | sv0->m_nature;
        } else {
          sv0 = new SVertex(vpositions[v0]);
          svertices_v[v0] = sv0;
          m_svertices.push_back(sv0);
        }
      }
      Vertex v1 = vertex(edge, (i0 + 1) % 2);
      SVertex *sv1;
      if (svertices_v[v1]) {
        sv1 = svertices_v[v1];
        if (sv1->m_fedges.size() >= 2) { // 3D intersection
          sv1->m_nature = VertexNature::INTERSECTION_3D | sv1->m_nature;
        }
      } else {
        sv1 = new SVertex(vpositions[v1]);
        svertices_v[v1] = sv1;
        m_svertices.push_back(sv1);
      }
      FEdge *fedge = new FEdgeSharp(sv0, sv1, edge, EdgeNature::SHARP_CONTOUR);
      m_fedges.push_back(fedge);
      fedges_e[edge] = fedge;
      sv0->m_fedges.push_back(fedge);
      sv1->m_fedges.push_back(fedge);

      if (!svertices_e[edge].empty()) { // 3D intersection
        for (auto v : svertices_e[edge]) {
          v->m_nature = VertexNature::INTERSECTION_3D | v->m_nature;
        }
      }
      return sv1;
    };

    m_chain_indices.resize(sum + m_chains.size() * 3);
    for (auto chain : m_chains) {
      Edge e = chain.front();
      SVertex *sv = nullptr;
      if (chain.size() == 1) {
        m_chain_indices[i++] = vertex(e, 0).idx();
        m_chain_indices[i++] = vertex(e, 0).idx();
        m_chain_indices[i++] = vertex(e, 1).idx();
        m_chain_indices[i++] = vertex(e, 1).idx();
        add_fedge(e, 0, sv);
      } else {
        auto it = chain.begin();
        it++;
        Edge ne = (*it);
        Vertex prev;
        bool prevConcave;
        if (vertex(e, 0) == vertex(ne, 0) || vertex(e, 0) == vertex(ne, 1)) {
          m_chain_indices[i++] = vertex(e, 1).idx();
          m_chain_indices[i++] = vertex(e, 1).idx();
          prev = vertex(e, 0);
          sv = add_fedge(e, 1, sv);
        } else if (vertex(e, 1) == vertex(ne, 0) ||
                   vertex(e, 1) == vertex(ne, 1)) {
          m_chain_indices[i++] = vertex(e, 0).idx();
          m_chain_indices[i++] = vertex(e, 0).idx();
          prev = vertex(e, 1);
          sv = add_fedge(e, 0, sv);
        } else {
          std::runtime_error("Invalid contour chain");
        }
        prevConcave = concave[e];
        m_chain_indices[i++] = prev.idx();
        while (it != chain.end()) {
          e = *it;
          if (prevConcave != concave[e]) { // curtain fold
            sv->m_nature = VertexNature::CONTOUR_CURTAIN_FOLD | sv->m_nature;
            m_debug_points.push_back(sv->m_pos3D);
            m_point_colors.push_back(curtainFoldColor.head<3>().cast<real_t>());
          }
          prevConcave = concave[e];
          if (prev != vertex(e, 0)) {
            prev = vertex(e, 0);
            sv = add_fedge(e, 1, sv);
          } else {
            prev = vertex(e, 1);
            sv = add_fedge(e, 0, sv);
          }
          m_chain_indices[i++] = prev.idx();
          it++;
        }
        m_chain_indices[i++] = prev.idx();
      }
    }
  }
}

void Mesh::computeSweepLineIntersections(ContourMode mode, const Camera &camera,
                                         const Vector2i &viewport) {
  typedef Segment<FEdge *> segment;
  typedef Intersection<segment> intersection;

  // we only want image-space intersections between any pair of edges except
  // for smooth-sharp pairs that intersect on the surface
  struct silhouette_binary_rule_no_same_face
      : public binary_rule<segment, segment> {
    silhouette_binary_rule_no_same_face() : binary_rule<segment, segment>() {}
    virtual bool operator()(segment &s1, segment &s2) {
      FEdge *e1 = s1.edge();
      FEdge *e2 = s2.edge();

      if ((e1->isSmooth() && e2->isSmooth()) ||
          (!e1->isSmooth() && !e2->isSmooth()))
        return true;

      if (e1->isSmooth()) {
        if (e1->m_vertexA->m_nature & VertexNature::INTERSECTION_3D ||
            e1->m_vertexB->m_nature & VertexNature::INTERSECTION_3D)
          return false;
      } else {
        if (e2->m_vertexA->m_nature & VertexNature::INTERSECTION_3D ||
            e2->m_vertexB->m_nature & VertexNature::INTERSECTION_3D)
          return false;
      }
      return true;
    }
  };

  vector<segment *> segments;

  sort(m_svertices.begin(), m_svertices.end(), VertexCompare(EPSILON));

  for (auto fedge : m_fedges) {
    segments.push_back(new segment(fedge, fedge->m_vertexA, fedge->m_vertexB));
    fedge->m_segment = segments.back();
  }

  SweepLine<FEdge *> SL;
  vector<segment *> vsegments;
  silhouette_binary_rule_no_same_face sbr;
  // Iterate over every sorted curve vertices
  for (auto v : m_svertices) {
    // Add new fedges
    for (auto fedge : v->m_fedges) {
      assert(fedge);
      vsegments.push_back(fedge->m_segment);
    }
    SL.process(v->m_pos2D.cast<double>().head<2>(), vsegments, sbr);
    vsegments.clear();
  }

  // retrieve the intersections
  vector<intersection *> &intersections = SL.intersections();
  m_2D_intersections.second = intersections.size();

  auto addIntersectionSVertex = [&](FEdge *e, real_t t) {
    SVertex *svA = e->m_vertexA;
    SVertex *svB = e->m_vertexB;
    Vector3f pos2D = svA->m_pos2D + t * (svB->m_pos2D - svA->m_pos2D);
    Vector3f pos3D = unproject(pos2D, camera.viewMatrix().matrix(),
                               camera.projectionMatrix(), viewport);
    // create new SVertex
    SVertex *sv_intersec = new SVertex(pos3D);
    sv_intersec->m_pos2D = pos2D;
    sv_intersec->m_nature =
        sv_intersec->m_nature | VertexNature::INTERSECTION_2D;
    m_svertices.push_back(sv_intersec);
    return sv_intersec;
  };

  // create new svertices
  for (auto intersect : intersections) {
    intersect->m_vertexA =
        addIntersectionSVertex(intersect->m_edgeA->edge(), intersect->m_tA);
    intersect->m_vertexB =
        addIntersectionSVertex(intersect->m_edgeB->edge(), intersect->m_tB);

    // create an FEdge between them
    FEdge *new_edge = new FEdge(intersect->m_vertexA, intersect->m_vertexB,
                                EdgeNature::IMAGE_INTERSECTION);
    intersect->m_vertexA->m_fedges.push_back(new_edge);
    intersect->m_vertexB->m_fedges.push_back(new_edge);
    m_fedges.push_back(new_edge);
  }

  auto addFEdge = [&](FEdge *e, SVertex *prev_sv, SVertex *curr_sv) {
    FEdge *new_edge;
    if (e->m_nature & EdgeNature::SHARP_CONTOUR ||
        e->m_nature & EdgeNature::BOUNDARY) {
      FEdgeSharp *fe_s = dynamic_cast<FEdgeSharp *>(e);
      new_edge = new FEdgeSharp(prev_sv, curr_sv, fe_s->m_edge, e->m_nature);
    } else if (e->m_nature & EdgeNature::SMOOTH_CONTOUR) {
      FEdgeSmooth *fe_s = dynamic_cast<FEdgeSmooth *>(e);
      new_edge = new FEdgeSmooth(prev_sv, curr_sv, fe_s->m_edgeA, fe_s->m_tA,
                                 fe_s->m_edgeB, fe_s->m_tB, e->m_nature);
    } else if (e->m_nature & EdgeNature::SURFACE_INTERSECTION) {
      FEdgeIntersection *fe_i = dynamic_cast<FEdgeIntersection *>(e);
      new_edge =
          new FEdgeIntersection(prev_sv, curr_sv, fe_i->m_f1, fe_i->m_f2);
    }
    curr_sv->m_fedges.push_back(new_edge);
    m_fedges.push_back(new_edge);
  };

  // retrieve the intersected edges
  std::vector<segment *> &iedges = SL.intersectedEdges();
  // split egdes at intersection
  for (auto s : iedges) {
    vector<intersection *> &eIntersections = s->intersections();
    // first need to sort these intersections along the edge
    std::sort(eIntersections.begin(), eIntersections.end(),
              less_Intersection<FEdge>(s));

    FEdge *e = s->edge();
    // retrieve the new svertices for all intersections on this edge
    vector<SVertex *> edgeSVertices;
    for (auto intersect : eIntersections) {
      if (e == intersect->m_edgeA->edge()) {
        edgeSVertices.push_back(intersect->m_vertexA);
      } else {
        assert(e == intersect->m_edgeB->edge());
        edgeSVertices.push_back(intersect->m_vertexB);
      }
    }

    // remove old FEdge from its extremities
    SVertex *svA = e->m_vertexA;
    SVertex *svB = e->m_vertexB;
    svA->m_fedges.remove(e);
    svB->m_fedges.remove(e);
    // chain new vertices
    SVertex *prev_sv = svA;
    for (auto curr_sv : edgeSVertices) {
      addFEdge(e, prev_sv, curr_sv);
      prev_sv = curr_sv;
    }
    addFEdge(e, prev_sv, svB);

    // delete the old FEdge 
    auto it = std::find(m_fedges.begin(), m_fedges.end(), e);
    if (it != m_fedges.end())
      m_fedges.erase(it);
    delete e;
  }

  // clean-up
  for (auto fedge : m_fedges)
    fedge->m_segment = nullptr;
  for (auto s : segments)
    delete s;
}

void Mesh::computeContourVisibility(const Vector3f &view_pos, ContourMode mode,
                                    const ProgressCallback &progress) {

  if (!m_bvh) {
    m_bvh = new BVH(
        &m_indices, &get_vertex_property<Vector3f>("v:point").vector(),
        &get_vertex_property<Vector3f>("v:normal").vector(), boundingBox());
    m_bvh->build(progress);
    if (progress)
      m_bvh->printStatistics();
  }

  auto concave = get_edge_property<bool>("e:concave");
  ray_intersections.clear();

  if (mode != INTERPOLATED_CONTOUR) {
    tbb::concurrent_vector<Vector2u> invisible_indices;
    invisible_indices.reserve(0.5 * m_mesh_contours_indices.size());
    tbb::concurrent_vector<Vector2u> visible_indices;
    visible_indices.reserve(0.5 * m_mesh_contours_indices.size());

    auto visibility = get_edge_property<bool>("e:visible");
    if (!visibility)
      visibility = add_edge_property<bool>("e:visible");

    auto vpositions = get_vertex_property<Vector3f>("v:point");
    parallel_for(tbb::blocked_range<uint32_t>(
                     0u, (uint32_t)m_contour_edges.size(), GRAIN_SIZE),
                 [&](const tbb::blocked_range<uint32_t> &range) {
                   for (uint32_t i = range.begin(); i != range.end(); ++i) {
                     Edge e = m_contour_edges[i];
                     const Vector3f &a = vpositions[vertex(e, 0)];
                     const Vector3f &b = vpositions[vertex(e, 1)];

                     if (concave[e]) {
                       //  ray_intersections.push_back(0.5f * (a + b));
                       invisible_indices.push_back(
                           Vector2u(vertex(e, 0).idx(), vertex(e, 1).idx()));
                       continue;
                     }

                     Vector3f dir = 0.5f * (a + b) - view_pos;
                     uint32_t idx;
                     real_t t;
                     visibility[e] = true;
                     Ray ray = Ray(view_pos, dir.normalized(), 0, dir.norm());
                     bool hit = m_bvh->rayIntersect(ray, idx, t);
                     if (hit && int(idx) != face(e, 0).idx() &&
                         int(idx) != face(e, 1).idx()) {
                       ray_intersections.push_back(ray(t));
                       visibility[e] = false;
                       invisible_indices.push_back(
                           Vector2u(vertex(e, 0).idx(), vertex(e, 1).idx()));
                     } else {
                       visible_indices.push_back(
                           Vector2u(vertex(e, 0).idx(), vertex(e, 1).idx()));
                     }
                   }
                 });

    m_mesh_contours_indices =
        MatrixXu(2, visible_indices.size() + invisible_indices.size());
    for (size_t i = 0; i < visible_indices.size(); ++i)
      m_mesh_contours_indices.col(i) = visible_indices.at(i);
    for (size_t i = 0; i < invisible_indices.size(); ++i)
      m_mesh_contours_indices.col(i + visible_indices.size()) =
          invisible_indices.at(i);
    m_visible = 2 * visible_indices.size();
  } else { // INTERPOLATED CONTOURS
    tbb::concurrent_vector<std::pair<Vector3f, Vector3f>> invisible_edges;
    invisible_edges.reserve(0.25 * m_interpolated_contours.cols());
    tbb::concurrent_vector<std::pair<Vector3f, Vector3f>> visible_edges;
    visible_edges.reserve(0.25 * m_interpolated_contours.cols());
    parallel_for(
        tbb::blocked_range<uint32_t>(
            0u, (uint32_t)m_interpolated_contours.cols() / 2, GRAIN_SIZE),
        [&](const tbb::blocked_range<uint32_t> &range) {
          for (uint32_t i = range.begin(); i != range.end(); ++i) {
            const Vector3f &p1 = m_interpolated_contours.col(2 * i);
            const Vector3f &p2 = m_interpolated_contours.col(2 * i + 1);
            Vector3f dir = 0.5f * (p1 + p2) - view_pos;
            uint32_t idx;
            real_t t;
            Ray ray = Ray(view_pos, dir.normalized(), 0, dir.norm());
            bool hit = m_bvh->rayIntersect(ray, idx, t);
            if (hit && int(idx) != m_interpolated_contour_faces[i]) {
              ray_intersections.push_back(ray(t));
              invisible_edges.push_back(std::make_pair(p1, p2));
            } else {
              visible_edges.push_back(std::make_pair(p1, p2));
            }
          }
        });
    for (size_t i = 0; i < visible_edges.size(); ++i) {
      m_interpolated_contours.col(i * 2) = visible_edges.at(i).first;
      m_interpolated_contours.col(i * 2 + 1) = visible_edges.at(i).second;
    }
    for (size_t i = 0; i < invisible_edges.size(); ++i) {
      m_interpolated_contours.col((i + visible_edges.size()) * 2) =
          invisible_edges.at(i).first;
      m_interpolated_contours.col((i + visible_edges.size()) * 2 + 1) =
          invisible_edges.at(i).second;
    }
    m_visible = 2 * visible_edges.size();
  }
}

void Mesh::floodFill() {
  // Flood fill ids to extract patches
  auto patchID = get_face_property<int>("f:patchID");
  if (!patchID) {
    patchID = add_face_property<int>("f:patchID", -1);
  } else {
    patchID.vector().assign(patchID.vector().size(), -1);
  }
  auto is_cut = edge_property<bool>("e:cut", false);
  auto is_contour = get_edge_property<real_t>("e:contour");

  int id = 0;
  m_patch_lengths.clear();
  m_patch_indices.clear();
  m_patch_indices.reserve(n_faces());
  for (Face f : faces()) {
    if (patchID[f] > -1) // already visited
      continue;
    int nb_faces = 0;
    std::queue<Face> queue;
    queue.push(f);
    while (!queue.empty()) {
      Face current_f = queue.front();
      queue.pop();
      if (patchID[current_f] > -1) // already processed, skip
        continue;
      nb_faces++;
      Vertex v[3];
      verticesOfFace(current_f, v);
      for (int i = 0; i < 3; i++)
        m_patch_indices.push_back(v[i].idx());
      patchID[current_f] = id;
      // add adjacent faces if no cut boundary or contour is crossed
      Halfedge_around_face_circulator hit, hend;
      hit = hend = halfedges(current_f);
      do {
        Edge e = edge(*hit);
        if (is_cut[e] || is_contour[e] >= 0.f || is_boundary(e))
          continue;
        // adjacent_f is in the same patch
        Face adjacent_f = face(opposite_halfedge(*hit));
        if (patchID[adjacent_f] == -1)
          queue.push(adjacent_f);
      } while (++hit != hend);
    }
    id++;
    m_patch_lengths.push_back(nb_faces);
  }
  assert(m_patch_indices.size() / 3 == n_faces());
}

void Mesh::markPatchBoundaryEdges(bool makeDisk) {
  auto is_boundary_edge = get_edge_property<bool>("e:is_boundary");
  if (!is_boundary_edge) {
    is_boundary_edge = add_edge_property<bool>("e:is_boundary", false);
  }
  is_boundary_edge = edge_property<bool>("e:is_boundary");

  auto patchID = get_face_property<int>("f:patchID");
  assert(patchID);

  // mark patch boundary edges, i.e., half-edges shared by two faces with
  // diffrent ids
  auto patchBoundary = get_halfedge_property<int>("h:patchBoundary");
  if (!patchBoundary) {
    add_halfedge_property<int>("h:patchBoundary", -1);
  }
  patchBoundary = halfedge_property<int>("h:patchBoundary");

  parallel_for(
      tbb::blocked_range<uint32_t>(0u, (uint32_t)n_edges(), GRAIN_SIZE),
      [&](const tbb::blocked_range<uint32_t> &range) {
        for (uint32_t i = range.begin(); i != range.end(); ++i) {
          Edge e(i);

          if (is_boundary(e))
            is_boundary_edge[e] = true;
          else
            is_boundary_edge[e] = false;

          Halfedge h0 = halfedge(e, 0);
          Halfedge h1 = halfedge(e, 1);
          // test if this edge is on an open boundary
          if (is_boundary(h0)) {
            patchBoundary[h1] = patchID[face(h1)];
          } else if (is_boundary(h1)) {
            patchBoundary[h0] = patchID[face(h0)];
          } else {
            // test if this edge is on a patch boundary
            int id0 = patchID[face(h0)];
            int id1 = patchID[face(h1)];
            if (id0 != id1) {
              patchBoundary[h0] = id0;
              patchBoundary[h1] = id1;
            }
          }
        }
      });

  if (!makeDisk)
    return;

  // make each patch homeomorphic to a disk by connecting boundary loops
  int start_index = 0;
  for (size_t id = 0; id < m_patch_lengths.size(); id++) {
    std::vector<std::vector<int32_t>> cuts;
    MatrixXi F;
    F.resize(m_patch_lengths[id], 3);
    for (int f = 0; f < m_patch_lengths[id]; f++) {
      F.row(f) = Vector3i(m_patch_indices[start_index + f * 3],
                          m_patch_indices[start_index + f * 3 + 1],
                          m_patch_indices[start_index + f * 3 + 2]);
    }
    igl::cut_to_disk(F, cuts);

    // tag edges accordingly
    for (auto &seam : cuts) {
      Vertex v1(seam.front()), v2;
      for (size_t i = 1; i < seam.size(); i++) {
        v2 = Vertex(seam[i]);
        Halfedge_around_vertex_circulator hit, hit_end;
        hit = hit_end = halfedges(v1);
        do {
          if (to_vertex(*hit) == v2) {
            break;
          }
        } while (++hit != hit_end);
        // skip previously tagged patch boundaries
        if (patchBoundary[*hit] == -1) {
          patchBoundary[*hit] = id;
          patchBoundary[opposite_halfedge(*hit)] = id;
        }
        v1 = v2;
      }
    }
    start_index += m_patch_lengths[id] * 3;
  }
}

void Mesh::extractPatchBoundaries(bool makeDisk) {
  markPatchBoundaryEdges(makeDisk);
  auto patchBoundary = get_halfedge_property<int>("h:patchBoundary");

  auto vpositions = get_vertex_property<Vector3f>("v:point");
  auto vnormals = get_vertex_property<Vector3f>("v:normal");
  auto ndotv = get_vertex_property<real_t>("v:ndotv");
  auto facing = get_vertex_property<FacingType>("v:facing");
  auto visited = add_halfedge_property<bool>("h:visited", false);
  auto vertex_mapping = add_vertex_property<Vertex>("v:mapping");
  auto VBO = get_face_property<FacingType>("f:VBO");
  m_patch_boundaries.resize(m_patch_lengths.size());
  m_patch_facing.resize(m_patch_lengths.size());

  // Create a new mesh to store the resuling tessellation
  res_mesh = new Mesh();
  auto new_normals = res_mesh->add_vertex_property<Vector3f>("v:normal");
  auto new_ndotv = res_mesh->add_vertex_property<real_t>("v:ndotv");
  auto new_facing = res_mesh->add_vertex_property<FacingType>("v:facing");
  res_mesh->add_face_property<int>("f:patchID");

  for (Halfedge fh : halfedges()) {
    if (patchBoundary[fh] < 0 || visited[fh])
      continue;
    Halfedge h = fh;
    int current_id = patchBoundary[h];
    m_patch_boundaries[current_id].resize(
        m_patch_boundaries[current_id].size() + 1);
    std::vector<Vertex> &polyline = m_patch_boundaries[current_id].back();
    if (face(h).idx() < 0)
      continue;
    m_patch_facing[current_id] = VBO[face(h)];
    do {
      if (visited[h])
        break;
      visited[h] = true;
      Vertex v = from_vertex(h);
      Vertex new_vertex = vertex_mapping[v];
      if (!new_vertex.is_valid()) {
        new_vertex = res_mesh->add_vertex(vpositions[v]);
        vertex_mapping[v] = new_vertex;
      }
      m_seams_indices.push_back(v.idx());
      new_normals[new_vertex] = vnormals[v];
      new_facing[new_vertex] = facing[v];
      new_ndotv[new_vertex] = ndotv[v];
      polyline.push_back(new_vertex);
      // search the next halfedge with the same patchID
      h = next_halfedge(h);
      while (patchBoundary[h] != current_id) {
        h = cw_rotated_halfedge(h);
      }
      assert(patchBoundary[h] == current_id);
      if (visited[h] == true && h != fh) {
        ray_intersections.push_back(vpositions[from_vertex(h)]);
        cout << "!!! Unexpected patch boundary edge !!!" << endl;
        break;
      }
    } while (h != fh);
    m_seams_lengths.push_back(polyline.size());
  }
  remove_halfedge_property(visited);
  remove_vertex_property(vertex_mapping);
}

Halfedge Mesh::get_patch_halfedge(Vertex const &v, size_t patch_id) const {
  auto patchBoundary = get_halfedge_property<int>("h:patchBoundary");
  Halfedge h = halfedge(v);
  auto hit = halfedges(v);
  auto hit_end = hit;

  do {
    if (patchBoundary[*hit] != (int)patch_id)
      continue;
    h = *hit;
    break;
  } while (++hit != hit_end);

  contess_assert_msg(patchBoundary[h] == (int)patch_id,
                     "Can't find a halfedge as the boundary of patch " +
                         std::to_string(patch_id));
  return h;
}

void Mesh::build_bvh() {
  if (m_bvh) {
    delete m_bvh;
    m_indices.reserve(3 * n_faces());
    // face iterator
    Face_iterator fit, fend = faces_end();
    // vertex circulator
    Vertex vertices[3];
    for (fit = faces_begin(); fit != fend; ++fit) {
      verticesOfFace(*fit, vertices);
      for (Vertex v : vertices)
        m_indices.push_back(v.idx());
    }
    m_sorted_indices = m_indices;
  }
  m_bvh = new BVH(
      &m_indices, &get_vertex_property<Vector3f>("v:point").vector(),
      &get_vertex_property<Vector3f>("v:normal").vector(), boundingBox());
  m_bvh->build(nullptr);
}

Mesh *Mesh::triangulateSimplePatches(std::vector<int> patchIDs,
                                     const Subdiv *subdiv,
                                     const Matrix4f &model,
                                     const Matrix4f &proj) {
  std::set<int> patchIDs_set(patchIDs.begin(), patchIDs.end());

  // Remesh each patch
  Matrix4f invC = model.matrix().inverse();
  Matrix3f proj3;
  proj3 << proj.topLeftCorner<2, 3>(), proj.bottomLeftCorner<1, 3>();
  Vector3f view_pos = invC.col(3).head<3>();

  build_bvh();

  auto new_positions = res_mesh->get_vertex_property<Vector3f>("v:point");
  auto new_patchID = res_mesh->get_face_property<int>("f:patchID");

  for (size_t p = 0; p < m_patch_boundaries.size(); ++p) {
    if (patchIDs.empty() ||
        (patchIDs[0] != -1 && patchIDs_set.count(int(p)) == 0))
      continue;

    // only one disk-like surface
    assert(m_patch_boundaries[p].size() == 1);

    auto &boundary_polygon = m_patch_boundaries[p].front();

    Eigen::MatrixXd P, V_out;
    Eigen::VectorXi R;
    Eigen::MatrixXi F_out;

    std::vector<Vector2f> currentPolygon;
    currentPolygon.resize(boundary_polygon.size());
    for (std::size_t j = 0; j < boundary_polygon.size(); j++) {
      Vector4f tmp;
      tmp << new_positions[boundary_polygon[j]], 1.f;
      tmp = model * tmp;
      tmp = proj * tmp;
      tmp = tmp.array() / tmp(3);
      currentPolygon[j] = tmp.head<2>();
    }

    {
      SVG svgWriter("test.svg", Vector2i(800, 600));
      for (size_t i = 0; i < m_patch_boundaries[p].size(); i++) {
        std::vector<Vector2f> pl;
        auto &polyline = m_patch_boundaries[p][i];
        for (size_t j = 0; j < polyline.size(); j++) {
          Vector4f tmp;
          tmp << new_positions[polyline[j]], 1.f;
          tmp = model * tmp;
          tmp = proj * tmp;
          tmp = tmp.array() / tmp(3);
          tmp = tmp.array() * 0.5f + 0.5f;
          tmp(0) = tmp(0) * 800;
          tmp(1) = tmp(1) * 600;

          pl.push_back(tmp.head<2>());
        }
        pl.push_back(pl.front());
        svgWriter.writePolyline(pl, 1, alphabetColors[i % 26], false);
      }
    }

    std::vector<Vertex> indexToVertex;
    P.resize(currentPolygon.size(), 2);
    if (m_patch_facing[p] == BACK) {
      std::reverse(std::begin(indexToVertex), std::end(indexToVertex));
      for (size_t i = 0; i < currentPolygon.size(); i++) {
        P.row(i) = currentPolygon[currentPolygon.size() - i - 1].cast<double>();
      }
    } else {
      for (size_t i = 0; i < currentPolygon.size(); i++) {
        P.row(i) = currentPolygon[i].cast<double>();
      }
    }
    R.setZero(P.rows());

    // Tessellate polygon
    if (!Shor_van_wyck(P, R, "", V_out, F_out, false)) {
      return this;
    }

    indexToVertex.resize(V_out.rows());

    // 2D tessellation
    parallel_for(
        tbb::blocked_range<uint32_t>(0, (uint32_t)V_out.rows(), GRAIN_SIZE),
        [&](const tbb::blocked_range<uint32_t> &range) {
          for (uint32_t i = range.begin(); i != range.end(); ++i) {
            Vector3f q;
            q << V_out.row(i).cast<real_t>().transpose(), 0.f;
            Vertex new_vertex = res_mesh->add_vertex(q);
            assert(new_vertex.is_valid());
            indexToVertex[i] = new_vertex;
          }
        });

    int mesh_faces_counter = 0;
    for (int i = 0; i < F_out.rows(); i++) {
      Vertex tris[3] = {indexToVertex[F_out.row(i)(0)],
                        indexToVertex[F_out.row(i)(1)],
                        indexToVertex[F_out.row(i)(2)]};
      if (m_patch_facing[p] == BACK)
        std::swap(tris[1], tris[2]);

      Face new_face = res_mesh->add_triangle(tris[0], tris[1], tris[2]);
      new_patchID[new_face] = p;
      res_mesh->m_patch_indices.push_back(tris[0].idx());
      res_mesh->m_patch_indices.push_back(tris[1].idx());
      res_mesh->m_patch_indices.push_back(tris[2].idx());
      mesh_faces_counter++;
    }
    res_mesh->m_patch_lengths.push_back(mesh_faces_counter);
  }

  res_mesh->init();
  res_mesh->add_face_property<FacingType>("f:VBO");
  updateVBO(res_mesh, view_pos);
  return res_mesh;
}

void Mesh::triangulateSimplePatches(ContourMode mode, const Subdiv *subdiv,
                                    const Matrix4f &model, const Matrix4f &proj,
                                    const Vector2i &viewportSize, MatrixXf &V2,
                                    MatrixXi &F2) {
  // Input polygon
  Eigen::MatrixXd V;
  MatrixXi E;
  Eigen::MatrixXd H;

  // boundary markers
  VectorXi VM, EM, VM2, EM2;

  Matrix3Xf contours = get_contours(mode, false);
  int offset = (mode == INTERPOLATED_CONTOUR ? 0 : 1);
  int start_index = 0;
  std::vector<Vector2f> polyline;
  std::vector<int> indices;
  for (size_t i = 0; i < m_chain_lengths.size(); i++) {
    for (int j = 1; j < m_chain_lengths[i] + offset; ++j) {
      Vector4f tmp;
      tmp << contours.col(j + start_index), 1;
      tmp = model * tmp;
      tmp = proj * tmp;
      tmp = tmp.array() / tmp(3);
      polyline.push_back(tmp.head<2>());
      indices.push_back(j + start_index);
    }
    start_index += m_chain_lengths[i] + offset + 2;
  }

  V.resize(polyline.size(), 2);
  E.resize(polyline.size(), 2);
  VM.resize(polyline.size());
  EM.setZero(E.rows());
  for (size_t i = 0; i < polyline.size(); i++) {
    V.row(i) = polyline[i].cast<double>();
    VM(i) = indices[i];
    if (int(i) == (m_chain_lengths[0] - 1))
      E.row(i) = Vector2i(i, 0);
    else
      E.row(i) = Vector2i(i, i + 1);
  }
  E.row(polyline.size() - 1) =
      Vector2i(polyline.size() - 1, m_chain_lengths[0]);

  H.resize(1, 2);
  H << 0, 0;

  Eigen::MatrixXd V_out;
  // Triangulate the interior
  igl::triangle::triangulate(V, E, H, VM, EM, "a0.001qQYY", V_out, F2, VM2,
                             EM2);

  // Find the 3D positions of the generated vertices by ray-casting
  V2.resize(3, V_out.rows());
  std::vector<FacingType> facing;
  facing.resize(V_out.rows());

  Matrix4f invC = model.matrix().inverse();
  Matrix3f proj3;
  proj3 << proj.topLeftCorner<2, 3>(), proj.bottomLeftCorner<1, 3>();
  Vector3f view_pos = invC.col(3).head<3>();

  if (!m_bvh) {
    m_bvh = new BVH(
        &m_indices, &get_vertex_property<Vector3f>("v:point").vector(),
        &get_vertex_property<Vector3f>("v:normal").vector(), boundingBox());
    m_bvh->build(nullptr);
  }

  auto param_loc = get_halfedge_property<Param_loc>("h:param_loc");
  parallel_for(
      tbb::blocked_range<uint32_t>(0u, (uint32_t)V_out.rows(), GRAIN_SIZE),
      [&](const tbb::blocked_range<uint32_t> &range) {
        for (uint32_t i = range.begin(); i != range.end(); ++i) {
          if (VM2(i) != 0) { // boundary vertex
            V2.col(i) = contours.col(VM2(i)).cast<real_t>();
            continue;
          }
          Vector3f q;
          q << V_out.row(i).cast<real_t>().transpose(), 1.f;
          Vector3f dir = invC.topLeftCorner<3, 3>() * (proj3.inverse() * q);
          Ray ray = Ray(view_pos, dir.normalized());
          uint32_t idx;
          real_t t;
          Vector2f uv;
          bool hit = m_bvh->rayIntersect(ray, idx, t, &uv);
          if (hit) {
            Vertex v[3];
            verticesOfFace(Face(idx), v);
            Halfedge h[3];
            halfedgesOfFace(Face(idx), h);
            Halfedge sh[3];
            // sort halfedge according to cw vertices
            for (int j = 0; j < 3; j++) {
              if (from_vertex(h[j]) == v[0])
                sh[0] = h[j];
              else if (from_vertex(h[j]) == v[1])
                sh[1] = h[j];
              else if (from_vertex(h[j]) == v[2])
                sh[2] = h[j];
              else
                assert(0);
            }
            // Interpolate parametric location using barycentric coordinates
            Param_loc param_res;
            param_res.ptexIndex = param_loc[h[0]].ptexIndex;
            assert(param_loc[h[0]].ptexIndex == param_loc[h[1]].ptexIndex &&
                   param_loc[h[0]].ptexIndex == param_loc[h[2]].ptexIndex);
            param_res.uv = param_loc[sh[1]].uv * uv(0) +
                           param_loc[sh[2]].uv * uv(1) +
                           param_loc[sh[0]].uv * (1.f - uv(0) - uv(1));
            Vector3f pos_res, normal_res;
            subdiv->evaluateLimit(param_res, pos_res, normal_res);
            V2.col(i) = pos_res;
            facing[i] = getFacing(view_pos - pos_res, normal_res);
          }
        }
      });

  parallel_for(
      tbb::blocked_range<uint32_t>(0u, (uint32_t)F2.rows(), GRAIN_SIZE),
      [&](const tbb::blocked_range<uint32_t> &range) {
        for (uint32_t i = range.begin(); i != range.end(); ++i) {
          const Vector3i &f = F2.row(i);
          FacingType VBO = UNDEFINED;
          for (int j = 0; j < 3; j++) {
            if (facing[f(j)] == UNDEFINED)
              continue;
            if (VBO != UNDEFINED && facing[f(j)] != VBO) {
              VBO = UNDEFINED;
              break;
            }
            VBO = facing[f(j)];
          }
          // compute consistency
          Vector3f a = V2.col(f(1)) - V2.col(f(0));
          Vector3f b = V2.col(f(2)) - V2.col(f(0));
          Vector3f normal = a.cross(b).normalized();
          real_t ndotv = normal.dot((view_pos - V2.col(f(0))).normalized());
          if ((ndotv > 0 && VBO != FRONT) || (ndotv < 0 && VBO != BACK)) {
            std::cout << "inconsistent" << std::endl;
            ray_intersections.push_back(
                (V2.col(f(0)) + V2.col(f(1)) + V2.col(f(2))) / 3.f);
          }
        }
      });
}

Matrix3Xf Mesh::get_contours(ContourMode mode, bool visible_only) {
  if (mode == INTERPOLATED_CONTOUR) {
    if (visible_only)
      return m_interpolated_contours.topLeftCorner(3, m_visible);
    return get_interpolated_contours();
  } else {
    auto vpositions = get_vertex_property<Vector3f>("v:point");

    Matrix3Xf contours = Matrix3Xf(3, m_chain_indices.size());
    for (size_t i = 0; i < m_chain_indices.size(); ++i) {
      Vertex v = Vertex(m_chain_indices[i]);
      contours.col(i) = vpositions[v];
    }

    return contours;
  }
  return Matrix3Xf();
}

Matrix3Xf Mesh::get_boundaries(bool visible_only) {
  auto vpositions = get_vertex_property<Vector3f>("v:point");

  Matrix3Xf boundaries = Matrix3Xf(3, m_boundary_indices.size());
  for (size_t i = 0; i < m_boundary_indices.size(); ++i) {
    Vertex v = Vertex(m_boundary_indices[i]);
    boundaries.col(i) = vpositions[v];
  }

  return boundaries;
}

MapXf Mesh::get_debug_points(nature_t nature) {
  using namespace VertexNature;

  m_debug_points.clear();
  m_point_colors.clear();
  m_2D_intersections.first = -1;

  // Points based on the view graph
  for (auto sv : m_svertices) {
    if (!sv)
      continue;
    nature_t sv_nature = sv->m_nature;
    if (nature & INTERSECTION_3D && sv_nature & INTERSECTION_3D) {
      m_debug_points.push_back(sv->m_pos3D);
      m_point_colors.push_back(intersection3DColor.head<3>().cast<real_t>());
    }
    if (nature & INTERSECTION_2D && sv_nature & INTERSECTION_2D) {
      if (m_2D_intersections.first == -1)
        m_2D_intersections.first = m_debug_points.size();
      m_debug_points.push_back(sv->m_pos3D);
      m_point_colors.push_back(intersection2DColor.head<3>().cast<real_t>());
    }
    if (nature & BOUNDARY_CURTAIN_FOLD && sv_nature & BOUNDARY_CURTAIN_FOLD) {
      m_debug_points.push_back(sv->m_pos3D);
      m_point_colors.push_back(curtainFoldColor.head<3>().cast<real_t>());
    }
    if (nature & CONTOUR_CURTAIN_FOLD && sv_nature & CONTOUR_CURTAIN_FOLD) {
      m_debug_points.push_back(sv->m_pos3D);
      m_point_colors.push_back(curtainFoldColor.head<3>().cast<real_t>());
    }
    if (nature & BIFURCATION && sv_nature & BIFURCATION) {
      m_debug_points.push_back(sv->m_pos3D);
      m_point_colors.push_back(bifurcationColor.head<3>().cast<real_t>());
    }
  }

  // 3D intersections
  if (nature & RAY_INTERSECTION) {
    m_debug_points.reserve(m_debug_points.size() + ray_intersections.size());
    m_point_colors.reserve(m_point_colors.size() + ray_intersections.size());
    for (size_t i = 0; i < ray_intersections.size(); i++) {
      m_debug_points.push_back(ray_intersections.at(i));
      m_point_colors.push_back(Vector3f(0.9, 0.0, 0.0));
    }
  }

  // Curtain folds in the interpolated contour mode
  if (nature & CONTOUR_CURTAIN_FOLD && !m_interpolated_contour_faces.empty()) {
    std::vector<Vector3f> curtain_folds;
    get_curtain_folds_interpolated(*this, curtain_folds);

    m_debug_points.insert(m_debug_points.end(), curtain_folds.begin(),
                          curtain_folds.end());
    m_point_colors.insert(m_point_colors.end(), curtain_folds.size(),
                          curtainFoldColor.head<3>().cast<real_t>());
  }

  // Vertices of the front portion adjacent to the contour edges
  if (nature & CONTOUR_FRONT) {
    auto is_contour_front_vertex = get_vertex_property<bool>("v:contour_front");
    auto vpositions = get_vertex_property<Vector3f>("v:point");

    if (!is_contour_front_vertex)
      return MapXf(m_debug_points[0].data(), 3, m_debug_points.size());

    m_debug_points.reserve(m_debug_points.size() + m_contour_edges.size());
    m_point_colors.reserve(m_point_colors.size() + m_contour_edges.size());
    for (Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit) {
      if (is_contour_front_vertex[*vit]) {
        m_debug_points.push_back(vpositions[*vit]);
        m_point_colors.push_back(frontVertexColor.head<3>().cast<real_t>());
      }
    }
  }

  return MapXf(m_debug_points[0].data(), 3, m_debug_points.size());
}

void Mesh::savePLYFile(const std::string &outputFilename, ContourMode mode,
                       int id) {
  // ---- output the PLY header ----
  FILE *fp = fopen(outputFilename.c_str(), "wt");

  if (fp == NULL) {
    printf("ERROR: CANNOT OPEN OUTPUT PLY FILE\n");
    exit(1);
  }

  fprintf(fp, "ply\n");
  fprintf(fp, "format ascii 1.0\n");
  fprintf(fp, "comment %s\n",
          (mode == MESH_CONTOUR || mode == VBO_CONTOUR) ? "mesh silhouettes"
                                                        : "smooth silhouettes");
  fprintf(fp, "element vertex %d\n", n_vertices());
  fprintf(fp, "property float x\n");
  fprintf(fp, "property float y\n");
  fprintf(fp, "property float z\n");
  fprintf(fp, "property float nx\n");
  fprintf(fp, "property float ny\n");
  fprintf(fp, "property float nz\n");
  fprintf(fp, "property float red\n");
  fprintf(fp, "property float green\n");
  fprintf(fp, "property float blue\n");
  fprintf(fp, "property float ndotv\n");

  auto patchID = get_face_property<int>("f:patchID");
  assert(patchID);
  if (id != -1) {
    fprintf(
        fp, "element face %d\n",
        int(std::count(patchID.vector().begin(), patchID.vector().end(), id)));
  } else {
    fprintf(fp, "element face %d\n", n_faces());
  }
  fprintf(fp, "property list uchar int vertex_index\n");
  fprintf(fp, "property uchar int\n"); // vbf
  fprintf(fp, "end_header\n");

  // ---- output all the vertices and save their IDs ----

  int nextVertID = 0;
  std::map<Vertex, int> vmapcc;
  auto vpositions = get_vertex_property<Vector3f>("v:point");
  auto vnormals = get_vertex_property<Vector3f>("v:normal");

  for (auto vit = vertices_begin(); vit != vertices_end(); ++vit) {
    vmapcc[*vit] = nextVertID;
    nextVertID++;
    const Vector3f &p = vpositions[*vit];
    fprintf(fp, "%.16f %.16f %.16f", p.x(), p.y(), p.z());

    Vector3f normal;
    if (mode == MESH_CONTOUR || mode == VBO_CONTOUR) {
      normal = vnormals[*vit];
    } else {
      normal = compute_vertex_normal(*vit);
    }
    fprintf(fp, " %.16f %.16f %.16f", double(normal[0]), double(normal[1]),
            double(normal[2]));
    fprintf(fp, " 0.0 0.0 0.0 0.0\n");
  }

  // --- output all the faces ----------
  auto VBO = get_face_property<FacingType>("f:VBO");
  Vertex v[3];
  for (auto fit = faces_begin(); fit != faces_end(); ++fit) {
    if (id != -1 && id != patchID[*fit])
      continue;

    FacingType vf = VBO[*fit];
    int vfint = (vf == FRONT ? 1 : (vf == BACK ? 2 : 3));

    verticesOfFace(*fit, v);

    fprintf(fp, "3 %d %d %d %d\n", vmapcc[v[0]], vmapcc[v[1]], vmapcc[v[2]],
            vfint);
  }

  // close output file
  fclose(fp);
}
