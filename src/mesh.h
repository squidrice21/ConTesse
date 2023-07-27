// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#ifndef _MESH_H
#define _MESH_H

#include <memory>
#include <vector>

#include <igl/predicates/predicates.h>
#include <surface_mesh/surface_mesh.h>

#include "bvh.h"
#include "camera.h"
#include "subdiv_common.h"
#include "viewgraph.h"

typedef std::vector<Halfedge> Chain;

struct VirtualHalfedge {
  Halfedge he;
  Vertex from_v, to_v;
};
typedef std::vector<VirtualHalfedge> SimplifiedChain;

class Mesh : public surface_mesh::Surface_mesh {

public:
  Mesh();
  Mesh(const std::string &filename);

  ~Mesh();

  bool load(const std::string &filename);
  void savePLYFile(const std::string &outputFilename, ContourMode mode,
                   int id = -1);

  void init();

  const Vector3f get_mesh_center() const { return m_bbox.center(); }

  real_t get_dist_max() const { return m_bbox.diagonal().norm(); }

  /// returns face indices as a 3xN matrix of integers
  MapXu get_indices() { return MapXu(m_indices.data(), 3, n_faces()); }

  MapXu get_sorted_indices() {
    return MapXu(m_sorted_indices.data(), 3, n_faces());
  }

  MapXu get_patch_indices() {
    return MapXu(m_patch_indices.data(), 3, n_faces());
  }

  MapConstXu get_const_patch_indices() const {
    return MapConstXu(m_patch_indices.data(), 3, n_faces());
  }

  void update_patch_facing(const Mesh &mesh);
  FacingType get_patch_facing(size_t patch_index) const;

  // Accessors to vertex attributes as Eigen's matrices:
  MapXf get_positions() {
    auto &vertices = get_vertex_property<Vector3f>("v:point").vector();
    return MapXf(vertices[0].data(), 3, vertices.size());
  }

  MapXf get_colors() {
    auto &vcolors = get_vertex_property<Vector3f>("v:color").vector();
    return MapXf(vcolors[0].data(), 3, vcolors.size());
  }

  MapXf get_normals() {
    auto &vnormals = get_vertex_property<Vector3f>("v:normal").vector();
    return MapXf(vnormals[0].data(), 3, vnormals.size());
  }

  const Matrix3Xf &get_normal_segments() { return m_normal_segments; }

  /// Re-compute the aligned bounding box (needs to be called after editing
  /// vertex positions)
  void updateBoundingBox();

  const AABB &boundingBox() const { return m_bbox; }

  const BVH *get_bvh() const { return m_bvh; }
  void build_bvh();

  Subdiv &subdivision() { return *m_subdiv; }
  Subdiv const &const_subdivision() const { return *m_subdiv; }

  void set_subdivision(std::shared_ptr<Subdiv> input) { m_subdiv = input; }

  /*============================================*/
  /* NPR algorithms */

  void tagConcaveEdges(ContourMode cont_mode, ConvexityMode mode,
                       Camera const &camera);

  void computeConsistency(const Vector3f &view_pos);
  void computeConsistencySubdivided(const Vector3f &view_pos);

  int splitEdgesWithZeroCrossing(const Vector3f &view_pos,
                                 const Subdiv *subdiv);

  int num_front_faces() const { return m_front / 3; }
  int num_back_faces() const { return m_back / 3; }

  std::vector<int> &get_patch_lengths() { return m_patch_lengths; }
  std::vector<int> const &get_const_patch_lengths() const {
    return m_patch_lengths;
  }
  std::vector<std::vector<std::vector<Vertex>>> get_patch_boundaries() const {
    return m_patch_boundaries;
  }

  MapXu get_inconsistent_indices() {
    return MapXu(m_inconsistent_indices.data(), 3,
                 m_inconsistent_indices.size() / 3);
  }

  void extractContours(ContourMode mode, Camera const &camera);
  bool is_contour_edge(Edge const &e) const;

  void extractBoundaries();

  void extractBoundaryCurtainFolds(const Vector3f &view_pos);

  void extractSurfaceIntersections();

  void projectSVertices(const Matrix4f &model, const Matrix4f &proj,
                        const Vector2i &viewportSize);

  MatrixXu &get_mesh_contours_indices() { return m_mesh_contours_indices; }

  MapXu get_boundary_indices() {
    return MapXu(m_boundary_indices.data(), 1, m_boundary_indices.size());
  }

  MapConstXu get_seams_indices() const {
    return MapConstXu(m_seams_indices.data(), 1, m_seams_indices.size());
  }

  Matrix3Xf &get_interpolated_contours() { return m_interpolated_contours; }
  std::vector<int> &get_interpolated_contour_faces() {
    return m_interpolated_contour_faces;
  }
  Matrix3Xf const &get_const_interpolated_contours() const {
    return m_interpolated_contours_orig;
  }
  std::vector<int> const &get_const_interpolated_contour_faces() const {
    return m_interpolated_contour_faces;
  }

  MapXf get_debug_points(nature_t nature);

  MapXf get_point_colors() {
    return MapXf(m_point_colors[0].data(), 3, m_point_colors.size());
  }

  const std::pair<int, int> &num2Dintersections() const {
    return m_2D_intersections;
  }

  void insertContours(const Vector3f &view_pos, const Subdiv *subdiv,
                      bool shift_allowed);

  void computeContourChains(ContourMode mode);

  void computeContourVisibility(const Vector3f &view_pos, ContourMode mode,
                                const ProgressCallback &progress);

  std::vector<int> &get_chain_lengths() { return m_chain_lengths; }

  MapXu get_chain_indices() {
    return MapXu(m_chain_indices.data(), 1, m_chain_indices.size());
  }

  int num_visible_contours() const { return m_visible / 2; }

  Matrix3Xf get_contours(ContourMode mode, bool visible_only);

  Matrix3Xf get_boundaries(bool visible_only);

  std::vector<int> &get_boundaries_lengths() { return m_boundaries_lengths; }

  std::vector<int> &get_seams_lengths() { return m_seams_lengths; }
  std::vector<int> const &get_seams_lengths_const() const {
    return m_seams_lengths;
  }

  Halfedge get_patch_halfedge(Vertex const &v, size_t patch_id) const;

  Matrix3Xf &get_surface_intersections() { return m_surface_intersections; }

  std::vector<int> &get_surface_intersections_lengths() {
    return m_surface_intersections_lengths;
  }

  std::vector<std::list<Edge>> &get_chains() { return m_chains; }

  std::vector<std::shared_ptr<Chain>> &get_oriented_chains() {
    return m_oriented_chains;
  }
  std::vector<std::shared_ptr<Chain>> const &get_const_oriented_chains() const {
    return m_oriented_chains;
  }
  std::multimap<int, std::shared_ptr<Chain>> &get_patch_chains() {
    return m_patch_chains;
  }
  std::multimap<int, std::shared_ptr<Chain>> const &
  get_const_patch_chains() const {
    return m_patch_chains;
  }

  void add_simplified_chain(int patch_idx,
                            std::shared_ptr<SimplifiedChain> const &chain) {
    m_patch_simplified_chains.insert(std::make_pair(patch_idx, chain));
    m_simplified_chains.emplace_back(chain);
  }
  std::vector<std::shared_ptr<SimplifiedChain>> &get_simplified_chains() {
    return m_simplified_chains;
  }
  std::vector<std::shared_ptr<SimplifiedChain>> const &
  get_const_simplified_chains() const {
    return m_simplified_chains;
  }
  std::multimap<int, std::shared_ptr<SimplifiedChain>> const &
  get_const_patch_simplified_chains() const {
    return m_patch_simplified_chains;
  }

  void computeSweepLineIntersections(ContourMode mode, const Camera &camera,
                                     const Vector2i &viewport);

  bool hasBoundaries() const { return m_boundary_indices.size() > 0; }
  bool hasSurfaceIntersections() const {
    return m_surface_intersections.size() > 0;
  }

  void floodFill();

  Mesh *triangulateSimplePatches(std::vector<int> patchIDs,
                                 const Subdiv *subdiv, const Matrix4f &model,
                                 const Matrix4f &proj);

  void triangulateSimplePatches(ContourMode mode, const Subdiv *subdiv,
                                const Matrix4f &model, const Matrix4f &proj,
                                const Vector2i &viewportSize, MatrixXf &V2,
                                MatrixXi &F2);

  void extractPatchBoundaries(bool makeDisk);
  void markPatchBoundaryEdges(bool makeDisk);

  static void updateVBO(Mesh *mesh, const Vector3f &view_pos);

  // Mesh data structure
  void verticesOfFace(Face f, Vertex v[3]) const;
  void halfedgesOfFace(Face f, Halfedge h[3]) const;
  bool are_adjacent(Face f1, Face f2) const;

  // External functions operating on the mesh class
  friend void insert_interpolated_contours(Mesh &mesh,
                                           Vector3f const &view_pos);
  friend bool chain_contour(Mesh &mesh, Camera const &camera,
                            bool to_cache_chain, bool ff_only);
  friend bool chain_cut_graph(Mesh &mesh, Camera const &camera, bool ff_only);
  friend void tag_non_fish_tail(Mesh &mesh, Camera const &camera, bool ff_only);
  friend void simplify_contours(Mesh &mesh);
  friend void subdivide_contour_edges(Mesh &mesh, Camera const &camera,
                                      std::vector<Edge> &subdiv_edges);
  friend void subdivide_contour_edges_even(Mesh &mesh, Camera const &camera,
                                           std::vector<Edge> &subdiv_edges);
  real_t bisect_search(const Vector3f &view_pos, const Subdiv *subdiv,
                       Surface_mesh::Edge e, Vector3f &pos_res,
                       Vector3f &normal_res, Param_loc &uv_res);
  real_t bisect_search(const Vector3f &view_pos, const Subdiv *subdiv,
                       Param_loc const &param_v1, Param_loc const &param_v2,
                       Vector3f &pos_res, Vector3f &normal_res,
                       Param_loc &uv_res);

  // View graph
  std::vector<SVertex *> m_svertices;
  std::vector<FEdge *> m_fedges;

private:
  // NOTE: Need to call this function after topological update to the mesh.
  void updateIndicesFromVBO(const Vector3f &view_pos);

  // \returns true if d and e are on the same side of the abc triangle
  inline bool sameSide(const Vector3f &a, const Vector3f &b, const Vector3f &c,
                       const Vector3f &d, const Vector3f &e) {
    return (igl::predicates::orient3d(a, b, c, d) ==
            igl::predicates::Orientation::POSITIVE) ==
           (igl::predicates::orient3d(a, b, c, e) ==
            igl::predicates::Orientation::POSITIVE);
  }

  // \returns returns a positive value if d is on the front side of triangle abc
  inline double frontSide(const Vector3f &a, const Vector3f &b,
                          const Vector3f &c, const Vector3f &d) {
    return -static_cast<double>(igl::predicates::orient3d(a, b, c, d));
  }

  /// \returns true if c and d are on the same side of the line ab
  inline bool sameSide(const Vector2f &a, const Vector2f &b, const Vector2f &c,
                       const Vector2f &d) {
    return (igl::predicates::orient2d(a, b, c) ==
            igl::predicates::Orientation::POSITIVE) ==
           (igl::predicates::orient2d(a, b, d) ==
            igl::predicates::Orientation::POSITIVE);
  }

  typedef std::tuple<Surface_mesh::Edge, Vector3f, Vector3f, Param_loc, real_t>
      EdgeTuple;
  bool canShift(Vertex vsource, const EdgeTuple &etuple);

  Vertex splitEdge(const EdgeTuple &etuple, bool shift_allowed);

  AABB m_bbox;
  std::vector<uint32_t> m_indices;
  std::vector<uint32_t> m_sorted_indices;
  std::vector<uint32_t> m_inconsistent_indices;
  Matrix3Xf m_normal_segments;
  int m_front, m_back;

  BVH *m_bvh;

  // Mesh contours
  MatrixXu m_mesh_contours_indices;
  tbb::concurrent_vector<Edge> m_contour_edges;
  std::vector<std::list<Edge>> m_chains;
  std::vector<int> m_chain_lengths;
  std::vector<uint32_t> m_chain_indices;

  // Interpolated contours
  // m_interpolated_contours: 3 x (2x#interpolated contour edge)
  // m_interpolated_contour_faces: indices of the faces passed by the
  // interpolated contours
  Matrix3Xf m_interpolated_contours, m_interpolated_contours_orig;
  std::vector<int> m_interpolated_contour_faces;

  // Boundaries
  std::vector<uint32_t> m_boundary_indices;
  std::vector<int> m_boundaries_lengths;

  // Surface-Surface intersections
  Matrix3Xf m_surface_intersections;
  std::vector<int> m_surface_intersections_lengths;

  // Patches
  std::vector<std::vector<std::vector<Vertex>>> m_patch_boundaries;
  std::vector<FacingType> m_patch_facing;
  Mesh *res_mesh;
  std::vector<uint32_t> m_patch_indices;
  std::vector<int> m_patch_lengths;
  std::vector<uint32_t> m_seams_indices;
  std::vector<int> m_seams_lengths;

  // Debug
  tbb::concurrent_vector<Vector3f> ray_intersections;

  std::vector<Vector3f> m_debug_points, m_point_colors;
  std::pair<int, int> m_2D_intersections;
  int m_visible;

  // Pack the subdivision here
  std::shared_ptr<Subdiv> m_subdiv;

  // The oriented chain
  std::multimap<int, std::shared_ptr<Chain>> m_patch_chains;
  std::vector<std::shared_ptr<Chain>> m_oriented_chains;
  std::multimap<int, std::shared_ptr<SimplifiedChain>>
      m_patch_simplified_chains;
  std::vector<std::shared_ptr<SimplifiedChain>> m_simplified_chains;
};

#endif //_MESH_H
