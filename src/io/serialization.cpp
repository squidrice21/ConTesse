// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include <assert.h>
#include <spdlog/fmt/ostr.h>
#include <surface_mesh/surface_mesh.h>

#include "camera.h"
#include "common.h"
#include "io/serialization.h"
#include "testing.h"

using namespace surface_mesh;

namespace serialization {

// ==================[ Basic Types ] ===================

void serialize(json_t &j, const double &v) { j = v; }
void deserialize(double &v, const json_t &j) { v = j; }
//
void serialize(json_t &j, const float &v) { j = v; }
void deserialize(float &v, const json_t &j) { v = j; }
//
void serialize(json_t &j, const int &v) { j = v; }
void deserialize(int &v, const json_t &j) { v = j; }
//
void serialize(json_t &j, const unsigned int &v) { j = v; }
void deserialize(unsigned int &v, const json_t &j) { v = j; }
//
void serialize(json_t &j, const std::string &v) { j = v; }
void deserialize(std::string &v, const json_t &j) { v = j; }

void serialize(json_t &j, const bool &v) { j = v; }
void deserialize(bool &v, const json_t &j) { v = j; }

void serialize(json_t &j, const std::filesystem::path &v) {
  j = v.generic_string();
}
void deserialize(std::filesystem::path &v, const json_t &j) {
  v = std::filesystem::path(std::string(j));
}

// ==================[ Surface_mesh ] ===================

// Helper for serializing properties of mesh
namespace {
template <typename T>
std::vector<T> get_vertex_property(const Surface_mesh &mesh,
                                   const std::string name) {
  std::vector<T> ans;
  Surface_mesh::Vertex_property<T> prop = mesh.get_vertex_property<T>(name);
  if (prop) {
    for (Surface_mesh::Vertex_iterator vit = mesh.vertices_begin();
         vit != mesh.vertices_end(); ++vit) {
      ans.push_back(prop[*vit]);
    }
  }
  return ans;
}

template <typename T>
std::vector<T> get_face_property(const Surface_mesh &mesh,
                                 const std::string name) {
  std::vector<T> ans;
  Surface_mesh::Face_property<T> prop = mesh.get_face_property<T>(name);
  if (prop) {
    for (Surface_mesh::Face_iterator fit = mesh.faces_begin();
         fit != mesh.faces_end(); ++fit) {
      ans.push_back(prop[*fit]);
    }
  }
  return ans;
}

template <typename T>
std::vector<T> get_edge_property(const Surface_mesh &mesh,
                                 const std::string name) {
  std::vector<T> ans;
  Surface_mesh::Edge_property<T> prop = mesh.get_edge_property<T>(name);
  if (prop) {
    for (Surface_mesh::Edge_iterator eit = mesh.edges_begin();
         eit != mesh.edges_end(); ++eit) {
      ans.push_back(prop[*eit]);
    }
  }
  return ans;
}

std::vector<Eigen::VectorXi> get_face_connectivity(const Surface_mesh &mesh) {
  std::vector<Eigen::VectorXi> ans;
  for (Surface_mesh::Face_iterator fit = mesh.faces_begin();
       fit != mesh.faces_end(); ++fit) {
    Surface_mesh::Vertex_around_face_circulator fvit = mesh.vertices(*fit),
                                                fvend = fvit;
    int num_verts = 0;
    do {
      ++num_verts;
    } while (++fvit != fvend);
    // logger().debug("{}", num_verts);
    //
    ans.emplace_back(num_verts);
    fvit = mesh.vertices(*fit), fvend = fvit;
    num_verts = 0;
    do {
      ans.back()[num_verts] = (*fvit).idx();
      ++num_verts;
    } while (++fvit != fvend);
    // logger().debug("{}",  ans.data.back());
  }
  return ans;
}

std::vector<Eigen::Vector2i> get_edge_connectivity(const Surface_mesh &mesh) {
  std::vector<Eigen::Vector2i> ans;
  for (Surface_mesh::Edge_iterator eit = mesh.edges_begin();
       eit != mesh.edges_end(); ++eit) {
    ans.emplace_back(mesh.vertex(*eit, 0).idx(), mesh.vertex(*eit, 1).idx());
  }
  return ans;
}

} // namespace

void serialize(json_t &json, const surface_mesh::Surface_mesh &mesh) {
  json["properties"] = json_t::object();
  serialize(json["properties"]["v:point"],
            get_vertex_property<Point>(mesh, "v:point"));
  serialize(json["properties"]["v:ndotv"],
            get_vertex_property<real_t>(mesh, "v:ndotv"));
  serialize(json["properties"]["v:normal"],
            get_vertex_property<Vector3f>(mesh, "v:normal"));
  serialize(json["properties"]["v:facing"],
            get_vertex_property<FacingType>(mesh, "v:facing"));
  if (mesh.get_vertex_property<real_t>("v:kappa"))
    serialize(json["properties"]["v:kappa"],
              get_vertex_property<real_t>(mesh, "v:kappa"));
  if (mesh.get_vertex_property<real_t>("v:kappa1"))
    serialize(json["properties"]["v:kappa1"],
              get_vertex_property<real_t>(mesh, "v:kappa1"));
  if (mesh.get_vertex_property<real_t>("v:kappa2"))
    serialize(json["properties"]["v:kappa2"],
              get_vertex_property<real_t>(mesh, "v:kappa2"));
  if (mesh.get_vertex_property<real_t>("v:kappa_fd"))
    serialize(json["properties"]["v:kappa_fd"],
              get_vertex_property<real_t>(mesh, "v:kappa_fd"));
  if (mesh.get_vertex_property<real_t>("v:kappa_ana"))
    serialize(json["properties"]["v:kappa_ana"],
              get_vertex_property<real_t>(mesh, "v:kappa_ana"));

  if (mesh.get_vertex_property<Vector3f>("v:d2x_dw2"))
    serialize(json["properties"]["v:d2x_dw2"],
              get_vertex_property<Vector3f>(mesh, "v:d2x_dw2"));
  if (mesh.get_vertex_property<Vector3f>("v:ds_fd"))
    serialize(json["properties"]["v:ds_fd"],
              get_vertex_property<Vector3f>(mesh, "v:ds_fd"));
  if (mesh.get_vertex_property<Vector3f>("v:dt_fd"))
    serialize(json["properties"]["v:dt_fd"],
              get_vertex_property<Vector3f>(mesh, "v:dt_fd"));
  if (mesh.get_vertex_property<Vector3f>("v:ds_osd"))
    serialize(json["properties"]["v:ds_osd"],
              get_vertex_property<Vector3f>(mesh, "v:ds_osd"));
  if (mesh.get_vertex_property<Vector3f>("v:dt_osd"))
    serialize(json["properties"]["v:dt_osd"],
              get_vertex_property<Vector3f>(mesh, "v:dt_osd"));

  if (mesh.get_vertex_property<Vector3f>("v:dsds_fd"))
    serialize(json["properties"]["v:dsds_fd"],
              get_vertex_property<Vector3f>(mesh, "v:dsds_fd"));
  if (mesh.get_vertex_property<Vector3f>("v:dsdt_fd"))
    serialize(json["properties"]["v:dsdt_fd"],
              get_vertex_property<Vector3f>(mesh, "v:dsdt_fd"));
  if (mesh.get_vertex_property<Vector3f>("v:dtdt_fd"))
    serialize(json["properties"]["v:dtdt_fd"],
              get_vertex_property<Vector3f>(mesh, "v:dtdt_fd"));
  if (mesh.get_vertex_property<Vector3f>("v:dsds_osd"))
    serialize(json["properties"]["v:dsds_osd"],
              get_vertex_property<Vector3f>(mesh, "v:dsds_osd"));
  if (mesh.get_vertex_property<Vector3f>("v:dsdt_osd"))
    serialize(json["properties"]["v:dsdt_osd"],
              get_vertex_property<Vector3f>(mesh, "v:dsdt_osd"));
  if (mesh.get_vertex_property<Vector3f>("v:dtdt_osd"))
    serialize(json["properties"]["v:dtdt_osd"],
              get_vertex_property<Vector3f>(mesh, "v:dtdt_osd"));

  if (mesh.get_vertex_property<bool>("v:cusp"))
    serialize(json["properties"]["v:cusp"],
              get_vertex_property<bool>(mesh, "v:cusp"));
  if (mesh.get_vertex_property<FacingType>("v:cusp_facing"))
    serialize(json["properties"]["v:cusp_facing"],
              get_vertex_property<FacingType>(mesh, "v:cusp_facing"));
  if (mesh.get_vertex_property<bool>("v:sim_cusp"))
    serialize(json["properties"]["v:sim_cusp"],
              get_vertex_property<bool>(mesh, "v:sim_cusp"));
  if (mesh.get_vertex_property<FacingType>("v:sim_cusp_facing"))
    serialize(json["properties"]["v:sim_cusp_facing"],
              get_vertex_property<FacingType>(mesh, "v:sim_cusp_facing"));
  if (mesh.get_vertex_property<Vector3f>("v:contour_laplacian"))
    serialize(json["properties"]["v:contour_laplacian"],
              get_vertex_property<Vector3f>(mesh, "v:contour_laplacian"));
  if (mesh.get_vertex_property<Vector3f>("v:view"))
    serialize(json["properties"]["v:view"],
              get_vertex_property<Vector3f>(mesh, "v:view"));
  if (mesh.get_vertex_property<bool>("v:extraordinary"))
    serialize(json["properties"]["v:extraordinary"],
              get_vertex_property<bool>(mesh, "v:extraordinary"));
  if (mesh.get_vertex_property<Vector3f>("v:ptex"))
    serialize(json["properties"]["v:ptex"],
              get_vertex_property<Vector3f>(mesh, "v:ptex"));
  if (mesh.get_vertex_property<bool>("v:is_inflated"))
    serialize(json["properties"]["v:is_inflated"],
              get_vertex_property<bool>(mesh, "v:is_inflated"));
  if (mesh.get_vertex_property<int>("v:orig_idx"))
    serialize(json["properties"]["v:orig_idx"],
              get_vertex_property<int>(mesh, "v:orig_idx"));
  if (mesh.get_vertex_property<int>("v:orig_f_idx"))
    serialize(json["properties"]["v:orig_f_idx"],
              get_vertex_property<int>(mesh, "v:orig_f_idx"));
  if (mesh.get_vertex_property<int>("v:intersection_2d"))
    serialize(json["properties"]["v:intersection_2d"],
              get_vertex_property<int>(mesh, "v:intersection_2d"));
  if (mesh.get_vertex_property<bool>("v:is_valid_intersection_2d"))
    serialize(json["properties"]["v:is_valid_intersection_2d"],
              get_vertex_property<bool>(mesh, "v:is_valid_intersection_2d"));
  if (mesh.get_vertex_property<bool>("v:stationary"))
    serialize(json["properties"]["v:stationary"],
              get_vertex_property<bool>(mesh, "v:stationary"));
  if (mesh.get_vertex_property<bool>("v:cut_branch"))
    serialize(json["properties"]["v:cut_branch"],
              get_vertex_property<bool>(mesh, "v:cut_branch"));
  if (mesh.get_vertex_property<ContourVertexType>("v:contour_type"))
    serialize(json["properties"]["v:contour_type"],
              get_vertex_property<ContourVertexType>(mesh, "v:contour_type"));
  if (mesh.get_vertex_property<int>("v:comp_idx"))
    serialize(json["properties"]["v:comp_idx"],
              get_vertex_property<int>(mesh, "v:comp_idx"));
  if (mesh.get_vertex_property<int>("v:sew_component"))
    serialize(json["properties"]["v:sew_component"],
              get_vertex_property<int>(mesh, "v:sew_component"));

  serialize(json["properties"]["e:connectivity"], get_edge_connectivity(mesh));
  if (mesh.get_edge_property<bool>("e:concave"))
    serialize(json["properties"]["e:concave"],
              get_edge_property<bool>(mesh, "e:concave"));
  if (mesh.get_edge_property<bool>("e:concave_kr"))
    serialize(json["properties"]["e:concave_kr"],
              get_edge_property<bool>(mesh, "e:concave_kr"));
  if (mesh.get_edge_property<bool>("e:concave_consist"))
    serialize(json["properties"]["e:concave_consist"],
              get_edge_property<bool>(mesh, "e:concave_consist"));
  serialize(json["properties"]["e:contour"],
            get_edge_property<real_t>(mesh, "e:contour"));
  if (mesh.get_edge_property<bool>("e:is_boundary"))
    serialize(json["properties"]["e:is_boundary"],
              get_edge_property<bool>(mesh, "e:is_boundary"));
  if (mesh.get_edge_property<bool>("e:vis_mesh_contour"))
    serialize(json["properties"]["e:vis_mesh_contour"],
              get_edge_property<bool>(mesh, "e:vis_mesh_contour"));
  if (mesh.get_edge_property<real_t>("e:cdott"))
    serialize(json["properties"]["e:cdott"],
              get_edge_property<real_t>(mesh, "e:cdott"));
  if (mesh.get_edge_property<Vector4f>("e:loop"))
    serialize(json["properties"]["e:loop"],
              get_edge_property<Vector4f>(mesh, "e:loop"));
  if (mesh.get_edge_property<int>("e:disk_cut"))
    serialize(json["properties"]["e:disk_cut"],
              get_edge_property<int>(mesh, "e:disk_cut"));
  if (mesh.get_edge_property<int>("e:disk_cut2"))
    serialize(json["properties"]["e:disk_cut2"],
              get_edge_property<int>(mesh, "e:disk_cut2"));
  if (mesh.get_edge_property<int>("e:disk_cut_debug"))
    serialize(json["properties"]["e:disk_cut_debug"],
              get_edge_property<int>(mesh, "e:disk_cut_debug"));
  if (mesh.get_edge_property<bool>("e:is_inflated_edge"))
    serialize(json["properties"]["e:is_inflated_edge"],
              get_edge_property<bool>(mesh, "e:is_inflated_edge"));
  if (mesh.get_edge_property<BoundaryType>("e:boundary_type"))
    serialize(json["properties"]["e:boundary_type"],
              get_edge_property<BoundaryType>(mesh, "e:boundary_type"));
  if (mesh.get_edge_property<int>("e:feasibility_collapsing"))
    serialize(json["properties"]["e:feasibility_collapsing"],
              get_edge_property<int>(mesh, "e:feasibility_collapsing"));
  if (mesh.get_edge_property<bool>("e:subdiv_contour"))
    serialize(json["properties"]["e:subdiv_contour"],
              get_edge_property<bool>(mesh, "e:subdiv_contour"));
  if (mesh.get_edge_property<bool>("e:cut_candidate"))
    serialize(json["properties"]["e:cut_candidate"],
              get_edge_property<bool>(mesh, "e:cut_candidate"));
  if (mesh.get_edge_property<int>("e:qi"))
    serialize(json["properties"]["e:qi"], get_edge_property<int>(mesh, "e:qi"));
  if (mesh.get_edge_property<bool>("e:patch_removed"))
    serialize(json["properties"]["e:patch_removed"],
              get_edge_property<bool>(mesh, "e:patch_removed"));

  serialize(json["properties"]["f:connectivity"], get_face_connectivity(mesh));
  serialize(json["properties"]["f:normal"],
            get_face_property<Vector3f>(mesh, "f:normal"));
  serialize(json["properties"]["f:patchID"],
            get_face_property<int>(mesh, "f:patchID"));
  serialize(json["properties"]["f:VBO"],
            get_face_property<FacingType>(mesh, "f:VBO"));
  serialize(json["properties"]["f:VBO_f"],
            get_face_property<FacingType>(mesh, "f:VBO_f"));
  serialize(json["properties"]["f:concave"],
            get_face_property<bool>(mesh, "f:concave"));
  if (mesh.get_face_property<bool>("f:cusp"))
    serialize(json["properties"]["f:cusp"],
              get_face_property<bool>(mesh, "f:cusp"));
  if (mesh.get_face_property<Vector3f>("f:ptex"))
    serialize(json["properties"]["f:ptex"],
              get_face_property<Vector3f>(mesh, "f:ptex"));
  if (mesh.get_face_property<real_t>("f:kappa_fd"))
    serialize(json["properties"]["f:kappa_fd"],
              get_face_property<real_t>(mesh, "f:kappa_fd"));
  if (mesh.get_face_property<Vector3f>("f:view"))
    serialize(json["properties"]["f:view"],
              get_face_property<Vector3f>(mesh, "f:view"));
  if (mesh.get_face_property<Vector3f>("f:view_tangent"))
    serialize(json["properties"]["f:view_tangent"],
              get_face_property<Vector3f>(mesh, "f:view_tangent"));
  if (mesh.get_face_property<Vector3f>("f:view_proj"))
    serialize(json["properties"]["f:view_proj"],
              get_face_property<Vector3f>(mesh, "f:view_proj"));
  if (mesh.get_face_property<Vector3f>("f:ds"))
    serialize(json["properties"]["f:ds"],
              get_face_property<Vector3f>(mesh, "f:ds"));
  if (mesh.get_face_property<Vector3f>("f:ds_orth"))
    serialize(json["properties"]["f:ds_orth"],
              get_face_property<Vector3f>(mesh, "f:ds_orth"));

  if (mesh.get_face_property<int>("f:near_contour"))
    serialize(json["properties"]["f:near_contour"],
              get_face_property<int>(mesh, "f:near_contour"));
  if (mesh.get_face_property<int>("f:disk"))
    serialize(json["properties"]["f:disk"],
              get_face_property<int>(mesh, "f:disk"));
  if (mesh.get_face_property<int>("f:patch_component"))
    serialize(json["properties"]["f:patch_component"],
              get_face_property<int>(mesh, "f:patch_component"));
  if (mesh.get_face_property<int>("f:sew_component"))
    serialize(json["properties"]["f:sew_component"],
              get_face_property<int>(mesh, "f:sew_component"));
}

// ==================[ Camera ] ===================

void serialize(json_t &j, const Camera &v) {
  // Ensure that it is up to date...
  v.updateProjectionMatrix();
  v.updateViewMatrix();

  // Now save
  j = json_t::object();

  // View port
  serialize(j["mVpX"], v.mVpX);
  serialize(j["mVpY"], v.mVpY);
  serialize(j["mVpWidth"], v.mVpWidth);
  serialize(j["mVpHeight"], v.mVpHeight);

  // Projection
  serialize(j["mFovY"], v.mFovY);
  serialize(j["mNearDist"], v.mNearDist);
  serialize(j["mFarDist"], v.mFarDist);

  // View Matrix
  serialize(j["mViewMatrix"], v.mViewMatrix.matrix());
}

// ==================[ SerializationObject ] ===================

void serialize(json_t &json, const SerializationObject &serialization_object) {
  json = json_t::object();

  if (serialization_object.mesh) {
    serialize(json["m_mesh"], *serialization_object.mesh);
  }

  if (serialization_object.camera) {
    serialize(json["m_contour_camera"], *serialization_object.camera);
    // Write an array of viewer cameras in case we want to add more camera
    // locations
    {
      std::vector<Camera> tmp(1, *serialization_object.camera);
      serialize(json["m_view_cameras"], tmp);
    }
  }

  // Pack rays in an array
  serialize(json["m_ray"], serialization_object.rays);
}

// ==================[ SerializationRay ] ===================

void serialize(json_t &json, const SerializationRay &serialization_ray) {
  json = json_t::object();

  serialize(json["origin"], serialization_ray.origin);
  serialize(json["dir"], serialization_ray.dir);
  serialize(json["query_point"], serialization_ray.query_point);
  serialize(json["query_point_epsilon"], serialization_ray.query_point_epsilon);
  serialize(json["query_edge"], serialization_ray.query_edge);
  serialize(json["intersection_points"], serialization_ray.intersection_points);
  serialize(json["intersection_faces"], serialization_ray.intersection_faces);
}

// ==================[ Testing ] ===================

void serialize(json_t &j, const BasicCameraCase &c) {
  j = json_t::object();
  serialize(j["camera_x"], c.camera_x);
  serialize(j["camera_y"], c.camera_y);
  serialize(j["camera_z"], c.camera_z);
  serialize(j["full_info"], c.full_info);
  serialize(j["lookat"], c.lookat);
}
void deserialize(BasicCameraCase &c, const json_t &j) {
  deserialize(c.camera_x, j["camera_x"]);
  deserialize(c.camera_y, j["camera_y"]);
  deserialize(c.camera_z, j["camera_z"]);
  deserialize(c.full_info, j["full_info"]);
  deserialize(c.lookat, j["lookat"]);
}

void serialize(json_t &j, const FileCameraCase &c) {
  j = json_t::object();
  serialize(j["path"], c.path);
  serialize(j["frame_idx"], static_cast<unsigned int>(c.frame_idx));
}
void deserialize(FileCameraCase &c, const json_t &j) {
  deserialize(c.path, j["path"]);
  unsigned int frame_idx;
  deserialize(frame_idx, j["frame_idx"]);
  c.frame_idx = frame_idx;
}

void serialize(json_t &j, const TestCase &c) {
  j = json_t::object();
  serialize(j["model_name"], c.model_name);
  serialize(j["model_path"], c.model_path);
  serialize(j["description"], c.description);
  serialize(j["camera_info"], c.camera_info);
  serialize(j["subdiv_level"], c.subdiv_level);
}

void deserialize(TestCase &c, const json_t &j) {
  deserialize(c.model_name, j["model_name"]);
  deserialize(c.model_path, j["model_path"]);
  deserialize(c.description, j["description"]);
  deserialize(c.camera_info, j["camera_info"]);
  deserialize(c.subdiv_level, j["subdiv_level"]);
}

void serialize(json_t &j, const TestingSingleton &c) {
  j = json_t::object();
  serialize(j["multicamera_test_cases"], c.multicamera_test_cases);
  serialize(j["pipeline_test_cases"], c.pipeline_test_cases);
}

void deserialize(TestingSingleton &c, const json_t &j) {
  deserialize(c.multicamera_test_cases, j["multicamera_test_cases"]);
  deserialize(c.pipeline_test_cases, j["pipeline_test_cases"]);
}

} // namespace serialization