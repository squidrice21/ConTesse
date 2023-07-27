#include <spdlog/fmt/ostr.h>

#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <unordered_set>

#include "foldertools.h"
#include "insert_interpolated_contours.h"
#include "insert_subdivision_cusps.h"
#include "logger.h"
#include "svg.h"
#include "testing.h"

std::string output_folder;

BasicCameraCase::BasicCameraCase(real_t a, real_t b)
    : camera_x(a), camera_y(b), camera_z(0), full_info(false),
      lookat(Vector3f::Zero()) {}

BasicCameraCase::BasicCameraCase(real_t a, real_t b, real_t c,
                                 Vector3f const &vec)
    : camera_x(a), camera_y(b), camera_z(c), full_info(true), lookat(vec) {}

FileCameraCase::FileCameraCase(std::filesystem::path p, std::size_t idx)
    : path(p), frame_idx(idx) {}

void TestingSingleton::load_cases(std::filesystem::path p) {
  if (!std::filesystem::exists(p)) {
    logger().error("File {} not found", p.string());
    return;
  }
  serialization::load(*this, p);
  spec_root_path = std::filesystem::absolute(p).parent_path();
}

std::filesystem::path TestingSingleton::repo_root_path() {
  return std::filesystem::path(__FILE__).parent_path().parent_path();
}

TestingSingleton &TestingSingleton::instance() {
  static TestingSingleton instance_;
  return instance_;
}

void setup_camera(Camera &camera, const Mesh &mesh, const TestCase &test_case) {
  std::visit(
      contess::overloaded{
          [&](const std::monostate &) -> void {
            contess_assert_msg(0, "Cameta is empty");
          },
          [&](const BasicCameraCase &camera_info) -> void {
            Vector3f cam_pos, lookat;
            if (!camera_info.full_info) {
              cam_pos = Vector3f(camera_info.camera_x, camera_info.camera_y,
                                 mesh.get_dist_max());
              lookat = mesh.get_mesh_center();
            } else {
              cam_pos = Vector3f(camera_info.camera_x, camera_info.camera_y,
                                 camera_info.camera_z);
              lookat = camera_info.lookat;
            }
            Vector3f up = Vector3f(0.0f, 1.0f, 0.0f);
            double cross = (lookat - cam_pos).cross(up).norm();
            up = (cross < std::numeric_limits<real_t>::epsilon())
                     ? Vector3f(0.0f, 0.0f, 1.0f)
                     : Vector3f(0.0f, 1.0f, 0.0f);
            camera.lookAt(cam_pos, lookat, up);
          },
          [&](const FileCameraCase &camera_info) -> void {
            int startframe, endframe;
            std::string cam_path = CONTOURS_TESSELATION_DATA_DIR "/cameras/" +
                                   camera_info.path.string() + ".cam";
            camera.load(cam_path, startframe, endframe);
            camera.getInterpolator().gotoPrevKeyFrame();
          }},
      test_case.camera_info);
}

std::string TestCase::get_short_name() const {
  std::string ans = this->model_name;

  ans += "_cam_";
  std::visit(
      contess::overloaded{
          [&](const std::monostate &) -> void {
            contess_assert_msg(0, "Camera is empty");
          },
          [&](const BasicCameraCase &camera_info) -> void {
            ans += to_string_with_precision<real_t>(camera_info.camera_x, 2) +
                   "_" +
                   to_string_with_precision<real_t>(camera_info.camera_y, 2) +
                   "_" +
                   to_string_with_precision<real_t>(camera_info.camera_z, 2);
          },
          [&](const FileCameraCase &camera_info) -> void {
            ans += camera_info.path.string() + "_" +
                   std::to_string(camera_info.frame_idx);
          }},
      this->camera_info);

  ans += "_subd_" + std::to_string(this->subdiv_level);

  return ans;
}

std::string TestCase::get_description() const {
  std::string enriched_description = this->description + ";";

  enriched_description += "Camera at (";

  std::visit(
      contess::overloaded{[&](const std::monostate &) -> void {
                            contess_assert_msg(0, "Camera is empty");
                          },
                          [&](const BasicCameraCase &camera_info) -> void {
                            enriched_description +=
                                std::to_string(camera_info.camera_x) + ", " +
                                std::to_string(camera_info.camera_y) + ", " +
                                std::to_string(camera_info.camera_z);
                          },
                          [&](const FileCameraCase &camera_info) -> void {
                            enriched_description +=
                                camera_info.path.string() + ", " +
                                std::to_string(camera_info.frame_idx);
                          }},
      this->camera_info);

  enriched_description += ");";

  enriched_description += "Subdiv: " + std::to_string(this->subdiv_level);

  return enriched_description;
}

std::string test_output_dir() {
  // Output to current dir.
  // Can be changed to something else like TEST_DATA_DIR+"/output"
  output_folder = foldertools::cwd() + "/test_out";
  return output_folder;
}

Mesh read_input(std::string model_filename, bool is_full_path, int subdiv_level,
                Subdiv::Backend sub_backend) {

  std::string input_filename =
      CONTOURS_TESSELATION_DATA_DIR "/meshes/" + model_filename;
  {
    std::ifstream f(input_filename);
    if (!f.is_open())
      input_filename = CONTOURS_TESSELATION_DATA_DIR "/" + model_filename;
  }
  if (is_full_path)
    input_filename = model_filename;

  // Check the file exists.
  {
    std::ifstream f(input_filename);
    if (!f.is_open())
      throw std::runtime_error(input_filename + " not found"); // catch macro
  }

  logger().info("Reading model: {} ", input_filename);
  Mesh mesh;

  if (subdiv_level <= 0) {
    mesh.load(input_filename);
    mesh.init();
  } else {
    std::shared_ptr<Subdiv> subdiv = std::make_shared<Subdiv>();
    subdiv->load(input_filename, sub_backend, subdiv_level);
    create_subdivided_mesh(subdiv, mesh);
    mesh.init();
  }

  return mesh;
}

void write_param_to_vertex(Mesh &mesh) {
  auto param_loc = mesh.get_halfedge_property<Param_loc>("h:param_loc");
  auto ptex = mesh.vertex_property<Vector3f>("v:ptex");
  if (!ptex) {
    ptex = mesh.add_vertex_property<Vector3f>("v:ptex");
  }
  ptex.vector().assign(ptex.vector().size(), -1 * Vector3f::Ones());

  for (uint32_t i = 0; i < mesh.n_vertices(); ++i) {
    Vertex v(i);
    auto hit = mesh.halfedges(v);

    // When the vertex is on the patch boundary, it takes an arbitrary parameter
    auto param = param_loc[*hit];
    ptex[v][0] = param.ptexIndex;

    ptex[v].tail(2) = param.uv;
  }
}

void write_facing_labelled_mesh(Mesh const &mesh, std::string name_pattern) {
  mesh.write(name_pattern + "_interpolated_contours.obj");
  logger().info("Wrote contour-inserted mesh to obj: {}",
                name_pattern + "_interpolated_contours.obj");

  Mesh mesh_front = mesh;
  mesh_front.m_fedges.clear();
  mesh_front.m_svertices.clear();
  auto VBO = mesh_front.face_property<FacingType>("f:VBO");
  for (size_t fi = 0; fi < mesh_front.faces_size(); fi++) {
    Face f(fi);
    if (VBO[f] != FacingType::FRONT)
      mesh_front.delete_face(f);
  }
  mesh_front.garbage_collection();
  mesh_front.write(name_pattern + "_interpolated_mesh_f.obj");
  logger().info("Wrote front-facing parts to obj: {}",
                name_pattern + "_interpolated_mesh_f.obj");

  Mesh mesh_back = mesh;
  mesh_back.m_fedges.clear();
  mesh_back.m_svertices.clear();
  VBO = mesh_back.face_property<FacingType>("f:VBO");
  for (size_t fi = 0; fi < mesh_back.faces_size(); fi++) {
    Face f(fi);
    if (VBO[f] != FacingType::BACK)
      mesh_back.delete_face(f);
  }
  mesh_back.garbage_collection();
  mesh_back.write(name_pattern + "_interpolated_mesh_b.obj");
  logger().info("Wrote front-facing parts to obj: {}",
                name_pattern + "_interpolated_mesh_b.obj");
}

void serialize_mesh(Mesh const &mesh, Camera const &camera,
                    std::vector<serialization::SerializationRay> const &rays,
                    std::string filename) {
  serialization::SerializationObject serial{&mesh, &camera, rays};

  logger().info("Serializating model to {}", filename);
  std::ofstream fstr(filename);
  serialization::save(fstr, serial, true, "json");
}

bool insert_contours(Mesh &mesh, Camera &camera, ContourMode contour_mode,
                     bool to_insert_cusps) {
  if (contour_mode == ContourMode::VBO_CONTOUR)
    mesh.computeConsistencySubdivided(camera.position());
  else
    mesh.computeConsistency(camera.position());

  if (contour_mode == ContourMode::VBO_CONTOUR) {
    int splits = 1, split_itr = 0;
    int max_split_itr = 5;
    do {
      logger().info("splitEdgesWithZeroCrossing round: {}", split_itr);
      splits = mesh.splitEdgesWithZeroCrossing(camera.position(),
                                               &mesh.subdivision());
      split_itr++;
    } while (splits > 0 && split_itr < max_split_itr);
    if (splits > 0 && split_itr == max_split_itr) {
      logger().warn("Warning: Mesh may miss zero-crossing vertices.");
    }
    mesh.computeConsistencySubdivided(camera.position());
    mesh.insertContours(camera.position(), &mesh.subdivision(), true);
  }
  mesh.extractBoundaries();
  // Throw unhandled exception when there's no intersection???
  // mesh.extractSurfaceIntersections();
  mesh.extractContours(contour_mode, camera);

  if (contour_mode == ContourMode::VBO_CONTOUR) {
    if (to_insert_cusps)
      insert_subdivision_cusps(mesh, mesh.subdivision(), camera);

    // Note that we don't update with computeConsistency even now we have
    // new topologies. Since we inserted new vertices to the limit
    // surfaces which may change the facing types of existing vertices.
    bool successful =
        consistently_label_interpolated_contours(mesh, camera.position());
    if (!successful)
      return false;
  }

  // The patch detection is required to relate the patch and the computed
  // rotation index in the output.
  mesh.floodFill();
  mesh.extractPatchBoundaries(false);

  return true;
}

bool verify_inconsistency(Mesh const &mesh, Camera const &camera,
                          FacingType patch_facing) {
  Vertex vs[3];
  auto positions = mesh.get_vertex_property<Vector3f>("v:point");
  size_t num_inconsist = 0;
  for (size_t i = 0; i < mesh.n_faces(); i++) {
    Face f = Face(i);
    // The facing type based on the face normal
    auto f_norm = mesh.compute_face_normal(f);

    mesh.verticesOfFace(f, vs);

    // Ignore extremely small triangles
    if ((positions[vs[0]] - positions[vs[1]]).norm() < 1e-5 &&
        (positions[vs[2]] - positions[vs[1]]).norm() < 1e-5 &&
        (positions[vs[0]] - positions[vs[2]]).norm() < 1e-5) {
      logger().warn("Ignore small triangle: {}", f);
      continue;
    }

    auto c_ray = positions[vs[0]] - camera.position();
    c_ray.normalized();
    FacingType f_ori =
        (f_norm.dot(c_ray) < 0) ? FacingType::FRONT : FacingType::BACK;
    if (f_ori != patch_facing) {
      num_inconsist++;
      logger().warn("Inconsistent triangle: {}", f);
    }
  }

  logger().info("\t=> #Inconsistent triangles: {} / {}", num_inconsist,
                mesh.n_faces());
  if (num_inconsist != 0) {
    return false;
  }
  return true;
}

void verify_qi_changes(Mesh const &mesh) {
  // Condition: change of QI must happen at a cusp or a 2D intersection
  auto is_cusp = mesh.get_vertex_property<bool>("v:cusp");
  auto intersection_2d = mesh.get_vertex_property<int>("v:intersection_2d");

  auto edge_qi = mesh.get_edge_property<int>("e:qi");
  contess_assert_msg(edge_qi && is_cusp && intersection_2d,
                     "verify_qi_changes: Missing contour information.");
  contess_assert_msg(!mesh.get_const_patch_chains().empty(),
                     "verify_qi_changes: Contour not chained.");

  std::unordered_set<int> visited_patch;
  size_t num_inconsist = 0;
  size_t num_qi_changes = 0;
  size_t num_inconsist_vis = 0;
  size_t num_qi_changes_vis = 0;
  std::unordered_set<int> seen_vertices;
  for (auto const &patch : mesh.get_const_patch_chains()) {
    if (visited_patch.count(patch.first))
      continue;

    visited_patch.insert(patch.first);

    // Set boundary cut directions
    auto chains = mesh.get_const_patch_chains().equal_range(patch.first);
    for (auto chain = chains.first; chain != chains.second; ++chain) {
      std::vector<Halfedge> chain_hes = *chain->second;
      for (size_t i = 0; i < chain_hes.size(); i++) {
        size_t j = (i + 1) % chain_hes.size();
        Halfedge h1 = chain_hes[i];
        Halfedge h2 = chain_hes[j];
        Vertex v = mesh.to_vertex(h1);

        if (seen_vertices.count(v.idx()))
          continue;
        seen_vertices.emplace(v.idx());

        // Determine the dead-end by checking QI
        if (edge_qi[mesh.edge(h1)] != edge_qi[mesh.edge(h2)]) {
          bool qi_change_condition = is_cusp[v] || intersection_2d[v] >= 0;
          num_qi_changes++;

          // If this change is visible
          if (std::min(edge_qi[mesh.edge(h1)], edge_qi[mesh.edge(h2)]) == 0)
            num_qi_changes_vis++;

          if (!qi_change_condition) {
            num_inconsist++;
            if (std::min(edge_qi[mesh.edge(h1)], edge_qi[mesh.edge(h2)]) == 0) {
              logger().error("Incorrect visible QI change at {}", v);
              num_inconsist_vis++;
            }
          }
        }
      }
    }
  }

  logger().info("\t=> #Inconsistent QI changes: {} / {}", num_inconsist,
                num_qi_changes);
  logger().info("\t=> #Inconsistent visible QI changes: {} / {}",
                num_inconsist_vis, num_qi_changes_vis);
}

void write_contours(Mesh const &mesh, Camera &camera, std::string filename,
                    int leq_qi) {
  Vector2i viewport = Vector2i(camera.vpWidth(), camera.vpHeight());
  SVG svgWriter(filename, viewport);
  auto is_contour = mesh.get_edge_property<real_t>("e:contour");
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto edge_qi = mesh.get_edge_property<int>("e:qi");

  for (auto eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit) {
    if (is_contour[*eit] < 0 && !mesh.is_boundary(*eit))
      continue;

    // QI threshold
    if (edge_qi[*eit] > leq_qi)
      continue;

    Vertex v = mesh.vertex(*eit, 0);
    Vertex v_to = mesh.vertex(*eit, 1);

    Vector3f v3d = vpositions[v];
    Vector2f v2d = project(v3d, camera.viewMatrix().matrix(),
                           camera.projectionMatrix(), viewport)
                       .head<2>();
    Vector3f v3d_to = vpositions[v_to];
    Vector2f v2d_to = project(v3d_to, camera.viewMatrix().matrix(),
                              camera.projectionMatrix(), viewport)
                          .head<2>();

    std::vector<Vector2f> poly;
    poly.emplace_back(v2d);
    poly.emplace_back(v2d_to);
    std::vector<SVG::PathVertexProperty> vertex_indices(
        {SVG::PathVertexProperty{v.idx()},
         SVG::PathVertexProperty{v_to.idx()}});

    svgWriter.writePolyline(poly, 1, Color(0, 0, 0, 255), false,
                            vertex_indices);
  }

  logger().info("Wrote to svg: {}", filename);
}

void write_chained_contours(Mesh const &mesh, Camera &camera,
                            std::string filename, int leq_qi) {
  Vector2i viewport = Vector2i(camera.vpWidth(), camera.vpHeight());
  SVG svgWriter(filename, viewport);
  auto vpositions = mesh.get_vertex_property<Vector3f>("v:point");
  auto edge_qi = mesh.get_edge_property<int>("e:qi");
  auto contour_type =
      mesh.get_vertex_property<ContourVertexType>("v:contour_type");
  auto intersection_2d = mesh.get_vertex_property<int>("v:intersection_2d");
  contess_assert_msg(edge_qi && contour_type && intersection_2d,
                     "write_chained_contours: Missing contour information.");
  contess_assert_msg(!mesh.get_const_patch_chains().empty(),
                     "write_chained_contours: Contour not chained.");

  auto get_type_str = [](ContourVertexType type,
                         int intersecting_index) -> std::string {
    switch (type) {
    case ContourVertexType::CONNECTOR:
      return "c";
    case ContourVertexType::DEAD_END:
      return "e";
    case ContourVertexType::JUNCTION:
    default:
      return "j_" + std::to_string(intersecting_index);
    }
    return "c";
  };

  std::unordered_set<int> visited_patch;
  for (auto const &patch : mesh.get_const_patch_chains()) {
    // Avoid duplicate contour output
    auto patch_facing = mesh.get_patch_facing(patch.first);

    if (visited_patch.count(patch.first))
      continue;

    visited_patch.insert(patch.first);

    // Set boundary cut directions
    auto chains = mesh.get_const_patch_chains().equal_range(patch.first);
    for (auto chain = chains.first; chain != chains.second; ++chain) {
      std::vector<Halfedge> chain_hes = *chain->second;
      bool first_visible = true;
      std::vector<Vector2f> poly;
      std::vector<SVG::PathVertexProperty> vertex_properties;

      // Find the start location (either an arbitrary vertex if all edges
      // are visible along the contour loop or the first edge adjacent to an
      // invisible one)
      size_t start_i = 0;
      for (size_t i = 0; i < chain_hes.size(); i++) {
        size_t j = (i + 1) % chain_hes.size();
        Halfedge h1 = chain_hes[i];
        Halfedge h2 = chain_hes[j];

        if (edge_qi[mesh.edge(h1)] != edge_qi[mesh.edge(h2)] &&
            edge_qi[mesh.edge(h1)] > leq_qi &&
            edge_qi[mesh.edge(h2)] <= leq_qi) {
          start_i = j;
        }
      }

      size_t i = start_i;
      do {
        Halfedge h1 = chain_hes[i];

        // QI threshold
        if (edge_qi[mesh.edge(h1)] > leq_qi ||
            (patch_facing != FacingType::FRONT &&
             !mesh.is_boundary(mesh.edge(h1)))) {
          first_visible = true;
          if (!poly.empty())
            svgWriter.writePolyline(poly, 1, Color(0, 0, 0, 255), false,
                                    vertex_properties);
          poly.clear();
          vertex_properties.clear();
        } else {
          Vertex v = mesh.from_vertex(h1);
          Vertex v_to = mesh.to_vertex(h1);

          Vector3f v3d = vpositions[v];
          Vector2f v2d = project(v3d, camera.viewMatrix().matrix(),
                                 camera.projectionMatrix(), viewport)
                             .head<2>();
          Vector3f v3d_to = vpositions[v_to];
          Vector2f v2d_to = project(v3d_to, camera.viewMatrix().matrix(),
                                    camera.projectionMatrix(), viewport)
                                .head<2>();

          if (first_visible) {
            poly.emplace_back(v2d);
            vertex_properties.emplace_back(SVG::PathVertexProperty{
                v.idx(), get_type_str(contour_type[v], intersection_2d[v]),
                (vpositions[v] - camera.position()).norm()});
          }
          poly.emplace_back(v2d_to);
          vertex_properties.emplace_back(SVG::PathVertexProperty{
              v_to.idx(),
              get_type_str(contour_type[v_to], intersection_2d[v_to]),
              (vpositions[v_to] - camera.position()).norm()});

          first_visible = false;
        }
        i = (i + 1) % chain_hes.size();
      } while (i != start_i);

      if (!poly.empty())
        svgWriter.writePolyline(poly, 1, Color(0, 0, 0, 255), false,
                                vertex_properties);
    }
  }
  logger().info("Wrote to svg: {}", filename);
}

Camera perturb_camera(int seed, real_t offset_scale,
                      Camera const &input_camera) {
  std::mt19937 eng(seed);
  std::uniform_real_distribution<real_t> urd(0, 1);

  // Perturb the three coordicates of the camera position
  Vector3f offset(offset_scale * urd(eng), offset_scale * urd(eng),
                  offset_scale * urd(eng));

  logger().info("Perturb camera with seed {}, offset: {}", seed,
                offset.transpose());

  Camera camera;
  Vector2i viewport({2048, 1600});
  camera.setViewport(viewport.x(), viewport.y());
  camera.lookAt(input_camera.position() + offset, input_camera.target(),
                input_camera.up());

  camera.projectionMatrix();
  camera.viewMatrix();

  return camera;
}

void perturb_mesh(Mesh &mesh, int seed, real_t offset_scale) {
  std::mt19937 eng(seed);
  std::uniform_real_distribution<real_t> urd(0, 1);
  auto vpositions = mesh.vertex_property<Vector3f>("v:point");

  logger().info("Perturb mesh with seed {}", seed);
  for (size_t i = 0; i < mesh.n_vertices(); i++) {
    Vertex v(i);

    Vector3f offset(offset_scale * urd(eng), offset_scale * urd(eng),
                    offset_scale * urd(eng));
    vpositions[v] += offset;
  }
}
