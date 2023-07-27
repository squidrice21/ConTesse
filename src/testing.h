#pragma once

#include <filesystem>
#include <string>
#include <variant>

#include "common.h"
#include "mesh.h"

extern std::string output_folder;

class BasicCameraCase {
public:
  real_t camera_x = 0;
  real_t camera_y = 0;
  real_t camera_z = 0;
  bool full_info = false;
  Vector3f lookat = Vector3f::Zero();
  BasicCameraCase() = default;
  BasicCameraCase(real_t, real_t);
  BasicCameraCase(real_t, real_t, real_t, Vector3f const &);
};

class FileCameraCase {
public:
  std::filesystem::path path = "";
  std::size_t frame_idx = 0;
  FileCameraCase() = default;
  FileCameraCase(std::filesystem::path, std::size_t);
};

struct TestCase {
  std::string model_name = "";
  std::string model_path = "";
  std::string description = "";
  std::variant<std::monostate, BasicCameraCase, FileCameraCase> camera_info;
  int subdiv_level = 2;

  TestCase() = default;
  std::string get_short_name() const;
  std::string get_description() const;
};

struct TestingSingleton {
  std::vector<std::filesystem::path> multicamera_test_cases;
  std::vector<std::filesystem::path> pipeline_test_cases;
  std::filesystem::path spec_root_path;

  void load_cases(std::filesystem::path);
  std::filesystem::path repo_root_path();

  static TestingSingleton &instance();
};

void setup_camera(Camera &, const Mesh &mesh, const TestCase &);

std::string test_output_dir();

// https://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 3) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

Mesh read_input(std::string model_filename, bool is_full_path = false,
                int subdiv_level = -1,
                Subdiv::Backend sub_backend = Subdiv::Backend::LACEWELL);
void write_param_to_vertex(Mesh &mesh);
void write_facing_labelled_mesh(Mesh const &mesh, std::string name_pattern);
void serialize_mesh(Mesh const &mesh, Camera const &camera,
                    std::vector<serialization::SerializationRay> const &rays,
                    std::string filename);

bool insert_contours(
    Mesh &mesh, Camera &camera,
    ContourMode contour_mode = ContourMode::INTERPOLATED_CONTOUR,
    bool to_insert_cusps = false);

bool verify_inconsistency(Mesh const &mesh, Camera const &camera,
                          FacingType patch_facing);
void verify_qi_changes(Mesh const &mesh);

void write_contours(Mesh const &mesh, Camera &camera, std::string filename,
                    int leq_qi = std::numeric_limits<int>::infinity());
void write_chained_contours(Mesh const &mesh, Camera &camera,
                            std::string filename,
                            int leq_qi = std::numeric_limits<int>::infinity());

Camera perturb_camera(int seed, real_t offset_scale,
                      Camera const &input_camera);

void perturb_mesh(Mesh &mesh, int seed, real_t offset_scale);
