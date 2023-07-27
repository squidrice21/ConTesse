// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#ifndef _COMMON_H
#define _COMMON_H

#include <algorithm>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <tbb/blocked_range.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/spin_mutex.h>
#include <tbb/task.h>
#include <tbb/task_scheduler_init.h>

#include <nanogui/common.h>
#include <nanogui/glutil.h>

// Some header is polluting macros by including windows.h
// I don't know which one it is, but these should take care of it.
#undef near
#undef far
#undef min
#undef max

#define GRAIN_SIZE 1024

#define contess_str(msg)                                                       \
  static_cast<std::ostringstream &>(std::ostringstream().flush() << msg).str()
#define contess_msg(msg)                                                       \
  contess_str("File: " << __FILE__ << " Line: " << __LINE__                    \
                       << "\n Message: " << msg)
#define contess_assert_msg(cond, msg)                                          \
  do {                                                                         \
    if (!(cond))                                                               \
      throw std::runtime_error(contess_msg(msg));                              \
  } while (0)
#define contess_assert(cond) contess_assert_msg(cond, "")

// #define SINGLE_PRECISION
/* Application precision -- can be set to single or double precision */
#if defined(SINGLE_PRECISION)
typedef float real_t;
#else
typedef double real_t;
#endif

/* Useful Eigen typedefs based on the current precision */
typedef Eigen::Matrix<int32_t, 2, 1> Vector2i;
typedef Eigen::Matrix<int32_t, 3, 1> Vector3i;
typedef Eigen::Matrix<int32_t, 4, 1> Vector4i;
typedef Eigen::Matrix<uint32_t, 2, 1> Vector2u;
typedef Eigen::Matrix<uint32_t, 3, 1> Vector3u;
typedef Eigen::Matrix<uint32_t, 4, 1> Vector4u;
typedef Eigen::Matrix<uint8_t, 4, 1> Vector4u8;
typedef Eigen::Matrix<real_t, 2, 1> Vector2f;
typedef Eigen::Matrix<real_t, 3, 1> Vector3f;
typedef Eigen::Matrix<real_t, 4, 1> Vector4f;
typedef Eigen::Matrix<int32_t, Eigen::Dynamic, 1> VectorXi;
typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, 1> VectorXu;
typedef Eigen::Matrix<uint8_t, Eigen::Dynamic, 1> VectorXu8;
typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> VectorXb;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, 1> VectorXf;
typedef Eigen::Matrix<int32_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXi;
typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXu;
typedef Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXu8;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXf;
typedef Eigen::Matrix<real_t, 3, Eigen::Dynamic> Matrix3Xf;
typedef Eigen::Matrix<real_t, 2, Eigen::Dynamic> Matrix2Xf;
typedef Eigen::Matrix<real_t, 2, 2> Matrix2f;
typedef Eigen::Matrix<real_t, 3, 3> Matrix3f;
typedef Eigen::Matrix<real_t, 4, 4> Matrix4f;
typedef Eigen::Quaternion<real_t> Quaternionf;
typedef Eigen::Transform<real_t, 3, Eigen::Affine> Affine3f;
typedef Eigen::Translation<real_t, 3> Translation3f;
typedef Eigen::AngleAxis<real_t> AngleAxisf;

typedef Eigen::Map<MatrixXu> MapXu;
typedef Eigen::Map<MatrixXf> MapXf;
typedef Eigen::Map<MatrixXu const> MapConstXu;
typedef Eigen::Map<MatrixXf const> MapConstXf;

using std::cerr;
using std::cout;
using std::endl;

// Helpers for std::variant
// https://dzone.com/articles/two-lines-of-code-and-three-c17-features-the-overl
namespace contess {
template <class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;
} // namespace contess

// TODO: Distinguish the use of UNDEFINED and CONTOUR
// It seems that the current code uses UNDEFINED to indicate interpolated
// contour faces.
typedef enum { FRONT, BACK, CONTOUR, UNDEFINED, NA } FacingType;
typedef enum { NEAR, FAR } NearFarType;
typedef enum { MESH_CONTOUR, INTERPOLATED_CONTOUR, VBO_CONTOUR } ContourMode;
typedef enum { MESH_CONVEXITY, CAMERA_CONVEXITY } ConvexityMode;
typedef enum { PATCH, CUT, COMPONENT, NONE } BoundaryType;
typedef enum { CONNECTOR, DEAD_END, JUNCTION } ContourVertexType;

typedef unsigned short nature_t;
namespace VertexNature {
static const nature_t S_VERTEX = 0;
static const nature_t INTERSECTION_3D = (1 << 1);
static const nature_t INTERSECTION_2D = (1 << 2);
static const nature_t BOUNDARY_CURTAIN_FOLD = (1 << 3);
static const nature_t CONTOUR_CURTAIN_FOLD = (1 << 4);
static const nature_t BIFURCATION = (1 << 5);
static const nature_t CURTAIN_FOLD =
    BOUNDARY_CURTAIN_FOLD | CONTOUR_CURTAIN_FOLD;
static const nature_t RAY_INTERSECTION = (1 << 6);
static const nature_t CONTOUR_FRONT = (1 << 7);
} // namespace VertexNature
namespace EdgeNature {
static const nature_t NO_FEATURE = 0;
static const nature_t SHARP_CONTOUR = (1 << 0);
static const nature_t SMOOTH_CONTOUR = (1 << 1);
static const nature_t CONTOUR = SHARP_CONTOUR | SMOOTH_CONTOUR;
static const nature_t BOUNDARY = (1 << 2);
static const nature_t SURFACE_INTERSECTION = (1 << 3);
static const nature_t IMAGE_INTERSECTION = (1 << 4);
} // namespace EdgeNature

const float CONTOUR_THRESHOLD = 1e-6;
// const float CONTOUR_THRESHOLD = 1e-8;
// const int MAX_ROOT_ITERATIONS = 100;
const int MAX_ROOT_ITERATIONS = 1000;
const float MAX_SHIFT_PERCENTAGE = 0.2;
const float EPSILON = 1e-8;
const float CUSP_ROOT_THRESHOLD = 1e-4;
const double EXTRAORDINARY_REGION_OFFSET = 1e-5;
const float MIN_EDGE_LENGTH = 1e-6;

// Epsilon constant used by the convexity tagging to ensure
// the front-portion vertex is sufficient far from the contour
const float FRONT_PORTION_VERTEX_EPSILON = 1e-3;

using Color = nanogui::Color;

const Color frontFaceColor = Color(254, 242, 192, 255);
const Color backFaceColor = Color(139, 191, 230, 255);
const Color undefinedFaceColor = Color(250, 50, 0, 255);
const Color inconsistentFaceColor = Color(0.8f, 0.1f, 0.8f, 1.f);
const Color wireframeColor = Color(0.3f, 0.3f, 0.3f, 1.f);
const Color normalsColor = Color(0.7f, 0.7f, 0.f, 1.f);
const Color controlPolygonColor = Color(0.f, 0.8f, 0.5f, 1.f);
const Color contourColor = Color(0.6f, 0.0f, 0.f, 1.f);
const Color hiddenContourColor = Color(0.8f, 0.8f, 0.8f, 1.f);
const Color convexContourColor = Color(252, 186, 3, 255);
const Color concaveContourColor = Color(5, 245, 245, 255);
const Color boundariesColor = Color(0.1f, 0.0f, 0.8f, 1.f);
const Color surfIntersectionsColor = Color(0.0f, 0.5f, 0.1f, 1.f);
const Color seamsColor = Color(1.0f, 0.0f, 0.0f, 1.f);
const Color intersection2DColor = Color(0.0f, 1.0f, 0.0f, 1.f);
const Color intersection3DColor = Color(0.0f, 1.0f, 1.0f, 1.f);
const Color curtainFoldColor = Color(1.0f, 0.6f, 0.0f, 1.f);
const Color bifurcationColor = Color(1.0f, 0.0f, 0.0f, 1.f);
const Color frontVertexColor = Color(227, 91, 0, 255);

// Paul Green-Armytage "A Colour Alphabet and the Limits of Colour Coding."
// (2010)
const Color alphabetColors[26] = {
    Color(240, 163, 255, 255), Color(0, 117, 220, 255),
    Color(153, 63, 0, 255),    Color(76, 0, 92, 255),
    Color(43, 206, 72, 255),   Color(255, 204, 153, 255),
    Color(128, 128, 128, 255), Color(148, 255, 181, 255),
    Color(143, 124, 0, 255),   Color(157, 204, 0, 255),
    Color(194, 0, 136, 255),   Color(0, 51, 128, 255),
    Color(25, 25, 25, 255),    Color(0, 92, 49, 255),
    Color(153, 0, 0, 255),     Color(255, 255, 128, 255),
    Color(255, 164, 5, 255),   Color(255, 168, 187, 255),
    Color(224, 255, 102, 255), Color(116, 10, 255, 255),
    Color(66, 102, 0, 255),    Color(255, 0, 16, 255),
    Color(94, 241, 242, 255),  Color(0, 153, 143, 255),
    Color(255, 255, 0, 255),   Color(255, 80, 5, 255)};

inline Vector3f project(const Vector3f &obj, const Matrix4f &model,
                        const Matrix4f &proj, const Vector2i &viewportSize) {
  Vector4f tmp;
  tmp << obj, 1;

  tmp = model * tmp;

  tmp = proj * tmp;

  tmp = tmp.array() / tmp(3);
  tmp = tmp.array() * 0.5f + 0.5f;
  tmp(0) = tmp(0) * viewportSize.x();
  tmp(1) = tmp(1) * viewportSize.y();

  return tmp.head(3);
}

inline Vector3f unproject(const Vector3f &win, const Matrix4f &model,
                          const Matrix4f &proj, const Vector2i &viewportSize) {
  Matrix4f Inverse = (proj * model).inverse();

  Vector4f tmp;
  tmp << win, 1;
  tmp(0) = tmp(0) / viewportSize.x();
  tmp(1) = tmp(1) / viewportSize.y();
  tmp = tmp.array() * 2.0f - 1.0f;

  Vector4f obj = Inverse * tmp;
  obj /= obj(3);

  return obj.head(3);
}

inline void projectToViewport(const Matrix3Xf &chain,
                              std::vector<Vector2f> &polyline,
                              const Matrix4f &projectionMatrix,
                              const Affine3f &viewMatrix,
                              const Vector2i &viewport) {
  polyline.reserve(polyline.size() + chain.cols());
  for (int j = 0; j < chain.cols(); ++j) {
    Vector3f pos2D =
        project(chain.col(j), viewMatrix.matrix(), projectionMatrix, viewport);
    polyline.push_back(pos2D.head<2>());
  }
}

inline void projectToViewportDedup(const Matrix3Xf &chain,
                                   std::vector<Vector2f> &polyline,
                                   const Matrix4f &projectionMatrix,
                                   const Affine3f &viewMatrix,
                                   const Vector2i &viewport) {
  projectToViewport(chain, polyline, projectionMatrix, viewMatrix, viewport);

  if (polyline.empty())
    return;

  // Clean the polygon vertices to remove duplicated ones
  std::vector<Vector2f> dedup_polyline;
  dedup_polyline.reserve(polyline.size());
  dedup_polyline.emplace_back(polyline[0]);
  for (size_t i = 1; i < polyline.size(); i++) {
    if ((dedup_polyline.back() - polyline[i]).norm() <
        std::numeric_limits<double>::epsilon())
      continue;
    dedup_polyline.emplace_back(polyline[i]);
  }
  // The polygon is closed by default. No need to duplicate.
  if (dedup_polyline.size() > 1 &&
      (dedup_polyline.front() - dedup_polyline.back()).norm() <
          std::numeric_limits<double>::epsilon())
    dedup_polyline.pop_back();
  polyline = dedup_polyline;
}

inline real_t vertex_to_edge_distance(Vector3f const &v1, Vector3f const &v2,
                                      Vector3f const &v) {
  Vector3f v12 = v2 - v1;
  Vector3f v1v = v - v1;

  // Unit vector orthogonal to the edge
  Vector3f up = (v12).cross(v1v);
  Vector3f unit_to_v = up.cross(v12);
  unit_to_v.normalize();

  // Measure distance
  real_t dist = std::abs(v1v.dot(unit_to_v));

  return dist;
}

template <typename T> void delete_pointed_to(T *const ptr) { delete ptr; }

#if defined(_WIN32)
#define RCPOVERFLOW_FLT 2.93873587705571876e-39f
#define RCPOVERFLOW_DBL 5.56268464626800345e-309
#else
#define RCPOVERFLOW_FLT 0x1p-128f
#define RCPOVERFLOW_DBL 0x1p-1024
#endif

#if defined(SINGLE_PRECISION)
#define RCPOVERFLOW RCPOVERFLOW_FLT
#else
#define RCPOVERFLOW RCPOVERFLOW_DBL
#endif

/* A callback to inform the GUI about progress of an operation */
typedef std::function<void(const std::string &, float)> ProgressCallback;

#define PROGRESS_BLKSIZE (1 << 18)
#define SHOW_PROGRESS(i, maxval, text)                                         \
  if (progress && (i % PROGRESS_BLKSIZE) == 0)                                 \
  progress(text, -PROGRESS_BLKSIZE / (float)maxval)

#define PROGRESS_SHIFT 18u
#define SHOW_PROGRESS_RANGE(range, maxval, text)                               \
  if (progress && range.begin() > 0) {                                         \
    uint32_t nUpdates = (range.end() >> PROGRESS_SHIFT) -                      \
                        ((range.begin() - 1) >> PROGRESS_SHIFT);               \
    if (nUpdates > 0) {                                                        \
      const uint32_t nUpdatesTotal =                                           \
          (uint32_t)(maxval) / (1 << PROGRESS_SHIFT);                          \
      progress(text, -(float)nUpdates / (float)nUpdatesTotal);                 \
    }                                                                          \
  }

template <typename TimeT = std::chrono::milliseconds> class Timer {
public:
  Timer() { start = std::chrono::system_clock::now(); }

  size_t value() const {
    auto now = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<TimeT>(now - start);
    return (size_t)duration.count();
  }

  size_t reset() {
    auto now = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<TimeT>(now - start);
    start = now;
    return (size_t)duration.count();
  }

private:
  std::chrono::system_clock::time_point start;
};

inline std::string timeString(double time, bool precise = false) {
  if (std::isnan(time) || std::isinf(time))
    return "inf";

  std::string suffix = "ms";
  if (time > 1000) {
    time /= 1000;
    suffix = "s";
    if (time > 60) {
      time /= 60;
      suffix = "m";
      if (time > 60) {
        time /= 60;
        suffix = "h";
        if (time > 12) {
          time /= 12;
          suffix = "d";
        }
      }
    }
  }

  std::ostringstream os;
  os << std::setprecision(precise ? 4 : 1) << std::fixed << time << suffix;

  return os.str();
}

inline std::string memString(size_t size, bool precise = false) {
  double value = (double)size;
  const char *suffixes[] = {"B", "KiB", "MiB", "GiB", "TiB", "PiB"};
  int suffix = 0;
  while (suffix < 5 && value > 1024.0f) {
    value /= 1024.0f;
    ++suffix;
  }

  std::ostringstream os;
  os << std::setprecision(suffix == 0 ? 0 : (precise ? 4 : 1)) << std::fixed
     << value << " " << suffixes[suffix];

  return os.str();
}

inline void coordinate_system(const Vector3f &a, Vector3f &b, Vector3f &c) {
  if (std::abs(a.x()) > std::abs(a.y())) {
    real_t invLen = 1.0f / std::sqrt(a.x() * a.x() + a.z() * a.z());
    c = Vector3f(a.z() * invLen, 0.0f, -a.x() * invLen);
  } else {
    real_t invLen = 1.0f / std::sqrt(a.y() * a.y() + a.z() * a.z());
    c = Vector3f(0.0f, a.z() * invLen, -a.y() * invLen);
  }
  b = c.cross(a);
}

inline std::string str_tolower(std::string str) {
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  return str;
}

inline uint32_t str_to_uint32_t(const std::string &str) {
  char *end_ptr = nullptr;
  uint32_t result = (uint32_t)strtoul(str.c_str(), &end_ptr, 10);
  if (*end_ptr != '\0')
    throw std::runtime_error("Could not parse unsigned integer \"" + str +
                             "\"");
  return result;
}

inline uint32_t str_to_int32_t(const std::string &str) {
  char *end_ptr = nullptr;
  int32_t result = (int32_t)strtol(str.c_str(), &end_ptr, 10);
  if (*end_ptr != '\0')
    throw std::runtime_error("Could not parse signed integer \"" + str + "\"");
  return result;
}

inline float str_to_float(const std::string &str) {
  char *end_ptr = nullptr;
  float result = (float)strtod(str.c_str(), &end_ptr);
  if (*end_ptr != '\0')
    throw std::runtime_error("Could not parse floating point value \"" + str +
                             "\"");
  return result;
}

#endif
