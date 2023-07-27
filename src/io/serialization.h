// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#pragma once

#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <tuple>
#include <variant>

#include <nlohmann/json.hpp>

#include <Eigen/Core>
#include <type_traits>
#include <vector>

#include "common.h"
#include "logger.h"

class Mesh;
class Camera;
struct TestCase;
class BasicCameraCase;
class FileCameraCase;
struct TestingSingleton;

namespace surface_mesh {
class Surface_mesh;
}

#define contess_befriend_serialization(class_)                                 \
  friend void serialization::serialize(serialization::json_t &,                \
                                       const class_ &);                        \
  friend void serialization::deserialize(class_ &,                             \
                                         const serialization::json_t &)

namespace serialization {

// A struct to pack the ray for visualization
struct SerializationRay {
  Vector3f origin, dir;
  Vector3f query_point, query_point_epsilon;
  int query_edge;
  std::vector<Vector3f, Eigen::aligned_allocator<Vector3f>> intersection_points;
  std::vector<int> intersection_faces;
};

// A struct to pack all the objects we want to serialize.
struct SerializationObject {
  const Mesh *mesh;
  const Camera *camera;
  std::vector<SerializationRay> rays;
};

// Helper for static asserts
template <typename... TypesL> constexpr bool always_false_v = false;

using json_t = nlohmann::json;

// ==================[Forward declare] ===================
template <typename Derived>
void serialize(json_t &j, const Eigen::DenseBase<Derived> &eig);
template <typename T, typename Allocator>
void serialize(json_t &j, const std::vector<T, Allocator> &vec);
template <typename T, typename GetterFunc>
void serialize_array(json_t &v, const GetterFunc &getter,
                     const std::size_t num_values);

// Variant
template <typename... Args>
void serialize(json_t &j, const std::variant<Args...> &v);
template <typename... Args>
void deserialize(std::variant<Args...> &v, const json_t &j);

// ==================[ Atomic types] ===================
void serialize(json_t &, const double &);
void deserialize(double &, const json_t &);
//
void serialize(json_t &, const float &);
void deserialize(float &, const json_t &);
//
void serialize(json_t &, const int &);
void deserialize(int &, const json_t &);
//
void serialize(json_t &j, const unsigned int &v);
void deserialize(unsigned int &v, const json_t &j);
//
void serialize(json_t &, const std::string &);
void deserialize(std::string &, const json_t &);

void serialize(json_t &, const bool &);
void deserialize(bool &, const json_t &);

void serialize(json_t &, const std::filesystem::path &);
void deserialize(std::filesystem::path &, const json_t &);

template <typename T>
std::enable_if_t<std::is_enum_v<T>> serialize(json_t &v, const T &enum_) {
  const int vint = static_cast<int>(enum_);
  v = vint;
}
template <typename T>
std::enable_if_t<std::is_enum_v<T>> deserialize(T &enum_, const json_t &v) {
  int vint = v;
  enum_ = static_cast<T>(vint);
}

// ==================[ Surface_mesh ] ===================
/**
Currently does not save everything.
Will add more as need be.
*/
void serialize(json_t &, const surface_mesh::Surface_mesh &);
inline void deserialize(surface_mesh::Surface_mesh &,
                        const json_t &) { /* not implemented */
}

// ==================[ Mesh ] ===================
/**
TODO: if we need to serialize things that are in Mesh but not surface_mesh, we
would need a separate serialization funtion.
*/
// void serialize(json_t&, const Mesh&);
// inline void deserialize(Mesh&, const json_t&) { /* not implemented */ }

// ==================[ Camera ] ===================
/**
 Note: does not support serializing animation or keyframes (can be added).
 Just saves the model view matrix, view port and the projection matrix.
*/
void serialize(json_t &, const Camera &);
inline void deserialize(Camera &, const json_t &) { /* not implemented */
}

// ==================[  SerializationObject  ] ===================
void serialize(json_t &, const SerializationObject &);
inline void deserialize(SerializationObject &,
                        const json_t &) { /* not implemented */
}

// ==================[  SerializationRay  ] ===================
void serialize(json_t &, const SerializationRay &);
inline void deserialize(SerializationRay &,
                        const json_t &) { /* not implemented */
}

// ==================[ TestCase ] ===================

void serialize(json_t &, const TestCase &);
void deserialize(TestCase &, const json_t &);

void serialize(json_t &, const BasicCameraCase &);
void deserialize(BasicCameraCase &, const json_t &);

void serialize(json_t &, const FileCameraCase &);
void deserialize(FileCameraCase &, const json_t &);

void serialize(json_t &, const TestingSingleton &);
void deserialize(TestingSingleton &, const json_t &);

// ==================[ Raw Array ] ===================
template <typename T, typename GetterFunc>
void serialize_array(json_t &v, const GetterFunc &getter,
                     const std::size_t num_values) {
  v = json_t::array();
  for (std::size_t i = 0; i < num_values; ++i) {
    json_t subv;
    T ith_value = getter(i);
    serialize(subv, ith_value);
    v.push_back(subv);
  }
}
template <typename T, typename ResizeFunc, typename SetValueFunc>
void deserialize_array(const json_t &v, const ResizeFunc &resize,
                       const SetValueFunc &set_value) {
  contess_assert(v.is_array());
  resize(v.size());
  for (std::size_t i = 0; i < v.size(); ++i) {
    T v_;
    deserialize(v_, v[i]);
    set_value(i, v_);
  }
}

// ==================[ Std vector ] ===================
template <typename T, typename Allocator>
void serialize(json_t &j, const std::vector<T, Allocator> &vec) {
  serialize_array<T>(
      j, [&vec](const std::size_t i) { return vec[i]; }, vec.size());
}
template <typename T, typename Allocator>
void deserialize(std::vector<T, Allocator> &vec, const json_t &v) {
  deserialize_array<T>(
      v, [&vec](const std::size_t sz) { vec.resize(sz); },
      [&vec](const std::size_t i, const T &v_) { vec[i] = v_; });
}

// ==================[ Eigen Types ] ===================
template <typename Derived>
void serialize(json_t &j, const Eigen::DenseBase<Derived> &eig) {
  using Scalar = typename Derived::Scalar;
  static constexpr bool is_rows_fixed =
      Derived::RowsAtCompileTime != Eigen::Dynamic;
  static constexpr bool is_cols_fixed =
      Derived::ColsAtCompileTime != Eigen::Dynamic;

  if constexpr (is_cols_fixed || is_rows_fixed) {
    // Always write as col major
    serialize_array<Scalar>(
        j,
        [&eig](const std::size_t i) {
          return eig(i % eig.rows(), i / eig.rows());
        },
        static_cast<std::size_t>(eig.size()));
  } else {
    // If you use false instead, the compiler will go crazy.
    static_assert(always_false_v<Derived>, "Not supported at the moment");
  }
}

template <typename Derived>
void deserialize(Eigen::DenseBase<Derived> &eig, const json_t &v) { /* todo */
}

// ==================[ Variant ] ===================
template <typename... Args>
void serialize(json_t &j, const std::variant<Args...> &v) {
  j = json_t::object();
  serialize(j["index"], (unsigned int)(v.index()));
  std::visit(
      contess::overloaded{[&](const std::monostate &v) {},
                          [&](const auto &v) { serialize(j["value"], v); }},
      v);
}

template <int N, typename... Ts>
using NthTypeOf = std::tuple_element_t<N, std::tuple<Ts...>>;

template <typename... Args>
void deserialize(std::variant<Args...> &v, const json_t &j) {
  static_assert(sizeof...(Args) <= 4, "Only 4 args supported at the moment");
  unsigned int index;
  deserialize(index, j["index"]);

#define IF_CONSTRUCTIBLE(N)                                                    \
  if constexpr (N < sizeof...(Args))                                           \
    if constexpr (std::is_default_constructible_v<NthTypeOf<N, Args...>>)      \
      if (index == N)

  IF_CONSTRUCTIBLE(0) { v = NthTypeOf<0, Args...>(); }
  IF_CONSTRUCTIBLE(1) { v = NthTypeOf<1, Args...>(); }
  IF_CONSTRUCTIBLE(2) { v = NthTypeOf<2, Args...>(); }
  IF_CONSTRUCTIBLE(3) { v = NthTypeOf<3, Args...>(); }
  IF_CONSTRUCTIBLE(4) { v = NthTypeOf<4, Args...>(); }

  std::visit(contess::overloaded{[&](std::monostate &v) {},
                                 [&](auto &v) { deserialize(v, j["value"]); }},
             v);

#undef IF_CONSTRUCTIBLE
}

// ==================[ Serialize and Save ] ===================
/**
 Note: this only works for types that have serialize() and deserize() defined
*/
template <typename T>
void save(std::ofstream &fstr, const T &obj, bool fancy = false,
          const std::string &format = "json") {
  contess_assert(fstr.is_open());
  json_t j;

  serialize(j, obj);
  {
    time_t rawtime;
    struct tm *timeinfo;
    char buffer[80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer, sizeof(buffer), "%d-%m-%Y %H:%M:%S", timeinfo);
    j["saved_at"] = buffer;
  }

  if (format == "json") {
    logger().debug("Writing json file");
    if (!fancy) {
      fstr << j.dump();
    } else {
      fstr << j.dump(1, '\t');
    }
  } else if (format == "cbor") {
    logger().debug("Writing cbor file");
    std::vector<std::uint8_t> v_cbor = json_t::to_cbor(j);
    fstr.write(reinterpret_cast<const char *>(v_cbor.data()),
               v_cbor.size() * sizeof(std::uint8_t));
  } else {
    contess_assert_msg(0, "format `" << format
                                     << "` not supported for serialization");
  }
}
template <typename T>
void load(T &v, const std::filesystem::path &fname) {
  std::ifstream fstr(fname);
  contess_assert_msg(fstr.is_open(), "Could not open " << fname.string());
  std::string format = fname.extension().string();
  json_t j;
  if (format == ".json") {
    fstr >> j;
  } else if (format == ".cbor") {
    contess_assert_msg(0, "I have not written reading cbor yet.");
  } else {
    contess_assert_msg(0, "format `" << format
                                     << "` not supported for serialization");
  }
  deserialize(v, j);
}

} // namespace serialization