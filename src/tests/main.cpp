// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "logger.h"
#include "spdlog/common.h"
#include "testing.h"
#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/sinks/stdout_sinks.h>

#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>

// for file sink of logging
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <igl/predicates/predicates.h>

int main(int argc, char *argv[]) {
  // Set up logger here
  logger().set_level(spdlog::level::info);
  logger().sinks().push_back(
      std::make_shared<spdlog::sinks::stdout_color_sink_st>());
  logger().sinks().back()->set_level(spdlog::level::info);

  igl::predicates::exactinit();

  // Parse argv and argc
  std::filesystem::path spec = TestingSingleton::instance().repo_root_path() /
                               "test-specs" / "spec.json";
  for (int i = argc - 1; i >= 0; --i) {
    auto arg = std::string(argv[i]);
    if (arg.find("spec=") == 0) {
      auto len = std::string("spec=").length();
      spec = arg.substr(len);
      argc--;
      for (int j = i; j < argc; ++j) {
        argv[j] = argv[j + 1];
      }
    }
    if (arg.find("out=") == 0) {
      auto len = std::string("out=").length();
      output_folder = arg.substr(len);
      argc--;
      for (int j = i; j < argc; ++j) {
        argv[j] = argv[j + 1];
      }
    }
    if (arg == "-o") {
      logger().sinks().push_back(
          std::make_shared<spdlog::sinks::ostream_sink_st>(std::cout));
      logger().sinks().back()->set_level(spdlog::level::debug);
    }
  }

  // Load the specs
  logger().info("Spec file is {}", spec.generic_string());
  TestingSingleton::instance().load_cases(spec);

  // This runs the unit tests.
  return Catch::Session().run(argc, argv);
}
