// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

#include "logger.h"

// Retrieve current logger
spdlog::logger &logger() {
  static auto default_logger = std::make_shared<spdlog::logger>("contess");
  return *default_logger;
}
