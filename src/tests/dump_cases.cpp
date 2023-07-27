#include <catch2/catch.hpp>
#include <foldertools.h>
#include <io/serialization.h>
#include <testing.h>
#include <testing_cases.h>

#include <string>

TEST_CASE("dump_cases", "[.]") {
  TestingSingleton singleton;

  auto dump_cases = [&](const std::vector<TestCase> &cases, std::string name) {
    auto spec_path = std::filesystem::path("test-specs");
    const std::filesystem::path outdir = spec_path / name;
    foldertools::makedir(outdir.generic_string().c_str());
    for (std::size_t i = 0; i < cases.size(); ++i) {
      const std::filesystem::path outfile =
          outdir / (cases[i].model_name + ".json");
      singleton.multicamera_test_cases.push_back(
          std::filesystem::relative(outfile, spec_path));
      singleton.pipeline_test_cases.push_back(
          std::filesystem::relative(outfile, spec_path));
      std::ofstream ofs(outfile);
      serialization::save(ofs, cases[i], true);
    }
  };

  dump_cases(bunny_test_cases_p1, "bunny-p1");
  dump_cases(bunny_test_cases_p2, "bunny-p2");
  dump_cases(Angela_test_cases, "Angela");
  dump_cases(walking_test_cases, "walking");
  dump_cases(seq_subdiv_test_cases, "seq_subdiv");
  dump_cases(perturb_test_cases, "perturb");
  dump_cases(test_pipeline_cases, "pipeline");
  std::ofstream ofs(std::filesystem::path("test-specs") / "spec.json");
  serialization::save(ofs, singleton, true);
}