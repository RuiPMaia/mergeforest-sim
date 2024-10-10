#include <mergeforest-sim/sparse_matrix.hpp>
#include <mergeforest-sim/matrix_IO.hpp>
#include <mergeforest-sim/gen_matrix.hpp>
#include <mergeforest-sim/simulator.hpp>

#include <CLI/CLI.hpp>
#include <spdlog/spdlog.h>

#include <string>
#include <filesystem>

using namespace mergeforest_sim;
namespace fs = std::filesystem;

int run_sim_app(CLI::App& app) {
  fs::path matrix_file1;
  fs::path matrix_file2;
  fs::path config_file;
  fs::path output_path;
  std::string out_filename;
  bool compute_result {true};

  app.add_option("-m,--matrix,--matrix1", matrix_file1, "matrix file")
    ->required()->check(CLI::ExistingFile);
  app.add_option("--matrix2", matrix_file2, "matrix file")
    ->check(CLI::ExistingFile);
  app.add_option("-c,--config", config_file, "config file")
    ->required()->check(CLI::ExistingFile);
  const auto sim_outdir_opt = app.add_option("-o,--outdir", output_path,
                                             "output directory");
  app.add_option("--outname", out_filename, "output filename")
    ->needs(sim_outdir_opt);
  app.add_flag("--compute-result,--no-compute-result{false}",
               compute_result, "compute result");

  try {
    app.parse(app.remaining_for_passthrough());
  } catch(const CLI::ParseError& e) { return app.exit(e); }

  if (!output_path.empty()) {
    fs::create_directories(output_path);
    if (out_filename.empty()) {
      out_filename = matrix_file1.stem();
      if (!matrix_file2.empty()) {
        out_filename += '_' + matrix_file2.stem().string();
      }
      out_filename += '_' + config_file.stem().string() + "_sim_results.txt";
    }
    output_path /= out_filename;
  }

  fmt::print("Loading matrix A: {}... ", matrix_file1.string());
  fflush(stdout);
  Spmat_Csr A(matrix_file1);
  fmt::print("Done\n");
  Spmat_Csr B;
  if (matrix_file2.empty()) {
    if (A.num_rows == A.num_cols) {
      fmt::print("Matrix B = A\n");
      B = A;
    } else {
      fmt::print("Matrix B = A^T\n");
      B = A.transpose();
    }
  } else {
    fmt::print("Loading matrix B: {}... ", matrix_file2.string());
    fflush(stdout);
    B = read_matrix_market_file(matrix_file2);
    fmt::print("Done\n");
  }
  Simulator simulator(config_file, output_path);
  simulator.set_mats(A, B);
  fmt::print("Starting simulation...\n");
  simulator.run_simulation(compute_result);
  if (!output_path.empty()) {
    fmt::print("Simulation results written to {}\n", output_path.c_str());
  }
  return 0;
}

int run_stats_app(CLI::App& app) {
  fs::path matrix_file1;
  fs::path matrix_file2;
  fs::path config_file;
  fs::path output_path;
  std::string out_filename;

  app.add_option("-m,--matrix,--matrix1", matrix_file1, "matrix file")->required()
    ->check(CLI::ExistingFile);
  app.add_option("--matrix2", matrix_file2, "matrix file")->check(CLI::ExistingFile);
  const auto outdir_opt = app.add_option("-o,--outdir", output_path,
                                               "output directory");
  app.add_option("--outname", out_filename, "output filename")->needs(outdir_opt);

  try {
    app.parse(app.remaining_for_passthrough());
  } catch(const CLI::ParseError& e) { return app.exit(e); }

  if (!output_path.empty()) {
    fs::create_directories(output_path);
    if (out_filename.empty()) {
      out_filename = matrix_file1.stem();
      if (!matrix_file2.empty()) {
        out_filename += '_' + matrix_file2.stem().string();
      }
      out_filename += "_spGEMM_stats.txt";
    }
    output_path /= out_filename;
  }

  fmt::print("Loading matrix A: {}... ", matrix_file1.string());
  fflush(stdout);
  Spmat_Csr A(matrix_file1);
  fmt::print("Done\n");
  Spmat_Csr B;
  if (matrix_file2.empty()) {
    if (A.num_rows == A.num_cols) {
      fmt::print("Matrix B = A\n");
      B = A;
    } else {
      fmt::print("Matrix B = A^T\n");
      B = A.transpose();
    }
  } else {
    fmt::print("Loading matrix B: {}... ", matrix_file2.string());
    fflush(stdout);
    B = read_matrix_market_file(matrix_file2);
    fmt::print("Done\n");
  }
  fmt::print("Computing spGEMM_stats...\n");
  print_spGEMM_stats(A, B, output_path.string());
  if (!output_path.empty()) {
    fmt::print("Stats written to {}\n", output_path.c_str());
  }
  return 0;
}

int run_gen_app(CLI::App& app) {
  unsigned num_nodes {};
  unsigned num_edges {};
  double A {};
  double B {};
  double C {};
  unsigned seed {};
  fs::path output_path;
  std::string out_filename;

  app.add_option("-n,--num-nodes", num_nodes, "number of nodes")->required();
  app.add_option("-e,--num-edges", num_edges, "number of edges")->required();
  app.add_option("-a", A, "a parameter")->required();
  app.add_option("-b", B, "b parameter")->required();
  app.add_option("-c", C, "c parameter")->required();
  app.add_option("--seed", seed, "seed");
  app.add_option("-o,--outdir", output_path, "output directory")->required();
  app.add_option("--outname", out_filename, "output filename")->required();

  try {
    app.parse(app.remaining_for_passthrough());
  } catch(const CLI::ParseError& e) { return app.exit(e); }

  if (A + B + C >= 1.0) {
    fmt::print("invalid parameters: A + B + C must be smaller than 1.0\n");
    return 0;
  }
  fs::create_directories(output_path);
  output_path /= out_filename;
  gen_RMat(output_path, num_nodes, num_edges, A, B, C, seed);
  return 0;
}

int main(int argc, char* argv[]) {
  try {
    const auto project_info = fmt::format("mergeforest-sim version 1.0.0");
    CLI::App app{project_info};
    app.set_version_flag("--version", project_info);
    app.require_subcommand(1);
    auto stats_app = app.add_subcommand("stats", "Print SpGEMM stats")->prefix_command();
    auto sim_app = app.add_subcommand("simulate", "Run simulation")->prefix_command();
    auto gen_app = app.add_subcommand("generate", "Generate random sparse matrix")
      ->prefix_command();

    CLI11_PARSE(app, argc, argv);

    if (*sim_app) {
      return run_sim_app(*sim_app);
    } else if (*stats_app) {
      return run_stats_app(*stats_app);
    } else if (*gen_app) {
      return run_gen_app(*gen_app);
    }
  } catch (const std::exception& e) {
    spdlog::error("Unhandled exception in main: {}", e.what());
  }
}
