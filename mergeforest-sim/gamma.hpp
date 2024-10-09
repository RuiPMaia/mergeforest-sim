#ifndef MERGEFOREST_SIM_GAMMA_HPP
#define MERGEFOREST_SIM_GAMMA_HPP

#include <mergeforest-sim/gamma/PE_manager.hpp>
#include <mergeforest-sim/gamma/fiber_cache.hpp>
#include <mergeforest-sim/main_memory.hpp>
#include <mergeforest-sim/sparse_matrix.hpp>
#include <mergeforest-sim/matrix_data.hpp>

#include <toml.hpp>

namespace mergeforest_sim {

class Gamma {
public:
  Gamma(const toml::value& parsed_config, Matrix_Data& matrix_data_,
	const std::string& out_path_);
  Spmat_Csr run_simulation(bool compute_result);
private:
  void reset();
  void print_progress();
  void check_valid_simulation();
  void print_stats();
  void print_stats_impl(std::ostream& os);

  const std::size_t progress_interval = 10000;
  const toml::value& parsed_config;
  Matrix_Data& matrix_data;
  const std::string& out_path;
  // system components
  gamma::PE_Manager PE_manager;
  gamma::Fiber_Cache fiber_cache;
  Main_Memory main_mem;
  std::size_t cycles {};
};

} // namespace mergeforest_sim

#endif // MERGEFOREST_SIM_GAMMA_HPP
