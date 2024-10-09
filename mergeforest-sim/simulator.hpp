#ifndef MERGEFOREST_SIM_SIMULATOR_HPP
#define MERGEFOREST_SIM_SIMULATOR_HPP

#include <mergeforest-sim/matrix_data.hpp>
#include <mergeforest-sim/sparse_matrix.hpp>
#include <mergeforest-sim/mergeforest.hpp>
#include <mergeforest-sim/gamma.hpp>

#include <toml.hpp>

#include <string>
#include <memory>
#include <variant>

namespace mergeforest_sim {

class Simulator {
public:
  Simulator(const std::string& config_file, const std::string& out_path_ = {});
  void set_mats(const Spmat_Csr& A, const Spmat_Csr& B);  
  Spmat_Csr run_simulation(bool compute_result = false);
private:
  const toml::value parsed_config;
  Matrix_Data matrix_data;
  std::string out_path; 

  std::variant<std::monostate, MergeForest, Gamma> arch;
};
  
} // namespace mergeforest_sim

#endif // MERGEFOREST_SIM_SIMULATOR_HPP
