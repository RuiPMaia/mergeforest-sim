#ifndef MERGEFOREST_SIM_MY_ARCH_HPP
#define MERGEFOREST_SIM_MY_ARCH_HPP

#include <mergeforest-sim/mergeforest/linked_list_cache.hpp>
#include <mergeforest-sim/mergeforest/merge_tree_manager.hpp>
#include <mergeforest-sim/main_memory.hpp>
#include <mergeforest-sim/matrix_data.hpp>
#include <mergeforest-sim/sparse_matrix.hpp>

#include <toml.hpp>

#include <string>
#include <iostream>

namespace mergeforest_sim {

class MergeForest {
public:
  MergeForest(const toml::value& parsed_config, Matrix_Data& matrix_data_,
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
  mergeforest::Merge_Tree_Manager merge_tree_manager;
  mergeforest::Linked_List_Cache linked_list_cache;
  Main_Memory main_mem;

  std::size_t cycles {};
};

} // namespace mergeforest_sim

#endif // MERGEFOREST_SIM_MY_ARCH_HPP
