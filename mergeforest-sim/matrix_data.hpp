#ifndef MERGEFOREST_SIM_MAT_DATA_HPP
#define MERGEFOREST_SIM_MAT_DATA_HPP

#include <cstddef>
#include <mergeforest-sim/port.hpp>
#include <mergeforest-sim/sparse_matrix.hpp>

#include <vector>
#include <utility>

namespace mergeforest_sim {

struct Matrix_Data {
  void preprocess_mats();
  void set_physical_addrs();
  bool spGEMM_check_result();
  // pointers to matrix objects
  const Spmat_Csr* A {nullptr};
  const Spmat_Csr* B = {nullptr};
  // result matrix
  Spmat_Csr C;
  bool compute_result {};
  // preprocessed arrays
  std::vector<uint32_t> preproc_A_row_ptr;
  std::vector<uint32_t> preproc_A_row_idx;
  std::vector<uint32_t> preproc_C_row_ptr;
  std::vector<double> preproc_A_values;
  std::vector<std::pair<uint32_t, uint32_t>> preproc_B_row_ptr_end;
  // physical addresses of the matrix arrays
  Address B_elements_addr {invalid_address};
  Address C_row_ptr_addr {invalid_address};
  Address C_row_end_addr {invalid_address};
  Address C_elements_addr {invalid_address};
  Address preproc_A_row_ptr_addr {invalid_address};
  Address preproc_A_row_idx_addr {invalid_address};
  Address preproc_A_values_addr {invalid_address};
  Address preproc_B_row_ptr_end_addr {invalid_address};
  Address C_partials_base_addr {invalid_address};
  //min and max B data bytes needed
  std::size_t B_data_min_reads {};
  std::size_t B_data_max_reads {};
  std::size_t B_data_min_reads_fiber_cache {};
  std::size_t B_data_max_reads_fiber_cache {};
  std::size_t min_bytes_B_data {};
  std::size_t max_bytes_B_data {};
  std::size_t num_mults {};
};

} // namespace mergeforest_sim

#endif // MERGEFOREST_SIM_MAT_DATA_HPP
