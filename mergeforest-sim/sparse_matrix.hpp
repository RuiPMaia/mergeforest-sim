#ifndef MERGEFOREST_SIM_SP_MAT_HPP
#define MERGEFOREST_SIM_SP_MAT_HPP

#include <string>
#include <vector>
#include <cstdint>

namespace mergeforest_sim {

struct Spmat_Csr {
  Spmat_Csr() = default;
  explicit Spmat_Csr(const std::string& filename);

  Spmat_Csr transpose();

  uint32_t num_rows {0};
  uint32_t num_cols {0};
  std::size_t nnz {0};
  std::vector<uint32_t> row_ptr;
  std::vector<uint32_t> row_end;
  std::vector<uint32_t> col_idx;
  std::vector<double> values;
};

struct Spmat_Packed {
  uint32_t n_rows;
  uint32_t num_sets;
  uint32_t* row_ptr;
  uint32_t* col_set_idx;
  uint64_t* col_set;
  Spmat_Packed();
  ~Spmat_Packed();
  void init(const Spmat_Csr& A);
};

void spGEMM_symbolic_phase(const Spmat_Csr& A, const Spmat_Csr& B, Spmat_Csr& C);

void print_spGEMM_stats(const Spmat_Csr& A, const Spmat_Csr& B, std::string_view out_path);
  
} // namespace mergeforest_sim

#endif // MERGEFOREST_SIM_SP_MAT_HPP
