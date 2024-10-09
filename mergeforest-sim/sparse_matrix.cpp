#include <mergeforest-sim/sparse_matrix.hpp>
#include <mergeforest-sim/matrix_IO.hpp>
#include <mergeforest-sim/math_utils.hpp>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include <queue>
#include <vector>
#include <unordered_set>
#include <tuple>
#include <algorithm>
#include <bit>
#include <climits>
#include <iostream>
#include <fstream>

namespace mergeforest_sim {

Spmat_Csr::Spmat_Csr(const std::string& filename)
  : Spmat_Csr{read_matrix_market_file(filename)}
{}

Spmat_Csr Spmat_Csr::transpose() {
  Spmat_Csr B;
  B.num_rows = num_cols;
  B.num_cols = num_rows;
  B.nnz = nnz;
  B.row_ptr = std::vector<uint32_t>(B.num_rows + 1);
  B.col_idx = std::vector<uint32_t>(B.nnz);
  B.values = std::vector<double>(B.nnz);

  std::vector<std::tuple<uint32_t, uint32_t, double>> coo;
  for (unsigned i = 0; i < num_rows; ++i) {
    for (unsigned j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
      coo.emplace_back(col_idx[j], i, values[j]);
    }
  }
  std::sort(coo.begin(), coo.end());
  B.row_ptr[0] = 0;
  uint32_t prev_row = 0;
  for (unsigned i = 0; i < B.nnz; ++i) {
    uint32_t row = std::get<0>(coo[i]);
    B.col_idx[i] = std::get<1>(coo[i]);
    B.values[i] = std::get<2>(coo[i]);
    for (unsigned j = prev_row; j < row; ++j) {
      B.row_ptr[j + 1] = i;
    }
    prev_row = row;
  }
  for (unsigned i = prev_row; i < B.num_rows; ++i) {
    B.row_ptr[i + 1] = static_cast<uint32_t>(B.nnz);
  }
  B.col_idx[B.nnz - 1] = std::get<0>(coo[B.nnz - 1]);
  B.values[B.nnz - 1] = std::get<2>(coo[B.nnz - 1]);
  return B;
}

Spmat_Packed::Spmat_Packed() {
  n_rows = num_sets = 0;
  row_ptr = col_set_idx = nullptr;
  col_set = nullptr;
}

Spmat_Packed::~Spmat_Packed() {
  if(row_ptr != nullptr) delete[] row_ptr;
  if(col_set_idx != nullptr) delete[] col_set_idx;
  if(col_set != nullptr) delete[] col_set;
}

void Spmat_Packed::init(const Spmat_Csr& A) {
  n_rows = A.num_rows;
  row_ptr = new uint32_t[n_rows + 1];
  row_ptr[0] = 0;
  // count number of sets per row
  for(unsigned i = 0; i < A.num_rows; ++i) {
    uint32_t counter = 0;
    unsigned idx = 0;
    for(unsigned j = A.row_ptr[i]; j < A.row_ptr[i+1]; ++j) {
      if(A.col_idx[j] >= idx) {
        ++counter;
        if(A.col_idx[j] % 64 == 0) {
          idx = A.col_idx[j] + 64;
        } else {
          idx = round_up_multiple(A.col_idx[j], 64U);
        }
      }
    }
    row_ptr[i + 1] = row_ptr[i] + counter;
  }
  num_sets = row_ptr[n_rows];
  col_set_idx = new uint32_t[num_sets];
  col_set = new uint64_t[num_sets];
  // calculate sets
  for (unsigned i = 0; i < n_rows; ++i) {
    uint32_t k = row_ptr[i];
    uint32_t idx = 0;
    for (unsigned j = A.row_ptr[i]; j < A.row_ptr[i+1]; ++j) {
      if (A.col_idx[j] >= idx) {
        if (A.col_idx[j] % 64 == 0) {
          idx = A.col_idx[j] + 64;
        } else {
          idx = round_up_multiple(A.col_idx[j], 64U);
        }
        col_set_idx[k] = idx / 64 - 1;
        col_set[k] = 0;
        ++k;
      }
      col_set[k - 1] |= uint64_t{1} << (A.col_idx[j] % 64);
    }
  }
}

void spGEMM_symbolic_phase(const Spmat_Csr& A, const Spmat_Csr& B, Spmat_Csr& C) {
  using pi = std::pair<uint32_t, uint32_t>;
  
  if (A.num_cols != B.num_rows) {
    throw std::runtime_error("matrices A and B don't have compatible dimensions");
  }
  Spmat_Packed B_packed;
  B_packed.init(B);
  C.num_rows = A.num_rows;
  C.num_cols = B.num_cols;
  C.row_ptr = std::vector<uint32_t>(C.num_rows + 1);
  unsigned max_row_size = 0;
  for (unsigned i = 0; i < A.num_rows; ++i) {
    if (max_row_size < A.row_ptr[i+1] - A.row_ptr[i])
      max_row_size = A.row_ptr[i+1] - A.row_ptr[i];
  }
  std::vector<uint32_t> row_idx(max_row_size);
  std::vector<uint32_t> row_end(max_row_size);

  std::priority_queue<pi, std::vector<pi>, std::greater<pi>> heap;
  C.row_ptr[0] = 0;
  for (unsigned i = 0; i < A.num_rows; ++i) {
    unsigned cur_idx = UINT_MAX;
    unsigned counter = 0;
    uint64_t cur_set = 0;
    for (unsigned j = 0; j < A.row_ptr[i+1] - A.row_ptr[i]; ++j) {
      row_idx[j] = B_packed.row_ptr[A.col_idx[A.row_ptr[i] + j]];
      row_end[j] = B_packed.row_ptr[A.col_idx[A.row_ptr[i] + j] + 1];
      if(row_idx[j] < row_end[j]) {
        heap.push(std::make_pair(B_packed.col_set_idx[row_idx[j]], j));
      }
    }
    while (!heap.empty()) {
      pi min = heap.top();
      heap.pop();
      if (min.first == cur_idx) {
        cur_set |= B_packed.col_set[row_idx[min.second]];
      } else {
        if (cur_idx != UINT_MAX) {
          counter += static_cast<unsigned>(std::popcount(cur_set));
        }
        cur_idx = min.first;
        cur_set = B_packed.col_set[row_idx[min.second]];
      }
      ++row_idx[min.second];
      if(row_idx[min.second] < row_end[min.second])
        heap.push(std::make_pair(B_packed.col_set_idx[row_idx[min.second]], min.second));
    }
    counter += static_cast<unsigned>(std::popcount(cur_set));
    C.row_ptr[i + 1] = C.row_ptr[i] + counter;
  }
  C.nnz = C.row_ptr.back();
}

void print_spGEMM_stats(const Spmat_Csr& A, const Spmat_Csr& B, std::string_view out_path) {
  Spmat_Csr C_symbolic_phase;
  spGEMM_symbolic_phase(A, B, C_symbolic_phase);
  std::size_t num_mults {0};
  std::size_t A_max_row_size {0};
  std::size_t A_min_row_size {A.num_rows};
  std::size_t B_max_row_size {0};
  std::size_t B_min_row_size {B.num_rows};
  std::size_t rows_to_process {0};
  std::size_t A_data_num_elements {0};
  std::size_t min_bytes_B_data {0};
  std::unordered_set<uint32_t> B_row_set;

  for (std::size_t i = 0; i < A.num_rows; ++i) {
    unsigned non_empty_rows {0};
    for (std::size_t j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
      const std::size_t B_row_size = B.row_ptr[A.col_idx[j] + 1] - B.row_ptr[A.col_idx[j]];
      if (B_row_size > 0) {
         if (!B_row_set.contains(A.col_idx[j])) {
          B_row_set.insert(A.col_idx[j]);
          min_bytes_B_data += B_row_size;
        }
        ++non_empty_rows;
        num_mults += B_row_size;
      }
      B_max_row_size = std::max(B_max_row_size, B_row_size);
      B_min_row_size = std::min(B_min_row_size, B_row_size);
    }
    const std::size_t A_row_size = A.row_ptr[i + 1] - A.row_ptr[i];
    A_max_row_size = std::max(A_max_row_size, A_row_size);
    A_min_row_size = std::min(A_min_row_size, A_row_size);
    if (non_empty_rows > 0) {
      ++rows_to_process;
      A_data_num_elements += non_empty_rows;
    }
  }

  std::ofstream of;
  std::ostream* os {};
  if (out_path.empty()) {
    os = &std::cout;
  } else {
    of.open(out_path.data());
    os = &of;
  }

  fmt::print(*os, "*---Matrix A---*\n");
  fmt::print(*os, "dimensions: {}x{}\n", A.num_rows, A.num_cols);
  fmt::print(*os, "nnz: {}\n", A.nnz);
  fmt::print(*os, "density: {:.4e}\n",
             (static_cast<double>(A.nnz) / static_cast<double>(A.num_rows))
             / static_cast<double>(A.num_cols));
  fmt::print(*os, "avg nnz per row: {:.4f}\n",
             static_cast<double>(A.nnz) / static_cast<double>(A.num_rows));
  fmt::print(*os, "max nnz per row: {}\n", A_max_row_size);
  fmt::print(*os, "min nnz per row: {}\n", A_min_row_size);
  fmt::print(*os, "*---Matrix B---*\n");
  fmt::print(*os, "dimensions: {}x{}\n", B.num_rows, B.num_cols);
  fmt::print(*os, "nnz: {}\n", B.nnz);
  fmt::print(*os, "density: {:.4e}\n",
             (static_cast<double>(B.nnz) / static_cast<double>(B.num_rows))
             / static_cast<double>(B.num_cols));
  fmt::print(*os, "avg nnz per row: {:.4f}\n",
             static_cast<double>(B.nnz) / static_cast<double>(B.num_rows));
  fmt::print(*os, "max nnz per row: {}\n", B_max_row_size);
  fmt::print(*os, "min nnz per row: {}\n", B_min_row_size);
  fmt::print(*os, "*---SpGEMM---*\n");
  fmt::print(*os, "number of mults: {}\n", num_mults);
  fmt::print(*os, "number of adds: {}\n", num_mults - C_symbolic_phase.nnz);
  fmt::print(*os, "nnz of result: {}\n", C_symbolic_phase.nnz);
  fmt::print(*os, "compression factor (n_mults/result nnz): {:.4f}\n",
             static_cast<double>(num_mults) / static_cast<double>(C_symbolic_phase.nnz));
  std::size_t A_bytes = rows_to_process * 3 * sizeof(uint32_t)
    + A_data_num_elements * (sizeof(double) + 2 * sizeof(uint32_t));
  std::size_t C_bytes = C_symbolic_phase.nnz * (sizeof(uint32_t) + sizeof(double));
  std::size_t B_max_bytes = num_mults * (sizeof(uint32_t) + sizeof(double));
  min_bytes_B_data *= (sizeof(uint32_t) + sizeof(double));
  fmt::print(*os, "A data bytes: {} ({:.4f} MB)\n", A_bytes, static_cast<double>(A_bytes) * 1E-6);
  fmt::print(*os, "C data bytes: {} ({:.4f} MB)\n", C_bytes, static_cast<double>(C_bytes) * 1E-6);
  fmt::print(*os, "B compulsory data bytes: {} ({:.4f} MB)\n",
             min_bytes_B_data, static_cast<double>(min_bytes_B_data) * 1E-6);
  fmt::print(*os, "B maximum data bytes: {} ({:.4f} MB)\n", B_max_bytes,
             static_cast<double>(B_max_bytes) * 1E-6);
  fmt::print(*os, "operational intensity (no B row reuse): {:.4f} flops/byte\n",
             static_cast<double>(num_mults) / static_cast<double>(A_bytes + B_max_bytes + C_bytes));
  fmt::print(*os, "operational intensity (full B row reuse): {:.4f} flops/byte\n",
             static_cast<double>(num_mults) / static_cast<double>(A_bytes + min_bytes_B_data + C_bytes));
}

} // namespace mergeforest_sim
