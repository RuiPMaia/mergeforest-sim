#include <cstddef>
#include <mergeforest-sim/matrix_data.hpp>
#include <mergeforest-sim/math_utils.hpp>
#include <mergeforest-sim/sparse_matrix.hpp>

#include <fmt/format.h>

#include <map>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <cmath>

namespace mergeforest_sim {

unsigned row_num_reads(unsigned B_row_ptr, unsigned B_row_end) {
  const auto begin_addr = round_down_multiple(B_row_ptr * element_size, mem_transaction_size);
  const auto end_addr = round_up_multiple(B_row_end * element_size, mem_transaction_size);
  return (end_addr - begin_addr) / mem_transaction_size; 
}

unsigned row_num_reads_fiber_cache(unsigned B_row_ptr, unsigned B_row_end) {
  const auto begin_addr = round_down_multiple(B_row_ptr, block_size);
  const auto end_addr = round_up_multiple(B_row_end, block_size);
  return (end_addr - begin_addr) / block_size; 
}

void Matrix_Data::preprocess_mats() {
  if (A->num_cols != B->num_rows) {
    throw std::runtime_error("matrices A and B don't have compatible dimensions");
  }
  fmt::print("Allocating space for result matrix using the upper-bound method... ");
  fflush(stdout);

  C.num_rows = A->num_rows;
  C.num_cols = B->num_cols;
  C.row_ptr = std::vector<uint32_t>(C.num_rows + 1);
  C.row_end = std::vector<uint32_t>(C.num_rows);
  C.row_ptr[0] = 0;

  B_data_max_reads = 0;
  B_data_min_reads = 0;
  B_data_min_reads_fiber_cache = 0;
  min_bytes_B_data = 0;
  max_bytes_B_data = 0;
  num_mults = 0;
  preproc_A_row_ptr.push_back(0);
  preproc_A_values.reserve(A->nnz);
  preproc_B_row_ptr_end.reserve(A->nnz);
  std::unordered_set<unsigned> B_row_set;
  std::unordered_set<unsigned> B_cache_block_set;
  bool C_row_ptr_overflow {false};

  for (unsigned i = 0; i < A->num_rows; ++i) {
    unsigned C_max_row_size {0};
    unsigned non_empty_rows {0};
    for (unsigned j = A->row_ptr[i]; j < A->row_ptr[i + 1]; ++j) {
      const unsigned B_row_ptr = B->row_ptr[A->col_idx[j]];
      const unsigned B_row_end = B->row_ptr[A->col_idx[j] + 1];
      const unsigned B_row_size = B_row_end - B_row_ptr;
      if (B_row_size == 0) continue;
      max_bytes_B_data += B_row_size;
      const auto B_row_num_reads = row_num_reads(B_row_ptr, B_row_end);
      B_data_max_reads += B_row_num_reads;
      B_data_max_reads_fiber_cache += row_num_reads_fiber_cache(B_row_ptr, B_row_end);
      if (!B_row_set.contains(A->col_idx[j])) {
        B_row_set.insert(A->col_idx[j]);
        min_bytes_B_data += B_row_size;
        B_data_min_reads += B_row_num_reads;
      }
      // iterate all blocks of a row to calculate min B reads with infinite cache
      for (unsigned idx = round_down_multiple(B_row_ptr, block_size); idx < B_row_end; idx += block_size) {
        if (!B_cache_block_set.contains(idx)) {
          B_cache_block_set.insert(idx);
          ++B_data_min_reads_fiber_cache;
        }
      }
      C_max_row_size += B_row_size;
      num_mults += B_row_size;
      ++non_empty_rows;
      preproc_A_values.emplace_back(A->values[j]);
      preproc_B_row_ptr_end.emplace_back(B_row_ptr, B_row_end);
    }
    C_max_row_size = std::min(C_max_row_size, B->num_cols);
    C.row_ptr[i + 1] = C.row_ptr[i] + C_max_row_size;
    C.row_end[i] =  C.row_ptr[i];
    if (C.row_ptr[i + 1] < C.row_ptr[i]) {
      C_row_ptr_overflow = true;
    }
    if (non_empty_rows > 0) {
      preproc_A_row_ptr.push_back(preproc_A_row_ptr.back() + non_empty_rows);
      preproc_A_row_idx.push_back(i);
      preproc_C_row_ptr.push_back(C.row_ptr[i]);
    }
  }
  B_data_min_reads_fiber_cache *= 3;
  B_data_max_reads_fiber_cache *= 3;
  fmt::print("Done\n");
  if (C_row_ptr_overflow) {
    fmt::print("Not enough space for the upper-bound method. Performing symbolic phase... ");
    fflush(stdout);
    spGEMM_symbolic_phase(*A, *B, C);
    C.row_end.clear();
    fmt::print("Done\n");
  }
  if (compute_result) {
    C.col_idx = std::vector<uint32_t>(C.row_ptr[C.num_rows]);
    C.values = std::vector<double>(C.row_ptr[C.num_rows]);
  }
  min_bytes_B_data *= (sizeof(int) + sizeof(double));
  max_bytes_B_data *= (sizeof(int) + sizeof(double));
}

void Matrix_Data::set_physical_addrs() {
  Address addr {0UL};
  B_elements_addr = addr;
  addr += round_up_multiple(B->nnz * (sizeof(int) + sizeof(double)), 96UL);
  C_row_ptr_addr = addr;
  addr += round_up_multiple((C.num_rows + 1) * sizeof(int), 32UL);
  C_row_end_addr = addr;
  addr += round_up_multiple(C.num_rows * sizeof(int), 32UL);
  C_elements_addr = addr;
  addr += round_up_multiple(C.row_ptr[C.num_rows] * (sizeof(int) + sizeof(double)), 96UL);
  preproc_A_row_ptr_addr = addr;
  addr += round_up_multiple(preproc_A_row_ptr.size() * sizeof(int), 32UL);
  preproc_A_row_idx_addr = addr;
  addr += round_up_multiple(preproc_A_row_idx.size() * sizeof(int), 32UL);
  preproc_A_values_addr = addr;
  addr += round_up_multiple(preproc_A_values.size() * sizeof(double), 32UL);
  preproc_B_row_ptr_end_addr = addr;
  addr += round_up_multiple(preproc_B_row_ptr_end.size() * 2 * sizeof(uint32_t), 32UL);
  C_partials_base_addr = round_up_multiple(addr, 96ul);
}

bool Matrix_Data::spGEMM_check_result() {
  using pi = std::pair<uint32_t, uint32_t>;

  fmt::print("Checking result... ");
  fflush(stdout);
  uint32_t max_row_size {0};
  for (unsigned i = 0; i < A->num_rows; ++i) {
    if (max_row_size < A->row_ptr[i+1] - A->row_ptr[i]) {
      max_row_size = A->row_ptr[i+1] - A->row_ptr[i];
    }
  }
  std::vector<uint32_t> B_row_addr(max_row_size);
  std::vector<uint32_t> B_row_end(max_row_size);
  std::vector<double> A_values(max_row_size);
  std::priority_queue<pi, std::vector<pi>, std::greater<pi>> heap;

  for (unsigned i = 0; i < A->num_rows; ++i) {
    uint32_t cur_idx = UINT32_MAX;
    double cur_value = 0.0;
    uint32_t offset = C.row_ptr[i];
    for (unsigned j = 0; j < A->row_ptr[i+1] - A->row_ptr[i]; ++j) {
      B_row_addr[j] = B->row_ptr[A->col_idx[A->row_ptr[i] + j]];
      B_row_end[j] = B->row_ptr[A->col_idx[A->row_ptr[i] + j] + 1];
      A_values[j] = A->values[A->row_ptr[i] + j];
      if (B_row_addr[j] < B_row_end[j]) {
        heap.push(std::make_pair(B->col_idx[B_row_addr[j]], j));
      }
    }
    while (!heap.empty()) {
      pi min = heap.top();
      heap.pop();
      if(min.first == cur_idx) {
        //cur_value += A_values[min.second] * B->values[B_row_addr[min.second]];
        cur_value = std::fma(A_values[min.second], B->values[B_row_addr[min.second]], cur_value);
      }
      else {
        if(cur_idx != UINT32_MAX) {
          if (C.col_idx[offset] != cur_idx || !almost_equal(C.values[offset], cur_value, 1e6)) {
            fmt::print("\nError in row {}: {}, {} should be {}, {}\n",
              i, C.col_idx[offset], C.values[offset], cur_idx, cur_value);
            return false;
          }
          ++offset;
        }
        cur_idx = min.first;
        cur_value = A_values[min.second] * B->values[B_row_addr[min.second]];
      }
      ++B_row_addr[min.second];
      if(B_row_addr[min.second] < B_row_end[min.second])
        heap.push(std::make_pair(B->col_idx[B_row_addr[min.second]], min.second));
    }
    if (cur_idx != UINT32_MAX) {
      if (C.col_idx[offset] != cur_idx || !almost_equal(C.values[offset], cur_value, 1e6)) {
        fmt::print("\nError in row {}: {}, {} should be {}, {}\n",
          i, C.col_idx[offset], C.values[offset], cur_idx, cur_value);
        return false;
      }
      ++offset;
    }
    if (offset != C.row_end[i]) {
      fmt::print("\nError in row end {}: {} should be {}\n", i, C.row_end[i], offset);
      return false;
    }
  }
  fmt::print("Correct!\n");
  return true;
}

} // namespace mergeforest_sim
