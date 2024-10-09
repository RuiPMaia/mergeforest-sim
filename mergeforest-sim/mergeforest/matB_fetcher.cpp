#include <mergeforest-sim/mergeforest/matB_fetcher.hpp>
#include <mergeforest-sim/math_utils.hpp>

#include <algorithm>
#include <cassert>

namespace mergeforest_sim {

namespace mergeforest {

std::tuple<unsigned, unsigned, bool> Row_Fetcher::get_data() {
  if (row_ptr_addr == invalid_address) return {0, UINT_MAX, false};
  bool last = row_ptr_addr == row_end_addr && pending_reqs.empty()
    && num_bytes_received <= block_size_bytes;
  if (last) {
    const unsigned num_elements = static_cast<unsigned>(num_bytes_received)
      / element_size;
    row_ptr_addr = invalid_address;
    row_end_addr = invalid_address;
    num_bytes_received = 0;
    return {num_elements, row_ptr, true};
  } else if (num_bytes_received >= block_size_bytes) {
    num_bytes_received -= block_size_bytes;
    return {block_size, row_ptr, false};
  }
  return {0, UINT_MAX, false};
}

void MatB_Fetcher::reset() {
  std::fill(row_fetchers.begin(), row_fetchers.end(), Row_Fetcher{});
  new_row_idx = 0;
  request_idx = 0;
  num_outstanding_reqs = 0;
  num_rows_fetch = 0;
  bytes_read_B_data = 0;
}

bool MatB_Fetcher::add_row(Address begin, Address end, unsigned row_ptr_cache) {
  if (!can_accept_row()) return false;
  for (unsigned i = 0; i < row_fetchers.size(); ++i) {
    new_row_idx = inc_mod(new_row_idx, row_fetchers.size());
    if (row_fetchers[new_row_idx].row_ptr_addr == invalid_address) {
      row_fetchers[new_row_idx].row_ptr = row_ptr_cache;
      row_fetchers[new_row_idx].row_ptr_addr = begin;
      row_fetchers[new_row_idx].row_end_addr = end;
      ++num_rows_fetch;
      return true;
    }
  }
  assert(false);
  return false;
}

bool MatB_Fetcher::can_accept_row() const {
  return num_rows_fetch < row_fetchers.size();
}

Mem_Request MatB_Fetcher::get_request() {
  if (num_outstanding_reqs == max_outstanding_reqs) { return Mem_Request{}; }
  for (unsigned i = 0; i < row_fetchers.size(); ++i) {
    request_idx = inc_mod(request_idx, row_fetchers.size());
    if (row_fetchers[request_idx].row_ptr_addr
        < row_fetchers[request_idx].row_end_addr)
    {
      Mem_Request request{.address = row_fetchers[request_idx].row_ptr_addr,
                          .id = static_cast<unsigned>(request_idx),
                          .is_write = false};
      row_fetchers[request_idx].pending_reqs.emplace_back(
        row_fetchers[request_idx].row_ptr_addr, false);
      const auto num_bytes =
        std::min(mem_transaction_size
                 - row_fetchers[request_idx].row_ptr_addr
                 % mem_transaction_size,
                 row_fetchers[request_idx].row_end_addr
                 - row_fetchers[request_idx].row_ptr_addr);
      row_fetchers[request_idx].row_ptr_addr += num_bytes;
      ++num_outstanding_reqs;
      bytes_read_B_data += num_bytes;
      return request;
    }
  }
  return Mem_Request{};
}

bool MatB_Fetcher::put_response(const Mem_Response& read_response) {
  if (!read_response.valid()) return false;
  assert(!row_fetchers[read_response.id].pending_reqs.empty());
  for (auto& req : row_fetchers[read_response.id].pending_reqs) {
    if (req.first == read_response.address) {
      req.second = true;
      break;
    }
  }
  while (!row_fetchers[read_response.id].pending_reqs.empty()) {
    if (row_fetchers[read_response.id].pending_reqs.front().second) {
      const Address address =
        row_fetchers[read_response.id].pending_reqs.front().first;
      row_fetchers[read_response.id].num_bytes_received +=
        static_cast<unsigned>(
          std::min(mem_transaction_size - address % mem_transaction_size,
                   row_fetchers[read_response.id].row_end_addr - address));
      row_fetchers[read_response.id].pending_reqs.pop_front();
    }
    else break;
  }
  --num_outstanding_reqs;
  return true;
}

} // namespace mergeforest

} // namespace mergeforest_sim
