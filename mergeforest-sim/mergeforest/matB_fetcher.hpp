#ifndef MERGEFOREST_SIM_MAT_B_FETCHER_HPP
#define MERGEFOREST_SIM_MAT_B_FETCHER_HPP

#include <mergeforest-sim/port.hpp>

#include <vector>
#include <deque>
#include <tuple>
#include <climits>

namespace mergeforest_sim {

namespace mergeforest {

struct Row_Fetcher {
  Address row_ptr_addr {invalid_address};
  Address row_end_addr {invalid_address};
  unsigned row_ptr {UINT_MAX};
  std::size_t num_bytes_received {};
  std::deque<std::pair<Address, bool>> pending_reqs;

  std::tuple<unsigned, unsigned, bool> get_data();
};

struct MatB_Fetcher {
  void reset();
  bool add_row(Address begin, Address end, unsigned row_ptr_cache);
  bool can_accept_row() const;
  Mem_Request get_request();
  bool put_response(const Mem_Response& read_response);

  std::vector<Row_Fetcher> row_fetchers;
  std::size_t new_row_idx {};
  std::size_t request_idx {};
  std::size_t num_outstanding_reqs {};
  std::size_t num_rows_fetch {};
  std::size_t max_outstanding_reqs {};
  // stats
  std::size_t bytes_read_B_data {};
};

} // namespace mergeforest

} // namespace mergeforest_sim

#endif // MERGEFOREST_SIM_MAT_B_FETCHER_HPP
