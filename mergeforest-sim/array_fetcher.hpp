#ifndef MERGEFOREST_SIM_ARRAY_FETCHER_HPP
#define MERGEFOREST_SIM_ARRAY_FETCHER_HPP

#include <mergeforest-sim/port.hpp>

#include <vector>
#include <deque>
#include <cassert>

namespace mergeforest_sim {

template<typename T>
class Array_Fetcher {
public:
  Array_Fetcher(const std::vector<T>& vec_)
    : vec {vec_}
  {}

  void reset() {
    idx = 0;
    idx_fetch = 0;
    num_elements = 0;
    pending_reqs.clear();
  }

  Address get_fetch_address() {
    if (idx_fetch >= vec.size()) return invalid_address;
    if (idx_fetch - idx > buffer_size - mem_transaction_size / sizeof(T)) return invalid_address;
    Address address = base_addr + idx_fetch * sizeof(T);
    pending_reqs.emplace_back(address, false);
    idx_fetch += mem_transaction_size / sizeof(T);
    return address;
  }

  std::size_t receive_data(Address address) {
    if (address == invalid_address) return 0;
    assert(!pending_reqs.empty());
    for (auto& req : pending_reqs) {
      if (req.first == address) {
        req.second = true;
        break;
      }
    }
    std::size_t total_elements_received {};
    while (!pending_reqs.empty()) {
      if (pending_reqs.front().second) {
        const auto num_elements_received = 
          std::min(mem_transaction_size / sizeof(T), vec.size() - num_elements - idx);
        num_elements += num_elements_received;
	total_elements_received += num_elements_received;
        assert(num_elements <= buffer_size);
        pending_reqs.pop_front();
      } else break;
    }
    return total_elements_received;
  }

  bool finished() const {
    return idx == vec.size();
  }

  const T& front() const {
    return vec[idx];
  }

  const T& at(std::size_t pos) const {
    return vec[idx + pos];
  }

  void pop() {
     if (num_elements == 0) return;
    ++idx;
    --num_elements;
  }

  std::size_t buffer_size {};
  Address base_addr {invalid_address};
  std::size_t num_elements {};
private:
  const std::vector<T>& vec;
  std::size_t idx {};
  std::size_t idx_fetch {};
  std::deque<std::pair<Address, bool>> pending_reqs;
};

} // namespace mergeforest_sim

#endif // MERGEFOREST_SIM_ARRAY_FETCHER_HPP
