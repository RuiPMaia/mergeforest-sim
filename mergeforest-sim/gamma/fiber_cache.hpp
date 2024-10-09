#ifndef MERGEFOREST_SIM_GAMMA_FIBER_CACHE_HPP
#define MERGEFOREST_SIM_GAMMA_FIBER_CACHE_HPP

#include <cstddef>
#include <cstdint>
#include <mergeforest-sim/matrix_data.hpp>
#include <mergeforest-sim/port.hpp>

#include <toml.hpp>

#include <unordered_set>
#include <vector>
#include <deque>
#include <unordered_map>

namespace mergeforest_sim {

namespace gamma {

struct Pending_Read {
  std::vector<std::pair<std::size_t, unsigned>> dest_ids;
  unsigned num_arrived_reqs {};
  unsigned num_uses {};
  bool C_partial {false};
};

struct Bank {
  std::deque<Mem_Request> mem_reqs;
  std::size_t read_arbiter {UINT64_MAX};
  std::size_t write_arbiter {UINT64_MAX};
};

struct Cache_Line {
  bool valid() const;
  
  Address address {invalid_address};
  unsigned num_uses {};
  bool C_partial {false};
};

class Fiber_Cache {
public:
  using Mem_Port = Port<Mem_Request, Mem_Response>;
  using Slave_Port = Port<Mem_Response, Mem_Request>;
  using Prefetch_Port = Port<Empty_Msg, std::size_t>;

  Fiber_Cache(const toml::value& parsed_config, const Matrix_Data& matrix_data_);
  void reset();
  void update();
  void apply();
  bool inactive();
  Mem_Port* get_mem_port(std::size_t id);
  Slave_Port* get_read_port(std::size_t id);
  Slave_Port* get_write_port(std::size_t id);
  Prefetch_Port* get_prefetch_port();
  // config params
  std::size_t num_blocks {};
  unsigned assoc {};
  unsigned sample_interval {};
  // stats
  std::size_t B_data_reads {};
  std::size_t C_partial_reads {};
  std::size_t C_partial_writes {};
  std::size_t reads {};
  std::size_t writes {};
  std::size_t read_hits {};
  std::size_t B_blocks_avg {};
  std::size_t C_partial_blocks_avg {};
  std::size_t num_samples {};
#ifndef NDEBUG  
  std::unordered_set<Address> C_addrs;
#endif
private:
  void get_config_params(const toml::value& parsed_config);
  void receive_mem_responses();
  void receive_read_requests();
  void process_read_request(std::size_t port); 
  void receive_write_requests();
  void receive_prefetch_data();
  std::size_t cache_search(Address address);
  void cache_insert(Address address, unsigned num_uses, bool C_partial);
  void cache_evict(Address address);
  std::size_t address_to_bank(Address address); 
  void sample_cache_utilization();
  
  const Matrix_Data& matrix_data;

  std::vector<Mem_Port> mem_ports;
  std::vector<Slave_Port> read_ports;
  std::vector<Slave_Port> write_ports;
  Prefetch_Port prefetch_port;
  
  std::size_t mem_arbiter {UINT64_MAX};
  std::size_t prefetch_idx {};
  std::deque<Mem_Request> prefetch_reqs;
  std::vector<Bank> banks;
  std::vector<Cache_Line> cache_lines;
  std::unordered_map<Address, Pending_Read> pending_reqs;
  std::vector<std::deque<Mem_Response>> finished_reqs;
  std::size_t num_B_blocks {};
  std::size_t num_C_partial_blocks {};
  unsigned cycles {};
}; 

} // namespace gamma

} // namespace mergeforest_sim

#endif // MERGEFOREST_SIM_GAMMA_FIBER_CACHE_HPP
