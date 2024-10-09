#ifndef MERGEFOREST_SIM_MAIN_MEM_HPP
#define MERGEFOREST_SIM_MAIN_MEM_HPP

#include <mergeforest-sim/port.hpp>
#include <mergeforest-sim/math_utils.hpp>

#include <toml.hpp>

#include <string>
#include <deque>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <memory>
#include <cstdint>

namespace mergeforest_sim {

class Main_Memory {
public:
  using Mem_Port = Port<Mem_Response, Mem_Request>;

  Main_Memory(const toml::value& parsed_config);
  void reset();
  void update();
  void set_num_ports(std::size_t num_ports);
  Mem_Port* get_port(std::size_t id);
  bool inactive() const;
  void print_dramsim3_stats() const;

  // stats
  std::size_t read_requests {};
  std::size_t write_requests {};
  std::size_t reads_completed {};
  std::size_t writes_completed {};
private:
  void get_config_params(const toml::value& parsed_config);

  std::vector<Mem_Port> slave_ports;
  std::deque<std::tuple<Mem_Response, std::size_t, std::size_t>> pending_reqs;
  std::size_t arbiter {UINT64_MAX};
  std::size_t cycle {};
  // config parameters
  unsigned latency {};
  unsigned requests_per_cycle {};
};

} // namespace mergeforest_sim

#endif // MERGEFOREST_SIM_MAIN_MEM_HPP
