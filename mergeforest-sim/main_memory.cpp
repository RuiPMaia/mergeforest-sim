#include <mergeforest-sim/main_memory.hpp>

#include <functional>
#include <cassert>

namespace mergeforest_sim {

Main_Memory::Main_Memory(const toml::value& parsed_config) {
  get_config_params(parsed_config);
}

void Main_Memory::reset() {
  for (auto& port : slave_ports) {
    port.reset();
  }
  pending_reqs.clear();
  arbiter = UINT64_MAX;
  cycle = 0;
  read_requests = 0;
  write_requests = 0;
  reads_completed = 0;
  writes_completed = 0;
}

void Main_Memory::update() {
  // add requests with round robin arbitration and respecting the maximum bandwidth
  unsigned count {0};
  for (unsigned i = 0; i < slave_ports.size(); ++i) {
    arbiter = inc_mod(arbiter, slave_ports.size());
    if (!slave_ports[arbiter].msg_received_valid()) continue;
    const auto request = slave_ports[arbiter].get_msg_received();
    assert(request.valid());
    if (request.is_write) {
      ++write_requests;
      ++writes_completed;
    } else {
      pending_reqs.emplace_back(Mem_Response{.address = request.address, .id = request.id},
        cycle + latency, arbiter);
      ++read_requests;
    }
    slave_ports[arbiter].clear_msg_received();
    ++count;
    if (count == requests_per_cycle) break;
  }
  // send responses that have waited the specified latency
  while (!pending_reqs.empty()) {
    auto [resp, req_cycle, idx] = pending_reqs.front();
    if (req_cycle <= cycle && !slave_ports[idx].has_msg_send()) {
      slave_ports[idx].add_msg_send(resp);
      pending_reqs.pop_front();
      ++reads_completed;
    } else break;
  }
  ++cycle;
  for (auto& port : slave_ports) {
    port.transfer();
  }
}

void Main_Memory::set_num_ports(std::size_t num_ports) {
  slave_ports = std::vector<Mem_Port>(num_ports);
}

Main_Memory::Mem_Port* Main_Memory::get_port(std::size_t id) {
  if ( id >= slave_ports.size()) return nullptr;
  return &slave_ports[id];
}

bool Main_Memory::inactive() const {
  return read_requests == reads_completed
    && write_requests == writes_completed;
};

void Main_Memory::get_config_params(const toml::value& parsed_config) {
  latency = toml::find_or(parsed_config, "mem", "latency", 80u);
  requests_per_cycle = toml::find_or(parsed_config, "mem", "bandwidth", 128u) / mem_transaction_size;
}

} // namespace mergeforest_sim
