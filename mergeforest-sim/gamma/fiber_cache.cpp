#include <mergeforest-sim/gamma/fiber_cache.hpp>
#include <mergeforest-sim/math_utils.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>

namespace mergeforest_sim {

namespace gamma {

bool Cache_Line::valid() const {
  return address != invalid_address;
}

Fiber_Cache::Fiber_Cache(const toml::value& parsed_config, const Matrix_Data& matrix_data_)
  : matrix_data{ matrix_data_ }
{
  get_config_params(parsed_config);    
  reset();
}

void Fiber_Cache::reset() {
  for (auto& i : mem_ports) { i.reset(); }
  for (auto& i : read_ports) { i.reset(); }
  for (auto& i : write_ports) { i.reset(); }
  prefetch_port.reset();
  mem_arbiter = UINT64_MAX;
  prefetch_idx = 0;
  prefetch_reqs.clear();
  std::ranges::fill(banks, Bank{});
  std::ranges::fill(cache_lines, Cache_Line{});
  pending_reqs.clear();
  for (auto& i : finished_reqs) { i.clear(); }
  num_B_blocks = 0;
  num_C_partial_blocks = 0;
  cycles = 0;
  B_data_reads = 0;
  C_partial_reads = 0;
  C_partial_writes = 0;
  reads = 0;
  writes = 0;
  read_hits = 0;
  B_blocks_avg = 0;
  C_partial_blocks_avg = 0;
  num_samples = 0;
  #ifndef NDEBUG  
  C_addrs.clear();
  #endif
}

void Fiber_Cache::update() {
  // send read responses
  for (unsigned i = 0; i < read_ports.size(); ++i) {
    if (read_ports[i].has_msg_send()) continue;
    if (finished_reqs[i].empty()) continue;
    read_ports[i].add_msg_send(finished_reqs[i].front());
    finished_reqs[i].pop_front();
  }
  // Send memory requests with bank misses having priority over prefetch 
  for (auto& p : mem_ports) {
    if (p.has_msg_send()) continue;
    for (unsigned i = 0; i < banks.size(); ++i) {
      mem_arbiter = inc_mod(mem_arbiter, banks.size());
      if (banks[mem_arbiter].mem_reqs.empty()) continue;
      p.add_msg_send(banks[mem_arbiter].mem_reqs.front());
      banks[mem_arbiter].mem_reqs.pop_front();
      break;
    }
    if (p.has_msg_send()) continue;
    // if no request from banks do a prefetch request
    if (prefetch_reqs.empty()) continue;
    p.add_msg_send(prefetch_reqs.front());
    prefetch_reqs.pop_front();
  }
  for (auto& p : read_ports) {
    p.transfer();
  }
  for (auto& p : mem_ports) {
    p.transfer();
  }
  cycles = inc_mod(cycles, sample_interval);
  if (cycles == 0) {
    sample_cache_utilization();
  } 
}

void Fiber_Cache::apply() {
  receive_mem_responses();
  receive_read_requests();
  receive_write_requests();
  receive_prefetch_data();
}

bool Fiber_Cache::inactive() {
  return true;
}

Fiber_Cache::Mem_Port* Fiber_Cache::get_mem_port(std::size_t id) {
  if (id >= mem_ports.size()) return nullptr;
  return &mem_ports[id];
}

Fiber_Cache::Slave_Port* Fiber_Cache::get_read_port(std::size_t id) {
  if (id >= read_ports.size()) return nullptr;
  return &read_ports[id];
}

Fiber_Cache::Slave_Port* Fiber_Cache::get_write_port(std::size_t id) {
  if (id >= write_ports.size()) return nullptr;
  return &write_ports[id];
}

Fiber_Cache::Prefetch_Port* Fiber_Cache::get_prefetch_port() {
  return &prefetch_port;
}

void Fiber_Cache::get_config_params(const toml::value& parsed_config) {
  const auto num_mem_ports = toml::find<std::size_t>(parsed_config, "fiber_cache", "num_mem_ports");
  mem_ports = std::vector<Mem_Port>(num_mem_ports);
  const auto num_slave_ports = toml::find<std::size_t>(parsed_config, "PE_manager", "num_PEs");
  read_ports = std::vector<Slave_Port>(num_slave_ports);
  write_ports = std::vector<Slave_Port>(num_slave_ports);
  finished_reqs = std::vector<std::deque<Mem_Response>>(num_slave_ports);
  num_blocks = toml::find<std::size_t>(parsed_config, "fiber_cache", "size") / block_size_bytes;
  cache_lines = std::vector<Cache_Line>(num_blocks);
  const auto num_banks = toml::find<std::size_t>(parsed_config, "fiber_cache", "num_banks"); 
  banks = std::vector<Bank>(num_banks);
  assoc = toml::find<unsigned>(parsed_config, "fiber_cache", "assoc");
  sample_interval = toml::find_or(parsed_config, "fiber_cache", "sample_interval", 10000U);
}

void Fiber_Cache::receive_mem_responses() {
  for (auto& p : mem_ports) {
    if (!p.msg_received_valid()) continue;
    const auto response = p.get_msg_received();
    const auto addr = round_down_multiple(response.address, static_cast<Address>(block_size_bytes));
    p.clear_msg_received();
    auto it = pending_reqs.find(addr);
    assert(it != pending_reqs.end());
    ++it->second.num_arrived_reqs;
    if (it->second.num_arrived_reqs == 3) {
      for (auto& i : it->second.dest_ids) {
	finished_reqs[i.first].push_back(Mem_Response{.address = addr, .id = i.second});
      }
      if (!it->second.C_partial) {
	cache_insert(addr, it->second.num_uses, false);
      }
      pending_reqs.erase(it);
    } 
  }
}

void Fiber_Cache::receive_read_requests() {
  for (unsigned i = 0; i < banks.size(); ++i) {
    for (unsigned j = 0; j < read_ports.size(); ++j) {
      banks[i].read_arbiter = inc_mod(banks[i].read_arbiter, read_ports.size());
      const auto p = banks[i].read_arbiter;
      if (!read_ports[p].msg_received_valid()) continue;
      auto req = read_ports[p].get_msg_received();
      assert(req.valid());
      if (address_to_bank(req.address) != i) continue;
      process_read_request(p);
      ++reads;
      read_ports[p].clear_msg_received();
    }
  }
}

void Fiber_Cache::process_read_request(std::size_t port) {
  auto req = read_ports[port].get_msg_received();
  // search write ports for the C partial data
  if (req.address >= matrix_data.C_partials_base_addr) {
    for (auto& p: write_ports) {
      if (!p.msg_received_valid()) continue;
      if (req.address == p.get_msg_received().address) {
        p.clear_msg_received();
        finished_reqs[port].push_back(Mem_Response{ .address = req.address, .id = req.id });
        ++read_hits;
        return;
      }
    }
    #ifndef NDEBUG
      assert(C_addrs.contains(req.address));
      C_addrs.erase(req.address);
    #endif
  }
  const auto idx = cache_search(req.address);
  if (idx != UINT_MAX) {
    if (cache_lines[idx].C_partial) {
      assert(req.address >= matrix_data.C_partials_base_addr);
      cache_lines[idx] = Cache_Line{};
      --num_C_partial_blocks;
    } else if (cache_lines[idx].num_uses > 0) {
      --cache_lines[idx].num_uses;
    }
    finished_reqs[port].push_back(Mem_Response{ .address = req.address, .id = req.id });
    ++read_hits;
    return;
  }
  auto it = pending_reqs.find(req.address);
  if (it != pending_reqs.end()) {
    it->second.dest_ids.emplace_back(port, req.id);
    if (it->second.num_uses > 0) --it->second.num_uses;
    return;
  }
  Pending_Read pending_read;
  pending_read.dest_ids.emplace_back(port, req.id);
  if (req.address >= matrix_data.C_partials_base_addr) {
    pending_read.C_partial = true;
    ++C_partial_reads;
  } else {
    ++B_data_reads;
  }
  pending_reqs.emplace(req.address, pending_read);
  for (unsigned k = 0; k < 3; ++k) {
    const auto b = address_to_bank(req.address);
    banks[b].mem_reqs.push_back(Mem_Request{ .address = req.address, .is_write = false });
    req.address += mem_transaction_size;
  }
}

void Fiber_Cache::receive_write_requests() {
  for (unsigned i = 0; i < banks.size(); ++i) {
    for (unsigned j = 0; j < write_ports.size(); ++j) {
      banks[i].write_arbiter = inc_mod(banks[i].write_arbiter, write_ports.size());
      const auto p = banks[i].write_arbiter;
      if (!write_ports[p].msg_received_valid()) continue;
      auto req = write_ports[p].get_msg_received();
      if (address_to_bank(req.address) != i) continue;
      #ifndef NDEBUG
        assert(!C_addrs.contains(req.address));
	C_addrs.emplace(req.address);
      #endif 
      cache_insert(req.address, 1, true);
      write_ports[p].clear_msg_received();
      ++writes;
      break;
    }
  }
}

void Fiber_Cache::receive_prefetch_data() {
  if (!prefetch_port.msg_received_valid()) return;
  auto prefetch_num_elements = prefetch_port.get_msg_received();
  prefetch_port.clear_msg_received();
  while (prefetch_num_elements > 0) {
    auto [B_row_ptr, B_row_end] = matrix_data.preproc_B_row_ptr_end[prefetch_idx];
    ++prefetch_idx;
    --prefetch_num_elements;
    // align B_row_ptr and B_row_end to block size
    B_row_ptr = round_down_multiple(B_row_ptr, block_size);
    B_row_end = round_up_multiple(B_row_end, block_size);
    // create memory requests for all blocks of the row
    while (B_row_ptr < B_row_end) {
      Address addr = matrix_data.B_elements_addr + B_row_ptr * element_size;
      B_row_ptr += block_size;
      const auto cache_idx = cache_search(addr);
      if (cache_idx != UINT_MAX) {
	++cache_lines[cache_idx].num_uses;
	continue;
      }
      const auto it = pending_reqs.find(addr);
      if (it != pending_reqs.end()) {
	++it->second.num_uses;
	continue;
      }
      Pending_Read pending_read;
      pending_read.num_uses = 1;
      pending_reqs.emplace(addr, pending_read);
      for (unsigned i = 0; i < 3; ++i) {
	prefetch_reqs.push_back(Mem_Request{.address = addr, .is_write = false});
	addr += mem_transaction_size;
      }
      ++B_data_reads;
    }
  }
}

std::size_t Fiber_Cache::cache_search(Address address) {
  address = round_down_multiple(address, static_cast<std::size_t>(block_size_bytes));
  const auto index = (address / block_size_bytes) % (cache_lines.size() / assoc);
  for (unsigned i = 0; i < assoc; ++i) {
    const auto idx = index * assoc + i;
    if (cache_lines[idx].address == address) return idx;
  }
  return UINT_MAX;
}

void Fiber_Cache::cache_insert(Address address, unsigned num_uses, bool C_partial) {
  const auto index = (address / block_size_bytes) % (cache_lines.size() / assoc);
  unsigned min_num_uses = UINT_MAX;
  std::size_t min_idx = 0;
  for (unsigned i = 0; i < assoc; ++i) {
    const auto idx = index * assoc + i;
    if (!cache_lines[idx].valid()) {
      cache_lines[idx].address = address;
      cache_lines[idx].num_uses = num_uses;
      cache_lines[idx].C_partial = C_partial;
      if (C_partial) {
        ++num_C_partial_blocks;
      } else {
        ++num_B_blocks;
      }
      return;
    }
    if (min_num_uses > cache_lines[idx].num_uses) {
      min_num_uses = cache_lines[idx].num_uses;
      min_idx = idx;
    }
  }
  if (num_uses > min_num_uses || (C_partial && min_num_uses <= 1)) {
    if (cache_lines[min_idx].C_partial) {
      cache_evict(cache_lines[min_idx].address);
      if (!C_partial) {
        ++num_B_blocks;
        --num_C_partial_blocks;
      }
    } else if (C_partial) {
      ++num_C_partial_blocks;
      --num_B_blocks;
    }
    cache_lines[min_idx].address = address;
    cache_lines[min_idx].num_uses = num_uses;
    cache_lines[min_idx].C_partial = C_partial;
  } else if (C_partial) {
    cache_evict(address);
  }
}

void Fiber_Cache::cache_evict(Address address) {
  const auto bank = address_to_bank(address);
  for (unsigned i = 0; i < 3; ++i) {
    banks[bank].mem_reqs.push_back(Mem_Request{.address = address, .is_write = true});
    address += mem_transaction_size;
  }
  ++C_partial_writes;
}

std::size_t Fiber_Cache::address_to_bank(Address address) {
  return (address / block_size_bytes) % banks.size();
}

void Fiber_Cache::sample_cache_utilization() {
  B_blocks_avg += num_B_blocks;
  C_partial_blocks_avg += num_C_partial_blocks;
  ++num_samples;
}

} // namespace gamma

} // namespace mergeforest_sim
