#include <mergeforest-sim/mergeforest/linked_list_cache.hpp>
#include <mergeforest-sim/math_utils.hpp>

#include <algorithm>
#include <stdexcept>
#include <climits>
#include <cassert>
#include <toml/get.hpp>
#include <vector>

namespace mergeforest_sim {

namespace mergeforest {

bool Active_Row::valid() const {
  return num_blocks > 0;
}

bool Inactive_Row::valid() const {
  return B_row_ptr != UINT_MAX;
}

Linked_List_Cache::Linked_List_Cache(const toml::value& parsed_config,
				     const Matrix_Data& matrix_data_)
  : matrix_data{matrix_data_}
  , B_row_ptr_end_fetcher{matrix_data_.preproc_B_row_ptr_end}
{
  get_config_params(parsed_config);
}

void Linked_List_Cache::reset() {
  for (auto& port : mem_ports) { port.reset(); }
  prefetch_port.reset();
  for (auto& port : read_ports) { port.reset(); }
  write_port.reset();
  arbiter = 0;

  B_row_ptr_end_fetcher.reset();
  B_row_ptr_end_fetcher.base_addr = matrix_data.preproc_B_row_ptr_end_addr;
  matB_fetcher.reset();
  pending_reqs.clear();
  for (auto& i : finished_reqs) { i.clear(); }

  active_rows.clear();
  std::ranges::fill(inactive_rows_cache, Inactive_Row{});
  for (unsigned i = 0; i < row_data_list.size() - 1; ++i) {
    row_data_list[i].num_elements = 0;
    row_data_list[i].next = i + 1;
    row_data_list[i].last = false;
  }
  row_data_list.back() = {0, UINT_MAX, true, false};
  free_list_heads = {0};
  C_partial_row_ptr = UINT_MAX;
  inactive_rows_list_head = UINT_MAX;
  inactive_rows_list_tail = UINT_MAX;
  num_inactive_rows = 0;
  num_active_blocks = 0;
  num_inactive_blocks = 0;
  num_C_partial_blocks = 0;
  num_free_blocks = row_data_list.size();
  num_fetching_blocks = 0;
  cycles = 0;

  reads = 0;
  writes = 0;
  B_reads = 0;
  B_elements_read = 0;
  reused_rows = 0;
  fetched_rows = 0;
  evictions = 0;
  num_active_blocks_avg = 0;
  num_inactive_blocks_avg = 0;
  num_C_partial_blocks_avg = 0;
  num_free_blocks_avg  = 0;
  max_free_lists = 0;
  num_samples = 0;
  stats_max_active_rows = 0;
  stats_max_inactive_rows = 0;
  stats_max_fetched_rows = 0;
  stats_max_outstanding_reqs = 0;
}

void Linked_List_Cache::update() {
  // send requests of B matrix data to main memory
  for (unsigned i = 0; i < mem_ports.size() - 1; ++i) {
    if (!mem_ports[i].has_msg_send()) {
      const auto request = matB_fetcher.get_request();
      if (request.valid()) {
        mem_ports[i].add_msg_send(request);
        stats_max_outstanding_reqs = std::max(matB_fetcher.num_outstanding_reqs,
					      stats_max_outstanding_reqs);
	++B_reads;
      }
    }
    mem_ports[i].transfer();
  }
  if (!mem_ports.back().has_msg_send()) {
    const Address addr = B_row_ptr_end_fetcher.get_fetch_address();
    if (addr != invalid_address) {
      mem_ports.back().add_msg_send({.address = addr, .is_write = false});
      ++preproc_A_reads;
    }
  }
  mem_ports.back().transfer();
  // send prefetched B_row_ptrs
  if (!prefetch_port.has_msg_send()) {
    std::vector<Prefetched_Row> prefetched_rows;
    for (unsigned i = 0; i < prefetched_rows_per_cycle; ++i) {
      if (B_row_ptr_end_fetcher.num_elements == 0) break;
      const auto B_row_ptr_end = B_row_ptr_end_fetcher.front();
      const unsigned B_row_head_ptr = add_new_row(B_row_ptr_end.first,
						  B_row_ptr_end.second);
      if (B_row_head_ptr == UINT_MAX) break;
      B_row_ptr_end_fetcher.pop();
      prefetched_rows.emplace_back(B_row_ptr_end.first, B_row_head_ptr);
    }
    if (!prefetched_rows.empty()) {
      prefetch_port.add_msg_send(prefetched_rows);
    }
  }
  prefetch_port.transfer();
  if (cycles == 0) {
    sample_cache_utilization();
  } 
  cycles = inc_mod(cycles, sample_interval);
}

void Linked_List_Cache::apply() {
  // put data from matB fetcher in linked list buffer
  write_B_row_data();
  // receive data from main memory
  for (unsigned i = 0; i < mem_ports.size() - 1; ++i) {
    if (!mem_ports[i].msg_received_valid()) continue;
    matB_fetcher.put_response(mem_ports[i].get_msg_received());
    mem_ports[i].clear_msg_received();
  }
  if (mem_ports.back().msg_received_valid()) {
    const auto mem_response = mem_ports.back().get_msg_received();
    B_row_ptr_end_fetcher.receive_data(mem_response.address);
    mem_ports.back().clear_msg_received();
  }
  receive_read_requests();
  send_read_responses();
  // write C partial data
  if (write_port.msg_received_valid()) {
    const unsigned response = write_C_partial_row(write_port.get_msg_received());
    write_port.clear_msg_received();
    if (response != UINT_MAX) {
      assert(!write_port.has_msg_send());
      write_port.add_msg_send(response);
    }
  }
  // update slave ports
  for (auto& port : read_ports) {
    port.transfer();
  }
  write_port.transfer();
}

Linked_List_Cache::Mem_Port* Linked_List_Cache::get_mem_port(std::size_t id) {
  if (id >= mem_ports.size()) return nullptr;
  return &mem_ports[id];
}

Linked_List_Cache::Prefetch_Port* Linked_List_Cache::get_prefetch_port() {
  return &prefetch_port;
}

Linked_List_Cache::Cache_Read_Port*
Linked_List_Cache::get_read_port(std::size_t id) {
  if (id >= read_ports.size()) return nullptr;
  return &read_ports[id];
}

Linked_List_Cache::Cache_Write_Port* Linked_List_Cache::get_write_port() {
  return &write_port;
}

std::size_t Linked_List_Cache::num_mem_ports() const {
  return mem_ports.size();
}

void Linked_List_Cache::get_config_params(const toml::value& parsed_config) {
  const auto num_mem_ports = toml::find_or(parsed_config, "linked_list_cache",
                                           "num_mem_ports", 4U);
  mem_ports = std::vector<Mem_Port>(num_mem_ports + 1);
  const auto num_cache_read_ports = toml::find<unsigned>(parsed_config,
                                               "merge_tree_manager",
                                               "num_merge_trees");
  read_ports = std::vector<Cache_Read_Port>(num_cache_read_ports);
  finished_reqs.assign(num_cache_read_ports, {});
  const auto max_rows_fetch = toml::find<std::size_t>(parsed_config,
                                                      "linked_list_cache",
                                                      "max_fetched_rows");
  B_row_ptr_end_fetcher.buffer_size = max_rows_fetch;
  matB_fetcher.row_fetchers = std::vector<Row_Fetcher>(max_rows_fetch);
  const auto max_inactive_rows = toml::find_or(parsed_config,
                                               "linked_list_cache",
                                               "max_inactive_rows", 32768U);
  inactive_rows_cache = std::vector<Inactive_Row>(max_inactive_rows);
  num_blocks = toml::find_or(parsed_config, "linked_list_cache",
                             "size",
                             3u * 1024u * 1024u) / block_size_bytes;
  row_data_list = std::vector<Linked_list_Node>(num_blocks);
  max_active_rows = toml::find_or(parsed_config, "linked_list_cache",
                                  "max_active_rows", 1024U);
  inactive_rows_assoc = toml::find_or(parsed_config, "linked_list_cache",
                                      "inactive_rows_assoc", 16u);
  inactive_rows_num_sets = max_inactive_rows / inactive_rows_assoc;
  num_banks = toml::find_or(parsed_config, "linked_list_cache", "num_banks", num_cache_read_ports);
  matB_fetcher.max_outstanding_reqs = toml::find_or(parsed_config,
                                                    "linked_list_cache",
                                                    "max_outstanding_reqs",
                                                    800U);
  prefetched_rows_per_cycle = toml::find_or(parsed_config, "linked_list_cache",
                                            "prefetched_rows_per_cycle", 4U);
  sample_interval = toml::find_or(parsed_config, "linked_list_cache",
                                            "sample_interval", 10000U);
}


unsigned Linked_List_Cache::add_new_row(uint32_t B_row_ptr, uint32_t B_row_end) {
  auto it = active_rows.find(B_row_ptr);
  // search row in the active rows hash table
  if (it != active_rows.end()) {
    ++(it->second.num_uses);
    ++reused_rows;
    return it->second.row_head;
  }
  if (active_rows.size() == max_active_rows) return UINT_MAX;
  // search row in inactive rows cache
  const unsigned index = B_row_ptr % inactive_rows_num_sets;
  for (unsigned i = 0; i < inactive_rows_assoc; ++i) {
    auto& inactive_row = inactive_rows_cache[index * inactive_rows_assoc + i];
    if (inactive_row.B_row_ptr != B_row_ptr) continue;
    // move row to active rows hash table
    active_rows[B_row_ptr] = {inactive_row.row_head, 1, inactive_row.num_blocks};
    stats_max_active_rows = std::max(stats_max_active_rows, active_rows.size());
    num_active_blocks += inactive_row.num_blocks;
    num_inactive_blocks -= inactive_row.num_blocks;
    assert(num_active_blocks + num_inactive_blocks + num_C_partial_blocks
           + num_free_blocks <= row_data_list.size());
    inactive_rows_list_remove(index * inactive_rows_assoc + i);
    ++reused_rows;
    return active_rows[B_row_ptr].row_head;
  }
  // add new row to the cache
  if (!matB_fetcher.can_accept_row()) return UINT_MAX;
  const unsigned row_num_blocks = div_ceil(B_row_end - B_row_ptr, block_size);
  assert(num_free_blocks + num_inactive_blocks >= num_fetching_blocks);
  if (row_num_blocks > num_free_blocks + num_inactive_blocks - num_fetching_blocks) {
    return UINT_MAX;
  }
  const unsigned ptr = allocate_block();
  assert(ptr != UINT_MAX);
  row_data_list[ptr].next = B_row_ptr;
  const Address begin = matrix_data.B_elements_addr + B_row_ptr * element_size;
  const Address end = matrix_data.B_elements_addr + B_row_end * element_size;
  matB_fetcher.add_row(begin, end, ptr);
  stats_max_fetched_rows = std::max(matB_fetcher.num_rows_fetch,
                                    stats_max_fetched_rows);
  active_rows[B_row_ptr] = {ptr, 1, row_num_blocks};
  stats_max_active_rows = std::max(stats_max_active_rows, active_rows.size());
  num_fetching_blocks += row_num_blocks;
  assert(num_free_blocks + num_inactive_blocks >= num_fetching_blocks);
  ++fetched_rows;
  return ptr;
}

void Linked_List_Cache::write_B_row_data() {
  for (auto& row_fetcher : matB_fetcher.row_fetchers) {
    auto [num_elements, ptr, last] = row_fetcher.get_data();
    if (num_elements == 0) continue;
    assert(num_fetching_blocks > 0);
    --num_fetching_blocks;
    assert(num_free_blocks > 0);
    --num_free_blocks;
    ++num_active_blocks;
    assert(num_active_blocks + num_inactive_blocks + num_C_partial_blocks
           + num_free_blocks <= row_data_list.size());
    row_data_list[ptr].num_elements = num_elements;
    B_elements_read += num_elements;
    if (last) {
      --matB_fetcher.num_rows_fetch;
    } else {
      row_data_list[ptr].last = false;
      const unsigned new_block_ptr = allocate_block();
      assert(new_block_ptr != UINT_MAX);
      row_data_list[new_block_ptr].next = row_data_list[ptr].next;
      row_data_list[ptr].next = new_block_ptr;
      row_fetcher.row_ptr = new_block_ptr;
    }
    finish_pending_reqs(ptr);
  }
}

void Linked_List_Cache::receive_read_requests() {
  for (unsigned i = 0; i != read_ports.size(); ++i) {
    if (!read_ports[i].msg_received_valid()) { continue; }
    const auto request = read_ports[i].get_msg_received();
    assert(request.valid());
    auto& row_block = row_data_list[request.row_ptr];
    if (row_block.num_elements == 0
        || (row_block.last == false && row_block.next == UINT_MAX))
    {
      pending_reqs.emplace(request.row_ptr, std::make_pair(i, request.id));
    } else {
      Cache_Response response = {.row_ptr = row_block.next,
                                 .num_elements = row_block.num_elements,
                                 .id = request.id};
      if (row_block.last) { response.row_ptr = UINT_MAX; }
      finished_reqs[i].push_back(response);
      update_cache_block(request.row_ptr);
    }
    read_ports[i].clear_msg_received();
    ++reads;
  }
}

void Linked_List_Cache::send_read_responses() {
  unsigned num_responses = 0;
  for (unsigned i = 0; i < read_ports.size(); ++i) {
    arbiter = inc_mod(arbiter, read_ports.size());
    if (!finished_reqs[arbiter].empty()) {
      assert(!read_ports[arbiter].has_msg_send());
      read_ports[arbiter].add_msg_send(finished_reqs[arbiter].front());
      finished_reqs[arbiter].pop_front(); 
    }
    ++num_responses;
    if (num_responses == num_banks) { break; }
  }
}

void Linked_List_Cache::finish_pending_reqs(unsigned ptr) {
  for (auto [it, end] = pending_reqs.equal_range(ptr); it != end; ++it) {
    Cache_Response response{.row_ptr = row_data_list[ptr].next,
                            .num_elements = row_data_list[ptr].num_elements,
                            .id = it->second.second};
    if (row_data_list[ptr].last) { response.row_ptr = UINT_MAX; }
    finished_reqs[it->second.first].push_back(response);
    update_cache_block(ptr);
  }
  pending_reqs.erase(ptr);
}

void Linked_List_Cache::update_cache_block(unsigned ptr) {
  if (row_data_list[ptr].C_partial_row) {
    // free C partial row block
    row_data_list[ptr].num_elements = 0;
    row_data_list[ptr].C_partial_row = false;
    assert(num_C_partial_blocks != 0);
    --num_C_partial_blocks;
    ++num_free_blocks;
    assert(num_active_blocks + num_inactive_blocks + num_C_partial_blocks
           + num_free_blocks <= row_data_list.size());
    if (free_list_heads.empty()) {
      row_data_list[ptr].next = UINT_MAX;
      row_data_list[ptr].last = true;
      free_list_heads.emplace_back(ptr);
    } else {
      row_data_list[ptr].next = free_list_heads.back();
      row_data_list[ptr].last = false;
      free_list_heads.back() = ptr;
    }
  } else if (row_data_list[ptr].last) {
    // decrease number of uses in active rows
    auto it = active_rows.find(row_data_list[ptr].next);
    assert(it != active_rows.end());
    --(it->second.num_uses);
    if (it->second.num_uses == 0) {
      add_to_inactive_rows(*it);
      active_rows.erase(it);
    }
  }
}

unsigned Linked_List_Cache::write_C_partial_row(Cache_Write request) {
  assert(request.type != Cache_Write::invalid);
  const unsigned new_block_ptr = allocate_block();
  if (new_block_ptr == UINT_MAX) {
    throw std::runtime_error("Linked list cache has no space for partial row\n");
  }
  ++writes;
  assert(num_free_blocks > 0);
  --num_free_blocks;
  ++num_C_partial_blocks;
  assert(num_active_blocks + num_inactive_blocks + num_C_partial_blocks
         + num_free_blocks <= row_data_list.size());
  assert(num_active_blocks + num_C_partial_blocks <= row_data_list.size());
  row_data_list[new_block_ptr].C_partial_row = true;
  row_data_list[new_block_ptr].num_elements = request.num_elements;
  unsigned response = UINT_MAX;
  // first block of partial row
  if (C_partial_row_ptr == UINT_MAX) {
    C_partial_row_ptr = new_block_ptr;
    response = new_block_ptr;
  } else {
    row_data_list[C_partial_row_ptr].next = new_block_ptr;
    finish_pending_reqs(C_partial_row_ptr);
    C_partial_row_ptr = new_block_ptr;
  }
  if (request.type != Cache_Write::write_last) {
    row_data_list[C_partial_row_ptr].last = false;
  } else {
    finish_pending_reqs(C_partial_row_ptr);
    C_partial_row_ptr = UINT_MAX;
  }
  return response;
}

unsigned Linked_List_Cache::allocate_block() {
  if (free_list_heads.empty()) {
    if (!free_inactive_row()) return UINT_MAX;
  }
  const unsigned ptr = free_list_heads.front();
  if (row_data_list[free_list_heads.front()].last) {
    free_list_heads.pop_front();
  } else {
    free_list_heads.front() = row_data_list[free_list_heads.front()].next;
  }
  row_data_list[ptr] = Linked_list_Node{};
  return ptr;
}

bool Linked_List_Cache::free_inactive_row() {
  if (inactive_rows_list_head == UINT_MAX) return false;
  assert(free_list_heads.empty());
  assert(inactive_rows_cache[inactive_rows_list_head].B_row_ptr != UINT_MAX);
  assert(num_inactive_blocks >= inactive_rows_cache[inactive_rows_list_head].num_blocks);
  num_inactive_blocks -= inactive_rows_cache[inactive_rows_list_head].num_blocks;
  num_free_blocks += inactive_rows_cache[inactive_rows_list_head].num_blocks;
  assert(num_active_blocks + num_inactive_blocks + num_C_partial_blocks
         + num_free_blocks <= row_data_list.size());
  assert(num_inactive_rows > 0);
  --num_inactive_rows;
  ++evictions;
  free_list_heads.emplace_back(
    inactive_rows_cache[inactive_rows_list_head].row_head);
  unsigned ptr = inactive_rows_list_head;
  inactive_rows_list_head = inactive_rows_cache[inactive_rows_list_head].next;
  inactive_rows_cache[ptr] = Inactive_Row{};
  if (inactive_rows_list_head == UINT_MAX) {
    assert(num_inactive_rows == 1);
    inactive_rows_list_tail = UINT_MAX;
  } else {
    inactive_rows_cache[inactive_rows_list_head].prev = UINT_MAX;
  }
  max_free_lists = std::max(max_free_lists, free_list_heads.size());
  return true;
}

void Linked_List_Cache::add_to_inactive_rows(
  const std::pair<uint32_t, Active_Row>& active_row)
{
  assert(num_active_blocks >= active_row.second.num_blocks);
  num_active_blocks -= active_row.second.num_blocks;
  num_inactive_blocks += active_row.second.num_blocks;
  assert(num_active_blocks + num_inactive_blocks + num_C_partial_blocks
         + num_free_blocks <= row_data_list.size());
  // search for an empty way or replace smallest inactive row
  const unsigned index = active_row.first % inactive_rows_num_sets;
  unsigned pos = 0;
  unsigned min_row_num_blocks = UINT_MAX;
  for (unsigned i = 0; i < inactive_rows_assoc; ++i) {
    auto& inactive_row = inactive_rows_cache[index * inactive_rows_assoc + i];
    if (!inactive_row.valid()) {
      pos = i;
      break;
    }
    if (inactive_row.num_blocks < min_row_num_blocks) {
      min_row_num_blocks = inactive_row.num_blocks;
      pos = i;
    }
  }
  pos += index * inactive_rows_assoc;
  // add inactive row to free lists queue
  if (inactive_rows_cache[pos].valid()) {
    assert(num_inactive_rows > inactive_rows_cache[pos].num_blocks);
    num_inactive_blocks -= inactive_rows_cache[pos].num_blocks;
    num_free_blocks += inactive_rows_cache[pos].num_blocks;
    assert(num_active_blocks + num_inactive_blocks + num_C_partial_blocks
           + num_free_blocks <= row_data_list.size());
    free_list_heads.emplace_back(inactive_rows_cache[pos].row_head);
    inactive_rows_list_remove(pos);
    ++evictions;
    max_free_lists = std::max(max_free_lists, free_list_heads.size());
  }
  inactive_rows_cache[pos].B_row_ptr = active_row.first;
  inactive_rows_cache[pos].row_head = active_row.second.row_head;
  inactive_rows_cache[pos].num_blocks = active_row.second.num_blocks;
  inactive_rows_cache[pos].prev = inactive_rows_list_tail;
  inactive_rows_cache[pos].next = UINT_MAX;
  if (inactive_rows_list_tail == UINT_MAX) {
    assert(num_inactive_rows == 0);
    assert(inactive_rows_list_head == UINT_MAX);
    inactive_rows_list_head = pos;
  } else {
    inactive_rows_cache[inactive_rows_cache[pos].prev].next = pos;
  }
  inactive_rows_list_tail = pos;
  ++num_inactive_rows;
  stats_max_inactive_rows =
    std::max(stats_max_inactive_rows, num_inactive_rows);
}

void Linked_List_Cache::inactive_rows_list_remove(unsigned ptr) {
  assert(num_inactive_rows > 0);
  if (inactive_rows_cache[ptr].next != UINT_MAX) {
    inactive_rows_cache[inactive_rows_cache[ptr].next].prev =
      inactive_rows_cache[ptr].prev;
  } else {
    inactive_rows_list_tail = inactive_rows_cache[ptr].prev;
  }
  if (inactive_rows_cache[ptr].prev != UINT_MAX) {
    inactive_rows_cache[inactive_rows_cache[ptr].prev].next =
      inactive_rows_cache[ptr].next;
  } else {
    inactive_rows_list_head = inactive_rows_cache[ptr].next;
  }
  inactive_rows_cache[ptr] = Inactive_Row{};
  --num_inactive_rows;
}

void Linked_List_Cache::sample_cache_utilization() {
  num_active_blocks_avg += num_active_blocks;
  num_inactive_blocks_avg += num_inactive_blocks;
  num_C_partial_blocks_avg += num_C_partial_blocks;
  num_free_blocks_avg += num_free_blocks;
  ++num_samples;
}

} // namespace mergeforest

} // namespace mergeforest_sim
