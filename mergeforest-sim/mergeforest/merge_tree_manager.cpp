#include <mergeforest-sim/mergeforest/merge_tree_manager.hpp>
#include <mergeforest-sim/math_utils.hpp>

#include <fmt/format.h>

#include <vector>
#include <algorithm>
#include <cassert>

namespace mergeforest_sim {

namespace mergeforest {


bool Fiber_Buffer::empty() const {
  return col_idx.empty();
}

bool Fiber_Buffer::finished() const {
  return col_idx.empty() && last;
}

std::size_t Fiber_Buffer::size() const {
  return col_idx.size();
}

bool Fiber_Buffer::ready_to_merge(unsigned size) const {
  return last || col_idx.size() >= size;
}

bool C_Partial_Fiber::finished() const {
  return data.finished();
}

bool Task_Output::valid() const {
  return C_partial || write_address != invalid_address;
}

Address Task_Output::get_C_write_address() {
  if (write_address == invalid_address || num_bytes_write == 0) { return invalid_address; }
  const auto ret_address = write_address;
  const auto write_size = mem_transaction_size - write_address % mem_transaction_size;
  if (num_bytes_write >= write_size) {
    num_bytes_write -= write_size;
    if (num_bytes_write == 0 && C_row_idx == UINT_MAX) {
      write_address = invalid_address;
      return ret_address;
    }
    write_address += write_size;
    return ret_address;
  }
  if (C_row_idx != UINT_MAX) { return invalid_address; }
  num_bytes_write = 0;
  write_address = invalid_address;
  return ret_address;
}

Cache_Write Task_Output::get_C_partial_write() {
  if (!C_partial) { return Cache_Write{}; }
  if (num_bytes_write >= block_size_bytes) {
    num_bytes_write -= block_size_bytes;
    if (num_bytes_write == 0 && C_partial->data.last) {
      C_partial = nullptr;
      return Cache_Write{.type = Cache_Write::write_last, .num_elements = block_size};
    }
    return Cache_Write{.type = Cache_Write::write, .num_elements = block_size};
  }
  if (num_bytes_write == 0 || !C_partial->data.last) { return Cache_Write{}; }
  const auto cache_write = Cache_Write{.type = Cache_Write::write_last,
				       .num_elements =
                                       static_cast<unsigned>(num_bytes_write / element_size)};
  num_bytes_write = 0;
  C_partial = nullptr;
  return cache_write;
}

bool Fiber_Source::valid() const {
  return index != UINT_MAX;
}

bool Fiber_Source::merge_tree_src() const {
  return task_idx != UINT_MAX;
}

bool Dynamic_Tree_Node::empty() const {
  return data.finished() && !output.valid();
}

bool Task_Allocator::empty() const {
  return num_B_rows == 0 && C_partial_fibers.empty() && allocated_sources.empty();
}

void Task_Allocator::reset() {
  num_B_rows = 0;
  C_partial_fibers.clear();
  trees_allocated.assign(trees_allocated.size(), false);
  allocated_sources.clear();
  output = Task_Output{};
}

bool Task_Allocator::all_rows_allocated() const {
  return num_B_rows == 0 && C_partial_fibers.empty();
}
bool Task_Allocator::last_merge() const {
  return allocated_sources.size() == 2 && all_rows_allocated(); 
}

bool Task_Tree::empty() const {
  return num_C_partials_level.empty();
}

void Task_Tree::reset() {
  tree_level = 0;
  B_rows_first_level = 0;
  B_rows_second_level = 0;
  C_row_idx = UINT_MAX;
  C_row_ptr = UINT_MAX;
  num_C_partials_level.clear();
  C_partial_fibers.clear();
}

Merge_Tree_Manager::Merge_Tree_Manager(const toml::value& parsed_config, Matrix_Data& matrix_data_)
  : matrix_data{matrix_data_}
  , A_row_ptr_fetcher(matrix_data.preproc_A_row_ptr)
  , A_row_idx_fetcher(matrix_data.preproc_A_row_idx)
  , C_row_ptr_fetcher(matrix_data.preproc_C_row_ptr)
  , A_values_fetcher(matrix_data.preproc_A_values)
{
  get_config_params(parsed_config);
}

void Merge_Tree_Manager::reset() {
  mem_read_port.reset();
  prefetch_port.reset();
  for (auto& port : cache_read_ports) { port.reset(); }
  cache_write_port.reset();
  for (auto& port : mem_write_ports) { port.reset(); }

  A_row_ptr_fetcher.reset();
  A_row_ptr_fetcher.base_addr = matrix_data.preproc_A_row_ptr_addr;
  A_row_idx_fetcher.reset();
  A_row_idx_fetcher.base_addr = matrix_data.preproc_A_row_idx_addr;
  C_row_ptr_fetcher.reset();
  C_row_ptr_fetcher.base_addr = matrix_data.C_row_ptr_addr;
  A_values_fetcher.reset();
  A_values_fetcher.base_addr = matrix_data.preproc_A_values_addr;
  read_arbiter = UINT_MAX;
  prefetched_B_rows.clear();

  for (auto& tree : merge_trees) { tree.reset(); }
  std::ranges::fill(dyn_nodes, Dynamic_Tree_Node{});
  std::ranges::fill(C_partial_fibers, C_Partial_Fiber{});
  task_allocator.reset();
  task_tree.reset();
  C_partial_write_idx = UINT_MAX;
  C_partial_head_ptr = nullptr;
  write_arbiter = UINT64_MAX;
  
  num_mults = 0;
  num_block_mults = 0;
  merge_tree_num_merges = 0;
  merge_tree_num_adds = 0;
  dyn_num_merges = 0;
  dyn_num_adds = 0;
  num_idle_cycles = 0;
  C_writes = 0;
  num_C_partial_rows = 0;
  num_C_partial_elements = 0;
  prefetch_stalls = 0;
  A_data_stalls = 0;
  max_write_bytes = 0;
}

void Merge_Tree_Manager::update() {
  write_C_data();
  write_C_partial_data();
  update_dynamic_nodes();
  update_merge_trees();
  allocate_task();
  get_new_task();
  send_A_data_request();
  send_cache_read_requests();
  
  mem_read_port.transfer();
  for (auto& p: cache_read_ports) { p.transfer(); }
  for (auto& p: mem_write_ports) { p.transfer(); }
  cache_write_port.transfer();
}

void Merge_Tree_Manager::apply() {
  receive_A_data();
  receive_prefetch_data();
  receive_cache_data();
}

bool Merge_Tree_Manager::finished() {
  if (num_mults != matrix_data.num_mults) { return false; }
  assert(A_row_idx_fetcher.finished());
  assert(C_row_ptr_fetcher.finished());
  assert(A_values_fetcher.finished());
  assert(prefetched_B_rows.empty());
  for (const auto& tree : merge_trees) {
    if (!tree.inactive()) { return false; }
  }
  for (const auto& node : dyn_nodes) {
    if (!node.empty()) { return false; }
  }
  return true;
}

Merge_Tree_Manager::Mem_Port* Merge_Tree_Manager::get_mem_read_port() {
  return &mem_read_port;
}

Merge_Tree_Manager::Prefetch_Port* Merge_Tree_Manager::get_prefetch_port() {
  return &prefetch_port;
}

Merge_Tree_Manager::Cache_Read_Port* Merge_Tree_Manager::get_cache_read_port(std::size_t id) {
  if (id >= cache_read_ports.size()) return nullptr;
  return &cache_read_ports[id];
}

Merge_Tree_Manager::Cache_Write_Port* Merge_Tree_Manager::get_cache_write_port() {
  return &cache_write_port;
}

Merge_Tree_Manager::Mem_Port* Merge_Tree_Manager::get_mem_write_port(std::size_t id) {
  if (id >= mem_write_ports.size()) return nullptr;
  return &mem_write_ports[id];
}

std::size_t Merge_Tree_Manager::num_mem_ports() const {
  return mem_write_ports.size();
}

std::size_t Merge_Tree_Manager::num_cache_read_ports() const {
  return cache_read_ports.size();
}

void Merge_Tree_Manager::get_config_params(const toml::value& parsed_config) {
  const auto A_row_ptr_buffer_size = toml::find_or(parsed_config, "merge_tree_manager",
                                                   "A_row_ptr_buffer_size", 16u);
  A_row_ptr_fetcher.buffer_size = A_row_ptr_buffer_size;
  A_row_idx_fetcher.buffer_size = A_row_ptr_buffer_size;
  C_row_ptr_fetcher.buffer_size = A_row_ptr_buffer_size;
  max_prefetched_rows = toml::find_or(parsed_config, "merge_tree_manager",
                                      "max_prefetched_rows", 1024u);
  A_values_fetcher.buffer_size = max_prefetched_rows;
  const auto num_merge_trees = toml::find<unsigned>(parsed_config,
                                                    "merge_tree_manager",
                                                    "num_merge_trees");
  merge_tree_size = toml::find<unsigned>(parsed_config, "merge_tree_manager",
                                         "merge_tree_size");
  merge_trees = std::vector<Merge_Tree>(num_merge_trees, Merge_Tree{*this});
  dyn_nodes = std::vector<Dynamic_Tree_Node>(num_merge_trees - 1);
  merge_tree_merger_width = toml::find<unsigned>(parsed_config,
                                                 "merge_tree_manager",
                                                 "merge_tree_merger_width");
  merge_tree_merger_num_adds = toml::find_or(parsed_config,
                                             "merge_tree_manager",
                                             "merge_tree_merger_num_adds",
                                             merge_tree_merger_width - 1);
  num_final_mergers = toml::find<unsigned>(parsed_config, "merge_tree_manager",
                                           "num_final_mergers");
  dyn_merger_width = toml::find<unsigned>(parsed_config, "merge_tree_manager",
                                          "final_merger_width");
  dyn_merger_num_adds = toml::find_or(parsed_config, "merge_tree_manager",
                                      "merge_tree_merger_num_adds",
                                      dyn_merger_width - 1);
  max_rows_merge = num_merge_trees * merge_tree_size;
  const auto task_tree_max_level = 32u / log2_ceil(max_rows_merge);
  const auto max_partial_rows = task_tree_max_level * max_rows_merge;
  C_partial_fibers = std::vector<C_Partial_Fiber>(max_partial_rows);
  task_allocator.C_partial_fibers.reserve(max_rows_merge);
  task_allocator.trees_allocated.assign(num_merge_trees, false);
  task_allocator.allocated_sources.reserve(num_merge_trees);
  task_tree.num_C_partials_level.reserve(task_tree_max_level);
  task_tree.C_partial_fibers.reserve(max_partial_rows);
  const auto num_mem_ports = toml::find<std::size_t>(parsed_config,
                                                     "merge_tree_manager",
                                                     "num_mem_ports");
  mem_write_ports = std::vector<Mem_Port>(num_mem_ports);
  cache_read_ports = std::vector<Cache_Read_Port>(num_merge_trees);
  input_buffer_size = toml::find_or(parsed_config, "merge_tree_manager",
                                    "input_buffer_size", 2 * block_size);
  output_buffer_size = toml::find_or(parsed_config, "merge_tree_manager",
                                    "output_buffer_size", 2 * dyn_merger_width);
}

void Merge_Tree_Manager::send_A_data_request() {
  if (mem_read_port.has_msg_send()) { return; }
  Mem_Request request{};
  for (unsigned i = 0; i < 4; ++i) {
    read_arbiter = inc_mod(read_arbiter, 4U);
    switch (read_arbiter) {
    case 0:
      request.address = A_row_ptr_fetcher.get_fetch_address();
      break;
    case 1:
      request.address = A_row_idx_fetcher.get_fetch_address();
      break;
    case 2:
      request.address = C_row_ptr_fetcher.get_fetch_address();
      break;
    case 3:
      request.address = A_values_fetcher.get_fetch_address();
      break;
    }
    if (request.valid()) {
      request.id = read_arbiter;
      mem_read_port.add_msg_send(request);
      ++preproc_A_reads;
      return;
    }
  }
}

void Merge_Tree_Manager::send_cache_read_requests() {
  for (unsigned i = 0; i < cache_read_ports.size(); ++i) {
    if (cache_read_ports[i].has_msg_send()) { continue; }
    const auto request = merge_trees[i].get_request();
    if (request.valid()) {
      cache_read_ports[i].add_msg_send(request);
    }
  }
}

void Merge_Tree_Manager::write_C_data() {
  const auto size = merge_trees.size() + dyn_nodes.size();
  for (auto& port : mem_write_ports) {
    if (port.has_msg_send()) { continue; }
    for (unsigned i = 0; i != size; ++i) {
      write_arbiter = inc_mod(write_arbiter, size);
      Address address {};
      if (write_arbiter < merge_trees.size()) {
        address = merge_trees[write_arbiter].get_C_write_address();
      } else {
        auto& output = dyn_nodes[write_arbiter - merge_trees.size()].output;
        address = output.get_C_write_address();
      }
      if (address == invalid_address) { continue; }
      ++C_writes;
      port.add_msg_send(Mem_Request{.address = address, .is_write = true});
      break;
    }
  }
}

void Merge_Tree_Manager::write_C_partial_data() {
  if (C_partial_write_idx == UINT_MAX || cache_write_port.has_msg_send()) { return; }
  Cache_Write cache_write {};
  if (C_partial_write_idx < merge_trees.size()) {
    cache_write = merge_trees[C_partial_write_idx].get_C_partial_write();
  } else {
    auto& output = dyn_nodes[C_partial_write_idx - merge_trees.size()].output;
    cache_write = output.get_C_partial_write();
  }
  if (!cache_write.valid()) { return; }
  if (cache_write.type == Cache_Write::write_last) {
    C_partial_write_idx = UINT_MAX;
    ++num_C_partial_rows;
  }
  cache_write_port.add_msg_send(cache_write);
}

void Merge_Tree_Manager::update_dynamic_nodes() {
  std::vector<unsigned> possible_merges;
  for (unsigned i = 0; i != dyn_nodes.size(); ++i) {
    if (dyn_nodes[i].data.size() > output_buffer_size - dyn_merger_width
	|| dyn_nodes[i].output.num_bytes_write >
        (output_buffer_size - dyn_merger_width) * element_size)
    {
      continue;
    }
    if (!dyn_nodes[i].src1.valid() && !dyn_nodes[i].src2.valid()) { continue; }
    if (dyn_nodes[i].src1.valid() && !fiber_source_ready(dyn_nodes[i].src1)) { continue; }
    if (dyn_nodes[i].src2.valid() && !fiber_source_ready(dyn_nodes[i].src2)) { continue; }    
    possible_merges.push_back(i);
  }
  for (unsigned i = 0, num_merges = 0; i != possible_merges.size(); ++i) {
    auto& node = dyn_nodes[possible_merges[i]];
    auto& node_dest = (node.output.C_partial) ? node.output.C_partial->data : node.data;
    unsigned num_elements_out = 0;
    if (node.src1.valid() && node.src2.valid()) {
      if (num_merges == num_final_mergers) { continue; }
      num_elements_out = do_merge_add(node_dest, fiber_source_node(node.src1),
                                      fiber_source_node(node.src2), false);
      if (fiber_source_node(node.src1).finished()) {
        fiber_source_reset(node.src1);
      }
      if (fiber_source_node(node.src2).finished()) {
        fiber_source_reset(node.src2);
      }
      ++num_merges;
    } else if (node.src1.valid()) {
      num_elements_out = fiber_buffer_transfer(fiber_source_node(node.src1), node_dest,
                                               dyn_merger_width);
      if (fiber_source_node(node.src1).finished()) {
        fiber_source_reset(node.src1);
      }
    } else {
      assert(node.src2.valid());
      num_elements_out = fiber_buffer_transfer(fiber_source_node(node.src2), node_dest,
                                               dyn_merger_width);
      if (fiber_source_node(node.src2).finished()) {
        fiber_source_reset(node.src2);
      }
    }
    node.data.last = node_dest.last;
    if (node.output.valid()) {
      write_C_output(node.output, node.data, num_elements_out);
    }
  //   if (node.src1.valid()) {
  //     auto& node_src1 = fiber_source_node(node.src1);
  //     if (node.src2.valid() && num_merges != num_final_mergers) {
  //       auto& node_src2 = fiber_source_node(node.src2);
  //       num_elements_out = do_merge_add(node_dest, node_src1, node_src2, false);
  //       if (node_src2.finished()) { node.src2 = Fiber_Source{}; }
  //       ++num_merges;
  //     } else {
  //       num_elements_out = fiber_buffer_transfer(node_src1, node_dest, dyn_merger_width);
  //     }
  //     if (node_src1.finished()) { node.src1 = Fiber_Source{}; }
  //   } else {
  //     auto& node_src2 = fiber_source_node(node.src2);
  //     num_elements_out = fiber_buffer_transfer(node_src2, node_dest, dyn_merger_width);
  //     if (node_src2.finished()) { node.src2 = Fiber_Source{}; }
  //   }
  //   node.data.last = node_dest.last;
  //   if (node.output.valid()) {
  //     node.output.num_bytes_write += num_elements_out * element_size;
  //     max_write_bytes = std::max(max_write_bytes, node.output.num_bytes_write);
  //     if (node.output.write_address != invalid_address) {
  //       write_C_output(node.output, node.data);
  //     }
  //   }
  }
}

bool Merge_Tree_Manager::fiber_source_ready(const Fiber_Source& src) const {
  assert(src.valid());
  if (src.merge_tree_src()) {
    auto& root_level = merge_trees[src.index].levels[0]; 
    if (root_level.task != src.task_idx) {
      return false;
    }
    return root_level.nodes[0].ready_to_merge(dyn_merger_width);
  }
  // dynamic node src
  return dyn_nodes[src.index].data.ready_to_merge(dyn_merger_width);
}

Fiber_Buffer& Merge_Tree_Manager::fiber_source_node(const Fiber_Source& src) {
  assert(src.valid());
  return (src.merge_tree_src()) ? merge_trees[src.index].levels[0].nodes[0]
    : dyn_nodes[src.index].data; 
}

void Merge_Tree_Manager::fiber_source_reset(Fiber_Source& src) {
  if (src.merge_tree_src()) {
    merge_trees[src.index].levels[0].task = UINT_MAX;
    merge_trees[src.index].levels[0].num_active_nodes = 0;
  }
  src = Fiber_Source{};
}

void Merge_Tree_Manager::write_C_output(Task_Output& output, Fiber_Buffer& node,
                                        unsigned num_elements_out)
{
  output.num_bytes_write += num_elements_out * element_size;
  max_write_bytes = std::max(max_write_bytes, output.num_bytes_write);
  if (output.write_address == invalid_address) { return; }
  assert(node.size() == num_elements_out);
  while (!node.empty()) {
    if (matrix_data.compute_result) {
      matrix_data.C.col_idx[output.C_row_ptr] = node.col_idx.front();
      matrix_data.C.values[output.C_row_ptr] = node.values.front();
      node.values.pop_front();
    }
    node.col_idx.pop_front();
    ++output.C_row_ptr;
    ++matrix_data.C.nnz;
  }
  if (node.finished()) {
    matrix_data.C.row_end[output.C_row_idx] = output.C_row_ptr;
    output.C_row_idx = UINT_MAX;
    output.C_row_ptr = UINT_MAX;
  }
}

void Merge_Tree_Manager::update_merge_trees() {
  for (auto& tree : merge_trees) { tree.update(); }
}

void Merge_Tree_Manager::allocate_task() {
  if (task_allocator.empty()) { return; }
  if (!task_allocator.all_rows_allocated()) {
    // find available tree to allocate
    for (unsigned i = 0; i != merge_trees.size(); ++i) {
      if (add_task_merge_tree(i)) { return; }
    }
  }
  // find two sources ready to merge
  if (task_allocator.allocated_sources.size() < 2) { return; }
  unsigned idx_merge = UINT_MAX;
  for (unsigned i = 1; i != task_allocator.allocated_sources.size(); ++i) {
    auto& prev_src = task_allocator.allocated_sources[i - 1];
    auto& cur_src = task_allocator.allocated_sources[i];
    if (prev_src.second == cur_src.second) {
      idx_merge = i;
      break;
    }
  }
  if (idx_merge == UINT_MAX) {
    if (task_allocator.all_rows_allocated()) {
      idx_merge = static_cast<unsigned>(task_allocator.allocated_sources.size()) - 1;
    } else {
      return;
    }
  }
  // find available node to allocate
  for (unsigned i = 0; i != dyn_nodes.size(); ++i) {
    if (dyn_nodes[i].empty()) {
       add_task_dyn_node(i, idx_merge);
       return;
    }
  }
}

bool Merge_Tree_Manager::task_allocator_single_subtask() const {
  return task_allocator.allocated_sources.empty()
    && task_allocator.num_B_rows + task_allocator.C_partial_fibers.size() <= merge_tree_size;
}

bool Merge_Tree_Manager::add_task_merge_tree(unsigned tree_idx) {
  if (task_allocator.trees_allocated[tree_idx]) { return false; }
  const auto B_rows_to_allocate = std::min(merge_tree_size, task_allocator.num_B_rows);
  if (A_values_fetcher.num_elements < B_rows_to_allocate
      || prefetched_B_rows.size() < B_rows_to_allocate)
  {
    ++A_data_stalls;
    return false;
  }
  auto& tree = merge_trees[tree_idx];
  if (tree.num_active_inputs > 0) { return false; }
  auto& output = tree.outputs[tree.input_task];
  if (output.valid()) { return false; }
  if (task_allocator_single_subtask()) {
    if (task_allocator.output.C_partial) {
      if (C_partial_write_idx != UINT_MAX || C_partial_head_ptr != nullptr) {
        ++C_partial_stalls;
        return false;
      }
      C_partial_write_idx = tree_idx;
      C_partial_head_ptr = task_allocator.output.C_partial;
    }
    output = task_allocator.output;
    task_allocator.output = Task_Output{};
  } else {
    task_allocator.allocated_sources.emplace_back(Fiber_Source{ .index = tree_idx,
                                                                .task_idx = tree.input_task}, 0);
    task_allocator.trees_allocated[tree_idx] = true;
  }
  // init B rows in inputs
  while (tree.num_active_inputs < B_rows_to_allocate) {
    tree.inputs[tree.num_active_inputs] =
      Input_Fiber{ .A_value = A_values_fetcher.front(),
                   .B_row_ptr = prefetched_B_rows.front().B_row_ptr,
                   .head_ptr = prefetched_B_rows.front().row_head_ptr };
    A_values_fetcher.pop();
    prefetched_B_rows.pop_front();
    ++tree.num_active_inputs;
  }
  task_allocator.num_B_rows -= B_rows_to_allocate;
  // init C partial rows in inputs
  while (tree.num_active_inputs < tree.inputs.size()
         && !task_allocator.C_partial_fibers.empty())
  {
    tree.inputs[tree.num_active_inputs] =
      Input_Fiber{ .C_partial_fiber =
                   task_allocator.C_partial_fibers.back() };
    task_allocator.C_partial_fibers.pop_back();
    ++tree.num_active_inputs;
  }
  return true;
}

void Merge_Tree_Manager::add_task_dyn_node(unsigned node_idx,
                                           unsigned idx_merge)
{
  auto& prev_src = task_allocator.allocated_sources[idx_merge - 1];
  auto& cur_src = task_allocator.allocated_sources[idx_merge];
  if (task_allocator.last_merge()) {
    if (task_allocator.output.C_partial) {
      if (C_partial_write_idx != UINT_MAX || C_partial_head_ptr != nullptr) {
        ++C_partial_stalls;
        return;
      }
      C_partial_write_idx = node_idx + static_cast<unsigned>(merge_trees.size());
      C_partial_head_ptr = task_allocator.output.C_partial;
    }
    dyn_nodes[node_idx].src1 = prev_src.first;
    dyn_nodes[node_idx].src2 = cur_src.first;
    dyn_nodes[node_idx].data.last = false;
    dyn_nodes[node_idx].output = task_allocator.output;
    task_allocator.reset();
  } else {
    dyn_nodes[node_idx].src1 = prev_src.first;
    dyn_nodes[node_idx].src2 = cur_src.first;
    dyn_nodes[node_idx].data.last = false;
    prev_src.first = Fiber_Source{ .index = node_idx,
                                   .task_idx = UINT_MAX };
    ++prev_src.second;
    auto it = task_allocator.allocated_sources.begin()
      + static_cast<std::ptrdiff_t>(idx_merge);
    task_allocator.allocated_sources.erase(it);
  }
}

void Merge_Tree_Manager::get_new_task() {
  if (!task_allocator.empty()) { return; }
  if (task_tree.empty()) {
    if (A_row_ptr_fetcher.num_elements < 2
        || A_row_idx_fetcher.num_elements == 0
        || C_row_ptr_fetcher.num_elements == 0)
    {
      return;
    }
    const unsigned A_row_idx = A_row_idx_fetcher.front();
    const unsigned C_row_ptr = C_row_ptr_fetcher.front();
    const unsigned num_rows_merge = A_row_ptr_fetcher.at(1)
      - A_row_ptr_fetcher.front();
    A_row_ptr_fetcher.pop();
    A_row_idx_fetcher.pop();
    C_row_ptr_fetcher.pop();
    if (num_rows_merge <= max_rows_merge) {
      task_allocator.output.C_row_idx = A_row_idx;
      task_allocator.output.C_row_ptr = C_row_ptr;
      task_allocator.output.write_address =
        matrix_data.C_elements_addr + C_row_ptr * element_size;
      task_allocator.num_B_rows = num_rows_merge;
      return;
    } else {
      init_task_tree(num_rows_merge, A_row_idx, C_row_ptr);
    }
  }
  assert(!task_tree.empty());
  // create task from task tree
  const auto last_level = task_tree.num_C_partials_level.size() - 1;
  if (task_tree.tree_level == 0) {
    assert(task_tree.B_rows_first_level > 0);
    const auto C_partial_ptr = get_C_partial_fiber();
    if (C_partial_ptr == nullptr) { return; }
    const auto B_rows_merge = std::min(task_tree.B_rows_first_level,
                                       max_rows_merge);
    task_tree.B_rows_first_level -= B_rows_merge;
    assert(task_tree.C_partial_fibers[task_tree.num_C_partials_level[0]]
           == nullptr);
    task_tree.C_partial_fibers[task_tree.num_C_partials_level[0]] =
      C_partial_ptr;
    task_allocator.output.C_partial = C_partial_ptr;
    task_allocator.num_B_rows = B_rows_merge;
    ++task_tree.num_C_partials_level[0];
    if (task_tree.num_C_partials_level[0] == max_rows_merge
        || task_tree.B_rows_first_level == 0)
    {
      task_tree.tree_level = 1;
    }
    return;
  }
  if (task_tree.tree_level == 1) {
    if (task_tree.tree_level == last_level) {
      assert(task_tree.B_rows_second_level + task_tree.num_C_partials_level[0]
             == max_rows_merge);
      task_allocator.output.C_row_idx = task_tree.C_row_idx;
      task_allocator.output.C_row_ptr = task_tree.C_row_ptr;
      task_allocator.output.write_address =
        matrix_data.C_elements_addr + task_tree.C_row_ptr * element_size;
      task_allocator.num_B_rows = task_tree.B_rows_second_level;
      for (unsigned i = 0; i < task_tree.num_C_partials_level[0]; ++i) {
        task_allocator.C_partial_fibers.push_back(
          task_tree.C_partial_fibers[i]);
      }
      task_tree.reset();
      return;
    }
    const auto C_partial_ptr = get_C_partial_fiber();
    if (C_partial_ptr == nullptr) { return; }
    const auto B_rows_merge = max_rows_merge
      - task_tree.num_C_partials_level[0];
    assert(task_tree.C_partial_fibers[max_rows_merge
                                      + task_tree.num_C_partials_level[1]]
           == nullptr);
    task_tree.C_partial_fibers[max_rows_merge
                               + task_tree.num_C_partials_level[1]] =
      C_partial_ptr;
    task_allocator.output.C_partial = C_partial_ptr;
    task_allocator.num_B_rows = B_rows_merge;
    for (unsigned i = 0; i < task_tree.num_C_partials_level[0]; ++i) {
      task_allocator.C_partial_fibers.push_back(task_tree.C_partial_fibers[i]);
      task_tree.C_partial_fibers[i] = nullptr;
    }
    task_tree.num_C_partials_level[0] = 0;
    ++task_tree.num_C_partials_level[1];
    if (task_tree.num_C_partials_level[1] == max_rows_merge) {
      ++task_tree.tree_level;
    } else if (task_tree.B_rows_first_level > 0) {
      task_tree.tree_level = 0;
    }
    return;
  }
  if (task_tree.tree_level < last_level) {
    assert(task_tree.num_C_partials_level[task_tree.tree_level - 1]
           == max_rows_merge);
    const auto C_partial_ptr = get_C_partial_fiber();
    if (C_partial_ptr == nullptr) { return; }
    const auto idx = max_rows_merge * task_tree.tree_level
      + task_tree.num_C_partials_level[task_tree.tree_level];
    assert(task_tree.C_partial_fibers[idx] == nullptr);
    task_tree.C_partial_fibers[idx] = C_partial_ptr;
    task_allocator.output.C_partial = C_partial_ptr;
    for (unsigned i = 0; i < max_rows_merge; ++i) {
      auto& C_partial =
        task_tree.C_partial_fibers[(task_tree.tree_level-1)
                                   * max_rows_merge + i];
      task_allocator.C_partial_fibers.push_back(C_partial);
      C_partial = nullptr;
    }
    task_tree.num_C_partials_level[task_tree.tree_level-1] = 0;
    ++task_tree.num_C_partials_level[task_tree.tree_level];
    if (task_tree.num_C_partials_level[task_tree.tree_level] == max_rows_merge) {
      ++task_tree.tree_level;
    } else if (task_tree.B_rows_first_level > 0) {
      task_tree.tree_level = 0;
    } else {
      task_tree.tree_level = 1;
    }
    return;
  }
  // last level
  assert(task_tree.num_C_partials_level[task_tree.tree_level - 1]
         == max_rows_merge);
  task_allocator.output.C_row_idx = task_tree.C_row_idx;
  task_allocator.output.C_row_ptr = task_tree.C_row_ptr;
  task_allocator.output.write_address =
    matrix_data.C_elements_addr + task_tree.C_row_ptr * element_size;
  for (unsigned i = 0; i < max_rows_merge; ++i) {
    auto& C_partial =
      task_tree.C_partial_fibers[(task_tree.tree_level - 1)
                                 * max_rows_merge + i];
    task_allocator.C_partial_fibers.push_back(C_partial);
    C_partial = nullptr;
  }
  task_tree.reset();
  return;
}

void Merge_Tree_Manager::init_task_tree(unsigned num_rows, unsigned C_row_idx_,
                                        unsigned C_row_ptr_)
{
  const auto second_level_num_rows = nearest_pow_floor(num_rows, max_rows_merge);
  task_tree.B_rows_first_level =
    div_ceil((num_rows - second_level_num_rows) * max_rows_merge,
                                          max_rows_merge - 1);
  task_tree.B_rows_second_level = num_rows - task_tree.B_rows_first_level;
  const auto num_levels = log_ceil(num_rows, max_rows_merge);
  task_tree.num_C_partials_level.assign(num_levels, 0);
  task_tree.C_partial_fibers.assign(num_levels * max_rows_merge, nullptr);
  task_tree.C_row_idx = C_row_idx_;
  task_tree.C_row_ptr = C_row_ptr_;
}

C_Partial_Fiber* Merge_Tree_Manager::get_C_partial_fiber() {
  for (auto& p : C_partial_fibers) {
    if (p.finished()) {
      p.data.last = false;
      return &p;
    }
  }
  return nullptr;
}

void Merge_Tree_Manager::receive_A_data() {
  if (!mem_read_port.msg_received_valid()) { return; }
  const auto mem_read_resp = mem_read_port.get_msg_received();
  assert(mem_read_resp.id < 4);
  switch (mem_read_resp.id) {
  case 0:
    A_row_ptr_fetcher.receive_data(mem_read_resp.address);
    break;
  case 1:
    A_row_idx_fetcher.receive_data(mem_read_resp.address);
    break;
  case 2:
    C_row_ptr_fetcher.receive_data(mem_read_resp.address);
    break;
  case 3:
    A_values_fetcher.receive_data(mem_read_resp.address);
    break;
  }
  mem_read_port.clear_msg_received();
}

void Merge_Tree_Manager::receive_prefetch_data() {
  if (prefetch_port.msg_received_valid()) {
    const auto& prefetch_resp = prefetch_port.get_msg_received();
    if (prefetched_B_rows.size() + prefetch_resp.size()
        <= max_prefetched_rows)
    {
      prefetched_B_rows.insert(prefetched_B_rows.end(), prefetch_resp.begin(),
                               prefetch_resp.end());
      prefetch_port.clear_msg_received();
    }
  }
}

void Merge_Tree_Manager::receive_cache_data() {
  for (unsigned i = 0; i < cache_read_ports.size(); ++i) {
    if (!cache_read_ports[i].msg_received_valid()) { continue; }
    merge_trees[i].receive_response(cache_read_ports[i].get_msg_received());
    cache_read_ports[i].clear_msg_received();
  }
  if (!cache_write_port.msg_received_valid()) { return; }
  assert(C_partial_head_ptr && C_partial_head_ptr->head_ptr == UINT_MAX);
  C_partial_head_ptr->head_ptr = cache_write_port.get_msg_received();
  C_partial_head_ptr = nullptr;
  cache_write_port.clear_msg_received();
}
unsigned Merge_Tree_Manager::do_merge_add(Fiber_Buffer& dest,
                                          Fiber_Buffer& src1,
                                          Fiber_Buffer& src2,
                                          bool is_merge_tree)
{
  assert(!src1.empty() && !src2.empty());
  const auto [merge_width, max_num_adds] = (is_merge_tree) ?
    std::tie(merge_tree_merger_width, merge_tree_merger_num_adds)
    : std::tie(dyn_merger_width, dyn_merger_num_adds);
  unsigned num_elements_output {};
  unsigned num_adds {};
  while (num_elements_output < merge_width && num_adds < max_num_adds) {
    if (src1.empty()) {
      num_elements_output +=
        fiber_buffer_transfer(src2, dest, merge_width - num_elements_output);
      break;
    } else if (src2.empty()) {
      num_elements_output +=
        fiber_buffer_transfer(src1, dest, merge_width - num_elements_output);
      break;
    }
    const auto res  = src1.col_idx.front() <=> src2.col_idx.front();
    if (res < 0) {
      dest.col_idx.push_back(src1.col_idx.front());
      src1.col_idx.pop_front();
      if (matrix_data.compute_result) {
        dest.values.push_back(src1.values.front());
        src1.values.pop_front();
      }
    } else if (res > 0) {
      dest.col_idx.push_back(src2.col_idx.front());
      src2.col_idx.pop_front();
      if (matrix_data.compute_result) {
        dest.values.push_back(src2.values.front());
        src2.values.pop_front();
      }
    } else {
      dest.col_idx.push_back(src1.col_idx.front());
      src1.col_idx.pop_front();
      src2.col_idx.pop_front();
      if (!src1.values.empty()) {
        const auto add_value = src1.values.front() + src2.values.front();
        dest.values.push_back(add_value);
        src1.values.pop_front();
        src2.values.pop_front();
      }
      ++num_adds;
    }
    ++num_elements_output;
  }
  if (src1.finished() && src2.finished()) {
    dest.last = true;
  }
  if (is_merge_tree) {
    ++merge_tree_num_merges;
    merge_tree_num_adds += num_adds;
  } else {
    ++dyn_num_merges;
    dyn_num_adds += num_adds;
  }
  return num_elements_output;
}

unsigned fiber_buffer_transfer(Fiber_Buffer& src, Fiber_Buffer& dest, std::size_t num_elements) {
  const auto n =
    static_cast<std::ptrdiff_t>(std::min(num_elements, src.size()));
  if (n == 0) { return 0; }
  dest.col_idx.insert(dest.col_idx.end(), src.col_idx.begin(),
                      src.col_idx.begin() + n);
  src.col_idx.erase(src.col_idx.begin(), src.col_idx.begin() + n);
  if (!src.values.empty()) {
    dest.values.insert(dest.values.end(), src.values.begin(),
                       src.values.begin() + n);
    src.values.erase(src.values.begin(), src.values.begin() + n);
  }
  if (src.finished()) {
    dest.last = true;
  }
  return static_cast<unsigned>(n);
}

} // namespace mergeforest

} // namespace mergeforest_sim
