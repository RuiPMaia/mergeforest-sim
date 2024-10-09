#include <mergeforest-sim/gamma/PE_manager.hpp>
#include <mergeforest-sim/math_utils.hpp>

namespace mergeforest_sim {

namespace gamma {

bool C_Partial_Fiber::empty() const {
  return begin == invalid_address;
}

bool C_Partial_Fiber::is_finished() const {
  return finished && col_idx.empty();
}

bool Input_Fiber::finished() const {
  if (C_partial_fiber) {
    return C_partial_fiber->is_finished();
  }
  return B_row_ptr == B_row_end;
}

bool Task::valid() const { return !inputs.empty(); }

PE::PE(Matrix_Data& matrix_data_) : matrix_data{ matrix_data_ } {
  reset();
}

void PE::reset() {
  cur_task = Task{};
  next_task = Task{};
  cur_task_finished = false;
  C_col_idx = UINT_MAX;
  C_value = 0.0;
  input_buffers = std::vector<Input_Buffer>(PE::radix);
  read_arbiter = UINT64_MAX;
  write_address = invalid_address;
  num_bytes_write = 0;
}

Mem_Request PE::get_cache_request() {
  if (!cur_task.valid()) return Mem_Request{};
  for (std::size_t i = 0; i < input_buffers.size(); ++i) {
    read_arbiter = inc_mod(read_arbiter, input_buffers.size());
    // check if fiber to fetch is from the current task or the next
    bool fetching_next_task {false};
    Input_Fiber* in_fiber {};
    if (cur_task.inputs.size() > read_arbiter && !cur_task.inputs[read_arbiter].finished()) {
      in_fiber = &cur_task.inputs[read_arbiter];
    } else if (next_task.inputs.size() > read_arbiter && !next_task.inputs[read_arbiter].finished()) {
      in_fiber = &next_task.inputs[read_arbiter];
      fetching_next_task = true;
    } else {
      continue;
    }
    auto& buffer = input_buffers[read_arbiter];
    unsigned num_elements_fetch {};
    if (in_fiber->C_partial_fiber) {
      const unsigned C_num_elements = static_cast<unsigned>(in_fiber->C_partial_fiber->col_idx.size());
      if (C_num_elements == 0) continue;
      if (C_num_elements >= block_size) {
	num_elements_fetch = block_size;
      } else if (in_fiber->C_partial_fiber->finished) {
	num_elements_fetch = C_num_elements;
      } else {
	continue;
      }
    } else {
      num_elements_fetch = std::min(in_fiber->B_row_end - in_fiber->B_row_ptr, block_size - in_fiber->B_row_ptr % block_size);
    }
    if (buffer.col_idx.size() + num_elements_fetch > Input_Buffer::buffer_size) continue;
    Mem_Request req{.id = static_cast<unsigned>(read_arbiter), .is_write = false};
    // put data in buffers and set request address
    if (in_fiber->C_partial_fiber) {
      if (in_fiber->C_partial_fiber->begin == in_fiber->C_partial_fiber->end) continue;
      req.address = in_fiber->C_partial_fiber->begin;
      in_fiber->C_partial_fiber->begin += num_elements_fetch * element_size;
      assert(in_fiber->C_partial_fiber->end >= in_fiber->C_partial_fiber->begin);
      // transfer data from C_partial fiber to buffer
      auto& C_partial_col_idx = in_fiber->C_partial_fiber->col_idx;
      buffer.col_idx.insert(buffer.col_idx.end(), C_partial_col_idx.begin(), C_partial_col_idx.begin() + num_elements_fetch); 
      C_partial_col_idx.erase(C_partial_col_idx.begin(), C_partial_col_idx.begin() + num_elements_fetch);
      if (matrix_data.compute_result) {
	auto& C_partial_values = in_fiber->C_partial_fiber->values;
	buffer.values.insert(buffer.values.end(), C_partial_values.begin(), C_partial_values.begin() + num_elements_fetch); 
	C_partial_values.erase(C_partial_values.begin(), C_partial_values.begin() + num_elements_fetch);
      }
      if (in_fiber->C_partial_fiber->is_finished()) {
	assert(C_Partial_Fiber::num_fibers > 0);
	--C_Partial_Fiber::num_fibers;
	*in_fiber->C_partial_fiber = C_Partial_Fiber{};
	in_fiber->C_partial_fiber = nullptr;
      }
    } else {
      req.address = matrix_data.B_elements_addr + in_fiber->B_row_ptr * element_size;
      for (unsigned j = 0; j < num_elements_fetch; ++j) {
	buffer.col_idx.push_back(matrix_data.B->col_idx[in_fiber->B_row_ptr + j]);
	if (matrix_data.compute_result) {
	  buffer.values.push_back(matrix_data.B->values[in_fiber->B_row_ptr + j]);
	}
      }
      in_fiber->B_row_ptr += num_elements_fetch;
      PE::num_mults += num_elements_fetch;
      assert(in_fiber->B_row_end >= in_fiber->B_row_ptr);
    }
    req.address = round_down_multiple(req.address, static_cast<std::size_t>(block_size_bytes));
    buffer.pending_reqs.emplace_back(req.address, num_elements_fetch, false);
    if (!fetching_next_task) {
      buffer.num_elems_fetched_cur_task += num_elements_fetch;
    }
    return req;
  }
  return Mem_Request{};
}

void PE::receive_cache_response(Mem_Response mem_response) {
  if (!mem_response.valid()) return;
  auto& buffer = input_buffers[mem_response.id];
  assert(!buffer.pending_reqs.empty());
  for (auto& req : buffer.pending_reqs) {
    if (std::get<0>(req) == mem_response.address) {
      std::get<2>(req) = true;
      break;
    }
  } 
  while (!buffer.pending_reqs.empty()) {
    if (!std::get<2>(buffer.pending_reqs.front())) break;
    buffer.num_elements_received += std::get<1>(buffer.pending_reqs.front());
    buffer.pending_reqs.pop_front();
  }     
  assert(buffer.num_elements_received <= buffer.col_idx.size());
}

void PE::update() {
  if (!cur_task.valid()) {
    ++PE::idle_cycles;
    return;
  }
  if (cur_task_finished) return;
  if (num_bytes_write + element_size > PE::output_buffer_size * element_size) {
    ++PE::write_stalls;
    return;
  }
  unsigned min_col_idx = {UINT_MAX};
  unsigned min_idx = {UINT_MAX};
  unsigned finished_inputs {0};
  bool stall {false};
  for (unsigned i = 0; i < cur_task.inputs.size(); ++i) {
    if (input_buffers[i].num_elems_fetched_cur_task == 0 && cur_task.inputs[i].finished()) {
      ++finished_inputs;
      continue;
    }
    if (input_buffers[i].num_elements_received == 0) {
      stall = true;
      continue;
    };
    if (input_buffers[i].col_idx.front() < min_col_idx) {
      min_col_idx = input_buffers[i].col_idx.front();
      min_idx = i;
    }
  }
  // finish task
  if (finished_inputs == cur_task.inputs.size()) {
    cur_task_finished = true;
    assert(C_col_idx != UINT_MAX);
    if (cur_task.C_partial_fiber) {
      cur_task.C_partial_fiber->col_idx.push_back(C_col_idx);
      if (matrix_data.compute_result) {
	cur_task.C_partial_fiber->values.push_back(C_value);
      }
      cur_task.C_partial_fiber->finished = true;
      ++PE::num_C_partial_elements;
      ++PE::num_C_partial_rows;
    } else {
      if (matrix_data.compute_result) {
	matrix_data.C.col_idx[cur_task.C_row_ptr] = C_col_idx;
	matrix_data.C.values[cur_task.C_row_ptr] = C_value;
      }
      ++cur_task.C_row_ptr;
      matrix_data.C.row_end[cur_task.C_row_idx] = cur_task.C_row_ptr;
      ++matrix_data.C.nnz;
      ++PE::num_finished_rows;
    }
    num_bytes_write += element_size;
    PE::max_bytes_write = std::max(PE::max_bytes_write, num_bytes_write);
    C_col_idx = UINT_MAX;
    C_value = 0.0;
    return;
  }
  if (stall) {
    ++PE::B_data_stalls;
    return;
  }
  assert(min_idx != UINT_MAX);
  // execute one multiply add
  if (C_col_idx == UINT_MAX) {
    C_col_idx = min_col_idx;
    if (matrix_data.compute_result) {
      C_value = cur_task.inputs[min_idx].A_value * input_buffers[min_idx].values.front();
    }
  } else if (min_col_idx > C_col_idx) {
    if (cur_task.C_partial_fiber) {
      ++PE::num_C_partial_elements;
      cur_task.C_partial_fiber->col_idx.push_back(C_col_idx);
      if (matrix_data.compute_result) {
	cur_task.C_partial_fiber->values.push_back(C_value);
      }
    } else {
      ++matrix_data.C.nnz;
      if (matrix_data.compute_result) {
	matrix_data.C.values[cur_task.C_row_ptr] = C_value;
	matrix_data.C.col_idx[cur_task.C_row_ptr] = C_col_idx;
      }
      ++cur_task.C_row_ptr;
    }
    num_bytes_write += element_size;
    PE::max_bytes_write = std::max(PE::max_bytes_write, num_bytes_write);
    C_col_idx = min_col_idx;
    if (matrix_data.compute_result) {
      C_value = cur_task.inputs[min_idx].A_value * input_buffers[min_idx].values.front();
    }
  } else {
    assert(min_col_idx == C_col_idx);
    ++PE::num_adds;
    if (matrix_data.compute_result) {
      C_value += cur_task.inputs[min_idx].A_value * input_buffers[min_idx].values.front();
    }
  }
  // pop element from input buffer
  --input_buffers[min_idx].num_elements_received;
  --input_buffers[min_idx].num_elems_fetched_cur_task;
  input_buffers[min_idx].col_idx.pop_front();
  if (matrix_data.compute_result) {
    input_buffers[min_idx].values.pop_front();
  }
}

void Task_Tree::reset() {
  tree_level = 0;
  B_rows_first_level = 0;
  B_rows_second_level = 0;
  C_row_ptr = UINT_MAX;
  C_row_idx = UINT_MAX;
  num_C_partials_level.clear();
  C_partial_fibers.clear();
}

void Task_Tree::init(unsigned num_rows, unsigned C_row_idx_, unsigned C_row_ptr_) {
  unsigned second_level_num_rows = nearest_pow_floor(num_rows, PE::radix);
  B_rows_first_level = div_ceil((num_rows - second_level_num_rows) * PE::radix, PE::radix - 1); 
  B_rows_second_level = num_rows - B_rows_first_level;
  unsigned num_levels = log_ceil(num_rows, PE::radix); 
  num_C_partials_level = std::vector<unsigned>(num_levels);
  C_partial_fibers = std::vector<C_Partial_Fiber*>(num_levels * PE::radix);
  C_row_idx = C_row_idx_;
  C_row_ptr = C_row_ptr_;
}

bool Task_Tree::valid() const {
  return !num_C_partials_level.empty();
}

PE_Manager::PE_Manager(const toml::value& parsed_config, Matrix_Data& matrix_data_)
  : matrix_data{matrix_data_}
  , mem_read_ports(2)
  , A_row_ptr_fetcher{matrix_data_.preproc_A_row_ptr}
  , A_row_idx_fetcher{matrix_data_.preproc_A_row_idx}
  , C_row_ptr_fetcher{matrix_data_.preproc_C_row_ptr}
  , A_values_fetcher{matrix_data_.preproc_A_values}
  , B_row_ptr_end_fetcher{matrix_data_.preproc_B_row_ptr_end}
{
  get_config_params(parsed_config);
  reset();
}

void PE_Manager::reset() {
  for (auto& i: mem_read_ports) {
    i.reset();
  }
  for (auto& i: mem_write_ports) {
    i.reset();
  }
  for (auto& i: cache_read_ports) {
    i.reset();
  }
  for (auto& i: cache_write_ports) {
    i.reset();
  }
  prefetch_port.reset();
  A_row_ptr_fetcher.reset();
  A_row_idx_fetcher.reset();
  C_row_ptr_fetcher.reset();
  A_values_fetcher.reset();
  B_row_ptr_end_fetcher.reset();
  A_row_ptr_fetcher.base_addr = matrix_data.preproc_A_row_ptr_addr;
  A_row_idx_fetcher.base_addr = matrix_data.preproc_A_row_idx_addr;
  C_row_ptr_fetcher.base_addr = matrix_data.C_row_ptr_addr;
  A_values_fetcher.base_addr = matrix_data.preproc_A_values_addr;
  B_row_ptr_end_fetcher.base_addr = matrix_data.preproc_B_row_ptr_end_addr;
  read_arbiter = UINT_MAX;
  num_elements_prefetch = 0;
  for (auto pe : PEs) { pe.reset(); };
  std::fill(C_partial_fibers.begin(), C_partial_fibers.end(), C_Partial_Fiber{});
  C_Partial_Fiber::num_fibers = 0;
  task_tree.reset();
  PE::num_mults = 0;
  PE::num_adds = 0;
  PE::num_finished_rows = 0;
  PE::num_C_partial_rows = 0;
  PE::num_C_partial_elements = 0;
  PE::idle_cycles = 0;
  PE::B_data_stalls = 0;
  PE::C_writes = 0;
  PE::max_bytes_write = 0;
  preproc_A_reads = 0;
}

void PE_Manager::update() {
  // send mem request of 1 of the arrays to main memory
  if (!mem_read_ports[0].has_msg_send()) {
    Mem_Request request {};
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
        mem_read_ports[0].add_msg_send(request);
	++preproc_A_reads;
        break;
      }
    }
  }
  // send request of B_row_ptr_end to main memory
  if (!mem_read_ports[1].has_msg_send()) {
    Mem_Request request {};
    request.address = B_row_ptr_end_fetcher.get_fetch_address();
    if (request.valid()) {
      mem_read_ports[1].add_msg_send(request);
      ++preproc_A_reads;
    }
  }
  // send prefetch request
  if (!prefetch_port.has_msg_send()) {
    const auto n = std::min(num_elements_prefetch, prefetched_rows_per_cycle);
    num_elements_prefetch -= n;
    prefetch_port.add_msg_send(n);
  }
  // send cache read requests
  for (unsigned i = 0; i < PEs.size(); ++i) {
    if (cache_read_ports[i].has_msg_send()) continue;
    const auto req = PEs[i].get_cache_request();
    if (req.valid()) { 
      cache_read_ports[i].add_msg_send(req);
    }
  }
  write_data();
  // update PEs
  for (auto& pe : PEs) {
    pe.update();
  }
  allocate_tasks();
  for (auto& p: mem_read_ports) {
    p.transfer();
  }
  for (auto& p: mem_write_ports) {
    p.transfer();
  }
  for (auto& p: cache_read_ports) {
    p.transfer();
  }
  for (auto& p: cache_write_ports) {
    p.transfer();
  }
  prefetch_port.transfer();
}

void PE_Manager::apply() {
  // receive mem responses
  if (mem_read_ports[0].msg_received_valid()) {
    const auto mem_read_resp = mem_read_ports[0].get_msg_received();
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
    mem_read_ports[0].clear_msg_received();
  }
  if (mem_read_ports[1].msg_received_valid()) {
    const auto mem_read_resp = mem_read_ports[1].get_msg_received();
    num_elements_prefetch += B_row_ptr_end_fetcher.receive_data(mem_read_resp.address);
    mem_read_ports[1].clear_msg_received();
  }
  // receive cache data
  for (unsigned i = 0; i < PEs.size(); ++i) {
    if (!cache_read_ports[i].msg_received_valid()) continue;
    PEs[i].receive_cache_response(cache_read_ports[i].get_msg_received());
    cache_read_ports[i].clear_msg_received();
  }
}

PE_Manager::Mem_Port* PE_Manager::get_mem_read_port(std::size_t id) {
  if (id >= mem_read_ports.size()) return nullptr;
  return &mem_read_ports[id];
}

PE_Manager::Mem_Port* PE_Manager::get_mem_write_port(std::size_t id){
  if (id >= mem_write_ports.size()) return nullptr;
  return &mem_write_ports[id];
}

PE_Manager::Mem_Port* PE_Manager::get_cache_read_port(std::size_t id){
  if (id >= cache_read_ports.size()) return nullptr;
  return &cache_read_ports[id];
}

PE_Manager::Mem_Port* PE_Manager::get_cache_write_port(std::size_t id){
  if (id >= cache_write_ports.size()) return nullptr;
  return &cache_write_ports[id];
}

PE_Manager::Prefetch_Port* PE_Manager::get_prefetch_port() {
  return &prefetch_port;
}

bool PE_Manager::finished() const {
  if (PE::num_finished_rows < matrix_data.preproc_A_row_idx.size()) return false;
  // check if all PEs finished writing
  for (const auto& pe : PEs) {
    if (pe.num_bytes_write > 0) return false;
  }
  return true;
}

void PE_Manager::get_config_params(const toml::value& parsed_config) {
  PE::radix = toml::find<unsigned>(parsed_config, "PE_manager", "PE_radix");
  Input_Buffer::buffer_size = toml::find_or(parsed_config, "PE_manager",
                                            "PE_input_buffer_size", 16ul);
  PE::output_buffer_size = toml::find_or(parsed_config, "PE_manager",
                                         "PE_output_buffer_size", 16u);
  const auto num_PEs = toml::find<unsigned>(parsed_config, "PE_manager", "num_PEs");
  mem_write_ports = std::vector<Mem_Port>(num_PEs);
  cache_read_ports = std::vector<Mem_Port>(num_PEs);
  cache_write_ports = std::vector<Mem_Port>(num_PEs);
  PEs = std::vector<PE>(num_PEs, PE{matrix_data});
  const auto task_tree_max_level = 32u / log2_ceil(PE::radix);
  const auto max_partial_fibers = std::max(task_tree_max_level * PE::radix, 2 * num_PEs); 
  C_partial_fibers = std::vector<C_Partial_Fiber>(max_partial_fibers);
  prefetched_rows_per_cycle = toml::find_or(parsed_config, "PE_manager", "prefetched_rows_per_cycle", 4u);
  A_row_ptr_fetcher.buffer_size = toml::find_or(parsed_config, "PE_manager", "A_row_ptr_buffer_size", 128u);
  A_row_idx_fetcher.buffer_size = A_row_ptr_fetcher.buffer_size;
  C_row_ptr_fetcher.buffer_size = A_row_ptr_fetcher.buffer_size;
  A_values_fetcher.buffer_size = toml::find_or(parsed_config, "PE_manager", "A_values_buffer_size", 1024u);
  B_row_ptr_end_fetcher.buffer_size = toml::find_or(parsed_config, "PE_manager", "B_row_ptr_end_buffer_size", 1024u);
}

void PE_Manager::write_data() {
  for (unsigned i = 0; i < PEs.size(); ++i) {
    if (!PEs[i].cur_task.valid()) continue;
    // set initial write address
    if (PEs[i].write_address == invalid_address) {
      if (PEs[i].cur_task.C_partial_fiber) {
	PEs[i].write_address = PEs[i].cur_task.C_partial_fiber->begin;
      } else {
	PEs[i].write_address = matrix_data.C_elements_addr + PEs[i].cur_task.C_row_ptr * element_size;
      }
    }
    if (PEs[i].cur_task.C_partial_fiber) {
      if (cache_write_ports[i].has_msg_send()) continue;
      unsigned num_bytes_write = static_cast<unsigned>(block_size_bytes - PEs[i].write_address % block_size_bytes);
      if (PEs[i].cur_task_finished) {
	assert(PEs[i].num_bytes_write > 0);
	num_bytes_write = std::min(num_bytes_write, PEs[i].num_bytes_write);
      }
      if (PEs[i].num_bytes_write < num_bytes_write) continue;
      Mem_Request req{.address = PEs[i].write_address, .is_write = true};
      cache_write_ports[i].add_msg_send(req);
      PEs[i].write_address += num_bytes_write;
      PEs[i].num_bytes_write -= num_bytes_write;
      PEs[i].cur_task.C_partial_fiber->end += num_bytes_write;
    } else {
      if (mem_write_ports[i].has_msg_send()) continue;
      unsigned num_bytes_write = static_cast<unsigned>(mem_transaction_size - PEs[i].write_address % mem_transaction_size);
      if (PEs[i].cur_task_finished) {
	assert(PEs[i].num_bytes_write > 0);
	num_bytes_write = std::min(num_bytes_write, PEs[i].num_bytes_write);
      }
      if (PEs[i].num_bytes_write < num_bytes_write) continue;
      Mem_Request req{.address = PEs[i].write_address, .is_write = true};
      mem_write_ports[i].add_msg_send(req);
      ++PE::C_writes;
      PEs[i].write_address += num_bytes_write;
      PEs[i].num_bytes_write -= num_bytes_write;
    }
    // switch to new task
    if (PEs[i].cur_task_finished && PEs[i].num_bytes_write == 0) {
      if (PEs[i].next_task.valid()) {
	PEs[i].cur_task = PEs[i].next_task;
	PEs[i].next_task = Task{};
	if (PEs[i].cur_task.C_partial_fiber) {
	  PEs[i].write_address = PEs[i].cur_task.C_partial_fiber->begin;
	} else {
	  PEs[i].write_address = matrix_data.C_elements_addr + PEs[i].cur_task.C_row_ptr * element_size;
	}
	for (auto& buffer : PEs[i].input_buffers) {
	  assert(buffer.num_elems_fetched_cur_task == 0);
	  buffer.num_elems_fetched_cur_task = buffer.col_idx.size();
	}
	PEs[i].cur_task_finished = false;
      } else {
	PEs[i].cur_task = Task{};
	PEs[i].cur_task_finished = false;
	PEs[i].write_address = invalid_address;
      }
    }
  }
}

void PE_Manager::allocate_tasks() {
  for (auto& pe : PEs) {
    if (!pe.cur_task.valid()) {
      pe.cur_task = get_new_task();
      if (!pe.cur_task.valid()) return;
    }
  }
  for (auto& pe : PEs) {
    if (pe.next_task.valid()) {
      pe.next_task = get_new_task();
      if (!pe.next_task.valid()) return;
    }
  }
}

Task PE_Manager::get_new_task() {
  if (!task_tree.valid()) {
    if (A_row_idx_fetcher.finished()) return Task{};
    if (A_row_ptr_fetcher.num_elements < 2 || A_row_idx_fetcher.num_elements == 0 || C_row_ptr_fetcher.num_elements == 0) {
      return Task{};
    }
    const unsigned A_row_idx = A_row_idx_fetcher.front();
    const unsigned C_row_ptr = C_row_ptr_fetcher.front();
    const unsigned num_rows_merge = A_row_ptr_fetcher.at(1) - A_row_ptr_fetcher.front();
    if (num_rows_merge <= PE::radix) {
      if (A_values_fetcher.num_elements < num_rows_merge || B_row_ptr_end_fetcher.num_elements < num_rows_merge) {
	return Task{};
      }
      Task task;
      task.C_row_idx = A_row_idx;
      task.C_row_ptr = C_row_ptr;
      for (unsigned i = 0; i < num_rows_merge; ++i) {
	task.inputs.push_back(get_B_input_fiber());
      }
      A_row_ptr_fetcher.pop();
      A_row_idx_fetcher.pop();
      C_row_ptr_fetcher.pop();
      return task;
    } else {
      A_row_ptr_fetcher.pop();
      A_row_idx_fetcher.pop();
      C_row_ptr_fetcher.pop();
      task_tree.init(num_rows_merge, A_row_idx, C_row_ptr);
    }
  }
  assert(task_tree.valid());
  // create task from task tree
  const auto last_level = task_tree.num_C_partials_level.size() - 1;
  if (task_tree.tree_level == 0) {
    assert(task_tree.B_rows_first_level > 0);
    if (C_Partial_Fiber::num_fibers == C_partial_fibers.size()) {
      return Task{};
    }
    unsigned B_rows_merge = std::min(task_tree.B_rows_first_level, PE::radix);
    if (A_values_fetcher.num_elements < B_rows_merge || B_row_ptr_end_fetcher.num_elements < B_rows_merge) {
      return Task{};
    }
    task_tree.B_rows_first_level -= B_rows_merge;
    const auto C_partial_ptr = get_C_partial_ptr();
    assert(C_partial_ptr != nullptr);
    assert(task_tree.C_partial_fibers[task_tree.num_C_partials_level[0]] == nullptr);
    task_tree.C_partial_fibers[task_tree.num_C_partials_level[0]] = C_partial_ptr;
    Task task;
    task.C_partial_fiber = C_partial_ptr;
    for (unsigned i = 0; i < B_rows_merge; ++i) {
      task.inputs.push_back(get_B_input_fiber());
    }
    ++task_tree.num_C_partials_level[0];
    if (task_tree.num_C_partials_level[0] == PE::radix || task_tree.B_rows_first_level == 0) {
      task_tree.tree_level = 1;
    }
    return task;
  }
  if (task_tree.tree_level == 1) {
    if (task_tree.tree_level == last_level) {
      assert(task_tree.B_rows_second_level + task_tree.num_C_partials_level[0] == PE::radix);
      if (A_values_fetcher.num_elements < task_tree.B_rows_second_level || B_row_ptr_end_fetcher.num_elements < task_tree.B_rows_second_level) {
	return Task{};
      }
      Task task;
      task.C_row_idx = task_tree.C_row_idx;
      task.C_row_ptr = task_tree.C_row_ptr;
      for (unsigned i = 0; i < task_tree.num_C_partials_level[0]; ++i) {
	task.inputs.push_back(Input_Fiber{.A_value = 1.0, .C_partial_fiber = task_tree.C_partial_fibers[i]});
	task_tree.C_partial_fibers[i] = nullptr;
      }
      for (unsigned i = 0; i < task_tree.B_rows_second_level; ++i) {
	task.inputs.push_back(get_B_input_fiber());
      }
      task_tree.reset();
      return task;
    }
    if (C_Partial_Fiber::num_fibers == C_partial_fibers.size()) {
      return Task{};
    }
    const unsigned B_rows_merge = PE::radix - task_tree.num_C_partials_level[0];
    if (A_values_fetcher.num_elements < B_rows_merge || B_row_ptr_end_fetcher.num_elements < B_rows_merge) {
      return Task{};
    }
    const auto C_partial_ptr = get_C_partial_ptr();
    assert(C_partial_ptr != nullptr);
    assert(task_tree.C_partial_fibers[PE::radix + task_tree.num_C_partials_level[1]] == nullptr);
    task_tree.C_partial_fibers[PE::radix + task_tree.num_C_partials_level[1]] = C_partial_ptr;
    Task task;
    task.C_partial_fiber = C_partial_ptr;
    for (unsigned i = 0; i < task_tree.num_C_partials_level[0]; ++i) {
      task.inputs.push_back(Input_Fiber{ .A_value = 1.0, .C_partial_fiber = task_tree.C_partial_fibers[i] });
      task_tree.C_partial_fibers[i] = nullptr;
    }
    for (unsigned i = 0; i < B_rows_merge; ++i) {
      task.inputs.push_back(get_B_input_fiber());
    }
    task_tree.num_C_partials_level[0] = 0;
    ++task_tree.num_C_partials_level[1];
    if (task_tree.num_C_partials_level[1] == PE::radix) {
      ++task_tree.tree_level;
    } else if (task_tree.B_rows_first_level > 0) {
      task_tree.tree_level = 0;
    }
    return task;
  }
  if (task_tree.tree_level < last_level) {
    assert(task_tree.num_C_partials_level[task_tree.tree_level-1] == PE::radix);
    if (C_Partial_Fiber::num_fibers == C_partial_fibers.size()) {
      return Task{};
    }
    const auto C_partial_ptr = get_C_partial_ptr();
    const auto idx = PE::radix * task_tree.tree_level + task_tree.num_C_partials_level[task_tree.tree_level];
    assert(C_partial_ptr != nullptr);
    assert(task_tree.C_partial_fibers[idx] == nullptr);
    task_tree.C_partial_fibers[idx] = C_partial_ptr;
    Task task;
    task.C_partial_fiber = C_partial_ptr;
    for (unsigned i = 0; i < PE::radix; ++i) {
      auto& C_partial = task_tree.C_partial_fibers[(task_tree.tree_level-1) * PE::radix + i];
      task.inputs.push_back(Input_Fiber{ .A_value = 1.0, .C_partial_fiber = C_partial });
      C_partial = nullptr;
    }
    task_tree.num_C_partials_level[task_tree.tree_level-1] = 0;
    ++task_tree.num_C_partials_level[task_tree.tree_level];
    if (task_tree.num_C_partials_level[task_tree.tree_level] == PE::radix) {
      ++task_tree.tree_level;
    } else if (task_tree.B_rows_first_level > 0) {
      task_tree.tree_level = 0;
    } else {
      task_tree.tree_level = 1;
    }
    return task;
  }
  // last level
  assert(task_tree.num_C_partials_level[task_tree.tree_level - 1] == PE::radix);
  Task task;
  task.C_row_idx = task_tree.C_row_idx;
  task.C_row_ptr = task_tree.C_row_ptr;
  for (unsigned i = 0; i < PE::radix; ++i) {
    auto& C_partial = task_tree.C_partial_fibers[(task_tree.tree_level - 1) * PE::radix + i];
    task.inputs.push_back(Input_Fiber{ .A_value = 1.0, .C_partial_fiber = C_partial });
    C_partial = nullptr;
  }
  task_tree.reset();
  return task;
}

Input_Fiber PE_Manager::get_B_input_fiber() {
  assert(B_row_ptr_end_fetcher.num_elements > 0);
  assert(A_values_fetcher.num_elements > 0);
  auto [B_row_ptr_, B_row_end_] = B_row_ptr_end_fetcher.front();
  Input_Fiber input_fiber{.A_value = A_values_fetcher.front(), .B_row_ptr = B_row_ptr_, .B_row_end = B_row_end_, .C_partial_fiber = nullptr};
  A_values_fetcher.pop();
  B_row_ptr_end_fetcher.pop();
  return input_fiber;
}

C_Partial_Fiber* PE_Manager::get_C_partial_ptr() {
  for (unsigned i = 0; i < C_partial_fibers.size(); ++i) {
    if (C_partial_fibers[i].empty()) {
      std::size_t C_region_size = (UINT64_MAX - matrix_data.C_partials_base_addr) / C_partial_fibers.size();
      C_region_size = round_up_multiple(C_region_size, 96ul);
      C_partial_fibers[i].begin = matrix_data.C_partials_base_addr + i * C_region_size;
      C_partial_fibers[i].end = C_partial_fibers[i].begin;
      C_partial_fibers[i].finished = false;
      ++C_Partial_Fiber::num_fibers;
      return &C_partial_fibers[i];
    }
  }
  return nullptr;
}

} // namespace gamma

} // namespace mergeforest_sim
