#include <mergeforest-sim/mergeforest/merge_tree_manager.hpp>
#include <mergeforest-sim/math_utils.hpp>

namespace mergeforest_sim {

namespace mergeforest {

bool Input_Fiber::finished() const {
  return !C_partial_fiber && head_ptr == UINT_MAX && B_num_elements == 0;
}

bool Tree_Level::empty() const {
  return nodes.empty();
}

void Tree_Level::init(unsigned new_task, unsigned num_nodes) {
  task = new_task;
  num_active_nodes = num_nodes;
  for (unsigned i = 0; i != num_nodes; ++i) {
    assert(nodes[i].finished());
    nodes[i].last = false;
  }
}

Merge_Tree::Merge_Tree(Merge_Tree_Manager& parent_)
  : parent{ parent_ }
{
  inputs.assign(parent.merge_tree_size, Input_Fiber{});
  levels.assign(log2_ceil(parent.merge_tree_size) + 1, Tree_Level{});
  unsigned size = 1;
  for (auto& level : levels) {
    level.nodes.assign(size, Fiber_Buffer{});
    size *= 2;
  }
  outputs.assign(levels.size(), Task_Output{});
}

void Merge_Tree::reset() {
  std::ranges::fill(inputs, Input_Fiber{});
  num_active_inputs = 0;
  input_task = 0;
  input_arbiter = UINT64_MAX;
  mult_arbiter = UINT64_MAX;
  for (auto& level : levels) {
    std::ranges::fill(level.nodes, Fiber_Buffer{});
    level.task = UINT_MAX;
    level.num_active_nodes = 0;
  }
  std::ranges::fill(outputs, Task_Output{});
}

bool Merge_Tree::inactive() const {
  return num_active_inputs == 0 && levels[0].task == UINT_MAX;
}

std::size_t Merge_Tree::input_buffer_size(std::size_t idx) const {
  assert(idx < parent.merge_tree_size);
  return levels.back().nodes[idx].size()
    + inputs[idx].B_num_elements + inputs[idx].next_data.size();
}

Cache_Read Merge_Tree::get_request() {
  for (unsigned i = 0; i != inputs.size(); ++i) {
    input_arbiter = inc_mod(input_arbiter, inputs.size());
    auto& input = inputs[input_arbiter];
    if (input.request_sent) { continue; }
    if (input.C_partial_fiber
        && input.head_ptr == UINT_MAX
        && !input.C_partial_fiber->finished()
        && input.C_partial_fiber->head_ptr != UINT_MAX)
    {
      input.head_ptr = input.C_partial_fiber->head_ptr;
    }
    if (input.head_ptr != UINT_MAX && 
        input_buffer_size(input_arbiter) + block_size <= parent.input_buffer_size)
    {
      input.request_sent = true;
      return Cache_Read{.row_ptr = input.head_ptr,
                        .id = static_cast<unsigned>(input_arbiter)};
    }
  }
  return Cache_Read{};
}

void Merge_Tree::receive_response(const Cache_Response& resp) {
  assert(resp.id < parent.merge_tree_size);
  auto& input = inputs[resp.id];
  assert(input.request_sent &&
         input_buffer_size(resp.id) + resp.num_elements <= parent.input_buffer_size);
  if (input.C_partial_fiber) {
    assert(input.C_partial_fiber->data.size() >= resp.num_elements);
    auto& base_node = levels.back().nodes[resp.id];
    auto& buffer = (input_task == levels.back().task) ? base_node : input.next_data;
    buffer.last = false;
    fiber_buffer_transfer(input.C_partial_fiber->data, buffer, resp.num_elements);
    if (input.C_partial_fiber->finished()) {
      assert(resp.row_ptr == UINT_MAX);
      input.C_partial_fiber->head_ptr = UINT_MAX;
      input.C_partial_fiber = nullptr;
    }
  } else {
    input.B_num_elements += resp.num_elements;
  }
  input.head_ptr = resp.row_ptr;
  input.request_sent = false;
  if (input.finished()) {
    assert(num_active_inputs > 0);
    --num_active_inputs;
    if (num_active_inputs == 0) {
      input_task = inc_mod(input_task, static_cast<unsigned>(levels.size()));
    }
  }
}

Address Merge_Tree::get_C_write_address() {
  auto& root_level = levels[0];
  if (root_level.task == UINT_MAX) { return invalid_address; }
  auto& output = outputs[root_level.task];
  if (!output.valid()) { return invalid_address; }
  const auto address = output.get_C_write_address();
  if (!output.valid()) {
    root_level.task = UINT_MAX;
  }
  return address;
}

Cache_Write Merge_Tree::get_C_partial_write() {
  auto& root_level = levels[0];
  if (root_level.task == UINT_MAX) { return Cache_Write{}; }
  auto& output = outputs[root_level.task];
  if (!output.valid()) { return Cache_Write{}; }
  const auto cache_write = output.get_C_partial_write();
  if (!output.valid()) {
    root_level.task = UINT_MAX;
  }
  return cache_write;
}

void Merge_Tree::update() {
  for (unsigned i = 0; i != levels.size() - 1; ++i) {
    update_level(i);
  }
  update_base();
}

void Merge_Tree::update_level(unsigned idx) {
  assert(idx < levels.size() - 1);
  auto& cur_level = levels[idx];
  auto& next_level = levels[idx + 1];
  if (cur_level.task == UINT_MAX) {
    if (next_level.task == UINT_MAX) { return; }
    cur_level.init(next_level.task, (next_level.num_active_nodes + 1) / 2);
  }
  if (cur_level.task != next_level.task) { return; }
  if (idx == 0) {
    update_root();
    return;
  }
  for (unsigned i = 0; i != cur_level.nodes.size(); ++i) {
    auto& src1 = next_level.nodes[2*i];
    auto& src2 = next_level.nodes[2*i + 1];
    auto& dest = cur_level.nodes[i];
    if (dest.size() > parent.merge_tree_merger_width) { continue; }
    if (src1.finished() && src2.finished()) { continue; }
    if (!src1.ready_to_merge(parent.merge_tree_merger_width)
        || !src2.ready_to_merge(parent.merge_tree_merger_width))
    {
      continue;
    }
    if (src1.finished()) {
      fiber_buffer_transfer(src2, dest, parent.merge_tree_merger_width);
      if (src2.finished()) {
        assert(next_level.num_active_nodes > 0);
        --next_level.num_active_nodes;
      }
    } else if (src2.finished()) {
      fiber_buffer_transfer(src1, dest, parent.merge_tree_merger_width);
      if (src1.finished()) {
        assert(next_level.num_active_nodes > 0);
        --next_level.num_active_nodes;
      }
    } else {
      parent.do_merge_add(dest, src1, src2, true);
      if (src1.finished()) {
        assert(next_level.num_active_nodes > 0);
        --next_level.num_active_nodes;
      }
      if (src2.finished()) {
        assert(next_level.num_active_nodes > 0);
        --next_level.num_active_nodes;
      }
    }
    if (next_level.num_active_nodes == 0) {
      next_level.task = UINT_MAX;
    }
    break;
  }
}

void Merge_Tree::update_root() {
  assert(levels[0].task != UINT_MAX);
  assert(levels[0].task == levels[1].task);
  auto& output = outputs[levels[0].task];
  if (output.num_bytes_write >
      (parent.output_buffer_size - parent.merge_tree_merger_width) * element_size)
  {
    return;
  }
  auto& src1 = levels[1].nodes[0];
  auto& src2 = levels[1].nodes[1];
  auto& dest = levels[0].nodes[0];
  assert(!src1.finished() || !src2.finished());
  if (dest.size() > std::max(parent.merge_tree_merger_width, parent.dyn_merger_width)
      || !src1.ready_to_merge(parent.merge_tree_merger_width)
      || !src2.ready_to_merge(parent.merge_tree_merger_width))
  { 
    return;
  }
  auto& buffer = output.C_partial ? output.C_partial->data : dest;
  unsigned num_elements_out = 0;
  if (src1.finished()) {
    num_elements_out = fiber_buffer_transfer(src2, buffer, parent.merge_tree_merger_width);
    if (src2.finished()) {
      assert(levels[1].num_active_nodes > 0);
      --levels[1].num_active_nodes;
    }
  } else if (src2.finished()) {
    num_elements_out = fiber_buffer_transfer(src1, buffer, parent.merge_tree_merger_width);
    if (src1.finished()) {
      assert(levels[1].num_active_nodes > 0);
      --levels[1].num_active_nodes;
    }
  } else {
    num_elements_out = parent.do_merge_add(buffer, src1, src2, true);
    if (src1.finished()) {
      assert(levels[1].num_active_nodes > 0);
      --levels[1].num_active_nodes;
    }
    if (src2.finished()) {
      assert(levels[1].num_active_nodes > 0);
      --levels[1].num_active_nodes;
    }
  }
  if (levels[1].num_active_nodes == 0) {
    levels[1].task = UINT_MAX;
  }
  dest.last = buffer.last;
  if (output.valid()) {
    parent.write_C_output(output, dest, num_elements_out);
  }
}

void Merge_Tree::update_base() {
  auto& base_level = levels.back();
  if (base_level.task == UINT_MAX) {
    // update merge tree base level with new task
    if (num_active_inputs == 0) { return; }
    base_level.task = input_task;
    base_level.num_active_nodes = num_active_inputs;
    for (unsigned i = 0; i != base_level.num_active_nodes; ++i) {
      auto& base_node = base_level.nodes[i];
      auto& input = inputs[i];
      assert(base_node.finished());
      if (input.next_data.finished()) {
        base_node.last = false;
      } else {
        std::swap(base_level.nodes[i], input.next_data);
        if (input.finished()) {
          assert(num_active_inputs > 0);
          --num_active_inputs;
          if (num_active_inputs == 0) {
            input_task = inc_mod(input_task, static_cast<unsigned>(levels.size()));
          }
        }
      }
    }
  }
  // select base node or input to do block mult
  for (unsigned i = 0; i != inputs.size(); ++i) {
    mult_arbiter = inc_mod(mult_arbiter, inputs.size());
    if (inputs[mult_arbiter].B_num_elements == 0) { continue; }
    auto& input = inputs[mult_arbiter];
    auto& buffer = (base_level.task == input_task) ? base_level.nodes[mult_arbiter]
      : input.next_data;
    assert(buffer.size() <= parent.input_buffer_size);
    // do block mult
    auto n = std::min(parent.merge_tree_merger_width, input.B_num_elements);
    input.B_num_elements -= n;
    parent.num_mults += n;
    ++parent.num_block_mults;
    const auto& mat_B = parent.matrix_data.B;
    while (n--) {
      buffer.col_idx.push_back(mat_B->col_idx[input.B_row_ptr]);
      if (parent.matrix_data.compute_result) {
        buffer.values.push_back(input.A_value * mat_B->values[input.B_row_ptr]);
      }
      ++input.B_row_ptr;
    }
    if (base_level.task != input_task) {
      buffer.last = false;
    }
    if (input.finished()) {
      buffer.last = true;
      if (input.next_data.finished()) {
        assert(num_active_inputs > 0);
        --num_active_inputs;
        if (num_active_inputs == 0) {
          input_task = inc_mod(input_task, static_cast<unsigned>(levels.size()));
        }
      }
    }
    break;
  }
}

} // namespace mergeforest

} // namespace mergeforest_sim
