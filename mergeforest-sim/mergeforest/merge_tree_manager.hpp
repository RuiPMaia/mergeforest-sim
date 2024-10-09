#ifndef MERGEFOREST_SIM_MERGE_TREE_MANAGER_HPP
#define MERGEFOREST_SIM_MERGE_TREE_MANAGER_HPP

#include <mergeforest-sim/array_fetcher.hpp>
#include <mergeforest-sim/port.hpp>
#include <mergeforest-sim/matrix_data.hpp>

#include <toml.hpp>

#include <vector>
#include <deque>
#include <climits>

namespace mergeforest_sim {

namespace mergeforest {

struct Fiber_Buffer {
  bool empty() const;
  bool finished() const;
  std::size_t size() const;
  bool ready_to_merge(unsigned size) const;

  std::deque<uint32_t> col_idx{};
  std::deque<double> values{};
  bool last{ true };
};

struct C_Partial_Fiber {
  bool finished() const;

  Fiber_Buffer data;
  unsigned head_ptr{ UINT_MAX };
};

struct Task_Output {
  bool valid() const;
  Address get_C_write_address();
  Cache_Write get_C_partial_write();

  C_Partial_Fiber* C_partial {};
  unsigned C_row_idx{ UINT_MAX };
  unsigned C_row_ptr{ UINT_MAX };
  std::size_t num_bytes_write {};
  Address write_address { invalid_address };
};

struct Fiber_Source {
  bool valid() const;
  bool merge_tree_src() const;
  
  unsigned index { UINT_MAX };
  // value of UINT_MAX -> dynamic node  
  unsigned task_idx { UINT_MAX };
};

struct Dynamic_Tree_Node {
  bool empty() const;

  Fiber_Buffer data;
  Fiber_Source src1;
  Fiber_Source src2;
  Task_Output output;
};

struct Task_Allocator {
  bool empty() const;
  void reset();
  bool all_rows_allocated() const;
  bool last_merge() const;

  unsigned num_B_rows{};
  std::vector<C_Partial_Fiber*> C_partial_fibers;
  std::vector<bool> trees_allocated;
  std::vector<std::pair<Fiber_Source, unsigned>> allocated_sources;
  Task_Output output;
};

struct Task_Tree {
  bool empty() const;
  void reset();

  unsigned tree_level{};
  unsigned B_rows_first_level{};
  unsigned B_rows_second_level{};
  unsigned C_row_idx{ UINT_MAX };
  unsigned C_row_ptr{ UINT_MAX };
  std::vector<unsigned> num_C_partials_level;
  std::vector<C_Partial_Fiber*> C_partial_fibers;
};

struct Merge_Tree;

class Merge_Tree_Manager {
public:
  friend Merge_Tree;
  
  using Mem_Port = Port<Mem_Request, Mem_Response>;
  using Prefetch_Port = Port<Empty_Msg, std::vector<Prefetched_Row>>;
  using Cache_Read_Port = Port<Cache_Read, Cache_Response>;
  using Cache_Write_Port = Port<Cache_Write, unsigned>;

  Merge_Tree_Manager(const toml::value& parsed_config, Matrix_Data& matrix_data_);
  void reset();
  void update();
  void apply();
  bool finished();
  Mem_Port* get_mem_read_port();
  Prefetch_Port* get_prefetch_port();
  Cache_Read_Port* get_cache_read_port(std::size_t id);
  Cache_Write_Port* get_cache_write_port();
  Mem_Port* get_mem_write_port(std::size_t id);
  std::size_t num_mem_ports() const;
  std::size_t num_cache_read_ports() const;
  // config parameters
  unsigned max_prefetched_rows {};
  unsigned merge_tree_size {};
  unsigned max_rows_merge {};
  unsigned merge_tree_merger_width {};
  unsigned merge_tree_merger_num_adds {};
  unsigned num_final_mergers {};
  unsigned dyn_merger_width {};
  unsigned dyn_merger_num_adds {};
  unsigned input_buffer_size {};
  unsigned output_buffer_size {};
  // stats
  std::size_t num_mults {};
  std::size_t num_block_mults {};
  std::size_t merge_tree_num_merges {};
  std::size_t dyn_num_merges {};
  std::size_t merge_tree_num_adds {};
  std::size_t dyn_num_adds {};
  std::size_t num_idle_cycles {};
  std::size_t C_writes {};
  std::size_t preproc_A_reads {};
  std::size_t num_C_partial_rows {};
  std::size_t num_C_partial_elements {};
  std::size_t prefetch_stalls {};
  std::size_t A_data_stalls {};
  std::size_t C_partial_stalls {};
  std::size_t max_write_bytes {};
private:
  void get_config_params(const toml::value& parsed_config);
  void send_A_data_request(); 
  void send_cache_read_requests();
  void write_C_data();
  void write_C_partial_data();
  void update_dynamic_nodes();
  bool fiber_source_ready(const Fiber_Source& src) const;
  Fiber_Buffer& fiber_source_node(const Fiber_Source& src); 
  void fiber_source_reset(Fiber_Source& src);
  void write_C_output(Task_Output& output, Fiber_Buffer& node,
                      unsigned num_elements_out);
  void update_merge_trees();
  void allocate_task();
  bool task_allocator_single_subtask() const;
  bool add_task_merge_tree(unsigned tree_idx); 
  void add_task_dyn_node(unsigned node_idx, unsigned idx_merge); 
  void get_new_task();
  void init_task_tree(unsigned num_rows, unsigned C_row_idx_,
                      unsigned C_row_ptr_);
  C_Partial_Fiber* get_C_partial_fiber();
  void receive_A_data();
  void receive_prefetch_data();
  void receive_cache_data();
  unsigned do_merge_add(Fiber_Buffer& dest, Fiber_Buffer& src1,
                        Fiber_Buffer& src2, bool is_merge_tree);
  
  Matrix_Data& matrix_data;

  Mem_Port mem_read_port;
  Prefetch_Port prefetch_port;
  std::vector<Cache_Read_Port> cache_read_ports;
  Cache_Write_Port cache_write_port;
  std::vector<Mem_Port> mem_write_ports;

  Array_Fetcher<uint32_t> A_row_ptr_fetcher;
  Array_Fetcher<uint32_t> A_row_idx_fetcher;
  Array_Fetcher<uint32_t> C_row_ptr_fetcher;
  Array_Fetcher<double> A_values_fetcher;
  unsigned read_arbiter {UINT_MAX};
  std::deque<Prefetched_Row> prefetched_B_rows;

  std::vector<Merge_Tree> merge_trees;
  std::vector<Dynamic_Tree_Node> dyn_nodes; 
  std::vector<C_Partial_Fiber> C_partial_fibers;
  Task_Allocator task_allocator;
  Task_Tree task_tree;
  unsigned C_partial_write_idx {UINT_MAX};
  C_Partial_Fiber* C_partial_head_ptr {nullptr};
  std::size_t write_arbiter {UINT64_MAX};
};

struct Input_Fiber {
  bool finished() const;

  C_Partial_Fiber* C_partial_fiber {};
  double A_value {};
  uint32_t B_row_ptr { UINT32_MAX };
  unsigned head_ptr { UINT_MAX };
  bool request_sent { false };
  unsigned B_num_elements {};
  Fiber_Buffer next_data {};
};

struct Tree_Level {
  bool empty() const;
  void init(unsigned new_task, unsigned num_nodes);

  std::vector<Fiber_Buffer> nodes;
  unsigned task { UINT_MAX };
  unsigned num_active_nodes {};
};

struct Merge_Tree {
  Merge_Tree(Merge_Tree_Manager& parent_);
  void reset();
  bool inactive() const;
  std::size_t input_buffer_size(std::size_t idx) const;
  Cache_Read get_request();
  void receive_response(const Cache_Response& resp);
  Address get_C_write_address();
  Cache_Write get_C_partial_write();
  void update();
  void update_level(unsigned idx);
  void update_root();
  void update_base();

  Merge_Tree_Manager& parent;

  std::vector<Input_Fiber> inputs;
  unsigned num_active_inputs {};
  unsigned input_task {};
  std::size_t input_arbiter {UINT64_MAX};
  std::size_t mult_arbiter {UINT64_MAX};
  std::vector<Tree_Level> levels;
  std::vector<Task_Output> outputs;
};

unsigned fiber_buffer_transfer(Fiber_Buffer& src, Fiber_Buffer& dest,
                               std::size_t num_elements);

} // namespace mergeforest

} // namespace mergeforest_sim

#endif // MERGEFOREST_SIM_MERGE_TREE_MANAGER_HPP
