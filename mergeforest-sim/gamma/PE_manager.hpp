#ifndef MERGEFOREST_SIM_GAMMA_PE_MANAGER_HPP
#define MERGEFOREST_SIM_GAMMA_PE_MANAGER_HPP

#include <mergeforest-sim/array_fetcher.hpp>
#include <mergeforest-sim/matrix_data.hpp>
#include <mergeforest-sim/port.hpp>

#include <toml.hpp>

#include <vector>
#include <utility>

namespace mergeforest_sim {

namespace gamma {

struct C_Partial_Fiber {
  bool empty() const;
  bool is_finished() const;
  
  inline static unsigned num_fibers;

  std::deque<uint32_t> col_idx;
  std::deque<double> values;
  Address begin {invalid_address};
  Address end {invalid_address};
  bool finished {};
};

struct Input_Fiber {
  bool finished() const;

  double A_value {};
  uint32_t B_row_ptr {};
  uint32_t B_row_end {};
  C_Partial_Fiber* C_partial_fiber {};
};

struct Task {
  bool valid () const;
  std::vector<Input_Fiber> inputs;
  uint32_t C_row_ptr {UINT32_MAX};
  uint32_t C_row_idx {UINT32_MAX};
  C_Partial_Fiber* C_partial_fiber {};
};

struct Input_Buffer {
  inline static std::size_t buffer_size;
  std::size_t num_elements_received {};
  std::size_t num_elems_fetched_cur_task {};
  std::deque<std::tuple<Address, unsigned, bool>> pending_reqs;
  std::deque<uint32_t> col_idx;
  std::deque<double> values;
};

struct PE {
  PE(Matrix_Data& matrix_data_);
  void reset();
  Mem_Request get_cache_request();
  void receive_cache_response(Mem_Response mem_response);
  void update();
  //config params
  inline static unsigned radix;
  inline static unsigned output_buffer_size;
  // stats
  inline static std::size_t num_mults;
  inline static std::size_t num_adds;
  inline static std::size_t num_finished_rows;
  inline static std::size_t num_C_partial_rows;
  inline static std::size_t num_C_partial_elements;
  inline static std::size_t idle_cycles;
  inline static std::size_t B_data_stalls;
  inline static std::size_t write_stalls;
  inline static std::size_t C_writes;
  inline static unsigned max_bytes_write;
  Matrix_Data& matrix_data;
  
  Task cur_task;
  Task next_task;
  bool cur_task_finished {};
  uint32_t C_col_idx {UINT32_MAX};
  double C_value {};
  std::vector<Input_Buffer> input_buffers;
  std::size_t read_arbiter {UINT64_MAX};
  Address write_address {invalid_address};
  unsigned num_bytes_write {};
};

struct Task_Tree {
  void reset();
  void init(unsigned num_rows, unsigned C_row_idx_, unsigned C_row_ptr_);
  bool valid() const;
  
  unsigned tree_level {};
  unsigned B_rows_first_level {};
  unsigned B_rows_second_level {};
  unsigned C_row_idx {UINT_MAX};
  unsigned C_row_ptr {UINT_MAX};
  std::vector<unsigned> num_C_partials_level;
  std::vector<C_Partial_Fiber*> C_partial_fibers;
};
  
class PE_Manager {
public:
  using Mem_Port = Port<Mem_Request, Mem_Response>;
  using Prefetch_Port = Port<std::size_t, Empty_Msg>;
  
  PE_Manager(const toml::value& parsed_config, Matrix_Data& matrix_data_); 
  void reset();
  void update();
  void apply();
  Mem_Port* get_mem_read_port(std::size_t id);
  Mem_Port* get_mem_write_port(std::size_t id);
  Mem_Port* get_cache_read_port(std::size_t id);
  Mem_Port* get_cache_write_port(std::size_t id);
  Prefetch_Port* get_prefetch_port();
  bool finished() const;
  // stats
  std::size_t preproc_A_reads {};
private:
  void get_config_params(const toml::value& parsed_config);
  void write_data();
  void allocate_tasks();
  Task get_new_task();
  Input_Fiber get_B_input_fiber();
  C_Partial_Fiber* get_C_partial_ptr();

  Matrix_Data& matrix_data;

  std::vector<Mem_Port> mem_read_ports;
  std::vector<Mem_Port> mem_write_ports;
  std::vector<Mem_Port> cache_read_ports; 
  std::vector<Mem_Port> cache_write_ports;
  Prefetch_Port prefetch_port; 

  Array_Fetcher<uint32_t> A_row_ptr_fetcher;
  Array_Fetcher<uint32_t> A_row_idx_fetcher;
  Array_Fetcher<uint32_t> C_row_ptr_fetcher;
  Array_Fetcher<double> A_values_fetcher;
  Array_Fetcher<std::pair<uint32_t,uint32_t>> B_row_ptr_end_fetcher;
  unsigned read_arbiter {UINT_MAX};
  std::size_t num_elements_prefetch {};
  
  std::vector<PE> PEs;
  std::vector<C_Partial_Fiber> C_partial_fibers;
  Task_Tree task_tree; 
  // config parameters
  std::size_t prefetched_rows_per_cycle {};
}; 
  
} // namespace gamma

} // namespace mergeforest_sim

#endif // MERGEFOREST_SIM_GAMMA_PE_MANAGER_HPP
