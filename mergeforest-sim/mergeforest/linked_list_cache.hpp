#ifndef MERGEFOREST_SIM_LINKED_LIST_CACHE_HPP
#define MERGEFOREST_SIM_LINKED_LIST_CACHE_HPP

#include <mergeforest-sim/array_fetcher.hpp>
#include <mergeforest-sim/mergeforest/matB_fetcher.hpp>
#include <mergeforest-sim/matrix_data.hpp>
#include <mergeforest-sim/port.hpp>

#include <toml.hpp>

#include <vector>
#include <deque>
#include <unordered_map>

namespace mergeforest_sim {

namespace mergeforest {

struct Linked_list_Node {
  unsigned num_elements {};
  unsigned next {UINT_MAX};
  bool last {true};
  bool C_partial_row {false};
};

struct Active_Row {
  bool valid() const;

  unsigned row_head {};
  unsigned num_uses {};
  unsigned num_blocks {};
};

struct Inactive_Row {
  bool valid() const;

  unsigned B_row_ptr {UINT_MAX};
  unsigned row_head {};
  unsigned num_blocks {};
  unsigned prev {};
  unsigned next {};
};

class Linked_List_Cache {
public:
  using Mem_Port = Port<Mem_Request, Mem_Response>;
  using Prefetch_Port = Port<std::vector<Prefetched_Row>, Empty_Msg>;
  using Cache_Read_Port = Port<Cache_Response, Cache_Read>;
  using Cache_Write_Port = Port<unsigned, Cache_Write>;

  Linked_List_Cache(const toml::value& parsed_config,
		    const Matrix_Data& matrix_data_);
  void reset();
  void update();
  void apply();
  Mem_Port* get_mem_port(std::size_t id);
  Prefetch_Port* get_prefetch_port();
  Cache_Read_Port* get_read_port(std::size_t id);
  Cache_Write_Port* get_write_port();
  std::size_t num_mem_ports() const;
  // config parameters
  std::size_t num_blocks {};
  unsigned max_active_rows {};
  unsigned inactive_rows_assoc {};
  unsigned inactive_rows_num_sets {};
  unsigned num_banks {};
  unsigned prefetched_rows_per_cycle {};
  unsigned sample_interval {};
  // stats
  std::size_t reads {};
  std::size_t writes {};
  std::size_t preproc_A_reads {};
  std::size_t B_reads {};
  std::size_t B_elements_read {};
  std::size_t C_partial_reads {};
  std::size_t C_partial_writes {};
  std::size_t reused_rows {};
  std::size_t fetched_rows {};
  std::size_t evictions {};
  std::size_t num_active_blocks_avg {};
  std::size_t num_inactive_blocks_avg {};
  std::size_t num_C_partial_blocks_avg {};
  std::size_t num_free_blocks_avg {};
  std::size_t num_samples {};
  std::size_t max_free_lists {};
  std::size_t stats_max_active_rows {};
  std::size_t stats_max_inactive_rows {};
  std::size_t stats_max_fetched_rows {};
  std::size_t stats_max_outstanding_reqs {};

private:
  void get_config_params(const toml::value& parsed_config);
  unsigned add_new_row(uint32_t B_row_ptr, uint32_t B_row_end);
  void write_B_row_data();
  void receive_read_requests();
  void send_read_responses();
  void finish_pending_reqs(unsigned ptr);
  void update_cache_block(unsigned ptr);
  unsigned write_C_partial_row(Cache_Write request);
  unsigned allocate_block();
  bool free_inactive_row();
  void add_to_inactive_rows(const std::pair<uint32_t, Active_Row>& active_row);
  void inactive_rows_list_remove(unsigned ptr);
  void sample_cache_utilization();

  const Matrix_Data& matrix_data;

  std::vector<Mem_Port> mem_ports;
  Prefetch_Port prefetch_port;
  std::vector<Cache_Read_Port> read_ports;
  Cache_Write_Port write_port;
  std::size_t arbiter {UINT64_MAX};

  Array_Fetcher<std::pair<uint32_t, uint32_t>> B_row_ptr_end_fetcher;
  MatB_Fetcher matB_fetcher;
  std::unordered_multimap<unsigned, std::pair<unsigned, unsigned>> pending_reqs;
  std::vector<std::deque<Cache_Response>> finished_reqs;

  std::unordered_map<uint32_t, Active_Row> active_rows;
  std::vector<Inactive_Row> inactive_rows_cache;
  std::vector<Linked_list_Node> row_data_list;
  std::deque<unsigned> free_list_heads;
  unsigned inactive_rows_list_head {UINT_MAX};
  unsigned inactive_rows_list_tail {UINT_MAX};
  std::size_t num_inactive_rows {};
  unsigned C_partial_row_ptr {UINT_MAX};
  std::size_t num_active_blocks {};
  std::size_t num_inactive_blocks {};
  std::size_t num_C_partial_blocks {};
  std::size_t num_free_blocks {};
  std::size_t num_fetching_blocks {};
  unsigned cycles {};
};

} // namespace mergeforest

} // namespace mergeforest_sim

#endif // MERGEFOREST_SIM_LINKED_LIST_CACHE_HPP
