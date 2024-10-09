#include <mergeforest-sim/mergeforest.hpp>
#include <mergeforest-sim/math_utils.hpp>

#include <spdlog/spdlog.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

namespace mergeforest_sim {

MergeForest::MergeForest(const toml::value& parsed_config_,
		 Matrix_Data& matrix_data_,
		 const std::string& out_path_)
  : parsed_config{parsed_config_}
  , matrix_data{matrix_data_}
  , out_path{out_path_}
  , merge_tree_manager{parsed_config, matrix_data_}
  , linked_list_cache{parsed_config, matrix_data_}
  , main_mem{parsed_config}
{
  main_mem.set_num_ports(1 + linked_list_cache.num_mem_ports()
                         + merge_tree_manager.num_mem_ports());
  // port connections
  merge_tree_manager.get_mem_read_port()->connect(main_mem.get_port(0));
  std::size_t port_idx = 1;
  for (std::size_t i = 0; i != linked_list_cache.num_mem_ports(); ++i) {
    linked_list_cache.get_mem_port(i)->connect(main_mem.get_port(port_idx));
    ++port_idx;
  }
  for (std::size_t i = 0; i != merge_tree_manager.num_mem_ports(); ++i) {
    merge_tree_manager.get_mem_write_port(i)->connect(
      main_mem.get_port(port_idx));
    ++port_idx;
  }
  merge_tree_manager.get_prefetch_port()->connect(
    linked_list_cache.get_prefetch_port());
  for (std::size_t i = 0; i != merge_tree_manager.num_cache_read_ports(); ++i) {
    merge_tree_manager.get_cache_read_port(i)->connect(
      linked_list_cache.get_read_port(i));
  }
  merge_tree_manager.get_cache_write_port()->connect(
    linked_list_cache.get_write_port());
}

void MergeForest::print_progress() {
  if (merge_tree_manager.num_mults == 0) {
    fmt::print("progress:   0.00%\r");
  } else {
    const auto progress = static_cast<double>(merge_tree_manager.num_mults)
      / static_cast<double>(matrix_data.num_mults) * 100.0;
    fmt::print("progress: {:6.2f}%\r", progress);
  }
  fflush(stdout);
}

Spmat_Csr MergeForest::run_simulation(bool compute_result) {
  matrix_data.compute_result = compute_result;
  matrix_data.preprocess_mats();
  matrix_data.set_physical_addrs();
  reset();
  // simulation loop
  for (;;) {
    merge_tree_manager.update();
    linked_list_cache.update();
    main_mem.update();
    linked_list_cache.apply();
    merge_tree_manager.apply();
    if (cycles % progress_interval == 0) {
      print_progress();
    }
    ++cycles;
    if (merge_tree_manager.finished() && main_mem.inactive()) {
      break;
    }
  }
  fmt::print("progress: 100.00%\n");
  check_valid_simulation();
  print_stats();
  if (compute_result) {
    matrix_data.spGEMM_check_result();
    return matrix_data.C;
  }
  return Spmat_Csr{};
}

void MergeForest::reset() {
  merge_tree_manager.reset();
  linked_list_cache.reset();
  main_mem.reset();
  cycles = 0;
}

void MergeForest::check_valid_simulation() {
  if (matrix_data.num_mults != merge_tree_manager.num_mults) {
    spdlog::error("Number of multiplications doesn't match the expected value");
  }
  const auto num_adds = merge_tree_manager.merge_tree_num_adds
    + merge_tree_manager.dyn_num_adds;
  if (matrix_data.num_mults != matrix_data.C.nnz + num_adds) {
    spdlog::error("Number of additions doesn't match the expected value");
  }
  const auto num_reads = merge_tree_manager.preproc_A_reads
    + linked_list_cache.preproc_A_reads + linked_list_cache.B_reads
    + linked_list_cache.C_partial_reads;
  if (main_mem.read_requests != num_reads) {
    spdlog::error("Number of reads in Main Memory doesn't match the rest of the system");
  }
  const auto num_writes = merge_tree_manager.C_writes
    + linked_list_cache.C_partial_writes;
  if (main_mem.write_requests != num_writes) {
    spdlog::error("Number of writes in Main Memory doesn't match the rest of the system");
  }
  if (linked_list_cache.C_partial_reads != linked_list_cache.C_partial_writes) {
    spdlog::error("Number of reads and writes of C partial data doesn't match");
  }
  const auto B_bytes_read = linked_list_cache.B_elements_read * element_size;
  if (B_bytes_read < matrix_data.min_bytes_B_data) {
    spdlog::error("Number of B bytes read too small");
  }
  if (B_bytes_read > matrix_data.max_bytes_B_data) {
    spdlog::error("Number of B bytes read too big");
  }
  if (linked_list_cache.B_reads < matrix_data.B_data_min_reads) {
    spdlog::error("Number of B reads read too small");
  }
  if (linked_list_cache.B_reads > matrix_data.B_data_max_reads) {
    spdlog::error("Number of B reads read too big");
  }
  if (linked_list_cache.fetched_rows + linked_list_cache.reused_rows
      != matrix_data.preproc_B_row_ptr_end.size())
  {
    spdlog::error("Number of fetched and reused B rows "
                  "doesn't match total number of B rows");
  }
}

void MergeForest::print_stats() {
  if (out_path.empty()) {
    print_stats_impl(std::cout);
  } else {
    std::ofstream of;
    of.open(out_path.data());
    print_stats_impl(of);
  }
}

void MergeForest::print_stats_impl(std::ostream& os) {
  const double period_ns = toml::find_or(parsed_config, "clock_period_ns", 1.0);
  const auto exec_time_ns = static_cast<double>(cycles) * period_ns;
  const auto exec_time_ms = exec_time_ns * 1e-6;
  const auto Gflops = static_cast<double>(matrix_data.num_mults) / exec_time_ns;
  const auto block_mults_ratio = ratio(matrix_data.num_mults,
                                       merge_tree_manager.num_block_mults
                                       * merge_tree_manager.merge_tree_merger_width) * 100.0;
  const auto num_adds = merge_tree_manager.merge_tree_num_adds + merge_tree_manager.dyn_num_adds;
  const auto merge_tree_adds_ratio = ratio(merge_tree_manager.merge_tree_num_adds,
                                           merge_tree_manager.merge_tree_num_merges
                                           * merge_tree_manager.merge_tree_merger_num_adds) * 100.0;
  const auto dyn_adds_ratio = ratio(merge_tree_manager.dyn_num_adds,
                                    merge_tree_manager.dyn_num_merges
                                    * merge_tree_manager.dyn_merger_num_adds) * 100.0;
  const auto dyn_merges_per_cycle = ratio(merge_tree_manager.dyn_num_merges, cycles);
  const auto num_trees = merge_tree_manager.num_cache_read_ports();
  const auto idle_cycles_ratio = ratio(merge_tree_manager.num_idle_cycles,
                                       cycles * num_trees) * 100.0;
  const auto A_data_stalls_ratio = ratio(merge_tree_manager.A_data_stalls,
                                         cycles) * 100.0;
  const auto C_partial_stalls_ratio = ratio(merge_tree_manager.C_partial_stalls,
                                         cycles) * 100.0;
  const auto cache_bandwidth = ratio(linked_list_cache.reads + linked_list_cache.writes, cycles);
  const auto active_blocks_avg = ratio(linked_list_cache.num_active_blocks_avg,
                                       linked_list_cache.num_samples);
  const auto inactive_blocks_avg = ratio(linked_list_cache.num_inactive_blocks_avg,
                                         linked_list_cache.num_samples);
  const auto C_partial_blocks_avg = ratio(linked_list_cache.num_C_partial_blocks_avg,
                                          linked_list_cache.num_samples);
  const auto free_blocks_avg = ratio(linked_list_cache.num_free_blocks_avg,
                                     linked_list_cache.num_samples);
  const auto active_blocks_ratio = ratio(active_blocks_avg, linked_list_cache.num_blocks) * 100.0;
  const auto inactive_blocks_ratio =
    ratio(inactive_blocks_avg, linked_list_cache.num_blocks) * 100.0;
  const auto C_partial_blocks_ratio =
    ratio(C_partial_blocks_avg, linked_list_cache.num_blocks) * 100.0;
  const auto free_blocks_ratio =
    ratio(free_blocks_avg, linked_list_cache.num_blocks) * 100.0;
  const auto mem_traffic = main_mem.read_requests + main_mem.write_requests;
  const auto mem_traffic_bytes = static_cast<double>(mem_traffic * mem_transaction_size);
  const auto bandwidth = mem_traffic_bytes / exec_time_ns;
  const auto op_intensity =
    static_cast<double>(matrix_data.num_mults) / mem_traffic_bytes;
  const auto preproc_A_reads = merge_tree_manager.preproc_A_reads
    + linked_list_cache.preproc_A_reads;
  const auto preproc_A_bytes_read = sizeof(uint32_t) * (
    matrix_data.preproc_A_row_ptr.size()
    + matrix_data.preproc_A_row_idx.size()
    + matrix_data.preproc_C_row_ptr.size()
    + 2 * matrix_data.preproc_B_row_ptr_end.size())
    + sizeof(double) * matrix_data.preproc_A_values.size();
  const auto B_bytes_read = linked_list_cache.B_elements_read * element_size;
  const auto C_partial_bytes_rw = linked_list_cache.C_partial_reads
    * mem_transaction_size;
  const auto mem_bytes_read = preproc_A_bytes_read + B_bytes_read + C_partial_bytes_rw;
  const auto unused_read_bytes_ratio =
    unused_bytes_ratio(main_mem.read_requests, mem_bytes_read);
  const auto C_bytes_write = matrix_data.C.nnz * element_size;
  const auto mem_bytes_write = C_bytes_write + C_partial_bytes_rw;
  const auto unused_write_bytes_ratio = unused_bytes_ratio(main_mem.write_requests,
                                                           mem_bytes_write);
  const auto unused_A_bytes_ratio = unused_bytes_ratio(preproc_A_reads,
                                                       preproc_A_bytes_read);
  const auto unused_B_bytes_ratio = unused_bytes_ratio(linked_list_cache.B_reads, B_bytes_read);
  const auto unused_C_bytes_ratio = unused_bytes_ratio(merge_tree_manager.C_writes,
                                                       C_bytes_write);
  const auto total_unused_bytes_ratio = unused_bytes_ratio(mem_traffic,
                                                           mem_bytes_read + mem_bytes_write);

  fmt::print(os, "*---Simulation Results---*\n");
  fmt::print(os, "Config file: {}\n", parsed_config.location().file_name());
  fmt::print(os, "Num cycles: {}\n", cycles);
  fmt::print(os, "Clock period: {} ns\n", period_ns);
  fmt::print(os, "Execution time: {:.4f} ms\n", exec_time_ms);
  fmt::print(os, "GFlops: {:.4f}\n", Gflops);
  fmt::print(os, "*---Merge_Tree_Manager---*\n");
  fmt::print(os, "Number flops (mults): {}\n", matrix_data.num_mults);
  fmt::print(os, "Number block mults: {} ({:.4f}%) utilization\n",
             merge_tree_manager.num_block_mults, block_mults_ratio);
  fmt::print(os, "Number adds : {}\n", num_adds);
  fmt::print(os, "Number merge tree merges : {} ({:.4f}% adder utilization)\n",
             merge_tree_manager.merge_tree_num_merges,
             merge_tree_adds_ratio);
  fmt::print(os, "Number dynamic merges : {} ({:.4f}% adder utilization)\n",
             merge_tree_manager.dyn_num_merges,
             dyn_adds_ratio);
  fmt::print(os, "Dynamic merges per cycle: {:.4f}\n", dyn_merges_per_cycle);
  fmt::print(os, "Idle cycles: {} ({:.4f}%)\n", merge_tree_manager.num_idle_cycles,
	     idle_cycles_ratio);
  fmt::print(os, "A data stalls: {} ({:.4f}%)\n", merge_tree_manager.A_data_stalls,
	     A_data_stalls_ratio);
  fmt::print(os, "C partial stalls: {} ({:.4f}%)\n", merge_tree_manager.C_partial_stalls,
	     C_partial_stalls_ratio);
  fmt::print(os, "C partial rows: {}\n", merge_tree_manager.num_C_partial_rows);
  fmt::print(os, "C partial elements: {}\n",
	     merge_tree_manager.num_C_partial_elements);
  fmt::print(os, "Max write bytes: {}\n", merge_tree_manager.max_write_bytes);
  fmt::print(os, "*---Linked List Cache---*\n");
  fmt::print(os, "Cache reads: {}\n", linked_list_cache.reads);
  fmt::print(os, "Cache writes: {}\n", linked_list_cache.writes);
  fmt::print(os, "Cache bandwidth: {:.4f} blocks/cycle\n",
	     cache_bandwidth);
  fmt::print(os, "Fetched rows: {}\n", linked_list_cache.fetched_rows);
  fmt::print(os, "Reused rows: {}\n", linked_list_cache.reused_rows);
  fmt::print(os, "Evicted rows: {}\n", linked_list_cache.evictions);
  fmt::print(os, "Max active rows: {}\n",
	     linked_list_cache.stats_max_active_rows);
  fmt::print(os, "Max inactive rows: {}\n",
	     linked_list_cache.stats_max_inactive_rows);
  fmt::print(os, "Average active blocks: {:.4f} ({:.4f}%)\n",
	     active_blocks_avg, active_blocks_ratio);
  fmt::print(os, "Average inactive blocks: {:.4f} ({:.4f}%)\n",
	     inactive_blocks_avg, inactive_blocks_ratio);
  fmt::print(os, "Average C_partial blocks: {:.4f} ({:.4f}%)\n",
	     C_partial_blocks_avg, C_partial_blocks_ratio);
  fmt::print(os, "Average free blocks: {:.4f} ({:.4f}%)\n",
	     free_blocks_avg, free_blocks_ratio);

  fmt::print(os, "Max free lists: {}\n",
	     linked_list_cache.max_free_lists);
  fmt::print(os, "Max fetched rows: {}\n",
	     linked_list_cache.stats_max_fetched_rows);
  fmt::print(os, "Max outstanding reqs: {}\n",
	     linked_list_cache.stats_max_outstanding_reqs);
  fmt::print(os, "*---Main Memory---*\n");
  fmt::print(os, "Memory bandwidth: {:.4f} GB/s\n", bandwidth);
  fmt::print(os, "Operational intensity: {:.4f} flop/byte\n", op_intensity);
  fmt::print(os,
	     "Memory traffic: {} transactions ({:.4f} MB) ({:.4f}% unused)\n",
	     mem_traffic, reqs_to_MB(mem_traffic), total_unused_bytes_ratio);
  fmt::print(os, "Memory reads: {} ({:.4f} MB) ({:.4f}% unused)\n",
	     main_mem.read_requests, reqs_to_MB(main_mem.read_requests),
	     unused_read_bytes_ratio);
  fmt::print(os, "Memory writes: {} ({:.4f} MB) ({:.4f}% unused)\n",
	     main_mem.write_requests, reqs_to_MB(main_mem.write_requests),
	     unused_write_bytes_ratio);
  fmt::print(os, "A data reads: {} ({:.4f} MB) ({:.4f}% unused)\n",
	     preproc_A_reads, reqs_to_MB(preproc_A_reads),
	     unused_A_bytes_ratio);
  fmt::print(os, "B data reads: {} ({:.4f} MB) ({:.4f}% unused)\n",
	     linked_list_cache.B_reads,
	     reqs_to_MB(linked_list_cache.B_reads),
	     unused_B_bytes_ratio);
  fmt::print(os, "B data min reads: {} ({:.4f} MB)\n",
	     matrix_data.B_data_min_reads,
	     reqs_to_MB(matrix_data.B_data_min_reads));
  fmt::print(os, "B data max reads: {} ({:.4f} MB)\n",
	     matrix_data.B_data_max_reads_fiber_cache,
	     reqs_to_MB(matrix_data.B_data_max_reads_fiber_cache));
  fmt::print(os, "C partial reads/writes: {} ({:.4f} MB) (0% unused)\n",
	     linked_list_cache.C_partial_reads,
	     reqs_to_MB(linked_list_cache.C_partial_reads));
  fmt::print(os, "C data writes: {} ({:.4f} MB) ({:.4f}% unused)\n",
	     merge_tree_manager.C_writes,
	     reqs_to_MB(merge_tree_manager.C_writes), unused_C_bytes_ratio);
  fmt::print(os, "A data bytes read: {}\n", preproc_A_bytes_read);
  fmt::print(os, "B data bytes read: {}\n", B_bytes_read);
  fmt::print(os, "C data bytes written: {}\n", C_bytes_write);
}

} // namespace mergeforest_sim
