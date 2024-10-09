#include <mergeforest-sim/gamma.hpp>
#include <mergeforest-sim/math_utils.hpp>

#include <spdlog/spdlog.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include <iostream>

namespace mergeforest_sim {

Gamma::Gamma(const toml::value& parsed_config_,
	     Matrix_Data& matrix_data_,
	     const std::string& out_path_)
  : parsed_config{parsed_config_}
  , matrix_data{ matrix_data_ }
  , out_path{out_path_}
  , PE_manager{parsed_config_, matrix_data_}
  , fiber_cache{parsed_config_, matrix_data_}
  , main_mem{parsed_config_}
{
  const auto fiber_cache_num_mem_ports = toml::find<std::size_t>(parsed_config,
								 "fiber_cache",
								 "num_mem_ports");
  const auto num_PEs = toml::find<std::size_t>(parsed_config, "PE_manager", "num_PEs");
  main_mem.set_num_ports(2 + fiber_cache_num_mem_ports + num_PEs);
  // port connections
  PE_manager.get_mem_read_port(0)->connect(main_mem.get_port(0));
  PE_manager.get_mem_read_port(1)->connect(main_mem.get_port(1));
  for (unsigned i = 0; i < fiber_cache_num_mem_ports; ++i) {
    fiber_cache.get_mem_port(i)->connect(main_mem.get_port(i + 2));
  }
  for (unsigned i = 0; i < num_PEs; ++i) {
    PE_manager.get_mem_write_port(i)->connect(main_mem.get_port(i + fiber_cache_num_mem_ports + 2));
  }
  PE_manager.get_prefetch_port()->connect(fiber_cache.get_prefetch_port());
  for (unsigned i = 0; i < num_PEs; ++i) {
    PE_manager.get_cache_read_port(i)->connect(fiber_cache.get_read_port(i));
    PE_manager.get_cache_write_port(i)->connect(fiber_cache.get_write_port(i));
  }
}

void Gamma::print_progress() {
  if (gamma::PE::num_mults == 0) {
    fmt::print("progress:   0.00%\r");
  } else {
    const auto progress = static_cast<double>(gamma::PE::num_mults) /
      static_cast<double>(matrix_data.num_mults) * 100.0;
    fmt::print("progress: {:6.2f}%\r", progress);
  }
  fflush(stdout);
} 

Spmat_Csr Gamma::run_simulation(bool compute_result) {
  matrix_data.compute_result = compute_result;
  matrix_data.preprocess_mats();
  matrix_data.set_physical_addrs();
  reset();
  // simulation loop
  for (;;) {
    PE_manager.update();
    fiber_cache.update();
    main_mem.update();
    fiber_cache.apply();
    PE_manager.apply();
    if (cycles % progress_interval == 0) {
      print_progress();
    }
    ++cycles;
    if (PE_manager.finished() && fiber_cache.inactive() && main_mem.inactive()) {
      break;
    }
  }
  fmt::print("progress: 100.00%\n");
  fiber_cache.B_data_reads *= 3;
  fiber_cache.C_partial_reads *= 3;
  fiber_cache.C_partial_writes *= 3;
  check_valid_simulation();
  if (compute_result) {
    matrix_data.spGEMM_check_result();
    print_stats();
    return matrix_data.C;
  }
  print_stats();
  return Spmat_Csr{};
}

void Gamma::reset() {
  PE_manager.reset();
  fiber_cache.reset();
  main_mem.reset();
  cycles = 0;
}

void Gamma::check_valid_simulation() {
  if (matrix_data.num_mults != gamma::PE::num_mults) {
    spdlog::error(R"(Error in simulation: number of multiplications doesn't
      match the expected value\n)");
  }
  if (gamma::PE::num_mults - gamma::PE::num_adds != matrix_data.C.nnz) {
    spdlog::error(R"(Error in simulation: number of multiplications and
      additions doesn'tmatch the nnz of the result\n)");
  }
  if (fiber_cache.B_data_reads < matrix_data.B_data_min_reads_fiber_cache) {
    spdlog::error("Error in simulation: number of B bytes read too small\n");
  }
  if (fiber_cache.B_data_reads > matrix_data.B_data_max_reads_fiber_cache) {
    spdlog::error("Error in simulation: number of B bytes read too big\n");
  }
  if (fiber_cache.C_partial_reads != fiber_cache.C_partial_writes) {
    spdlog::error(R"(Error in simulation: number of C bytes read doesn't match
      the number of C bytes written\n)");
  }
  if (main_mem.read_requests !=
      PE_manager.preproc_A_reads
      + fiber_cache.B_data_reads
      + fiber_cache.C_partial_reads)
  {
    spdlog::error(R"(Error in simulation: memory reads don't match PE manager
      and fiber cache reads\n)");
  }
  if (main_mem.write_requests !=
      gamma::PE::C_writes + fiber_cache.C_partial_writes)
  {
    spdlog::error(R"(Error in simulation: memory reads don't match PE manager
      and fiber cache reads\n)");
  }
}

void Gamma::print_stats() {
  if (out_path.empty()) {
    print_stats_impl(std::cout);
  } else {
    std::ofstream of;
    of.open(out_path.data());
    print_stats_impl(of);
  }
}

void Gamma::print_stats_impl(std::ostream& os) {
  const double period_ns = toml::find_or(parsed_config, "clock_period_ns", 1.0);
  const auto exec_time_ns = static_cast<double>(cycles) * period_ns;
  const auto exec_time_ms = exec_time_ns * 1e-6;
  const auto Gflops = static_cast<double>(matrix_data.num_mults) / exec_time_ns;
  const auto num_PEs = toml::find<std::size_t>(parsed_config,
					       "PE_manager",
					       "num_PEs");
  const auto idle_cycles_ratio = ratio(gamma::PE::idle_cycles, cycles * num_PEs) * 100.0;
  const auto B_data_stalls_ratio = ratio(gamma::PE::B_data_stalls, cycles * num_PEs) * 100.0;
  const auto write_stalls_ratio = ratio(gamma::PE::write_stalls, cycles * num_PEs) * 100.0;

  const auto mem_traffic = main_mem.read_requests + main_mem.write_requests;
  const auto mem_traffic_bytes = static_cast<double>(mem_traffic * mem_transaction_size);
  const auto bandwidth = mem_traffic_bytes / exec_time_ns;
  const auto op_intensity =
    static_cast<double>(matrix_data.num_mults) / mem_traffic_bytes;
  const auto cache_bandwidth =
    static_cast<double>(fiber_cache.reads + fiber_cache.writes)
    / static_cast<double>(cycles);
  const auto B_blocks_avg = ratio(fiber_cache.B_blocks_avg, fiber_cache.num_samples);
  const auto C_partial_blocks_avg = ratio(fiber_cache.C_partial_blocks_avg, fiber_cache.num_samples);
  const auto free_blocks_avg = static_cast<double>(fiber_cache.num_blocks)
    - B_blocks_avg - C_partial_blocks_avg;
  const auto B_blocks_ratio = ratio(B_blocks_avg, fiber_cache.num_blocks) * 100.0;
  const auto C_partial_blocks_ratio = ratio(C_partial_blocks_avg, fiber_cache.num_blocks) * 100.0;
  const auto free_blocks_ratio = ratio(free_blocks_avg, fiber_cache.num_blocks) * 100.0;
  const auto preproc_A_bytes_read = sizeof(uint32_t) * (
    matrix_data.preproc_A_row_ptr.size()
    + matrix_data.preproc_A_row_idx.size()
    + matrix_data.preproc_C_row_ptr.size()
    + 2 * matrix_data.preproc_B_row_ptr_end.size())
    + sizeof(double) * matrix_data.preproc_A_values.size();
  const auto mem_bytes_read = preproc_A_bytes_read
    + (fiber_cache.B_data_reads + fiber_cache.C_partial_reads)
    * mem_transaction_size;
  const auto unused_read_bytes_ratio = unused_bytes_ratio(main_mem.read_requests,
						  mem_bytes_read);
  const auto C_data_bytes_write = matrix_data.C.nnz * element_size;
  const auto mem_bytes_write = C_data_bytes_write
    + fiber_cache.C_partial_writes * mem_transaction_size;
  const auto unused_write_bytes_ratio = unused_bytes_ratio(main_mem.write_requests,
						   mem_bytes_write);
  const auto unused_A_bytes_ratio = unused_bytes_ratio(PE_manager.preproc_A_reads,
					       preproc_A_bytes_read);
  const auto unused_C_bytes_ratio = unused_bytes_ratio(gamma::PE::C_writes,
					       C_data_bytes_write);
  const auto total_unused_bytes_ratio = unused_bytes_ratio(mem_traffic,
					     mem_bytes_read + mem_bytes_write);
  const auto cache_hit_rate =
    static_cast<double>(fiber_cache.read_hits)
    / static_cast<double>(fiber_cache.reads) * 100.0;

  fmt::print(os, "*---Simulation Results---*\n");
  fmt::print(os, "Config file: {}\n", parsed_config.location().file_name());
  fmt::print(os, "Num cycles: {}\n", cycles);
  fmt::print(os, "Clock period: {} ns\n", period_ns);
  fmt::print(os, "Execution time: {:.4f} ms\n", exec_time_ms);
  fmt::print(os, "GFlops: {:.4f}\n", Gflops);
  fmt::print(os, "*---Processing Elements---*\n");
  fmt::print(os, "Number flops (mults): {}\n", matrix_data.num_mults);
  fmt::print(os, "Number adds : {}\n", gamma::PE::num_adds);
  fmt::print(os, "Idle cycles: {} ({:.4f}%)\n", gamma::PE::idle_cycles,
	     idle_cycles_ratio);
  fmt::print(os, "B data stalls: {} ({:.4f}%)\n", gamma::PE::B_data_stalls,
	     B_data_stalls_ratio);
  fmt::print(os, "Write stalls: {} ({:.4f}%)\n", gamma::PE::write_stalls,
	     write_stalls_ratio);
  fmt::print(os, "C partial rows: {}\n", gamma::PE::num_C_partial_rows);
  fmt::print(os, "C partial elements: {}\n",
	     gamma::PE::num_C_partial_elements);
  fmt::print(os, "Max bytes write: {}\n", gamma::PE::max_bytes_write);
  fmt::print(os, "*---Fiber Cache---*\n");
  fmt::print(os, "Fiber cache reads: {}\n", fiber_cache.reads);
  fmt::print(os, "Fiber cache writes: {}\n", fiber_cache.writes);
  fmt::print(os, "Fiber cache read hits: {} ({:.4f}% hit rate)\n",
	     fiber_cache.read_hits, cache_hit_rate);
  fmt::print(os, "Fiber cache bandwidth: {:.4f} blocks/cycle\n",
	     cache_bandwidth);
  fmt::print(os, "Average B blocks: {:.4f} ({:.4f}%)\n", B_blocks_avg, B_blocks_ratio);
  fmt::print(os, "Average C partial blocks: {:.4f} ({:.4f}%)\n",
             C_partial_blocks_avg, C_partial_blocks_ratio);
  fmt::print(os, "Average free blocks: {:.4f} ({:.4f}%)\n",
	     free_blocks_avg, free_blocks_ratio);
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
	     PE_manager.preproc_A_reads, reqs_to_MB(PE_manager.preproc_A_reads),
	     unused_A_bytes_ratio);
  fmt::print(os, "B data reads: {} ({:.4f} MB (0% unused)\n",
	     fiber_cache.B_data_reads, reqs_to_MB(fiber_cache.B_data_reads));
  fmt::print(os, "B data min reads: {} ({:.4f} MB)\n",
	     matrix_data.B_data_min_reads_fiber_cache,
	     reqs_to_MB(matrix_data.B_data_min_reads_fiber_cache));
  fmt::print(os, "B data max reads: {} ({:.4f} MB)\n",
	     matrix_data.B_data_max_reads_fiber_cache,
	     reqs_to_MB(matrix_data.B_data_max_reads_fiber_cache));
  fmt::print(os, "C partial reads/writes: {} ({:.4f} MB) (0% unused)\n",
	     fiber_cache.C_partial_reads,
	     reqs_to_MB(fiber_cache.C_partial_reads));
  fmt::print(os, "C data writes: {} ({:.4f} MB) ({:.4f}% unused)\n",
	     gamma::PE::C_writes,
	     reqs_to_MB(gamma::PE::C_writes), unused_C_bytes_ratio);
  fmt::print(os, "A data bytes read: {}\n", preproc_A_bytes_read);
  fmt::print(os, "C data bytes written: {}\n", C_data_bytes_write);
}

} // namespace mergeforest_sim
