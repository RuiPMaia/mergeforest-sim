#include <mergeforest-sim/simulator.hpp>

namespace mergeforest_sim {

template<typename T>
concept Arch = requires(T a, bool compute_res) {
  {a.run_simulation(compute_res)} -> std::same_as<Spmat_Csr>;
};

struct Arch_Visitor {
  Arch_Visitor(bool compute_result) : compute_res{compute_result} {} 

  Spmat_Csr operator()(Arch auto& a) {
    return a.run_simulation(compute_res);
  }

  Spmat_Csr operator()([[maybe_unused]] auto& a) {
    return Spmat_Csr{};
  }

  bool compute_res {};
};

Simulator::Simulator(const std::string& config_file,
		     const std::string& out_path_)
  : parsed_config(toml::parse(config_file))
  , out_path{out_path_}
{
  const auto arch_str = toml::find<std::string>(parsed_config, "arch");
  if (arch_str == "my_arch") {
    arch.emplace<MergeForest>(parsed_config, matrix_data, out_path);
  }
  else if (arch_str == "gamma") {
    arch.emplace<Gamma>(parsed_config, matrix_data, out_path);
  }
  else { 
    throw std::runtime_error("Error: architecture \""
			     + arch_str + "\" not implemented");
  }
}

void Simulator::set_mats(const Spmat_Csr& A, const Spmat_Csr& B) {
  matrix_data.A = &A;
  matrix_data.B = &B;
}

Spmat_Csr Simulator::run_simulation(bool compute_result) {
  return std::visit(Arch_Visitor{compute_result}, arch);
}

} // namespace mergeforest_sim
