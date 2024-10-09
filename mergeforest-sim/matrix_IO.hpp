#ifndef MERGEFOREST_SIM_MAT_IO_HPP
#define MERGEFOREST_SIM_MAT_IO_HPP

#include <mergeforest-sim/sparse_matrix.hpp>

#include <string>

namespace mergeforest_sim {

Spmat_Csr read_matrix_market_file(const std::string& filename);

} // namespace mergeforest_sim

#endif // MERGEFOREST_SIM_MAT_IO_HPP
