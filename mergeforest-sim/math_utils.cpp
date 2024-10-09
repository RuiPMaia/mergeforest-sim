#include <mergeforest-sim/math_utils.hpp>
#include <mergeforest-sim/port.hpp>

namespace mergeforest_sim {

double reqs_to_MB(std::size_t reqs) {
  return static_cast<double>(reqs * mem_transaction_size) * 1E-6;
}

double unused_bytes_ratio(std::size_t reqs, std::size_t bytes) {
  const auto reqs_bytes = reqs * mem_transaction_size;
  return static_cast<double>(reqs_bytes - bytes)
    / static_cast<double>(reqs_bytes) * 100.0;
}

double ratio(std::size_t a, std::size_t b) {
  if (b == 0) { return 0.0; }
  return static_cast<double>(a) / static_cast<double>(b);
}

double ratio(double a, std::size_t b) {
  if (b == 0) { return 0.0; }
  return a / static_cast<double>(b);
}

} // namespace mergeforest_sim
