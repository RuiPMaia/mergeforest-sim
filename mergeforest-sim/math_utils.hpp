#ifndef MERGEFOREST_SIM_COMMON_H
#define MERGEFOREST_SIM_COMMON_H

#include <concepts>
#include <bit>
#include <limits>
#include <cmath>

namespace mergeforest_sim {

template<std::unsigned_integral T>
constexpr T log2_ceil(T number) {
  return T(std::bit_width(number)) - T(std::has_single_bit(number));
}

template<std::integral T>
constexpr T round_up_multiple(T number, T multiple) {
  if (number % multiple == 0) return number;
  else return number + multiple - number % multiple;
}

template<std::integral T>
constexpr T round_down_multiple(T number, T multiple) {
  return number - number % multiple;
}

template<std::integral T>
constexpr T inc_mod(T number, T divisor) {
  ++number;
  return (number >= divisor) ? static_cast<T>(0) : number;
}

template<std::integral T>
constexpr T div_ceil(T number, T divisor) {
  return (number - 1) / divisor + 1;
}

template<std::unsigned_integral T>
constexpr T log_ceil(T number, T base) {
  T result = 0;
  T aux = 1;
  while (aux < number) {
    aux *= base;
    ++result;
  }
  return result;
}

template<std::integral T>
T pow_2(T exp) {
    return T(1) << exp;
}

template<std::unsigned_integral T>
constexpr T nearest_pow_floor(T number, T base) {
  T result = T(1);
  for (;;) {
    T aux = result * base;
    if (aux > number) break;
    result = aux;
  }
  return result;
}

template<std::floating_point T>
constexpr bool almost_equal(T a, T b, T c = 1.0) {
  return std::fabs(a - b) <= std::numeric_limits<T>::epsilon()
    * std::fmax(std::fabs(a), std::fabs(b)) * c;
}

double reqs_to_MB(std::size_t reqs);

double unused_bytes_ratio(std::size_t reqs, std::size_t bytes);

double ratio(std::size_t a, std::size_t b);

double ratio(double a, std::size_t b);

} // namespace mergeforest_sim

#endif //MERGEFOREST_SIM_COMMON_H
