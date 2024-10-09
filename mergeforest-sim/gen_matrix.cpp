#include <mergeforest-sim/gen_matrix.hpp>

#include <fmt/os.h>

#include <vector>
#include <utility>
#include <string>
#include <random>
#include <cassert>

namespace mergeforest_sim {

void gen_RMat(const std::string& out_path, unsigned num_nodes, unsigned num_edges,
              double A, double B, double C, unsigned seed) {
  using Edge = std::pair<unsigned, unsigned>;

  std::mt19937 rnd(seed);
  assert (A + B + C < 1.0);
  std::vector<Edge> edge_list;
  edge_list.reserve(num_edges);
  // sum of parameters (probabilities)
  std::vector<double> sumA, sumAB, sumAC, sumABC;  // up to 2^128 vertices ~ 3.4e38
  sumA.reserve(128);
  sumAB.reserve(128);
  sumAC.reserve(128);
  sumABC.reserve(128);
  for (unsigned i = 0; i < 128; ++i) {
    std::uniform_real_distribution dist{0.5, 1.5};
    const double a = A * (dist(rnd));
    const double b = B * (dist(rnd));
    const double c = C * (dist(rnd));
    const double d = (1.0 - (A + B + C)) * (dist(rnd));
    const double abcd = a + b + c + d;
    sumA.emplace_back(a / abcd);
    sumAB.emplace_back((a + b) / abcd);
    sumAC.emplace_back((a + c) / abcd);
    sumABC.emplace_back((a + b + c) / abcd);
  }
  std::uniform_real_distribution uni_dist{0.0, 1.0};
  unsigned depth = 0, collisions = 0, cnt = 0, pct_done = 0;
  const unsigned edge_gap = num_edges / 100 + 1;
  for (unsigned edge = 0; edge < num_edges; ) {
    unsigned rng_X = num_nodes, rng_Y = num_nodes, off_X = 0, off_Y = 0;
    depth = 0;
    // recurse the matrix
    while (rng_X > 1 || rng_Y > 1) {
      const double rnd_prob = uni_dist(rnd);
      if (rng_X > 1 && rng_Y > 1) {
        if (rnd_prob < sumA[depth]) {
          rng_X /= 2;
          rng_Y /= 2;
        } else if (rnd_prob < sumAB[depth]) {
          off_X += rng_X / 2;
          rng_X -= rng_X / 2;
          rng_Y /= 2;
        } else if (rnd_prob < sumABC[depth]) {
          off_Y += rng_Y / 2;
          rng_X /= 2;
          rng_Y -= rng_Y / 2;
        } else {
          off_X += rng_X / 2;
          off_Y += rng_Y / 2;
          rng_X -= rng_X / 2;
          rng_Y -= rng_Y / 2;
        }
      } else if (rng_X > 1) { // row vector
        if (rnd_prob < sumAC[depth]) {
          rng_X /= 2;
          rng_Y /= 2;
        } else {
          off_X += rng_X / 2;
          rng_X -= rng_X / 2;
          rng_Y /= 2;
        }
      } else if (rng_Y > 1) { // column vector
        if (rnd_prob < sumAB[depth]) {
          rng_X /= 2;
          rng_Y /= 2;
        } else {
          off_Y += rng_Y / 2;
          rng_X /= 2;
          rng_Y -= rng_Y / 2;
        }
      } else { assert(false); }
      ++depth;
    }
    // add edge
    if (off_X == off_Y) {
      ++collisions;
    } else {
      Edge new_edge {off_X, off_Y};
      const auto it = std::lower_bound(edge_list.begin(), edge_list.end(), new_edge);
      if (*it == new_edge) {
        ++collisions;
      } else {
        edge_list.insert(it, new_edge);
        ++cnt;
        if (cnt > edge_gap) {
          cnt = 0;
          fmt::print("{}% edges\r", ++pct_done);
          fflush(stdout);
        }
        ++edge;
      }
    }
  }
  fmt::print("RMat: nodes:{}, edges:{}, Iterations:{}, Collisions:{} ({:.1f}%).\n",
             num_nodes, num_edges, num_edges + collisions, collisions,
             100.0 * static_cast<double>(collisions)
             / static_cast<double>(num_edges + collisions));
  auto out = fmt::output_file(out_path);
  out.print("%%MatrixMarket matrix coordinate pattern general\n");
  out.print("%seed: {}\n", seed);
  out.print("{} {} {}\n", num_nodes, num_nodes, num_edges);
  for (const auto& edge : edge_list) {
    out.print("{} {}\n", edge.first + 1, edge.second + 1);
  }
}

} // namespace mergeforest_sim
