#ifndef MERGEFOREST_SIM_GEN_MAT_HPP
#define MERGEFOREST_SIM_GEN_MAT_HPP

#include <string>

namespace mergeforest_sim {

// R-MAT Generator. The modes is based on the recursive descent into a 2x2
// matrix [A,B; C, 1-(A+B+C)].
// See: R-MAT Generator: A Recursive Model for Graph Mining. 
// D. Chakrabarti, Y. Zhan and C. Faloutsos, in SIAM Data Mining 2004. 
// URL: http://www.cs.cmu.edu/~deepay/mywww/papers/siam04.pdf
void gen_RMat(const std::string& out_path, unsigned num_nodes, unsigned num_edges,
              double A, double B, double C, unsigned seed = 0);

} // namespace mergeforest_sim

#endif // MERGEFOREST_SIM_GEN_MAT_HPP
