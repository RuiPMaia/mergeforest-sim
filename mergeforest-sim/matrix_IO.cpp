#include <mergeforest-sim/matrix_IO.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <cstdint>

namespace mergeforest_sim {

enum class Format {
  array,
  coordinate
};

enum class Symmetry {
  general,
  symmetric,
  hermitian,
  skew_symmetric
};

enum class Type {
  pattern,
  real,
  integer,
  complex
};

struct Matrix_Market_Header {
  Format format;
  Type type;
  Symmetry symmetry;
};

Matrix_Market_Header read_matrix_market_header(std::ifstream& input_stream) {
  Matrix_Market_Header header;
  std::stringstream stream;
  std::string buffer;
  // tokenize header
  std::getline(input_stream, buffer);
  stream << buffer;
  std::string identifier, object, format, type, symmetry;
  stream >> identifier >> object >> format >> type >> symmetry;
  // check if header is valid
  if (stream.fail() || identifier != "%%MatrixMarket") {
    throw std::runtime_error("invalid MatrixMarket header");
  }
  if (object != "matrix") {
    throw std::runtime_error("invalid MatrixMarket object type [" + object + "]");
  }
  if (format != "coordinate") {
    throw std::runtime_error("invalid MatrixMarket storage format [" + format + "]");
  } 
  header.format = Format::coordinate;
  
  if (type == "pattern") {
    header.type = Type::pattern;
  } else if (type == "real") {
    header.type = Type::real;
  } else if (type == "integer") {
    header.type = Type::integer;
  } else if (type == "complex") {
    header.type = Type::complex;
  } else {
    throw std::runtime_error("invalid MatrixMarket data type [" + type + "]");
  }
  if (symmetry == "general") {
    header.symmetry = Symmetry::general;
  } else if (symmetry == "symmetric") {
    header.symmetry = Symmetry::symmetric;
  } else if (symmetry == "hermitian") {
    header.symmetry = Symmetry::hermitian;
  } else if (symmetry == "skew-symmetric") {
    header.symmetry = Symmetry::skew_symmetric;
  } else {
    throw std::runtime_error("invalid MatrixMarket symmetry type [" + symmetry + "]");
  }
  if (header.type != Type::complex && header.symmetry == Symmetry::hermitian) {
    throw std::runtime_error("invalid MatrixMarket combination [" + type + ", hermitian]");
  }
  if (header.type == Type::pattern && header.symmetry == Symmetry::skew_symmetric) {
    throw std::runtime_error("invalid MatrixMarket combination [pattern, skew-symmetric]");
  }
  if (header.type == Type::complex) {
    throw std::runtime_error("MatrixMarket data type [complex] not supported");
  }
  return header;
}

Spmat_Csr read_matrix_market_file(const std::string& filename) {
  Spmat_Csr mtx;
  std::ifstream input_stream(filename);
  if (!input_stream)
    throw std::runtime_error("unable to open file \"" + filename + "\" for reading");
  // read header
  Matrix_Market_Header header = read_matrix_market_header(input_stream);
  std::string line;
  // skip over comments
  do {
    std::getline(input_stream, line);
  } while (line[0] == '%');
  // tokenize size line
  std::istringstream(line) >> mtx.num_rows >> mtx.num_cols >> mtx.nnz;
  
  std::vector<std::tuple<uint32_t, uint32_t, double>> coo;
  
  for (unsigned i = 0; i < mtx.nnz; ++i) {
    // tokenize each data line
    std::getline(input_stream, line);
    std::stringstream stream(line);
    unsigned row_idx = 0;
    unsigned col_idx = 0;
    double value = 1.0;
    stream >> row_idx >> col_idx;
    if (row_idx < 1 || col_idx < 1 || row_idx > mtx.num_rows || col_idx > mtx.num_cols) {
      throw std::runtime_error("MatrixMarket invalid index");
    }
    if (header.type == Type::real || header.type == Type::integer) {
      stream >> value;
    }
    if (stream.fail()) {
      throw std::runtime_error("MatrixMarket invalid data");
    }
    coo.push_back(std::make_tuple(row_idx - 1, col_idx - 1, value));
    if (row_idx != col_idx) {
      if (header.symmetry == Symmetry::symmetric) {
        coo.push_back(std::make_tuple(col_idx - 1, row_idx - 1, value));
      }
      else if (header.symmetry == Symmetry::skew_symmetric) {
        coo.push_back(std::make_tuple(col_idx - 1, row_idx - 1, -value));
      }
    }
  }
  std::sort(coo.begin(), coo.end());
  // nnz might change with symmetry
  mtx.nnz = coo.size();
  mtx.row_ptr = std::vector<uint32_t>(mtx.num_rows + 1);
  mtx.col_idx = std::vector<uint32_t>(mtx.nnz);
  mtx.values = std::vector<double>(mtx.nnz);
  mtx.row_ptr[0] = 0;
  uint32_t prev_row = 0;
  for (unsigned i = 0; i < mtx.nnz; ++i) {
    const auto [row, col_idx, value] = coo[i];
    mtx.col_idx[i] = col_idx;
    mtx.values[i] = value;
    for (unsigned j = prev_row; j < row; ++j) {
      mtx.row_ptr[j + 1] = i;
    }
    prev_row = row;
  }
  for (unsigned i = prev_row; i < mtx.num_rows; ++i) {
    mtx.row_ptr[i + 1] = static_cast<unsigned>(mtx.nnz);
  }
  mtx.col_idx[mtx.nnz - 1] = std::get<1>(coo[mtx.nnz - 1]);
  mtx.values[mtx.nnz - 1] = std::get<2>(coo[mtx.nnz - 1]);
  return mtx;
}

} // namespace mergeforest_sim
