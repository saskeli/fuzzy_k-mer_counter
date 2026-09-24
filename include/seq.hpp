#pragma once

#include <array>

#include "util.hpp"

namespace fkc {
using namespace fkc::internal;
template <uint8_t w_count>
class sequence {
 private:
  std::array<uint64_t, w_count> words_;
  uint64_t size_;

 public:
  sequence() : words_(), size_() {}

  void append(char c) {
    uint8_t v = nuc_to_v[c];
    for (size_t i = w_count - ONE; i > ZERO; --i) {
      words_[i] = (words_[i] << TWO) | (words_[i - ONE] >> SIXTYTWO);
    }
    words_[0] = (words_[0] << 2) | v;
    ++size_;
  }

  const uint64_t& operator[](size_t i) const { return words_[i]; }

  size_t length() const { return size_; }
};

};  // namespace fkc