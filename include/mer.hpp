#include <bitset>
#include <cstdint>
#include <cstring>
#include <iostream>

#include "util.hpp"

namespace fkc {

class mer {
 private:
  uint64_t symbols_;
  uint64_t mask_;
  uint8_t size_;

 public:
  mer(const char* chars) : symbols_(0), mask_(0), size_(std::strlen(chars)) {
    if (size_ > 32) {
      std::cerr << "Maximum k-mer length is 32. " << chars << " is too long"
                << std::endl;
      exit(1);
    }

    for (uint8_t i = 0; i < size_; ++i) {
      char c = chars[i];
      uint64_t v = internal::nuc_to_v[c];
      if (internal::v_to_nuc[v] != c && c != 'N') {
        std::cerr << "Invalid mer: " << chars << std::endl;
        exit(1);
      }
      symbols_ = (symbols_ << 2) | v;
      mask_ <<= 2;
      mask_ |= 0b11 * (c != 'N');
    }
  }

  mer(uint64_t symbols, uint64_t mask, uint8_t size)
      : symbols_(symbols), mask_(mask), size_(size) {}

  bool operator==(const mer& rhs) const {
    return (symbols_ == rhs.symbols_) && (mask_ == rhs.mask_) &&
           (size_ == rhs.size_);
  }

  bool operator!=(const mer& rhs) const { return not(*this == rhs); }

  mer rc() const {
    uint64_t rc_sym = 0;
    uint64_t rc_mask = 0;
    uint8_t* rc_bytes = reinterpret_cast<uint8_t*>(&rc_sym);
    uint8_t* rc_m_bytes = reinterpret_cast<uint8_t*>(&rc_mask);
    const uint8_t* symbol_bytes = reinterpret_cast<const uint8_t*>(&symbols_);
    const uint8_t* mask_bytes = reinterpret_cast<const uint8_t*>(&mask_);
    for (uint8_t i = 0; i < 7; ++i) {
      rc_bytes[7 - i] = internal::byte_rc[symbol_bytes[i]];
      rc_m_bytes[7 - i] = internal::byte_rev[mask_bytes[i]];
    }
    uint8_t un_pad = 64 - (size_ * 2);
    rc_sym >>= un_pad;
    rc_mask >>= un_pad;
    return {rc_sym & rc_mask, rc_mask, size_};
  }

  bool match(uint64_t seq) const {
    return symbols_ == (seq & mask_);
  }

  uint8_t length() const { return size_; }

  template <class OS_T>
  OS_T& print(OS_T& out) const {
    for (size_t i = size_ * 2 - 2; i < 64; i -= 2) {
      char c = internal::v_to_nuc[(symbols_ >> i) & 0b11];
      if (((mask_ >> i) & 0b11) == 0) {
        c = 'N';
      }
      out << c;
    }
    return out;
  }
};

}  // namespace fkc
