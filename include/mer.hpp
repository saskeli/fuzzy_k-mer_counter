#pragma once

#include <bitset>
#include <cstdint>
#include <cstring>
#include <iostream>

#include "util.hpp"

namespace fkc {
using namespace fkc::internal;
class mer {
 private:
  uint64_t symbols_;
  uint64_t mask_;
  uint8_t size_;

  friend class long_mer;

 public:
  mer(const char* chars)
      : symbols_(ZERO), mask_(ZERO), size_(std::strlen(chars)) {
    if (size_ > THIRTYTWO) {
      std::cerr << "Maximum k-mer length for 'short' mers is 32. " << chars
                << " is too long." << std::endl;
      exit(1);
    }

    for (uint8_t i = ZERO; i < size_; ++i) {
      char c = chars[i];
      uint64_t v = internal::nuc_to_v[c];
      if (internal::v_to_nuc[v] != c && c != 'N') {
        std::cerr << "Invalid mer: " << chars << std::endl;
        exit(1);
      }
      symbols_ = (symbols_ << TWO) | v;
      mask_ <<= TWO;
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
    uint64_t rc_sym = ZERO;
    uint64_t rc_mask = ZERO;
    uint8_t* rc_bytes = reinterpret_cast<uint8_t*>(&rc_sym);
    uint8_t* rc_m_bytes = reinterpret_cast<uint8_t*>(&rc_mask);
    const uint8_t* symbol_bytes = reinterpret_cast<const uint8_t*>(&symbols_);
    const uint8_t* mask_bytes = reinterpret_cast<const uint8_t*>(&mask_);
    for (uint8_t i = 0; i < 7; ++i) {
      rc_bytes[7 - i] = internal::byte_rc[symbol_bytes[i]];
      rc_m_bytes[7 - i] = internal::byte_rev[mask_bytes[i]];
    }
    uint8_t un_pad = W - (size_ * TWO);
    rc_sym >>= un_pad;
    rc_mask >>= un_pad;
    return {rc_sym & rc_mask, rc_mask, size_};
  }

  template <class seq_t>
  bool match(const seq_t& seq) const {
    return symbols_ == (seq[0] & mask_);
  }

  uint8_t length() const { return size_; }

  template <class OS_T>
  OS_T& print(OS_T& out) const {
    for (size_t i = size_ * TWO - TWO; i < W; i -= TWO) {
      char c = internal::v_to_nuc[(symbols_ >> i) & 0b11];
      if (((mask_ >> i) & 0b11) == ZERO) {
        c = 'N';
      }
      out << c;
    }
    return out;
  }
};

class long_mer {
 private:
  uint64_t high_bits_;
  uint64_t low_bits_;
  uint64_t high_mask_;
  uint64_t low_mask_;

  uint8_t size_;

 public:
  long_mer(const char* chars)
      : high_bits_(ZERO),
        low_bits_(ZERO),
        high_mask_(ZERO),
        low_mask_(ZERO),
        size_(std::strlen(chars)) {
    if (size_ > W) {
      std::cerr << "Maximum k-mer length for 'long' mers is 64. " << chars
                << " is too long." << std::endl;
      exit(1);
    }

    for (uint8_t i = ZERO; i < size_; ++i) {
      char c = chars[i];
      uint64_t v = internal::nuc_to_v[c];
      if (internal::v_to_nuc[v] != c && c != 'N') {
        std::cerr << "Invalid mer: " << chars << std::endl;
        exit(1);
      }
      high_bits_ = (high_bits_ << TWO) | (low_bits_ >> SIXTYTWO);
      low_bits_ = (low_bits_ << TWO) | v;
      high_mask_ = (high_mask_ << TWO) | (low_mask_ >> SIXTYTWO);
      low_mask_ <<= TWO;
      low_mask_ |= 0b11 * (c != 'N');
    }
  }

  long_mer(uint64_t high_bits, uint64_t low_bits, uint64_t high_mask,
           uint64_t low_mask, uint8_t size)
      : high_bits_(high_bits),
        low_bits_(low_bits),
        high_mask_(high_mask),
        low_mask_(low_mask),
        size_(size) {}

  bool operator==(const long_mer& rhs) const {
    return (high_bits_ == rhs.high_bits_) && (low_bits_ == rhs.low_bits_) &&
           (high_mask_ == rhs.high_mask_) && (low_mask_ == rhs.low_mask_) &&
           (size_ == rhs.size_);
  }

  bool operator!=(const long_mer& rhs) const { return not(*this == rhs); }

  long_mer rc() const {
    if (size_ <= 32) {
      mer m = {low_bits_, low_mask_, size_};
      m = m.rc();
      return {ZERO, m.symbols_, ZERO, m.mask_, size_};
    }
    uint64_t rc_h_sym = ZERO;
    uint64_t rc_l_sym = ZERO;
    uint64_t rc_h_mask = ZERO;
    uint64_t rc_l_mask = ZERO;
    uint8_t* rc_h_bytes = reinterpret_cast<uint8_t*>(&rc_h_sym);
    uint8_t* rc_l_bytes = reinterpret_cast<uint8_t*>(&rc_l_sym);
    uint8_t* rc_hm_bytes = reinterpret_cast<uint8_t*>(&rc_h_mask);
    uint8_t* rc_lm_bytes = reinterpret_cast<uint8_t*>(&rc_l_mask);
    const uint8_t* h_bytes = reinterpret_cast<const uint8_t*>(&high_bits_);
    const uint8_t* l_bytes = reinterpret_cast<const uint8_t*>(&low_bits_);
    const uint8_t* h_mask_bytes = reinterpret_cast<const uint8_t*>(&high_mask_);
    const uint8_t* l_mask_bytes = reinterpret_cast<const uint8_t*>(&low_mask_);
    for (uint8_t i = 0; i < 7; ++i) {
      rc_h_bytes[7 - i] = byte_rc[l_bytes[i]];
      rc_l_bytes[7 - i] = byte_rc[h_bytes[i]];
      rc_hm_bytes[7 - i] = byte_rev[h_mask_bytes[i]];
      rc_lm_bytes[7 - i] = byte_rev[l_mask_bytes[i]];
    }
    uint8_t un_pad = TWO * W - (size_ * TWO);
    rc_l_sym = (rc_l_sym >> un_pad) | (rc_h_sym << (W - un_pad));
    rc_l_mask = (rc_l_mask >> un_pad) | (rc_h_mask << (W - un_pad));
    rc_h_sym = rc_h_sym >> un_pad;
    rc_h_mask = rc_h_mask >> un_pad;

    return {high_bits_ & high_mask_, low_bits_ & low_mask_, high_mask_,
            low_mask_, size_};
  }

  template <class seq_t>
  bool match(const seq_t& seq) const {
    return (low_bits_ == (seq[0] & low_mask_)) &&
           (high_bits_ == (seq[1] & high_mask_));
  }

  uint8_t length() const { return size_; }

  template <class OS_T>
  OS_T& print(OS_T& out) const {
    if (size_ <= 32) {
      return mer(low_bits_, low_mask_, size_).print(out);
    }
    for (size_t i = (size_ - THIRTYTWO) * TWO - TWO; i < W; i -= TWO) {
      char c = internal::v_to_nuc[(high_bits_ >> i) & 0b11];
      if (((high_mask_ >> i) & 0b11) == ZERO) {
        c = 'N';
      }
      out << c;
    }
    for (size_t i = W - TWO; i < W; i -= TWO) {
      char c = internal::v_to_nuc[(low_bits_ >> i) & 0b11];
      if (((low_mask_ >> i) & 0b11) == ZERO) {
        c = 'N';
      }
      out << c;
    }
    return out;
  }
};

}  // namespace fkc
