#pragma once

#include <murk2/aa/group.hpp>

#include <iostream>

#include <array>
#include <variant>

#include <cstddef>

namespace murk2::geo {
  /// Represents a finite dimensional vector
//  template<typename Ring, size_t Dimension, bool is_field = true>
//  class vector;

//  template<typename Ring, size_t Dimension>
//  class vector<Ring, Dimension, false> : public std::array<typename Ring::elem_t, Dimension> {
//  public:
//    using ring_t = Ring;
//    constexpr static size_t dimension = Dimension;

//  protected:
//    c3lt::managed<const Ring> ring;

//  public:
//    inline Ring operator+(vector const& other) const {
//      std::array<typename Ring::elem_t, Dimension> arr;
//      std::transform(this->begin(), this->end(), other.begin(), arr.begin(), [this](auto& a, auto& b) { return ring->add()->op(a, b); });
//      return arr;
//    }
//    inline Ring operator+=(vector const& other) {
//      std::transform(this->begin(), this->end(), other.begin(), this->begin());
//    }
//    inline bool operator==(const vector& other) const { return std::equal(this->begin(), this->end(), other.begin()); }
//    inline bool is_zero_vector() const { return std::all_of(this->begin(), this->end(), [this](auto& i) { return ring->add()->is_identity(i); }); }

//  public:
//    vector(c3lt::managed<const Ring> ring_, std::array<typename Ring::elem_t, Dimension> coords) : std::array<typename Ring::elem_t, Dimension>{coords}, ring{ring_} {}
//  };

  template<typename Ring, size_t Dimension>
  class vector/*<Ring, Dimension, true>*/ : public std::array<typename Ring::elem_t, Dimension> {
  public:
    using ring_t = Ring;
    constexpr static size_t dimension = Dimension;

  protected:
    c3lt::safe_ptr <const Ring> ring;

  public:
    constexpr c3lt::safe_ptr <const Ring> get_ring() const noexcept { return ring; }
    bool is_linearly_dependent(vector const& other) {
      std::optional<bigint> last_ratio;
      for (size_t i = 0; i < Dimension + 1; ++i) {
        auto& our_val = this->operator[](i);
        auto& other_val = other[i];
        if (ring->add()->is_identity(our_val)) {
          if (!ring->add()->is_identity(other_val))
            return false;
          else
            continue;
        }
        else if (!last_ratio) {
          last_ratio = other_val / our_val;
          continue;
        }
        else if (other_val / our_val != *last_ratio)
          return false;
      }
      return true;
    }

  public:
    inline Ring operator+(vector const& other) const {
      std::array<typename Ring::elem_t, Dimension> arr;
      std::transform(this->begin(), this->end(), other.begin(), arr.begin(), [this](auto& a, auto& b) { return ring->add()->op(a, b); });
      std::cout << *this << std::endl;
      return arr;
    }
    inline Ring operator+=(vector const& other) {
      std::transform(this->begin(), this->end(), other.begin(), this->begin());
    }
    inline bool operator==(const vector& other) const { return std::equal(this->begin(), this->end(), other.begin()); }
    inline virtual bool is_zero_vector() const { return std::all_of(this->begin(), this->end(), [this](auto& i) { return ring->add()->is_identity(i); }); }

  public:
    vector(c3lt::safe_ptr <const Ring> ring_, std::array<typename Ring::elem_t, Dimension> coords) : std::array<typename Ring::elem_t, Dimension>{coords}, ring{ring_} {}
  };

  template<typename Ring, size_t Dimension>
  vector(c3lt::safe_ptr <const Ring>, std::array<typename Ring::elem_t, Dimension>) -> vector<Ring, Dimension>;

  template<typename Ring, size_t Dimension>
  std::ostream& operator<<(std::ostream& os, vector<Ring, Dimension> const& coord) {
    os << '(';
    for (size_t i = 0; i < Dimension - 1; ++i) {
      os << coord.at(i) << ", ";
    }
    os << coord.back();
    return os << ')';
  }

  template<typename Ring>
  std::ostream& operator<<(std::ostream& os, vector<Ring, 0> const&) { return os << "()"; }
}
