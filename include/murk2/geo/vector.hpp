#pragma once

#include <murk2/aa/group.hpp>

#include <array>
#include <variant>

#include <cstddef>

namespace murk2::geo {
  /// Represents a finite dimensional vector
  template<typename Field, size_t Dimension>
  class vector : public std::array<typename Field::elem_t, Dimension> {
  public:
    using field_t = Field;
    constexpr static size_t dimension = Dimension;

  private:
    c3lt::managed<const Field> field;

  public:
    bool is_linearly_dependent(vector const& other) {
      std::optional<bigint> last_ratio;
      for (size_t i = 0; i < Dimension + 1; ++i) {
        auto& our_val = this->operator[](i);
        auto& other_val = other[i];
        if (field->add()->is_identity(our_val)) {
          if (!field->add()->is_identity(other_val))
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
    inline Field operator+(vector const& other) const {
      std::array<typename Field::elem_t, Dimension> arr;
      std::transform(this->begin(), this->end(), other.begin(), arr.begin(), [this](auto& a, auto& b) { return field->add()->op(a, b); });
      return arr;
    }
    inline Field operator+=(vector const& other) {
      std::transform(this->begin(), this->end(), other.begin(), this->begin());
    }
    inline bool operator==(const vector& other) const { return std::equal(this->begin(), this->end(), other.begin()); }

  public:
    vector(c3lt::managed<const Field> field_, std::array<typename Field::elem_t, Dimension> coords) : std::array<typename Field::elem_t, Dimension>{coords}, field{field_} {}
  };

  template<typename Field, size_t Dimension>
  vector(c3lt::managed<const Field>, std::array<typename Field::elem_t, Dimension>) -> vector<Field, Dimension>;

  template<typename Field, size_t Dimension>
  std::ostream& operator<<(std::ostream& os, vector<Field, Dimension> const& coord) {
    os << '(';
    for (size_t i = 0; i < Dimension - 1; ++i) {
      os << coord.at(i) << ", ";
    }
    os << coord.back();
    return os << ')';
  }

  template<typename Field>
  std::ostream& operator<<(std::ostream& os, vector<Field, 0> const&) { return os << "()"; }
}
