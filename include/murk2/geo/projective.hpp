#pragma once

#include <murk2/geo/vector.hpp>

#include <concepts>
#include <iostream>

namespace murk2::geo {
//  struct point_at_infinity_t {};
//  constexpr point_at_infinity_t point_at_infinity;

  template<typename Field, size_t Dimension>
  class affinisation_vector;

  template<typename Field, size_t Dimension>
  class point_at_infinity : public geo::vector<Field, Dimension> {
  public:
    point_at_infinity(c3lt::managed<const Field> field_, std::array<typename Field::elem_t, Dimension> coords) : geo::vector<Field, Dimension>{field_, coords} {}
    point_at_infinity(point_at_infinity const&) = default;
  };

  template<typename Field, size_t Dimension>
  std::ostream& operator<<(std::ostream& os, point_at_infinity<Field, Dimension> const& coord) {
    return os << static_cast<geo::vector<Field, Dimension> const&>(coord) << "@inf";
  }

  template<typename Field, size_t Dimension>
  class projective_coordinate : public geo::vector<Field, Dimension + 1> {
  public:
    using field_t = Field;
    constexpr static size_t dimension = Dimension;

  private:
    c3lt::managed<const field_t> field;

  public:
    affinisation_vector<Field, Dimension> affinise(size_t with_respect_to) {
      if (this->at(with_respect_to) == 0) {
        std::array<Field, Dimension> coords;
        auto iter = std::copy(this->begin(), this->begin() + with_respect_to, coords.begin());
        std::copy(this->begin() + with_respect_to + 1, this->end(), iter);
        return point_at_infinity{field, coords};
      }
      else {
        std::array<Field, Dimension> coords;
        for (size_t i = 0; i < with_respect_to; ++i)
          coords[i] = this->operator[](i) / this->operator[](with_respect_to);
        for (size_t i = with_respect_to + 1; i <= Dimension; ++i)
          coords[i - 1] = this->operator[](i) / this->operator[](with_respect_to);
        return geo::vector{field, coords};
      }
    }

  public:
    bool operator==(projective_coordinate const& other) {
      return is_linearly_dependent(other);
    }

  public:
    projective_coordinate(c3lt::managed<const field_t> field_, std::array<typename Field::elem_t, Dimension + 1> coords) : std::array<typename Field::elem_t, Dimension + 1>{coords}, field{field_} {}
  };

  template<typename Field, size_t Dimension>
  std::ostream& operator<<(std::ostream& os, projective_coordinate<Field, Dimension> const& coord) {
    os << '(';
    for (size_t i = 0; i < Dimension; ++i) {
      os << coord.at(i) << ':';
    }
    os << coord.back();
    return os << ')';
  }

  // Represents the (potentially infinite) coordinates in the affinisation of a projective space
  template<typename Field, size_t Dimension>
  class affinisation_vector : public std::variant<geo::vector<Field, Dimension>, point_at_infinity<Field, Dimension>> {
  public:
    inline bool is_point_at_infinity() const noexcept { return std::holds_alternative<point_at_infinity<Field, Dimension>>(*this); }
    inline geo::vector<Field, Dimension>& finite_point() { return std::get<geo::vector<Field, Dimension>>(*this); }
    inline geo::vector<Field, Dimension> const& finite_point() const { return std::get<geo::vector<Field, Dimension>>(*this); }
    inline geo::vector<Field, Dimension>& infinite_point() { return std::get<point_at_infinity<Field, Dimension>>(*this); }
    inline geo::vector<Field, Dimension> const& infinite_point() const { return std::get<point_at_infinity<Field, Dimension>>(*this); }

    projective_coordinate<Field, Dimension> projectivise(size_t with_respect_to) {
      std::array<Field, Dimension + 1> coords;
      if (with_respect_to >= coords.size())
        throw std::out_of_range{"Attempted to projectivise with respect to out of range coordinate"};

      std::visit([this, with_respect_to, &coords](auto const& x) {
        using T = std::remove_cvref_t<decltype(x)>;
        if (with_respect_to == x.size())
          std::copy(x.begin(), x.end(), coords.begin());
        else {
          auto iter = std::copy(x.begin(), x.begin() + with_respect_to, coords.begin());
          std::copy(x.begin() + with_respect_to + 1, x.end(), ++iter);
        }
        if constexpr (std::same_as<T, point_at_infinity<Field, Dimension>>)
          coords[with_respect_to] = this->field->add()->identity();
        else
          coords[with_respect_to] = this->field->mul()->identity();
      }, *this);

      return {this->field, coords};
    }

  public:
    affinisation_vector(c3lt::managed<const Field> field, std::array<typename Field::elem_t, Dimension> coords): std::variant<geo::vector<Field, Dimension>, point_at_infinity<Field, Dimension>>{field, coords} {}
    affinisation_vector(geo::vector<Field, Dimension> vec): std::variant<geo::vector<Field, Dimension>, point_at_infinity<Field, Dimension>>{vec} {}
    affinisation_vector(point_at_infinity<Field, Dimension> point): std::variant<geo::vector<Field, Dimension>, point_at_infinity<Field, Dimension>>{point} {}
  };

  template<typename Field, size_t Dimension>
  affinisation_vector(geo::vector<Field, Dimension> vec) -> affinisation_vector<Field, Dimension>;
  template<typename Field, size_t Dimension>
  affinisation_vector(c3lt::managed<const Field>, std::array<typename Field::elem_t, Dimension>) -> affinisation_vector<Field, Dimension>;

  template<typename Field, size_t Dimension>
  std::ostream& operator<<(std::ostream& os, affinisation_vector<Field, Dimension> const& coord) {
    std::visit([&os](auto& x) { os << x; }, coord);
    return os;
  }
}
