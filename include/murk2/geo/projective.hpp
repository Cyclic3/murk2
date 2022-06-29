#pragma once

#include <murk2/geo/vector.hpp>

#include <concepts>
#include <iostream>

namespace murk2::geo {
//  struct point_at_infinity_t {};
//  constexpr point_at_infinity_t point_at_infinity;

  template<typename Ring, size_t Dimension>
  class affinisation_coordinate;

  template<typename Ring, size_t Dimension>
  class point_at_infinity : public geo::vector<Ring, Dimension> {
  public:
    point_at_infinity(c3lt::managed<const Ring> ring_, std::array<typename Ring::elem_t, Dimension> coords) : geo::vector<Ring, Dimension>{ring_, coords} {}
    point_at_infinity(point_at_infinity const&) = default;
  };

  template<typename Ring, size_t Dimension>
  std::ostream& operator<<(std::ostream& os, point_at_infinity<Ring, Dimension> const& coord) {
    return os << static_cast<geo::vector<Ring, Dimension> const&>(coord) << "@inf";
  }

  template<typename Ring, size_t Dimension>
  class projective_coordinate : public geo::vector<Ring, Dimension + 1> {
  public:
    using ring_t = Ring;
    constexpr static size_t dimension = Dimension;

  public:
    affinisation_coordinate<Ring, Dimension> affinise(size_t with_respect_to) {
      std::array<typename Ring::elem_t, Dimension> coords;

      if (this->ring->add()->is_identity(this->at(with_respect_to))) {
        auto iter = std::copy(this->begin(), this->begin() + with_respect_to, coords.begin());
        std::copy(this->begin() + with_respect_to + 1, this->end(), iter);
        return point_at_infinity{this->ring, coords};
      }
      else {
        auto mul = this->ring->ring_mul();
        auto inverse_opt = mul->try_invert(this->operator[](with_respect_to));
        if (!inverse_opt)
          throw aa::missing_structure{"Tried to affinise with respect to non-invertible element"};

        for (size_t i = 0; i < with_respect_to; ++i)
          coords[i] = mul->op(this->operator[](i), *inverse_opt);
        for (size_t i = with_respect_to + 1; i <= Dimension; ++i)
          coords[i - 1] = mul->op(this->operator[](i), *inverse_opt);
        return geo::vector{this->ring, coords};
      }
    }

  public:
    bool operator==(projective_coordinate const& other) {
      return this->is_linearly_dependent(other);
    }

  public:
    projective_coordinate(c3lt::managed<const ring_t> ring_, std::array<typename Ring::elem_t, Dimension + 1> coords) : geo::vector<Ring, Dimension + 1>{ring_, std::move(coords)} {
      if (this->is_zero_vector())
        throw std::invalid_argument{"Zero vector given as projective coordinate"};
    }
    projective_coordinate(geo::vector<Ring, Dimension + 1> vec) : geo::vector<Ring, Dimension + 1>{std::move(vec)} {
      if (this->is_zero_vector())
        throw std::invalid_argument{"Zero vector given as projective coordinate"};
    }
  };

  template<typename Ring, size_t Dimension>
  std::ostream& operator<<(std::ostream& os, projective_coordinate<Ring, Dimension> const& coord) {
    os << '(';
    for (size_t i = 0; i < Dimension; ++i) {
      os << coord.at(i) << ':';
    }
    os << coord.back();
    return os << ')';
  }

  // Represents the (potentially infinite) coordinates in the affinisation of a projective space
  template<typename Ring, size_t Dimension>
  class affinisation_coordinate : public std::variant<geo::vector<Ring, Dimension>, point_at_infinity<Ring, Dimension>> {
  public:
    inline bool is_point_at_infinity() const noexcept { return std::holds_alternative<point_at_infinity<Ring, Dimension>>(*this); }
    inline geo::vector<Ring, Dimension>& finite_point() { return std::get<geo::vector<Ring, Dimension>>(*this); }
    inline geo::vector<Ring, Dimension> const& finite_point() const { return std::get<geo::vector<Ring, Dimension>>(*this); }
    inline geo::vector<Ring, Dimension>& infinite_point() { return std::get<point_at_infinity<Ring, Dimension>>(*this); }
    inline geo::vector<Ring, Dimension> const& infinite_point() const { return std::get<point_at_infinity<Ring, Dimension>>(*this); }

    projective_coordinate<Ring, Dimension> projectivise(size_t with_respect_to) {
      std::array<Ring, Dimension + 1> coords;
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
        if constexpr (std::same_as<T, point_at_infinity<Ring, Dimension>>)
          coords[with_respect_to] = this->ring->add()->identity();
        else
          coords[with_respect_to] = this->ring->mul()->identity();
      }, *this);

      return {this->ring, coords};
    }

  public:
    affinisation_coordinate(c3lt::managed<const Ring> ring, std::array<typename Ring::elem_t, Dimension> coords): std::variant<geo::vector<Ring, Dimension>, point_at_infinity<Ring, Dimension>>{ring, coords} {}
    affinisation_coordinate(geo::vector<Ring, Dimension> vec): std::variant<geo::vector<Ring, Dimension>, point_at_infinity<Ring, Dimension>>{vec} {}
    affinisation_coordinate(point_at_infinity<Ring, Dimension> point): std::variant<geo::vector<Ring, Dimension>, point_at_infinity<Ring, Dimension>>{point} {}
  };

  template<typename Ring, size_t Dimension>
  affinisation_coordinate(geo::vector<Ring, Dimension> vec) -> affinisation_coordinate<Ring, Dimension>;
  template<typename Ring, size_t Dimension>
  affinisation_coordinate(c3lt::managed<const Ring>, std::array<typename Ring::elem_t, Dimension>) -> affinisation_coordinate<Ring, Dimension>;

  template<typename Ring, size_t Dimension>
  std::ostream& operator<<(std::ostream& os, affinisation_coordinate<Ring, Dimension> const& coord) {
    std::visit([&os](auto& x) { os << x; }, coord);
    return os;
  }
}
