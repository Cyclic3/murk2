#pragma once

#include <murk2/aa/group.hpp>
#include <murk2/geo/projective.hpp>
#include <murk2/common/err.hpp>

#include <array>

namespace murk2::aa {
  enum class elliptic_group_type {
    Generic,
    Finite,
    FiniteCyclic
  };

  bigint compute_elliptic_curve_group_order(bigint A, bigint B, bigint order);


  template<typename GroundField, elliptic_group_type Type = elliptic_group_type::Generic>
  class elliptic_curve_group;

  template<typename GroundField>
  elliptic_curve_group(c3lt::managed<const GroundField> ground_field_, typename GroundField::elem_t a, typename GroundField::elem_t b) -> elliptic_curve_group<std::enable_if_t<!std::derived_from<GroundField, finite_ring<typename GroundField::elem_t>>, GroundField>>;
  template<typename GroundField, int Dummy1 = 0>
  elliptic_curve_group(c3lt::managed<const GroundField> ground_field_, typename GroundField::elem_t a, typename GroundField::elem_t b) -> elliptic_curve_group<std::enable_if_t<std::derived_from<GroundField, finite_field<typename GroundField::elem_t>>, GroundField>, elliptic_group_type::Finite>;

  template<typename GroundField>
  class elliptic_curve_group<GroundField, elliptic_group_type::Generic>: public virtual group<geo::affinisation_vector<GroundField, 2>> {
  public:
    constexpr static elliptic_group_type type = elliptic_group_type::Generic;

  private:
    c3lt::managed<const GroundField> ground_field;
    typename GroundField::elem_t A;
    typename GroundField::elem_t B;
//    geo::affinisation_vector<GroundField, 2> gen;
    geo::affinisation_vector<GroundField, 2> id;

  public:
    constexpr typename GroundField::elem_t const& get_A() const noexcept { return A; }
    constexpr typename GroundField::elem_t const& get_B() const noexcept { return B; }
    constexpr GroundField const& get_ground_field() const noexcept { return *ground_field; }

  private:
    geo::affinisation_vector<GroundField, 2> op_internal(geo::vector<GroundField, 2> const& P, geo::vector<GroundField, 2> const& Q, typename GroundField::elem_t const& grad) const {
      std::array<typename GroundField::elem_t, 2> coords;
      coords[0] = ground_field->mul()->op_iter(grad, 2); // grad^2
      ground_field->add()->op_mut(coords[0], ground_field->add()->invert(ground_field->add()->op(P[0], Q[0]))); // grad^2 + -(x_P + x_Q)

      coords[1] = ground_field->add()->invert(P[0]); // -x_P
      ground_field->add()->op_mut(coords[1], coords[0]); // x_R - x_P
      ground_field->mul()->op_mut(coords[1], grad); // s(x_R - x_P)
      ground_field->add()->op_mut(coords[1], P[1]); // y_P + s(x_R - x_P)
      ground_field->add()->invert_mut(coords[1]); // -(y_P + s(x_R - x_P))

      return geo::vector<GroundField, 2>{ground_field, coords};
    }

  public:
    typename geo::affinisation_vector<GroundField, 2> op(geo::affinisation_vector<GroundField, 2> const& a, geo::affinisation_vector<GroundField, 2> const& b) const override {
      return std::visit([this](auto& a, auto& b) -> geo::affinisation_vector<GroundField, 2> {
        // Handle the identity
        if constexpr (std::same_as<std::remove_cvref_t<decltype(a)>, geo::point_at_infinity<GroundField, 2>>)
          return b;
        else if constexpr (std::same_as<std::remove_cvref_t<decltype(b)>, geo::point_at_infinity<GroundField, 2>>)
          return a;
        else {
          if (a[0] == b[0]) {
            // Check if adding to own inverse
            if (ground_field->add()->is_inverse(a[1], b[1]))
              return this->identity();

            // We therefore have that a[1] == b[1] != 0, and so this is a simple duplication
            auto denom = ground_field->add()->op(a[1], a[1]); // y^2
            ground_field->mul()->invert_mut(denom); // y^(-2)

            auto grad = ground_field->mul()->op_iter(a[0], 2); // x^2
            ground_field->add()->op_iter_mut(grad, 3); // 3 * (x^2)
            ground_field->add()->op_mut(grad, this->A); // 3 * (x^2) + A
            ground_field->mul()->op_mut(grad, denom); // (3 * (x^2) + A) / y^2

            return op_internal(a, b, grad);
          }
          else {
            auto grad = ground_field->add()->op(a[1], ground_field->add()->invert(b[1])); // y_P - y_Q
            auto grad_denom = ground_field->add()->op(a[0], ground_field->add()->invert(b[0])); // x_P - x_Q
            ground_field->mul()->op_mut(grad, ground_field->mul()->invert(grad_denom));
            return op_internal(a, b, grad);
          }
        }
      }, a, b);
    }

    geo::affinisation_vector<GroundField, 2> identity() const noexcept override { return id; }
    bool is_identity(geo::affinisation_vector<GroundField, 2> const& i) const noexcept override { return i.is_point_at_infinity(); }

    geo::affinisation_vector<GroundField, 2> invert(geo::affinisation_vector<GroundField, 2> const& a) const override {
      return std::visit([this](auto& a) -> geo::affinisation_vector<GroundField, 2> {
        if constexpr (std::same_as<std::remove_cvref_t<decltype(a)>, geo::point_at_infinity<GroundField, 2>>)
          return a;
        else
          return geo::vector{ground_field, std::array{a[0], ground_field->add()->invert(a[1])}};
      }, a);
    }

    inline auto element(typename GroundField::elem_t x, typename GroundField::elem_t y) const {
      return element(geo::vector{ground_field, std::array<decltype(x), 2>{std::move(x), std::move(y)}});
    }
    using group<geo::affinisation_vector<GroundField, 2>>::element;

    inline auto operator()(typename GroundField::elem_t x, typename GroundField::elem_t y) const { return element(std::move(x), std::move(y)); }
    using group<geo::affinisation_vector<GroundField, 2>>::operator();

//    geo::affinisation_vector<GroundField, 2> generator() const override;

//    bigint order() const override;

  public:
    elliptic_curve_group(c3lt::managed<const GroundField> ground_field_, typename GroundField::elem_t a, typename GroundField::elem_t b) :
        ground_field{ground_field_},
        A{std::move(a)},
        B{std::move(b)},
        id{geo::point_at_infinity<GroundField, 2>(ground_field_, std::array{ground_field_->add()->identity(), ground_field_->mul()->identity()})} {}
    virtual ~elliptic_curve_group() = default;
  };

  template<typename GroundField>
  class elliptic_curve_group<GroundField, elliptic_group_type::Finite> : public virtual elliptic_curve_group<GroundField, elliptic_group_type::Generic>, public virtual finite_magma<typename elliptic_curve_group<GroundField>::elem_t> {
  public:
    constexpr static elliptic_group_type type = elliptic_group_type::Finite;

  public:
    bigint order() const override { return compute_elliptic_curve_group_order(this->get_A(), this->get_B(), this->get_ground_field().order()); }

  public:
    elliptic_curve_group(c3lt::managed<const GroundField> ground_field, typename GroundField::elem_t a, typename GroundField::elem_t b) : elliptic_curve_group<GroundField, elliptic_group_type::Generic>{ground_field, std::move(a), std::move(b)} {}
    virtual ~elliptic_curve_group() = default;
  };

  // TODO: work out how to fix the inheritence
  template<typename GroundField>
  class elliptic_curve_group<GroundField, elliptic_group_type::FiniteCyclic> : public elliptic_curve_group<GroundField, elliptic_group_type::Generic>, public virtual finite_magma<typename elliptic_curve_group<GroundField>::elem_t> , public virtual cyclic_group<typename elliptic_curve_group<GroundField>::elem_t> {
  public:
    constexpr static elliptic_group_type type = elliptic_group_type::FiniteCyclic;

  private:
    std::optional<geo::affinisation_vector<GroundField, 2>> gen;

  private:
    geo::affinisation_vector<GroundField, 2> find_generator() const MURK2_UNIMPLEMENTED;

  public:
    bigint order() const override { return compute_elliptic_curve_group_order(this->get_A(), this->get_B(), this->get_ground_field().order()); }
    geo::affinisation_vector<GroundField, 2> generator() const override { return gen ? *gen : find_generator(); }

  public:
    elliptic_curve_group(c3lt::managed<const GroundField> ground_field, typename GroundField::elem_t a, typename GroundField::elem_t b) :
        elliptic_curve_group<GroundField, elliptic_group_type::Generic>{ground_field, std::move(a), std::move(b)} {}
    elliptic_curve_group(c3lt::managed<const GroundField> ground_field, typename GroundField::elem_t a, typename GroundField::elem_t b, geo::affinisation_vector<GroundField, 2> gen_) :
        elliptic_curve_group<GroundField, elliptic_group_type::Generic>{ground_field, std::move(a), std::move(b)},
        gen{gen_} {}
    elliptic_curve_group(c3lt::managed<const GroundField> ground_field, typename GroundField::elem_t a, typename GroundField::elem_t b, std::array<typename GroundField::elem_t, 2> gen_) :
        elliptic_curve_group<GroundField, elliptic_group_type::Generic>{ground_field, std::move(a), std::move(b)},
        gen{geo::vector<GroundField, 2>{ground_field, gen_}} {}
    virtual ~elliptic_curve_group() = default;
  };
}
