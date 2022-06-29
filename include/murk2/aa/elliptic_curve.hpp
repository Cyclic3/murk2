#pragma once

#include <murk2/aa/group.hpp>
#include <murk2/geo/projective.hpp>
#include <murk2/common/err.hpp>

#include <array>

namespace murk2::aa {
  enum class elliptic_group_type {
    Generic,
    Finite
  };

  bigint compute_elliptic_curve_group_order(bigint A, bigint B, bigint order);
  template<typename GroundRing>
  bool check_elliptic_curve_nonsingular(ring_element<GroundRing> const& A, ring_element<GroundRing> const& B) {
    return !A.context->add()->is_identity(((((A^3)&4) + (B*B)&27)&16).elem);
  }


  template<typename GroundField, elliptic_group_type Type = elliptic_group_type::Generic>
  class elliptic_curve_group;

  template<typename GroundField>
  elliptic_curve_group(c3lt::managed<const GroundField> ground_field_, typename GroundField::elem_t a, typename GroundField::elem_t b) -> elliptic_curve_group<std::enable_if_t<!std::derived_from<GroundField, finite_ring<typename GroundField::elem_t>>, GroundField>>;
  template<typename GroundField>
  elliptic_curve_group(c3lt::managed<GroundField> ground_field_, typename GroundField::elem_t a, typename GroundField::elem_t b) -> elliptic_curve_group<std::enable_if_t<!std::derived_from<GroundField, finite_ring<typename GroundField::elem_t>>, GroundField>>;
  template<typename GroundField, int Dummy1 = 0>
  elliptic_curve_group(c3lt::managed<const GroundField> ground_field_, typename GroundField::elem_t a, typename GroundField::elem_t b) -> elliptic_curve_group<std::enable_if_t<std::derived_from<GroundField, finite_field<typename GroundField::elem_t>>, GroundField>, elliptic_group_type::Finite>;
  template<typename GroundField, int Dummy1 = 0>
  elliptic_curve_group(c3lt::managed<GroundField> ground_field_, typename GroundField::elem_t a, typename GroundField::elem_t b) -> elliptic_curve_group<std::enable_if_t<std::derived_from<GroundField, finite_field<typename GroundField::elem_t>>, GroundField>, elliptic_group_type::Finite>;

  template<typename GroundField>
  class elliptic_curve_group<GroundField, elliptic_group_type::Generic>: public virtual group<geo::affinisation_coordinate<GroundField, 2>> {
  public:
    constexpr static elliptic_group_type type = elliptic_group_type::Generic;

  private:
    c3lt::managed<const GroundField> ground_field;
    typename GroundField::elem_t A;
    typename GroundField::elem_t B;
    geo::affinisation_coordinate<GroundField, 2> id;

  public:
    constexpr typename GroundField::elem_t const& get_A() const noexcept { return A; }
    constexpr typename GroundField::elem_t const& get_B() const noexcept { return B; }
    constexpr GroundField const& get_ground_field() const noexcept { return *ground_field; }

  private:
    geo::affinisation_coordinate<GroundField, 2> op_internal(geo::vector<GroundField, 2> const& P, geo::vector<GroundField, 2> const& Q, typename GroundField::elem_t const& grad) const {
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
    typename geo::affinisation_coordinate<GroundField, 2> op(geo::affinisation_coordinate<GroundField, 2> const& a, geo::affinisation_coordinate<GroundField, 2> const& b) const override {
      return std::visit([this](auto& a, auto& b) -> geo::affinisation_coordinate<GroundField, 2> {
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

    geo::affinisation_coordinate<GroundField, 2> identity() const noexcept override { return id; }
    bool is_identity(geo::affinisation_coordinate<GroundField, 2> const& i) const noexcept override { return i.is_point_at_infinity(); }

    geo::affinisation_coordinate<GroundField, 2> invert(geo::affinisation_coordinate<GroundField, 2> const& a) const override {
      return std::visit([this](auto& a) -> geo::affinisation_coordinate<GroundField, 2> {
        if constexpr (std::same_as<std::remove_cvref_t<decltype(a)>, geo::point_at_infinity<GroundField, 2>>)
          return a;
        else
          return geo::vector{ground_field, std::array{a[0], ground_field->add()->invert(a[1])}};
      }, a);
    }

    inline auto element(typename GroundField::elem_t x, typename GroundField::elem_t y) const {
      return element(geo::vector{ground_field, std::array<decltype(x), 2>{std::move(x), std::move(y)}});
    }
    using group<geo::affinisation_coordinate<GroundField, 2>>::element;

    inline auto operator()(typename GroundField::elem_t x, typename GroundField::elem_t y) const { return element(std::move(x), std::move(y)); }
    using group<geo::affinisation_coordinate<GroundField, 2>>::operator();

  public:
    elliptic_curve_group(c3lt::managed<const GroundField> ground_field_, typename GroundField::elem_t a, typename GroundField::elem_t b) :
        ground_field{ground_field_},
        A{std::move(a)},
        B{std::move(b)},
        id{geo::point_at_infinity<GroundField, 2>(ground_field, std::array{ground_field->add()->identity(), ground_field->mul()->identity()})} {}
    virtual ~elliptic_curve_group() = default;
  };

  template<typename GroundField>
  class elliptic_curve_group<GroundField, elliptic_group_type::Finite> : public virtual elliptic_curve_group<GroundField, elliptic_group_type::Generic>, public virtual finite_monoid<typename elliptic_curve_group<GroundField>::elem_t> {
  public:
    constexpr static elliptic_group_type type = elliptic_group_type::Finite;

  public:
    bigint order() const override { return compute_elliptic_curve_group_order(this->get_A(), this->get_B(), this->get_ground_field().order()); }

  public:
    elliptic_curve_group(c3lt::managed<const GroundField> ground_field, typename GroundField::elem_t a, typename GroundField::elem_t b) : elliptic_curve_group<GroundField, elliptic_group_type::Generic>{ground_field, std::move(a), std::move(b)} {}
    virtual ~elliptic_curve_group() = default;
  };


  // SOURCE: https://arxiv.org/pdf/2010.15543.pdf
  template<typename GroundRing>
  class elliptic_curve_proj_group : public group<geo::projective_coordinate<GroundRing, 2>> {
  private:
    using ground_elem_t = typename GroundRing::elem_t;

  private:
    c3lt::managed<const GroundRing> ground_ring;
    geo::projective_coordinate<GroundRing, 2> id;
    ground_elem_t A;
    ground_elem_t B;

  public:
    bool check_point(geo::projective_coordinate<GroundRing, 2> const& P) const {
      auto& X = P[0];
      auto& Y = P[1];
      auto& Z = P[2];

      auto z2 = ground_ring->ring_mul()->op(Z, Z); // z^2
      auto x3_part = ground_ring->ring_mul()->op_iter(X, 3); // x^3
      auto xz2_part = ground_ring->ring_mul()->op(A, ground_ring->ring_mul()->op(X, z2)); // A * X.Z^2
      auto z3_part = ground_ring->ring_mul()->op(B, ground_ring->ring_mul()->op(z2, Z)); // B * Z^3
      auto y2z_part = ground_ring->ring_mul()->op(ground_ring->ring_mul()->op(Y, Y), Z); // Y^2.Z

      auto res = ground_ring->add()->op(ground_ring->add()->op(ground_ring->add()->op(x3_part, xz2_part), z3_part), ground_ring->add()->invert(y2z_part));

      return ground_ring->add()->is_identity(std::move(res));
    }

    geo::projective_coordinate<GroundRing, 2> canonicalise(geo::projective_coordinate<GroundRing, 2>&& elem) const override {
//      if (!check_point(elem))
//        throw std::invalid_argument{"Given point off curve"};
      return std::move(elem);
    }

    geo::projective_coordinate<GroundRing, 2> identity() const noexcept override { return id; }
    // If the last coord is zero, then it is the point at infinity (assuming that the point is actually on the curve)
    bool is_identity(geo::projective_coordinate<GroundRing, 2> const& P) const override {
      return !ground_ring->ring_mul()->is_invertible(P[0]) && !ground_ring->ring_mul()->is_invertible(P[2]);
    }

    // SOURCE: https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Standard_Projective_Coordinates
    geo::projective_coordinate<GroundRing, 2> op(geo::projective_coordinate<GroundRing, 2> const& P, geo::projective_coordinate<GroundRing, 2> const& Q) const override {
      if (is_identity(P))
        return Q;
      else if (is_identity(Q))
        return P;

      auto const& X1 = P[0];
      auto const& Y1 = P[1];
      auto const& Z1 = P[2];

      auto const& X2 = Q[0];
      auto const& Y2 = Q[1];
      auto const& Z2 = Q[2];

      ring_element<GroundRing> U1{ground_ring, ground_ring->ring_mul()->op(Y2, Z1)};
      ring_element<GroundRing> U2{ground_ring, ground_ring->ring_mul()->op(Y1, Z2)};
      ring_element<GroundRing> V1{ground_ring, ground_ring->ring_mul()->op(X2, Z1)};
      ring_element<GroundRing> V2{ground_ring, ground_ring->ring_mul()->op(X1, Z2)};

      if (V1 == V2) {
        if (U1 != U2)
          return identity();
        else {
          auto const& X = P[0];
          auto const& Y = P[1];
          auto const& Z = P[2];

          ring_element<GroundRing> W{ground_ring, A * ground_ring->ring_mul()->op(Z, Z) + ground_ring->add()->op_iter(ground_ring->ring_mul()->op(X, X), 3)}; // a*Z^2 + 3*X^2
          ring_element<GroundRing> S{ground_ring, ground_ring->ring_mul()->op(Y, Z)}; // Y*Z
          ring_element<GroundRing> B = S * ground_ring->ring_mul()->op(X, Y); // X.Y.S
          ring_element<GroundRing> H = (W*W) - (B&8); // W^2 - 8*B
          ring_element<GroundRing> X3 = (H * S)&2; // 2*H.S
          ring_element<GroundRing> Y3 = W*((B&4) - H) - (((S*S) * ground_ring->ring_mul()->op(Y, Y)) & 8); // W*(4*B - H) - 8*Y^2*S^2
          ring_element<GroundRing> Z3 = (S^3) & 8; // 8*S^3

//          std::cerr << this->A << ", " << this->B << ", " << ground_ring->order() << std::endl;

//          std::cerr << P << " * 2 == " << geo::projective_coordinate<GroundRing, 2>{ground_ring, {X3.elem, Y3.elem, Z3.elem}} << std::endl;

          return {ground_ring, {std::move(X3.elem), std::move(Y3.elem), std::move(Z3.elem)}};
        }
      }

      ring_element<GroundRing> U = U1 - U2;
      ring_element<GroundRing> V = V1 - V2;
      ring_element<GroundRing> W{ground_ring, Z1 * Z2};
      ring_element<GroundRing> V3 = (V*V) * V2;
      ring_element<GroundRing> A = (U*U)*W - (V^3) - (V3 & 2);

      ring_element<GroundRing> X3 = V * A;
      ring_element<GroundRing> Y3 = U*(V3 - A) - (U2 * (V^3));
      ring_element<GroundRing> Z3 = W * (V^3);

      return {ground_ring, {std::move(X3.elem), std::move(Y3.elem), std::move(Z3.elem)}};
    }

    geo::projective_coordinate<GroundRing, 2> invert(geo::projective_coordinate<GroundRing, 2> const& a) const override {
      if (is_identity(a))
        return a;
      else
        return geo::projective_coordinate<GroundRing, 2>{ground_ring, {a[0], ground_ring->add()->invert(a[1]), a[2]}};
    }

    inline auto element(typename GroundRing::elem_t x, typename GroundRing::elem_t y, typename GroundRing::elem_t z) const {
      return element(geo::vector{ground_ring, std::array<decltype(x), 3>{std::move(x), std::move(y), std::move(z)}});
    }
    using group<geo::projective_coordinate<GroundRing, 2>>::element;

    inline auto operator()(typename GroundRing::elem_t x, typename GroundRing::elem_t y, typename GroundRing::elem_t z) const { return element(std::move(x), std::move(y), std::move(z)); }
    using group<geo::projective_coordinate<GroundRing, 2>>::operator();

  public:
    elliptic_curve_proj_group(c3lt::managed<const GroundRing> ground_ring_, typename GroundRing::elem_t a, typename GroundRing::elem_t b) :
      ground_ring{ground_ring_},
      A{std::move(a)},
      B{std::move(b)},
      id{geo::projective_coordinate<GroundRing, 2>(ground_ring, std::array{ground_ring->add()->identity(), ground_ring->ring_mul()->identity(), ground_ring->add()->identity()})} {}
    virtual ~elliptic_curve_proj_group() = default;
  };
}
