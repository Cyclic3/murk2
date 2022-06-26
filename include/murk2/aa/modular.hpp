#pragma once

#include <murk2/aa/group.hpp>

namespace murk2::aa {
  class mod_add_group final : public cyclic_group<bigint>, public finite_magma<bigint> {
  public:
    bigint modulus;

  public:
    bigint op(bigint const& a, bigint const& b) const override;
    void op_mut(bigint& a, bigint const& b) const override;
    bigint op_iter(bigint const& a, bigint const& reps = 2) const override;
    void op_iter_mut(bigint& a, bigint const& reps = 2) const override;

    bigint identity() const noexcept override;
    bool is_identity(bigint const& i) const noexcept override;

    bigint invert(bigint const& a) const override;

    bigint generator() const override;

    bigint order() const override;

  public:
    mod_add_group(bigint modulus);
  };

  class mod_mul_base : public virtual monoid<bigint>, public finite_magma<bigint> {
  public:
    bigint modulus;

  public:
    bigint op(bigint const& a, bigint const& b) const override;
    void op_mut(bigint& a, bigint const& b) const override;

    bigint identity() const noexcept override;
    bool is_identity(bigint const& i) const noexcept override;

    bigint order() const override;

  public:
    mod_mul_base(bigint modulus);
  };

  class mod_mul_monoid : public mod_mul_base {
  public:
    bigint op_iter(bigint const& a, bigint const& reps = 2) const override;
    bigint op_iter(bigint&& a, bigint const& reps = 2) const override;
    void op_iter_mut(bigint& a, bigint const& reps = 2) const override;

    bool is_invertible(bigint const&) const override;
    std::optional<bigint> try_invert(bigint const&) const override;
    std::optional<bigint> try_invert(bigint&&) const override;
    bool try_invert_mut(bigint&) const override;

  public:
    mod_mul_monoid(bigint modulus);
  };

  class mod_mul_prime_group: public virtual cyclic_group<bigint>, public mod_mul_base  {
  private:
    std::optional<bigint> gen;

  public:
    bigint op_iter(bigint const& a, bigint const& reps = 2) const override;
    bigint op_iter(bigint&& a, bigint const& reps = 2) const override;
    void op_iter_mut(bigint& a, bigint const& reps = 2) const override;

    bigint invert(bigint const& a) const override;
    void invert_mut(bigint& a) const override;

    bigint generator() const override;

  public:
    mod_mul_prime_group(bigint modulus);
    mod_mul_prime_group(bigint modulus, bigint generator);
  };

  class mod_ring final : public finite_ring<bigint> {
  private:
    mod_add_group add_;
    mod_mul_monoid mul_;

  public:
    bigint const& get_modulus() const;
    void set_modulus(bigint n);

  public:
    c3lt::managed<const group<bigint>> add() const noexcept override;
    c3lt::managed<const monoid<bigint>> ring_mul() const noexcept override;

  public:
    mod_ring(bigint modulus);
  };

  class mod_prime_field final : public virtual finite_field<bigint> {
  private:
    mod_add_group add_;
    mod_mul_prime_group mul_;

  public:
    bigint const& get_modulus() const;
    void set_modulus(bigint n);

  public:
    c3lt::managed<const group<bigint>> add() const noexcept override;
    c3lt::managed<const group<bigint>> mul() const noexcept override;

    bigint order() const override final;
    bigint order_prime() const override final;
    bigint order_exponent() const override final;

  public:
    mod_prime_field(bigint modulus);
    mod_prime_field(bigint modulus, bigint mul_generator);
  };

  /// Puts the element with the range [0, modulus)
  void canonicalise_mod_mut(bigint& x, bigint const& modulus);
  /// Puts the element with the range [0, modulus)
  [[nodiscard]]
  bigint canonicalise_mod(bigint const& x, bigint const& modulus);
  /// Puts the element with the range [0, modulus)
  [[nodiscard]]
  bigint canonicalise_mod(bigint&& x, bigint const& modulus);

  /// Puts the element with the range [0, modulus)
  inline void canonicalise_mod_mut(ring_element<mod_ring, bigint>& x) {
    canonicalise_mod_mut(x.elem, x.context->get_modulus());
  }
  /// Puts the element with the range [0, modulus)
  [[nodiscard]]
  inline ring_element<mod_ring, bigint> canonicalise_mod(ring_element<mod_ring, bigint> const& x) {
    return {x.context, canonicalise_mod(x.elem, x.context->get_modulus())};
  }
  /// Puts the element with the range [0, modulus)
  [[nodiscard]]
  inline ring_element<mod_ring, bigint> canonicalise_mod(ring_element<mod_ring, bigint>&& x) {
    return {x.context, canonicalise_mod(std::move(x.elem), x.context->get_modulus())};
  }

  bigint find_primitive_root(bigint const& modulus);
}
