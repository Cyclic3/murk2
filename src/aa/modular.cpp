#include <murk2/aa/modular.hpp>

#include <murk2/common/err.hpp>

#include <boost/integer/mod_inverse.hpp>

namespace murk2::aa {
  bigint mod_add_group::op(bigint const& a, bigint const& b) const { return (a + b) % modulus; }
  void mod_add_group::op_mut(bigint& a, bigint const& b) const { a += b; a %= modulus; }
  bigint mod_add_group::op_iter(bigint const& a, bigint const& reps) const { return (a * reps) % modulus; }
  void mod_add_group::op_iter_mut(bigint& a, bigint const& reps) const { a *= reps; a %= modulus; }
  bigint mod_add_group::identity() const noexcept { return 0; }
  bigint mod_add_group::invert(bigint const& a) const { return a.is_zero() ? a : modulus - a; }
  bool mod_add_group::is_identity(bigint const& a) const noexcept { return a.is_zero() || ((a % modulus) == 0); }
  bigint mod_add_group::generator() const { return 1; }
  bigint mod_add_group::order() const { return modulus; }
  mod_add_group::mod_add_group(bigint modulus_) : modulus{modulus_} {}

  // TODO: optimise multiplication
  bigint mod_mul_base::op(bigint const& a, bigint const& b) const { return (a * b) % modulus; }
  void mod_mul_base::op_mut(bigint& a, bigint const& b) const { a *= b; a %= modulus;}
  bigint mod_mul_base::identity() const noexcept { return 1; }
  bool mod_mul_base::is_identity(bigint const& a) const noexcept { return a == 1 || canonicalise_mod(a, modulus) == 1; }
  bigint mod_mul_base::order() const { return modulus - 1; }
  mod_mul_base::mod_mul_base(bigint modulus_) : modulus{std::move(modulus_)} {}

  bigint mod_mul_monoid::op_iter(bigint const& a, bigint const& reps) const {
    if (reps.sign() == -1) {
      if (auto new_a = try_invert(a))
        return boost::multiprecision::powm(*new_a, -reps, modulus);
      else
        throw missing_structure{"Tried to invert non-invertible element"};
    }
    return boost::multiprecision::powm(a, reps, modulus);
  }
  void mod_mul_monoid::op_iter_mut(bigint& a, bigint const& reps) const {
    if (reps.sign() == -1) {
      if (!try_invert_mut(a))
        throw missing_structure{"Tried to invert non-invertible element"};
      a = boost::multiprecision::powm(a, -reps, modulus);
    }
    else
      a = boost::multiprecision::powm(a, reps, modulus);
  }
  bool mod_mul_monoid::is_invertible(bigint const& a) const { return boost::multiprecision::gcd(a, modulus) == 1; }
  std::optional<bigint> mod_mul_monoid::try_invert(bigint const& a) const {
    std::optional<bigint> ret = boost::integer::mod_inverse(a, modulus);
    if (ret->is_zero())
      ret.reset();
    return ret;
  }
  bool mod_mul_monoid::try_invert_mut(bigint& a) const {
    auto res = boost::integer::mod_inverse(a, modulus);
    if (!a)
      return false;
    a = std::move(res);
    return true;
  }
  std::optional<bigint> mod_mul_monoid::try_invert(bigint&& a) const {
    a = boost::integer::mod_inverse(a, modulus);
    return a.is_zero() ? std::nullopt : std::optional{std::move(a)};
  }
  mod_mul_monoid::mod_mul_monoid(bigint modulus) : mod_mul_base{std::move(modulus)} {}
  bigint const& mod_ring::get_modulus() const { return add_.modulus; }
  void mod_ring::set_modulus(bigint n) {
    add_.modulus = n;
    mul_.modulus = std::move(n);
  }
  c3lt::managed<const group<bigint>> mod_ring::add() const noexcept { return c3lt::managed<const group<bigint>>{&add_}; }
  c3lt::managed<const monoid<bigint>> mod_ring::ring_mul() const noexcept { return c3lt::managed<const monoid<bigint>>{&mul_}; }

  mod_ring::mod_ring(bigint modulus) : add_{modulus}, mul_{std::move(modulus)} {}

  void canonicalise_mod_mut(bigint& x, bigint const& modulus) {
    x %= modulus;
    if (x.sign() == -1)
      x += modulus;
  }

  bigint canonicalise_mod(bigint const& x, bigint const& modulus) {
    bigint ret = x % modulus;
    if (ret.sign() == -1)
      ret += modulus;
    return ret;
  }

  bigint canonicalise_mod(bigint&& x, bigint const& modulus) {
    canonicalise_mod_mut(x, modulus);
    return std::move(x);
  }

  bigint mod_mul_prime_group::op_iter(bigint const& a, bigint const& reps) const {
    bigint ret;
    mpz_powm(ret.backend().data(), a.backend().data(), reps.backend().data(), modulus.backend().data());
    return ret;
//    if (reps.sign() == -1) {
//      auto new_a = invert(a);
//      return boost::multiprecision::powm(new_a, -reps, modulus);
//    }
//    return boost::multiprecision::powm(a, reps, modulus);
  }
  void mod_mul_prime_group::op_iter_mut(bigint& a, bigint const& reps) const {
    mpz_powm(a.backend().data(), a.backend().data(), reps.backend().data(), modulus.backend().data());
//    if (reps.sign() == -1) {
//      invert_mut(a);
//      a = boost::multiprecision::powm(a, -reps, modulus);
//    }
//    else
//      a = boost::multiprecision::powm(a, reps, modulus);
  }
  bigint mod_mul_prime_group::invert(bigint const& a) const {
    bigint ret;
    // Get rekt boost this is 10x quicker
    mpz_invert(ret.backend().data(), a.backend().data(), modulus.backend().data());
    return ret;
  }
  void mod_mul_prime_group::invert_mut(bigint& a) const {
    // Get rekt boost this is 10x quicker
    mpz_invert(a.backend().data(), a.backend().data(), modulus.backend().data());
  }
  bigint mod_mul_prime_group::generator() const { return gen ? *gen : find_primitive_root(modulus); }
  mod_mul_prime_group::mod_mul_prime_group(bigint modulus) : mod_mul_base{std::move(modulus)} {}
  mod_mul_prime_group::mod_mul_prime_group(bigint modulus, bigint mul_generator) : mod_mul_base{std::move(modulus)}, gen{std::move(mul_generator)} {}

  c3lt::managed<const group<bigint>> mod_prime_field::add() const noexcept { return c3lt::managed<const group<bigint>>{&add_}; }
  c3lt::managed<const group<bigint>> mod_prime_field::mul() const noexcept { return c3lt::managed{&mul_}; }
  bigint const& mod_prime_field::get_modulus() const { return add_.modulus; }
  bigint mod_prime_field::order_prime() const { return get_modulus(); }
  bigint mod_prime_field::order_exponent() const { return 1; }
  bigint mod_prime_field::order() const { return get_modulus(); }
  mod_prime_field::mod_prime_field(bigint modulus) : add_{modulus}, mul_{std::move(modulus)} {}
  mod_prime_field::mod_prime_field(bigint modulus, bigint mul_generator) : add_{modulus}, mul_{std::move(modulus), std::move(mul_generator)} {}

  bigint find_primitive_root(bigint const& modulus) MURK2_UNIMPLEMENTED;
}
