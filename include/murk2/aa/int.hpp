#pragma once

#include <murk2/aa/group.hpp>

namespace murk2::aa {
  struct int_add_group final : public cyclic_group<bigint> {
    bigint op(bigint const& a, bigint const& b) const override { return a + b; }
    void op_mut(bigint& a, bigint const& b) const override { a += b; }
    bigint op_iter(bigint const& a, bigint const& reps = 2) const override { return a * static_cast<bigint>(reps); }
    void op_iter_mut(bigint& a, bigint const& reps = 2) const override { a *= static_cast<bigint>(reps); }

    bigint identity() const noexcept override { return 0; }

    bigint invert(bigint const& a) const override { return -a; }

    bigint generator() const override { return 1; }

    int_add_group() = default;
    virtual ~int_add_group() = default;
  };

  struct int_mul_monoid final : public monoid<bigint> {
    bigint op(bigint const& a, bigint const& b) const override { return a * b; }
    void op_mut(bigint& a, bigint const& b) const override { a *= b; }

    bigint identity() const noexcept override { return 1; }

    int_mul_monoid() = default;
    virtual ~int_mul_monoid() = default;
  };

  class int_ring : public ring<bigint> {
  private:
    static int_add_group add_;
    static int_mul_monoid mul_;

  public:
    c3lt::safe_ptr<const group<bigint>> add() const noexcept override { return c3lt::managed(&add_); }
    c3lt::safe_ptr<const monoid<bigint>> ring_mul() const noexcept override { return c3lt::managed(&mul_); }
  };
}
