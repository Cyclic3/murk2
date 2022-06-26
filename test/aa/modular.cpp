#include <gtest/gtest.h>

#include <murk2/aa/modular.hpp>

TEST(murk2aa, ModSignFix) {
  static murk2::aa::mod_ring ring(21);

  murk2::aa::ring_element minus_one{c3lt::managed(&ring), -1};
  murk2::aa::ring_element twenty{c3lt::managed(&ring), 20};

  EXPECT_EQ(murk2::aa::canonicalise_mod(minus_one), twenty);
}

TEST(murk2aa, ModWrapping) {
  static murk2::aa::mod_ring ring(21);

  auto zero = ring(0);
  auto one = ring(1);
  auto twenty = ring(20);

  EXPECT_EQ(one + twenty, zero);
}

TEST(murk2aa, ModMulInv) {
  static murk2::aa::mod_ring ring(21);

  auto one = ring(1);
  auto two = ring(2);
  auto eleven = ring(11);

  ASSERT_NO_THROW(two ^ -1);
  EXPECT_EQ(two ^ -1, eleven);
  EXPECT_EQ(two.try_invert().value(), eleven) << two.try_invert().value().elem << " != " << eleven.elem;
  EXPECT_EQ(two * (two^(-1)), one);
}

TEST(murk2aa, ModMulInvFail) {
  static murk2::aa::mod_ring ring(21);

  auto three = ring(3);

  EXPECT_FALSE(three.try_invert());
  EXPECT_THROW(three ^ -1, murk2::aa::missing_structure);
}

//TEST(murk2aa, ModMulFieldBruteForce) {
//  static murk2::aa::mod_prime_field field(7);

//  auto absorber = field(field.add()->identity());
//  EXPECT_EQ(absorber.elem, 0);
//  EXPECT_THROW(absorber.invert(), murk2::aa::missing_structure);

//  murk2::bigint generator_val = murk2::aa::require_structure<murk2::aa::cyclic_group>(*field.add()).generator();
//  auto generator = murk2::aa::ring_element{c3lt::managed(&field), generator_val};
//  EXPECT_EQ(generator.elem, 1);
//  auto identity = murk2::aa::ring_element{c3lt::managed(&field), field.mul()->identity()};
//  EXPECT_EQ(generator, identity);

//  for (auto i = generator; i != absorber; i += generator)
//    EXPECT_EQ(i.invert() * i, identity);
//}

TEST(murk2aa, ModMulInvBig) {
  murk2::bigint modulus = 1;
  modulus <<= 110503;
  modulus -= 1;

  static murk2::aa::mod_prime_field ring(modulus);

  auto x = ring(42069);
  x.elem <<= 1024;
  auto x_inv = x.invert();

  auto res = x * x_inv;

  EXPECT_EQ(res, ring(ring.mul()->identity())) << x.elem << " * " << x_inv.elem << " == " << res.elem << " != " << ring.mul()->identity();
}
