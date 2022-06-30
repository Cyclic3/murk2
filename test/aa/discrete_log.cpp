#include <gtest/gtest.h>

#include <murk2/aa/discrete_log.hpp>


TEST(murk2aa, ModAddPollardRhoDL) {
  murk2::aa::mod_add_group group(101);
  murk2::aa::group_element gen{c3lt::managed(&group), group.generator()};
  murk2::bigint exponent = (group.order() / 2) - 1;

  auto target = gen * exponent;

  auto res = murk2::aa::pollard_rho_dl<decltype(group)>(target);

  EXPECT_EQ(res, exponent);
  EXPECT_EQ(gen*res, target);
}

TEST(murk2aa, ModMulPollardRhoDL) {
  murk2::aa::mod_mul_prime_group parent_group(4941056177);
  murk2::bigint generator{"2853960960"};
  murk2::bigint order{"139169"};
  murk2::aa::cyclic_subgroup<murk2::aa::mod_mul_prime_group, true> group(c3lt::managed(&parent_group), generator, order);
  murk2::aa::group_element gen{c3lt::managed(&group), group.generator()};
  murk2::bigint exponent = (group.order() / 2) - 1;

  auto target = gen * exponent;

  auto res = murk2::aa::pollard_rho_dl<decltype(group)>(target);

  EXPECT_EQ(res, exponent);
  EXPECT_EQ(gen*res, target);
}

TEST(murk2aa, ModMulPollardRhoDL2) {
  murk2::aa::mod_mul_prime_group group(4941056177, 3);
  murk2::aa::group_element gen{c3lt::managed(&group), group.generator()};
  murk2::bigint exponent = (group.order() / 2) - 1;

  auto target = gen * exponent;

  auto res = murk2::aa::pollard_rho_dl<decltype(group)>(target);

  EXPECT_EQ(res, exponent);
  EXPECT_EQ(gen*res, target);
}

TEST(murk2aa, ModMulDLNonSmooth) {
  murk2::bigint p = 4941056177;
  // Find non-smooth value
  for (; !murk2::number::is_prime((p - 1)/2); p = murk2::number::next_prime(p));

  murk2::aa::mod_mul_prime_group group(p, 3);
  murk2::aa::group_element gen{c3lt::managed(&group), group.generator()};
  murk2::bigint exponent = (group.order() / 2) - 1;

  auto target = gen * exponent;

  auto res = murk2::aa::discrete_log(target);

  EXPECT_TRUE(res);

  EXPECT_EQ(*res, exponent);
  EXPECT_EQ(gen * *res, target);
}

TEST(murk2aa, ModMulDLSometimesEdgeCaseIdk) {
  murk2::aa::mod_mul_prime_group group(13, 2);
  murk2::aa::group_element gen{c3lt::managed(&group), group.generator()};
  murk2::bigint exponent = 9;

  auto target = gen * exponent;

  auto res = murk2::aa::discrete_log(target);

  EXPECT_TRUE(res);

  EXPECT_EQ(*res, exponent);
  EXPECT_EQ(gen * *res, target);
}


// TODO: elliptic curve test
