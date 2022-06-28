#include <gtest/gtest.h>

#include <murk2/aa/discrete_log.hpp>

template<typename Group>
uint64_t f(murk2::aa::group_element<Group> const& i) {
  return murk2::simulation_hash(i.elem);
}


TEST(murk2aa, ModAddPollardRhoDL) {
  murk2::aa::mod_add_group group(101);
  murk2::aa::group_element gen{c3lt::managed(&group), group.generator()};
  murk2::bigint exponent = group.order() / 2;

  auto target = gen * exponent;

  auto res = murk2::aa::pollard_rho_dl<decltype(group)>(target, f<decltype(group)>);

  EXPECT_EQ(res, exponent);
  EXPECT_EQ(gen*res, target);
}

TEST(murk2aa, ModMulPollardRhoDL) {
  murk2::aa::mod_mul_prime_group parent_group(4941056177);
  murk2::bigint generator{"2853960960"};
  murk2::bigint order{"139169"};
  murk2::aa::cyclic_subgroup<murk2::aa::mod_mul_prime_group, true> group(c3lt::managed(&parent_group), generator, order);
  murk2::aa::group_element gen{c3lt::managed(&group), group.generator()};
  murk2::bigint exponent = group.order() / 2;

  auto target = gen * exponent;

  auto res = murk2::aa::pollard_rho_dl<decltype(group)>(target, f<decltype(group)>);

  EXPECT_EQ(res, exponent);
  EXPECT_EQ(gen*res, target);
}

// TODO: elliptic curve test
