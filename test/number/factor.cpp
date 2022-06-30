#include <gtest/gtest.h>

#include <murk2/number/factor.hpp>


TEST(murk2number, PollardRhoEz) {
  auto factor_opt = murk2::number::pollard_rho(15);
  ASSERT_TRUE(factor_opt);

  auto& factor = *factor_opt;
  EXPECT_TRUE(factor == 3 || factor == 5) << factor;
}

TEST(murk2number, PollardRho30) {
  murk2::bigint p{"49417"};
  murk2::bigint q{"50789"};

  auto factor_opt = murk2::number::pollard_rho(p * q);
  ASSERT_TRUE(factor_opt);

  auto& factor = *factor_opt;
  EXPECT_TRUE(factor == p || factor == q) << factor;
}

TEST(murk2number, PollardRho70) {
  murk2::bigint p{"8447004649"};
  murk2::bigint q{"2458956821"};

  auto factor_opt = murk2::number::pollard_rho(p * q);
  ASSERT_TRUE(factor_opt);

  auto& factor = *factor_opt;
  EXPECT_TRUE(factor == p || factor == q) << factor;
}

TEST(murk2number, PollardRhoStopCheck) {
  murk2::bigint p{"8447004649"};

  auto factor_opt = murk2::number::pollard_rho(p, 128);

  EXPECT_FALSE(factor_opt);
}

TEST(murk2number, LenstraECMEz) {
  auto factor_opt = murk2::number::lenstra_ecm(15);
  ASSERT_TRUE(factor_opt);

  auto& factor = *factor_opt;
  EXPECT_TRUE(factor == 3 || factor == 5) << factor;
}


TEST(murk2number, LenstraECM13) {
  murk2::bigint p{"17"};
  murk2::bigint q{"59"};

  auto factor_opt = murk2::number::lenstra_ecm(p * q);
  ASSERT_TRUE(factor_opt);

  auto& factor = *factor_opt;
  EXPECT_TRUE(factor == p || factor == q) << factor;
}

TEST(murk2number, LenstraECM30) {
  murk2::bigint p{"49417"};
  murk2::bigint q{"50789"};

  auto factor_opt = murk2::number::lenstra_ecm(p * q);
  ASSERT_TRUE(factor_opt);

  auto& factor = *factor_opt;
  EXPECT_TRUE(factor == p || factor == q) << factor;
}

TEST(murk2number, LenstraECM70) {
  murk2::bigint p{"8447004649"};
  murk2::bigint q{"2458956821"};

  auto factor_opt = murk2::number::lenstra_ecm(p * q);
  ASSERT_TRUE(factor_opt);

  auto& factor = *factor_opt;
  EXPECT_TRUE(factor == p || factor == q) << factor;
}

TEST(murk2number, LenstraECM100) {
  murk2::bigint p{"563594796096133"};
  murk2::bigint q{"122904209049907"};

  auto factor_opt = murk2::number::lenstra_ecm(p * q);
  ASSERT_TRUE(factor_opt);

  auto& factor = *factor_opt;
  EXPECT_TRUE(factor == p || factor == q) << factor;
}

//TEST(murk2number, LenstraECM128) {
//  murk2::bigint p{"39986299716732592789"};
//  murk2::bigint q{"33333395883311093771"};

//  auto factor_opt = murk2::number::lenstra_ecm(p * q);
//  ASSERT_TRUE(factor_opt);

//  auto& factor = *factor_opt;
//  EXPECT_TRUE(factor == p || factor == q) << factor;
//}


TEST(murk2number, LenstraECMStop) {
  murk2::bigint p{"8447004649"};

  auto factor_opt = murk2::number::lenstra_ecm(p, 8);
  EXPECT_FALSE(factor_opt);
}

TEST(murk2number, Factor) {
  murk2::bigint p{"17"};
  murk2::bigint q{"50789"};
  murk2::bigint r{"563594796096133"};

  murk2::bigint n = (p * p * p * p * p) * (q * r) * 2;
  murk2::number::factorisation_t expected {
    {2, 1},
    {p, 5},
    {q, 1},
    {r, 1}
  };

  auto factor_res = murk2::number::factor(n);
  EXPECT_TRUE(factor_res.failed.empty());

  EXPECT_EQ(factor_res.factorisation, expected);
}
