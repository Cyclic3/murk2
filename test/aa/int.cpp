#include <gtest/gtest.h>

#include <murk2/aa/int.hpp>

TEST(murk2aa, IntRing1Plus1) {
  murk2::aa::int_ring ring;

  auto one = ring(1);
  auto two = ring(2);

  EXPECT_EQ(one + one, two);
}

TEST(murk2aa, IntRing2Times3) {
  murk2::aa::int_ring ring;

  auto two = ring(2);
  auto three = ring(3);
  auto six = ring(6);

  EXPECT_EQ(two * three, six);
}
