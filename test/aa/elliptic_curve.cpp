#include <gtest/gtest.h>

#include <murk2/aa/modular.hpp>
#include <murk2/aa/elliptic_curve.hpp>

static murk2::aa::mod_prime_field const& get_field() {
  static thread_local murk2::aa::mod_prime_field ret{11};
  return ret;
}
static auto& get_curve() {
  // Hand checked non-singular
  static thread_local murk2::aa::elliptic_curve_group ret{c3lt::managed(&get_field()), -7, 6};
  return ret;
}

TEST(murk2aa, EllipticPrimeFiniteFieldIdentity) {
  auto& group = get_curve();

  auto id = group(group.identity());
  auto root1 = group(1, 0);

  EXPECT_EQ(id + id, id);
  EXPECT_EQ(root1 + id, root1);
  EXPECT_EQ(id + root1, root1);
}

TEST(murk2aa, EllipticFiniteField) {
  auto& group = get_curve();

  auto root1 = group(1, 0);
  auto root2 = group(2, 0);
  auto root3 = group(8, 0);
  auto point1 = group(4, 3);
  auto point1_inv = group(4, 8);
  auto point2 = group(7, 5);
  auto point3 = group(9, 1);
  auto id = group(group.identity());

  EXPECT_EQ(root1 + root1, id) << "2 * " << root1.elem << " != " << id.elem;
  EXPECT_EQ(root1 * 2, id);
  EXPECT_EQ(root1 * 3, root1);
  EXPECT_EQ(root1 + root2, root3) << root1.elem << " + " << root2.elem << " != " << root3.elem;
  EXPECT_TRUE(group.is_inverse(point1.elem, point1_inv.elem));
  EXPECT_EQ(point1 + id, point1);
  EXPECT_EQ(point1 + root1, point2);
  EXPECT_EQ(point1 + point2, point3);
}

static_assert(std::remove_cvref_t<decltype(get_curve())>::type != murk2::aa::elliptic_group_type::Generic, "Elliptic curve type was incorrectly derived as Generic");
static_assert(std::remove_cvref_t<decltype(get_curve())>::type == murk2::aa::elliptic_group_type::Finite, "Elliptic curve type was incorrectly derived");
