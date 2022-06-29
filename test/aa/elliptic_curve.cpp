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
static auto& get_proj_curve() {
  // Hand checked non-singular
  static thread_local murk2::aa::elliptic_curve_proj_group ret{c3lt::managed(&get_field()), -7, 6};
  return ret;
}

const murk2::bigint curve_order = 16;

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
  EXPECT_EQ(point1 * curve_order, id);
  EXPECT_EQ(root1 * 2, id);
  EXPECT_EQ(root1 * 3, root1);
  EXPECT_EQ(root1 + root2, root3) << root1.elem << " + " << root2.elem << " != " << root3.elem;
  EXPECT_TRUE(group.is_inverse(point1.elem, point1_inv.elem));
  EXPECT_EQ(point1 + id, point1);
  EXPECT_EQ(point1 + root1, point2);
  EXPECT_EQ(point1 + point2, point3);
}

TEST(murk2aa, EllipticProjFiniteFieldMatch) {
  auto& group = get_proj_curve();
  auto& affine_group = get_curve();

  auto root1 = group(1, 0, 1);
  auto point1 = group(4, 3, 1);

  auto root1_afn = affine_group(1, 0);
  auto point1_afn = affine_group(4, 3);

  EXPECT_TRUE(group.is_identity((root1 + root1).elem));
  EXPECT_TRUE(affine_group.is_identity((root1 + root1).elem.affinise(2)));
  EXPECT_EQ((root1 + point1).elem.affinise(2), (root1_afn + point1_afn).elem);
  EXPECT_EQ((point1 + point1 + point1).elem.affinise(2), (point1_afn + point1_afn + point1_afn).elem);
  EXPECT_EQ((root1 + point1 + point1).elem.affinise(2), (root1_afn + point1_afn + point1_afn).elem);
  EXPECT_TRUE(group.is_identity((point1 * curve_order).elem));
}

static_assert(std::remove_cvref_t<decltype(get_curve())>::type != murk2::aa::elliptic_group_type::Generic, "Elliptic curve type was incorrectly derived as Generic");
static_assert(std::remove_cvref_t<decltype(get_curve())>::type == murk2::aa::elliptic_group_type::Finite, "Elliptic curve type was incorrectly derived");
