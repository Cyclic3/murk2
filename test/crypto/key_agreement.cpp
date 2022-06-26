#include <gtest/gtest.h>

#include <murk2/crypto/key_agreement/diffie_hellman.hpp>

#include <murk2/aa/modular.hpp>
#include <murk2/aa/elliptic_curve.hpp>

TEST(murk2crypto, Dh) {
  static murk2::aa::mod_mul_prime_group group(13, 2);
  static murk2::crypto::diffie_hellman kex{c3lt::managed(&group)};

  murk2::bigint priv1 = 3;
  murk2::bigint priv2 = 9;

  auto pub1 = kex.make_pubkey(priv1);
  auto pub2 = kex.make_pubkey(priv2);

  auto shared1 = kex.derive_shared(priv1, pub2);
  auto shared2 = kex.derive_shared(priv2, pub1);

  EXPECT_EQ(shared1, shared2);
  EXPECT_EQ(shared1.elem, 8);
}

TEST(murk2crypto, Ecdh) {
  static murk2::aa::mod_prime_field ground_field(13);
  // Checked by hand to by non-singular and prime order
  static murk2::aa::elliptic_curve_group<murk2::aa::mod_prime_field, murk2::aa::elliptic_group_type::FiniteCyclic> group(c3lt::managed(&ground_field), 0, -2, std::array<murk2::bigint, 2>{1, 5});
  static murk2::crypto::diffie_hellman kex{c3lt::managed(&group)};

  murk2::bigint priv1 = 3;
  murk2::bigint priv2 = 9;

  auto pub1 = kex.make_pubkey(priv1);
  auto pub2 = kex.make_pubkey(priv2);

  auto shared1 = kex.derive_shared(priv1, pub2);
  auto shared2 = kex.derive_shared(priv2, pub1);

  EXPECT_EQ(shared1, shared2);
}

TEST(murk2crypto, EcdhSecp256k1) {
  static murk2::aa::mod_prime_field ground_field(murk2::bigint("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F"));
  std::array<murk2::bigint, 2> generator{murk2::bigint{"0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798"},
                                         murk2::bigint{"0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8"}};
  static murk2::aa::elliptic_curve_group<murk2::aa::mod_prime_field, murk2::aa::elliptic_group_type::FiniteCyclic> group(c3lt::managed(&ground_field), 0, 7, generator);
  static murk2::crypto::diffie_hellman kex{c3lt::managed(&group)};

  murk2::bigint priv1{"0x1109078cda9dfc2fac17b817ec4f601d2437ebbab3d15d544551221fb8ceb313"};
  murk2::bigint priv2{"0x25eccca92a77a89b3b303ce3ef604a286853c39664a7de8ab45f1baab07b1048"};

  auto expected_pub1 = group(murk2::bigint{"0xf6dd6ef5c2f63c9223fc8d2f5f420a4137fc403dc92d6b7e09fd24b748c48803"},
                             murk2::bigint{"0xbec6bf235e11510a44c340477cce7104a0c0a0c01c4e91d3f60afea2efdc8539"});
  auto expected_pub2 = group(murk2::bigint{"0x1a537513bae37825b8dfa62375dfbbc54a8215f72ceda56c626ba1e70e94de95"},
                             murk2::bigint{"0x65203bf187dc0f1bb025fa47233e60adcb5b8d3c93b8df682ceab2614114a43f"});

  auto pub1 = kex.make_pubkey(priv1);
  auto pub2 = kex.make_pubkey(priv2);

  EXPECT_EQ(pub1, expected_pub1);
  EXPECT_EQ(pub2, expected_pub2);

  auto shared1 = kex.derive_shared(priv1, pub2);
  auto shared2 = kex.derive_shared(priv2, pub1);
  murk2::bigint expected_shared_x{"0x9fe36641fa24aa5f0f248e9e7ac97a6801b0aac1386acd39def0faf6a531824d"};

  EXPECT_EQ(shared1.elem.finite_point()[0], expected_shared_x);
}
