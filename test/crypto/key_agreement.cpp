#include <gtest/gtest.h>

#include <murk2/crypto/key_agreement/diffie_hellman.hpp>

#include <murk2/aa/modular.hpp>
#include <murk2/aa/elliptic_curve.hpp>

TEST(murk2crypto, Dh) {
  murk2::aa::mod_mul_prime_group group{13, 2};
  murk2::crypto::diffie_hellman kex{c3lt::managed(&group)};

  murk2::bigint priv1 = 3;
  murk2::bigint priv2 = 9;

  auto pub1 = kex.make_pubkey(priv1);
  auto pub2 = kex.make_pubkey(priv2);

  auto shared1 = kex.derive_shared(priv1, pub2);
  auto shared2 = kex.derive_shared(priv2, pub1);

  EXPECT_EQ(shared1, shared2);
  EXPECT_EQ(shared1.elem, 8);
}

TEST(murk2crypto, DhMODP2048With256BitSubgroup) {
  murk2::aa::mod_mul_prime_group parent_group{murk2::bigint{"0x87A8E61DB4B6663CFFBBD19C651959998CEEF608660DD0F25D2CEED4435E3B00E00DF8F1D61957D4FAF7DF4561B2AA3016C3D91134096FAA3BF4296D830E9A7C209E0C6497517ABD5A8A9D306BCF67ED91F9E6725B4758C022E0B1EF4275BF7B6C5BFC11D45F9088B941F54EB1E59BB8BC39A0BF12307F5C4FDB70C581B23F76B63ACAE1CAA6B7902D52526735488A0EF13C6D9A51BFA4AB3AD8347796524D8EF6A167B5A41825D967E144E5140564251CCACB83E6B486F6B3CA3F7971506026C0B857F689962856DED4010ABD0BE621C3A3960A54E710C375F26375D7014103A4B54330C198AF126116D2276E11715F693877FAD7EF09CADB094AE91E1A1597"}};
  murk2::aa::cyclic_subgroup group{c3lt::managed(&parent_group), murk2::bigint{"0x3FB32C9B73134D0B2E77506660EDBD484CA7B18F21EF205407F4793A1A0BA12510DBC15077BE463FFF4FED4AAC0BB555BE3A6C1B0C6B47B1BC3773BF7E8C6F62901228F8C28CBB18A55AE31341000A650196F931C77A57F2DDF463E5E9EC144B777DE62AAAB8A8628AC376D282D6ED3864E67982428EBC831D14348F6F2F9193B5045AF2767164E1DFC967C1FB3F2E55A4BD1BFFE83B9C80D052B985D182EA0ADB2A3B7313D3FE14C8484B1E052588B9B7D2BBD2DF016199ECD06E1557CD0915B3353BBB64E0EC377FD028370DF92B52C7891428CDC67EB6184B523D1DB246C32F63078490F00EF8D647D148D47954515E2327CFEF98C582664B4C0F6CC41659"}};
  murk2::crypto::diffie_hellman kex{c3lt::managed(&group)};

  murk2::bigint order{"0x8CF83642A709A097B447997640129DA299B1A47D1EB3750BA308B0FE64F5FBD3"};

  murk2::bigint priv1{"15465687299313203123107744777962903225766691206903108469843505699655530827446"};
  murk2::bigint priv2{"491469383105845453898377478124598850203215052608712351998708989026032062936"};

  auto pub1 = kex.make_pubkey(priv1);
  auto pub2 = kex.make_pubkey(priv2);

  auto shared1 = kex.derive_shared(priv1, pub2);
  auto shared2 = kex.derive_shared(priv2, pub1);
  murk2::bigint expected_shared{"7794030862298244700678946066007909732991935227992811667202939070523499707297860119448877384062369422942882086543740016432195518533249456783210885893455651417150073632928998996579981747903602434660408241795052106559124827721020203795901941376206173663651264853948350569884908015640610242223791067463693954346999616891742411972504967194618691490332712082669623188544837819574891321575795986654588599776627876409610924570031286641220029348431935042678677342282289229362139472152362163718474039807532822230065347656364434473300121401876354029201165351271474186903349154977262821899041146106361172623968389877922949225104"};

  EXPECT_EQ(shared1, shared2);
  EXPECT_EQ(shared1.elem, expected_shared);
}

TEST(murk2crypto, Ecdh) {
  murk2::aa::mod_prime_field ground_field(13);
  // Checked by hand to by non-singular and prime order
  murk2::aa::elliptic_curve_group parent_group(c3lt::managed(&ground_field), 0, -2);
  murk2::aa::cyclic_subgroup group{parent_group(1, 5)};
  murk2::crypto::diffie_hellman kex{c3lt::managed(&group)};

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
  static murk2::aa::elliptic_curve_group parent_group(c3lt::managed(&ground_field), 0, 7);
  murk2::aa::cyclic_subgroup group{parent_group(generator[0], generator[1])};
  static murk2::crypto::diffie_hellman kex{c3lt::managed(&group)};

  murk2::bigint priv1{"0x1109078cda9dfc2fac17b817ec4f601d2437ebbab3d15d544551221fb8ceb313"};
  murk2::bigint priv2{"0x25eccca92a77a89b3b303ce3ef604a286853c39664a7de8ab45f1baab07b1048"};

  auto expected_pub1 = group(parent_group(murk2::bigint{"0xf6dd6ef5c2f63c9223fc8d2f5f420a4137fc403dc92d6b7e09fd24b748c48803"},
                                          murk2::bigint{"0xbec6bf235e11510a44c340477cce7104a0c0a0c01c4e91d3f60afea2efdc8539"}));
  auto expected_pub2 = group(parent_group(murk2::bigint{"0x1a537513bae37825b8dfa62375dfbbc54a8215f72ceda56c626ba1e70e94de95"},
                                          murk2::bigint{"0x65203bf187dc0f1bb025fa47233e60adcb5b8d3c93b8df682ceab2614114a43f"}));

  auto pub1 = kex.make_pubkey(priv1);
  auto pub2 = kex.make_pubkey(priv2);

  EXPECT_EQ(pub1, expected_pub1);
  EXPECT_EQ(pub2, expected_pub2);

  auto shared1 = kex.derive_shared(priv1, pub2);
  auto shared2 = kex.derive_shared(priv2, pub1);

  EXPECT_EQ(shared1, shared2);

  murk2::bigint expected_shared_x{"0x9fe36641fa24aa5f0f248e9e7ac97a6801b0aac1386acd39def0faf6a531824d"};

  EXPECT_EQ(shared1.elem.finite_point()[0], expected_shared_x);
}

TEST(murk2crypto, DhCrack) {
  murk2::aa::mod_mul_prime_group group{13, 2};
  murk2::crypto::diffie_hellman kex{c3lt::managed(&group)};

  murk2::bigint priv1 = 3;
  murk2::bigint priv2 = 9;

  auto pub1 = kex.make_pubkey(priv1);
  auto pub2 = kex.make_pubkey(priv2);


  auto shared1 = kex.derive_shared(priv1, pub2);
  auto priv1_cracked = murk2::crypto::diffie_hellman_crack(pub1);

  EXPECT_TRUE(priv1_cracked);
  EXPECT_EQ(priv1, *priv1_cracked);
}
