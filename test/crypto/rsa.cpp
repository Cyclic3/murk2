#include <gtest/gtest.h>

#include <murk2/crypto/rsa.hpp>

TEST(murk2crypto, RSA) {
  murk2::bigint p = 3;
  murk2::bigint q = 11;

  murk2::crypto::rsa rsa;

  auto [pubkey, privkey] = murk2::crypto::derive_rsa_key(p, q, 3, true);

  EXPECT_EQ(pubkey.modulus, privkey.modulus);
  EXPECT_EQ(pubkey.modulus, murk2::bigint{p * q});
  EXPECT_EQ(pubkey.exponent, 3);
  EXPECT_EQ(privkey.exponent, 7);
}


TEST(murk2crypto, RSA2048) {
  murk2::bigint p{"183041381519527082618197165590323097476437239654939254242079007499210106447748113496411210777382744061596816679719488257688320965228815068013703023249603718796353233515160683446887568146688074082930772971214458181926622109390522626258667801652308159687531728999128073927326044750301847660571434646563"};
  murk2::bigint q{"120495867045013207486889392880656905118052635563201700030073704728271406249043689811921305579486300409178611321991332074520961646519446581376165377678935931099543956263054192693153010331037424085025251248273345031981325257913208222002586041107884617532877325770469415638595964935684941974145533483149"};

  std::cout << q << std::endl;

  murk2::crypto::rsa rsa;

  auto [pubkey, privkey] = murk2::crypto::derive_rsa_key(p, q, 3, true);

  murk2::aa::mod_ring ring{pubkey.modulus};

  EXPECT_EQ(pubkey.modulus, privkey.modulus);
  EXPECT_EQ(pubkey.modulus, murk2::bigint{p * q});
  EXPECT_EQ(pubkey.exponent, 3);
  EXPECT_EQ((ring(2)^(privkey.exponent *pubkey.exponent)).elem, 2);
}
