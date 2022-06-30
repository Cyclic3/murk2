#pragma once

#include <murk2/crypto/asymmetric.hpp>

#include <murk2/aa/modular.hpp>

#include <murk2/number/factor.hpp>

namespace murk2::crypto {
  struct rsa_pubkey {
    bigint exponent;
    bigint modulus;
  };

  struct rsa_privkey {
    bigint exponent;
    bigint modulus;
  };

  std::pair<rsa_pubkey, rsa_privkey> derive_rsa_key(bigint const& p, bigint const& q, bigint const& exponent = 3, bool set_pub_exponent = true);
  /// For multiprime keys
  std::pair<rsa_pubkey, rsa_privkey> derive_rsa_key(number::factorisation_t const&, bigint const& exponent = 3, bool set_pub_exponent = true);

  class rsa : public signature_alg<rsa_pubkey, rsa_privkey, bigint>, public asymmetric_encryption_alg<rsa_pubkey, rsa_privkey, bigint> {
  public:
    bigint sign(rsa_privkey const& key, bigint const& input) const;
    bool verify(rsa_pubkey const& key, bigint const& input, bigint const& signature) const;

    bigint encrypt(rsa_pubkey const& key, bigint const& input) const;
    bigint decrypt(rsa_privkey const& key, bigint const& output) const;
  };

  std::optional<rsa_privkey> crack_key(rsa_pubkey const&);
}
