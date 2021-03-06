#include <murk2/crypto/rsa.hpp>

#include <murk2/number/factor.hpp>
#include <murk2/number/modular.hpp>

namespace murk2::crypto {
  bigint rsa::sign(rsa_privkey const& key, bigint const& input) const {
    aa::mod_ring ring(key.modulus);
    return (ring(input) ^ key.exponent).elem;
  }
  bool rsa::verify(rsa_pubkey const& key, bigint const& input, bigint const& signature) const {
    aa::mod_ring ring(key.modulus);
    return (ring(input) ^ key.exponent).elem == signature;
  }

  bigint rsa::encrypt(rsa_pubkey const& key, bigint const& input) const {
    aa::mod_ring ring(key.modulus);
    return (ring(input) ^ key.exponent).elem;
  }
  bigint rsa::decrypt(rsa_privkey const& key, bigint const& output) const {
    aa::mod_ring ring(key.modulus);
    return (ring(output) ^ key.exponent).elem;
  }

  std::pair<rsa_pubkey, rsa_privkey> derive_rsa_key(bigint const& p, bigint const& q, bigint const& exponent, bool set_pub_exponent) {
    std::pair<rsa_pubkey, rsa_privkey> ret;
    ret.first.modulus = ret.second.modulus = p * q;

    bigint* known_exponent;
    bigint* unknown_exponent;
    if (set_pub_exponent) {
      known_exponent = &ret.first.exponent;
      unknown_exponent = &ret.second.exponent;
    }
    else {
      known_exponent = &ret.second.exponent;
      unknown_exponent = &ret.first.exponent;
    }

    *known_exponent = exponent;

    aa::mod_mul_monoid monoid(number::carmichael_semiprime(p, q));
    if (auto res = monoid.try_invert(*known_exponent))
      *unknown_exponent = std::move(*res);
    else
      throw std::invalid_argument{"Non-invertible exponent"};

    return ret;
  }


  std::pair<rsa_pubkey, rsa_privkey> derive_rsa_key(number::factorisation_t const& factors, bigint const& exponent, bool set_pub_exponent) {
    std::pair<rsa_pubkey, rsa_privkey> ret;
    ret.first.modulus = ret.second.modulus = number::unfactor(factors);

    bigint* known_exponent;
    bigint* unknown_exponent;
    if (set_pub_exponent) {
      known_exponent = &ret.first.exponent;
      unknown_exponent = &ret.second.exponent;
    }
    else {
      known_exponent = &ret.second.exponent;
      unknown_exponent = &ret.first.exponent;
    }

    *known_exponent = exponent;

    aa::mod_mul_monoid monoid(number::carmichael(factors));
    if (auto res = monoid.try_invert(*known_exponent))
      *unknown_exponent = std::move(*res);
    else
      throw std::invalid_argument{"Non-invertible exponent"};

    return ret;
  }

  std::optional<rsa_privkey> crack_key(rsa_pubkey const& key) {
    auto factor = number::factor(key.modulus);
    if (!factor.failed.empty())
      return std::nullopt;

    return derive_rsa_key(factor.factorisation, key.exponent, true).second;
  }
}
