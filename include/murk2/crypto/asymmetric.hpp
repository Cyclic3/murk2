#pragma once

namespace murk2::crypto {
  template<typename PublicKey, typename PrivateKey, typename Input, typename Signature = Input>
  struct signature_alg {
    using public_key_t = PublicKey;
    using private_key_t = PrivateKey;
    using input_t = Input;
    using signature_t = Signature;

    virtual Signature sign(PrivateKey const& key, Input const& input) const = 0;
    virtual bool verify(PublicKey const& key, Input const& input, Signature const& signature) const = 0;

    virtual ~signature_alg() = default;
  };

  template<typename PublicKey, typename PrivateKey, typename PlainText, typename CypherText = PlainText>
  struct asymmetric_encryption_alg {
    using public_key_t = PublicKey;
    using private_key_t = PrivateKey;
    using plaintext_t = PlainText;
    using cyphertext_t = CypherText;

    virtual CypherText encrypt(PublicKey const& key, PlainText const& ptext) const = 0;
    virtual PlainText decrypt(PrivateKey const& key, CypherText const& ctext) const = 0;

    virtual ~asymmetric_encryption_alg() = default;
  };

  template<typename PublicKey, typename PrivateKey, typename Domain>
  struct asymmetric_encryption_alg<PublicKey, PrivateKey, Domain, Domain> {
    using public_key_t = PublicKey;
    using private_key_t = PrivateKey;
    using plaintext_t = Domain;
    using cyphertext_t = Domain;

    virtual Domain encrypt(PublicKey const& key, Domain const& ptext) const = 0;
    virtual void encrypt_mut(PublicKey const& key, Domain& x) const { x = encrypt(key, x); }
    virtual Domain decrypt(PrivateKey const& key, Domain const& ctext) const = 0;
    virtual void decrypt_mut(PrivateKey const& key, Domain& x) const { x = decrypt(key, x); }

    virtual ~asymmetric_encryption_alg() = default;
  };
}
