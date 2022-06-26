#pragma once

namespace murk2::crypto {
  template<typename SharedKey, typename PubKey, typename PrivKey>
  struct key_agreement {
    using shared_key_t = SharedKey;
    using public_key_t = PubKey;
    using private_key_t = PrivKey;

    virtual PubKey make_pubkey(PrivKey const&) const = 0;
    virtual SharedKey derive_shared(PrivKey const& our, PubKey const& other) const = 0;
    virtual ~key_agreement() = default;
  };
}
