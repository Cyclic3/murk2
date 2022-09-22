#pragma once

#include <murk2/crypto/key_agreement/common.hpp>

#include <murk2/aa/discrete_log.hpp>
#include <murk2/aa/group.hpp>

namespace murk2::crypto {
  /// Be aware that ECDH uses the abscissa (x-value) of the resultant coordinate
  template<typename CyclicGroup, typename = std::enable_if_t<std::derived_from<CyclicGroup, aa::cyclic_group<typename CyclicGroup::elem_t>>>>
  class diffie_hellman : public key_agreement<aa::group_element<CyclicGroup>, aa::group_element<CyclicGroup>, bigint> {
  public:
    c3lt::safe_ptr<const CyclicGroup> field;

  public:
    aa::group_element<CyclicGroup> make_pubkey(bigint const& priv) const override {
      return {field, field->op_iter(field->generator(), priv)};
    }
    aa::group_element<CyclicGroup> derive_shared(bigint const& priv, aa::group_element<CyclicGroup> const& pub) const override {
      return {field, field->op_iter(pub.elem, priv)};
    }

  public:
    diffie_hellman(c3lt::safe_ptr<const CyclicGroup> field_) : field{field_} {}
  };

  template<typename CyclicGroup>
  diffie_hellman(c3lt::safe_ptr<CyclicGroup>) -> diffie_hellman<CyclicGroup>;

  template<typename CyclicGroup>
  std::optional<bigint> diffie_hellman_crack(aa::group_element<CyclicGroup> pubkey) {
    return discrete_log(pubkey);
  }
}
