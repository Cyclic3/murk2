#pragma once

#include <murk2/number/factor.hpp>

#include <murk2/common/bigint.hpp>

namespace murk2::number {
  bigint euler_totient_inner(bigint const& prime, uint64_t power);
  inline bigint euler_totient() { return 1; }
  bigint euler_totient(const factorisation_t& primes_and_exps);

  // For a single prime power
  bigint carmichael_inner(bigint const& prime, uint64_t power);

  inline bigint carmichael() { return 1; }
  template<typename... Remaining>
  bigint carmichael(bigint const& prime, uint64_t power, Remaining&&... remaining) {
    return lcm(carmichael_inner(prime, power), carmichael(std::forward<Remaining>(remaining)...));
  }
  // Special case
  bigint carmichael_semiprime(bigint const& prime1, bigint const& prime2);
  // For all remaining cases
  bigint carmichael(factorisation_t const& primes_and_exps);

  /// Chinese remainder theorem solver
  std::pair<bigint, bigint> crt(std::vector<std::pair<bigint, bigint>> residues_and_moduli);
  /// Chinese remainder theorem solver for coprime case
  std::pair<bigint, bigint> crt_one(std::pair<bigint, bigint> const& a, std::pair<bigint, bigint> const& b);
  /// Chinese remainder theorem solver for coprime case
  std::pair<bigint, bigint> crt_coprime(std::vector<std::pair<bigint, bigint>> residues_and_moduli);
  /// Chinese remainder theorem solver for coprime case
  std::pair<bigint, bigint> crt_coprime_one(std::pair<bigint, bigint> const& a, std::pair<bigint, bigint> const& b);
}
