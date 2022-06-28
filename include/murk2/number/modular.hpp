#pragma once

#include <murk2/common/bigint.hpp>

namespace murk2::number {
  bigint euler_totient_inner(bigint const& prime, uint64_t power);
  inline bigint euler_totient() { return 1; }
  template<typename... Remaining>
  bigint euler_totient(bigint const& prime, uint64_t power, Remaining&&... remaining) {
    return lcm(euler_totient_inner(prime, power), euler_totient(std::forward<Remaining>(remaining)...));
  }
  bigint euler_totient(std::vector<std::pair<bigint, uint64_t>> primes_and_exps);

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
  bigint carmichael(std::vector<std::pair<bigint, uint64_t>> primes_and_exps);
}
