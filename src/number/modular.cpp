#include <murk2/number/modular.hpp>

namespace murk2::number {
  bigint euler_totient_inner(bigint const& prime, uint64_t power) {
    return do_gmp(mpz_pow_ui, prime, power);
  }
  bigint euler_totient(std::vector<std::pair<bigint, uint64_t>> primes_and_exps) {
    bigint ret = 1;
    for (auto& [p, e] : primes_and_exps)
      ret *= euler_totient_inner(p, e);
    return ret;
  }

  // For a single prime power
  bigint carmichael_inner(bigint const& prime, uint64_t power) {
    if (prime == 2 || power >= 3)
      return euler_totient(prime, power) / 2;
    else
      return euler_totient(prime, power);
  }

  // Special case
  bigint carmichael_semiprime(bigint const& prime1, bigint const& prime2) {
    return do_gmp(mpz_lcm, bigint{prime1 - 1}, bigint{prime2 - 1});
  }
  // For all remaining cases
  bigint carmichael(std::vector<std::pair<bigint, uint64_t>> primes_and_exps) {
    bigint ret = 1;
    for (auto& [p, e] : primes_and_exps)
      mpz_lcm(ret, ret, carmichael_inner(p, e));
    return ret;
  }
}
