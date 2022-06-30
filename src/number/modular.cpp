#include <murk2/number/modular.hpp>

#include <murk2/aa/modular.hpp>

#include <numeric>

namespace murk2::number {
  bigint euler_totient_inner(bigint const& prime, uint64_t power) {
    return do_gmp(mpz_pow_ui, prime, power - 1) * (prime - 1);
  }
  bigint euler_totient(factorisation_t const& primes_and_exps) {
    bigint ret = 1;
    for (auto& [p, e] : primes_and_exps)
      ret *= euler_totient_inner(p, e);
    return ret;
  }

  // For a single prime power
  bigint carmichael_inner(bigint const& prime, uint64_t power) {
    if (prime == 2 || power >= 3)
      return euler_totient_inner(prime, power) / 2;
    else
      return euler_totient_inner(prime, power);
  }

  // Special case
  bigint carmichael_semiprime(bigint const& prime1, bigint const& prime2) {
    return do_gmp(mpz_lcm, bigint{prime1 - 1}, bigint{prime2 - 1});
  }
  // For all remaining cases
  bigint carmichael(factorisation_t const& primes_and_exps) {
    bigint ret = 1;
    for (auto& [p, e] : primes_and_exps)
      mpz_lcm(ret, ret, carmichael_inner(p, e));
    return ret;
  }

  std::pair<bigint, bigint> crt_coprime_one(std::pair<bigint, bigint> const& a, std::pair<bigint, bigint> const& b) {
    bigint m_1, m_2, one;
    mpz_gcdext(one, m_1, m_2, a.second, b.second);
    bigint modulus = a.second * b.second;
    return {murk2::aa::canonicalise_mod(bigint{(a.first * b.second * m_2) + (b.first * a.second * m_1)}, modulus), modulus};
  }
  std::pair<bigint, bigint> crt_coprime(std::vector<std::pair<bigint, bigint>> residues_and_moduli) {
    // `at` checks index for us
    return std::accumulate(residues_and_moduli.begin() + 1, residues_and_moduli.end(), residues_and_moduli.at(0), &crt_coprime_one);
  }
}
