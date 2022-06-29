#pragma once

#include <murk2/common/bigint.hpp>
#include <murk2/common/err.hpp>

namespace murk2::number {
  inline bool is_prime(bigint const& n) { return mpz_probab_prime_p(n, 32); }
  inline bigint next_prime(bigint const& after) { return do_gmp(mpz_nextprime, after); }
}
