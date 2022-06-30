#pragma once

#include <murk2/common/bigint.hpp>
#include <murk2/number/prime.hpp>

#include <map>
#include <optional>

namespace murk2::number {
  using factorisation_t = std::map<bigint, uint64_t>;

  /// @returns A non-trivial factor of num
  std::optional<bigint> pollard_rho(bigint const& num, uint64_t max_iters = -1);

  /// @returns A non-trivial factor of num
  /// @arg smoothness_prime a prime-valued smoothness bound
  std::optional<bigint> lenstra_ecm(bigint const& num, std::pair<uint64_t, uint64_t> bounds, uint64_t max_iters = -1);

  std::optional<std::pair<bigint, uint64_t>> check_if_power(bigint const& num);

  std::pair<uint64_t, uint64_t> lenstra_ecm_estimate_bounds(bigint const& num);
  inline std::optional<bigint> lenstra_ecm(bigint const& num, uint64_t max_iters = -1) {
    auto bounds = lenstra_ecm_estimate_bounds(num);
    return lenstra_ecm(num, bounds, max_iters);
  }

  factorisation_t trial_division(bigint num);
  factorisation_t trial_division(bigint& num, bigint const& bound);

  struct factor_results {
    factorisation_t factorisation;
    factorisation_t failed;
  };

  factor_results factor(bigint num);
  bigint unfactor(factorisation_t const& factors);
}

inline std::ostream& operator<<(std::ostream& os, murk2::number::factorisation_t const& fact) {
  if (fact.empty())
    return os << "1^1";

  auto fake_end = fact.end();
  --fake_end;

  for (auto iter = fact.begin(); iter != fake_end; ++iter) {
    os << iter->first << "^" << iter->second << " * ";
  }
  os << fake_end->first << "^" << fake_end->second;
  return os;
}
