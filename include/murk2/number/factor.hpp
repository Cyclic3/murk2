#pragma once

#include <murk2/common/bigint.hpp>
#include <murk2/number/prime.hpp>

#include <optional>

namespace murk2::number {
  /// @returns A non-trivial factor of num
  std::optional<bigint> pollard_rho(bigint const& num, uint64_t max_iters = -1);

  /// @returns A non-trivial factor of num
  /// @arg smoothness_prime a prime-valued smoothness bound
  std::optional<bigint> lenstra_ecm(bigint const& num, std::pair<uint64_t, uint64_t> bounds, uint64_t max_iters = -1);

  std::pair<uint64_t, uint64_t> lenstra_ecm_estimate_bounds(bigint const& num);
  inline std::optional<bigint> lenstra_ecm(bigint const& num, uint64_t max_iters = -1) {
    auto bounds = lenstra_ecm_estimate_bounds(num);
    return lenstra_ecm(num, bounds, max_iters);
  }
}
