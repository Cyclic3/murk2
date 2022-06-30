#pragma once

#include <murk2/common/bigint.hpp>

namespace murk2 {
  class simulation_rng {
  private:
    gmp_randstate_t state;

  public:
    inline bigint operator()(bigint const& upper_bound) {
      return do_gmp(mpz_urandomm, state, upper_bound);
    }
    inline bigint operator()(bigint const& lower_bound, bigint const& upper_bound) {
      return do_gmp(mpz_urandomm, state, bigint{upper_bound - lower_bound}) + lower_bound;
    }

  public:
    inline simulation_rng() { gmp_randinit_default(state); }
    inline simulation_rng(simulation_rng const&) = delete;
    inline simulation_rng(simulation_rng&&) = delete;

    inline ~simulation_rng() { gmp_randclear(state); }
  };

  uint64_t simulation_hash(bigint const&);
  template<std::integral T>
  uint64_t simulation_hash(T const& t) { return simulation_hash(bigint{t}); }
  template<typename T>
  uint64_t simulation_hash(T const* t) { return simulation_hash(reinterpret_cast<uintptr_t>(t)); }
}
