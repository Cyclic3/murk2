#pragma once

#include <murk2/number/factor.hpp>
#include <murk2/number/modular.hpp>

#include <murk2/aa/group.hpp>
#include <murk2/aa/modular.hpp>

#include <murk2/common/simulation.hpp>

namespace murk2::aa {
  // TODO: check out anomalous and supersingular dlp

  // Pollard's rho algorithm for discrete logarithms with Teske's Adding-walk
  //
  //
  // SOURCE: https://maths-people.anu.edu.au/~brent/pd/rpb231.pdf
  // SOURCE: https://www.ams.org/journals/mcom/2001-70-234/S0025-5718-00-01213-8/S0025-5718-00-01213-8.pdf
  // Algorithm generic outline from SOURCE: https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm_for_logarithms
  template<typename Group, typename = std::enable_if_t<std::derived_from<Group, cyclic_group<typename Group::elem_t>>>>
  bigint pollard_rho_dl(group_element<Group> target, std::vector<bigint> m, std::vector<bigint> n) {
    Group const& group = *target.context;
    group_element<Group> generator{target.context, group.generator()};

    if (m.size() != n.size())
      throw std::invalid_argument{"There must be the same number of m values as n values"};

    size_t r = m.size();

    std::vector<group_element<Group>> M;
    M.reserve(r);
    for (size_t i = 0; i < r; ++i)
      M.emplace_back((generator * m[i]) + (target * n[i]));

    group_element<Group> x_0{target.context, group.identity()};
    bigint a_0 = 0;
    bigint b_0 = 0;

    group_element<Group> x_i = x_0;
    bigint a_i;
    bigint b_i;

    group_element<Group> x_j = x_0;
    bigint a_j;
    bigint b_j;

    auto do_iter = [&M, &m, &n, &r](auto& x, bigint& a, bigint& b) {
      auto x_class = simulation_hash(x) % r;
      x += M.at(x_class);
      a += m[x_class];
      b += n[x_class];
    };

    mod_ring ring(group.order());

    std::optional<ring_element<aa::ring<bigint>>> b_mod;

    for (bigint n_resets = 0; true; ++n_resets) {
      x_i = x_j = x_0;
      a_i = a_j = a_0;
      b_i = b_j = b_0;

      do {
        do_iter(x_i, a_i, b_i);
        do_iter(x_j, a_j, b_j);
        do_iter(x_j, a_j, b_j);
      } while (x_i != x_j);

      if ((b_mod = ring(b_i - b_j).try_invert()))
        break;
      x_0 += generator;
      a_0 += 1;
    }

    auto a_mod = ring(a_j - a_i);

    return (a_mod * *b_mod).elem;
  }

  template<typename Group, typename = std::enable_if_t<std::derived_from<Group, cyclic_group<typename Group::elem_t>>>>
  bigint pollard_rho_dl(group_element<Group> target, size_t n_gen = 20) {
    simulation_rng rng;
    bigint order = target.context->order();

    std::vector<bigint> m(n_gen);
    std::vector<bigint> n(n_gen);
    std::generate(m.begin(), m.end(), std::bind(std::ref(rng), 1, order));
    std::generate(n.begin(), n.end(), std::bind(std::ref(rng), 1, order));

    return pollard_rho_dl<Group>(std::move(target), std::move(m), std::move(n));
  }

  template<typename Group, typename = std::enable_if_t<std::derived_from<Group, cyclic_group<typename Group::elem_t>>>>
  std::optional<bigint> discrete_log_prime_power(group_element<Group> const& target) {
    auto order = target.context->order();
    if (order < 20) {
      group_element<Group> g{target.context, target.context->generator()};
      group_element<Group> trial{target.context, target.context->identity()};
      for (bigint i = 0; i <= order; ++i, trial += g) {
        if (trial == target)
          return i;
      }
      throw std::logic_error{"discrete_log_prime got incorrect group order"};
    }

    return pollard_rho_dl(target);
  }

  template<typename Group, typename = std::enable_if_t<std::derived_from<Group, cyclic_group<typename Group::elem_t>>>>
  std::optional<bigint> discrete_log(group_element<Group> const& target) {
    // Pohlig–Hellman algorithm implementation
    auto n = target.context->order();
    group_element<Group> g{target.context, target.context->generator()};
    auto factors = number::factor(n);

    if (!factors.failed.empty())
      return std::nullopt;

    std::vector<std::pair<bigint, bigint>> congs;

    for (auto& [p, e] : factors.factorisation) {
      bigint subgroup_order = do_gmp(mpz_pow_ui, p, e);

      bigint remaining = n/subgroup_order;
      group_element<Group> g_i = g * remaining;
      group_element<Group> h_i = target * remaining;
      cyclic_subgroup<Group, true> subgroup{g_i, subgroup_order};

      auto res = discrete_log_prime_power(subgroup(h_i));
      if (!res)
        return std::nullopt;

      congs.emplace_back(*res, subgroup_order);
    }

    return number::crt_coprime(congs).first;
  }
}
