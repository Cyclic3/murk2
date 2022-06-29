#include <murk2/number/factor.hpp>

#include <murk2/aa/elliptic_curve.hpp>
#include <murk2/aa/modular.hpp>

#include <murk2/common/simulation.hpp>

#include <cmath>
#include <chrono>
#include <map>

namespace murk2::number {
  static void pollard_rho_iter(auto& elem)  {
    elem = (elem * elem) + 1;
  }

  std::optional<bigint> pollard_rho(bigint const& num, uint64_t max_iters) {
    aa::mod_ring ring(num);

    auto tortoise = ring(2);
    auto hare = tortoise;
    pollard_rho_iter(hare);
    bigint ret;

    uint64_t iter = 0;
    for (; iter < max_iters && (ret = gcd(hare.elem - tortoise.elem, num)) == 1; pollard_rho_iter(tortoise), pollard_rho_iter(hare), pollard_rho_iter(hare));

    if (iter == max_iters || ret == num)
      return std::nullopt;

    return ret;
  }

  // https://hal.inria.fr/inria-00070192v1/document
  std::optional<bigint> lenstra_ecm(bigint const& num, std::pair<uint64_t, uint64_t> bounds, uint64_t max_iters) {
    simulation_rng rng;
    aa::mod_ring ring(num);

    std::vector<std::pair<bigint, uint64_t>> primes_and_exps;

    uint64_t iters = 0;

    bigint k = 1;
    double log_bound1 = std::log(bounds.first);
    for (bigint p = 2; p < bounds.first; p = next_prime(p)) {
      auto log_p = std::log(p.get_d());
      primes_and_exps.emplace_back(p, static_cast<uint64_t>(log_bound1 / log_p));
    }

    uint64_t C = (bounds.first / 6) * 6;

//    double sqrt_b2 = std::sqrt(bounds.second);
//    bigint d = 1;
//    for (bigint p = 2; d < sqrt_b2; p = next_prime(p))
//      d *= p;

    for (bigint curve_no = 0; iters < max_iters; ++iters, ++curve_no) {
      auto a   = ring(rng(num));
      auto x_P = ring(rng(num));
      auto y_P = ring(rng(num));

      auto b = (y_P * y_P) - (x_P ^ 3) - (a * x_P);
      if (!aa::check_elliptic_curve_nonsingular(a, b)) {
        continue;
      }

      aa::elliptic_curve_proj_group<murk2::aa::mod_ring> group{c3lt::managed(&ring), a.elem, b.elem};
      auto P = group(x_P.elem, y_P.elem, 1);

      auto Q = P;
      for (auto& [p, k]: primes_and_exps)
        for (uint64_t i = 0; i < k; ++i)
          Q *= p;

      if (group.is_identity(Q.elem)) {
        // catch zero divisors
        if (Q.elem[2])
          return gcd(Q.elem[2], num);
        else
          continue; // TODO: work out what to do here (but it never happens so we don't really care???)
      }

      // https://www.rieselprime.de/ziki/Elliptic_curve_method#Step_2
//      auto Q6 = Q * 6;
//      auto Q_acc = Q * C;
//      auto prod = ring(1);

//      for (uint64_t acc = 0; acc < bounds.second; acc += 6, Q_acc += Q6)
//        prod *= (Q_acc.elem[0] - Q.elem[0]);

//      if (bigint res = gcd(prod.elem, num); res != 1) {
//        std::cout << "FOUND" << std::endl;
//        return res;
//      }


//      std::cout << "Missing out!" << std::endl;

      // TODO: impl continuation
      continue;

//      for (bigint p = bounds.first; p <= bounds.second; p = next_prime(p), ++iters) { //
//        auto Q = P * p;
//        if (auto res = gcd(bigint{Q.elem[0] - P.elem[0]}, num); res != 1) {
//          if (auto ret = gcd(res, num); ret != num)
//            return ret;
//          else
//            // We found a dud
//            goto next_curve;
//        }
//        if (iters >= max_iters)
//          return std::nullopt;
//        auto Q = P * p;
//        auto res = gcd(Q.elem[2], num);
//        if (res != 1)
//          return res;
//      }

//      for (bigint p = next_prime(bounds.first); p <= bounds.second; p = next_prime(p), ++iters) { //
//        if (iters >= max_iters)
//          return std::nullopt;
//        auto Q = P * p;
//        auto res = gcd(Q.elem[2], num);
//        if (res != 1)
//          return res;
//      }

      // We now know that one step along the multiplication would have failed in the affine case, so let's find it
//      auto acc = P;
//      for (bigint i = 1; i < k; ++i, ++iters) {
//      if (iters >= max_iters)
//        return std::nullopt;

//        if (auto res = gcd(bigint{acc.elem[0] - P.elem[0]}, num); res != 1) {
//          if (auto ret = gcd(res, num); ret != num)
//            return ret;
//          else
//            // We found a dud
//            goto next_curve;
//        }
//        acc += P;
//      }

//      throw std::logic_error{"Somehow missed factor in lenstra_ecm"};


next_curve:{}
    }

    return std::nullopt;
  }

  // old SOURCE: https://www.ams.org/journals/mcom/1993-61-203/S0025-5718-1993-1122078-7/S0025-5718-1993-1122078-7.pdf
  // old SOURCE: https://members.loria.fr/PZimmermann/records/ecm/params.html
  //
  // // Proper evaluation would be better, but for now I use Big Table (TM) (Table 3) from the above source
  std::pair<uint64_t, uint64_t> lenstra_ecm_estimate_bounds(bigint const& n) {
    ssize_t n_bits = mpz_sizeinbase(n, 2);

    static const constexpr std::array<uint64_t, 8> b1_tab {{
        256,
        512,
        1024,
        2048,
        4096,
        6000,
        6000,
        65536
    }};

    ssize_t index = n_bits / 16;
    uint64_t b1 = b1_tab[std::clamp<ssize_t>(index, 0, b1_tab.size() - 1)];
    uint64_t b2 = b1 * 100;
//    std::cout << n_bits << ": " << b1 << std::endl;
    return {b1, b2};
  }
}
