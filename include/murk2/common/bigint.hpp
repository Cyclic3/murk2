#pragma once

#include <gmpxx.h>
#include <iostream>

namespace murk2 {
  class bigint : public mpz_class {
  public:
    inline int sign() const noexcept { return get_mpz_t()->_mp_size; }

  public:
    constexpr operator bool() const noexcept {
      return get_mpz_t()->_mp_size != 0;
    }
    constexpr operator bool() noexcept {
      return get_mpz_t()->_mp_size != 0;
    }
    inline operator mpz_ptr() { return get_mpz_t(); }
    inline operator mpz_srcptr() const { return get_mpz_t(); }


  public:
    bigint() = default;
    template<typename... Args>
    inline bigint(Args&&... args) : mpz_class{std::forward<Args>(args)...} {}
  };

  inline std::ostream& operator<<(std::ostream& os, bigint const& i) {
    int base;

    switch (os.flags() & std::ios_base::basefield) {
      case std::ios_base::hex: base = 16; break;
      case std::ios_base::oct: base = 8; break;
      default: base = 10;
    }

    return os << i.get_str(base);
  }

  template<typename F, typename... Args>
  inline bigint do_gmp(F f, Args&&... args) {
    bigint ret;
    f(ret, std::forward<Args>(args)...);
    return ret;
  }

//  using bigint = mpz_class;//boost::multiprecision::mpz_int;
}
