#include <murk2/common/simulation.hpp>

#include <xxhash.h>

namespace murk2 {
  uint64_t simulation_hash(bigint const& x) {
    std::string str = x.get_str();

    return XXH64(str.data(), str.size(), 0);
  }
}
