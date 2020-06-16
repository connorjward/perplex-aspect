#include <cassert>

#include <perplex/state.h>
#include <perplex/utils.h>

namespace perplex
{
  namespace utils
  {
    double convert_pascals_to_bar(const double pressure_in_pascals) {
      return pressure_in_pascals / 1e6;
    }

    double convert_bar_to_pascals(const double pressure_in_bar) {
      return pressure_in_bar * 1e6;
    }
  }
}
