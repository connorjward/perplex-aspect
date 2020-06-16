#include <string>
#include <vector>

namespace perplex 
{
  /**
   * Initialize Perple_X.
   *
   * @param filename The Perple_X problem definition file.
   */
  void init(const std::string filename);

  /**
   * Perform the minimization.
   *
   * @param pressure    The pressure (Pa).
   * @param temperature The temperature (K).
   */
  void minimize(const double pressure, const double temperature);
} 
