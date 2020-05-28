#include <cstddef>
#include <string>
#include <vector>

class MeemumWrapper {
  public:
    MeemumWrapper(const std::string);
    void minimize(double, double, std::vector<double>);
};
