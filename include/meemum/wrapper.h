#include <cstddef>
#include <string>
#include <vector>

class MeemumWrapper {
  public:
    void init(const std::string);
    void minimize(double, double, const std::vector<double>);
};
