#include <iostream>
#include "meemum_wrapper.hpp"

int main() {
  c_init();
  c_minimize();

  double components [12];

  update_components(components);

  for (int i = 0; i < 12; i++) {
    //std::cout << get_component_amount(&i) << std::endl;
    std::cout << components[i] << std::endl;
  }
}
