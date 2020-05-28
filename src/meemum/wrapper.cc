#include <meemum/wrapper.h>
#include "ftoc.h"

void MeemumWrapper::init(const std::string filename) {
  // TODO: set filename before c_init
  c_init(filename.c_str());
}

void MeemumWrapper::minimize(const double pressure, 
                             const double temperature,
			     const std::vector<double> composition) {
  // set the temperature and pressure
  c_set_temperature(temperature);
  c_set_pressure(pressure);

  // set the composition
  for (size_t i = 0; i < composition.size(); i++)
    c_set_composition_component(i, composition[i]);

  // run the minimization
  c_minimize();
}
