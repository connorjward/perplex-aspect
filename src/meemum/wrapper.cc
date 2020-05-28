#include <meemum/wrapper.h>

extern "C" {
  void c_init(const char*);
  void c_minimize();

  void c_set_pressure(double);
  void c_set_temperature(double);

  double c_get_composition_component(size_t);
  void c_set_composition_component(size_t, double);
}

MeemumWrapper::MeemumWrapper(const std::string filename) {
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
