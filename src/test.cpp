#include <iostream>
#include "meemum_wrapper.hpp"

int main() {
    double pressure = 100;
    double temperature = 500;

    init();
    minimize(&pressure, &temperature);

    int comp_id = 2;
    char comp_name [6];
    double components [12];

    get_component_name(&comp_id, comp_name);
    get_component_amount(&comp_id, components);

    std::cout << get_n_components() << std::endl;
    std::cout << comp_name << std::endl;

    for (int i = 0; i < 12; i++) {
	std::cout << components[i] << std::endl;
    }
}
