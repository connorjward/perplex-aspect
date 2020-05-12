#include <iostream>
#include "meemum/calcs.hpp"
#include "meemum/props.hpp"

int main() {
    meemum::init("khgp");

    size_t id = 2;
    char* name = meemum::props::abbr_soln_name(id);
    std::cout << name << std::endl;

    for (size_t id = 1; id <= meemum::props::n_soln_models(); id++) {
	std::cout << meemum::props::abbr_soln_name(id) << std::endl;
    }
	std::cout << name << std::endl;
}
