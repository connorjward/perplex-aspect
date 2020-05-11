#include <iostream>
#include "meemum/calcs.hpp"
#include "meemum/props.hpp"

int main() {
    const size_t id { 1 };

    meemum::init("khgp");

    char* name { new char[20] };

    for (size_t id = 1; id <= meemum::props::n_soln_models(); id++) {
	meemum::props::abbr_soln_name(&id, name);
	std::cout << name << std::endl;
    }
}
