#include <iostream>
#include "meemum/interface.hpp"

int main() {
    const size_t id { 1 };

    meemum::init("khgp");

    char* name { new char[20] };
    meemum::load_abbr_soln_name(&id, name);

    std::cout << "Start" << std::endl;
    std::cout << name << std::endl;
    std::cout << meemum::abbr_soln_name(&id) << std::endl;
    std::cout << "End" << std::endl;
}
