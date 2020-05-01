#include <iostream>
#include "meemum_wrapper.hpp"

int main()
{
  std::cout << "Starting c_init()..." << std::endl;
  c_init();

  std::cout << "Starting c_minimize()..." << std::endl;
  c_minimize();
}
