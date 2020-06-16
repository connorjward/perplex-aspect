#include <fcntl.h>
#include <unistd.h>

#include <perplex/main.h>
#include <perplex/utils.h>
#include "c_interface.h"

namespace 
{
  const int disable_stdout() {
    // flush stdout
    fflush(stdout);

    // get file descriptors
    const int stdout_descriptor = dup(1);
    const int null_descriptor = open("/dev/null", O_WRONLY);

    // reassign stdout to /dev/null
    dup2(null_descriptor, 1);
    close(null_descriptor);

    return stdout_descriptor;
  }

  void enable_stdout(const int stdout_descriptor) {
    // flush stdout
    fflush(stdout);

    // reassign descriptor
    dup2(stdout_descriptor, 1);
    close(stdout_descriptor);
  }
}

namespace perplex 
{
  void init(const std::string filename) 
  {
    // disable stdout to prevent Perple_X dominating stdout
    const int fd = disable_stdout();
    solver_init(filename.c_str());
    enable_stdout(fd);
  }

  void minimize(const double pressure, const double temperature)
  {
    solver_set_pressure(utils::convert_pascals_to_bar(pressure));
    solver_set_temperature(temperature);

    // disable stdout to prevent Perple_X dominating stdout
    /* const int fd = disable_stdout(); */
    solver_minimize();
    /* enable_stdout(fd); */
  }
}
