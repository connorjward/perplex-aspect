add_custom_command(OUTPUT src/parameters.h
                   COMMAND $(python?) build_parameters_header.py?
