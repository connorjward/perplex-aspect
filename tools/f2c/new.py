fortran_header = read_fortran_header()
c_header = []
dtype_lookup = {}

for line in fortran_header:
    if is_dtype_defn(line):
        dtypes = parse_dtype_defn(line)
        for name, dtype in dtypes:
            dtype_lookup[name] = dtype
    elif is_parameter_defn(line):
        line = convert_parameter_defn(line)
        c_header.append(line)
    elif is_common_defn(line):
        line = convert_common_defn(line, dtype_lookup)
        c_header.append(line)

write_c_header(c_header)

