'''
docstring
'''

import os
import re

F2C_DTYPES = {'integer': 'int',
              'double precision': 'double',
              'logical': 'int',
              'character': 'int'}

def read_fortran_header():
    # get Perple_X directory
    try:
        perplex_dir = os.environ['PERPLEX_DIR']
    except KeyError:
        print('ERROR: PERPLEX_DIR environment variable not defined.')
        sys.exit(1)

    # read file
    fortran_header = []
    with open('{}/perplex_parameters.h'.format(perplex_dir)) as f:
        for line in f:
            # ignore comments and blank lines
            if ( line.isspace() 
                    or line.startswith('c') 
                    or line.startswith('!') ):
                continue

            # remove leading and trailing whitespace
            line = line.strip()

            # handle line continuations
            if line.startswith('*'):
                # remove asterisk and whitespace
                line = line.lstrip('* ')

                fortran_header[-1] += line
            else:
                fortran_header.append(line) 

    return fortran_header

def write_c_header(c_header):
    with open('../src/_parameters.h', 'w') as f:
        f.write('\n'.join(c_header))

def is_dtype_defns(line):
    return re.match('integer|double precision|logical|character', line)

def parse_dtype_defns(line):
    # extract dtype
    dtype, variables = re.match('((?:integer'
                                '|double precision'
                                '|logical'
                                '|character)(?:\*\d+)?)'
                                '\s+'
                                '(.+)'
                                , line).groups()

    names = []
    for variable in re.split('\s*,\s*', variables):
        name, shape = re.match('(\w+)(\*\d+)?', variable).groups()

        if shape is not None:
            dtype += shape
        names.append(name)

    return dtype, names

def is_parameter_defns(line):
    return line.startswith('parameter')

def convert_parameter_defns(line):
    c_header_lines = []

    #remove parameter keyword and outer brackets
    line = re.sub('parameter\s*\(', '', line)
    line = re.sub('\)$', '', line)
            
    # identify parameters
    name_value_pairs = re.split('\s*,\s*', line)

    for pair in name_value_pairs:
        name, value = re.split('\s*=\s*', pair)
        c_header_lines.append('#define {} {}'.format(name, value))

    return c_header_lines

def is_common_defn(line):
    return line.startswith('common')

def convert_common_defn(line, dtype_lookup):
    c_header_line = ['extern struct {']

    #remove common keyword and forward slash
    line = re.sub('common\s*/', '', line)

    # identify name and variables 
    name, variables_string = re.split('\s*/\s*', line)
    entries = re.findall('[a-z0-9]+\([a-z0-9,:]+\)|[a-z0-9]+',
                           variables_string)

    for entry in entries:
        name, *shape = re.findall('[0-9a-z]+', entry)
        dtype = dtype_lookup[name]
        match = re.fullmatch('(\w+)\*(\d+)', dtype)
        if match:
            dtype = match.group(1)
            shape.append(match.group(2))

        # convert dtype to c
        dtype = F2C_DTYPES[dtype]

        if shape:
            c_header_line.append('\t{} {}[{}];'.format(dtype, name,
                                 ']['.join(shape)))
        else:
            c_header_line.append('\t{} {};'.format(dtype, name))

    c_header_line.append('}} {}_'.format(name))
            
    return '\n'.join(c_header_line)


fortran_header = read_fortran_header()
c_header = []
dtype_lookup = {}

for line in fortran_header:
    if is_dtype_defns(line):
        dtype, names = parse_dtype_defns(line)
        for name in names:
            dtype_lookup[name] = dtype
    elif is_parameter_defns(line):
        line = convert_parameter_defns(line)
        c_header.extend(line)
    elif is_common_defn(line):
        line = convert_common_defn(line, dtype_lookup)
        c_header.append(line)

write_c_header(c_header)

