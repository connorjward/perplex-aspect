'''
docstring
'''

import re

F2C_DTYPES = {'integer': 'int',
              'double precision': 'double',
              'logical': 'int',
              'character': 'char'}

class Compiler:
    '''docstring'''

    def __init__(self, source):
        self.source = source

        self.variables = {}
        self.parameters = []
        self.commons = []

    def parse(self):
        def _parse_variables(line, dtype):
            # identify variables
            variable_names = re.findall('[0-9a-z]+', line)

            # store variables
            for name in variable_names:
                variable = Variable(name, dtype)

                self.variables[name] = variable

        def _parse_parameters(line):
            # identify parameters
            name_value_pairs = re.split('\s*,\s*', line)

            # store parameters
            for pair in name_value_pairs:
                name, value = re.split('\s*=\s*', pair)
            variable = self.variables[name]

            self.parameters.append(Parameter(variable, value))

        def _parse_common(line):
            name, data_string = re.split('\s*/\s*', line)
            variables = re.findall('[a-z0-9]+\([a-z0-9,:]+\)|[a-z0-9]+',
                                   data_string)

            scalars = []
            arrays = []
            for var_string in variables:
                # check if var is a scalar
                if re.fullmatch('[0-9a-z]+', var_string):
                    var = self.variables[var_string]
                    scalars.append(var)
                # else var is an array
                else:
                    var_name, *shape = re.findall('[0-9a-z]+', var_string)
                    var = self.variables[var_name]
                    arrays.append(Array(var, shape))

            self.commons.append(Common(name, scalars, arrays))

        for line in self.source:
            for dtype in ['integer', 'double precision', 
                    'logical', 'character']:
                if line.startswith(dtype):
                    line = re.sub('{}\s*'.format(dtype), '', line) 
                    _parse_variables(line, dtype)
                    break

            if line.startswith('parameter'):
                line = re.sub('parameter\s*\(', '', line)
                line = re.sub('\)$', '', line)
                _parse_parameters(line)
            elif line.startswith('common'):
                line = re.sub('common\s*/', '', line)
                _parse_common(line)

    def compile(self, filename):
        compiled = []
        for parameter in self.parameters:
            compiled.append(parameter.header_string)
        for common in self.commons:
            compiled.append(common.header_string)

        print('\n'.join(compiled))


    @classmethod
    def from_file(cls, filename):
        source = []

        with open(filename) as f:
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

                    source[-1] += line
                else:
                    source.append(line) 

        return cls(source)


class Parameter:
    '''docstring'''

    def __init__(self, variable, value):
        self.variable = variable
        self.value = value

    @property
    def header_string(self):
        return '#define {} {}'.format(self.variable.name, self.value)


class Common:
    '''docstring'''

    def __init__(self, name, scalars, arrays):
        self.name = name
        self.scalars = scalars
        self.arrays = arrays

    @property
    def header_string(self):
        header = ['extern struct {']
        for scalar in self.scalars:
            header.append(scalar.header_string)
        for array in self.arrays:
            header.append(array.header_string)
        header.append('}} {}_'.format(self.name))

        return '\n'.join(header)


class Array:
    '''docstring'''

    def __init__(self, variable, shape):
        self.variable = variable
        self.shape = shape

    @property
    def header_string(self):
        return '{} {}[{}];'.format(self.variable.cdtype, self.variable.name, ']['.join(self.shape))


class Variable:
    '''docstring'''

    def __init__(self, name, dtype):
        self.name = name
        self.dtype = dtype
        self.cdtype = F2C_DTYPES[dtype]

    @property
    def header_string(self):
        return '{} {};'.format(self.dtype, self.name)

