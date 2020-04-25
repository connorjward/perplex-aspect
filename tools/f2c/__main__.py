'''
TODO write description here
TODO include copyright info about Bob Myhill
'''

import os
import sys

from compiler import Compiler

def to_header(parameters, commons):
    for parameter in parameters:
        print(parameter.header_string)

    for common in commons:
        print(common.header_string)


DTYPES = {'integer': 'int', 
          'double precision': 'double',
          'logical': 'bool',
          'character': 'char'}

# get Perple_X directory
try:
    perplex_dir = os.environ['PERPLEX_DIR']
except KeyError:
    print('ERROR: PERPLEX_DIR environment variable not defined.')
    sys.exit(1)

compiler = Compiler.from_file('{}/perplex_parameters.h'.format(perplex_dir))

compiler.parse() 
compiler.compile(None) 

