import shutil

# This python code takes perplex_parameters.h
# (the Fortran header file for the PerpleX executables)
# wrapper_definitions.txt and a change_log.txt file containing a list of changes
# to the perplex_parameters.h file, and creates a C header file which can be
# used to create a wrapper for meemum

# Bob Myhill, 2018/06/22

# Comments at start of file
string=("/*\n" 
        " * Header file for C/C++ wrapper for Perple_X-meemum.\n"
        " * - Perple_X-meemum by Jamie Connolly (see http://www.perplex.ethz.ch/)\n"
        " * - This header was created automatically using create_C_header.py.\n"
        " * - It is designed for use with the C wrapper written by Lars Kaislaniemi (lars.kaislaniemi@iki.fi)\n"
        " */\n"
        "\n"
        "/* Changes:\n")

# Write change log
with open("change_log.txt") as f:
    changes = f.readlines()

for line in changes:
    string+=" * - "+line

# Add comment about parameters
string+=("*/\n"
         "\n"
         "/* PerpleX constant parameters (see explanation below or in Perple_X\n"
         " * sources perplex_parameters.h. These MUST HAVE the same values as\n"
         " * in Perple_X code (perplex_parameters.h), otherwise resulting\n"
         " * in array size mismatch and overall havoc.\n"
         " */\n")

# Copy fortran header into include directory
shutil.copy2('../perplex_parameters.h', './perplex_parameters.h')

# Read parameters from fortran header
with open("perplex_parameters.h") as f:
    params = f.readlines()

# Process parameters
params=[x.strip() for x in params
        if x[0]!="!"
        and x[0]!="c"
        and x.strip()[0:7]!="integer"
        and len(x.strip()) != 0]

merged_params = []
for i in range(len(params)):
    if params[i][0]=="*":
        merged_params[-1] += params[i][1:].replace(" ", "")
    else:
        merged_params.append(params[i])

merged_params = [p[:-1].replace("parameter (", "").split(",") for p in merged_params]
merged_params = [item for sublist in merged_params for item in sublist]

param_list = [p.split("=")[0] for p in merged_params]
for param in param_list:
    merged_params = [p.replace("p_"+param, param).replace(param, "p_"+param)
                     for p in merged_params]

for p in merged_params:
    string+="#define "+p.replace("=", " ")+"\n"
    
string+="\n"

# Add wrapper definitions
with open("wrapper_definitions.txt") as f:
    defs = f.readlines()

for definition in defs:
    string+=definition

# Write string to C header file
with open("perplex_c.h", "w") as f:
    f.write(string)
