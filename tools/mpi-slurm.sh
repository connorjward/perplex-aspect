#!/bin/csh

# simple template job script for submission of an MPI job 
# on Hamilton

# sbatch -p parallel_queue  name_of_job_script

# or with overriding the number of tasks as option

# sbatch -p parallel_queue -n number_of_tasks  name_of_job_script

# If successful SLURM will return a jobID for this job which can be
# used to query its status.

#############################################################################

## All lines that start with #SBATCH will be processed by SLURM.
## Lines in this template script that have white space between # and SBATCH 
## will be ignored. They provide examples of further options that can be
## activated by deleting the white space and replacing any text after the 
## option.

## By default SLURM uses as working directory the directory from where the
## job script is submitted. To change this the standard Linux cd command has
## to be used.

## Name of the job as it will appear when querying jobs with squeue (the
## default is the name of the job script)

#SBATCH  -J  aspect

## By default SLURM combines the standard output and error streams in a single
## file based on the jobID and with extension .out
## These streams can be directed to separate files with these two directives

##SBATCH  -o  file_name.o%j
##SBATCH  -e  err_file_name.e%j

## where SLURM will expand %j to the jobID of the job.

## Request email to the user when certain type of events occur to 
## the job

##SBATCH  --mail-type=ALL

## where <type> can be one of BEGIN, END, FAIL, REQUEUE or ALL,
## and send to email address

##SBATCH  --mail-user  jeroen.van-hunen@durham.ac.uk

## The default email name is that of the submitting user as known to the system.

## Specify project or account name (currently not required).
##SBATCH -A ITSall

#############################################################################

## This job requests number_of_tasks MPI tasks (without OpenMP)

#SBATCH  -n  8 

# Request submission to a queue (partition) for parallel jobs

#SBATCH  -p  par6.q

module purge
module load slurm/current
## Load any modules required here

module load module-git gcc/4.9.1 cmake/3.6.2 lapack/gcc/3.5.0 zlib/gcc/1.2.7 sge/current openmpi/gcc/2.1.1 gsl/gcc/64/1.15

## Execute the MPI program

#mpirun -np 8 ./aspect $1
mpirun -np 8 aspect $1
