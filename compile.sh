#!/bin/sh

# please execute this file with "./compile.sh" to move to bin folder
# after compilation has been succeed (please note there is extra dot before file name)

module purge
module load intel/2018-mpi
module load mpi/intel/mpich-1.2.7p1
make -f MakefileEdison clean
make -f MakefileEdison
