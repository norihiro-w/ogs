
#!/bin/bash

ogs="~/ogs/ogs6-thmc/BuildPETScGccDebug/bin/ogs6 -i line --ls petsc"

#normal run
mpirun -np 2 $ogs

#total view
#mpirun -tv -np 2 $ogs

# valgrid
#mpirun -np 2 valgrind --tool=memcheck -q --log-file=valgrind.log.%p $ogs

