#!/bin/bash                                                                                         
casename='turbPipe'
noPrcs=8    #number of processors

echo $casename > SESSION.NAME
echo $PWD/ >> SESSION.NAME

if [ -f "logfile" ]; then
   rm logfile
fi
mpirun -np $noPrcs nek5000 >>logfile&
