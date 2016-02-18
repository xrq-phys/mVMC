#!/bin/bash
#------ pjsub option --------#
#PJM -L "rscgrp=F12"
#PJM -L "node=2x3x2"
#PJM --mpi "proc=192"
#PJM -L "elapse=1440:00"
#PJM -j
#------- Program execution -------# 
export OMP_NUM_THREADS=1
fipp -C -I call,hwm -d test mpiexec ./vmc.out xnamelist.def zqp_opt.dat
