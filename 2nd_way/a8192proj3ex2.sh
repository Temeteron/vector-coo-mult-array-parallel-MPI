#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
mpirun prog.out bcs8192.mtx
