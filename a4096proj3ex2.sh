#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
mpirun prog.out bcs4096.mtx
