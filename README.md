# vector-coo-mult-array-parallel-MPI
Multiplication of a vecotor n*n in Coordinate Format(COO) with an array in Parallel with MPI

Seperation of jobs to N/p block of lines, N=size of vector(N*N), p=number of processes

Makefile commands

	make: compiles code
	make clean: deletes files
	make run: executes prog without args
	make ARGS=my__vector.mtx run: executes prog with arg the file which contains a vector

Run with more processes

	mpirun -np "numberOfProcesses" ./prog.out bcs4096.mtx

Use of mtx vectors as input, sparse/diagonial
