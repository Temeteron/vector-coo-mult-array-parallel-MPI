#include "libraries.h"
/*Greasidis Dimitrios AEM : 1624
 *Kontogiannis Dimitrios AEM : 1565
 *
 *Parallel implementation of sparse matrix-vector multiplication using MPI in COO format.
 *Steps of Algorithm:
 *1)Root process reads matrix from the .mtx file using the Matrix Market library, then 
 * the multiplication vector is generated and preprocessed to match the COO format in 
 * an attempt to avoid excessive all-to-all communication
 *2)Node Setup : Data distribution by row-block of N/p rows to each computation node
 * in MPI_COMM_WORLD is done with one-to-all personalized communication,with the use 
 * of MPI_Scatterv.Each node should have a part of I_index ,J_index,A_values and 
 * b_vector at the end of the setup.
 *3)Multiplication is done locally on each without without any ongoing communication 
 * thanks to the initial preprocessing of b.At the end of their computation each 
 * processing node should have N/p elements of the resulting y vector.
 *4)Result reconstruction is done by using all-to-one personalized communication of 
 * fixed message length with the source as destination.This is achieved with the usage
 * of MPI_Gather.
 *5)The timing of the program consists of the computation and the reconstruction 
 * (without considering the setup and the initialization).
 * 
 * Compile: make
 * Run: mpirun -np numberofprocessors ./prog.out matrixmarketfilename.mtx
 */
int main(int argc, char *argv[])
{
    int ret_code,rank,size;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz,sum,id;   
    int i, *I_index, *J_index,*sendcounts,*displs,*rec_buf1,*rec_buf2;
    double *A_values,*b_vector,*y_vector,*btsend,*rec_buf3,*rec_buf4,*result,tstart,tend,time_spent;
	
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	
	
	sendcounts = malloc(sizeof(int)*size);
    displs = malloc(sizeof(int)*size);
	
	if (rank==0){

    if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}
    else    
    { 
        if ((f = fopen(argv[1], "r")) == NULL) 
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    /* reseve memory for matrices */
	//i_index array
    I_index = (int *) malloc(nz * sizeof(int));
    J_index = (int *) malloc(nz * sizeof(int));
    A_values = (double *) malloc(nz * sizeof(double));
	
	b_vector = (double *) malloc(N * sizeof(double));
	btsend = (double *) malloc(nz * sizeof(double));
	result = (double *) malloc(N * sizeof(double));
    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
//We work as example with the transpose of the MatrixMarket Sparse Matrix
    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n" , &J_index[i],&I_index[i], &A_values[i]);
        I_index[i]--;  /* adjust from 1-based to 0-based */
        J_index[i]--;
    }

    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    //for (i=0; i<10; i++) {
    //    fprintf(stdout, "%d %d %20.19g\n", I_index[i], J_index[i], A_values[i]);
	//}
	//Fill b_vector with random numbers
	for (i=0; i<N; i++)
    {
        b_vector[i]=(double)rand();
    }
    
    for(i=0; i<size; i++){
		sendcounts[i]=0;
		displs[i]=0;
	}
	//Preparation for data distribution and preprocessing
	sum=0;
	for(id=0; id<size; id++){
		for(i=0; i<nz; i++){
			if(I_index[i]>=id*N/size && I_index[i]<(id*N/size + N/size)){
				sendcounts[id]++;
				btsend[i]=b_vector[J_index[i]];
			}
		}
		
		displs[id]=displs[id]+sum;
		sum += sendcounts[id];
		
	}
    
	}
	//Prepare Collective Communication
	MPI_Bcast(sendcounts,size,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(displs,size,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&nz,1,MPI_INT,0,MPI_COMM_WORLD);
	
	y_vector = (double *) malloc(N * sizeof(double));
	
	rec_buf1 = (int *) malloc(sendcounts[rank] * sizeof(int));
    rec_buf2 = (int *) malloc(sendcounts[rank] * sizeof(int));
    rec_buf3 = (double *) malloc(sendcounts[rank] * sizeof(double));
	rec_buf4 = (double *) malloc(sendcounts[rank] * sizeof(double));
	
		
	MPI_Scatterv(I_index, sendcounts, displs, MPI_INT, rec_buf1, sendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(J_index, sendcounts, displs, MPI_INT, rec_buf2, sendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(A_values, sendcounts, displs, MPI_DOUBLE, rec_buf3, sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(btsend, sendcounts, displs, MPI_DOUBLE, rec_buf4, sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	//Make sure the setup has finished
	MPI_Barrier(MPI_COMM_WORLD);
	
    //Do the multiplication
	
	tstart=MPI_Wtime();
    for(i=rank*N/size; i<(rank*N/size+N/size); i++){
		y_vector[i]=0;
	}
    for ( i = 0; i < sendcounts[rank]; i++) {
        y_vector[rec_buf1[i]]=y_vector[rec_buf1[i]]+rec_buf3[i]*rec_buf4[i];
	}

	
	tend=MPI_Wtime();
	time_spent=tend-tstart;
	
	printf("Time spent computing rank %d COO : %lf \n",rank,time_spent);
	
	tstart=MPI_Wtime();
	//Result reconstruction
	MPI_Gather(y_vector+rank*N/size,N/size,MPI_DOUBLE,result,N/size,MPI_DOUBLE,0,MPI_COMM_WORLD);
	tend=MPI_Wtime();
	time_spent=tend-tstart;
	printf("Time spent sending rank %d COO : %lf \n",rank,time_spent);
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	MPI_Finalize();
	free(sendcounts);
    free(displs);
	free(rec_buf1);
    free(rec_buf2);
	free(rec_buf3);
    free(rec_buf4);
	free(y_vector);
	

	return 0;
}

