#include "libraries.h"
/*Greasidis Dimitrios AEM : 1624
 *Kontogiannis Dimitrios AEM: 1565
 * 
 *Load balanced parallel implementation of sparse matrix-vector multiplication 
 *using MPI in COO format.
 *Steps of Algorithm:
 *1)Root process reads matrix from the .mtx file using the Matrix Market library, then 
 * the multiplication vector is generated and preprocessed to match the COO format in 
 * an attempt to avoid excessive all-to-all communication
 *2)Node Setup : Data distribution by row-blocks of approximately NZ/p elements to each computation
 * node in MPI_COMM_WORLD is done with one-to-all personalized communication,with the usage 
 * of MPI_Scatterv.Each node should have a part of I_index ,J_index,A_values and 
 * b_vector at the end of the setup.
 *3)Multiplication is done locally on each without any ongoing communication 
 * thanks to the initial preprocessing of b.At the end of their computation each 
 * processing node should have N/p elements of the resulting y vector.
 *4)Result reconstruction is done by using all-to-one personalized communication of 
 * variable message length with the source as destination.This is achieved with the usage
 * of MPI_Gatherv.
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
    int M, N, nz,sum,id,k,count,z,resumed,total,sum2;   
    int i, *I_index, *J_index,*sendcounts,*displs,*rec_buf1,*rec_buf2,*grcount,*recvcounts,*recdispls;
    double *A_values,*b_vector,*y_vector,*btsend,*rec_buf3,*rec_buf4,*result,tstart,tend,time_spent,ii;
	
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
	grcount = (int *) malloc(N * sizeof(int));
	for(i=0; i<N; i++){
			grcount[i]=0;
	}
	
    I_index = (int *) malloc(nz * sizeof(int));
    J_index = (int *) malloc(nz * sizeof(int));
    A_values = (double *) malloc(nz * sizeof(double));
	
	b_vector = (double *) malloc(N * sizeof(double));
	btsend = (double *) malloc(nz * sizeof(double));
	result = (double *) malloc(N * sizeof(double));
	recvcounts = (int *) malloc(size * sizeof(int));
	recdispls = (int *) malloc(size * sizeof(int));
	for(i=0; i<N; i++){
		result[i]=0.0;
	}
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
		recdispls[i]=0;
	}
	
	i=0;
	for(z=0; z<N; z++){
		 while(I_index[i]==z){
			grcount[z]=grcount[z]+1;
			i++;
		}
	//	printf("Row %d HAS %d elements \n",z,grcount[z]);
	}
	//Partition from COO format in approximately NZ/p parts so that each process is load balanced
	//The preparation for the distribution is done here.
	i=0;
	sum=0;
	sum2=0;
	resumed=0;
	z=0;
	for(id=0; id<size; id++){
		count=(int)ceil((float)nz/(float)size);
		//printf("COUNT : %d \n",count);
		while(count>=grcount[z]){
				sendcounts[id]=sendcounts[id]+grcount[z];
				count=count-grcount[z];
				while (I_index[i]==z){
					btsend[i]=b_vector[J_index[i]];
					//printf( "Index [ %d ] gets btsend [ %d ] = %lf ",i,i,b_vector[J_index[i]]);
					i++;
				}
				z=z+1;
				
			}
		recvcounts[id]=z-1-resumed;
		recdispls[id]=recdispls[id]+sum2;
		sum2+=recvcounts[id];
		resumed=z-1;
		displs[id]=displs[id]+sum;
		sum += sendcounts[id];
		}
		//printf(" %d , %d \n",sum ,nz );
		if (sum!=nz){
			for(i=sum; i<nz; i++){
				sendcounts[size-1]++;
				btsend[i]=b_vector[J_index[i]];
			}
		}
	total=0;
	 for ( i = 0; i < size; i++) {
		 total=total+recvcounts[i];
         //  printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
       }
	recvcounts[size-1]=recvcounts[size-1] + (N-total);
	for ( i = 0; i < size; i++) {
	//printf("Lines at rank %d  = %d | %d \n",i,recvcounts[i],recdispls[i]);
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
	
	
	//Distribute Data	
	MPI_Scatterv(I_index, sendcounts, displs, MPI_INT, rec_buf1, sendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(J_index, sendcounts, displs, MPI_INT, rec_buf2, sendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(A_values, sendcounts, displs, MPI_DOUBLE, rec_buf3, sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(btsend, sendcounts, displs, MPI_DOUBLE, rec_buf4, sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	//Make sure the setup has finished
	MPI_Barrier(MPI_COMM_WORLD);
    //Do the multiplication
	
	//Computation After Setup
	tstart=MPI_Wtime();
	for (i=rec_buf1[0]; i<=rec_buf1[sendcounts[rank]-1]; i++){
		y_vector[i]=0.0;
		
	}
    
    for ( i = 0; i < sendcounts[rank]; i++) {
        y_vector[rec_buf1[i]]=y_vector[rec_buf1[i]]+rec_buf3[i]*rec_buf4[i];
		
	}
	for ( i=rec_buf1[0]; i<=rec_buf1[sendcounts[rank]-1]; i++) {
	//printf("RANK %d has y[%d] = %lf \n",rank,i,y_vector[i]);
	}
	
	tend=MPI_Wtime();
	time_spent=tend-tstart;
	//clock_t end=clock();
	//double time_spent=(double)(end-begin) / CLOCKS_PER_SEC;
	printf("Time spent computing rank %d COO : %lf \n",rank,time_spent);
	//Result Retrieval
	tstart=MPI_Wtime();
	//Result reconstruction
	MPI_Gatherv(&y_vector[rec_buf1[0]],rec_buf1[sendcounts[rank]-1]-rec_buf1[0],MPI_DOUBLE,result,recvcounts,recdispls,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
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

