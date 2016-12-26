#include "libraries.h"

void COOSpMV(int *I_index ,int *J_index,double *A_values,double *b_vector,double *y_vector,int N,int NZ){
	int i;
	
	for (i=0; i<N; i++){
		y_vector[i]=0;
	}
	for (i=0; i<NZ; i++){
		y_vector[I_index[i]]=y_vector[I_index[i]] + A_values[i]*b_vector[J_index[i]];
	}
	
}

int main(int argc, char *argv[]){
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;   
    int i, *I_index, *J_index;
    double *A_values,*b_vector,*y_vector;

    // MPI VARS
    int id, ierr, p, wtime;


    if (argc < 2){
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	} else {
        if ((f = fopen(argv[1], "r")) == NULL){
            exit(1);
        }
    }

    if (mm_read_banner(f, &matcode) != 0){
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) ){
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
	y_vector = (double *) malloc(N * sizeof(double));

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (i=0; i<nz; i++){
        fscanf(f, "%d %d %lg\n", &I_index[i], &J_index[i], &A_values[i]);
        I_index[i]--;  /* adjust from 1-based to 0-based */
        J_index[i]--;
    }

    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
   // for (i=0; i<nz; i++) {
   //     fprintf(stdout, "%d %d %20.19g\n", I_index[i]+1, J_index[i]+1, A_values[i]);
	//}
	//Fill b_vector with random numbers
	for (i=0; i<N; i++){
        b_vector[i]=(double)rand();
    }

    /* Initialize MPI. */
    ierr = MPI_Init ( &argc, &argv );
    /* Get the number of processes. */
    ierr = MPI_Comm_size ( MPI_COMM_WORLD, &p );
    /* Get the individual process ID. */
    ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id );

	// clock_t begin=clock();
    wtime = MPI_Wtime ( );

	COOSpMV(I_index,J_index,A_values,b_vector,y_vector,N,nz);
	
    wtime = MPI_Wtime ( ) - wtime;
    // clock_t end=clock();
	
	// double time_spent=(double)(end-begin) / CLOCKS_PER_SEC;
    printf("Time spent in process: %d is: %f\n", id, wtime);
	// printf("Time spent in serial Sparse Matrix-Vector Multiplication using COO : %lf \n",time_spent);
	/*Debug
		for(i=0;i<N;i++) {	
		printf("Position: %d   %lf \n",i,y_vector[i]); }
	*/

    ierr = MPI_Finalize ( );
	return 0;
}

