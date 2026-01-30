/**
 * This version is stamped on May 10, 2016
 *
 * Contact:
 *   Louis-Noel Pouchet <pouchet.ohio-state.edu>
 *   Tomofumi Yuki <tomofumi.yuki.fr>
 *
 * Web address: http://polybench.sourceforge.net
 */
/* gesummv.c: this file is part of PolyBench/C */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* SuiteSparse file to read Rutherford-Boeing format */
#include <suitesparse/RBio.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
#define DATA_TYPE float
#define DATA_PRINTF_MODIFIER "%0.2f "

/* Array initialization. */
static
void init_array(char * rbpath,
                SuiteSparse_long * nrow,
		SuiteSparse_long * ncol,
                SuiteSparse_long ** Ai,
                SuiteSparse_long ** Ap,
                double ** A,
                double ** x,
                double ** y )
{
  int i, j;
  SuiteSparse_long ret;
  char title[73];
  char key[9];
  char mtype[4];
  SuiteSparse_long mkind;
  SuiteSparse_long skind;
  SuiteSparse_long asize;
  SuiteSparse_long znz;
  double * Az;
  SuiteSparse_long *Zp;
  SuiteSparse_long *Zi;

  ret = RBread( rbpath, 1, 0, title, key, mtype, nrow, ncol, &mkind, &skind, &asize, &znz, Ap, Ai, A, &Az, &Zp, &Zi );
//  printf( "Returned value: %d\n", ret );
  int nnz = (*Ap)[*ncol];
//  printf( "Matrix is %d X %d with %d NNZs\n", *nrow, ncol, nnz );
  SuiteSparse_long * CSRAi;
  posix_memalign(&CSRAi, 2*1024*1024, sizeof(SuiteSparse_long)*nnz);
  SuiteSparse_long * CSRAp;
  posix_memalign(&CSRAp, 2*1024*1024, sizeof(SuiteSparse_long)*(*nrow+1));
  double * CSRdata;
  posix_memalign(&CSRdata, 2*1024*1024, sizeof(double)*nnz);
  // Initialize pointer array to zero
  for( int i = 0; i < *nrow; ++i ) {
      CSRAp[i] = 0;
  }
  // Add one to each entry in pointer array for each value of CSC index array
  for( int i = 0; i < nnz; ++i ) {
      CSRAp[(*Ai)[i]]++;
  }
  // Cumsum to obtain Ap
  for( int i = 0, cumsum = 0; i < *nrow; ++i ) {
      int temp = CSRAp[i];
      CSRAp[i] = cumsum;
      cumsum += temp;
  }
  CSRAp[*nrow] = nnz;

  for( int col = 0; col < *ncol; ++col ) {
      for( int jj = (*Ap)[col]; jj < (*Ap)[col+1]; jj++ ) {
          int row = (*Ai)[jj];
          int dest = CSRAp[row];

          CSRAi[dest] = col;
          CSRdata[dest] = (*A)[jj];

          CSRAp[row]++;
      }
  }

  for( int row= 0, last = 0; row<= *nrow; ++row ) {
      int temp = CSRAp[row];
      CSRAp[row] = last;
      last = temp;
  }

  *x = (double *)malloc(sizeof(double)*(*ncol));
  *y = (double *)malloc(sizeof(double)*(*nrow));

  for (i = 0; i < *ncol; i++)
    {
      (*x)[i] = (double)( i % *ncol) / *ncol;
    }

  for (i = 0; i < *nrow; i++)
    {
      (*y)[i] = 0.0;
    }

  free(*A);
  *A = CSRdata;
  free(*Ai);
  *Ai=CSRAi;
  free(*Ap);
  *Ap=CSRAp;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(y,n,n))

{
  int i;

  POLYBENCH_DUMP_START;
  POLYBENCH_DUMP_BEGIN("y");
  for (i = 0; i < n; i++) {
    if (i % 20 == 0) fprintf (POLYBENCH_DUMP_TARGET, "\n");
    fprintf (POLYBENCH_DUMP_TARGET, DATA_PRINTF_MODIFIER, y[i]);
  }
  POLYBENCH_DUMP_END("y");
  POLYBENCH_DUMP_FINISH;
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static __attribute__((noinline))
void kernel_gesummv(SuiteSparse_long n,
		    SuiteSparse_long * indices,
                    SuiteSparse_long * indptr,
                    double * A,
                    double * x,
		    double * y )
{
  int i, j,t;

#pragma omp parallel for private(j)
  for( i = 0; i < n; ++i ) {
      y[i] = 0.0;
      for( j = indptr[i]; j < indptr[i+1]; ++j )
	{
          y[i] = A[j] * x[indices[j]] + y[i];
#ifdef  TEST_DEBUG
          printf( "%d:%d:%d:%lf:%lf:%lf\n", i, j, indices[j], A[j], x[indices[j]], y[i] );
#endif
	}
    }

}


int polybench_ntasks = 1;
int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  SuiteSparse_long nrow, ncol;

  /* Variable declaration/allocation. */
  double * A;
  double * x;
  double * y;
  SuiteSparse_long * indices;
  SuiteSparse_long * indptr;
  char * rbpath;

  if( argc < 2 ) {
      printf( "Usage: spmv_regular <RB file>\n" );
      exit(0);
  }

  rbpath = argv[1];

  if( argc == 3 ) {
    polybench_ntasks = atoi( argv[2] );
  }

  /* Initialize array(s). */
  init_array(rbpath,
              &nrow, &ncol,
              &indices, 
              &indptr,
	      &A,
	      &x,
              &y );

  // Float conversion

#ifdef POLYBENCH_GFLOPS
    polybench_program_total_flops = 2L * indptr[nrow] * REPS;
#endif
  polybench_start_instruments;

  /* Run kernel. */
  for( int t = 0; t < REPS; ++t )
  kernel_gesummv (nrow, 
                  indices, 
                  indptr,
		  A,
		  x,
		  y);

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(nrow, y));

#ifdef TEST_RESULTS
for( int i = 0; i < nrow; ++i ) {
    printf( "%d %lf\n", i, y[i] );
}
#endif
  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);

  return 0;
}
