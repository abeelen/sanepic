#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <time.h>
#include <cmath>
#include <vector>
#include <string>

#include "cholesky.h"
#include "imageIO.h"
#include "temporary_IO.h"
#include "mpi_architecture_builder.h"
#include "todprocess.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>


extern "C" {
#include "nrutil.h"
}


using namespace std;




void system_triang(double **A, long n, long m, double *x, double *b, bool sup_inf){
	// we solve A.x=b with A : n*m triangular matrix
	// we supose that the system has a single solution
	// sup_inf bool: 0 if the matrix A is upper triangular, 1 for lower triangular
	int i, j;



	// upper case

	if(sup_inf==0){

		i=n-1;
		j=m-1;


		// find first non-null element starting with right-bottom element and raising the lines
		while (A[i][j]==0 && i>=0){
			i--;
		}


		for(int k=i; k>=0; k--){
			// solution is unique so if a_kj is null, we don't need to modify the results
			if(A[k][j]!=0){
				// each component is computed starting from the botom of the matrix
				double element= b[j];
				for(int l=j+1; l<m; l++){
					element=element-x[l]*A[k][l];
				}
				element=element/A[k][j];
				x[j]=element;
				j--;
			}else{
				cerr << "warning a[i][j]=0 in system_triang\n";
			}
		}
	}
	// the same steps are done if the matrix is lower-triang but we run the matrix from the top to the bottom
	else{

		i=0;
		j=0;

		// first non-null element starting from left-top element and running down the lines
		while (A[i][j]==0 && i<n){
			i++;
		}

		// each component is computing running the matrix from top to bottom
		for(int k=i; k<n; k++){
			if(A[k][j]!=0){
				double element= b[k];
				for(int l=0; l<j; l++){
					element=element-x[l]*A[k][l];
				}
				element=element/A[k][j];
				x[j]=element;
				j++;
			}else{
				cerr << "warning a[i][j]=0 in system_triang\n";
			}
		}
	}

}




void cholesky(long n, double **a, double **l){

	double element;

	for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
			l[i][j]=0;

	// a : matrix symetric, positive-definite

	l[0][0]=sqrt(a[0][0]);



	for(int i=1; i<n; i++){

		// the matrix is fill from left to right until diagonal
		for(int j=0; j<i; j++){
			element = a[i][j];
			for(int k=0; k<j; k++){
				element=element-(l[i][k]*l[j][k]);
			}
			l[i][j]=element/l[j][j];
		}

		// diagonal coefficient are computed
		element=a[i][i];
		for(int k=0; k<i; k++){
			element=element-l[i][k]*l[i][k];
		}
		l[i][i]=sqrt(element);
	}
}





void solve_cholesky(double **a, double *b, double **l, double *x, long n){

	//we resolve ax=b, a : matrix symetric n*n, positive-definite
	// l is the matrix computed by cholesky

	double **l_transp, *y;

	// malloc
	y=new double[n];
	l_transp=new double*[n];
	for(int i=0; i<n; i++)
		l_transp[i]=new double[n];

	// we resolve l.y=b
	system_triang(l,n,n,y,b,1);

	// l is transposed
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++)
			l_transp[i][j]=l[j][i];
	}

	// we resolve l_transp*x=y
	system_triang(l_transp,n,n,x,y,0);

	// memory cleanup
	delete [] y;
    for (int i = n; i > 0; --i)
		delete[] l_transp[i-1];
	delete[] l_transp;


}

