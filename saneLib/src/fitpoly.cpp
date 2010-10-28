#include "fitpoly.h"
#include <iostream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

using namespace std;



void fitpoly(int norder, long taille, double *sx, double *sy, double *a){

	// norder: output poly order, taille: data size, sx et sy: data, a: output poly coefficients
	// WARNING  : taille must be > norder
	// sx is normalised between -1 and 1

	gsl_matrix *V;
	gsl_vector *Y, *resultat, *tau, *residual;

	// malloc
	V=gsl_matrix_calloc ((size_t) taille, (size_t) norder+1);
	Y = gsl_vector_calloc ((size_t) taille);
	residual=gsl_vector_calloc((size_t) taille);
	resultat = gsl_vector_calloc ((size_t) norder+1);
	tau = gsl_vector_calloc ((size_t) norder+1);


	// sy is converted to a gsl vector
	for(int i=0; i<taille; i++){
		gsl_vector_set(Y, i, sy[i]);
	}



	// vandermonde matrix computation

	for(int i=0; i<taille; i++){
		double element=1.0;
		gsl_matrix_set(V,i,norder,element);
		for(int j=norder-1; j>= 0; j--){
			element=sx[i]*element;
			gsl_matrix_set(V,i,j,element);
		}
	}

	// QR-decomposition and least-squares computation
	gsl_linalg_QR_decomp (V, tau);
	gsl_linalg_QR_lssolve (V, tau, Y, resultat,residual);

	// converting a gls vector to the expected result
	for(int i=0; i<norder+1; i++)
		a[i]=gsl_vector_get(resultat,norder-i);

	// memory cleanup
	gsl_matrix_free (V);
	gsl_vector_free (Y);
	gsl_vector_free(residual);
	gsl_vector_free (resultat);
	gsl_vector_free (tau);

}
