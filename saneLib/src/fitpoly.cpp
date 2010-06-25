/*
 * fitpoly.cpp
 *
 *  Created on: 31 mai 2010
 *      Author: stage
 */

#include "fitpoly.h"
#include <iostream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

using namespace std;



void fitpoly(int norder, long taille, double *sx, double *sy, double *a){

	//norder: ordre demande, taille: taille des donnees, sx et sy: donnees, a: coefficients du polynome resultat
	//taille doit etre plus grand que norder
	//sx doit etre normalise entre -1 et 1

	gsl_matrix *V;
	gsl_vector *Y, *res, *tau, *residual;


	if(taille<norder+1){
		cout << "Warning in fitpoly : Polynomial is not unique; degree >= number of data points.\n";
		cout << "data points : " << taille << endl;
		cout << "degree : " << norder +1 << endl;
		cout << "Returning degree 0 polynomia instead ...\n";
		double mean=0.0;
		for(long jj=0; jj<taille; jj++)
			mean += sy[jj];
		a[0]=mean/(double)taille;
		for(long jj=1; jj<norder+1; jj++)
			a[jj]=0.0;
		return;
	}

	//declaration de la taille des matrices et vecteurs utilisés
	V=gsl_matrix_calloc ((size_t) taille, (size_t) norder+1);
	Y = gsl_vector_calloc ((size_t) taille);
	residual=gsl_vector_calloc((size_t) taille);
	res = gsl_vector_calloc ((size_t) norder+1);
	tau = gsl_vector_calloc ((size_t) norder+1);


	// conversion de sy en gsl_vector
	for(int i=0; i<taille; i++){
		gsl_vector_set(Y, i, sy[i]);
	}



	// ecriture de la matrice de vandermonde

	for(int i=0; i<taille; i++){
		double element=1.0;
		gsl_matrix_set(V,i,norder,element);
		for(int j=norder-1; j>= 0; j--){
			element=sx[i]*element;
			gsl_matrix_set(V,i,j,element);
		}
	}

	//	gsl_matrix_fprintf(stdout,V, " %lf");

	//	cout << "avant decomp \n";
	//	gsl_vector_fprintf(stdout,tau, " %lf");

	// decomposition QR et résolution par methode des moindres carres
	gsl_linalg_QR_decomp (V, tau);

	//	gsl_vector_fprintf(stdout,tau, " %lf");

	//	cout << "avant solve \n";

	gsl_linalg_QR_lssolve (V, tau, Y, res,residual);

	//conversion de la matrice res en tableau a
	for(int i=0; i<norder+1; i++)
		a[i]=gsl_vector_get(res,norder-i);

	//on libere la memoire
	gsl_matrix_free (V);
	gsl_vector_free (Y);
	gsl_vector_free(residual);
	gsl_vector_free (res);
	gsl_vector_free (tau);

}
