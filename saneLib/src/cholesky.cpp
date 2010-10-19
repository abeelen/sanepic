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
	// solution de Ax=b avec A triangulaire de taille n*m
	// on suppose que le systeme possede une unique solution
	// sup_inf booléen: vaut 0 si la matrice est triangulaire supérieure, 1 si elle est triangulaire inférieure
	int i, j;



// cas d'une matrice triangulaire superieure

if(sup_inf==0){

	i=n-1;
	j=m-1;


	//on recherche le premier element non nul en partant de l'element en bas a droite de la matrice et on remonte les lignes
	while (A[i][j]==0 && i>=0){
		i--;
	}


	for(int k=i; k>=0; k--){
		//on  suppose que le systeme a une unique solution donc si a_kj est nul, on n'a pas besoin de modifier le resultat
		if(A[k][j]!=0){
			// on calcul les valeurs de chaque composante du resultat en partant du bas de la matrice
			double element= b[j];
			for(int l=j+1; l<m; l++){
				element=element-x[l]*A[k][l];
			}
			element=element/A[k][j];
			x[j]=element;
			j--;
		}else{
			cerr << "attention a[i][j]=0 dans system_triang\n";
		}
	}
}
//on fait de meme si la matrice est triangulaire inferieur, mais on parcourt la matrice dans l'autre sens:
else{

	i=0;
	j=0;

// on cherche de meme le premier coefficient non nul en partant en haut à gauche de la matrice
	while (A[i][j]==0 && i<n){
		i++;
	}

// on calcul chaque coefficient en parcourant la matrice de haut en bas
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
					cerr << "attention a[i][j]=0 dans system_triang\n";
		}
	}
}

}




void cholesky(long n, double **a, double **l){

	double element;

	for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
			l[i][j]=0;

	//on suppose A symetrique, definie positive

// en parcourant la matrice dans le bon sens, on deduit chaque coefficient de L d'apres ceux calcules precedement

	l[0][0]=sqrt(a[0][0]);



	for(int i=1; i<n; i++){

//on rempli la matrice de gauche a droite, jusqu'a la diagonale
		for(int j=0; j<i; j++){
			element = a[i][j];
			for(int k=0; k<j; k++){
				element=element-(l[i][k]*l[j][k]);
			}
			l[i][j]=element/l[j][j];
		}
//on calcul le coefficient sur la diagonale
		element=a[i][i];
		for(int k=0; k<i; k++){
			element=element-l[i][k]*l[i][k];
		}
		l[i][i]=sqrt(element);
	}
}





void solve_cholesky(double **a, double *b, double **l, double *x, long n){

//on resoud ax=b, a matrice n*n symetrique, definie positive
// l est la matrice calculee par cholesky

	double **l_transp, *y;

// declaration memoire
	y=new double[n];
	l_transp=new double*[n];
	for(int i=0; i<n; i++)
		l_transp[i]=new double[n];

// on resoud ly=b
	system_triang(l,n,n,y,b,1);

// on transpose l
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++)
			l_transp[i][j]=l[j][i];
	}

//on resoud l_transp*x=y
	system_triang(l_transp,n,n,x,y,0);

//on libere la memoire
	delete [] y;
    for (int i = n; i > 0; --i)
      delete[] l_transp[i-1];
    delete[] l_transp;


}

