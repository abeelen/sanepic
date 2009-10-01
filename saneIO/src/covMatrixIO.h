

#ifndef COVMATRIX_H_
#define COVMATRIX_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;


//void read_noisefile(string fname, string bolo1bolo2, double *ell, double *SPN, long *nbins);

/*! Return the longest string length of a string vector */
long maxStringLength(std::vector<string> strings);


/*!  This function writes the Inverse Covariance Matrices in binary format */
void write_InvNoisePowerSpectra(std::vector<string> bolos, long nbins, double * ell, double **Rellth, string ouputDir, string suffix);

/*! This function reads the Inverse Covariance Matrices in binary format */
void read_InvNoisePowerSpectra(string prefix, string boloName, string suffix, long * nbins, long * ndet, double ** ell, double *** SpN_all);

/*! This function writes the NoiseNoise Matrices in a fits file (also Writes the mixing matrices) */
void write_CovMatrix(string fname, std::vector<string> bolos, long nbins, double *ell, double **Rellth);
void write_CovMatrix2(string fname, std::vector<string> bolos, long nbins, double *ell, double **Rellth);

/*! Writes the reduced mixing matrix in a binary file */
void write_ReducedMixingMatrix(double **mixmat,long ndet,int ncomp, string outputDir);

/*! Writes the reduced mixing matrix from a binary file */
void read_ReducedMixingMatrix(double **&mixmat,long &ndet,int &ncomp, string dir);

/*! This function reads the NoiseNoise Matrices in a fits file (also reads the mixing matrices) */
void read_CovMatrix(string fname, std::vector<string> &bolos, long &nbins, double *&ell, double **&Rellth);

/*! Returns the table Format for the given vector of string*/
char* tableFormat(std::vector<string> strings);

/*! Transform a vector of string into a array of char*/
char** vString2carray(std::vector<string> strings);

#endif /* COVMATRIX_H_ */
