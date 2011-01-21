
#ifndef COVMATRIX_IO_H_
#define COVMATRIX_IO_H_

#include <string>
#include <vector>

extern "C" {
#include "getdata.h"
}

/*!  This function writes the Inverse Covariance Matrices in binary format */
int write_InvNoisePowerSpectra(DIRFILE* D, std::vector<std::string> bolos, long nbins, double * ell, double **Rellth, std::string suffix);

/*! This function reads the Inverse Covariance Matrices in binary format */
int read_InvNoisePowerSpectra(DIRFILE* D, std::string outputDir, std::string boloName, std::string suffix, long nbins, long ndet, double ** ell, double *** SpN_all);

/*! This function writes the NoiseNoise Matrices in a fits file (also Writes the mixing matrices) */
int write_CovMatrix(std::string fname, std::vector<std::string> bolos, long nbins, double *ell, double **Rellth);

/*! Writes the reduced mixing matrix from a binary file */
//void read_ReducedMixingMatrix(double **&mixmat,long &ndet,int &ncomp, std::string dir);

/*! This function reads the NoiseNoise Matrices in a fits file (also reads the mixing matrices) */
int read_CovMatrix(std::string fname, std::vector<std::string> &bolos, long &nbins, double *&ell, double **&Rellth);

/*! Returns the table Format for the given vector of string*/
char* tableFormat(std::vector<std::string> strings);

/*! Return the longest string length of a string vector */
long maxStringLength(std::vector<std::string> strings);

/*! Transform a vector of string into a array of char*/
char** vString2carray(std::vector<std::string> strings);

int write_psd_tofits(std::string fname, long nx, long ny,char dtype, void * psd1d);


#endif /* COVMATRIX_IO_H_ */
