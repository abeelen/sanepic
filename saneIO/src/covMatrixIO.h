

#ifndef COVMATRIX_H_
#define COVMATRIX_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

void read_bolofile(string fname, std::vector<string> &bolos);
void read_noisefile(string fname, string bolo1bolo2, double *ell, double *SPN, long *nbins);
long maxStringLength(std::vector<string> strings);
void write_InvNoisePowerSpectra(std::vector<string> bolos, long nbins, double * ell, double **Rellth, string ouputDir, string suffix);
void read_InvNoisePowerSpectra(string prefix, string boloName, string suffix, long * nbins, long * ndet, double ** ell, double *** SpN_all);
void write_CovMatrix(string fname, std::vector<string> &bolos, long nbins, double *ell, double **Rellth);
void read_CovMatrix(string fname, std::vector<string> &bolos, long *nbins, double **ell, double ***Rellth);
char* tableFormat(std::vector<string> strings);
char** vString2carray(std::vector<string> strings);

#endif /* COVMATRIX_H_ */
