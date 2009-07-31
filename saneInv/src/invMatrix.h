
#ifndef INVMATRIX_H_
#define INVMATRIX_H_

#include <cmath>

#include <string>
#include <vector>


void reorderMatrix(long nbins, std::vector<string> listIn, double **MatrixIn, std::vector<string> listOut, double ***MatrixOut,double **mixmatOrig, int ncomp, double **&mixmat);
void inverseCovMatrixByMode(long nbins, long ndet, double **MatrixIn, double ***MatrixOut);



#endif /* INVMATRIX_H_ */
