
#ifndef INVMATRIX_H_
#define INVMATRIX_H_


#include <string>
#include <vector>

void reorderMatrix(long nbins, std::vector<std::string> listIn, double **MatrixIn, std::vector<std::string> listOut, double ***MatrixOut);
void inverseCovMatrixByMode(long nbins, long ndet, double **MatrixIn, double ***MatrixOut);
int who_do_it(int size, int rank, int ii);


#endif /* INVMATRIX_H_ */
