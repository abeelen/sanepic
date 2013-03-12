#ifndef INVMATRIX_H_
#define INVMATRIX_H_


#include <string>
#include <vector>
#include <gsl/gsl_matrix.h>

//! Resizes the covariance matrix with only needed detectors
/*!
 * \param nbins Number of power spectrum bins
 * \param listIn Input list of detectors
 * \param MatrixIn Input Covariance Matrix
 * \param listOut Output list of desired detectors
 * \param MatrixOut Filled Covariance Matrix with only desired detectors informations
 * \return An integer specifying if there were an error (>0) or not (=0)
 */
int reorderMatrix(long nbins, std::vector<std::string> listIn, gsl_matrix * MatrixIn, std::vector<std::string> listOut, gsl_matrix *& MatrixOut);

//! Inverses the Covariance PowerSpectrum by mode
/*!
 * \param nbins Number of power spectrum bins
 * \param ndet Number of detectors
 * \param MatrixIn The Matrix that will be inverted
 * \param MatrixOut The Inverted Matrix (all modes)
 */
void inverseCovMatrixByMode(long nbins, long ndet, gsl_matrix *MatrixIn, gsl_matrix *& MatrixOut);


#endif /* INVMATRIX_H_ */
