#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>

#include "SaneInvTools.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>




using namespace std;

int reorderMatrix(long nbins, std::vector<string> listIn, gsl_matrix *MatrixIn, std::vector<string> listOut, gsl_matrix * & MatrixOut)
/* Resizes the covariance matrix with only needed detectors */
{
	std::vector<int> indexIn; /* Used to match input and output channels */
	long ndetIn = listIn.size(); /* Input number of channels*/
	long ndetOut = listOut.size();/* Output number of channels */

	// output number of detector cannot be larger than input number
	if(ndetOut>ndetIn){
		cerr << "Input Noise Power Spectra must include all requested channels"	<< endl;
		return 1;
	}

	// find indexes of input bolo file corresponding to the output bolo file
	indexIn.resize(listOut.size(), -1);
	for (int idetOut = 0; idetOut < ndetOut; idetOut++) {
		for (int idetIn = 0; idetIn < ndetIn; idetIn++) {
			if (listOut[idetOut] == listIn[idetIn]){
				indexIn[idetOut] = idetIn;
				break; }
		}
	}

	// Check for missing data
	if ( *std::min_element( indexIn.begin(), indexIn.end()) == -1 ) {
		cerr << "EE - Input Noise Power Spectra must include all requested channels (missing : ";
				for (int idetOut = 0; idetOut < ndetOut; idetOut++) {
					if (indexIn[idetOut] == -1) {
						cerr << listOut[idetOut] << ", ";
					}
				}
		cerr << " )" << endl;
		return 1;

	}

	MatrixOut = gsl_matrix_alloc(ndetOut, ndetOut*nbins);

	for (int idetOut1 = 0; idetOut1 < ndetOut; idetOut1++) {
		for (int idetOut2 = 0; idetOut2 < ndetOut; idetOut2++) {
			for (int ii = 0; ii < nbins; ii++) {
				// MatrixIn/Out should be symetric definite positive, but well... a mean does not hurt...
				gsl_matrix_set(MatrixOut, idetOut1, idetOut2 * nbins + ii, \
						( gsl_matrix_get(MatrixIn,indexIn[idetOut1] * ndetIn + indexIn[idetOut2],ii) \
								+ gsl_matrix_get(MatrixIn,indexIn[idetOut2] * ndetIn + indexIn[idetOut1],ii)) / 2);
			}
		}
	}

	return 0;
}




void inverseCovMatrixByMode(long nbins, long ndet, gsl_matrix * MatrixIn, gsl_matrix *& MatrixOut)
/* Inverses the Covariance PowerSpectrum by mode */
{
	gsl_matrix * Mat_k, * iMat_k;
	gsl_vector * uvec, * ivec;


	Mat_k     = gsl_matrix_alloc(ndet, ndet);
	iMat_k    = gsl_matrix_alloc(ndet, ndet);
	MatrixOut = gsl_matrix_alloc(ndet, ndet*nbins);

	uvec = gsl_vector_alloc(ndet);
	ivec = gsl_vector_alloc(ndet);

	//TODO: Parallelize this loop... -> Means we need to MPI_Reduce gsl_matrix
	for (int ibin = 0; ibin < nbins; ibin++) {

		//		cout << "Progress : " << ibin * 100. / nbins << "% \r" << flush;

		gsl_matrix_set_zero(iMat_k);
		// Matrix preparation
		for (int idet1 = 0; idet1 < ndet; idet1++)
			for (int idet2 = 0; idet2 < ndet; idet2++)
				gsl_matrix_set(Mat_k, idet1, idet2, gsl_matrix_get(MatrixIn,idet1, idet2*nbins+ibin));



		// Check for definite positive
		for (int idet1 = 0; idet1 < ndet; idet1++)
			for (int idet2 = 0; idet2 < ndet; idet2++)
				if (gsl_matrix_get(Mat_k,idet1,idet2) > sqrt(gsl_matrix_get(Mat_k,idet1,idet1) * gsl_matrix_get(Mat_k,idet2,idet2)))
					cerr << "EE - Matrix not definite positive, can not inverse... " << gsl_matrix_get(Mat_k,idet1,idet2) << endl;

		gsl_linalg_cholesky_decomp (Mat_k);

		for (int idet1 = 0; idet1 < ndet; idet1++) {

			gsl_vector_set_basis(uvec,idet1);
			gsl_linalg_cholesky_solve(Mat_k, uvec, ivec);

			for (int idet2 = 0; idet2 < ndet; idet2++)
				gsl_matrix_set(iMat_k,idet1,idet2, gsl_vector_get(ivec,idet2));
		}

		for (int idet1 = 0; idet1 < ndet; idet1++)
			for (int idet2 = 0; idet2 < ndet; idet2++)
				gsl_matrix_set(MatrixOut,idet1,idet2 * nbins + ibin, gsl_matrix_get(iMat_k,idet1,idet2));

	}

	// just to get a 100% value printed on screen
	//	cout << "Progress : 100.00% \r" << flush;

	gsl_matrix_free(Mat_k);
	gsl_matrix_free(iMat_k);
	gsl_vector_free(ivec);
	gsl_vector_free(uvec);
}
