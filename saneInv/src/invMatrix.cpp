#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>
#include "invMatrix.h"
#include "cholesky.h"

extern "C" {
#include "nrutil.h"
}

using namespace std;

int reorderMatrix(long nbins, std::vector<string> listIn, double **MatrixIn,
		std::vector<string> listOut, double ***MatrixOut)
/*!\brief Resizes the covariance matrix with only needed detectors
 * \param nbins Number of power spectrum bins
 * \param listIn Input list of detectors
 * \param MatrixIn Input Covariance Matrix
 * \param listOut Output list of desired detectors
 * \param MatrixOut Filled Covariance Matrix with only desired detectors informations
 */
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
	for (int idetOut = 0; idetOut < ndetOut; idetOut++) {
		if (indexIn[idetOut] == -1) {
			cerr << "Input Noise Power Spectra must include all requested channels (" << listOut[idetOut] << ")" << endl;
			return 1;
		}
	}

	// memory allocation
	*MatrixOut = dmatrix(0, ndetOut - 1, 0, ndetOut * nbins - 1);

	for (int idetOut1 = 0; idetOut1 < ndetOut; idetOut1++) {
		for (int idetOut2 = 0; idetOut2 < ndetOut; idetOut2++) {
			for (int ii = 0; ii < nbins; ii++) {
				// MatrixIn/Out should be symetric definite positive, but well... a mean does not hurt...
				(*MatrixOut)[idetOut1][idetOut2 * nbins + ii]
				                       = (MatrixIn[indexIn[idetOut1] * ndetIn
				                                   + indexIn[idetOut2]][ii]
				                                                        + MatrixIn[indexIn[idetOut2] * ndetIn
				                                                                   + indexIn[idetOut1]][ii]) / 2;
			}
		}
	}

	return 0;
}




void inverseCovMatrixByMode(long nbins, long ndet, double **MatrixIn,
		double ***MatrixOut)
/*!\brief Inverse the Covariance PowerSpectrum by mode
 * \param nbins Number of power spectrum bins
 * \param ndet Number of detectors
 * \param MatrixIn The Matrix that will be inverted
 * \param MatrixOut The Inverted Matrix (all modes)
 */
{
	double **Mat_k, **iMat_k;
	double *uvec, *ivec;

	double **l;


	l=new double*[ndet];
	for(long ii=0;ii<ndet;ii++)
		l[ii]=new double[ndet];


	Mat_k = dmatrix(0, ndet - 1, 0, ndet - 1); // k-Mode Matrix
	iMat_k = dmatrix(0, ndet - 1, 0, ndet - 1); // inverted k-Mode Matrix
	*MatrixOut = dmatrix(0, ndet - 1, 0, ndet * nbins - 1); // whole Mode inverted Matrix

	uvec = new double[ndet];
	ivec = new double[ndet];

	for (int ibin = 0; ibin < nbins; ibin++) {

//		cout << "Progress : " << ibin * 100. / nbins << "% \r" << flush;

		// Matrix preparation
		for (int idet1 = 0; idet1 < ndet; idet1++) {
			for (int idet2 = 0; idet2 < ndet; idet2++) {

				Mat_k[idet1][idet2] = MatrixIn[idet1][idet2 * nbins + ibin];
				iMat_k[idet1][idet2] = 0.0;
			}
		}

		// Check for definit positive
		for (int idet1 = 0; idet1 < ndet; idet1++) {
			for (int idet2 = 0; idet2 < ndet; idet2++) {
				if (Mat_k[idet1][idet2] > sqrt(Mat_k[idet1][idet1]
				                                            * Mat_k[idet2][idet2]))
					cerr << "ERROR " << Mat_k[idet1][idet2] << endl;
			}
		}

		// invert the covariance matrix per mode
		cholesky(ndet, Mat_k, l);

		for (int idet1 = 0; idet1 < ndet; idet1++) {
			for (int idet2 = 0; idet2 < ndet; idet2++)
				uvec[idet2] = 0.0;

			uvec[idet1] = 1.0;
			solve_cholesky(Mat_k, uvec, l, ivec, ndet);

			for (int idet2 = 0; idet2 < ndet; idet2++)
				iMat_k[idet1][idet2] = ivec[idet2];
		}

		for (int idet1 = 0; idet1 < ndet; idet1++)
			for (int idet2 = 0; idet2 < ndet; idet2++)
				(*MatrixOut)[idet1][idet2 * nbins + ibin] = iMat_k[idet1][idet2];

	}

	// just to get a 100% value printed on screen
//	cout << "Progress : 100.00% \r" << flush;

	// clean up
	free_dmatrix(Mat_k,0, ndet - 1, 0, ndet - 1);
	free_dmatrix(iMat_k,0, ndet - 1, 0, ndet - 1);
	delete [] ivec;
	delete [] uvec;
	for(long ii=0;ii<ndet;ii++)
		delete [] l[ii];
	delete [] l;
}

