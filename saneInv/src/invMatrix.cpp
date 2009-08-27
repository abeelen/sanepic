#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

extern "C" {
#include "nrutil.h"
#include "nrcode.h"
}

void reorderMatrix(long nbins, std::vector<string> listIn, double **MatrixIn,
		std::vector<string> listOut, double ***MatrixOut,double **mixmatOrig, int ncomp, double **&mixmat)

/*
 * Take a subsample of MatrixIn[ndet^2,nbins] by matching listIn & listOut
 * ouput MatrixOut[ndet,ndet*nbins]
 *  &
 *  Take a subsample of Mixing_matrix[ndet,ncomp] by matching listIn & listOut
 * ouput Mixmat[ndet,ncomp]
 */
{
	std::vector<int> indexIn; /*! Used to match input and output channels */
	long ndetIn = listIn.size(); /*! Input number of channels*/
	long ndetOut = listOut.size();/*! Output number of channels */

	// output number of detector cannot be larger than input number
	if(ndetOut>ndetIn){ // ajout mat 27/07
		cerr << "Input Noise Power Spectra must include all requested channels"	<< endl;
		exit(1);
	}

	// find indexes of input bolo file corresponding to the output bolo file
	indexIn.resize(listOut.size(), -1);
	for (int idetOut = 0; idetOut < ndetOut; idetOut++) {
		for (int idetIn = 0; idetIn < ndetIn; idetIn++) {
			if (listOut[idetOut] == listIn[idetIn]){
				indexIn[idetOut] = idetIn;
				break; } // ajout mat 27/07 TODO : Check que ca fait la meme chose
			}
	}

	// Check for missing data
	for (int idetOut = 0; idetOut < ndetOut; idetOut++) {
		if (indexIn[idetOut] == -1) {
			cerr << "Input Noise Power Spectra must include all requested channels"	<< endl;
			exit(1);
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

	// memory allocation
	mixmat=dmatrix(0,ndetOut - 1, 0, ncomp-1);

	// match input and output channels in the mixing matrix
	for(long ii=0; ii < ndetOut;ii++)
		for(int jj=0; jj < ncomp; jj++)
			mixmat[ii][jj]=mixmatOrig[indexIn[ii]][jj];
}




void inverseCovMatrixByMode(long nbins, long ndet, double **MatrixIn,
		double ***MatrixOut)

/*
 * Inverse the Covariance PowerSpectrum by mode
 */
{
	double **Mat_k, **iMat_k;
	double *p, *uvec, *ivec;

	Mat_k = dmatrix(0, ndet - 1, 0, ndet - 1);
	iMat_k = dmatrix(0, ndet - 1, 0, ndet - 1);
	*MatrixOut = dmatrix(0, ndet - 1, 0, ndet * nbins - 1);

	p = new double[ndet];
	uvec = new double[ndet];
	ivec = new double[ndet];

	for (int ibin = 0; ibin < nbins; ibin++) {

		cout << "Progress : " << ibin * 100. / nbins << " %\r" << flush;

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
		dcholdc(Mat_k, ndet, p);

		for (int idet1 = 0; idet1 < ndet; idet1++) {
			for (int idet2 = 0; idet2 < ndet; idet2++)
				uvec[idet2] = 0.0;

			uvec[idet1] = 1.0;
			dcholsl(Mat_k, ndet, p, uvec, ivec);

			for (int idet2 = 0; idet2 < ndet; idet2++)
				iMat_k[idet1][idet2] = ivec[idet2];
		}

		for (int idet1 = 0; idet1 < ndet; idet1++)
			for (int idet2 = 0; idet2 < ndet; idet2++)
				(*MatrixOut)[idet1][idet2 * nbins + ibin] = iMat_k[idet1][idet2];

	}

}


