#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <cmath>

#include <string>
#include <vector>
#include <list>

#include "covMatrixIO.h"
#include "invMatrix.h"
#include "parseInv.h"
#include "inputFileIO.h"

using namespace std;

extern "C" {
#include "nrutil.h"
#include "nrcode.h"
#include <fitsio.h>
}



int main(int argc, char *argv[]) {

	// data parameters
	/*!
	 * -ndet = number of detectors to output
	 * -ndetOrig = number of detectors in the NoiseNoise matrix
	 * -nbins = number of bins (Ell)
	 */
	long ndet, ndetOrig, nbins;

	double *ell; /*! bins values */
	/*!
	 * -Rellth : Reduced NoiseNoise matrix
	 * -RellthOrig : Original NoiseNoise matrix
	 * -iRellth : Inverted reduced NoiseNoise matrix
	 * -mixmatOrig : original mixing matrix
	 * -mixmat : Reduced mixing matrix
	 */
	double **Rellth, **RellthOrig, **iRellth/*,**mixmatOrig,**mixmat*/;

	string noiseSp_dir_output;/*! output directory */
	string extentnoiseSp;/*! output noise file suffix */
	string fname; /*! covariance matrix fits filename */
	string boloname;/*! channels list file */

	std::vector<string> channelIn; /*! Covariance matrix channel vector*/
	std::vector<string> channelOut; /*! bolometer reduction : Reduced vector of output channel */
	std::vector<int> indexIn; /*! bolometer index, used to determine which intput detector corresponds to which output detector*/

	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {
		int parsed=1;
		parsed=parse_saneInv_ini_file(argv[1], fname, boloname, noiseSp_dir_output, extentnoiseSp);


		if (parsed==-1){
#ifdef USE_MPI
			MPI_Finalize();
#endif
			exit(1);
		}
	}


	// read covariance matrix in a fits file named fname
	// returns : -the bins => Ell
	// -the input channel list => channelIn
	// -The number of bins (size of Ell) => nbins
	// -The original NoiseNoise covariance matrix => RellthOrig
	read_CovMatrix(fname, channelIn, nbins, ell, RellthOrig);
	// total number of detectors in the covmatrix fits file
	ndetOrig = channelIn.size();
	printf("TOTAL NUMBER OF DETECTORS IN PS file: %d\n", (int) channelIn.size());

	// Input argument for output : channellist
	read_strings(boloname, channelOut);

	//Total number of detectors to ouput (if ndet< ndetOrig : bolometer reduction)
	ndet = channelOut.size();
	printf("TOTAL NUMBER OF DETECTORS TO OUTPUT : %d\n", (int) ndet);

	//Deal with bolometer reduction and fill Rellth and mixmat
	reorderMatrix(nbins, channelIn, RellthOrig, channelOut, &Rellth);

	// Inverse reduced covariance Matrix : Returns iRellth
	inverseCovMatrixByMode(nbins, ndet, Rellth, &iRellth);

	// write inversed noisePS in a binary file for each detector
	write_InvNoisePowerSpectra(channelOut, nbins, ell, iRellth, noiseSp_dir_output, extentnoiseSp);

	// write Reduced mixing matrix in a binary file
	//write_ReducedMixingMatrix(mixmat,ndet,ncomp,noiseSp_dir_output);

	return 0;
}

