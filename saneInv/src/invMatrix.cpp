#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <cmath>

#include <string>
#include <vector>
#include <list>

#include "covMatrixIO.h"

using namespace std;

extern "C" {
#include "nrutil.h"
#include "nrcode.h"
#include <fitsio.h>
}

//**********************************************************************************//
//**********************************************************************************//
//*************************** Beginning of main program ****************************//
//**********************************************************************************//
//**********************************************************************************//


int main(int argc, char *argv[]) {

	long ii;
	long idet1, idet2, ibin;

	// data parameters
	long ndet, ndetOrig, nbins;

	double *ell1, *p, *uvec, *ivec, *temparray;
	double **Mat_k, **iMat_k, **Rellth, **RellthOrig, **iRellth;

	string field;
	string field1;
	string field2;
	string bolonamesIn;
	string bolofield;
	string bolofield1;
	string bolofield2;
	string file_offsets;
	string termin;
	string noiseSpfile;
	string noiseSp_dir_output;
	string extentnoiseSp;

	std::vector<string> channelIn;
	std::vector<string> channelOut;
	std::vector<int> indexIn;

	// Read input spectra

	//  // Input argument c: hannelList ....
	//  read_bolofile(argv[1], channelIn);
	//  //  cerr << "num ch: "<< channel.size() << endl;
	//  if (channelIn.size() == 0) {
	//    cerr << "Must provide a valid input channel list.\n\n";
	//  }
	//  printf("TOTAL NUMBER OF DETECTORS IN PS file: %d\n",(int)channelIn.size());
	//
	//  // ... Corresponding NoisePS file
	//  noiseSpfile = argv[2];
	//
	//  // Read input PS file
	//  ifstream Spfile(argv[2], ios::binary);
	//  if (!Spfile.is_open()){
	//    cerr << "Could not open PS file" << endl;
	//    exit(1);
	//  }
	//
	//  Spfile.read(reinterpret_cast<char *>(&ndetOrig),sizeof(long));
	//  Spfile.read(reinterpret_cast<char *>(&nbins),sizeof(long));
	//
	//  nbins = (int) nbins;
	//
	//  if (ndetOrig != channelIn.size()) {
	//    cerr << "Input Channel List does not correspond to PS file\n";
	//  }
	//
	//  //  printf("%d bins\n",(int)nbins);
	//
	//  ell1 = new double[nbins+1];
	//  RellthOrig  = dmatrix(0,ndetOrig*ndetOrig-1,0,nbins-1);
	//
	//  for (ii=0; ii<nbins+1; ii++)
	//    Spfile.read(reinterpret_cast<char *>(&ell1[ii]), sizeof(double));
	//
	//  for (idet1 = 0;   idet1 < ndetOrig; idet1++)
	//    for (idet2 = 0; idet2 < ndetOrig; idet2++)
	//      for ( ii=0; ii<nbins; ii++)
	//	Spfile.read(reinterpret_cast<char *>(&RellthOrig[idet1*ndetOrig+idet2][ii]),sizeof(double));
	//
	//  Spfile.close();
	//
	//
	//  write_CovMatrix("!essai.fits", channelIn, nbins, ell1, RellthOrig);
	read_CovMatrix(argv[1], channelIn, &nbins, &ell1, &RellthOrig);

	ndetOrig = channelIn.size();
	printf("TOTAL NUMBER OF DETECTORS IN PS file: %d\n", (int) channelIn.size());

	// Input argument for output : channellist
	read_bolofile(argv[2], channelOut);
	if (channelOut.size() == 0) {
		cerr << "Must provide a valid output channel list.\n\n";
	}
	ndet = channelOut.size();
	printf("TOTAL NUMBER OF DETECTORS TO OUTPUT : %d\n", (int) ndet);

	/*
	 * TODO: Check Channel match
	 */

	noiseSp_dir_output = argv[3];
	extentnoiseSp = argv[4];


	//************************************************************************//
	//************************************************************************//
	//program starts here
	//************************************************************************//
	//************************************************************************//


	// Take subsample of the RellthOrig
	Rellth = dmatrix(0, ndet - 1, 0, ndet * nbins - 1);

	// find indexes of input bolo file corresponding to the output bolo file
	indexIn.resize(channelOut.size(), -1);
	for (idet1 = 0; idet1 < ndet; idet1++) {
		for (idet2 = 0; idet2 < ndetOrig; idet2++) {
			if (channelOut[idet1] == channelIn[idet2])
				indexIn[idet1] = idet2;
		}
	}
	//   for (idet1=0; idet1<ndet; idet1++)
	//     cout << idet1 << " " << channelOut[idet1] << " " << indexIn[idet1] << endl;

	for (idet1 = 0; idet1 < ndet; idet1++) {
		for (idet2 = 0; idet2 < ndet; idet2++) {
			for (ii = 0; ii < nbins; ii++) {
				Rellth[idet1][idet2 * nbins + ii]
						= (RellthOrig[indexIn[idet1] * ndetOrig
								+ indexIn[idet2]][ii]
								+ RellthOrig[indexIn[idet2] * ndetOrig
										+ indexIn[idet1]][ii]) / 2;
			}
		}
	}

	Mat_k = dmatrix(0, ndet - 1, 0, ndet - 1);
	iMat_k = dmatrix(0, ndet - 1, 0, ndet - 1);
	iRellth = dmatrix(0, ndet - 1, 0, ndet * nbins - 1);
	p = new double[ndet];
	uvec = new double[ndet];
	ivec = new double[ndet];
	temparray = new double[(ndet + 1) * nbins + 2];
	//  SPN = new double[nbins];


	for (ibin = 0; ibin < nbins; ibin++) {

		cout << "Progress : " << ibin * 100. / nbins << " %\r" << flush;

		for (idet1 = 0; idet1 < ndet; idet1++) {
			for (idet2 = 0; idet2 < ndet; idet2++) {
				Mat_k[idet1][idet2] = Rellth[idet1][idet2 * nbins + ibin];
				iMat_k[idet1][idet2] = 0.0;
			}
		}

		for (idet1 = 0; idet1 < ndet; idet1++) {
			for (idet2 = 0; idet2 < ndet; idet2++) {
				if (Mat_k[idet1][idet2] > sqrt(Mat_k[idet1][idet1]
						* Mat_k[idet2][idet2]))
					printf("ERROR %10.15g\t", Mat_k[idet1][idet2]);
			}
		}

		// invert the covariance matrix per mode
		dcholdc(Mat_k, ndet, p);
		for (idet1 = 0; idet1 < ndet; idet1++) {
			for (idet2 = 0; idet2 < ndet; idet2++)
				uvec[idet2] = 0.0;
			uvec[idet1] = 1.0;
			dcholsl(Mat_k, ndet, p, uvec, ivec);
			for (idet2 = 0; idet2 < ndet; idet2++)
				iMat_k[idet1][idet2] = ivec[idet2];
		}

		for (idet1 = 0; idet1 < ndet; idet1++)
			for (idet2 = 0; idet2 < ndet; idet2++)
				iRellth[idet1][idet2 * nbins + ibin] = iMat_k[idet1][idet2];

	}

	/// write inverse covariance matrix to disk

	write_InvNoisePowerSpectra(channelOut, nbins, ell1, iRellth, noiseSp_dir_output, extentnoiseSp);


	//******************************************************************//
	//******************************************************************//
	//************************  End of loop ****************************//
	//******************************************************************//
	//******************************************************************//


	return 0;
}

