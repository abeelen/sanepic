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

using namespace std;

extern "C" {
#include "nrutil.h"
#include "nrcode.h"
#include <fitsio.h>
}



int main(int argc, char *argv[]) {

	// data parameters
	long ndet, ndetOrig, nbins;

	double *ell1;
	double **Rellth, **RellthOrig, **iRellth;

	string noiseSp_dir_output;
	string extentnoiseSp;
	string fname;
	string boloname;

	std::vector<string> channelIn;
	std::vector<string> channelOut;
	std::vector<int> indexIn;

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



	read_CovMatrix(fname/*argv[1]*/, channelIn, &nbins, &ell1, &RellthOrig);

	ndetOrig = channelIn.size();
	printf("TOTAL NUMBER OF DETECTORS IN PS file: %d\n", (int) channelIn.size());

	// Input argument for output : channellist
	read_bolofile(boloname/*argv[2]*/, channelOut);
	ndet = channelOut.size();
	printf("TOTAL NUMBER OF DETECTORS TO OUTPUT : %d\n", (int) ndet);


	//noiseSp_dir_output = string(argv[3]);
	//extentnoiseSp = string(argv[4]);

	reorderMatrix(nbins, channelIn, RellthOrig, channelOut, &Rellth);
	inverseCovMatrixByMode(nbins, ndet, Rellth, &iRellth);

	write_InvNoisePowerSpectra(channelOut, nbins, ell1, iRellth, noiseSp_dir_output, extentnoiseSp);

	return 0;
}

