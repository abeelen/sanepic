
#include "covMatrixIO.h"
#include "invMatrix.h"
#include "parseInv.h"
#include "inputFileIO.h"
#include "mpi_architecture_builder.h"

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()
#include <algorithm>

extern "C" {
#include "nrutil.h"
#include "nrcode.h"
}

using namespace std;


int main(int argc, char *argv[]) {

	printf("Begin of saneInv\n\n");

	// data parameters
	/*!
	 * -ndet = number of detectors to output
	 * -ndetOrig = number of detectors in the NoiseNoise matrix
	 * -nbins = number of bins (Ell)
	 */
	long ndet, ndetOrig, nbins;
	int nbolos;
	int n_iter;

	double *ell; /*! bins values */

	struct directories dir;
	struct samples samples_struct;
	/*!
	 * -Rellth : Reduced NoiseNoise matrix
	 * -RellthOrig : Original NoiseNoise matrix
	 * -iRellth : Inverted reduced NoiseNoise matrix
	 * -mixmatOrig : original mixing matrix
	 * -mixmat : Reduced mixing matrix
	 */
	double **Rellth, **RellthOrig, **iRellth/*,**mixmatOrig,**mixmat*/;

	//	string noiseSp_dir_output;/*! output directory */
	string base_name="";/*! output noise file suffix */
	string fname; /*! covariance matrix fits filename */
	string fname2;
	//	char* fname_temp;
	string boloname;/*! channels list file */
	string extname = "_InvNoisePS";

	std::vector<string> channelIn; /*! Covariance matrix channel vector*/
	std::vector<string> channelOut; /*! bolometer reduction : Reduced vector of output channel */
	std::vector<int> indexIn; /*! bolometer index, used to determine which intput detector corresponds to which output detector*/

	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {
		int parsed=1;
		parsed=parse_saneInv_ini_file(argv[1], fname,samples_struct,dir, boloname, base_name);


		if (parsed==-1){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			exit(1);
		}
	}

	std::vector<string>::iterator it;
	long size_tmp;

//	cout << dir.tmp_dir + "bolonum_" + base_name + extname << endl;
//	cout << (int)samples_struct.noisevect.size() << endl;


	it = unique(samples_struct.noisevect.begin(), samples_struct.noisevect.end());
	size_tmp = it - samples_struct.noisevect.begin();
//	cout << size_tmp << endl;
	if(size_tmp==1){
		n_iter=1;
		cout << "The same covariance Matrix will be inverted for all the scans" << endl;
	}else{

		n_iter = (int)samples_struct.noisevect.size();
		if(n_iter==0){
			cerr << "Warning. You have forgot to mention covariance matrix in ini file or fits_filelist\n";
			exit(1);
		}
		cout << n_iter << " covariance Matrix will be inverted" << endl;
	}



	// Input argument for output : channellist
	read_strings(boloname, channelOut);

	//Total number of detectors to ouput (if ndet< ndetOrig : bolometer reduction)
	ndet = channelOut.size();
	printf("TOTAL NUMBER OF DETECTORS TO OUTPUT : %d\n", (int) ndet);



	for(int ii=0; ii<n_iter; ii++){
		fname="";
		fname+=(string)samples_struct.noisevect[ii];
		base_name=FitsBasename(fname);
//		cout << base_name << endl;

		fname2=(string)samples_struct.noisevect[ii];
		//		cout << fname2 << endl;
		//		cout << samples_struct.noisevect[ii] << endl;

		//		getchar();

		// read covariance matrix in a fits file named fname
		// returns : -the bins => Ell
		// -the input channel list => channelIn
		// -The number of bins (size of Ell) => nbins
		// -The original NoiseNoise covariance matrix => RellthOrig
		read_CovMatrix(fname2, channelIn, nbins, ell, RellthOrig);
		// total number of detectors in the covmatrix fits file
		ndetOrig = channelIn.size();
		printf("TOTAL NUMBER OF DETECTORS IN PS file: %d\n", (int) channelIn.size());
		//		getchar();

		//Deal with bolometer reduction and fill Rellth and mixmat
		reorderMatrix(nbins, channelIn, RellthOrig, channelOut, &Rellth);

		// Inverse reduced covariance Matrix : Returns iRellth
		inverseCovMatrixByMode(nbins, ndet, Rellth, &iRellth);

		//		cout << dir.tmp_dir + base_name + extname << endl;
		// write inversed noisePS in a binary file for each detector
		write_InvNoisePowerSpectra(channelOut, nbins, ell, iRellth, dir.tmp_dir, base_name + extname);

	}

	printf("\nEND OF SANEINV \n");

	nbolos = (int) channelIn.size();

	free_dmatrix(Rellth,0, ndet - 1, 0, ndet * nbins - 1);
	free_dmatrix(iRellth,0, ndet - 1, 0, ndet * nbins - 1);
	free_dmatrix(RellthOrig,0, nbolos * nbolos - 1, 0, nbins - 1);
	delete [] ell;

	return 0;
}

