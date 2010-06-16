
#include "covMatrix_IO.h"
#include "invMatrix.h"
#include "inline_IO2.h"
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
}

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;


int main(int argc, char *argv[]) {

	int size; // MPI number of procs
	int rank; // My proc number
#ifdef USE_MPI

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if(rank==0)
		printf("Begin of saneInv\n\n");
#else
	size = 1;
	rank = 0;
	printf("Begin of saneInv\n\n");
#endif


	// data parameters
	/*!
	 * -ndet = number of detectors to output
	 * -ndetOrig = number of detectors in the NoiseNoise matrix
	 * -nbins = number of bins (Ell)
	 */
	long ndet, ndetOrig, nbins;
	int nbolos;
	long n_iter;

	double *ell; /*! bins values */

	struct common dir;
	struct samples samples_struct;
	/*!
	 * -Rellth : Reduced NoiseNoise matrix
	 * -RellthOrig : Original NoiseNoise matrix
	 * -iRellth : Inverted reduced NoiseNoise matrix
	 * -mixmatOrig : original mixing matrix
	 * -mixmat : Reduced mixing matrix
	 */
	double **Rellth, **RellthOrig, **iRellth;

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

	int parsed=1;

	if (argc<2) { // wrong number of arguments
		if(rank==0)
			printf("Please run %s using a *.ini file\n",argv[0]);
		parsed = -2;
	} else {
		parsed=parse_saneInv_ini_file(argv[1],samples_struct,dir, boloname, rank);
	}


	if (parsed<0){
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		exit(1);
	}

	std::vector<string>::iterator it;
	long size_tmp;

	it = unique(samples_struct.noisevect.begin(), samples_struct.noisevect.end());
	size_tmp = it - samples_struct.noisevect.begin();

	if(size_tmp==1){
		n_iter=1;
		if(rank==0)
			cout << "The same covariance Matrix will be inverted for all the scans\n" << endl;
	}else{
		n_iter = size_tmp;
//		n_iter = (int)samples_struct.noisevect.size();
		if(n_iter==0){
			if(rank==0)
				cerr << "WARNING. You have forgotten to mention covariance matrix in ini file or fits_filelist\n";
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			exit(1);
		}
		if(rank==0)
			cout << n_iter << " covariance Matrix will be inverted\n" << endl;
	}

	// Input argument for output : channellist
	read_strings(boloname, channelOut);

	//Total number of detectors to ouput (if ndet< ndetOrig : bolometer reduction)
	ndet = channelOut.size();

	if(rank==0)
		printf("TOTAL NUMBER OF DETECTORS TO OUTPUT : %d\n", (int) ndet);

	for(int ii=0; ii<n_iter; ii++){
		// select which proc number compute this loop
		if(rank==who_do_it(size,rank,ii)){
			fname="";
			fname+=(string)samples_struct.noisevect[ii];
			base_name=FitsBasename(fname);
			cout << base_name << endl;
			base_name=FitsBasename(fname);

			// get input covariance matrix file name
			fname2=dir.noise_dir + (string)samples_struct.noisevect[ii];

			// read covariance matrix in a fits file named fname
			// returns : -the bins => Ell
			// -the input channel list => channelIn
			// -The number of bins (size of Ell) => nbins
			// -The original NoiseNoise covariance matrix => RellthOrig
			read_CovMatrix(fname2, channelIn, nbins, ell, RellthOrig);

			// total number of detectors in the covmatrix fits file
			ndetOrig = channelIn.size();

			printf("TOTAL NUMBER OF DETECTORS IN PS file: %d\n", (int) channelIn.size());

			//Deal with bolometer reduction and fill Rellth and mixmat
			reorderMatrix(nbins, channelIn, RellthOrig, channelOut, &Rellth);

			// Inverse reduced covariance Matrix : Returns iRellth
			inverseCovMatrixByMode(nbins, ndet, Rellth, &iRellth);

			// write inversed noisePS in a binary file for each detector
			write_InvNoisePowerSpectra(channelOut, nbins, ell, iRellth, dir.tmp_dir, base_name + extname);

			// MAJ format file
			compute_dirfile_format_noisePS(dir.tmp_dir, channelOut, base_name + extname);

			// number of detector in the input channel list
			nbolos = (int) channelIn.size();

			// clean up
			free_dmatrix(Rellth,0, ndet - 1, 0, ndet * nbins - 1);
			free_dmatrix(iRellth,0, ndet - 1, 0, ndet * nbins - 1);
			free_dmatrix(RellthOrig,0, nbolos * nbolos - 1, 0, nbins - 1);
			delete [] ell;
		}
	}

	// clean up
	delete [] samples_struct.nsamples;

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

	printf("\nEND OF SANEINV \n");
	return 0;
}

