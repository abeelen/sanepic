#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()
#include <algorithm>
#include <sysexits.h>

#include "covMatrix_IO.h"
#include "invMatrix.h"
#include "temporary_IO.h"
#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "parser_functions.h"
#include "struct_definition.h"
#include "error_code.h"

#include <gsl/gsl_matrix.h>




#ifdef PARA_FRAME
#include "mpi.h"
#endif

using namespace std;

/*!
 *  This is organized as :
 *
 *  - parse the input ini file and verify his validity
 *  - check for existence of directory/files pointed from the ini file
 *  - Print parser output to screen
 *
 *  - for each file :
 *      - Generate or clear the dirfile parts that will be filled : Noise and Noise/ell
 *
 *  - Read all channel files, store it into a vector<vector> (and commit to other ranks if needed)
 *
 *	For each Covariance Matrix to invert :
 *		- Read Covariance Matrix from disk (fits file)
 *		- Deal with bolometer reduction if needed (ReorderMatrix)
 *		- Invert Matrix
 *		- Write Inverted Power Spectra to disk in dirfiles
 */

int main(int argc, char *argv[]) {

	int size; // MPI number of procs
	int rank; // this proc number

#ifdef PARA_FRAME

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
	size = 1;
	rank = 0;
#endif

	if(rank==0)
	  cout << endl << "saneInv : inversion of the Noise-Noise PowerSpectra" << endl;

	// data parameters
	/*
	 * -ndet = number of detectors to output
	 * -ndetOrig = number of detectors in the NoiseNoise matrix
	 * -nbins = number of bins (Ell)
	 */
	long ndet, nbins;
	long n_iter;

	double *ell; /* bins values */

	struct param_common dir;
	struct samples samples_struct;

	/*
	 * -Rellth : Reduced NoiseNoise matrix
	 * -RellthOrig : Original NoiseNoise matrix
	 * -iRellth : Inverted reduced NoiseNoise matrix
	 * -mixmatOrig : original mixing matrix
	 * -mixmat : Reduced mixing matrix
	 */
	gsl_matrix *Rellth, *RellthOrig, *iRellth;

	//	string noiseSp_dir_output;/* output directory */
	string boloname;/* channels list file */
	string noise_suffix = "_InvNoisePS";
	string output = "";

	struct param_sanePos pos_param;
	struct param_saneProc proc_param;
	struct param_sanePS structPS;
	struct param_saneInv saneInv_struct;
	struct param_sanePic struct_sanePic;

	std::vector<std::vector<std::string> > bolo_list; // this vector contains all bolonames for all the scans

	std::vector<int> indexIn; /* bolometer index, used to determine which intput detector corresponds to which output detector*/

	uint16_t mask_saneInv = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | TMP_PATH_PROBLEM |
			BOLOFILE_NOT_FOUND | SANEINV_INPUT_ERROR | FITS_FILELIST_NOT_FOUND; // 0x601b

	uint16_t parsed=0x0000; // parser error status
	uint16_t compare_to_mask; // parser error status

	//	if (rank==0){ // root parse ini file and fill the structures. Also print warnings or errors

	// Parse ini file
	if (argc<2) {
		compare_to_mask=0x001;
	} else {
		parsed=parser_function(argv[1], output, dir, samples_struct, pos_param, proc_param,
				structPS, saneInv_struct, struct_sanePic, size, rank);

		compare_to_mask = parsed & mask_saneInv;

		// print parser warning and/or errors
		if (rank == 0)
			cout << endl << output << endl;
	}

	if(compare_to_mask>0x0000){

		switch (compare_to_mask){/* error during parsing phase */

		case 0x0001: printf("Please run %s using a correct *.ini file\n",argv[0]);
		break;

		default : printf("Wrong program options or argument. Exiting !\n");
		break;


		}

#ifdef PARA_FRAME
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return EX_CONFIG;
	}
	//	}


#ifdef PARA_FRAME

	MPI_Barrier(MPI_COMM_WORLD);

	if(configure_PARA_FRAME_samples_struct(dir.tmp_dir, samples_struct, rank, size)){
		MPI_Abort(MPI_COMM_WORLD, 1);
		return EX_IOERR;
	}

	MPI_Barrier(MPI_COMM_WORLD);

#endif


	/* ------------------------------------- READ bolo list ----------------------------*/

	if(channel_list_to_vect_list(samples_struct, bolo_list, rank)){
		cout << "error in channel_list_to_vect_list" << endl;
		return EX_CONFIG;
	}

	/* ------------------------------------------------------------------------------------*/

	if(rank==0)
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				structPS, struct_sanePic, saneInv_struct);


	// START OF saneInv

	//	std::vector<string>::iterator it;
	//	long size_tmp;
	//
	//	it = unique(samples_struct.noisevect.begin(), samples_struct.noisevect.end());
	//	size_tmp = it - samples_struct.noisevect.begin();
	//
	//	if(size_tmp==1){
	//		n_iter=1;
	//		if(rank==0)
	//			cout << "The same covariance Matrix will be inverted for all the scans\n" << endl;
	//	}else{
	n_iter = (long)samples_struct.noisevect.size();
	if(n_iter==0){
		if(rank==0)
			cerr << "WARNING. You have forgotten to mention covariance matrix in ini file\n";
#ifdef PARA_FRAME
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_CONFIG;
	}
	if(rank==0)
		cout << n_iter << " covariance Matrix will be inverted\n" << endl;
	//	}

	if(rank==0)
		cleanup_dirfile_saneInv(dir.tmp_dir, samples_struct, n_iter, noise_suffix, bolo_list);

#ifdef PARA_FRAME
	MPI_Barrier(MPI_COMM_WORLD); // other procs wait untill rank 0 has created dirfile architecture.
#endif

	// Open the dirfile to read temporary files
	string filedir = dir.tmp_dir + "dirfile";
	samples_struct.dirfile_pointer = gd_open((char *) filedir.c_str(), GD_RDWR | GD_VERBOSE
			| GD_UNENCODED);

	if (gd_error(samples_struct.dirfile_pointer) != 0) {
		cout << "error opening dirfile : " << filedir << endl;
#ifdef PARA_FRAME
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return 1;
	}

	if(rank==0)
		cout << "Inverting Covariance Matrices..." << endl;

	for (long iframe=samples_struct.iframe_min;iframe<samples_struct.iframe_max;iframe++){

		std::vector<string> channelIn; /* Covariance matrix channel vector*/

		// read covariance matrix in a fits file named fname
		// returns : -the bins => Ell
		// -the input channel list => channelIn
		// -The number of bins (size of Ell) => nbins
		// -The original NoiseNoise covariance matrix => RellthOrig
		read_CovMatrix(samples_struct.noisevect[iframe], channelIn, nbins, ell, RellthOrig);

		std::vector<string> channelOut; /* bolometer reduction : Reduced vector of output channel */

		channelOut = bolo_list[iframe];

		//Total number of detectors to ouput (if ndet< ndetOrig : bolometer reduction)
		ndet = channelOut.size();

		//Deal with bolometer reduction and fill Rellth and mixmat
		reorderMatrix(nbins, channelIn, RellthOrig, channelOut, Rellth);

		// Inverse reduced covariance Matrix : Returns iRellth
		inverseCovMatrixByMode(nbins, ndet, Rellth, iRellth);

		// write inversed noisePS in a binary file for each detector
		if(write_InvNoisePowerSpectra(samples_struct.dirfile_pointer, channelOut, nbins, ell, iRellth, samples_struct.basevect[iframe] + noise_suffix)){
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return EX_CANTCREAT;
		}

		// clean up
		gsl_matrix_free(Rellth);
		gsl_matrix_free(iRellth);
		gsl_matrix_free(RellthOrig);
		delete [] ell;


	} // n_iter loop

	if(rank==0)
		printf("done. \n");

#ifdef PARA_FRAME
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if (gd_close(samples_struct.dirfile_pointer))
		cout << "error closing dirfile : " << filedir << endl;

#ifdef PARA_FRAME
	MPI_Finalize();
#endif

	if(rank==0)
		printf("\nEnd of saneInv\n");

	return EXIT_SUCCESS;
}

