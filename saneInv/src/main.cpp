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

#include "CovMatrixIO.h"
#include "SaneInvTools.h"
#include "TemporaryIO.h"
#include "InputFileIO.h"
#include "MPIConfiguration.h"
#include "ParserFunctions.h"
#include "StructDefinition.h"
#include "ErrorCode.h"

#include <gsl/gsl_matrix.h>


#ifdef USE_MPI
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

	int rank, size; /* MPI processor rank and MPI total number of used processors */
	int  bolo_rank, bolo_size;  /* As for parallel scheme */
	int node_rank, node_size; /* On a node basis, same as *sub* but for frame scheme */

#ifdef USE_MPI

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	MPI_Comm MPI_COMM_NODE, MPI_COMM_MASTER_NODE;
#else
	size = 1;
	rank = 0;
	bolo_size  = 1;
	bolo_rank  = 0;
	node_size  = 1;
	node_rank  = 0;
#endif

	if(rank==0)
		cout << endl << "saneInv : inversion of the noise power spectra" << endl;

	struct param_common dir;
	struct samples samples_struct;
	struct param_saneInv Inv_Param;

	struct param_sanePos Pos_param;
	struct param_saneProc Proc_param;
	struct param_sanePic Pic_param;
	struct param_sanePS PS_param;
	string parser_output = "";


	// data parameters
	/*
	 * -ndet = number of detectors to parser_output
	 * -ndetOrig = number of detectors in the NoiseNoise matrix
	 * -nbins = number of bins (Ell)
	 */
	long ndet, nbins;
	double *ell; /* bins values */

	/*
	 * -Rellth : Reduced NoiseNoise matrix
	 * -RellthOrig : Original NoiseNoise matrix
	 * -iRellth : Inverted reduced NoiseNoise matrix
	 * -mixmatOrig : original mixing matrix
	 * -mixmat : Reduced mixing matrix
	 */
	gsl_matrix *Rellth, *RellthOrig, *iRellth;

	string boloname;/* channels list file */

	std::vector<int> indexIn; /* bolometer index, used to determine which intput detector corresponds to which output detector*/

	uint16_t mask_saneInv = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | TMP_PATH_PROBLEM |
			BOLOFILE_NOT_FOUND | SANEINV_INPUT_ERROR | FITS_FILELIST_NOT_FOUND; // 0x601b

	uint16_t parsed=0x0000; // parser error status
	uint16_t compare_to_mask; // parser error status

	//	if (rank==0){ // root parse ini file and fill the structures. Also print warnings or errors

	// Parse ini file
	if (argc<2) {
		if (rank == 0)
			cerr << "EE - Please run  " << StringOf(argv[0]) << " with a .ini file" << endl;
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD );
		MPI_Finalize();
#endif
		exit(EXIT_FAILURE);

	} else {
		parsed=parser_function(argv[1], parser_output, dir, samples_struct, Pos_param, Proc_param,
				PS_param, Inv_Param, Pic_param, size, rank);

		compare_to_mask = parsed & mask_saneInv;

		// print parser warning and/or errors
		if (rank == 0)
			cout << endl << parser_output << endl;

		if(compare_to_mask>0x0000){

			switch (compare_to_mask){/* error during parsing phase */

			case 0x0001:
				if (rank==0)
					cerr << " EE - Please run " << StringOf(argv[0]) << " using a correct *.ini file" << endl;
				break;

			default :
				if (rank==0)
					cerr << "EE - Wrong program options or argument. Exiting ! " <<  "("<< hex << compare_to_mask << ")" << endl;
				break;
			}

#ifdef USE_MPI
			MPI_Finalize();
#endif
			return EX_CONFIG;
		}
	}

	/* ------------------------------------------------------------------------------------*/
	// Start...

	if(rank==0) {
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, Pos_param,  Proc_param,
				PS_param, Pic_param, Inv_Param);

	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);

	if(configureMPI(dir.output_dir, samples_struct, rank, size,
			bolo_rank,  bolo_size, node_rank, node_size,
			MPI_COMM_NODE, MPI_COMM_MASTER_NODE)){
		if (rank==0)
			cerr << endl << endl << "Exiting..." << endl;

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return EX_CONFIG;
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif


	if ( cleanup_dirfile_saneInv(dir.tmp_dir, samples_struct, bolo_rank) ) {
		cerr << "EE - Error in initializing dirfile - Did you run sanePre ?" << endl;
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_CONFIG;
	}


#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD); // other procs wait untill rank 0 has created dirfile architecture.
#endif

	if(rank==0)
		cout << endl << "Inverting Covariance Matrices..." << endl;

	if (bolo_rank == 0) {
		for (long iframe=samples_struct.iframe_min;iframe<samples_struct.iframe_max;iframe++){

			std::vector<string> channelIn; /* Covariance matrix channel vector*/

			// read covariance matrix in a fits file named fname
			// returns : -the bins => Ell
			// -the input channel list => channelIn
			// -The number of bins (size of Ell) => nbins
			// -The original NoiseNoise covariance matrix => RellthOrig
			read_CovMatrix(samples_struct.noisevect[iframe], channelIn, nbins, ell, RellthOrig);

			std::vector<string> channelOut; /* bolometer reduction : Reduced vector of parser_output channel */

			channelOut = samples_struct.bolo_list[iframe];

			//Total number of detectors to ouput (if ndet< ndetOrig : bolometer reduction)
			ndet = channelOut.size();

			//Deal with bolometer reduction and fill Rellth and mixmat
			reorderMatrix(nbins, channelIn, RellthOrig, channelOut, Rellth);

			// Inverse reduced covariance Matrix : Returns iRellth
			inverseCovMatrixByMode(nbins, ndet, Rellth, iRellth);

			// write inversed noisePS in a binary file for each detector
			if(write_InvNoisePowerSpectra(samples_struct.dirfile_pointers[iframe], channelOut, samples_struct.basevect[iframe], nbins, ell, iRellth)){
#ifdef USE_MPI
				MPI_Abort(MPI_COMM_WORLD, 1);
#endif
				return EX_CANTCREAT;
			}

			// clean up
			gsl_matrix_free(Rellth);
			gsl_matrix_free(iRellth);
			gsl_matrix_free(RellthOrig);
			delete [] ell;


		} // iframe loop
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// Close previously openened dirfile
	for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++){
		if (samples_struct.dirfile_pointers[iframe]) {
			if (gd_close(samples_struct.dirfile_pointers[iframe])){
				cerr << "EE - error closing dirfile...";
			} else {
				samples_struct.dirfile_pointers[iframe] = NULL;
			}
		}
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Comm_free(&MPI_COMM_NODE);
	MPI_Comm_free(&MPI_COMM_MASTER_NODE);

	MPI_Finalize();
#endif

	if(rank==0)
		cout << endl << "End of "<< StringOf(argv[0]) << endl;

	return EXIT_SUCCESS;
}

