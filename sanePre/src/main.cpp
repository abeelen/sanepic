#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <time.h>
#include <cmath>
#include <vector>
#include <string>
#include <fftw3.h>
#include <sysexits.h>


#include "imageIO.h"
#include "temporary_IO.h"
#include "mpi_architecture_builder.h"
#include "parser_functions.h"
#include "struct_definition.h"
#include "inputFileIO.h"
#include "error_code.h"


extern "C" {
#include "nrutil.h"
#include "wcslib/wcs.h"     // to export wcs header keys
#include "wcslib/wcshdr.h"
#include <fftw3.h>          // for wisdom
}


#ifdef USE_MPI
#include "mpi.h"
#include <algorithm>
#include <fstream>
#endif

using namespace std;


//**********************************************************************************//
//**********************************************************************************//
//*************************** Beginning of main program ****************************//
//**********************************************************************************//
//**********************************************************************************//


/*! \mainpage Sanepic for SPIRE
 *
 * \section intro_sec Matthieu HUSSON & Alexandre Beelen
 *
 * Sanepic for SPIRE Manual
 */

/*!
 *  This is organized as :
 *
 *  - parse the input ini file and verify his validity
 *  - check for existence of directory/files pointed from the ini file
 *  - Print parser output to screen
 *
 *  - Generate, or clear dirfile tree (folders and format files)
 *  - for each file :
 *  	- Generate or clear the dirfile parts that will be filled : data, flag, LON, LAT
 *
 *  - Read all channel files, store it into a vector<vector> (and commit to other ranks if needed)
 *
 *  - Copy from scans (fits) to disk (getdata binaries in a dirfile tree) : data, flag, LON, LAT
 *
 */

int main(int argc, char *argv[])
/* Sanepic preprocess main function */
{

	int size; /* number of processors */
	int rank; /* rank = processor MPI rank*/

#ifdef USE_MPI

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
	size = 1;
	rank = 0;
#endif

	if(rank==0)
		cout << endl << "sanePre :  pre processing of the data" << endl;

	struct param_saneProc Proc_param; /* contains user options about preprocessing properties */
	struct samples samples_struct;  /*  everything about frames, noise files and frame processing order */
	struct param_sanePos Pos_param; /* contains user options about map projection and properties */
	struct param_common dir; /* contains output input temp directories */

	// those variables will not be used by sanePre but they are read in ini file (to check his conformity)
	struct param_sanePS PS_param;
	struct param_saneInv Inv_param;
	struct param_sanePic Pic_param;
	string parser_output = "";

	//	uint16_t mask_sanePre = 0x405f;
	uint16_t mask_sanePre = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | OUPUT_PATH_PROBLEM | TMP_PATH_PROBLEM |
			BOLOFILE_NOT_FOUND | FILEFORMAT_NOT_FOUND | FITS_FILELIST_NOT_FOUND; // 0x405f

#ifdef DEBUG
	// Processing time estimation
	time_t t2, t3;
#endif

	//	if (rank==0){ // root parse ini file and fill the structures. Also print warnings or errors

	uint16_t parsed=0x0000; // parser error status
	uint16_t compare_to_mask; // parser error status


	if (argc<2)/* not enough argument */
		compare_to_mask=0x001;
	else {
		/* parse ini file and fill structures */
		parsed=parser_function(argv[1], parser_output, dir, samples_struct, Pos_param, Proc_param, PS_param, Inv_param, Pic_param, size, rank);

		compare_to_mask = parsed & mask_sanePre;


		// print parser warning and/or errors
		if (rank==0)
			cout << endl << parser_output << endl;

	}


	if(compare_to_mask>0x0000){

		switch (compare_to_mask){/* error during parsing phase */

		case 0x0001: printf("Please run %s using a correct *.ini file\n",argv[0]);
		break;

		default : cout << "Wrong program options or argument. Exiting ! " <<  "("<< hex << compare_to_mask << ")" << endl;
		break;


		}

#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return EX_CONFIG;
	}
	//	}

#ifdef DEBUG
	// processing begins here
	t2=time(NULL);
#endif


#ifdef USE_MPI

	MPI_Barrier(MPI_COMM_WORLD);


	if(configure_PARA_FRAME_samples_struct(dir.tmp_dir, samples_struct, rank, size)){
		MPI_Abort(MPI_COMM_WORLD, 1);
		exit(EX_IOERR);
	}

	MPI_Barrier(MPI_COMM_WORLD);

#endif

	/* ------------------------------------------------------------------------------------*/

	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, Pos_param,  Proc_param,
				PS_param, Pic_param, Inv_param);

		compute_dirfile_format_file(dir.tmp_dir, samples_struct, Pos_param.fileFormat);
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(rank==0)
		cout << endl << "Exporting signal and flags... " << endl;

	if(write_data_flag_to_dirfile(dir, samples_struct)){
		cerr << "EE - write_data_flag_to_dirfile !! Exiting ..." << endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		exit(EXIT_FAILURE);
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif


	if(rank==0)
		cout << "Exporting positions... " << endl;


	//TODO: What if Pos_param = 0 ?????
	switch (Pos_param.fileFormat) {
	case 0:
		if(export_LON_LAT_to_dirfile(dir, samples_struct)){
			cerr << "EE - write_LON_LAT_to_dirfile !! Exiting ..." << endl;
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			exit(EXIT_FAILURE);
		}
		break;
	case 1:
		if(write_LON_LAT_to_dirfile(dir, samples_struct)){
			cerr << "EE - write_LON_LAT_to_dirfile !! Exiting ..." << endl;
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			exit(EXIT_FAILURE);
		}
		break;
	}


	// Create a fake WCS image header and populate it with info from the first fits file
	if (rank == 0){

		cout << "Exporting header keys... " << endl;

		int nwcs=1;                      // We will only deal with one wcs....
		struct wcsprm * wcs;             // wcs structure of the image

		char * subheader;               // Additionnal header keywords
		int nsubkeys;                   //


		if (get_fits_META(samples_struct.fitsvect[0], wcs, &subheader, &nsubkeys, rank)){
			cout << "pb getting fits META\n";
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		}

		if(save_keyrec(dir.tmp_dir,wcs, 0, 0, subheader, nsubkeys)){
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return(EX_CANTCREAT);
		}

		wcsvfree(&nwcs, &wcs);

	}


	if (rank == 0 && Proc_param.wisdom ){
		cout << "Building some wisdom..." << endl;

		fftw_complex * fdata;
		double *        data;
		long ns;
		fftw_plan plan;

		for (long iframe=0 ;iframe<samples_struct.ntotscan;iframe++){

			ns = samples_struct.nsamples[iframe];

			data  = (double *) fftw_malloc(sizeof(double)*ns);
			fdata = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (ns/2+1));

			plan = fftw_plan_dft_r2c_1d(ns, data, fdata, FFTW_PATIENT);
			//			fftw_print_plan(plan);
			fftw_destroy_plan(plan);

			plan = fftw_plan_dft_c2r_1d(ns, fdata, data, FFTW_PATIENT);
			//			fftw_print_plan(plan);
			fftw_destroy_plan(plan);

			fftw_free(data);
			fftw_free(fdata);
		}

		string filename = dir.tmp_dir + "fftw.wisdom";
		FILE * pFilename;

		pFilename = fopen((const char*) filename.c_str(), "w+");

		if ( pFilename != NULL ){
			fftw_export_wisdom_to_file( pFilename );
		} else {
			cerr << "EE - Problem while writing FFTW wisdom in " << filename << endl;
		}

		fclose(pFilename);

	}


#ifdef DEBUG
	//Get processing time
	t3=time(NULL);

	if(rank==0)
		cout << "Total Time : " << t3-t2 << " sec\n\n";

#endif

#ifdef USE_MPI
	MPI_Finalize();
#endif

	if(rank==0)
		cout << endl << "end of sanePre" << endl;

	return EXIT_SUCCESS;

}
