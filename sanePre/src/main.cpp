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


extern "C" {
#include "nrutil.h"
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
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
		cout << endl << "sanePre :  Pre Processing of the data" << endl;

	struct param_sanePre proc_param; /* contains user options about preprocessing properties */
	struct samples samples_struct;  /*  everything about frames, noise files and frame processing order */
	struct param_sanePos pos_param; /* contains user options about map projection and properties */
	struct param_common dir; /* contains output input temp directories */

	// default parameters
	long iframe_min=0, iframe_max=0; /* min and max number of frame (used with mpi) */

	// those variables will not be used by sanePre but they are read in ini file (to check his conformity)
	struct param_sanePS structPS;
	struct param_saneInv saneInv_struct;
	struct param_sanePic struct_sanePic;
	string parser_output = "";

	std::vector<std::vector<std::string> > bolo_list; // this vector contains all bolonames for all the scans

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
		parsed=parser_function(argv[1], parser_output, dir, samples_struct, pos_param, proc_param,
				structPS, saneInv_struct, struct_sanePic, size, rank);

		compare_to_mask = parsed & mask_sanePre;


		// print parser warning and/or errors
		if (rank==0)
			cout << endl << parser_output << endl;

	}


	if(compare_to_mask>0x0000){

		switch (compare_to_mask){/* error during parsing phase */

		case 0x0001: printf("Please run %s using a correct *.ini file\n",argv[0]);
		break;

		default : printf("Wrong program options or argument. Exiting !\n");
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

	if(configure_PARA_FRAME_samples_struct(dir.output_dir, samples_struct, rank, size, iframe_min, iframe_max)){
		MPI_Abort(MPI_COMM_WORLD, 1);
		exit(EX_IOERR);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (iframe_max==iframe_min){ // ifram_min=iframe_max => This processor will not do anything
		cout << "Warning. Rank " << rank << " will not do anything ! please run saneFrameorder\n";
	}
#else

	iframe_min = 0; // single processor, will compute all the scans
	iframe_max = samples_struct.ntotscan;

#endif


//	if (rank == 0){
//		cout << rank << " : " << iframe_min << " " << iframe_max << endl;
//		for (long iframe = iframe_min; iframe < iframe_max; iframe++){
//			cout << iframe << " " << samples_struct.basevect[iframe] << endl;
//		}
//
//	}
	/* ------------------------------------- READ bolo list ----------------------------*/

	if(channel_list_to_vect_list(samples_struct, bolo_list, rank)){
		cout << "error in channel_list_to_vect_list" << endl;
		return EX_CONFIG;
	}

	/* ------------------------------------------------------------------------------------*/

	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				structPS, struct_sanePic, saneInv_struct);

		compute_dirfile_format_file(dir.tmp_dir, samples_struct, pos_param.fileFormat);
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(write_data_flag_to_dirfile(dir, samples_struct, iframe_min, iframe_max, bolo_list)){
		cout << "error write_data_flag_to_dirfile !! Exiting ...\n";
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		exit(EXIT_FAILURE);
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(pos_param.fileFormat==1)
		if(write_LON_LAT_to_dirfile(dir, samples_struct, iframe_min, iframe_max, bolo_list)){
			cout << "error write_LON_LAT_to_dirfile !! Exiting ...\n";
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			exit(EXIT_FAILURE);
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
