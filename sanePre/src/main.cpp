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



int main(int argc, char *argv[])
/*! Sanepic preprocess main function */
{

	int size;/*!< number of processors */
	int rank;

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

	struct param_sanePre proc_param; /*! contains user options about preprocessing properties */
	struct samples samples_struct;  /*  everything about frames, noise files and frame processing order */
	struct param_sanePos pos_param; /*! contains user options about map projection and properties */
	struct param_common dir; /*! contains output input temp directories */

	// default parameters
	long iframe_min=0, iframe_max=0; /*!  min and max number of frame (used with mpi) */

	// those variables will not be used by sanePre but they are read in ini file (to check his conformity)
	struct param_sanePS structPS;
	struct param_saneInv saneInv_struct;
	struct param_sanePic struct_sanePic;
	string parser_output = "";

	// Processing time estimation
	time_t t2, t3;


	if (rank==0){ // root parse ini file and fill the structures. Also print warnings or errors

		int parsed=0; // parser error status


		if (argc<2)/* not enough argument */
			parsed=1;
		else {
			/* parse ini file and fill structures */
			parsed=parser_function(argv[1], parser_output, dir, samples_struct, pos_param, proc_param,
					structPS, saneInv_struct, struct_sanePic, size, rank);

		}


		// print parser warning and/or errors
		cout << endl << parser_output << endl;
		switch (parsed){/* error during parsing phase */

		case 1: printf("Please run %s using a *.ini file\n",argv[0]);
		break;

		case 2 : printf("Wrong program options or argument. Exiting !\n");
		break;

		case 3 : printf("Exiting...\n");
		break;

		default :;
		}


		// in case there is a parsing error or the dirfile format file was not created correctly
		if (parsed>0){
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return EX_CONFIG;
		}
	}


	// processing begins here
	t2=time(NULL);

#ifdef USE_MPI

	MPI_Datatype message_type;
	struct ini_var_strings ini_v;
	int ntotscan;

	if(rank==0){
		fill_var_sizes_struct(dir, pos_param, proc_param,
				saneInv_struct, structPS, samples_struct, ini_v);

		ntotscan = ini_v.ntotscan;
	}

	MPI_Bcast(&ntotscan, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(rank!=0){
		ini_v.fitsvect=new int[ntotscan];
		ini_v.noisevect=new int[ntotscan];
		ini_v.bolovect=new int[ntotscan];
	}

	ini_v.ntotscan=ntotscan;

	Build_derived_type_ini_var (&ini_v,	&message_type);

	MPI_Bcast(&ini_v, 1, message_type, 0, MPI_COMM_WORLD);

	commit_struct_from_root(dir, pos_param, proc_param, saneInv_struct, struct_sanePic, structPS, samples_struct, ini_v, rank);

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

	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				structPS, struct_sanePic, saneInv_struct);

		compute_dirfile_format_file(dir.tmp_dir, samples_struct, pos_param.fileFormat);
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(write_data_flag_to_dirfile(dir, samples_struct, iframe_min, iframe_max)){
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
		if(write_RA_DEC_to_dirfile(dir, samples_struct, iframe_min, iframe_max)){
			cout << "error write_RA_DEC_to_dirfile !! Exiting ...\n";
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			exit(EXIT_FAILURE);
		}

	//Get processing time
	t3=time(NULL);

	if(rank==0)
		cout << "Total Time : " << t3-t2 << " sec\n";

#ifdef USE_MPI

	delete [] ini_v.fitsvect;
	delete [] ini_v.noisevect;
	delete [] ini_v.bolovect;

	MPI_Finalize();
#endif

	if (rank == 0)
		cout << endl << "done." << endl;

	return EXIT_SUCCESS;

}
