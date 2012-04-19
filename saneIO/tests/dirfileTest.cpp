#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <list>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <sysexits.h>

#include "mpi_architecture_builder.h"
#include "struct_definition.h"

#include "dataIO.h"
#include "imageIO.h"
#include "temporary_IO.h"
#include "inputFileIO.h"
#include "parser_functions.h"

extern "C" {
#include "nrutil.h"
#include "getdata.h"
}

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;

// template<typename T, size_t N>
// T * end(T (&ra)[N]) {
//     return ra + N;
// }
#define Elements_in(arrayname) (sizeof arrayname/sizeof *arrayname)

void print_struct(struct param_saneProc proc_param, struct samples samples_struct, struct param_sanePos pos_param, struct param_common dir,
		struct param_saneInv inv_param, struct param_sanePic pic_param, struct param_sanePS ps_param);

int main(int argc, char *argv[])
{


	int size;
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


	struct param_saneProc proc_param; /*! A structure that contains user options about preprocessing properties */
	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_sanePos pos_param; /*! A structure that contains user options about map projection and properties */
	struct param_common dir; /*! structure that contains output input temp directories */

	string field; /*! actual boloname in the bolo loop */
	struct param_sanePic pic_param;

	// those variables will not be used by sanePic but they are read in ini file (to check his conformity)
	struct param_sanePS ps_param;
	struct param_saneInv inv_param;
	string parser_output = "";


	// parser variables
	int parsed = 0;

	if (argc < 2) // no enough arguments
		parsed = 1;
	else {

		/* parse ini file and fill structures */
		parsed = parser_function(argv[1], parser_output, dir,
				samples_struct, pos_param, proc_param, ps_param, inv_param,
				pic_param, size, rank);


		if(rank==0)
			// print parser warning and/or errors
			cout << endl << parser_output << endl;


	}

	if (parsed > 0) { // error during parser phase
		if (rank == 0)
			switch (parsed) {

			case 1:
				printf(
						"Please run %s using the following options : sanepic_ini.ini \n",
						argv[0]);
				break;

			case 2:
				printf("Wrong program options or argument. Exiting !\n");
				break;

			case 3:
				printf("Exiting...\n");
				break;

			default:
				;
			}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_CONFIG;
	}


	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				ps_param, pic_param, inv_param);

	}

	cout << rank << " here" << endl;

#ifdef PARA_FRAME
	//	MPI_Barrier(MPI_COMM_WORLD);

	if(configure_PARA_FRAME_samples_struct(dir.tmp_dir, samples_struct, rank, size)){
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return EX_IOERR;
	}

	MPI_Barrier(MPI_COMM_WORLD);

#endif
	
	std::vector<DIRFILE *> dirfile_pointers;
	dirfile_pointers.clear();
	dirfile_pointers.resize(samples_struct.ntotscan, NULL);

	for (long iframe=samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++){
	  string dirfile = dir.tmp_dir+"dirfile/"+samples_struct.basevect[iframe];
	  dirfile_pointers[iframe] = gd_open((char *) dirfile.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED);
	  
	  const char * subdirs[] = { "fData", "Indexes", "data", "flag", "LON", "LAT", "Noise_data", "Noise_data/ell" };
	  for (unsigned long ii=0; ii< Elements_in(subdirs); ii++)
	    cout << rank << " " << iframe << " " << subdirs[ii] << endl;
	}

	
	


#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
	return 0;

}
