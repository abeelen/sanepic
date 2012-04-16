#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <cstdlib>

#include <fcntl.h>
#include <unistd.h>
#include <sysexits.h>
#include <vector>
#include <stdio.h>
#include <algorithm>

#include "mpi_architecture_builder.h"
#include "struct_definition.h"
#include "parser_functions.h"
#include "error_code.h"

#ifdef PARA_BOLO
#define PARA_FRAME
#endif

#ifdef PARA_FRAME
#include "mpi.h"
#endif

using namespace std;


int main(int argc, char *argv[])
{

	int size=1;
	int rank=0;

#ifdef PARA_FRAME
	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);


	if(rank==0)
		cout << endl << "saneFrameOrder :  distributing the data frames" << endl;

	if(size==1) {cerr << "Please run mpirun -n# with # > 1\n"; MPI_Barrier(MPI_COMM_WORLD); MPI_Finalize(); exit(1);}


	struct samples samples_struct;
	struct param_common dir;

	// struct used in the parser
	struct param_sanePS structPS;
	struct param_saneInv saneInv_struct;
	struct param_sanePic struct_sanePic;
	struct param_saneProc proc_param;
	struct param_sanePos pos_param;

	long *ruleorder ;
	long *frnum ;

	string parser_output = "";

	uint16_t parsed=0x0000; // parser error status
	uint16_t compare_to_mask; // parser error status

	//	uint16_t mask_sanePre = 0x405f;
	uint16_t mask_saneFrameOrder = INI_NOT_FOUND | ARG_INPUT_PROBLEM | DATA_INPUT_PATHS_PROBLEM | TMP_PATH_PROBLEM | FITS_FILELIST_NOT_FOUND; // 0x405f

	if (argc<2)/* not enough argument */
		compare_to_mask=0x001;
	else {
		/* parse ini file and fill structures */
		parsed=parser_function(argv[1], parser_output, dir, samples_struct, pos_param, proc_param, structPS, saneInv_struct, struct_sanePic, size, rank);

		if (rank == 0) {
			if(samples_struct.ntotscan<size){
				cerr << "WW - You are using more processors ("<< size << ") than avaible scans (" << samples_struct.ntotscan << ")" << endl;
				cerr << "WW - This will return non-optimal use of your processors" << endl;
			}


			if(samples_struct.framegiven){
				parser_output += "EE - You have already given processors order in the fits file list.\n";
				parser_output += "EE - Exiting\n";
				parsed |= ARG_INPUT_PROBLEM;
			}

		}
		compare_to_mask = parsed & mask_saneFrameOrder;

		// print parser warning and/or errors
		if (rank==0)
			cout << endl << parser_output << endl;

	}

	if (rank == 0) {

	if(compare_to_mask>OK){

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

	}

	if (rank == 0){

		ruleorder     = new long[samples_struct.ntotscan];
		frnum         = new long[size+1];

		if(samples_struct.ntotscan==size){ // special case : number of proc = number of scan

			for(long hh=0; hh<samples_struct.ntotscan;hh++){ // one scan per proc ...
				ruleorder[hh]=hh;
			}
			for(long hh=0; hh<size+1; hh++){
				frnum[hh] = hh;
			}

			frnum[samples_struct.ntotscan]=frnum[samples_struct.ntotscan-1]+1;

			//write parallel schema in a file
			parsed=write_ParallelizationScheme(dir.tmp_dir, ruleorder, frnum, size,samples_struct);


		}else{ // less procs than number of scans

			/********************* Define parallelization scheme   *******/
			find_best_order_frames(ruleorder, frnum, samples_struct.nsamples, samples_struct.ntotscan, size);

			//write parallel schema in a file
			parsed=write_ParallelizationScheme(dir.tmp_dir, ruleorder, frnum, size,samples_struct);
		}

		if(parsed==-1)
			cerr << "Write parallelization Error !" << endl;

		delete [] frnum;
		delete [] ruleorder;

		cout << endl << "End of saneFrameOrder" << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Close MPI process
	MPI_Finalize();
	return EXIT_SUCCESS;

#else
	cout << "Mpi is not used for this step. Exiting" << endl;
	return EXIT_SUCCESS;
#endif


}
