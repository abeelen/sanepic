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

	if(size==1) {cerr << "Please run mpirun -n# with # > 1\n"; MPI_Barrier(MPI_COMM_WORLD); MPI_Finalize(); exit(1);}

	if (rank == 0){

		struct samples samples_struct;
		struct param_common dir;

		// struct used in the parser
		struct param_sanePS structPS;
		struct param_saneInv saneInv_struct;
		struct param_sanePic struct_sanePic;
		struct param_saneProc proc_param;
		struct param_sanePos pos_param;

		sortclass_int sortobject;
		long *ruleorder ;
		long *frnum ;
		int parsed;
		int error_code = 0;

		string fname, parser_output = "";

		parsed=parser_function(argv[1], parser_output, dir, samples_struct, pos_param, proc_param,
				structPS, saneInv_struct, struct_sanePic, size, rank);


		// print parser warning and/or errors
		cout << endl << parser_output << endl;

		if (parsed==-1){
			cerr << "Problem in the parse function. Exiting\n";
		}



		if(samples_struct.ntotscan<size){
			cerr << "WW - You are using more processors ("<< size << ") than avaible scans (" << samples_struct.ntotscan << ")" << endl;
			cerr << "WW - This will return non-optimal use of your processors" << endl;
		}


		if(samples_struct.framegiven){
			cerr << "EE - You have already given processors order in the fits file list." << endl;
			cerr << "EE - Exiting" << endl;
			error_code=1;
		}

		if(error_code==0){

			ruleorder     = new long[samples_struct.ntotscan];
			frnum         = new long[size+1];


			fname = dir.output_dir + parallel_scheme_filename;

			if(samples_struct.ntotscan==size){ // special case : number of proc = number of scan

				for(long hh=0; hh<samples_struct.ntotscan;hh++){ // one scan per proc ...
					ruleorder[hh]=hh;
				}
				for(long hh=0; hh<size+1; hh++){
					frnum[hh] = hh;
				}

				frnum[samples_struct.ntotscan]=frnum[samples_struct.ntotscan-1]+1;

				//write parallel schema in a file
				parsed=write_ParallelizationScheme(fname, ruleorder, frnum, size,samples_struct);


			}else{ // less procs than number of scans

				/********************* Define parallelization scheme   *******/
				find_best_order_frames(ruleorder, frnum, samples_struct.nsamples, samples_struct.ntotscan, size);

				//write parallel schema in a file
				parsed=write_ParallelizationScheme(fname, ruleorder, frnum, size,samples_struct);
			}

			if(parsed==-1)
				cerr << "Write parallelization Error !" << endl;

			delete [] frnum;
			delete [] ruleorder;

		} // if error = 0

		cout << endl << "End of saneFrameOrder" << endl;
	} else {
		MPI_Barrier(MPI_COMM_WORLD);
	}

	// Close MPI process
	MPI_Finalize();
	return EXIT_SUCCESS;

#else
	cout << "Mpi is not used for this step. Exiting" << endl;
	return EXIT_SUCCESS;
#endif


}
