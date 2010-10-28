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

#include "parse_FBFO.h"
#include "mpi_architecture_builder.h"
#include "struct_definition.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;



struct sortclass {
	bool operator() (int i,int j) { return (i<j);}
} sortobject;



int main(int argc, char *argv[])
{


	int size=1;
	int rank=0;

#ifdef USE_MPI
	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if(rank==0)
		cout << "size : " << size << " rank : " << rank << endl;

	if(size==1) {cerr << "Please run mpirun -n# with # > 1\n"; MPI_Barrier(MPI_COMM_WORLD); MPI_Finalize(); exit(1);}

	if (rank == 0){

		struct samples samples_struct;
		struct param_common dir;

		long *ruleorder ;
		long *frnum ;
		int parsed;

		string fname, parser_output = "";

		parsed=parse_FBFO(argv[1], parser_output, samples_struct,dir);

		// print parser warning and/or errors
		cout << endl << parser_output << endl;

		if (parsed==-1){
			if(rank==0)
				cerr << "Problem in the parse function. Exiting\n";
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			return EX_CONFIG;
		}


		if(samples_struct.ntotscan<size){
			cerr << "You must use at least " << size << " scans or reduce the number of processor to be used.\n";
			cerr << "This number must be at least equal to the number of scans you are using : ie. " << samples_struct.ntotscan << endl;
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			return EX_CONFIG;
		}


		if(samples_struct.scans_index.size()!=0){
			cerr << "You have already given processors order in the fits filelist. Exiting\n";
			sort (samples_struct.scans_index.begin(), samples_struct.scans_index.end(), sortobject);

			if (size!=(samples_struct.scans_index[samples_struct.scans_index.size()-1]+1))
				cout << "Warning, you have to run MPI with " << samples_struct.scans_index[samples_struct.scans_index.size()-1]+1 << " processors and you are currently running MPI with " << size << " processors\n";
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			return EX_CONFIG;

		}


		ruleorder     = new long[samples_struct.ntotscan];
		frnum         = new long[samples_struct.ntotscan+1];


		fname = dir.output_dir + parallel_scheme_filename;

		if(samples_struct.ntotscan==size){

			for(long hh=0; hh<samples_struct.ntotscan;hh++){
				ruleorder[hh]=hh;
				frnum[hh]=hh;
			}
			frnum[samples_struct.ntotscan]=frnum[samples_struct.ntotscan-1]+1;

			//write parallel schema in a file
			parsed=write_ParallelizationScheme(fname, ruleorder, frnum, size,samples_struct);
		}else{


			/********************* Define parallelization scheme   *******/
			find_best_order_frames(ruleorder, frnum, samples_struct.nsamples, samples_struct.ntotscan, size);

			//write parallel schema in a file
			parsed=write_ParallelizationScheme(fname, ruleorder, frnum, size,samples_struct);
		}

		if(parsed==-1){
			if(rank==0)
				cerr << "Write parallelization Error !" << endl;
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			return EX_CANTCREAT;
		}

		delete [] frnum;
		delete [] ruleorder;


	}


#else
	cout << "Mpi is not used for this step. Exiting" << endl;
	return EXIT_SUCCESS;
#endif

#ifdef USE_MPI
	//wait for the other processors
	MPI_Barrier(MPI_COMM_WORLD);
	// Close MPI process
	MPI_Finalize();
#endif

	if(rank==0)
		cout << "\nEnd of saneFrameOrder\n";
	return EXIT_SUCCESS;
}
