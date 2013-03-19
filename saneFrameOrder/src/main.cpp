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

#include "MPIConfiguration.h"
#include "StructDefinition.h"
#include "ParserFunctions.h"
#include "ErrorCode.h"
#include "FrameOrder.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;

int main(int argc, char *argv[]) {

	int size = 1;
	int rank = 0;

#ifdef USE_MPI
	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0)
		cout << endl << "saneFrameOrder : distributing the frames over nodes" << endl << endl;

//	if (size == 1) {
//		cerr << "EE - using saneFrameOrder with only 1 CPU is meaningless" << endl;
//		cerr << "EE - Please run mpirun -n# with # > 1" << endl << endl;
//		MPI_Barrier(MPI_COMM_WORLD );
//		MPI_Finalize();
//		exit(EXIT_FAILURE);
//	}

	// To determine node ID based on processor names
	char *proc_names = NULL;
	int *nodeID = NULL;
	vector<long> nodeSizes;
	vector<string> nodeName;

	// struct needed
	string parser_output = "";
	struct param_common dir; /* structure that contains output input temp directories */
	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */

	// those variables will not be used by sanePic but they are read in ini file (to check his conformity)
	struct param_saneProc Proc_param; /* A structure that contains user options about preprocessing properties */
	struct param_sanePos Pos_param; /* A structure that contains user options about map projection and properties */
	struct param_sanePS PS_param;
	struct param_sanePic Pic_param;
	struct param_saneInv Inv_param;

	struct param_sanePS structPS;
	struct param_saneInv saneInv_struct;
	struct param_sanePic struct_sanePic;
	struct param_saneProc proc_param;
	struct param_sanePos pos_param;

	vector<long> order;


	uint32_t parsed = 0; // parser error status
	uint32_t compare_to_mask; // parser error status

	//	uint32_t mask_sanePre = 0x405f;
	uint32_t mask_saneFrameOrder = INI_NOT_FOUND | ARG_INPUT_PROBLEM
			| DATA_INPUT_PATHS_PROBLEM | TMP_PATH_PROBLEM
			| FITS_FILELIST_NOT_FOUND; // 0x405f

	if (argc < 2){
		if (rank == 0)
			cerr << "EE - Please run  " << argv[0] << " with a .ini file" << endl;
		MPI_Barrier(MPI_COMM_WORLD );
		MPI_Finalize();
		exit(EXIT_FAILURE);
	} else {
		/* parse ini file and fill structures */
		parsed = parser_function(argv[1], parser_output, dir,
				samples_struct, Pos_param, Proc_param, PS_param, Inv_param,
				Pic_param, size, rank);

		if (rank == 0)
			cout << parser_output << endl;

		if (rank == 0 && parsed == 0) {
			if (samples_struct.framegiven) {
				parser_output +=
						"EE - You have already given processors order in the fits file list.\n";
				parser_output += "EE - Exiting\n";
				parsed |= ARG_INPUT_PROBLEM;
			}

		}
		compare_to_mask = parsed & mask_saneFrameOrder;

		if (compare_to_mask > OK) {

			switch (compare_to_mask) {/* error during parsing phase */

			case 0x0001:
				if (rank == 0 )
					cerr << "Please run with a correct ini file" << endl;
				break;

			default:
				if (rank == 0)
					cout << "Wrong program options or argument. Exiting ! " << "("
					<< hex << compare_to_mask << ")" << endl;
				break;
			}

			MPI_Finalize();
			return EX_CONFIG;
		}
	}

	// Start...

	if (rank == 0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, Pos_param,  Proc_param,
				PS_param, Pic_param, Inv_param);

		removeProcName(dir.output_dir);

	}
	if (rank == 0){
		cout << endl << "II - Reading file sizes" << endl;
		readFramesFromFits(samples_struct);
	}


	if ( checkProcName( rank, size, dir.output_dir) )
		MPI_Abort(MPI_COMM_WORLD, -1);

	if (rank == 0)
		proc_names = new char[size*MPI_MAX_PROCESSOR_NAME];

	gatherProcName(rank, size, proc_names);

	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0){
		cout << "II - Distributing files" << endl;

		string output;

		vector<float> nodeWeight;
		map<string, float> map_nodeWeight;

		if (readNodeWeight(output, dir.input_dir, map_nodeWeight))
			cout << output << endl;

		// Find all computers
		nodeID = new int[size];
		AssignNodeByProcname(proc_names, size, nodeID, nodeName, nodeSizes);
		delete [] nodeID;

		switch(samples_struct.parallel_scheme){

		case 0:{
					// ... in the bolo/mixed case, each computer become a node with a nodeSizes cpus ...

					// ... assign the previously read weights to the nodeName....
					AssignNodeFloatWeight(nodeName, map_nodeWeight, nodeWeight);

					// ... and distribute the nsamples !
					distributeFrames(samples_struct.nsamples, nodeSizes, nodeWeight, order);
					break;
				}

//		case 0: {
//			// if parallel_scheme is bolo, when need to insure a few things...
//
//			// ... assign the previously read weights to the nodeName....
//			AssignNodeFloatWeight(nodeName, map_nodeWeight, nodeWeight);
//
//			// ... find the maximum weight ...
//			vector<float>::iterator it = find( nodeWeight.begin(), nodeWeight.end(), *( max_element(nodeWeight.begin(), nodeWeight.end()) ) );
//
//			size_t idNode = distance(nodeWeight.begin(), it);
//
//			order.clear();
//			order.assign(samples_struct.nsamples.size(), idNode);
//			break;
//		}
		case 1: {
			//  .. in the para case, each proc become a node of size 1 (1 proc) ...
			nodeName.clear();
			nodeName.resize(size);
			for (int ii=0; ii< size; ii++)
				nodeName[ii] = proc_names+(ii*MPI_MAX_PROCESSOR_NAME);

			nodeSizes.clear();
			nodeSizes.assign(size, 1);

			// ... assign the previously read weights to the nodeName....
			AssignNodeFloatWeight(nodeName, map_nodeWeight, nodeWeight);

			// ... and distribute the nsamples !
			distributeFrames(samples_struct.nsamples, nodeSizes, nodeWeight, order);
			break;
		}

		}

		delete [] proc_names;

		printNodeUsage(nodeName, nodeSizes, order);

		// used_nodeSizes will contains the given nodeSize or 0 if one node is not used...
		vector<long> used_nodeSizes;
		used_nodeSizes.assign(nodeSizes.size(), 0);
		for (size_t ii=0; ii < order.size(); ii++)
			used_nodeSizes[order[ii]] = nodeSizes[order[ii]];


		if ( *min_element(used_nodeSizes.begin(), used_nodeSizes.end()) == 0 ){
			cerr << "EE - Exiting..." << endl;
			MPI_Abort(MPI_COMM_WORLD, -1);
		}


		if ( writeParallelScheme(dir.output_dir, order, samples_struct) )
			MPI_Abort(MPI_COMM_WORLD, -1);
	}

	MPI_Barrier(MPI_COMM_WORLD );

	// Close MPI process
	MPI_Finalize();

	if(rank==0)
		cout << endl << "End of "<< StringOf(argv[0]) << endl;

	return EXIT_SUCCESS;

#else
	cerr << "WW -saneFrameOrder is only useful when used with mpi" << endl;
	cerr << "WW - Exiting" << endl;
	return EXIT_SUCCESS;
#endif

}
