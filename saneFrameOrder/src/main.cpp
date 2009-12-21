/*
 * main.cpp
 *
 *  Created on: 15 juin 2009
 *      Author: matthieu
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <cstdlib>

#include <fcntl.h>
#include <unistd.h>

#include <vector>
#include <stdio.h>
#include <algorithm>

#include "parse_FBFO.h"
#include "mpi_architecture_builder.h"

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
	// int tag = 10;
//	MPI_Status status;

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	cout << "size : " << size << " rank : " << rank << endl;

	if(size==1) {cerr << "Please run mpirun -n# with # > 1\n"; MPI_Barrier(MPI_COMM_WORLD); MPI_Finalize(); exit(1);}

	if (rank == 0){

		long ntotscan;
		long *nsamples; // number of samples table nsamples_vec -> nsamples
		long *ruleorder ;
		long *frnum ;
		int parsed;

		std::vector<string> fitsvect;
		std::vector<string> noisevect;
		std::vector<int> scans_index;
		//	vector<long>::iterator it;

		//std::vector<long> fframes_vec, nsamples_vec;
		string fname,outdir;

		parsed=parse_FBFO(argv[1], outdir, ntotscan, nsamples, fitsvect, noisevect, scans_index);
		if (parsed==-1){
			cerr << "Problem in the parse function. Exiting\n";
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			exit(0);
		}

		//cout << fitsvect[0] << " "  << fitsvect[1] << " "  << fitsvect[2] << " "  << fitsvect[3] << endl;
		//cout << noisevect[0] << " "  << noisevect[1] << " "  << noisevect[2] << " "  << noisevect[3] << endl;

		if(scans_index.size()!=0){
		  cerr << "You have already given processors order in the fits filelist. Exiting\n";
		  sort (scans_index.begin(), scans_index.end(), sortobject);
		  cout << scans_index[0] << " " <<  scans_index[1] << " " << scans_index[2] << " " << scans_index[3] << endl;

		  //scans_index.unique();
		  if (size!=(scans_index[scans_index.size()-1]+1))
		    cout << "Warning, you have to run MPI with " << scans_index[scans_index.size()-1]+1 << " processors and you are currently running MPI with " << size << " processors\n";
		  MPI_Barrier(MPI_COMM_WORLD);
		  MPI_Finalize();
		  exit(0);

		}

		//cout << ntotscan << endl;

		//nsamples      = new long[ntotscan];
		ruleorder     = new long[ntotscan];
		frnum         = new long[ntotscan+1];

		//vector2array(nsamples_vec,nsamples);

		fname = outdir + parallel_scheme_filename;

		/********************* Define parallelization scheme   *******/
		//		cout << ntotscan << endl;
		//		cout << nsamples[0] << " " << nsamples[1] << " " << nsamples[2] << " " << nsamples[3] <<  endl;
		//getchar();
		/*for (int ii=0; ii<ntotscan; ii++) {
			nsamples[ii] *= 20;      // convert nframes to nsamples
		}*/
		//cout << nsamples[0] << " " << nsamples[1] << " " << nsamples[2] << " " << endl;
		// reorder nsamples
		find_best_order_frames(ruleorder, frnum, nsamples, ntotscan, size);
		//cout << "ruleorder : " << ruleorder[0] << " " << ruleorder[1] << " " << ruleorder[2] << " \n";

		//		cout << size << endl;
		//		cout << ntotscan << endl;
		//		cout << nsamples[0] << " " << nsamples[1] << " " << nsamples[2] << " " << nsamples[3] << endl;
		//		cout << ruleorder[0] << " " << ruleorder[1] << " " << ruleorder[2] << " " << ruleorder[3] << " \n";
		//		cout << frnum[0] << " " << frnum[1] << " " << frnum[2] << " " << frnum[3] << " " << frnum[4] << endl;
		//cout << fitsvect[0] << " "  << fitsvect[1] << " "  << fitsvect[2] << " "  << fitsvect[3] << endl;
		//cout << noisevect[0] << " "  << noisevect[1] << " "  << noisevect[2] << " "  << noisevect[3] << endl;


		//write parallel schema in a file
		  parsed=write_ParallelizationScheme(fname, ruleorder, frnum, nsamples, ntotscan, size, fitsvect, noisevect, scans_index);
		if(parsed==-1){
			cerr << "merde" << endl;
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			exit(0);
		}
		delete [] frnum;
		delete [] ruleorder;


	}


#else
	cout << "Mpi is not used for this step. Exiting" << endl;
	exit(0);
#endif

#ifdef USE_MPI
	//wait for the other processors
	MPI_Barrier(MPI_COMM_WORLD);
	// Close MPI process
	MPI_Finalize();
#endif
}
