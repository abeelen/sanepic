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


template<class T> void vector2array(std::vector<T> l, T* a)
{
	// copy list of type T to array of type T
	typename std::vector<T>::iterator iter;
	int i;

	for (iter=l.begin(), i=0; iter != l.end(); iter++, i++) {
		a[i] = *iter;
	}
}




int main(int argc, char *argv[])
{



	int size;
	int rank;

#ifdef USE_MPI
	// int tag = 10;
	//MPI_Status status;

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//cout << size << endl;
	//cout << rank << endl;


	//test params
	//size=2;
	//rank=0;

	if(size==1) {cerr << "Please run mpirun -n# first\n"; MPI_Finalize(); exit(1);}

	if (rank == 0){

		char * pPath;
		pPath = getenv ("TMPBATCH");
		if (pPath!=NULL)
			printf ("The current path is: %s\n",pPath);

		long ntotscan;
		long *nsamples ; // number of samples table nsamples_vec -> nsamples
		long *ruleorder ;
		long *frnum ;
		int parsed;
		std::vector<long> fframes_vec, nsamples_vec;
		string fname;

		parsed=parse_FBFO(argv[1], fname, ntotscan, fframes_vec, nsamples_vec);
		if (parsed==-1){
			cerr << "Problem in the parse function. Exiting\n";
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			exit(0);
		}

		nsamples      = new long[ntotscan];
		ruleorder     = new long[ntotscan];
		frnum         = new long[ntotscan+1];

		vector2array(nsamples_vec,nsamples);



		/********************* Define parallelization scheme   *******/
		cout << ntotscan << endl;
		cout << nsamples[0] << " " << nsamples[1] << " " << nsamples[2] << " " << endl;
		//getchar();
		for (int ii=0; ii<ntotscan; ii++) {
			nsamples[ii] *= 20;      // convert nframes to nsamples
		}
		//cout << nsamples[0] << " " << nsamples[1] << " " << nsamples[2] << " " << endl;
		// reorder nsamples
		find_best_order_frames(ruleorder, frnum, nsamples, ntotscan, size);
		//cout << "ruleorder : " << ruleorder[0] << " " << ruleorder[1] << " " << ruleorder[2] << " \n";

		cout << size << endl;
		cout << ntotscan << endl;
		cout << nsamples[0] << " " << nsamples[1] << " " << nsamples[2] << " " << endl;
		cout << ruleorder[0] << " " << ruleorder[1] << " " << ruleorder[2] << " \n";
		cout << frnum[0] << " " << frnum[1] << " " << frnum[2] << " " << frnum[3] << " " << endl;



		// write parallel schema in a file
		parsed=write_ParallelizationScheme(fname, ruleorder, frnum, nsamples, ntotscan, size);
		if(parsed==-1){
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			exit(0);
		}


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
