/*
 * mpi_architecture_builder.h
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */

#ifndef MPI_ARCHITECTURE_BUILDER_H_
#define MPI_ARCHITECTURE_BUILDER_H_

#include <vector>
#include <cstdlib>
#include <string>

#include "stdio.h"

//#include "mpi.h"

using namespace std;

// scans are distributed over processors
void find_best_order_frames(long *position, long *frnum, long *ns, long ntotscan, int size);
//double* randg(long nombre, int seedpass); // on garde

int compare_array_double (const void *array_1, const void *array_2);
double randg_archi(long nombre, int seedpass);

int write_ParallelizationScheme(string fname, long  *position, long  *frnum, long  *ns,  long  ntotscan, int  size);
void read_ParallelizationScheme(string fname,  long **position, long **frnum, long **ns,  long *ntotscan, int *size);
int check_ParallelizationScheme(string fname, long *ns, long ntotscan, int size, long **position, long **frnum);
int define_parallelization_scheme(int rank,string fname,long **frnum,long ntotscan,int size,long *nsamples,long *fframes);

long readFitsLength(string filename);
void readFrames(long * nScan, std::vector<string> &inputFiles, long *& fframes, long *& nsamples);
void read_fits_list(string fname, std::vector<string> &fitsfiles, std::vector<string> &noisefiles, std::vector<long> &frameorder, bool &framegiven);

#define parallel_scheme_filename  "parallel_scheme.bin";

#endif /* MPI_ARCHITECTURE_BUILDER_H_ */
