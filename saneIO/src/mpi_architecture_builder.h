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

#include "stdio.h"

using namespace std;

// scans are distributed over processors
void find_best_order_frames(long *pos, long *frnum, long *ns, long ntotscan, int size);
double* randg(long nombre, int seedpass); // on garde

void write_ParallelizationScheme(string fname, long  *pos, long  *frnum, long  *ns,  long  ntotscan, int  size);
void read_ParallelizationScheme(string fname,  long **pos, long **frnum, long **ns,  long *ntotscan, int *size);
void check_ParallelizationScheme(string fname, long *ns, long ntotscan, int size, long **pos, long **frnum);

#endif /* MPI_ARCHITECTURE_BUILDER_H_ */
