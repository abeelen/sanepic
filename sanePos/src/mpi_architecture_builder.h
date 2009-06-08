/*
 * mpi_architecture_builder.h
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */

#ifndef MPI_ARCHITECTURE_BUILDER_H_
#define MPI_ARCHITECTURE_BUILDER_H_


#include <stdio.h>
#include <stdlib.h>

// scans are distributed over processors
void find_best_order_frames(long *pos, long *frnum, long *ns, int ntotscan, int size);
double* randg(long nombre, int seedpass); // on garde


#endif /* MPI_ARCHITECTURE_BUILDER_H_ */
