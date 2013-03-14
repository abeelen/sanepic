#ifndef ESTIMPS_STEPS_H_
#define ESTIMPS_STEPS_H_

#ifdef USE_MPI
#include "mpi.h"
#endif

double fdsf(double **Rellexp, double *w, double **A, double **P, double **N, long ndet, long ncomp, long nbins);

#ifdef USE_MPI
double fdsf_MPI(double **Rellexp, double *w, double **A, double **P, double **N, long ndet, long ncomp, long nbins, int bolo_rank, int bolo_size, MPI_Comm Comm);
#endif

void rndInitMixmat(long ndet, long ncomp, double ** &mixmat);

void rescaleAP(double **A, double **P, long ndet, long ncomp, long nbins);

double crudeSigma(double * data, long ns);

#endif /* ESTIMPS_STEPS_H_ */
