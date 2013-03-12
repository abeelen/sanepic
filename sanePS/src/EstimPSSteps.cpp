
#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

#include <vector>
#include <sstream>
#include <cmath>
#include <string>

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sysexits.h>

#include "TodProcess.h"
#include "TemporaryIO.h"
#include "InputFileIO.h"
#include "DataIO.h"
#include "CovMatrixIO.h"
#include "EstimPSSteps.h"
#include "sanePSIO.h"
#include "ErrorCode.h"
#include "Utilities.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>

extern "C" {
#include <fitsio.h>
#include "nrutil.h"
}


#ifdef USE_MPI
#include "mpi.h"
#endif

// #include <omp.h>

using namespace std;


double fdsf(double **Rellexp, double *w, double **A, double **P, double **N, long ndet, long ncomp, long nbins){

	double f;

	double triRhR, logdetiR;

	gsl_vector *uvec, *ivec;
	gsl_matrix *R;

	double **hR, **eR, **iR;

	hR = dmatrix(0,ndet-1,0,ndet-1);
	eR = dmatrix(0,ndet-1,0,ndet-1);
	iR = dmatrix(0,ndet-1,0,ndet-1);

	uvec = gsl_vector_alloc(ndet);
	ivec = gsl_vector_alloc(ndet);
	R    = gsl_matrix_alloc(ndet, ndet);

	// init
	f   = 0.0 ;
	long idet1, idet2, iComp;

	for (long iBin=0;iBin<nbins;iBin++) {

		for (idet1=0;idet1<ndet;idet1++)
			for (idet2=0;idet2<ndet;idet2++)
				hR[idet1][idet2] = Rellexp[idet1+idet2*ndet][iBin];

		/// building Rth from A, P, N

		gsl_matrix_set_zero(R);
		for(idet1=0; idet1< ndet; idet1++)
			gsl_matrix_set(R,idet1,idet1,N[idet1][iBin]);

		for (idet1=0;idet1<ndet;idet1++)
			for (idet2=0;idet2<ndet;idet2++)
				for (iComp=0;iComp<ncomp;iComp++)
					gsl_matrix_set(R, idet1, idet2, gsl_matrix_get(R, idet1, idet2) + A[idet1][iComp] * P[iComp][iBin] * A[idet2][iComp]);

		gsl_linalg_cholesky_decomp(R);

		for (idet1=0;idet1<ndet;idet1++){

			gsl_vector_set_basis(uvec,idet1);
			gsl_linalg_cholesky_solve(R, uvec, ivec);

			for (idet2=0;idet2<ndet;idet2++)
				iR[idet1][idet2] = gsl_vector_get(ivec,idet2);
		}

		triRhR = 0.0;
		for (idet1=0;idet1<ndet;idet1++)
			for (idet2=0;idet2<ndet;idet2++)
				triRhR += iR[idet1][idet2]*hR[idet2][idet1] ;

		logdetiR = 0;
		for (idet1=0;idet1<ndet;idet1++)
			logdetiR -= log(gsl_matrix_get(R,idet1,idet1)*gsl_matrix_get(R,idet1,idet1));

		f   +=  w[iBin] * (triRhR - logdetiR - ndet ) ; //pb when hR non inversible

	}

	gsl_vector_free(uvec);
	gsl_vector_free(ivec);
	gsl_matrix_free(R);

	free_dmatrix(hR,0,ndet-1,0,ndet-1);
	free_dmatrix(eR,0,ndet-1,0,ndet-1);
	free_dmatrix(iR,0,ndet-1,0,ndet-1);

	return f;

}

#ifdef USE_MPI

double fdsf_MPI(double **Rellexp, double *w, double **A, double **P, double **N, long ndet, long ncomp, long nbins, int bolo_rank, int bolo_size, MPI_Comm Comm){

	double f;

	double triRhR, logdetiR;

	gsl_vector *uvec, *ivec;
	gsl_matrix *R;

	double **hR, **eR, **iR;

	hR = dmatrix(0,ndet-1,0,ndet-1);
	eR = dmatrix(0,ndet-1,0,ndet-1);
	iR = dmatrix(0,ndet-1,0,ndet-1);

	uvec = gsl_vector_alloc(ndet);
	ivec = gsl_vector_alloc(ndet);
	R    = gsl_matrix_alloc(ndet, ndet);

	// init
	f   = 0.0;
	long idet1, idet2, iComp;

	for (long iBin = floor(bolo_rank*nbins*1.0/bolo_size); iBin<floor((bolo_rank+1)*nbins*1.0/bolo_size); iBin++) {

		for (idet1=0;idet1<ndet;idet1++)
			for (idet2=0;idet2<ndet;idet2++)
				hR[idet1][idet2] = Rellexp[idet1+idet2*ndet][iBin];

		/// building Rth from A, P, N

		gsl_matrix_set_zero(R);
		for(idet1=0; idet1< ndet; idet1++)
			gsl_matrix_set(R,idet1,idet1,N[idet1][iBin]);

		for (idet1=0;idet1<ndet;idet1++)
			for (idet2=0;idet2<ndet;idet2++)
				for (iComp=0;iComp<ncomp;iComp++)
					gsl_matrix_set(R, idet1, idet2, gsl_matrix_get(R, idet1, idet2) + A[idet1][iComp] * P[iComp][iBin] * A[idet2][iComp]);

		gsl_linalg_cholesky_decomp(R);

		for (idet1=0;idet1<ndet;idet1++){

			gsl_vector_set_basis(uvec,idet1);
			gsl_linalg_cholesky_solve(R, uvec, ivec);

			for (idet2=0;idet2<ndet;idet2++)
				iR[idet1][idet2] = gsl_vector_get(ivec,idet2);
		}

		triRhR = 0.0;
		for (idet1=0;idet1<ndet;idet1++)
			for (idet2=0;idet2<ndet;idet2++)
				triRhR += iR[idet1][idet2]*hR[idet2][idet1] ;

		logdetiR = 0;
		for (idet1=0;idet1<ndet;idet1++)
			logdetiR -= log(gsl_matrix_get(R,idet1,idet1)*gsl_matrix_get(R,idet1,idet1));

		f   +=  w[iBin] * (triRhR - logdetiR - ndet ) ; //pb when hR non inversible


	}
	MPI_Barrier(Comm);
	MPI_Allreduce(MPI_IN_PLACE, &f,  1, MPI_DOUBLE, MPI_SUM, Comm );

	gsl_vector_free(uvec);
	gsl_vector_free(ivec);
	gsl_matrix_free(R);

	free_dmatrix(hR,0,ndet-1,0,ndet-1);
	free_dmatrix(eR,0,ndet-1,0,ndet-1);
	free_dmatrix(iR,0,ndet-1,0,ndet-1);

	return f;

}

#endif



void rescaleAP(double **A, double **P, long ndet, long ncomp, long nbins){

	double *norm2ratio;

	norm2ratio = new double[ncomp];

	fill(norm2ratio,norm2ratio+ncomp,0.0);

	for (long ii=0;ii<ncomp;ii++){
		for (long jj=0;jj<ndet;jj++)
			norm2ratio[ii] += A[jj][ii] * A[jj][ii];
		norm2ratio[ii] = 1.0/norm2ratio[ii];
	}

	for (long ii=0;ii<ncomp;ii++){
		for (long jj=0;jj<ndet;jj++)
			A[jj][ii] = A[jj][ii] * sqrt(norm2ratio[ii]) ;
		for (long ib=0;ib<nbins;ib++)
			P[ii][ib] = P[ii][ib] / norm2ratio[ii] ;
	}

	delete [] norm2ratio;

}

// TODO: on 500 samples only in the middle of the timeline... is there a better/simpler/faster way...
// TODO: Stupid way of computing sigma of the noise
double crudeSigma(double * data, long ns){

	double mm      = 0.0;
	double tmpsign = 0.0;
	for (long ii=ns/2;ii<ns/2+500;ii++) mm += data[ii];
	mm = mm/501.0;
	for (long ii=ns/2;ii<ns/2+500;ii++) tmpsign += gsl_pow_2(data[ii]-mm);
	return sqrt(tmpsign/500.0);

}
