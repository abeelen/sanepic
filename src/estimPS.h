/*
 * estimPS.h
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */

#ifndef ESTIMPS_H_
#define ESTIMPS_H_


#include "todprocess.h"
#include "map_making.h"
#include "sane_io.h"



void EstimPowerSpectra(double fsamp, long ns, long ff, long ndet, int nn, long npix, long napod,
		long iframe, bool flgdupl, int factdupl, long *indpix,
		double *S, string MixMatfile, string *bolonames, string dirfile, string bextension,
		string fextension, /*string cextension,*/ int shift_data_to_point, string dir,
		string termin, bool NORMLIN, bool NOFILLGAP, string noiseSppreffile,
		string extentnoiseSp, string outdirSpN);

double fdsf(double **Rellexp, double *w, double **A, double **P, double **N, long ndet, long ncomp, long nbins);

void rescaleAP(double **A, double **P, long ndet, long ncomp, long nbins);

#endif /* ESTIMPS_H_ */
