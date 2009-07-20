/*
 * estimPS.h
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */

#ifndef ESTIMPS_H_
#define ESTIMPS_H_


#include <vector>
#include <algorithm>
#include "todprocess.h"
#include "map_making.h"
#include "inline_IO2.h"
//#include "sane_io.h"

extern "C" {
#include <fitsio.h>
}

void EstimPowerSpectra(double fsamp, long ns, long ff, long ndet, int nn, long npix, long napod,
		long iframe, bool flgdupl, int factdupl, long *indpix,
		double *S, string MixMatfile, std::vector<string> bolonames, string dirfile, string bextension,
		string fextension, int shift_data_to_point, string dir,
		string termin, bool NORMLIN, bool NOFILLGAP, bool remove_polynomia, string noiseSppreffile,
		string extentnoiseSp, string outdirSpN);

double fdsf(double **Rellexp, double *w, double **A, double **P, double **N, long ndet, long ncomp, long nbins);

void rescaleAP(double **A, double **P, long ndet, long ncomp, long nbins);

#endif /* ESTIMPS_H_ */
