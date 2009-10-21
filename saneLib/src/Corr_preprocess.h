/*
 * Corr_preprocess.h
 *
 *  Created on: 20 juil. 2009
 *      Author: matthieu
 */

#ifndef CORR_PREPROCESS_H_
#define CORR_PREPROCESS_H_

#include "todprocess.h"
#include "map_making.h"
#include "inline_IO2.h"
#include <gsl/gsl_math.h>


extern "C" {
#include "nrutil.h"
}

void do_PtNd(double *PNd, string *extentnoiseSp_all, string noiseSppreffile,
		string dir, string prefixe, std::vector<string> bolonames,
		double f_lppix, double fsamp,  long ns, long ndet,
		long long *indpix, long NAXIS1, long NAXIS2, long long npix, long iframe,
		double *Mp, long *hits);

void write_ftrProcesdata(double *S, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2, long long npix,
		long long npixsrc, long ntotscan, long long addnpix, bool flgdupl, int factdupl,
		int fillg, string dir, string dirfile,
		std::vector<string> bolonames,string *fits_table, double f_lppix,  long ns,
		long napod, long ndet, bool NORMLIN, bool NOFILLGAP,bool remove_polynomia,
		long iframe);



void write_tfAS(double *S, long long *indpix, long NAXIS1, long NAXIS2, long long npix,
			bool flgdupl, int factdupl,
		string dir, long ns, long ndet, long iframe);


#endif /* CORR_PREPROCESS_H_ */
