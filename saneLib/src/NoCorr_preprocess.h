/*
 * NoCorr_preprocess.h
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */

#ifndef NOCORR_PREPROCESS_H_
#define NOCORR_PREPROCESS_H_

#include <vector>
#include <string>

#include "todprocess.h"
#include "map_making.h"
#include "inline_IO2.h"

void do_PtNd_nocorr(double *PNd, string *extentnoiseSp_all, string noiseSppreffile,
		string dir,  string dirfile,
		 std::vector<string> bolonames, string *fits_table,
		double f_lppix, double f_lppix_Nk,
		double fsamp, long ntotscan, long addnpix, bool flgdupl, int factdupl,
		int fillg, long ns, long napod, long ndet,
		long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2, long long npix, long long npixsrc,
		bool NORMLIN, bool NOFILLGAP,bool remove_polynomia,  long iframe, double *S);



void do_PtNPS_nocorr(double *S, string *extentnoiseSp_all, string noiseSppreffile, string dir,
		string dirfile, std::vector<string> bolonames, double f_lppix,
		double fsamp, bool flgdupl, int factdupl, long ns,
		long ndet, long long *indpix, long NAXIS1, long NAXIS2, long long npix,
		long iframe, double *PtNPmatS, double *Mp, long *hits);

#endif /* NOCORR_PREPROCESS_H_ */
