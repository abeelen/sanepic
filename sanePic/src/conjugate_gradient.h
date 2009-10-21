/*
 * conjugate graident.h
 *
 *  Created on: 3 juil. 2009
 *      Author: matthieu
 */

#ifndef CONJUGATEGRAIDENT_H_
#define CONJUGATEGRAIDENT_H_

#include <vector>
#include <string>

extern "C" {
#include <fftw3.h>
#include "wcslib/wcs.h"
}



using namespace std;

void sanepic_conjugate_gradient(bool flgdupl, long long npix, double* &S, long iframe_min, long iframe_max,
		long *nsamples, std::vector<double> fcut,double f_lp,double fsamp,
		long long *indpix,
		struct wcsprm * wcs, long NAXIS1, long NAXIS2,
		int factdupl, string poutdir, long ndet,
		string *extentnoiseSp_all,string noiseSppreffile, std::vector<string> bolonames, int iterw,
		long long *indpsrc, long long npixsrc, int flagon, bool projgaps, int rank, bool CORRon,
		string dirfile, double *&PNdtot, long ntotscan,long long addnpix,bool NORMLIN,bool NOFILLGAP,
		long napod,bool remove_polynomia, string outdir, string *fits_table);

#endif /* CONJUGATEGRAIDENT_H_ */
