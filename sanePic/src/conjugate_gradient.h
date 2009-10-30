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
#include "wcslib/wcs.h"
}


//void sanepic_conjugate_gradient(bool flgdupl, long long npix, double* &S, long iframe_min, long iframe_max,
//		long *nsamples, std::vector<double> fcut,double f_lp,double fsamp,
//		long long *indpix,struct wcsprm * wcs, long NAXIS1, long NAXIS2,
//		int factdupl, std::string poutdir, long ndet,
//		std::string *extentnoiseSp_all,std::string noiseSppreffile, std::vector<std::string> bolonames, int iterw,
//		long long *indpsrc, long long npixsrc, int flagon, bool projgaps, int rank, bool CORRon,
//		std::string dirfile, double *&PNdtot, long ntotscan,long long addnpix,bool NORMLIN,bool NOFILLGAP,
//		long napod,bool remove_polynomia, std::string outdir, std::string *fits_table);


void sanepic_conjugate_gradient(struct samples samples_struct,struct input_commons com,struct detectors det,
		struct directories dir,struct user_options u_opt, long long npix, double* &S,long iframe_min, long iframe_max,
		std::vector<double> fcut, long long *indpix, struct wcsprm * wcs, long NAXIS1, long NAXIS2,
		int iterw, long long *indpsrc, long long npixsrc, int flagon, int rank,
	    double *&PNdtot, long long addnpix);

#endif /* CONJUGATEGRAIDENT_H_ */
