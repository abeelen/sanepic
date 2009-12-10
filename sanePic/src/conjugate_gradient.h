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

void sanepic_conjugate_gradient(struct samples samples_struct,struct detectors det,
		struct directories dir,struct param_process proc_param, long long npix, double* &S,long iframe_min, long iframe_max,
		std::vector<double> fcut, long long *indpix, struct wcsprm * wcs, long NAXIS1, long NAXIS2,
		int iterw, long long *indpsrc, long long npixsrc, int flagon, int rank,
	    double *&PNdtot, long long addnpix);

#endif /* CONJUGATEGRAIDENT_H_ */
