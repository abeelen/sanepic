/*
 * Corr_preprocess.h
 *
 *  Created on: 20 juil. 2009
 *      Author: matthieu
 */

#ifndef CORR_PREPROCESS_H_
#define CORR_PREPROCESS_H_

#include <string>
#include <vector>
#include "struct_definition.h"

extern "C" {
#include "nrutil.h"
}

void write_ftrProcesdata(double *S, struct param_process proc_param, struct samples samples_struct, struct param_positions pos_param,
		std::string tmp_dir,	struct detectors det, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2,
		long long npix,	long long npixsrc, long long addnpix, double f_lppix, long ns, long iframe, int rank, int size);

void do_PtNd(double *PNd, std::string *extentnoiseSp_all, std::string dir, std::string prefixe,
		struct detectors det, double f_lppix, double fsamp, long ns, int rank, int size,
		long long *indpix, long NAXIS1, long NAXIS2, long long npix, long iframe,
		double *Mp, long *hits);


//void write_tfAS(double *S, long long *indpix, long NAXIS1, long NAXIS2, long long npix,
//		bool flgdupl, int factdupl,
//		std::string dir, long ns, long ndet, long iframe, std::vector<std::string> bolonames);

void write_tfAS(double *S, struct detectors det,long long *indpix, long NAXIS1, long NAXIS2, long long npix,
		bool flgdupl, std::string dir, long ns, long iframe, int rank, int size);

#endif /* CORR_PREPROCESS_H_ */
