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
#include "mpi_architecture_builder.h"

extern "C" {
#include "nrutil.h"
}

//void do_PtNd(double *PNd, std::string *extentnoiseSp_all, std::string noiseSppreffile,
//		std::string dir, std::string prefixe, std::vector<std::string> bolonames,
//		double f_lppix, double fsamp,  long ns, long ndet,
//		long long *indpix, long NAXIS1, long NAXIS2, long long npix, long iframe,
//		double *Mp, long *hits);

//void write_ftrProcesdata(double *S, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2, long long npix,
//		long long npixsrc, long ntotscan, long long addnpix, bool flgdupl, int factdupl,
//		int fillg, std::string dir, std::string dirfile,
//		std::vector<std::string> bolonames,std::string *fits_table, double f_lppix,  long ns,
//		long napod, long ndet, bool NORMLIN, bool NOFILLGAP,bool remove_polynomia,
//		long iframe);

void write_ftrProcesdata(double *S, struct user_options u_opt, struct samples samples_struct, struct input_commons com,
		std::string tmp_dir,	struct detectors det, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2,
		long long npix,	long long npixsrc, long long addnpix, double f_lppix, long ns, long iframe);

void do_PtNd(double *PNd, std::string *extentnoiseSp_all, std::string dir, std::string prefixe,
		struct detectors det, double f_lppix, double fsamp, long ns,
		long long *indpix, long NAXIS1, long NAXIS2, long long npix, long iframe,
		double *Mp, long *hits);


//void write_tfAS(double *S, long long *indpix, long NAXIS1, long NAXIS2, long long npix,
//		bool flgdupl, int factdupl,
//		std::string dir, long ns, long ndet, long iframe, std::vector<std::string> bolonames);

void write_tfAS(double *S, struct detectors det,long long *indpix, long NAXIS1, long NAXIS2, long long npix,
		bool flgdupl, std::string dir, long ns, long iframe);

#endif /* CORR_PREPROCESS_H_ */
