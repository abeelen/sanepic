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
#include "mpi_architecture_builder.h"



//void do_PtNd_nocorr(double *PNd, std::string *extentnoiseSp_all, std::string noiseSppreffile,
//		std::string dir,  std::string dirfile,
//		 std::vector<std::string> bolonames, std::string *fits_table,
//		double f_lppix, double f_lppix_Nk,
//		double fsamp, long ntotscan, long addnpix, bool flgdupl, int factdupl,
//		int fillg, long ns, long napod, long ndet,
//		long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2, long long npix, long long npixsrc,
//		bool NORMLIN, bool NOFILLGAP,bool remove_polynomia,  long iframe, double *S);

void do_PtNd_nocorr(double *PNd,std::string tmp_dir, struct user_options u_opt,
		struct samples samples_struct, struct input_commons com, struct detectors det, double f_lppix, double f_lppix_Nk,
		long addnpix, long ns, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2, long long npix,
		long long npixsrc, long iframe, double *S);


//void do_PtNPS_nocorr(double *S, std::string *extentnoiseSp_all, std::string noiseSppreffile, std::string dir,
//		std::string dirfile, std::vector<std::string> bolonames, double f_lppix,
//		double fsamp, bool flgdupl, int factdupl, long ns,
//		long ndet, long long *indpix, long NAXIS1, long NAXIS2, long long npix,
//		long iframe, double *PtNPmatS, double *Mp, long *hits);

void do_PtNPS_nocorr(double *S, std::string *extentnoiseSp_all, struct directories dir,
		struct detectors det,double f_lppix,double fsamp, bool flgdupl, long ns,
		long long *indpix, long NAXIS1, long NAXIS2, long long npix,
		long iframe, double *PtNPmatS, double *Mp, long *hits);

#endif /* NOCORR_PREPROCESS_H_ */
