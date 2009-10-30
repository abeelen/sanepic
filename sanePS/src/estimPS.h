/*
 * estimPS.h
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */

#ifndef ESTIMPS_H_
#define ESTIMPS_H_

#include <vector>
#include <string>
#include "mpi_architecture_builder.h"


//void EstimPowerSpectra(double fsamp, long ns, long ff,  long ndet, long NAXIS1, long NAXIS2, long long npix, long napod,
//		long iframe, bool flgdupl, int factdupl, long long *indpix,
//		double *S, std::string MixMatfile, std::vector<std::string> bolonames, std::string dirfile, std::string ellFile, /*string bextension,
//		string fextension,*/  std::string dir,
//		bool NORMLIN, bool NOFILLGAP, bool remove_polynomia, std::string noiseSppreffile,
//		std::string extentnoiseSp, std::string outdirSpN, std::string fits_filename);

void EstimPowerSpectra(struct user_options u_opt,struct detectors det,struct directories dir,struct input_commons com,
		long ns, long ff, long NAXIS1, long NAXIS2, long long npix, long iframe,
		long long *indpix, double *S, std::string MixMatfile, std::string ellFile, std::string extentnoiseSp, std::string fits_filename);

#endif /* ESTIMPS_H_ */
