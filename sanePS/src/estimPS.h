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


void EstimPowerSpectra(struct param_process proc_param,struct detectors det,struct directories dir, struct param_positions pos_param,
		long ns, long ff, long NAXIS1, long NAXIS2, long long npix, long iframe,
		long long *indpix, double *S, std::string MixMatfile, std::string ellFile, std::string extentnoiseSp, std::string fits_filename);

#endif /* ESTIMPS_H_ */
