#ifndef ESTIMPS_H_
#define ESTIMPS_H_

#include <vector>
#include <string>
#include "mpi_architecture_builder.h"


int EstimPowerSpectra(struct param_process proc_param,struct detectors det,struct common dir, struct param_positions pos_param,
		long ns, long NAXIS1, long NAXIS2, long long npix, long iframe,
		long long *indpix, double *S, std::string MixMatfile, std::string ellFile, std::string fits_filename, long ncomp, double fcut, int rank);

#endif /* ESTIMPS_H_ */
