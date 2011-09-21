#ifndef ESTIMPS_H_
#define ESTIMPS_H_

#include <vector>
#include <string>
#include "mpi_architecture_builder.h"

int EstimPowerSpectra(std::vector<std::string> det, struct param_saneProc proc_param, struct param_common dir, struct param_sanePos pos_param, struct param_sanePS structPS, struct samples samples_struct,
		long NAXIS1, long NAXIS2, long long npix, long iframe, long long *indpix, double *S, int rank);

#endif /* ESTIMPS_H_ */
