#ifndef SANEPOS_PREPROCESS_H_
#define SANEPOS_PREPROCESS_H_

#include <string>
#include <vector>
#include "mpi_architecture_builder.h"


extern "C" {
	#include <wcslib/wcs.h>
}


int computePixelIndex(std::string outdir, std::vector<detectors> det_vect,
		struct samples samples_struct, struct param_process proc_param, struct param_positions pos_param, long iframe_min, long iframe_max,
		struct wcsprm * wcs, long NAXIS1, long NAXIS2, short *&mask,
		int factdupl,long long addnpix, long long *&pixon, int rank,
		long long *indpsrc, long long npixsrc, int &flagon, bool &pixout);

int computePixelIndex_HIPE(std::string outdir, std::vector<detectors> det_vect,
		struct samples samples_struct, struct param_process proc_param, struct param_positions pos_param, long iframe_min, long iframe_max,
		struct wcsprm * wcs, long NAXIS1, long NAXIS2, short *&mask,
		int factdupl,long long addnpix, long long *&pixon, int rank,
		long long *indpsrc, long long npixsrc, int &flagon, bool &pixout);

#endif /* SANEPOS_PREPROCESS_H_ */
