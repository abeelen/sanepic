#ifndef SANEPOS_PREPROCESS_H_
#define SANEPOS_PREPROCESS_H_

#include <string>
#include <vector>
#include "mpi_architecture_builder.h"


extern "C" {
	#include <wcslib/wcs.h>
}


int modify_mask_flag_in_dirfile(std::string tmp_dir, struct samples samples_struct, long long *indpsrc,
		long NAXIS1, long NAXIS2, long iframe_min, long iframe_max);

int computePixelIndex(std::string outdir,
		struct samples samples_struct, struct param_sanePre proc_param, struct param_sanePos pos_param, long iframe_min, long iframe_max,
		struct wcsprm * wcs, long NAXIS1, long NAXIS2, short *&mask,
		int factdupl,long long addnpix, long long *&pixon, int rank,
		long long *indpsrc, long long npixsrc, int &flagon, bool &pixout, std::vector<std::vector<std::string> > bolo_vect);

int computePixelIndex_HIPE(std::string outdir,
		struct samples samples_struct, struct param_sanePre proc_param, struct param_sanePos pos_param, long iframe_min, long iframe_max,
		struct wcsprm * wcs, long NAXIS1, long NAXIS2, short *&mask,
		int factdupl,long long addnpix, long long *&pixon, int rank,
		long long *indpsrc, long long npixsrc, int &flagon, bool &pixout, std::vector<std::vector<std::string> > bolo_vect);

#endif /* SANEPOS_PREPROCESS_H_ */
