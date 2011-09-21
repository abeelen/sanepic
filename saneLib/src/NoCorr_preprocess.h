#ifndef NOCORR_PREPROCESS_H_
#define NOCORR_PREPROCESS_H_

#include <vector>
#include <string>
#include "mpi_architecture_builder.h"


int do_PtNd_nocorr(double *PNd,std::string tmp_dir, struct param_saneProc proc_param, struct param_sanePos pos_param,
		struct samples samples_struct, std::vector<std::string> det, long ndet, double f_lppix, double f_lppix_Nk,
		long addnpix, long ns, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2, long long npix,
		long long npixsrc, long iframe, double *S, int para_bolo_indice, int para_bolo_size);


void do_PtNPS_nocorr(struct samples samples_struct, double *S, std::vector<std::string> noisevect, struct param_common dir,
		std::vector<std::string> det, long ndet, double f_lppix,double fsamp, bool flgdupl, long ns,
		long long *indpix, long NAXIS1, long NAXIS2, long long npix,
		long iframe, std::string fname, double *PtNPmatS, double *Mp, long *hits, int para_bolo_indice, int para_bolo_size);

#endif /* NOCORR_PREPROCESS_H_ */
