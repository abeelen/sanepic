#ifndef CORR_PREPROCESS_H_
#define CORR_PREPROCESS_H_

#include <string>
#include <fstream>
#include <vector>
#include "struct_definition.h"
#include <fftw3.h>

extern "C" {
#include "nrutil.h"
}

int write_ftrProcesdata(double *S, struct param_sanePre proc_param, struct samples samples_struct, struct param_sanePos pos_param,
		std::string tmp_dir, std::vector<std::string> det, long ndet, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2,
		long long npix,	long long npixsrc, long long addnpix, double f_lppix, long ns, long iframe, int para_bolo_indice, int para_bolo_size, std::string fname);

int do_PtNd(struct samples samples_struct, double *PNd, std::string dir, std::string prefixe,
		std::vector<std::string> det, long ndet, double f_lppix, double fsamp, long ns, int para_bolo_indice, int para_bolo_size,
		long long *indpix, long NAXIS1, long NAXIS2, long long npix, long iframe, std::string filename,
		double *Mp, long *hits,std::string fname);

int write_tfAS(struct samples samples_struct, double *S, std::vector<std::string> det, long ndet, long long *indpix, long NAXIS1, long NAXIS2, long long npix,
		bool flgdupl, std::string dir, long ns, std::string filename, int para_bolo_indice, int para_bolo_size);

#endif /* CORR_PREPROCESS_H_ */
