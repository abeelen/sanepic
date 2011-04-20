#ifndef MAP_MAKING_H_
#define MAP_MAKING_H_

//#include <cstdlib>
#include <string>
#include <vector>
#include "mpi_architecture_builder.h"


extern "C" {
#include <wcslib/wcs.h>
}



int computeMapMinima(struct samples samples_struct,
		long iframe_min, long iframe_max,
		double &ra_min,double &ra_max,double &dec_min,double &dec_max, std::vector<std::vector<std::string> > bolo_vect);

int computeMapMinima_HIPE(std::string tmp_dir, struct samples samples_struct,
		long iframe_min, long iframe_max,
		double &ra_min,double &ra_max,double &dec_min,double &dec_max, std::vector<std::vector<std::string> > bolo_vect);

int minmax_flag(double  *& array, int *& flag, long size, double & min_array, double &  max_array);

void computeMapHeader(double pixdeg, char *ctype, char* prjcode, double * coordscorner,
		struct wcsprm *& wcs, long &NAXIS1, long &NAXIS2);

int do_PtNd_Naiv(struct samples samples_struct, double *PNd, std::string dir, std::vector<std::string> file, std::vector<std::string> det, long ndet, int orderpoly, int napod, double f_lppix, long ns,
		long long *indpix, long iframe, long *hits);

#endif /* MAP_MAKING_H_ */
