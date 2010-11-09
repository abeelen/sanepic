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
		double &ra_min,double &ra_max,double &dec_min,double &dec_max);

int computeMapMinima_HIPE(struct samples samples_struct,
		long iframe_min, long iframe_max,
		double &ra_min,double &ra_max,double &dec_min,double &dec_max);

int minmax_flag(double  *& array, int *& flag, long size, double & min_array, double &  max_array);

void computeMapHeader(double pixdeg, char *ctype, char* prjcode, double * coordscorner,
		struct wcsprm *& wcs, long &NAXIS1, long &NAXIS2);

#endif /* MAP_MAKING_H_ */
