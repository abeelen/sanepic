#ifndef MAP_MAKING_H_
#define MAP_MAKING_H_

//#include <cstdlib>
#include <string>
#include <vector>
#include "mpi_architecture_builder.h"


extern "C" {
#include <wcslib/wcs.h>
}



#define D2PI 6.2831853071795864769252867665590057683943387987502
/* dmod(A,B) - A modulo B (double) */
#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
	:(A)+(B)*floor(-(A)/(B))):(A))


//void computeMapMinima(std::vector<std::string> bolonames, std::string *fits_table,
//		long iframe_min, long iframe_max, long *nsamples,double pixdeg,
//		double &ra_min,double &ra_max,double &dec_min,double &dec_max);

int computeMapMinima(std::vector<detectors> det_vect, struct samples samples_struct,
		long iframe_min, long iframe_max,
		double &ra_min,double &ra_max,double &dec_min,double &dec_max);

int computeMapMinima_HIPE(std::vector<detectors> det_vect, struct samples samples_struct,
		long iframe_min, long iframe_max,
		double &ra_min,double &ra_max,double &dec_min,double &dec_max);

int minmax_flag(double  *& array, int *& flag, long size, double & min_array, double &  max_array);

void computeMapHeader(double pixdeg, char *ctype, char* prjcode, double * coordscorner,
		struct wcsprm *& wcs, long &NAXIS1, long &NAXIS2);

#endif /* MAP_MAKING_H_ */
