#ifndef MAP_MAKING_H_
#define MAP_MAKING_H_

#include <cstdlib>
#include <string>
#include <vector>

extern "C" {
	#include <wcslib/wcs.h>
}

using namespace std;

#define D2PI 6.2831853071795864769252867665590057683943387987502
/* dmod(A,B) - A modulo B (double) */
#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
	:(A)+(B)*floor(-(A)/(B))):(A))


void computeMapMinima(std::vector<string> bolonames, string *fits_table,
		long iframe_min, long iframe_max, long *nsamples,double pixdeg,
		double &ra_min,double &ra_max,double &dec_min,double &dec_max);

void computeMapMinima_HIPE(std::vector<string> bolonames, string *fits_table,
		long iframe_min, long iframe_max, long *nsamples,double pixdeg,
		double &ra_min,double &ra_max,double &dec_min,double &dec_max);

void minmax_flag(double  *& array, short *& flag, long size, double & min_array, double &  max_array);

void computeMapHeader(double pixdeg, char *ctype, char* prjcode, double * coordscorner,
		struct wcsprm &wcs, long &NAXIS1, long &NAXIS2);

////TODO: slalib routine, will disapear
//double slaDranrm ( double angle );
//
//void slaDs2tp ( double ra, double dec, double raz, double decz,
//		double *xi, double *eta, int *j );
//
//void slaDtp2s ( double xi, double eta, double raz, double decz,
//		double *ra, double *dec );
//
//void sph_coord_to_sqrmap(double pixdeg, double *ra, double *dec, double *phi,
//		double *offsets, int ns, int *xx, int *yy, int *nn,
//		double *coordscorner, double *tancoord, double *tanpix,
//		bool fixcoord, double radius, /*double *offmap,*/ double *radecsrc = NULL,bool compute_xx_yy=0);
//
///*void reproj_to_map( double *data, int *xx, int *yy, int ns, double **map,
//		double **count, int nn, short *flag,
//		double **map_f, double **count_f );*/
//
//
//void flag_conditions(short *flag,/* double *scerr,*/ short *flpoint,
//		long ns, long napod, int *xx, int *yy, int nn, double errarcsec,
//		bool NOFILLGAP, unsigned char *rejectsamp);
//


#endif /* MAP_MAKING_H_ */
