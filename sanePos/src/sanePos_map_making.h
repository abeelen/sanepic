#ifndef MAP_MAKING_H_
#define MAP_MAKING_H_

#include <cstdlib>

#define D2PI 6.2831853071795864769252867665590057683943387987502
/* dmod(A,B) - A modulo B (double) */
#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
	:(A)+(B)*floor(-(A)/(B))):(A))

double slaDranrm ( double angle );

void slaDs2tp ( double ra, double dec, double raz, double decz,
		double *xi, double *eta, int *j );

void slaDtp2s ( double xi, double eta, double raz, double decz,
		double *ra, double *dec );

void sph_coord_to_sqrmap(double pixdeg, double *ra, double *dec, double *phi,
		double *offsets, int ns, int *xx, int *yy, int *nn,
		double *coordscorner, double *tancoord, double *tanpix,
		bool fixcoord, double radius, double *offmap, double *radecsrc = NULL,bool compute_xx_yy=0);

/*void reproj_to_map( double *data, int *xx, int *yy, int ns, double **map,
		double **count, int nn, short *flag,
		double **map_f, double **count_f );*/


void flag_conditions(short *flag,/* double *scerr,*/ short *flpoint,
		long ns, long napod, int *xx, int *yy, int nn, double errarcsec,
		bool NOFILLGAP, unsigned char *rejectsamp);



#endif /* MAP_MAKING_H_ */
