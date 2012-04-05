
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>


#include "struct_definition.h"
#include "temporary_IO.h"

extern "C" {
#include "getdata.h"
}

#include "coord.h"

using namespace std;


const double D2PI = 2*M_PI;
const double DEG2RAD = M_PI/180 ;

/* dsign(A,B) - magnitude of A with sign of B (double) */
#define dsign(A,B) ((B)<0.0?-(A):(A))

/* dmod(A,B) - A modulo B (double) */
#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
		:(A)+(B)*floor(-(A)/(B))):(A))


void equatorial2galactic (long nx, double *RA, double *DEC, double **glon, double **glat)
// Transformation from J200 alpha,delta to Galactic
{

	double inCart[3], outCart[3];
	double tmp1, tmp2;

	static double RotMat[3][3] = { { -0.054875539726, -0.873437108010, -0.483834985808 } , \
			{  0.494109453312, -0.444829589425,  0.746982251810 } , \
			{ -0.867666135858, -0.198076386122,  0.455983795705 } };



	for (long ii=0; ii< nx; ii++) {
		Spherical2Cartesian (RA[ii]*DEG2RAD, DEC[ii]*DEG2RAD, inCart);
		RotVect (RotMat, inCart, outCart);
		Cartesian2Spherical (outCart, &tmp1, &tmp2);
		tmp1 = ranrm (tmp1);
		tmp2 = range (tmp2);
		(*glon)[ii] = tmp1/DEG2RAD;
		(*glat)[ii] = tmp2/DEG2RAD;
	}

}

void galactic2equatorial (long nx, double *glon, double *glat, double **RA, double **DEC)
// Transformation from J200 alpha,delta to Galactic
{

	double inCart[3], outCart[3];
	double tmp1, tmp2;

	// http://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3.C3.973_matrices
/* //  Has to be inverted
 * 	static double RotMat[3][3] = { {-0.054875539726, -0.873437108010, -0.483834985808},  \
 *			{  0.494109453312, -0.444829589425, 0.746982251810}, \
 *			{-0.867666135858, -0.198076386122, 0.455983795705 } };
*/
	static double RotMat[3][3] = { {-0.054875539726, 0.494109453312, -0.867666135858},  \
			{ -0.873437108010 , -0.444829589425, -0.198076386122}, \
			{ -0.483834985808 ,  0.746982251810 , 0.455983795705 } };


	for (long ii=0; ii< nx; ii++) {
		Spherical2Cartesian (glon[ii]*DEG2RAD, glat[ii]*DEG2RAD, inCart);
		RotVect (RotMat, inCart, outCart);
		Cartesian2Spherical (outCart, &tmp1, &tmp2);
		tmp1 = ranrm (tmp1);
		tmp2 = range (tmp2);
		(*RA)[ii]  = tmp1/DEG2RAD;
		(*DEC)[ii] = tmp2/DEG2RAD;
	}

}


void Spherical2Cartesian (double RA, double DEC, double Cart[3])
// Spherical to cartesian
{
	// http://en.wikipedia.org/wiki/Spherical_coordinate_system
	// with theta in [0, pi[ and phi in [0, 2 pi[
	//   x = r sin (theta) * cos (phi)
	//   y = r sin (theta) * sin (phi)
	//   z = r cos (theta)
	// as DEC is in [-pi/2, pi/2] and RA in [-pi, pi ]
	// changing coordinates leave
	//  x  = r sin (DEC+pi/2) cos (RA + pi) = r cos(DEC) * cos(RA)
	//  y  = r sin (DEC+pi/2) sin (RA + pi) = r cos(DEC) * sin(RA)
	//  z  = r cos (DEC+pi/2)               = r sin(DEC)
	double cosDEC;
	cosDEC = cos (DEC);
	Cart[0] = cos (RA) * cosDEC;
	Cart[1] = sin (RA) * cosDEC;
	Cart[2] = sin (DEC);
}


void RotVect (double Mat[3][3], double IN[3], double OUT[3])
// Apply 3D rotation
// IN  rotation matrix,
//  q1 vector to be rotated
// OUT result
{
	//Hardcoded for speed....
	OUT[0] = Mat[0][0]*IN[0] + Mat[0][1]*IN[1] + Mat[0][2]*IN[2] ;
	OUT[1] = Mat[1][0]*IN[0] + Mat[1][1]*IN[1] + Mat[1][2]*IN[2] ;
	OUT[2] = Mat[2][0]*IN[0] + Mat[2][1]*IN[1] + Mat[2][2]*IN[2] ;

}

void Cartesian2Spherical (double Cart[3], double *RA, double *DEC)
// Cartesian to Spherical
{
	// http://en.wikipedia.org/wiki/Spherical_coordinate_system#Cartesian_coordinates
	// r = sqrt(x^2+y^2+z^2)
	// theta = acos(z/r)
	// phi   = atan2(y/x)
	double X, Y, Z, R;
	X = Cart[0];
	Y = Cart[1];
	Z = Cart[2];
	R = sqrt (X * X + Y * Y);
	*RA   = (R != 0.0) ? atan2 (Y, X) : 0.0;
	*DEC  = (Z != 0.0) ? atan2 (Z, R) : 0.0;

}

double ranrm (double angle)
// Put Angle in  Range 0 - 2 pi
{
	double tmp;
	tmp = dmod (angle, D2PI);
	return (tmp >= 0.0) ? tmp : tmp + D2PI;
}

double range (double angle)
// Put Angle in Range +-pi
{
	double tmp;
	tmp = dmod (angle, D2PI);
	return (fabs (tmp) < M_PI) ? tmp : tmp - dsign (D2PI, angle);
}



//using namespace std;
//int main(int argc, char *argv[])
//{
//
//  double *ra, *dec;
//  double *glon, *glat;
//  int nx = 2;
//
//
//  ra  = new double[nx];
//  dec = new double[nx];
//
//  glon = new double[nx];
//  glat = new double[nx];
//
//  ra[0]  = 24.5;
//  dec[0] = 29.4;
//
//  ra[1]  = 12;
//  dec[1] = -2.4;
//
//
//  cout << "Hello World" << endl;
//
//  equatorial2galactic(nx, ra,dec, &glon, &glat);
//
//  for (long ii=0; ii< nx; ii++){
//    cout << "(ra, dec)    = "  << setprecision(12) << ra[ii]   << " " << dec[ii]  << endl;
//    cout << "(glon, glat) = "  << setprecision(12) << glon[ii] << " " << glat[ii] << endl;
//  }
//
//}


// TODO : test this...
int convert_Dirfile_LON_LAT(struct samples samples_struct, struct param_sanePos pos_param)
{

	double *lon, *lat;
	double *new_lon, *new_lat;
	DIRFILE* D;

	long ns;
	int n_get;
	D = samples_struct.dirfile_pointer;


	for (long iframe=samples_struct.iframe_min;iframe<samples_struct.iframe_max;iframe++){

		ns = samples_struct.nsamples[iframe];

		lon = new double[ns];
		lat = new double[ns];

		new_lon = new double[ns];
		new_lat = new double[ns];

		for (unsigned long idet = 0; idet < samples_struct.bolo_list[iframe].size() ; idet++) {

			string lon_file = "LON_" + samples_struct.basevect[iframe] + "_" + samples_struct.bolo_list[iframe][idet];
			string lat_file = "LAT_" + samples_struct.basevect[iframe] + "_" + samples_struct.bolo_list[iframe][idet];

			const char * field_lon, *field_lat;

			field_lon = lon_file.c_str();
			field_lat = lat_file.c_str();

			// Read lon and lat

			n_get = gd_getdata(D, field_lon, 0, 0, 0, ns, GD_DOUBLE, lon);
			if (gd_error(D) != 0) {
				cout << "error getdata in Convert : read " << n_get << endl;
				return 1;
			}

			n_get = gd_getdata(D, field_lat, 0, 0, 0, ns, GD_DOUBLE, lat);
			if (gd_error(D) != 0) {
				cout << "error getdata in Convert : read " << n_get << endl;
				return 1;
			}

			//convert here

			if ( pos_param.eq2gal )
				equatorial2galactic(ns, lon, lat, &new_lon, &new_lat);

			if ( pos_param.gal2eq)
				galactic2equatorial(ns, lon, lat, &new_lon, &new_lat);



			// Write results
			n_get = gd_putdata(D, field_lon, 0, 0, 0, ns, GD_DOUBLE, new_lon);
			if (gd_error(D) != 0) {
				cout << "error putdata in Convert_Dirfile : wrote " << n_get << " and expected " << ns << endl;
				return 1;
			}
			n_get = gd_putdata(D, field_lat, 0, 0, 0, ns, GD_DOUBLE, new_lat);
			if (gd_error(D) != 0) {
				cout << "error putdata in write_lon : wrote " << n_get << " and expected " << ns << endl;
				return 1;
			}

			gd_flush(D, field_lon);
			gd_flush(D, field_lat);
		}


		delete [] lon;
		delete [] lat;
		delete [] new_lon;
		delete [] new_lat;
	}

	return 0;

}
