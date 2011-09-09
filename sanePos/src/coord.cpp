
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
		Spherical2Cartesian (RA[ii], DEC[ii], inCart);
		RotVect (RotMat, inCart, outCart);
		Cartesian2Spherical (outCart, &tmp1, &tmp2);
		tmp1 = ranrm (tmp1);
		tmp2 = range (tmp2);
		(*glon)[ii] = tmp1;
		(*glat)[ii] = tmp2;
	}

}

void galactic2equatorial (long nx, double *glon, double *glat, double **RA, double **DEC)
// Transformation from J200 alpha,delta to Galactic
{

	double inCart[3], outCart[3];
	double tmp1, tmp2;

	static double RotMat[3][3] = { {-0.054875539726, -0.873437108010, -0.483834985808},  \
			{  0.494109453312, -0.444829589425, 0.746982251810}, \
			{-0.867666135858, -0.198076386122, 0.455983795705 } };

	for (long ii=0; ii< nx; ii++) {
		Spherical2Cartesian (glon[ii], glat[ii], inCart);
		RotVect (RotMat, inCart, outCart);
		Cartesian2Spherical (outCart, &tmp1, &tmp2);
		tmp1 = ranrm (tmp1);
		tmp2 = range (tmp2);
		(*RA)[ii]  = tmp1;
		(*DEC)[ii] = tmp2;
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
//  ra[0]  = 24.5*DEG2RAD;
//  dec[0] = 29.4*DEG2RAD;
//
//  ra[1]  = 12*DEG2RAD;
//  dec[1] = -2.4*DEG2RAD;
//
//
//  cout << "Hello World" << endl;
//
//  equatorial2galactic(nx, ra,dec, &glon, &glat);
//
//  for (long ii=0; ii< nx; ii++){
//    cout << "(ra, dec)    = "  << setprecision(12) << ra[ii]/DEG2RAD   << " " << dec[ii]/DEG2RAD << endl;
//    cout << "(glon, glat) = "  << setprecision(12) << glon[ii]/DEG2RAD << " " << glat[ii]/DEG2RAD << endl;
//  }
//
//}


// TODO : fix this...

int convert_Dirfile_LON_LAT(struct param_common dir, struct samples samples_struct, struct param_sanePos pos_param, long iframe_min, long iframe_max, std::vector<std::vector<std::string> > bolo_vect)
{

	double *lon;
	double *lat;
	long ns;
	DIRFILE* D, *H;

	std::vector<string> det_vect;

	for (long iframe=iframe_min;iframe<iframe_max;iframe++){

		det_vect=bolo_vect[iframe];

		string base_name = samples_struct.basevect[iframe];
		string LONdir = dir.tmp_dir + "dirfile/" + base_name + "/LON";
		string LATdir = dir.tmp_dir + "dirfile/" + base_name + "/LAT";

		D = gd_open((char *) LONdir.c_str(), GD_RDWR | GD_CREAT |
				GD_VERBOSE | GD_UNENCODED); // | GD_TRUNC |
		H = gd_open((char *) LATdir.c_str(), GD_RDWR | GD_CREAT |
				GD_VERBOSE | GD_UNENCODED); // | GD_TRUNC |


		for(long idet=0; idet < (long)det_vect.size(); idet ++){

			string field = det_vect[idet];
			string lon_outfile = "LON_" + base_name + "_" + field;
			string lat_outfile = "LAT_" + base_name + "_" + field;

			// from dirfile...
			//			read_LON_LAT_from_fits(dir.data_dir + samples_struct.fitsvect[iframe], field, lon, lat, ns);

			//configure dirfile field
			gd_entry_t E;
			E.field = (char*) lon_outfile.c_str();
			E.field_type = GD_RAW_ENTRY;
			E.fragment_index = 0;
			E.spf = 1;
			E.data_type = GD_DOUBLE;
			E.scalar[0] = NULL;

			// add to the dirfile
			gd_add(D, &E);

			// write binary file on disk
			int n_write = gd_putdata(D, (char*) lon_outfile.c_str(), 0, 0, 0, ns, GD_DOUBLE, lon);
			if (gd_error(D) != 0) {
				cout << "error putdata in write_lon : wrote " << n_write
						<< " and expected " << ns << endl;
				return 1;
			}

			//configure dirfile field
			gd_entry_t F;
			F.field = (char*) lat_outfile.c_str();
			F.field_type = GD_RAW_ENTRY;
			F.fragment_index = 0;
			F.spf = 1;
			F.data_type = GD_DOUBLE;
			F.scalar[0] = NULL;

			// add to the dirfile
			gd_add(H, &F);

			// write binary file on disk
			n_write = gd_putdata(H, (char*) lat_outfile.c_str(), 0, 0, 0, ns, GD_DOUBLE, lat);
			if (gd_error(H) != 0) {
				cout << "error putdata in write_lat : wrote " << n_write
						<< " and expected " << ns << endl;
				return 1;
			}

			delete [] lon;
			delete [] lat;

			gd_flush(D,NULL);
			gd_flush(H,NULL);
		}

		// close dirfile
		if (gd_close(D)) {
			cout << "Dirfile gd_close error  " << LONdir
					<< endl;
			return 1;
		}

		// close dirfile
		if (gd_close(H)) {
			cout << "Dirfile gd_close error for : " << LATdir
					<< endl;
			return 1;
		}

	}


	return 0;
}
