#include <iostream>
#include <string>


// OK

#include "imageIO.h"

using namespace std;

extern "C" {
#include <fitsio.h>
}


void print_fits_error(int status){
  if(status){
    fits_report_error(stderr, status); /* print error report */
    exit(status);    /* terminate the program, returning error status */
  }
  return;
}



void write_fits(string fname, double pixsize, long nx, long ny,
		double *tancoord, double *tanpix, int coordsyst, char dtype, void *data)
{
  // all angles in degrees
  // coordcenter is a 2-element array containing RA/DEC (or l/b) of the central pixel

  fitsfile *fp;
  int fits_status = 0;

  long naxis = 2;           // number of dimensions
  long naxes[] = {nx, ny};  // size of dimensions
  long fpixel[] = {1, 1};   // index for write_pix
  long ndata = nx * ny;     // number of data points

  double dtmp;
  char *strx, *stry;


  // create fits file
  if ( fits_create_file(&fp, fname.c_str(), &fits_status) )
    print_fits_error(fits_status);

  // create fits image (switch on data type)
  switch (dtype) {
  case 'd':    // double
    if ( fits_create_img(fp, DOUBLE_IMG, naxis, naxes, &fits_status) )
      print_fits_error(fits_status);
    break;
  case 'l':    // long
    if ( fits_create_img(fp, LONG_IMG, naxis, naxes, &fits_status) )
      print_fits_error(fits_status);
    break;
  default:
	  printf("write_fits: data type %c not supported. Exiting.\n",dtype);
    exit(1);
  }

  // write date to file
  if ( fits_write_date(fp, &fits_status) )
    print_fits_error(fits_status);

  // write map parameters (keywords)
  if ( fits_write_key(fp, TLONG, (char*)"NROW", &nx, (char*)"Number of rows", &fits_status) )
    print_fits_error(fits_status);

  if ( fits_write_key(fp, TLONG, (char*)"NCOL", &ny, (char*)"Number of columns", &fits_status) )
    print_fits_error(fits_status);

  if ( fits_write_key(fp, TDOUBLE, (char*)"PIXSIZE", &pixsize, (char*)"Size of pixels (deg)", &fits_status) )
    print_fits_error(fits_status);

  if ( fits_write_comment(fp, "Galactic coordinates",  &fits_status) )
    print_fits_error(fits_status);

  dtmp = (tanpix[0]); // 0-based index to 1
  if ( fits_write_key(fp, TDOUBLE, (char*)"CRPIX1", &dtmp, (char*)"X PIXEL OF TANGENT POINT", &fits_status) )
    print_fits_error(fits_status);

  dtmp = (tanpix[1]); // 0-based index to 1
  if ( fits_write_key(fp, TDOUBLE, (char*)"CRPIX2", &dtmp, (char*)"Y PIXEL OF TANGENT POINT", &fits_status) )
    print_fits_error(fits_status);

  dtmp = -pixsize;
  if ( fits_write_key(fp, TDOUBLE, (char*)"CDELT1", &dtmp, (char*)"COORD VALUE INCR DEG/PIXEL AT ORIGIN ON LINE AXIS",
		      &fits_status) )
    print_fits_error(fits_status);

  if ( fits_write_key(fp, TDOUBLE, (char*)"CDELT2", &pixsize, (char*)"COORD VALUE INCR DEG/PIXEL AT ORIGIN ON LINE AXIS",
		      &fits_status) )
    print_fits_error(fits_status);

  if (coordsyst == 2){
    if ( fits_write_key(fp, TDOUBLE, (char*)"CRVAL1", tancoord, (char*)"GLON AT TANGENT POINT (DEG)", &fits_status) )
      print_fits_error(fits_status);

    if ( fits_write_key(fp, TDOUBLE, (char*)"CRVAL2", tancoord+1, (char*)"GLAT AT TANGENT POINT (DEG)", &fits_status) )
      print_fits_error(fits_status);

    strx = (char *)"GLON-TAN";
    stry = (char *)"GLAT-TAN";

  } else {
    if ( fits_write_key(fp, TDOUBLE, (char*)"CRVAL1", tancoord, (char*)"RA AT TANGENT POINT (DEG)", &fits_status) )
      print_fits_error(fits_status);

    if ( fits_write_key(fp, TDOUBLE, (char*)"CRVAL2", tancoord+1, (char*)"DEC AT TANGENT POINT (DEG)", &fits_status) )
      print_fits_error(fits_status);

    strx = (char *)"RA---TAN";
    stry = (char *)"DEC--TAN";

  }

  if ( fits_write_key(fp, TSTRING, (char*)"CTYPE1", strx, (char*)"TANGENT PLANE PROJECTION", &fits_status) )
    print_fits_error(fits_status);

  if ( fits_write_key(fp, TSTRING, (char*)"CTYPE2", stry, (char*)"TANGENT PLANE PROJECTION", &fits_status) )
    print_fits_error(fits_status);


  // write map data
  switch (dtype) {
  case 'd':    // double
    if ( fits_write_pix(fp, TDOUBLE, fpixel, ndata, (double*) data, &fits_status) )
      print_fits_error(fits_status);
    break;
  case 'l':    // long
    if ( fits_write_pix(fp, TLONG, fpixel, ndata, (long*) data, &fits_status) )
      print_fits_error(fits_status);
    break;
  }

  // close file
  if(fits_close_file(fp, &fits_status))
    print_fits_error(fits_status);

}

