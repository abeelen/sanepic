/*
 * psdIO.cpp
 *
 *  Created on: 16 juin 2009
 *      Author: matthieu
 */

#include "psdIO.h"
#include "imageIO.h"

void write_psd_tofits(string fname, long nx, long ny,
		char dtype, void * psd1d) {

	fitsfile *fp;
	int fits_status = 0;

	long naxis = 2;           // number of dimensions
	long naxes[] = {nx, ny};  // size of dimensions
	long fpixel[] = {1, 1};   // index for write_pix
	long ndata = nx * ny;     // number of data points

	// create fits file
	if ( fits_create_file(&fp, fname.c_str(), &fits_status) )
		print_fits_error(fits_status);

	// create fits image (switch on data type)
	switch (dtype) {
	case 'd':    // double
		if ( fits_create_img(fp, DOUBLE_IMG, naxis, naxes, &fits_status) )
			print_fits_error(fits_status);
		if ( fits_write_pix(fp, TDOUBLE, fpixel, ndata, (double*) psd1d, &fits_status) )
			print_fits_error(fits_status);
		break;
	case 'l':    // long
		if ( fits_create_img(fp, LONG_IMG, naxis, naxes, &fits_status) )
			print_fits_error(fits_status);
		if ( fits_write_pix(fp, TLONG, fpixel, ndata, (long*) psd1d, &fits_status) )
			print_fits_error(fits_status);
		break;
	default:
		cerr << "write_fits: data type '" << dtype << "' not supported. Exiting.\n";
		exit(1);
	}

	// close file
	if (fits_close_file(fp, &fits_status))
		fits_report_error(stderr, fits_status);

}
