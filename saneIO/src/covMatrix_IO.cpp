#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <algorithm>


#include "covMatrix_IO.h"
#include "utilities.h"

extern "C" {
#include <nrutil.h>
#include <fitsio.h>
#include "getdata.h"
}


using namespace std;


////sanePS
int write_CovMatrix(string fname, std::vector<string> bolos, long nbins, double *ell, double **Rellth)
/*
 * This function write the NoiseNoise Matrices in a fits file.
 */
{
	fitsfile *fptr;
	int status = 0;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	long nBolos = bolos.size();

	if (fits_create_file(&fptr, fname.c_str(), &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// write the Channel List

	char *ttype[] = { (char*) "NAME" };
	char *tform[] = { tableFormat(bolos) };
	char *tunit[] = { (char*) "None" };
	char **data;
	data = vString2carray(bolos);


	if (fits_create_tbl(fptr, BINARY_TBL, nBolos, 1, ttype, tform, tunit,
			(char*)"Channel List", &status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 1, 1, 1, nBolos, data, &status))
		return 1;
	if (fits_write_key(fptr, TSTRING, (char *) "TUNIT1", (char *) "NONE",
			(char *) "physical unit of the field", &status))
		return 1;

	// ---------------------------------------------
	// write the Ells
	naxes[0] = nbins + 1;
	if (fits_create_img(fptr, FLOAT_IMG, 1, naxes, &status))
		return 1;
	if (fits_write_pix(fptr, TDOUBLE, fpixel, naxes[0], ell, &status))
		return 1;
	if (fits_write_key(fptr, TSTRING, (char *) "TUNIT1", (char *) "Hz",
			(char *) "physical unit of the field", &status))
		return 1;
	if (fits_write_key(fptr, TSTRING, (char *) "EXTNAME", (char *) "Frequency",
			(char *) "name of this binary table extension", &status))
		return 1;

	if (fits_write_chksum(fptr, &status)){
		cout << "error checksum !\n";
		return 1;
	}

	// ---------------------------------------------
	// write the spectras
	naxes[0] = nbins;
	naxes[1] = nBolos * nBolos;
	if (fits_create_img(fptr, DOUBLE_IMG, 2, naxes, &status))
		return 1;

	// since Rellth is a NR matrix, one has to write it line by line :
	for (long i = 0; i < nBolos * nBolos; i++) {
		fpixel[1] = i + 1;
		if (fits_write_pix(fptr, TDOUBLE, fpixel, nbins, Rellth[i], &status))
			return 1;
	}
	if (fits_write_key(fptr, TSTRING, (char *) "EXTNAME",
			(char *) "Covariance Matrices",
			(char *) "name of this binary table extension", &status))
		return 1;
	if (fits_write_comment(
			fptr,
			(char *) "This contains the Fourrier transform of the covariance matrices",
			&status))
		return 1;
	if (fits_write_comment(
			fptr,
			(char *) "Each line contains a couple of detector (NAXIS1) vs Frequency (NAXIS2)",
			&status))
		return 1;

	if (fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}


	for(long ii=0;ii<nBolos;ii++)
		delete [] data[ii];
	delete [] data;
	delete [] *tform;

	return 0;
}


//saneInv
int read_CovMatrix(string fname, std::vector<string> &bolos, long &nbins, double *&ell, double **&Rellth)
/*
 * This function read the NoiseNoise Matrices.
 */
{
	fitsfile *fptr;
	int status = 0;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	long nBolos, repeat, width;
	int colnum, typecode;

	//	cout << fname << endl;

	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// read the Channel List
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "Channel List", NULL, &status)){
		fits_report_error(stderr, status);
		return 1;
	}



	fits_get_num_rows(fptr, &nBolos, &status);
	fits_get_colnum(fptr, CASEINSEN, (char*) "NAME", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);

	// Initialize the data container
	char ** data;
	data = new char*[nBolos];
	for (int i = 0; i < nBolos; i++) {
		data[i] = new char[repeat];
	}

	fits_read_col(fptr, TSTRING, colnum, 1, 1, nBolos, NULL, data, 0, &status);

	// convert to string vector and free the container
	bolos.resize(nBolos);
	for (int i = 0; i < nBolos; i++)
		bolos[i] = data[i];
	//		free(data[i]);

	for (int i = 0; i < nBolos; i++)
		delete [] data[i];
	//	free(data);
	delete [] data;

	// ---------------------------------------------
	// read the Ell
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "Frequency", NULL, &status))
		fits_report_error(stderr, status);

	fits_get_img_size(fptr, 1, naxes, &status);
	nbins = naxes[0] - 1;
	ell = new double[naxes[0]];
	fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL, ell, NULL, &status);

	// ---------------------------------------------
	// read the spectras
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "Covariance Matrices", NULL, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	fits_get_img_size(fptr, 2, naxes, &status);
	if (naxes[1] != nBolos*nBolos || naxes[0] != nbins){
		fits_report_error(stderr,213);
		return 1;
	}

	Rellth = dmatrix(0, nBolos * nBolos - 1, 0, nbins - 1);

	for (int i = 0; i < nBolos * nBolos; i++) {
		fpixel[1] = i + 1;
		fits_read_pix(fptr, TDOUBLE, fpixel, nbins, NULL, (Rellth)[i], NULL, &status);
	}

	if (fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	return 0;
}

//saneInv
int write_InvNoisePowerSpectra(DIRFILE* D, std::vector<string> bolos, long nbins, double * ell,
		double **Rellth, string suffix)
/*
 * This function writes the Inverse Covariance Matrices in binary format
 */
{
	long ndet = bolos.size();

	for (int idet = 0; idet < ndet; idet++) {

		// ell binary filename
		string outfile = bolos[idet] + "_" + suffix + "_ell"; // suffix is now made with fitsbasename of the scan file + noise suffix from inifile !

		// write binary file on disk
		int n_write = gd_putdata(D, (char*)outfile.c_str(), 0, 0, 0, nbins+1, GD_DOUBLE, ell);
		//		cout << "n_write : " << n_write << endl;
		if(gd_error(D)!=0){
			cout << "error gd_putdata : wrote " << n_write << " and expected " << nbins+1 << " for " << outfile <<  endl;
			return 1;
		}

		// spectra filename
		outfile = bolos[idet] + "_" + suffix;

		// write binary file on disk
		n_write = gd_putdata(D, (char*)outfile.c_str(), 0, 0, 0, ndet*nbins, GD_DOUBLE, Rellth[idet]);
		if(gd_error(D)!=0){
			cout << "error gd_putdata : wrote " << n_write << " and expected " << nbins * ndet << " for " << outfile << endl;
			return 1;
		}

	}

	// flush dirfile
	gd_flush(D,NULL);

	return 0;
}

int read_InvNoisePowerSpectra(DIRFILE* D, string outputDir, string boloName, string suffix,
		long nbins, long ndet, double ** ell, double *** SpN_all)
/*
 * This function reads the Inverse Covariance Matrices in binary format
 */
{
	//binary name
	string outfile = boloName + "_" + suffix + "_ell";

	// temp 1D array
	double *Rellth_full;

	// alloc ell
	*ell=new double[nbins+1];

	// fill ell with binary
	int nget = gd_getdata(D, (char*)outfile.c_str(), 0, 0, 0, nbins+1, GD_DOUBLE, *ell);
	if(gd_error(D)!=0){
		cout << "error getdata in read_InvNoisePowerSpectra : reading " << outfile << endl;
		return 1;
	}

	// Power spectra binary file name
	outfile = boloName + "_" + suffix;

	//alloc temp 1D array
	Rellth_full = new double[nbins*ndet];

	// read whole 1 D array
	nget = gd_getdata(D, (char*)outfile.c_str(), 0, 0, 0, nbins*ndet, GD_DOUBLE, Rellth_full);
	if(gd_error(D)!=0){
		cout << "error getdata in read_InvNoisePowerSpectra : reading " << outfile << endl;
		return 1;
	}

	// alloc spectra 2D array
	*SpN_all = dmatrix(0, (ndet) - 1, 0, (nbins) - 1);

	// reorganize as a 2D array
	for (long i=0; i<(ndet); i++)
		for (long ibin=0; ibin<(nbins); ibin++)
			(*SpN_all)[i][ibin] = Rellth_full[i*(nbins) + ibin];

	// clear temp 1D array
	delete [] Rellth_full;

	return 0;

}

int write_psd_tofits(string fname, long nx, long ny,
		char dtype, void * psd1d) {

	fitsfile *fp;
	int fits_status = 0;

	long naxis = 2;           // number of dimensions
	long naxes[] = {nx, ny};  // size of dimensions
	long fpixel[] = {1, 1};   // index for write_pix
	long ndata = nx * ny;     // number of data points

	// create fits file
	if ( fits_create_file(&fp, fname.c_str(), &fits_status) ){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	// create fits image (switch on data type)
	switch (dtype) {
	case 'd':    // double
		if ( fits_create_img(fp, DOUBLE_IMG, naxis, naxes, &fits_status) )
		{
			fits_report_error(stderr, fits_status);
			return 1;
		}
		if ( fits_write_pix(fp, TDOUBLE, fpixel, ndata, (double*) psd1d, &fits_status) ){
			fits_report_error(stderr, fits_status);
			return 1;
		}
		break;
	case 'l':    // long
		if ( fits_create_img(fp, LONG_IMG, naxis, naxes, &fits_status) ){
			fits_report_error(stderr, fits_status);
			return 1;
		}
		if ( fits_write_pix(fp, TLONG, fpixel, ndata, (long*) psd1d, &fits_status) ){
			fits_report_error(stderr, fits_status);
			return 1;
		}
		break;
	default:
		cerr << "write_fits: data type '" << dtype << "' not supported. Exiting.\n";
		return 1;
	}

	// close file
	if (fits_close_file(fp, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	return 0;

}


