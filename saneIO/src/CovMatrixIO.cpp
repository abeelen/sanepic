
#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include <gsl/gsl_matrix.h>


#include "CovMatrixIO.h"
#include "Utilities.h"


extern "C" {
#include <fitsio.h>
#include "nrutil.h"
#include "getdata.h"
}


using namespace std;

int write_CovMatrix(string fname, std::vector<string> bolos, long nbins, double *ell, double **Rellth)
{
	fitsfile *fptr;
	int status = 0;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	long ndet = bolos.size();

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


	if (fits_create_tbl(fptr, BINARY_TBL, ndet, 1, ttype, tform, tunit,
			(char*)"channels", &status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 1, 1, 1, ndet, data, &status))
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
	naxes[1] = ndet * ndet;
	if (fits_create_img(fptr, DOUBLE_IMG, 2, naxes, &status))
		return 1;

	// since Rellth is a NR matrix, one has to write it line by line :
	for (long i = 0; i < ndet * ndet; i++) {
		fpixel[1] = i + 1;
		if (fits_write_pix(fptr, TDOUBLE, fpixel, nbins, Rellth[i], &status))
			return 1;
	}
	if (fits_write_key(fptr, TSTRING, (char *) "EXTNAME",
			(char *) "PowerSpectra",
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


	for(long ii=0;ii<ndet;ii++)
		delete [] data[ii];
	delete [] data;
	delete [] *tform;

	return 0;
}


int read_CovMatrix(string fname, std::vector<string> &det_vect, long &nbins, double *&ell, gsl_matrix *& Rellth)
/*
 * This function read the NoiseNoise Matrices.
 */
{
	fitsfile *fptr;
	int status = 0;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	long ndet, repeat, width;
	int colnum, typecode;

	//	cout << fname << endl;

	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// read the Channel List
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "channels", 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	fits_get_num_rows(fptr, &ndet, &status);
	fits_get_colnum(fptr, CASEINSEN, (char*) "name", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);

	// Initialize the data container
	char ** data;
	data = new char*[ndet];
	for (int i = 0; i < ndet; i++) {
		data[i] = new char[repeat];
	}

	fits_read_col(fptr, TSTRING, colnum, 1, 1, ndet, NULL, data, 0, &status);

	// convert to string vector and free the container
	det_vect.resize(ndet);
	for (int i = 0; i < ndet; i++)
		det_vect[i] = data[i];
	//		free(data[i]);

	for (int i = 0; i < ndet; i++)
		delete [] data[i];
	//	free(data);
	delete [] data;

	// ---------------------------------------------
	// read the Ell
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "Frequency", 0, &status))
		fits_report_error(stderr, status);

	fits_get_img_size(fptr, 1, naxes, &status);
	nbins = naxes[0] - 1;
	ell = new double[naxes[0]];
	fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL, ell, NULL, &status);

	// ---------------------------------------------
	// read the spectras
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "PowerSpectra", 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	fits_get_img_size(fptr, 2, naxes, &status);
	if (naxes[1] != ndet*ndet || naxes[0] != nbins){
		fits_report_error(stderr,213);
		return 1;
	}

	Rellth = gsl_matrix_alloc(ndet*ndet, nbins);

	for (int i = 0; i < ndet * ndet; i++) {
		double * matrix_ptr = gsl_matrix_ptr(Rellth, i, 0);
		fpixel[1] = i + 1;
		fits_read_pix(fptr, TDOUBLE, fpixel, nbins, NULL, matrix_ptr, NULL, &status);
	}

	if (fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	return 0;
}

//saneInv
int write_InvNoisePowerSpectra(DIRFILE* D, std::vector<string> det_vect, string scan_name, long nbins, double * ell, gsl_matrix *Rellth)
/*
 * This function writes the Inverse Covariance Matrices in binary format
 */
{
	long ndet = det_vect.size();

	double * vector_ptr;
	string outfile;

	// ell binary filename
	outfile = "Ell_" + scan_name;

	// write binary file on disk
	int n_write = gd_putdata(D, (char*)outfile.c_str(), 0, 0, 0, nbins+1, GD_DOUBLE, ell);
	//		cout << "n_write : " << n_write << endl;
	if(gd_error(D)!=0){
		cout << "error gd_putdata : wrote " << n_write << " and expected " << nbins+1 << " for " << outfile <<  endl;
		return 1;
	}

	// flush dirfile
	gd_flush(D,(char*)outfile.c_str());


	for (int idet = 0; idet < ndet; idet++) {

		// spectra filename
		outfile = "InvNoisePS_" + scan_name + "_" + det_vect[idet];

		vector_ptr = gsl_matrix_ptr(Rellth, idet, 0);

		// write binary file on disk
		n_write = gd_putdata(D, (char*)outfile.c_str(), 0, 0, 0, ndet*nbins, GD_DOUBLE, vector_ptr);
		if(gd_error(D)!=0){
			cout << "error gd_putdata : wrote " << n_write << " and expected " << nbins * ndet << " for " << outfile << endl;
			return 1;
		}

		// flush dirfile
		gd_flush(D,(char*)outfile.c_str());

	}


	return 0;
}

uint32_t read_Ell(DIRFILE* D, string scan_name, long nbins, double * ell){

	//binary name
	string outfile;
	int nget;
	outfile = "Ell_" + scan_name;

	// fill ell with binary
	nget = gd_getdata(D, (char*)outfile.c_str(), 0, 0, 0, nbins+1, GD_DOUBLE, ell);
	if(gd_error(D)!=0 or nget != (nbins+1)){
		cout << "error getdata in read_InvNoisePowerSpectra : reading " << outfile << endl;
		return 1;
	}

	return 0;

}


//TODO: Separate into Read ONE Ell and SpN_all
uint32_t read_InvNoisePowerSpectra(DIRFILE* D, string scan_name, string boloName, long nbins, long ndet, double ** SpN_all)
/*
 * This function reads the Inverse Covariance Matrices in binary format
 */
{
	//binary name
	string outfile;
	int nget;

	// temp 1D array
	double *Rellth_full;

	// Power spectra binary file name
	outfile = "InvNoisePS_" + scan_name + "_" + boloName;

	//alloc temp 1D array
	Rellth_full = new double[nbins*ndet];

	// read whole 1 D array
	nget = gd_getdata(D, (char*)outfile.c_str(), 0, 0, 0, nbins*ndet, GD_DOUBLE, Rellth_full);
	if(gd_error(D)!=0 or nget != (nbins*ndet)){
		cout << "error getdata in read_InvNoisePowerSpectra : reading " << outfile << endl;
		return 1;
	}

	gd_flush(D, (char *) outfile.c_str());

	// reorganize as a 2D array
	for (long idet=0; idet<(ndet); idet++)
		for (long ibin=0; ibin<(nbins); ibin++)
			SpN_all[idet][ibin] = Rellth_full[idet*(nbins) + ibin];

	// clear temp 1D array
	delete [] Rellth_full;

	return 0;

}
