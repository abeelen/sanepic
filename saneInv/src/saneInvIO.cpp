#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <vector>


extern "C" {
#include <nrutil.h>
#include <fitsio.h>
}


using namespace std;


void read_CovMatrix(string fname, std::vector<string> &bolos, long &nbins, double *&ell, double **&Rellth)
/*
 * This function read the NoiseNoise Matrices.
 */
{
	fitsfile *fptr;
	int status = 0;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	long nBolos, repeat, width;
	int colnum, typecode;


	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// read the Channel List
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "Channel List", NULL, &status))
		fits_report_error(stderr, status);



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
	for (int i = 0; i < nBolos; i++) {
		bolos[i] = data[i];
	}
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
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "Covariance Matrices", NULL, &status))
		fits_report_error(stderr, status);

	fits_get_img_size(fptr, 2, naxes, &status);
	if (naxes[1] != nBolos*nBolos || naxes[0] != nbins)
		fits_report_error(stderr,213);

	Rellth = dmatrix(0, nBolos * nBolos - 1, 0, nbins - 1);

	for (int i = 0; i < nBolos * nBolos; i++) {
		fpixel[1] = i + 1;
		fits_read_pix(fptr, TDOUBLE, fpixel, nbins, NULL, (Rellth)[i], NULL, &status);
	}

	if (fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

}


void write_InvNoisePowerSpectra(std::vector<string> bolos, long nbins, double * ell,
		double **Rellth, string outputDir, string suffix)
/*
 * This function writes the Inverse Covariance Matrices in binary format
 */
{

	string filename;
	FILE *fpw;
	long ndet = bolos.size();

	for (int idet = 0; idet < ndet; idet++) {

		// open file
		filename = outputDir + bolos[idet] + "_" + suffix;
		if ((fpw = fopen(filename.c_str(),"w")) == NULL){
			cerr << "ERROR: Can't write noise power spectra file" << filename << endl;
			exit(1);
		}

		// write sizes
		fwrite(&nbins, sizeof(long), 1, fpw);
		fwrite(&ndet, sizeof(long), 1, fpw);

		// write arrays
		fwrite(ell, sizeof(double), nbins + 1, fpw);
		fwrite(Rellth[idet], sizeof(double), ndet * nbins, fpw);

		// close file
		fclose(fpw);
	}

}
