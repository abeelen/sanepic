#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <vector>

#include "covMatrixIO.h"

extern "C" {
#include <nrutil.h>
#include <fitsio.h>
}

using namespace std;
/*
void read_noisefile(string fname, string bolo1bolo2, double *ell, double *SPN,
		long *nbins) {
	string line;

	//long ii;
	double dummy1, dummy2;

	ifstream Spfile(fname.c_str());
	if (!Spfile.is_open()) {
		cerr << "Error opening Noise power spectra file '" << fname
		<< "'Exiting.\n";
		exit(1);
	}

	while (!Spfile.eof()) {
		getline(Spfile, line);

		line.erase(0, line.find_first_not_of(" \t")); // remove leading white space
		if (line.empty() || line[0] == '#') // skip if empty or commented
			continue;
		line = line.substr(0, line.find_first_of(" \t")); // pick out first word

		if (line == bolo1bolo2) {
			getline(Spfile, line);
			*nbins = atoi(line.c_str());
			for (long ii = 0; ii < *nbins; ii++) {
				getline(Spfile, line);
				sscanf(line.c_str(), "%lf%lf", &dummy1, &dummy2);
				ell[ii] = dummy1;
				SPN[ii] = dummy2;
			}
			getline(Spfile, line);
			ell[*nbins] = atof(line.c_str());

		}

	}

	Spfile.close();

}
*/

void write_CovMatrix(string fname, std::vector<string> bolos, long nbins, double *ell, double **Rellth, double** mixmat, int ncomp)
/*
 * This function write the NoiseNoise Matrices in a fits file.
 */
{
	fitsfile *fptr;
	int status = 0;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	long nBolos = bolos.size();

	if (fits_create_file(&fptr, fname.c_str(), &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// write the Channel List

	char *ttype[] = { (char*) "NAME" };
	char *tform[] = { tableFormat(bolos) };
	char *tunit[] = { (char*) "NONE" };
	char **data;
	data = vString2carray(bolos);

	fits_create_tbl(fptr, BINARY_TBL, nBolos, 1, ttype, tform, tunit,
			(char*)"Channel List", &status);
	fits_write_col(fptr, TSTRING, 1, 1, 1, nBolos, data, &status);
	fits_write_key(fptr, TSTRING, (char *) "TUNIT1", (char *) "NONE",
			(char *) "physical unit of the field", &status);

	// ---------------------------------------------
	// write the Ells
	naxes[0] = nbins + 1;
	fits_create_img(fptr, FLOAT_IMG, 1, naxes, &status);
	fits_write_pix(fptr, TDOUBLE, fpixel, naxes[0], ell, &status);
	fits_write_key(fptr, TSTRING, (char *) "TUNIT1", (char *) "Hz",
			(char *) "physical unit of the field", &status);
	fits_write_key(fptr, TSTRING, (char *) "EXTNAME", (char *) "Frequency",
			(char *) "name of this binary table extension", &status);

	// ---------------------------------------------
	// write the spectras
	naxes[0] = nbins;
	naxes[1] = nBolos * nBolos;
	fits_create_img(fptr, DOUBLE_IMG, 2, naxes, &status);

	// since Rellth is a NR matrix, one has to write it line by line :
	for (long i = 0; i < nBolos * nBolos; i++) {
		fpixel[1] = i + 1;
		fits_write_pix(fptr, TDOUBLE, fpixel, nbins, Rellth[i], &status);
	}
	fits_write_key(fptr, TSTRING, (char *) "EXTNAME",
			(char *) "Covariance Matrices",
			(char *) "name of this binary table extension", &status);
	fits_write_comment(
			fptr,
			(char *) "This contains the Fourrier transform of the covariance matrices",
			&status);
	fits_write_comment(
			fptr,
			(char *) "Each line contains a couple of detector (NAXIS1) vs Frequency (NAXIS2)",
			&status);


	//-------------------------------------------------
	// write the mixing matrix
	naxes[0] = ncomp;
	naxes[1] = nBolos;
	fits_create_img(fptr, DOUBLE_IMG, 2, naxes, &status);

	// since mixmat is a NR matrix, one has to write it line by line :
	for (long i = 0; i < nBolos; i++) {
		fpixel[1] = i + 1;
		fits_write_pix(fptr, TDOUBLE, fpixel, ncomp, mixmat[i], &status);
	}
	fits_write_key(fptr, TSTRING, (char *) "EXTNAME",
			(char *) "mixing Matrices",
			(char *) "name of this binary table extension", &status);
	fits_write_comment(
			fptr,
			(char *) "This contains the mixing matrix",
			&status);
	fits_write_comment(
			fptr,
			(char *) "Each line contains a detector (NAXIS1) for each component (NAXIS2)",
			&status);

	if (fits_close_file(fptr, &status))
		fits_report_error(stderr, status);


}

void read_CovMatrix(string fname, std::vector<string> &bolos, long *nbins,
		double **ell, double ***Rellth,double ***mixmat, int *ncomp)
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
		free(data[i]);
	}
	free(data);

	// ---------------------------------------------
	// read the Ell
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "Frequency", NULL, &status))
		fits_report_error(stderr, status);

	fits_get_img_size(fptr, 1, naxes, &status);
	*nbins = naxes[0] ;//- 1; // TODO : verifier que c'est OK sans ce -1 !
	*ell = new double[naxes[0]];
	fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL, *ell, NULL, &status);

	// ---------------------------------------------
	// read the spectras
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "Covariance Matrices", NULL, &status))
		fits_report_error(stderr, status);

	fits_get_img_size(fptr, 2, naxes, &status);
	if (naxes[1] != nBolos*nBolos || naxes[0] != *nbins)
		fits_report_error(stderr,213);

	*Rellth = dmatrix(0, nBolos * nBolos - 1, 0, *nbins - 1);

	for (int i = 0; i < nBolos * nBolos; i++) {
		fpixel[1] = i + 1;
		fits_read_pix(fptr, TDOUBLE, fpixel, *nbins, NULL, (*Rellth)[i], NULL, &status);
	}

	//------------------------------------------------
	// read the mixmat
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mixing Matrices", NULL, &status))
		fits_report_error(stderr, status);

	fits_get_img_size(fptr, 2, naxes, &status);
	if (naxes[1] != nBolos)// || naxes[0] != *ncomp)
		fits_report_error(stderr,213);

	*ncomp=naxes[0];
	*mixmat = dmatrix(0, nBolos - 1, 0, *ncomp - 1);

	for (int i = 0; i < nBolos; i++) {
		fpixel[1] = i + 1;
		fits_read_pix(fptr, TDOUBLE, fpixel, *ncomp, NULL, (*mixmat)[i], NULL, &status);
	}

	*ncomp=(int)*ncomp;

	if (fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

}

char** vString2carray(std::vector<string> strings) {
	// Transform a vector of string into a array of char

	int stringLength = maxStringLength(strings);
	int nBolos = strings.size();

	char **data;

	data = new char*[nBolos];

	for (int i = 0; i < nBolos; i++) {
		data[i] = new char[stringLength];
		strcpy(data[i], strings[i].c_str());
	}

	return data;

}

char* tableFormat(std::vector<string> strings) {
	// Return the table Format for the given vector of string
	stringstream stream_tform;
	string string_tform;
	char* p_tform;

	stream_tform << maxStringLength(strings) << "A";
	string_tform = stream_tform.str();

	p_tform = new char[string_tform.size() + 1];
	strcpy(p_tform, string_tform.c_str());
	return p_tform;
}

long maxStringLength(std::vector<string> strings) {
	// Return the longest string length of a string vector
	unsigned long maxSize = 0;
	std::vector<string>::iterator itString;

	for (itString = strings.begin(); itString != strings.end(); itString++) {
		string iString = *(itString);
		if (iString.size() > maxSize)
			maxSize = iString.size();

	}
	return maxSize;
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
		filename = outputDir + bolos[idet] + "-all_Inv" + suffix;
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

void write_ReducedMixingMatrix(double **mixmat,long ndet,int ncomp, string outputDir)
// Writes the reduced mixing matrix in a binary file
{

	string filename; /*! Reduced mixing matrix internal filename (fixed by us, not modifiable by users)*/
	FILE *fp;

	// open file
	filename = outputDir + "Reduced_MixingMatrix_internal_data.bin"; //Reduced mixing matrix binary file
	if((fp=fopen(filename.c_str(),"w"))== NULL){
		cerr << "ERROR: Can't write Reduced MixingMatrix file" << filename << endl;
		exit(1);
	}

	//Read sizes
	fwrite(&ndet,sizeof(long),1,fp); // number of detectors in the mixmat
	fwrite(&ncomp,sizeof(int),1,fp); // number of noise component in the mixmat

	for (long idet=0;idet<ndet;idet++)
		for (int icomp=0;icomp<ncomp;icomp++)
			fwrite(&mixmat[idet][icomp],sizeof(double),1,fp); // writes the mixmat element by element
			// TODO : verify if it is faster to read line by line due to dmatrix allocation

	//close file
	fclose(fp);

	//------------------------------DEBUG MODE------------------------------------
	// open file
	filename = outputDir + "Reduced_MixingMatrix_internal_data_test.txt";
	if((fp=fopen(filename.c_str(),"w"))== NULL){
		cerr << "ERROR: Can't write Reduced MixingMatrix file" << filename << endl;
		exit(1);
	}

	//read sizes
	fprintf(fp,"%ld",ndet);
	fprintf(fp,"%d",ncomp);

	for (long idet=0;idet<ndet;idet++)
		for (int icomp=0;icomp<ncomp;icomp++)
			fprintf(fp,"%lf",mixmat[idet][icomp]);



	//close file
	fclose(fp);



}

void read_ReducedMixingMatrix(double **&mixmat,long &ndet,int &ncomp, string dir){

	string filename; /*! Reduced mixing matrix internal filename (fixed by us, not modifiable by users)*/
	FILE *fp;

	// open file
	filename = dir + "Reduced_MixingMatrix_internal_data.bin"; //Reduced mixing matrix binary file
	if((fp=fopen(filename.c_str(),"r"))== NULL){
		cerr << "ERROR: Can't find Reduced MixingMatrix file" << filename << endl;
		exit(1);
	}

	//Read sizes
	fread(&ndet,sizeof(long),1,fp); // number of detector in the mixmat
	fread(&ncomp,sizeof(int),1,fp); // number of noise component

	//allocate memory considering readed sizes
	mixmat=dmatrix(0, ndet - 1, 0, ncomp - 1);


	for (long idet=0;idet<ndet;idet++)
		for (int icomp=0;icomp<ncomp;icomp++)
			fread(&mixmat[idet][icomp],sizeof(double),1,fp); // reads mixmat element by element


	//close file
	fclose(fp);

}

// TODO : la fonction fait doublon avec read_noise_file dans inline_IO2.cpp, celle ci permet de read ndet en plus (depend du format d'ecriture)
void read_InvNoisePowerSpectra(string prefix, string boloName, string suffix,
		long * nbins, long * ndet, double ** ell, double *** SpN_all)
/*
 * This function reads the Inverse Covariance Matrices in binary format
 */
{

	string filename;
	FILE *fp;

	filename = prefix + boloName + "-all2" + suffix; // TODO : (reminder) remplacé -all2 par -all : mat 28_07
	cout << filename << endl;
	if ((fp = fopen(filename.c_str(), "r")) == NULL) {
		cerr << "ERROR: Can't read noise power spectra file" << filename
		<< endl;
		exit(1);
	}
	// Read sizes
	fread(nbins, sizeof(long), 1, fp);
	fread(ndet, sizeof(long), 1, fp);


	// Allocate memory
	*ell = new double[(*nbins) + 1];
	*SpN_all = dmatrix(0, (*ndet) - 1, 0, (*nbins) - 1);

	// Read arrays
	fread(*ell,     sizeof(double), (*nbins) + 1, fp);
	for (int i=0; i<(*ndet); i++)
		fread((*SpN_all)[i], sizeof(double), (*nbins), fp);

	for (int i=0; i< *nbins; i++)
		cout << (*SpN_all)[0][i] << " ";
	cout << endl;

	cout << "here final" << endl;

	fclose(fp);

}
