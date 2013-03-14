#ifdef HAVE_CONFIG_H
#include "../../config.h"
#else
#define PACKAGE_VERSION "Unknown"
#endif


#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "ErrorCode.h"
#include "CovMatrixIO.h"
#include "StructDefinition.h"
#include "InputFileIO.h"
#include "ErrorCode.h"
#include "Utilities.h"


extern "C" {
#include "nrutil.h"
#include <fitsio.h>
}

#include "sanePSIO.h"


using namespace std;

int restore_session(string tmp_dir, string filename, int &goto_step, double **commonm2, double **N,
		double **P, double **Rellth, double **Rellexp, double *SPref, long ndet, int ncomp, long nbins, long ns)
{

	cout << "Restoring : " << filename << endl;
	FILE* fp;
	string file = tmp_dir + "data_saved_sanePS_" + filename + ".bin";
	size_t len=0;

	if (NULL == (fp = fopen(file.c_str(), "rb")))
	{
		cout << "Unable to open " << file << " for reading\n";
		return 0;
	}

	len+=fread(&goto_step,sizeof(int),1,fp);

	switch(goto_step)
	{
	case 2:
		// load commonm2
		for (int ii = 0; ii < ncomp; ii++)
			for(long in=0; in<ns; in++)
				len+=fread(&commonm2[ii][in], sizeof(double), 1, fp);
		break;

	case 3:
		for (int idet = 0; idet < ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				//				fprintf(fp,"%10.15g\n",N[idet][ibin]);
				len+=fread(&N[idet][ibin], sizeof(double), 1, fp);

		for (int ii = 0; ii < ncomp; ii++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&P[ii][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&Rellth[idet][ibin], sizeof(double), 1, fp);
		break;

	case 4:
	case 5:
		// load N, P, Rellexp, Rellth, SPref
		for (int idet = 0; idet < ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&N[idet][ibin], sizeof(double), 1, fp);

		for (int ii = 0; ii < ncomp; ii++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&P[ii][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&Rellth[idet][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&Rellexp[idet][ibin], sizeof(double), 1, fp);

		len+=fread(SPref, sizeof(double), nbins, fp);
		break;

	case 6:
		return 0;

	default :
		cerr << "EE - completed_step has an incorrect value in : " << file << ". Exiting ...\n" << endl;
		return 1;

	}

	fclose(fp);

	return 0;
}


int save_session(string tmp_dir, string filename, int goto_step, double **commonm2, double **N,
		double **P, double **Rellexp, double **Rellth, double *SPref, long ndet, int ncomp, long nbins, long ns){

	FILE* fp;
	string file = tmp_dir + "data_saved_sanePS_" + filename + ".bin";
	size_t len=0;

	if (NULL == (fp = fopen(file.c_str(), "w+")))
	{
		cout << "Unable to open " << file << " for writing\n";
		return 1;
	}

	len+=fwrite(&goto_step,sizeof(int),1,fp);

	switch(goto_step)
	{
	case 2:
		// load commonm2
		for (int ii = 0; ii < ncomp; ii++)
			for(long in=0; in<ns; in++)
				len+=fwrite(&commonm2[ii][in], sizeof(double), 1, fp);
		break;

	case 3:
		// write N, P, Rellth
		for (int idet = 0; idet < ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&N[idet][ibin], sizeof(double), 1, fp);

		for (int ii = 0; ii < ncomp; ii++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&P[ii][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&Rellth[idet][ibin], sizeof(double), 1, fp);

		break;

	case 4:
	case 5:
		// write N, P, Rellexp, SPref
		for (int idet = 0; idet < ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&N[idet][ibin], sizeof(double), 1, fp);

		for (int ii = 0; ii < ncomp; ii++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&P[ii][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&Rellexp[idet][ibin], sizeof(double), 1, fp);

		len+=fwrite(SPref, sizeof(double), nbins, fp);
		break;

	default :
		cerr << " EE - completed_step has an incorrect value in : " << file << ". Exiting ...\n" << endl;
		return 1;

	}

	fclose(fp);

	return 0;
}


int write_to_disk(string outdir, string fits_filename, struct param_sanePS structPS, std::vector<std::string> det,  long nbins, double *ell, double **mixmat,
		double **Rellth, double **Rellexp, double **N, double *SPref, double **P)
{


	ostringstream temp_stream;
	string outFilename;
	string basename;

	long ndet = (long)det.size();

	double *data1d;

	basename = FitsBasename(fits_filename);

	temp_stream << "!" << outdir << basename << structPS.cov_matrix_suffix;
	outFilename= temp_stream.str();
	temp_stream.str("");

	if(write_CovMatrix(outFilename, det, nbins, ell, Rellth))
		return FILE_PROBLEM;


	temp_stream << "!" << outdir << basename << "_exp" << structPS.cov_matrix_suffix;
	outFilename= temp_stream.str();
	temp_stream.str("");

	if(write_CovMatrix(outFilename, det, nbins, ell, Rellexp))
		return FILE_PROBLEM;


#ifdef DEBUG
	string tempstr1, tempstr2;
	FILE *fp;

	temp_stream << outdir << basename << "_exp.psd";
	outFilename= temp_stream.str();
	temp_stream.str("");
	fp = fopen(outFilename.c_str(),"w");

	for (long idet1=0;idet1<ndet;idet1++){
		for (long idet2=0;idet2<ndet;idet2++){

			///// write power spectrum to disk
			tempstr1 = det[idet1];
			tempstr2 = det[idet2];
			fprintf(fp,"%s%s%s\n",tempstr1.c_str(),"-",tempstr2.c_str());
			fprintf(fp,"%d\n",(int)nbins);
			for (long ii=0;ii<nbins;ii++){
				fprintf(fp,"%g\t",ell[ii]);
				fprintf(fp,"%10.15g\n",(Rellexp[idet1*(ndet)+idet2][ii]+Rellexp[idet2*(ndet)+idet1][ii])/2.0*SPref[ii]);
			}
			fprintf(fp,"%g\n",ell[nbins]);
		}
	}
	fclose(fp);

#endif
	//
	//	temp_stream << outdirSpN  << basename << structPS.mix_suffix;
	//	testfile= temp_stream.str();
	//	temp_stream.str("");
	//
	//	fp = fopen(testfile.c_str(),"w");
	//	fprintf(fp,"%i\n",structPS.ncomp);
	//	for (long idet=0;idet<ndet;idet++){
	//		for (long iComp=0;iComp<structPS.ncomp;iComp++){
	//			fprintf(fp,"%10.15g ",mixmat[idet][iComp]);
	//		}
	//		fprintf(fp,"\n");
	//	}
	//	fclose(fp);

	temp_stream << "!" << outdir << basename << structPS.mix_suffix ;
	outFilename= temp_stream.str();
	temp_stream.str("");
	write_MixMatrix(outFilename, det, structPS.ncomp, mixmat);

	data1d = new double[structPS.ncomp*nbins];
	for (long i=0; i< structPS.ncomp; i++)
		for (long j=0; j<nbins; j++)
			data1d[i*nbins+j] = P[i][j]*SPref[j];

	temp_stream << outdir + "Nfinal_" << basename << ".fits";
	outFilename= temp_stream.str();
	temp_stream.str("");
	if (write_psd_tofits(outFilename,nbins,structPS.ncomp,'d',data1d, (char *) "CorrelatedNoise", false))
		cerr << "WW - Could not write correlated noise power spectra";
	delete [] data1d;

	data1d = new double[ndet*nbins];
	for (long i=0; i< ndet; i++)
		for (long j=0; j<nbins; j++)
			data1d[i*nbins+j] = N[i][j]*SPref[j];

	if (write_psd_tofits(outFilename, nbins, ndet, 'd', data1d, (char *) "UncorrelatedNoise", true))
		cerr << "WW - Could not write uncorrelated noise power spectra";
	delete [] data1d;

	write_psd_tofits(outFilename,    1, nbins+1,'d',    ell, (char *) "Frequency", true);


	return EXIT_SUCCESS;

}

int write_MixMatrix(string fname, std::vector<string> bolos, long ncomp, double **mixmat)
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
	// write the mixing matrix
	naxes[0] = ncomp;
	naxes[1] = ndet;
	if (fits_create_img(fptr, DOUBLE_IMG, 2, naxes, &status))
		return 1;

	// since mixmat is a NR matrix, one has to write it line by line :
	for (long idet = 0; idet < ndet; idet++) {
		fpixel[1] = idet + 1;
		if (fits_write_pix(fptr, TDOUBLE, fpixel, ncomp, mixmat[idet], &status))
			return 1;
	}
	if (fits_write_key(fptr, TSTRING, (char *) "EXTNAME",
			(char *) "MixingMatrix",
			(char *) "name of this binary table extension", &status))
		return 1;
	if (fits_write_comment(
			fptr,
			(char *) "This contains the mixing coefficient for the common modes",
			&status))
		return 1;
	if (fits_write_comment(
			fptr,
			(char *) "Each line contains a detector (NAXIS1) vs mode (NAXIS2)",
			&status))
		return 1;


	if (fits_write_chksum(fptr, &status)){
		cerr << "error checksum !\n";
		return 1;
	}

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


int read_MixMatrix(string fname, std::vector<string> &det_vect, long &ncomp, double **& mixmat)
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
		//		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// read the Channel List
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "channels", 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	fits_get_num_rows(fptr, &nBolos, &status);
	fits_get_colnum(fptr, CASEINSEN, (char*) "name", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);

	// Initialize the data container
	char ** data;
	data = new char*[nBolos];
	for (int i = 0; i < nBolos; i++) {
		data[i] = new char[repeat];
	}

	fits_read_col(fptr, TSTRING, colnum, 1, 1, nBolos, NULL, data, 0, &status);

	// convert to string vector and free the container
	det_vect.resize(nBolos);
	for (int i = 0; i < nBolos; i++)
		det_vect[i] = data[i];
	//		free(data[i]);

	for (int i = 0; i < nBolos; i++)
		delete [] data[i];
	//	free(data);
	delete [] data;

	// ---------------------------------------------
	// read the spectras
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "MixingMatrix", 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	fits_get_img_size(fptr, 2, naxes, &status);
	ncomp = naxes[0];
	if (naxes[1] != nBolos ){
		fits_report_error(stderr,213);
		return 1;
	}

	mixmat = dmatrix(0, nBolos - 1, 0, ncomp - 1);

	for (int idet = 0; idet < nBolos; idet++) {
		fpixel[1] = idet + 1;
		fits_read_pix(fptr, TDOUBLE, fpixel, ncomp, NULL, mixmat[idet], NULL, &status);
	}

	if (fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	return 0;
}

uint16_t assignMixMat(string fname, std::vector<string> det, long ncomp, double **& mixmat){

	double ** mixmatIn;
	long ndetOut, ndetIn, nCompIn;
	std::vector<string> detIn;
	std::vector<int> indexIn;

	if ( read_MixMatrix(fname, detIn, nCompIn, mixmatIn) ){
		cerr << "WW - Could not find input Mixing Matrix " << fname << endl;
		return 0;
	}

	ndetOut =   det.size();
	ndetIn  = detIn.size();

	// assign detectors
	indexIn.resize(ndetOut, -1);
	for (long idetOut = 0; idetOut < ndetOut; idetOut++) {
		for (long idetIn = 0; idetIn < ndetIn; idetIn++) {
			if (det[idetOut] == detIn[idetIn]){
				indexIn[idetOut] = idetIn;
				break;
			}
		}
	}

	// Check for missing data
	if ( *std::min_element( indexIn.begin(), indexIn.end()) == -1 ) {
		cerr << "EE - Input MixingMatrix must include all requested channels (missing : ";
		for (int idetOut = 0; idetOut < ndetOut; idetOut++) {
			if (indexIn[idetOut] == -1) {
				cerr << det[idetOut] << ", ";
			}
		}
		cerr << " )" << endl;
		return 1;
	}

	for (int idet = 0; idet < ndetOut; idet ++) {
		for (int iComp = 0; iComp < nCompIn; iComp++) {
			if (iComp < ncomp)
				mixmat[idet][iComp] = mixmatIn[indexIn[idet]][iComp];
		}
	}

	free_dmatrix(mixmatIn,0,ndetIn-1,0,nCompIn-1);

	return 0;
}

int write_psd_tofits(string fname, long nx, long ny, char dtype, void * psd1d, string ext_name, bool extend) {


	fitsfile *fp;
	int fits_status = 0;

	long naxis = 2;           // number of dimensions
	long naxes[] = {nx, ny};  // size of dimensions
	long fpixel[] = {1, 1};   // index for write_pix
	long ndata = nx * ny;     // number of data points

	fname = ( extend==false ? (std::string) "!": (std::string) "" ) + fname;

	if (extend) {
		if (fits_open_file(&fp, fname.c_str(), READWRITE, &fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}
	} else {
		// create fits file
		if (fits_create_file(&fp, fname.c_str(), &fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}

		// Dummy image
		long dummy_naxes[] = {0, 0};
		if (fits_create_img(fp, DOUBLE_IMG, 0, dummy_naxes, &fits_status)) {
			fits_report_error(stderr, fits_status);
			return 1;
		}

		// Add sanepic Version as a key
		if (fits_write_key(fp, TSTRING, (const char*) "SANEVER", (char *) PACKAGE_VERSION, (char *) "Sanepic Version", &fits_status)) {
			fits_report_error(stderr, fits_status); return 1;
		}

		// Add Comments to the header
		std::vector<std::string> comments;
		comments.push_back(" ");
		comments.push_back("This fits file was generated by SANEPIC version "+ (string) PACKAGE_VERSION);
		comments.push_back("For more informations about SANEPIC and for SANEPIC updates, please");
		comments.push_back("check our website at http://www.ias.u-psud.fr/sanepic ");
		comments.push_back(" ");

		for (u_int ii = 0; ii < comments.size(); ii++) {
			if (fits_write_comment(fp, (char *) comments[ii].c_str(), &fits_status))
				return 1;
		}


		if (fits_write_chksum(fp, &fits_status)) {
			cout << "error checksum !\n";
			return 1;
		}

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

	// write date to file
	if (fits_write_date(fp, &fits_status)) {
		fits_report_error(stderr, fits_status);
		return 1;
	}


	if (fits_update_key(fp, TSTRING, (char *) "EXTNAME", (void*) (ext_name.c_str()), (char *) "table name", &fits_status))
		return 1;


	// Add sanepic Version as a key
	if (fits_write_key(fp, TSTRING, (const char*) "SANEVER", (char *) PACKAGE_VERSION, (char *) "Sanepic Version", &fits_status)) {
		fits_report_error(stderr, fits_status); return 1;
	}

	// Add Comments to the header
	std::vector<std::string> comments;
	comments.push_back(" ");
	comments.push_back("This fits file was generated by SANEPIC version "+ (string) PACKAGE_VERSION);
	comments.push_back("For more informations about SANEPIC and for SANEPIC updates, please");
	comments.push_back("check our website at http://www.ias.u-psud.fr/sanepic ");
	comments.push_back(" ");

	for (u_int ii = 0; ii < comments.size(); ii++) {
		if (fits_write_comment(fp, (char *) comments[ii].c_str(), &fits_status))
			return 1;
	}

	if (fits_write_chksum(fp, &fits_status)) {
		cout << "error checksum !\n";
		return 1;
	}

	// close file
	if (fits_close_file(fp, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	return 0;

}

//int readMixmatTxt(string MixMatfile, long ndet, long ncomp, double **&mixmat) {
//	FILE *fp;
//	long ncomp2;
//	double dummy1; // used to read mixing matrix
//
//	if ((fp = fopen(MixMatfile.c_str(), "r")) == NULL) {
//		cerr << "ERROR: Can't find Mixing Matrix file. Exiting. \n";
//		cout
//		<< "Advice : verify the file is in your noise directory and that his name is : "
//		<< MixMatfile << endl;
//		return FILE_PROBLEM;
//	}
//	if (1 != fscanf(fp, "%ld", &ncomp2))
//		cerr << "EE - failed reading ncomp2 from file" << MixMatfile << endl;
//
//
//	for (long idet = 0; idet < ndet; idet++) {
//		for (long iComp = 0; iComp < ncomp2; iComp++) {
//			if (1 != fscanf(fp, "%lf", &dummy1))
//				cerr << "EE - failed reading element from file " << MixMatfile << endl;
//			if (iComp < ncomp)
//				mixmat[idet][iComp] = dummy1;
//		}
//	}
//	fclose(fp);
//	return 0;
//
//}
