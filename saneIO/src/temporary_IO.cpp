
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "struct_definition.h"
#include "inputFileIO.h"
#include "temporary_IO.h"
#include "dataIO.h"


extern "C" {
#include "getdata.h"
#include <fitsio.h>
}

using namespace std;

//TODO: Absolutely inefficient.... Need to rewrite that, open the fits file ONCE !
int write_data_flag_to_dirfile(struct param_common dir, struct samples samples_struct, long iframe_min, long iframe_max, std::vector<std::vector<std::string> > bolo_vect)
{

	double *signal;
	int *flag;
	long ns;
	DIRFILE* D, *H;
	std::vector<string> det_vect;


	fitsfile *fptr;
	int status = 0;
	int rowIndex = 0, naxis = 0, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };


	for (long iframe=iframe_min;iframe<iframe_max;iframe++){

		det_vect=bolo_vect[iframe];

		string base_name = samples_struct.basevect[iframe];
		string fits_name = dir.data_dir + samples_struct.fitsvect[iframe];

		string datadir = dir.tmp_dir + "dirfile/" + base_name + "/data";
		string flagdir = dir.tmp_dir + "dirfile/" + base_name + "/flag";

		D = gd_open((char *) datadir.c_str(), GD_RDWR | GD_CREAT |
				GD_VERBOSE | GD_UNENCODED);
		H = gd_open((char *) flagdir.c_str(), GD_RDWR | GD_CREAT |
				GD_VERBOSE | GD_UNENCODED);


		base_name = FitsBasename(samples_struct.basevect[iframe]); //TODO: needed ???


		if (fits_open_file(&fptr, fits_name.c_str(), READONLY, &status)){
			fits_report_error(stderr, status);
			return 1;
		}

		for(long idet=0; idet < (long)det_vect.size(); idet ++){


			string field = det_vect[idet];
			string data_outfile = "data_" + base_name + "_" + field;
			string flag_outfile = "flag_" + base_name + "_" + field;

			// Retrieve the row of the specified channel
			rowIndex = find_channel_index(fptr, field.c_str());

			// Move ptr to signal hdu
			if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", 0, &status)){
				fits_report_error(stderr, status);
				return 1;
			}

			// Retrieve the size of the signal
			if (fits_get_img_dim(fptr, &naxis, &status)){
				fits_report_error(stderr, status);
				return 1;
			}
			if(naxis != 2){
				fits_report_error(stderr, status);
				return 1;
			}
			if (fits_get_img_size(fptr, 2, naxes, &status)){
				fits_report_error(stderr, status);
				return 1;
			}

			// Allocate Memory
			ns = naxes[0];
			signal = new double[ns];
			flag   = new int[ns];

			// Retrieve the corresponding row
			fpixel[0] = 1;
			fpixel[1] = rowIndex;
			if (fits_read_pix(fptr, TDOUBLE, fpixel, ns, 0, signal, &anynul, &status)){
				fits_report_error(stderr, status);
				return 1;
			}

			// Move ptr to mask hdu ...
			if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", 0, &status)){
				fits_report_error(stderr, status);
				return 1;
			}

			// Retrieve the size of the mask
			if (fits_get_img_dim(fptr, &naxis, &status)){
				fits_report_error(stderr, status);
				return 1;
			}
			if(naxis != 2){
				fits_report_error(stderr, status);
				return 1;
			}
			if (fits_get_img_size(fptr, 2, naxes, &status)){
				fits_report_error(stderr, status);
				return 1;
			}
			if (ns != naxes[0])
				return 1;

			if (fits_read_pix(fptr, TINT, fpixel, ns, 0, flag, &anynul, &status)){
				fits_report_error(stderr, status);
				return 1;
			}

			//configure dirfile field
			gd_entry_t E;
			E.field = (char*) data_outfile.c_str();
			E.field_type = GD_RAW_ENTRY;
			E.fragment_index = 0;
			E.spf = 1;
			E.data_type = GD_DOUBLE;
			E.scalar[0] = NULL;

			// add to the dirfile
			gd_add(D, &E);
			//			gd_flush(D,NULL);

			// write binary file on disk
			int n_write = gd_putdata(D, (char*) data_outfile.c_str(), 0, 0, 0, ns, GD_DOUBLE, signal);
			if (gd_error(D) != 0) {
				cout << "error putdata in write_data_flag : wrote " << n_write
						<< " and expected " << ns << endl;
				return 1;
			}

			//configure dirfile field
			gd_entry_t F;
			F.field = (char*) flag_outfile.c_str();
			F.field_type = GD_RAW_ENTRY;
			F.fragment_index = 0;
			F.spf = 1;
			F.data_type = GD_INT32;
			F.scalar[0] = NULL;

			// add to the dirfile
			gd_add(H, &F);

			// write binary file on disk
			n_write = gd_putdata(H, (char*) flag_outfile.c_str(), 0, 0, 0, ns, GD_INT32, flag);
			if (gd_error(H) != 0) {
				cout << "error putdata in write_data_flag : wrote " << n_write
						<< " and expected " << ns << endl;
				return 1;
			}

			delete [] flag;
			delete [] signal;

			gd_flush(D, data_outfile.c_str());
			gd_flush(H, flag_outfile.c_str());
		}


		// close fits file
		if(fits_close_file(fptr, &status)){
			fits_report_error(stderr, status);
			return 1;
		}

		// close dirfile
		if (gd_close(D)) {
			cout << "Dirfile gd_close error in write_data_flag for : " << datadir
					<< endl;
			return 1;
		}

		// close dirfile
		if (gd_close(H)) {
			cout << "Dirfile gd_close error in write_data_flag for : " << flagdir
					<< endl;
			return 1;
		}

	}

	return 0;
}

//TODO: What about the other format ???
int write_LON_LAT_to_dirfile(struct param_common dir, struct samples samples_struct, long iframe_min, long iframe_max, std::vector<std::vector<std::string> > bolo_vect)
{

	double *lon;
	double *lat;
	long ns;
	DIRFILE* D, *H;


	fitsfile *fptr;
	int status = 0;
	int rowIndex = 0, naxis = 0, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };


	std::vector<string> det_vect;

	for (long iframe=iframe_min;iframe<iframe_max;iframe++){

		det_vect=bolo_vect[iframe];

		string base_name = samples_struct.basevect[iframe];
		string fits_name = dir.data_dir + samples_struct.fitsvect[iframe];

		string LONdir = dir.tmp_dir + "dirfile/" + base_name + "/LON";
		string LATdir = dir.tmp_dir + "dirfile/" + base_name + "/LAT";

		D = gd_open((char *) LONdir.c_str(), GD_RDWR | GD_CREAT |
				GD_VERBOSE | GD_UNENCODED); // | GD_TRUNC |
		H = gd_open((char *) LATdir.c_str(), GD_RDWR | GD_CREAT |
				GD_VERBOSE | GD_UNENCODED); // | GD_TRUNC |


		base_name = FitsBasename(samples_struct.basevect[iframe]); //TODO: needed ???

		if (fits_open_file(&fptr, fits_name.c_str(), READONLY, &status)){
			fits_report_error(stderr, status);
			return 1;
		}

		for(long idet=0; idet < (long)det_vect.size(); idet ++){

			string field = det_vect[idet];
			string lon_outfile = "LON_" + base_name + "_" + field;
			string lat_outfile = "LAT_" + base_name + "_" + field;

			// Retrieve the row of the specified channel
			rowIndex = find_channel_index(fptr, field.c_str());

			// Move ptr to lon hdu
			if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "lon", 0, &status)){
				fits_report_error(stderr, status);
				return 1;
			}

			// Retrieve the size of the signal
			if (fits_get_img_dim(fptr, &naxis, &status)){
				fits_report_error(stderr, status);
				return 1;
			}
			if(naxis != 2){
				fits_report_error(stderr, status);
				return 1;
			}
			if (fits_get_img_size(fptr, 2, naxes, &status)){
				fits_report_error(stderr, status);
				return 1;
			}

			// Allocate Memory
			ns = naxes[0];
			lon = new double[ns];
			lat = new double[ns];

			// ---------------------------------------------
			// Retrieve the corresponding row
			fpixel[0] = 1;
			fpixel[1] = rowIndex;
			if (fits_read_pix(fptr, TDOUBLE, fpixel, ns, 0, lon, &anynul, &status)){
				fits_report_error(stderr, status);
				return 1;
			}


			// ---------------------------------------------
			// Move ptr to lat hdu
			if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "lat", 0, &status)){
				fits_report_error(stderr, status);
				return 1;
			}

			// Retrieve the size of the signal
			if (fits_get_img_dim(fptr, &naxis, &status)){
				fits_report_error(stderr, status);
				return 1;
			}
			if(naxis != 2){
				fits_report_error(stderr, status);
				return 1;
			}
			if (fits_get_img_size(fptr, 2, naxes, &status)){
				fits_report_error(stderr, status);
				return 1;
			}
			if (ns != naxes[0])
				return 1;

			// ---------------------------------------------
			// Retrieve the corresponding row
			if (fits_read_pix(fptr, TDOUBLE, fpixel, ns, 0, lat, &anynul, &status)){
				fits_report_error(stderr, status);
				return 1;
			}


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

			gd_flush(D, lon_outfile.c_str());
			gd_flush(H, lat_outfile.c_str());
		}


		// close fits file
		if(fits_close_file(fptr, &status)){
			fits_report_error(stderr, status);
			return 1;
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


int read_data_from_dirfile(DIRFILE* D, string filename, string field, double *data, long ns){

	// set dirfile name and binary name
	string outfile = "data_" + filename + "_" + field;
	const char * field_code;
	field_code = outfile.c_str();

	// fill data array
	int nget = gd_getdata(D, field_code, 0, 0, 0, ns, GD_DOUBLE,data);
	if (gd_error(D) != 0) {
		cout << "error getdata in read_data_from_dirfile : read " << nget << endl;
		return 1;
	}

	gd_flush(D,field_code);

	return 0;
}

int read_flag_from_dirfile(DIRFILE* H, string filename, string field, int *mask, long ns){

	// set dirfile name and binary name
	string outfile = "flag_" + filename + "_" + field;
	const char * field_code;
	field_code = outfile.c_str();

	// fill mask array
	int nget = gd_getdata(H, field_code, 0, 0, 0, ns, GD_INT32,	mask);
	if (gd_error(H) != 0) {
		cout << "error getdata in read_flag_from_dirfile : read " << nget << endl;
		return 1;
	}

	gd_flush(H,field_code);

	return 0;
}


int read_LON_from_dirfile(DIRFILE* D, string filename, string field, double *lon, long ns){

	// set dirfile name and binary name
	string outfile = "LON_" + filename + "_" + field;
	const char * field_code;
	field_code = outfile.c_str();

	// fill ra array
	int nget = gd_getdata(D, field_code, 0, 0, 0, ns, GD_DOUBLE, lon);
	if (gd_error(D) != 0) {
		cout << "error getdata in read_LON_from_dirfile : read " << nget << endl;
		return 1;
	}

	gd_flush(D,field_code);

	return 0;
}

int read_LAT_from_dirfile(DIRFILE* D, string filename, string field, double *lat, long ns){

	// set dirfile name and binary name
	string outfile = "LAT_" + filename + "_" + field;
	const char * field_code;
	field_code = outfile.c_str();

	// fill lat array
	int nget = gd_getdata(D, field_code, 0, 0, 0, ns, GD_DOUBLE, lat);
	if (gd_error(D) != 0) {
		cout << "error getdata in read_LAT_from_dirfile : read " << nget << endl;
		return 1;
	}

	gd_flush(D,field_code);

	return 0;
}


int write_samptopix(DIRFILE *D, long ns, long long *&samptopix,
		string filename, std::string boloname)
/*!  write a sample to pixel vector to disk  */
{
	string outfile = filename + "_" + boloname;

	// write binary file on disk
	int n_write = gd_putdata(D, (char*) outfile.c_str(), 0, 0, 0, ns, GD_INT64,
			samptopix);
	if (gd_error(D) != 0) {
		cout << "error putdata in write_samptopix : wrote " << n_write
				<< " and expected " << ns << endl;
		return 1;
	}

	gd_flush(D,NULL);

	return 0;
}

int read_samptopix(DIRFILE* D, long ns, long long *samptopix,
		string filename, std::string boloname)
/*!  read a sample to pixel vector from disk  */
{
	// set binary file name
	string outfile = filename + "_" + boloname;
	const char * field_code;
	field_code = outfile.c_str();

	// fill samptopix array
	int nget = gd_getdata(D, field_code, 0, 0, 0, ns, GD_INT64,	samptopix);
	if (gd_error(D) != 0 || nget != ns) {
		cout << "error getdata in read_samptopix : read " << nget << endl;
		return 1;
	}

	gd_flush(D,field_code);

	return 0;
}

int write_indpsrc(long long map_size, long long npixsrc, long long * indpsrc,
		std::string outdir)
/*!  write the bright sources index to disk  */
{
	FILE *fp;
	string testfile = outdir + "indpsrc.bin";

	if ((fp = fopen(testfile.c_str(), "w")) != NULL) {
		fwrite(&map_size, sizeof(long long), 1, fp);
		fwrite(&npixsrc, sizeof(long long), 1, fp);
		fwrite(indpsrc, sizeof(long long), map_size, fp);
		fclose(fp);
	} else {
		cerr << "ERROR : Could not open " << testfile << endl;
		return 1;
	}

	return 0;
}

int read_indpsrc(long long &map_size, long long &npixsrc, long long *&indpsrc,
		std::string outdir)
/*!  read the bright sources index from disk  */
{
	FILE *fp;
	string testfile;
	size_t result;

	testfile = outdir + "indpsrc.bin";
	if ((fp = fopen(testfile.c_str(), "r")) != NULL) {
		result = fread(&map_size, sizeof(long long), 1, fp);
		result = fread(&npixsrc, sizeof(long long), 1, fp);
		indpsrc = new long long[map_size];
		result = fread(indpsrc, sizeof(long long), map_size, fp);
		fclose(fp);
	} else {
		cerr << "Error : cannot find indpsrc.bin file at " << testfile << endl;
		return 1;
	}

	return 0;
}

int write_indpix(long long ind_size, long long npix, long long *indpix,
		string outdir, int flagon)
/*!  write the map index to disk  */
{
	FILE *fp;
	string testfile2;

	testfile2 = outdir + "Indpix_for_conj_grad.bi";

	if ((fp = fopen(testfile2.c_str(), "w")) != NULL) {
		fwrite(&flagon, sizeof(int), 1, fp); // mat 04/06
		fwrite(&npix, sizeof(long long), 1, fp);
		fwrite(&ind_size, sizeof(long long), 1, fp);
		fwrite(indpix, sizeof(long long), ind_size, fp);
		fclose(fp);
	} else {
		cerr << "ERROR : Could not open " << testfile2 << endl;
		return 1;
	}

#ifdef DEBUG
	testfile2 = outdir + "Indpix_for_conj_grad.txt";
	if((fp = fopen(testfile2.c_str(),"w"))) {
		fprintf(fp,"%d\n",flagon);
		fprintf(fp,"%lld\n",npix);
		fprintf(fp,"%lld\n",ind_size);
		for(long ii =0;ii<ind_size;ii++)
			fprintf(fp,"%lld ",indpix[ii]);
		fclose(fp);
	} else {
		cerr << "ERROR : Could not find " << testfile2 << endl;
		return 1;
	}

#endif

	return 0;
}

int read_indpix(long long &ind_size, long long &npix, long long *&indpix,
		string outdir, int &flagon)
/*!  read the map index from disk  */
{
	FILE *fp;
	string testfile2;
	size_t result;

	testfile2 = outdir + "Indpix_for_conj_grad.bi";
	if ((fp = fopen(testfile2.c_str(), "r")) != NULL) {
		result = fread(&flagon, sizeof(int), 1, fp);
		result = fread(&npix, sizeof(long long), 1, fp);
		result = fread(&ind_size, sizeof(long long), 1, fp);
		indpix = new long long[ind_size];
		result = fread(indpix, sizeof(long long), ind_size, fp);
		fclose(fp);
	} else {
		cerr << "Error : cannot find Indpix file " << testfile2 << endl;
		return 1;
	}

	return 0;
}

int write_PNd(double *PNd, long long npix, string outdir, std::string filename)
/*!  write the map preconditioner to disk  */
{
	FILE *fp;
	string testfile2;

	testfile2 = outdir + filename;

	if ((fp = fopen(testfile2.c_str(), "w")) != NULL) {
		fwrite(&npix, sizeof(long long), 1, fp);
		fwrite(PNd, sizeof(double), npix, fp);
		fclose(fp);
	} else {
		cerr << "ERROR : Could not find " << testfile2 << endl;
		return 1;
	}

#ifdef DEBUG
	testfile2 = outdir + "PNdCorr.txt";

	if((fp = fopen(testfile2.c_str(),"w"))) {
		fprintf(fp,"%lld\n",npix);
		for(int ii= 0; ii<npix;ii++)
			fprintf(fp,"%lf ",PNd[ii]);
		fclose(fp);
	} else {
		cerr << "ERROR : Could not find " << testfile2 << endl;
		return 1;
	}

#endif
	return 0;
}

// correlation between npix here and npix in Indpix file is done in sanePic (main.cpp)
int read_PNd(double *&PNdtot, long long &npix, string outdir, std::string filename)
/*!  read the map preconditioner from disk  */
{
	FILE *fp;
	string testfile2;
	size_t result;

	testfile2 = outdir + filename;

	if ((fp = fopen(testfile2.c_str(), "r")) != NULL) {
		result = fread(&npix, sizeof(long long), 1, fp);
		PNdtot = new double[npix];
		result = fread(PNdtot, sizeof(double), npix, fp);
		fclose(fp);
	} else {
		cerr << "Error. Unable to read file : " << testfile2 << endl;
		return 1;
	}

	return 0;
}

int write_fdata(DIRFILE *D, long ns, fftw_complex *fdata, string prefixe,
		long idet, string filename, std::vector<std::string> bolonames)
/*! write Fourier data file to disk */
{

	// set dirfile and binary names
	string outfile = prefixe + filename + "_" + bolonames[idet];
	const char * field_code;
	field_code = outfile.c_str();

	// write binary file on disk
	int n_write = gd_putdata(D, field_code, 0, 0, 0, ns/2+1, GD_COMPLEX128,
			fdata);
	if (gd_error(D) != 0) {
		cout << "error putdata in write_fdata : wrote " << n_write
				<< " and expected " << ns << endl;
		return 1;
	}

	gd_flush(D,field_code);

	return 0;
}

int read_fdata(DIRFILE* D, long ns, fftw_complex *fdata, string prefixe,
		long idet, string filename, std::vector<std::string> bolonames)
/*!  read the map preconditioner from disk  */
{

	// set dirfile and bin names
	string outfile = prefixe + filename + "_" + bolonames[idet];
	const char * field_code;
	field_code = outfile.c_str();

	// fill fdata with binary
	int nget = gd_getdata(D, field_code, 0, 0, 0, ns/2+1, GD_COMPLEX128, fdata);
	if (gd_error(D) != 0) {
		cout << "error getdata in read_fdata : read " << nget << endl;
		return 1;
	}

	gd_flush(D,field_code);

	return 0;
}

int read_mixmat_txt(string MixMatfile, long ndet, long ncomp, double **&mixmat) {
	FILE *fp;
	int result;
	long ncomp2;
	double dummy1; // used to read mixing matrix

	if ((fp = fopen(MixMatfile.c_str(), "r")) == NULL) {
		cerr << "ERROR: Can't find Mixing Matrix file. Exiting. \n";
		cout
		<< "Advice : verify the file is in your noise directory and that his name is : "
		<< MixMatfile << endl;
		return 1;
	}
	result = fscanf(fp, "%ld", &ncomp2);

	mixmat = dmatrix(0, ndet - 1, 0, ncomp - 1);

	for (long ii = 0; ii < ndet; ii++) {
		for (long jj = 0; jj < ncomp2; jj++) {
			result = fscanf(fp, "%lf", &dummy1);
			if (jj < ncomp)
				mixmat[ii][jj] = dummy1;
		}
	}
	fclose(fp);
	return 0;

}

