
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "struct_definition.h"
#include "inputFileIO.h"
#include <iostream>

#include "temporary_IO.h"
#include "dataIO.h"

extern "C" {
#include "nrutil.h"
#include "getdata.h"
}

using namespace std;

int write_data_flag_to_dirfile(struct param_common dir, struct samples samples_struct) {

	string datadir = dir.tmp_dir + "dirfile/data/";
	string flagdir = dir.tmp_dir + "dirfile/flag/";

	DIRFILE* D = gd_open((char *) datadir.c_str(), GD_RDWR | GD_CREAT
			| GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
	DIRFILE* H = gd_open((char *) flagdir.c_str(), GD_RDWR | GD_CREAT
			| GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);

	double *d;
	int *flag;
	long ns;

	//	// get format number of fields
	//	unsigned int nfields =  gd_nfields(D);
	//
	//	// get list of fields
	//	const char** field_list = gd_field_list(D);
	//
	//	// fill a vector with fields names
	//	std::vector<std::string> fields((char**)field_list, (char**)(field_list+nfields));
	//
	//	// seek whether field already exists
	//	int mycount = (int) count(fields.begin(), fields.end(), outfile);
	//
	//	if(mycount==1){ // delete field + bin in case it already exists
	//		gd_delete(D, (char*)outfile.c_str(), GD_DEL_DATA);
	//		if(gd_error(D)){
	//			cout << "error putdata in write_fdata : gd_delete " << outfile << " failed" << endl;
	//			return 1;
	//		}
	//	}

	for(long iframe = 0; iframe < samples_struct.ntotscan; iframe++){

		string base_name = FitsBasename(samples_struct.fitsvect[iframe]);
		string output_read = "";
		std::vector<string> bolonames;
		if(read_channel_list(output_read, samples_struct.bolovect[iframe], bolonames)){
			cout << output_read << endl;
			return 1;
		}
		long ndet = (long)bolonames.size();

		for(long idet=0; idet < ndet; idet ++){

			string field = bolonames[idet];
			string data_outfile = "data_" + base_name + "_" + field;
			string flag_outfile = "flag_" + base_name + "_" + field;

			read_signal_from_fits(samples_struct.fitsvect[iframe], field, d, ns);
			read_flag_from_fits(samples_struct.fitsvect[iframe], field, flag, ns);

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

			// write binary file on disk
			int n_write = gd_putdata(D, (char*) data_outfile.c_str(), 0, 0, 0, ns, GD_DOUBLE,
					d);
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
			n_write = gd_putdata(H, (char*) flag_outfile.c_str(), 0, 0, 0, ns, GD_INT32,
					flag);
			if (gd_error(H) != 0) {
				cout << "error putdata in write_data_flag : wrote " << n_write
						<< " and expected " << ns << endl;
				return 1;
			}

			delete [] flag;
			delete [] d;
		}
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

	return 0;
}


int read_data_flag_from_dirfile(string tmp_dir, string filename, string field, double *&data, int *&mask){ // TODO :  add ns for size verification ??

	// set dirfile name and binary name
	string datadir = tmp_dir + "dirfile/data/";
	string flagdir = tmp_dir + "dirfile/flag/";
	string data_outfile = "data_" + FitsBasename(filename) + "_" + field;
	string flag_outfile = "flag_" + FitsBasename(filename) + "_" + field;


	// open dirfile
	DIRFILE* D = gd_open((char *) datadir.c_str(), GD_RDWR | GD_VERBOSE
			| GD_UNENCODED | GD_BIG_ENDIAN);
	DIRFILE* H = gd_open((char *) flagdir.c_str(), GD_RDWR | GD_VERBOSE
			| GD_UNENCODED | GD_BIG_ENDIAN);

	// get ns from format
	unsigned int nframe = gd_nframes(D);
	//	cout << " nframesD : " << nframe << endl;
	//	cout << " nframesH : " << nframe << endl;

	data = new double[nframe];

	// fill data array
	int nget = gd_getdata(D, (char*) data_outfile.c_str(), 0, 0, 0, nframe, GD_DOUBLE,
			data);
	if (gd_error(D) != 0) {
		cout << "error getdata in read_data_flag_from_dirfile : read " << nget << endl;
		return 1;
	}


	nframe = gd_nframes(H);

	mask = new int[nframe];

	// fill mask array
	nget = gd_getdata(H, (char*) flag_outfile.c_str(), 0, 0, 0, nframe, GD_INT32,
			mask);
	if (gd_error(H) != 0) {
		cout << "error getdata in read_data_flag_from_dirfile : read " << nget << endl;
		return 1;
	}

	// close dirfile
	if (gd_close(D)) {
		cout << "Dirfile gd_close error in read_data_flag_from_dirfile for : " << datadir
				<< endl;
		return 1;
	}

	// close dirfile
	if (gd_close(H)) {
		cout << "Dirfile gd_close error in read_data_flag_from_dirfile for : " << flagdir
				<< endl;
		return 1;
	}

	return 0;
}

int write_samptopix(long ns, long long *&samptopix, string tmpdir,
		string filename, std::string boloname)
/*!  write a sample to pixel vector to disk  */
{
	// set dirfile path name
	string filedir = tmpdir + "dirfile/Indexes/";
	// set binary file name
	string outfile = FitsBasename(filename) + "_" + boloname;

	// open dirfile
	DIRFILE* D = gd_open((char *) filedir.c_str(), GD_RDWR | GD_CREAT
			| GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);

	//configure dirfile field
	gd_entry_t E;
	E.field = (char*) outfile.c_str();
	E.field_type = GD_RAW_ENTRY;
	E.fragment_index = 0;
	E.spf = 1;
	E.data_type = GD_INT64;
	E.scalar[0] = NULL;

	// add to the dirfile
	gd_add(D, &E);
	if (gd_error(D) != 0) {
		cout << "gd_add error in write_samptopix for file : " << filedir
				+ outfile << endl;
		return 1;
	}

	// write binary file on disk
	int n_write = gd_putdata(D, (char*) outfile.c_str(), 0, 0, 0, ns, GD_INT64,
			samptopix);
	if (gd_error(D) != 0) {
		cout << "error putdata in write_samptopix : wrote " << n_write
				<< " and expected " << ns << endl;
		return 1;
	}

	// close dirfile
	if (gd_close(D)) {
		cout << "Dirfile gd_close error in write_samptopix for : " << filedir
				<< endl;
		return 1;
	}

#ifdef DEBUG_PRINT

	outfile=tmpdir + "/Indexes/samptopix_" + FitsBasename(filename) + "_" + boloname + ".txt";
	FILE *fp;

	if((fp = fopen(outfile.c_str(),"w"))) {
		fprintf(fp,"%ld ",ns);
		for(long ii = 0; ii< ns; ii++)
			fprintf(fp,"%lld ",samptopix[ii]);
		fclose(fp);
	} else {
		cerr << "ERROR : Could not find " << outfile << endl;
		return 1;
	}
#endif

	return 0;
}

int read_samptopix(long ns, long long *&samptopix, string tmpdir,
		string filename, std::string boloname)
/*!  read a sample to pixel vector from disk  */
{
	// set dirfile name and binary name
	string filedir = tmpdir + "dirfile/Indexes/";
	string outfile = FitsBasename(filename) + "_" + boloname;

	// open dirfile
	DIRFILE* D = gd_open((char *) filedir.c_str(), GD_RDWR | GD_VERBOSE
			| GD_UNENCODED | GD_BIG_ENDIAN);

	unsigned int nframe = gd_nframes(D);

	// fill samptopix array
	int nget = gd_getdata(D, (char*) outfile.c_str(), 0, 0, 0, nframe, GD_INT64,
			samptopix);
	if (gd_error(D) != 0) {
		cout << "error getdata in read_samptopix : read " << nget << endl;
		return 1;
	}

	// close dirfile
	if (gd_close(D)) {
		cout << "Dirfile gd_close error in read_samptopix for : " << filedir
				<< endl;
		return 1;
	}

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

#ifdef DEBUG_PRINT
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

int write_PNd(double *PNd, long long npix, string outdir)
/*!  write the map preconditioner to disk  */
{
	FILE *fp;
	string testfile2;

	testfile2 = outdir + "PNdCorr.bi";

	if ((fp = fopen(testfile2.c_str(), "w")) != NULL) {
		fwrite(&npix, sizeof(long long), 1, fp);
		fwrite(PNd, sizeof(double), npix, fp);
		fclose(fp);
	} else {
		cerr << "ERROR : Could not find " << testfile2 << endl;
		return 1;
	}

#ifdef DEBUG_PRINT
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
int read_PNd(double *&PNdtot, long long &npix, string outdir)
/*!  read the map preconditioner from disk  */
{
	FILE *fp;
	string testfile2;
	size_t result;

	testfile2 = outdir + "PNdCorr.bi";

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

int write_fdata(long ns, fftw_complex *fdata, string prefixe, string tmpdir,
		long idet, string filename, std::vector<std::string> bolonames)
/*! write Fourier data file to disk */
{
	// set dirfile and binary names
	string filedir = tmpdir + "dirfile/Fourier_transform/";
	string outfile = prefixe + FitsBasename(filename) + "_" + bolonames[idet];

	//open dirfile
	DIRFILE* D = gd_open((char *) filedir.c_str(), GD_RDWR | GD_VERBOSE
			| GD_UNENCODED | GD_BIG_ENDIAN);

	// get format number of fields
	unsigned int nfields = gd_nfields(D);

	// get list of fields
	const char** field_list = gd_field_list(D);

	// fill a vector with fields names
	std::vector<std::string> fields((char**) field_list, (char**) (field_list
			+ nfields));

	// seek whether field already exists
	int mycount = (int) count(fields.begin(), fields.end(), outfile);

	if (mycount == 1) { // delete field + bin in case it already exists
		gd_delete(D, (char*) outfile.c_str(), GD_DEL_DATA);
		if (gd_error(D)) {
			cout << "error putdata in write_fdata : gd_delete " << outfile
					<< " failed" << endl;
			return 1;
		}
	}

	//configure dirfile field
	gd_entry_t E;
	E.field = (char*) outfile.c_str();
	E.field_type = GD_RAW_ENTRY;
	E.fragment_index = 0;
	E.spf = 1;
	E.data_type = GD_DOUBLE;
	E.scalar[0] = NULL;

	// add to the dirfile
	gd_add(D, &E);

	// write binary file on disk
	int n_write = gd_putdata(D, (char*) outfile.c_str(), 0, 0, 0, ns, GD_DOUBLE,
			fdata);
	if (gd_error(D) != 0) {
		cout << "error putdata in write_fdata : wrote " << n_write
				<< " and expected " << ns << endl;
		return 1;
	}

	// close dirfile
	if (gd_close(D)) {
		cout << "Dirfile gd_close error in write_fdata for : " << filedir
				<< endl;
		return 1;
	}

#ifdef DEBUG_PRINT
	FILE *fp;
	outfile=tmpdir + "/Fourier_data/" + prefixe + FitsBasename(filename) + "_" + bolonames[idet] + ".txt";
	if((fp = fopen(outfile.c_str(),"w"))) {
		data_size = (ns/2+1)*2;
		fprintf(fp,"%ld ",(long)data_size);
		for(int ii=0;ii<(ns/2+1);ii++) {
			fprintf(fp,"%lf ",fdata[ii][0]);
			fprintf(fp,"%lf \n",fdata[ii][1]);}
		fclose(fp);
	} else {
		cerr << "ERROR : Could not open " << outfile << endl;
		return 1;
	}
#endif

	return 0;
}

int read_fdata(long ns, fftw_complex *&fdata, string prefixe, string tmpdir,
		long idet, string filename, std::vector<std::string> bolonames)
/*!  read the map preconditioner from disk  */
{

	// set dirfile and bin names
	string filedir = tmpdir + "dirfile/Fourier_transform/";
	string outfile = prefixe + FitsBasename(filename) + "_" + bolonames[idet];

	// open dirfile
	DIRFILE* D = gd_open((char *) filedir.c_str(), GD_RDWR | GD_VERBOSE
			| GD_UNENCODED | GD_BIG_ENDIAN);

	unsigned int nframe = gd_nframes(D);

	// fill fdata with binary
	int nget = gd_getdata(D, (char*) outfile.c_str(), 0, 0, 0, nframe, GD_DOUBLE,
			fdata);
	if (gd_error(D) != 0) {
		cout << "error getdata in read_fdata : read " << nget << endl;
		return 1;
	}

	// close dirfile
	if (gd_close(D)) {
		cout << "Dirfile gd_close error in read_fdata for : " << filedir
				<< endl;
		return 1;
	}

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

