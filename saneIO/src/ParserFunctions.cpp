#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype> // for toupper
#include <typeinfo>
#include <sys/types.h>  // For stat()
#include <sys/stat.h>   // For stat()
#include <typeinfo>

#include <wordexp.h>     // To perform word expansion in path (system variable)

extern "C" {
#include "iniparser.h"
#include "dictionary.h"
#include "getdata.h"
}

#include "Utilities.h"
#include "InputFileIO.h"
#include "ParserFunctions.h"
#include "MPIConfiguration.h"
#include "Crc.h"
#include "ErrorCode.h"


#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;

string checkTrailingDir(string str) {
	if (str[str.length() - 1] != '/')
		str = str + '/';

	return str;
}

string expandDir(string str){
	string output;

	wordexp_t p;
	wordexp (str.c_str(), &p, 0);
	if (p.we_wordc != 1)
		cerr << "EE - Problem with directory expansion " +str << endl;
	output = p.we_wordv[0];

	//	  for (size_t i = 0; i < p.we_wordc; i++)
	//	    cout << w[i] << endl;

	wordfree (&p);
	return output;

}

string checkDir(string str){
	return checkTrailingDir(expandDir(str));
}

uint16_t check_path(string &output, string strPath, bool create, int rank) {

	if (access(strPath.c_str(), 0) == 0) {
		struct stat status;
		stat(strPath.c_str(), &status);

		if (status.st_mode & S_IFDIR) {
			return OK;
		} else {
			output += "EE - " + strPath + " is a file.\n";
			return DIR_PROBLEM;
		}
	} else {
		if (! create) {
			output += "EE - " + strPath + " does not exist\n";
			return DIR_PROBLEM;
		}

		if (rank == 0) {
			int status;
			status = mkdir(strPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			if (status == 0)
				output += "WW - " + strPath + " created\n";
			else {
				output += "EE - " + strPath + " failed to create\n";
				return DIR_PROBLEM;
			}
		}
	}
	return OK;
}


uint16_t init_tmpdir(string &output, struct samples &samples_struct, std::string dir, int rank) {

	uint16_t returnCode = 0;

	if (check_path(output, dir, true, rank))
		returnCode |= TMP_PATH_PROBLEM;
	if (check_path(output, dir + "dirfile", true, rank))
		returnCode |= TMP_PATH_PROBLEM;

	return returnCode;

}


uint16_t check_file(string strPath) {

	if (access(strPath.c_str(), 0) == 0) {
		struct stat status;
		stat(strPath.c_str(), &status);

		if ( (status.st_mode & S_IFREG) && (status.st_mode & S_IFLNK) )
			return 0;
		else
			return FILE_PROBLEM;

	} else {
		return FILE_PROBLEM;
	}
}

uint16_t fillvect_double(string & output, double value, string file, string dir, long ntotscan, vector<double> &outputVector){

	vector<double> dummy;
	outputVector.resize(ntotscan);

	if (file != "" && value <= 0.0){

		if ( read_file(output, dir + file, dummy) )
			return FILE_PROBLEM;
		if ( dummy.size() != (long unsigned) ntotscan )
			return FILE_SIZE_PROBLEM;

		outputVector = dummy;
	}

	if ( value > 0.0 ) {
		for (long ii=0; ii < ntotscan; ii++)
			outputVector[ii] = value;
	}

	return OK;
}

void fillvect_strings(string commonFile, vector<string> FitsFilename, string suffix, string dir, vector<string> & outputVector){

	unsigned long size;

	size = FitsFilename.size();
	outputVector.resize(size);

	if (commonFile != ""){

		for (unsigned long ii=0; ii < size; ii++)
			outputVector[ii] = dir + commonFile;

	} else {

		for (unsigned long ii=0; ii < size; ii++)
			outputVector[ii] = dir + FitsBasename(FitsFilename[ii])+ suffix;

	}

}

uint16_t fill_samples_struct(string &output, struct samples &samples_struct,
		struct param_common &dir, struct param_saneInv &Inv_param,
		struct param_saneProc &Proc_param, int rank, int size) {

	uint16_t returnCode = 0;

	if (rank == 0)
		returnCode |= readFitsList(output, dir.input_dir + dir.fits_filelist, samples_struct);


#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD); // other procs wait untill rank 0 has read the fits_list
	// Then we need to BCast the read vectors
	returnCode |=  MPI_Bcast_vector_string(samples_struct.fitsvect,  0,  MPI_COMM_WORLD);
	returnCode |=  MPI_Bcast_vector_long(samples_struct.scans_index, 0,  MPI_COMM_WORLD);

	returnCode |=  MPI_Bcast(&(samples_struct.framegiven), 1, MPI_INT,  0, MPI_COMM_WORLD);
	returnCode |=  MPI_Bcast(&(samples_struct.ntotscan),   1, MPI_LONG, 0, MPI_COMM_WORLD);

#endif


#ifdef USE_MPI
	// Default dummy value...
	samples_struct.parallel_scheme = -1;
	if (dir.parallel_scheme.compare("bolo")  == 0) samples_struct.parallel_scheme = 0;
	if (dir.parallel_scheme.compare("frame") == 0) samples_struct.parallel_scheme = 1;
#endif


	// Default values
	samples_struct.iframe_min = 0;
	samples_struct.iframe_max = samples_struct.ntotscan;

	samples_struct.nsamples.clear();
	samples_struct.nsamples.assign(samples_struct.ntotscan, 0);

	samples_struct.nbins.clear();
	samples_struct.nbins.assign(samples_struct.ntotscan, 0);

	samples_struct.ndet.clear();
	samples_struct.ndet.assign(samples_struct.ntotscan, 0);


	samples_struct.dirfile_pointers.clear();
	samples_struct.dirfile_pointers.resize(samples_struct.ntotscan, NULL);

	// Add data directory to fitsvect
	for (long iframe = 0; iframe < samples_struct.ntotscan; iframe++)
		samples_struct.fitsvect[iframe] = dir.data_dir + samples_struct.fitsvect[iframe];

	// Fill basevect
	for (long iframe = 0; iframe < samples_struct.ntotscan; iframe++) {
		samples_struct.basevect.push_back(
				dirfile_Basename(samples_struct.fitsvect[iframe]));
	}

	vector<string> dummy_string;
	vector<double> dummy_double;

	fillvect_strings(dir.bolo, samples_struct.fitsvect, dir.bolo_suffix, dir.input_dir, dummy_string);
	samples_struct.bolovect = dummy_string;
	dummy_string.clear();

	fillvect_strings(Inv_param.cov_matrix, samples_struct.fitsvect, Inv_param.cov_matrix_suffix, Inv_param.noise_dir, dummy_string);
	samples_struct.noisevect = dummy_string;
	dummy_string.clear();

	returnCode |= fillvect_double(output, Proc_param.fcut, Proc_param.fcut_file, dir.input_dir, samples_struct.ntotscan, dummy_double);
	samples_struct.fcut = dummy_double;
	dummy_double.clear();

	returnCode |= fillvect_double(output, Proc_param.fsamp,Proc_param.fsamp_file, dir.input_dir, samples_struct.ntotscan, dummy_double);
	samples_struct.fsamp = dummy_double;
	dummy_double.clear();

	returnCode |= fillvect_double(output, Proc_param.fhp, Proc_param.fhp_file, dir.input_dir, samples_struct.ntotscan, dummy_double);
	samples_struct.fhp = dummy_double;
	dummy_double.clear();

	// Read all channels list once for all
	returnCode |= fill_channel_list(output, samples_struct, rank, size);

	return returnCode;

}

uint16_t fill_channel_list(std::string &output, struct samples &samples_struct, int rank, int size){
	uint16_t returnCode = 0;

	std::vector<string> det_vect;

	for (long iframe = 0; iframe < samples_struct.ntotscan; iframe++) {

		if (rank == 0)
			returnCode |= readChannelList(output, samples_struct.bolovect[iframe], det_vect);

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD); // other procs wait untill rank 0 has read the fits_list
		// Then we need to BCast the read vectors
		returnCode |=  MPI_Bcast_vector_string(det_vect, 0, MPI_COMM_WORLD);
#endif

		samples_struct.bolo_list.push_back(det_vect);

		det_vect.clear();
	}

	return returnCode;
}



// Structure reading

void read_common(string &output, dictionary *ini, struct param_common &common) {

	char *s;
	string output2 = "";

	s = iniparser_getstring(ini, "common:data_directory", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "common:data_directory : default value [" + StringOf(
				common.data_dir) + "]\n";
	else
		common.data_dir = checkDir(StringOf(s));

	s = iniparser_getstring(ini, "common:input_directory", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "common:input_directory : default value [" + StringOf(
				common.input_dir) + "]\n";
	else
		common.input_dir = checkDir(StringOf(s));

	s = iniparser_getstring(ini, "common:output_dir", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "common:output_dir : default value [" + StringOf(
				common.output_dir) + "]\n";
	else
		common.output_dir = checkDir(StringOf(s));

	s = iniparser_getstring(ini, "common:temp_dir", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "common:temp_dir : default value [" + StringOf(
				common.tmp_dir) + "]\n";
	else
		common.tmp_dir = checkDir(StringOf(s));

	s = iniparser_getstring(ini, "common:fits_filelist", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "common:fits_filelist : default value [" + StringOf(
				common.fits_filelist) + "]\n";
	else
		common.fits_filelist = StringOf(s);

	s = iniparser_getstring(ini, "common:bolo_suffix", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "common:bolos_suffix : default value [" + StringOf(
				common.bolo_suffix) + "]\n";
	else
		common.bolo_suffix = StringOf(s);

	s = iniparser_getstring(ini, "common:bolo", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "common:bolos : default value [" + StringOf(
				common.bolo) + "]\n";
	else
		common.bolo = StringOf(s);

	s = iniparser_getstring(ini, "common:para", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "common:para: default value [" + StringOf(
				common.parallel_scheme ) + "]\n";
	else {
		string s_lower = StringOf(s);
		transform(s_lower.begin(), s_lower.end(), s_lower.begin(), (int(*)(int)) tolower);
		common.parallel_scheme = s_lower;
	}

#ifdef DEBUG
	output += output2;
#endif

}

void read_param_sanePos(string &output, dictionary *ini, struct param_sanePos &Pos_param) {

	char *s;
	double d;
	int i;
	string output2 = "";

	d = iniparser_getdouble(ini, (char*) "sanePos:lon", -1.0);
	if (d == -1.0)
		output2 += "sanePos:lon : default value [" + StringOf(
				Pos_param.lon) + "]\n";
	else
		Pos_param.lon = d;

	d = iniparser_getdouble(ini, (char*) "sanePos:lat", -1.0);
	if (d == -1.0)
		output2 += "sanePos:lat : default value [" + StringOf(
				Pos_param.lat) + "]\n";
	else
		Pos_param.lat = d;

	d = iniparser_getdouble(ini, (char*) "sanePos:pixsize", -1.0);
	if (d == -1.0)
		output2 += "sanePos:pixsize : default value [" + StringOf(
				Pos_param.pixdeg) + "]\n";
	else
		Pos_param.pixdeg = d;

	s = iniparser_getstring(ini, (char*) "sanePos:proj_code", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePos:proj_code : default value [" + StringOf(
				Pos_param.projcode) + "]\n";
	else
		Pos_param.projcode = StringOf(s);

	s = iniparser_getstring(ini, (char*) "sanePos:axis_type", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePos:axis_type : default value [" + StringOf(
				Pos_param.axistype) + "]\n";
	else
		Pos_param.axistype = StringOf(s);

	s = iniparser_getstring(ini, "sanePos:mask_file", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePos:mask_file : default value [" + string(
				Pos_param.maskfile) + "]\n";
	else
		Pos_param.maskfile = StringOf(s);

	i = iniparser_getint(ini, (char*) "sanePos:file_format", -1);
	if (i == -1)
		output2 += "sanePos:file_format : default value [" + StringOf(
				Pos_param.fileFormat) + "]\n";
	else
		Pos_param.fileFormat = i;

	i = iniparser_getboolean(ini, "sanePos:map_flagged_data", -1);
	if (i == -1)
		output2 += "sanePos:flgdupl : default value [" + StringOf(
				Pos_param.flgdupl) + "]\n";
	else
		Pos_param.flgdupl = (bool) i;

	i = iniparser_getboolean(ini, "sanePos:project_gaps", -1);
	if (i == -1)
		output2 += "sanePos:project_gaps : default value [" + StringOf(
				Pos_param.projgaps) + "]\n";
	else
		Pos_param.projgaps = (bool) i;


	i = iniparser_getboolean(ini, "sanePos:eq2gal", -1);
	if (i == -1)
		output2 += "sanePos:eq2gal : default value [" + StringOf(
				Pos_param.eq2gal) + "]\n";
	else
		Pos_param.eq2gal = (bool) i;

	i = iniparser_getboolean(ini, "sanePos:gal2eq", -1);
	if (i == -1)
		output2 += "sanePos:gal2eq : default value [" + StringOf(
				Pos_param.gal2eq) + "]\n";
	else
		Pos_param.gal2eq = (bool) i;

	s = iniparser_getstring(ini, "sanePos:radesys", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePos:radesys : default value [" + string(
				Pos_param.radesys) + "]\n";
	else
		Pos_param.radesys = StringOf(s);

	d = iniparser_getdouble(ini, (char*) "sanePos:equinox", -1.0);
	if (d == -1.0)
		output2 += "sanePos:equinox : default value [" + StringOf(
				Pos_param.equinox) + "]\n";
	else
		Pos_param.equinox = d;

	d = iniparser_getdouble(ini, (char*) "sanePos:restwav", -1.0);
	if (d == -1.0)
		output2 += "sanePos:restwav : default value [" + StringOf(
				Pos_param.restwav) + "]\n";
	else
		Pos_param.restwav = d;

#ifdef DEBUG
	output += output2;
#endif

}

void read_param_saneProc(string &output, dictionary *ini, struct param_saneProc &Proc_param) {

	int i;
	double d;
	char *s;
	string output2 = "";

	i = iniparser_getint(ini, (char*) "saneProc:apodize_Nsamples", -1);
	if (i == -1)
		output2 += "saneProc:apodize_Nsamples : default value [" + StringOf(
				Proc_param.napod) + "]\n";
	else
		Proc_param.napod = i;


	i = iniparser_getboolean(ini, (char*) "saneProc:wisdom", -1);
	if (i == -1)
		output2 += "saneProc:wisdom : default value [" + StringOf(
				Proc_param.wisdom) + "]\n";
	else
		Proc_param.wisdom = (bool) i;


	i = iniparser_getboolean(ini, "saneProc:fill_gap", -1);
	if (i == -1)
		output2 += "saneProc:fill_gap : default value [" + StringOf(
				Proc_param.fill_gap) + "]\n";
	else
		Proc_param.fill_gap = (bool) i;

	d = iniparser_getdouble(ini, (char*) "saneProc:fsamp", -1.0);
	if (d == -1.0)
		output2 += "saneProc:fsamp: default value [" + StringOf(
				Proc_param.fsamp) + "]\n";
	else
		Proc_param.fsamp = d;


	s = iniparser_getstring(ini, "saneProc:fsamp_file", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "saneProc:fsamp_file : default value [" + StringOf(
				Proc_param.fsamp_file) + "]\n";
	else
		Proc_param.fsamp_file = StringOf(s);


	d = iniparser_getdouble(ini, (char*) "saneProc:fhp", -1.0);
	if (d == -1.0)
		output2 += "saneProc:fhp: default value [" + StringOf(
				Proc_param.fhp) + "]\n";
	else
		Proc_param.fhp = d;


	s = iniparser_getstring(ini, "saneProc:fhp_file", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "saneProc:fhp_file : default value [" + StringOf(
				Proc_param.fhp_file) + "]\n";
	else
		Proc_param.fhp_file = StringOf(s);


	i = iniparser_getboolean(ini, "saneProc:linear_baseline", -1);
	if (i == -1)
		output2 += "saneProc:no_baseline: default value [" + StringOf(
				Proc_param.remove_linear) + "]\n";
	else
		Proc_param.remove_linear = (bool) i;

	i = iniparser_getboolean(ini, "saneProc:correlation", -1);
	if (i == -1)
		output2 += "saneProc:correlation: default value [" + StringOf(
				Proc_param.CORRon) + "]\n";
	else
		Proc_param.CORRon = (bool) i;

	i = iniparser_getint(ini, "saneProc:poly_order", -42);
	if (i == -42)
		output2 += "saneProc:poly_order: default value [" + StringOf(
				Proc_param.poly_order) + "]\n";
	else
		Proc_param.poly_order = i;



	d = iniparser_getdouble(ini, (char *) "saneProc:fcut", -1.0);
	if (d == -1.0)
		output2 += "saneProc:fcut : default value [" + StringOf(
				Proc_param.fcut) + "]\n";
	else
		Proc_param.fcut = d;

	s = iniparser_getstring(ini, "saneProc:fcut_file", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "saneProc:fcut_file : default value [" + StringOf(
				Proc_param.fcut_file) + "]\n";
	else
		Proc_param.fcut_file = StringOf(s);

#ifdef DEBUG
	output += output2;
#endif
}

void read_param_saneInv(std::string &output, dictionary *ini, struct param_saneInv &Inv_param) {

	char *s;
	string output2 = "";

	s = iniparser_getstring(ini, "saneInv:cov_matrix", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "saneInv:cov_matrix : default value [" + StringOf(
				Inv_param.cov_matrix) + "]\n";
	else
		Inv_param.cov_matrix = StringOf(s);

	s = iniparser_getstring(ini, "saneInv:cov_matrix_suffix", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "saneInv:cov_matrix_suffix : default value [" + StringOf(
				Inv_param.cov_matrix_suffix) + "]\n";
	else
		Inv_param.cov_matrix_suffix = StringOf(s);

	s = iniparser_getstring(ini, "saneInv:noise_dir", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "saneInv:noise_dir : default value [" + StringOf(
				Inv_param.noise_dir) + "]\n";
	else
		Inv_param.noise_dir = checkDir(StringOf(s));

#ifdef DEBUG
	output += output2;
#endif

}

void read_param_sanePS(std::string &output, dictionary *ini, struct param_sanePS &PS_param) {

	int i;
	char *s;
	string output2 = "";

	i = iniparser_getint(ini, "sanePS:niter", -1);
	if (i == -1)
		output2 += "sanePS:niter: default value [" + StringOf(
				PS_param.niter) + "]\n";
	else
		PS_param.niter= i;


	i = iniparser_getint(ini, "sanePS:ncomp", -1);
	if (i == -1)
		output2 += "sanePS:ncomp : default value [" + StringOf(
				PS_param.ncomp) + "]\n";
	else
		PS_param.ncomp = i;

	s = iniparser_getstring(ini, "sanePS:map_file", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePS:map_file : default value [" + StringOf(
				PS_param.signame) + "]\n";
	else
		PS_param.signame = StringOf(s);

	s = iniparser_getstring(ini, "sanePS:MixingMatrix_Suffix", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePS:MixingMatrix_Suffix : default value [" + StringOf(
				PS_param.mix_suffix) + "]\n";
	else
		PS_param.mix_suffix = StringOf(s);

	s = iniparser_getstring(ini, "sanePS:ell_suffix", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePS:ell_suffix : default value [" + StringOf(
				PS_param.ell_suffix) + "]\n";
	else
		PS_param.ell_suffix = StringOf(s);

	s = iniparser_getstring(ini, "sanePS:ell", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePS:ell : default value [" + StringOf(
				PS_param.ell) + "]\n";
	else
		PS_param.ell = StringOf(s);

	s = iniparser_getstring(ini, "sanePS:MixingMatrix",
			(char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePS:MixingMatrix : default value ["
				+ StringOf(PS_param.mix) + "]\n";
	else
		PS_param.mix = StringOf(s);

	i = iniparser_getboolean(ini, "sanePS:save_data", -1);
	if (i == -1)
		output2 += "sanePS:save_data: default value [" + StringOf(
				PS_param.save_data) + "]\n";
	else
		PS_param.save_data = (bool) i;

	//TODO: Ugly turnaround until sanePS is released;
	s = iniparser_getstring(ini, "saneInv:cov_matrix", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "saneInv:cov_matrix : default value [" + StringOf(
				PS_param.cov_matrix) + "]\n";
	else
		PS_param.cov_matrix = StringOf(s);

	s = iniparser_getstring(ini, "saneInv:cov_matrix_suffix", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "saneInv:cov_matrix_suffix : default value [" + StringOf(
				PS_param.cov_matrix_suffix) + "]\n";
	else
		PS_param.cov_matrix_suffix = StringOf(s);

#ifdef DEBUG
	output += output2;
#endif

}

void read_param_sanePic(std::string &output, dictionary *ini, struct param_sanePic &Pic_param){

	int i;
	double d;
	char *s;
	string output2 = "";

	i = iniparser_getint(ini, "sanePic:iterW", -1);
	if (i == -1)
		output2 += "sanePic:iterW : default value [" + StringOf(
				Pic_param.iterw) + "]\n";
	else
		Pic_param.iterw = i;

	d = iniparser_getdouble(ini, (char *) "sanePic:tolerance", -1.0);
	if (d == -1)
		output2 += "sanePic:tolerance : default value [" + StringOf(
				Pic_param.tolerance) + "]\n";
	else
		Pic_param.tolerance = d;

	d = iniparser_getdouble(ini, (char *) "sanePic:subtolerance", -1.0);
	if (d == -1)
		output2 += "sanePic:subtolerance : default value [" + StringOf(
				Pic_param.subtolerance) + "]\n";
	else
		Pic_param.subtolerance = d;

	if (Pic_param.iterw == 0) {
		Pic_param.save_data = 0;
		Pic_param.iterw = 10;
	} else {
		if (Pic_param.iterw < 0)
			Pic_param.save_data = 0;
		else
			Pic_param.save_data = 1;
	}

	i = iniparser_getint(ini, "sanePic:iterMAX", -1);
	if (i <= 0)
		output2 += "sanePic:iterMAX : default value [" + StringOf(
				Pic_param.itermax) + "]\n";
	else
		Pic_param.itermax = i;

	s = iniparser_getstring(ini, "sanePic:map_prefix", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePic:map_prefix : default value [" + StringOf(
				Pic_param.map_prefix) + "]\n";
	else
		Pic_param.map_prefix = StringOf(s);

#ifdef DEBUG
	output += output2;
#endif

}

void read_wisdom(std::string &output, struct param_common &dir, struct param_saneProc Proc_param, int rank){
	// Retrieve wisdom if asked / possible

	if (rank == 0 && Proc_param.wisdom) {
		string filename = dir.tmp_dir + "fftw.wisdom";
		FILE * pFilename;

		pFilename = fopen((const char*) filename.c_str(), "r");

		if ( pFilename != NULL ){
			fftw_import_wisdom_from_file( pFilename );
			fclose(pFilename);
		} else {
#ifdef DEBUG
			output += "WW - no FFTW wisdom imported from fftw.wisdom\n";
#endif
		}
	}

#ifdef USE_MPI
	if ( Proc_param.wisdom ) {

		char * wisdom = NULL;
		int size_wisdom;

		if (rank == 0){
			wisdom = fftw_export_wisdom_to_string();
			size_wisdom = strlen(wisdom);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&size_wisdom, 1, MPI_INT, 0, MPI_COMM_WORLD);

		if (rank != 0){
			wisdom = (char *)calloc(size_wisdom+1, sizeof(char));
			fill(wisdom, wisdom + size_wisdom+1, '\0');
		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(wisdom, size_wisdom, MPI_CHAR, 0, MPI_COMM_WORLD);

		if (rank != 0)
			if ( fftw_import_wisdom_from_string(wisdom) != 0 )
				output += "WW - Problem propagating FFTW wisdom\n";
	}
#endif

}


// Structure checks

uint16_t check_common(string &output, struct param_common &dir, int rank) {

	uint16_t returnCode = 0;

	if (check_path(output, dir.data_dir, false, rank))
		returnCode |= DATA_INPUT_PATHS_PROBLEM;
	if (check_path(output, dir.input_dir, false, rank))
		returnCode |= DATA_INPUT_PATHS_PROBLEM;
	if (check_path(output, dir.output_dir, true, rank))
		returnCode |= OUPUT_PATH_PROBLEM;


	if ((dir.bolo == "") && (dir.bolo_suffix == "")) {
		output += "EE - You must mention one of those parameters :\n";
		output += "     param_common:bolo_suffix or param_common:bolo\n";
		returnCode |= BOLOFILE_NOT_FOUND;
	}

	if ((dir.bolo != "") && (check_file(dir.input_dir+dir.bolo) != 0)){
		output += "EE - " + dir.bolo + " not found\n";
		returnCode |= FILE_PROBLEM;
	}

	if (check_file(dir.input_dir+dir.fits_filelist) != 0){
		output += "EE - " + dir.fits_filelist + " not found\n";
		returnCode |= FILE_PROBLEM;
	}

#ifdef USE_MPI

	if (    (dir.parallel_scheme.compare("bolo")  != 0) &&
			(dir.parallel_scheme.compare("frame") != 0) ) {
		output += "EE - parallel scheme unknown ["+dir.parallel_scheme +"]\n";
		output += "     valid values are (bolo, frame)\n";
		returnCode |= PARA_PROBLEM;
	}
#endif

	return returnCode;
}

uint16_t check_param_sanePos(string &output, struct param_sanePos &Pos_param) {

	uint16_t returnCode = 0;

	if (Pos_param.pixdeg <= 0 && Pos_param.maskfile == "") {
		output += "EE - Pixsize cannot be negative ! \n";
		returnCode |= PIXDEG_WRONG_VALUE;
	}

	if ((Pos_param.fileFormat != 0) && (Pos_param.fileFormat != 1)) {
		output += "EE - Fileformat must be 0 (SANEPIC) or 1 (HIPE) \n";
		returnCode |= FILEFORMAT_NOT_FOUND;
	}

	// Force axis type to GAL if converting from EQ to GAL
	if (Pos_param.eq2gal)
		Pos_param.axistype = "GAL";

	// Force axis type to EQ if converting from GAL to EQ
	if (Pos_param.gal2eq)
		Pos_param.axistype = "EQ";

	return returnCode;
}

uint16_t check_param_saneProc(string &output, struct param_common dir, struct param_saneProc &Proc_param) {

	uint16_t returnCode = 0;

	if (Proc_param.napod < 0) {
		output
		+= "EE - You must choose a positive number of samples to apodize\n";
		returnCode |= NAPOD_WRONG_VALUE;
	}

	if ( ( Proc_param.fsamp <= 0.0 ) && ( Proc_param.fsamp_file == "" ) ) {
		output += "EE - You must mention one of those parameters :\n";
		output += "     saneProc:fsamp or saneProc:fsamp_file\n";
		returnCode |= FSAMP_PROBLEM;
	}

	if ( (Proc_param.fsamp_file != ""  ) && (check_file(dir.input_dir + Proc_param.fsamp_file) != 0) ){
		output += "EE - " + Proc_param.fsamp_file + " not found\n";
		returnCode |= FSAMP_PROBLEM;
	}

	if ( ( Proc_param.fhp < 0.0 ) && ( Proc_param.fhp_file == "" ) ) {
		output += "EE - You must mention one of those parameters :\n";
		output += "     saneProc:fhp or saneProc:fhp_file\n";
		returnCode |= FHP_PROBLEM;
	}

	if ( (Proc_param.fhp_file != ""  ) && (check_file(dir.input_dir + Proc_param.fhp_file) != 0) ){
		output += "EE - " + Proc_param.fhp_file + " not found\n";
		returnCode |= FHP_PROBLEM;
	}

	if ( ( Proc_param.fcut <= 0.0 ) && ( Proc_param.fcut_file == "" ) ) {
		output += "EE - You must mention one of those parameters :\n";
		output += "     saneProc:fcut or saneProc:fcut_file\n";
		returnCode |= FCUT_PROBLEM;
	}

	if ( (Proc_param.fcut_file != ""  ) && (check_file(dir.input_dir + Proc_param.fcut_file) != 0) ){
		output += "EE - " + Proc_param.fcut_file + " not found\n";
		returnCode |= FCUT_PROBLEM;
	}

	if (Proc_param.poly_order >= 0)
		Proc_param.remove_polynomia = true;
	else
		Proc_param.remove_polynomia = false;

	if (Proc_param.fhp > 0 || Proc_param.fhp_file != "" )
		Proc_param.highpass_filter = true;
	else
		Proc_param.highpass_filter = false;

	return returnCode;
}

uint16_t check_param_sanePS(string &output, struct param_common dir, struct param_sanePS &PS_param) {

	uint16_t returnCode = 0;

	if (PS_param.ncomp <= 0) {
		output += "EE - Number of component ncomp cannot be negative or zero ! \n";
		returnCode |= NCOMP_WRONG_VALUE;
	}
	if ((PS_param.ell == "") && (PS_param.ell_suffix == "")) {
		output += "EE - You must mention one of those parameters :\n";
		output += "     sanePS:ell or sanePS:ell_suffix\n";
		returnCode |= ELL_FILE_NOT_FOUND;
	}

	if ( (PS_param.ell != ""  ) && (check_file(dir.input_dir + PS_param.ell) != 0) ){
		output += "EE - " + PS_param.ell + " not found\n";
		returnCode |= ELL_FILE_NOT_FOUND;
	}

	if ((PS_param.mix == "") && (PS_param.mix_suffix == "")) {
		output += "EE - You must mention one of those parameters :\n";
		output += "     sanePS:mix or sanePS:mix_suffix\n";
		returnCode |= MIX_FILE_NOT_FOUND;
	}

	if ( (PS_param.mix != ""  ) && (check_file(dir.input_dir + PS_param.mix) != 0) ){
		output += "EE - " + PS_param.mix + " not found\n";
		returnCode |= MIX_FILE_NOT_FOUND;
	}


	return returnCode;
}

uint16_t check_param_saneInv(string &output, struct param_common dir, struct param_saneInv &Inv_param) {

	uint16_t returnCode = 0;

	if (check_path(output, Inv_param.noise_dir, false, -1))
		returnCode |= SANEINV_INPUT_ERROR;

	if ((Inv_param.cov_matrix == "") && (Inv_param.cov_matrix_suffix == "")) {
		output += "EE - You must mention one of those parameters :\n";
		output += "     saneInv:cov_matrix_suffix or saneInv:cov_matrix\n";
		returnCode |= SANEINV_INPUT_ERROR;
	}

	if ( (Inv_param.cov_matrix != ""  ) && (check_file(dir.input_dir + Inv_param.cov_matrix) != 0) ){
		output += "EE - " + Inv_param.cov_matrix + " not found\n";
		returnCode |= SANEINV_INPUT_ERROR;
	}


	return returnCode;
}

uint16_t check_param_samples(string &output, struct samples &samples_struct){

	uint16_t returnCode = 0;

	for (long iframe=0; iframe < samples_struct.ntotscan; iframe++) {
		if ( (check_file(samples_struct.fitsvect[iframe])  & check_file(samples_struct.fitsvect[iframe] +".gz")) != 0 ){
			output += "EE - " + samples_struct.fitsvect[iframe] + " not found\n";
			returnCode |= FILE_PROBLEM;
		}
		if ( (check_file(samples_struct.noisevect[iframe])  & check_file(samples_struct.noisevect[iframe] +".gz")) != 0 ){
			output += "WW - " + samples_struct.noisevect[iframe] + " not found\n";
		}
		if ( check_file(samples_struct.bolovect[iframe]) != 0 ){
			output += "EE - " + samples_struct.bolovect[iframe] + " not found\n";
			returnCode |= FILE_PROBLEM;
		}

	}
	return returnCode;
}


// Structure defaults value

void default_param(struct param_common &dir,
		struct param_sanePos &Pos_param, struct param_saneProc &Proc_param,
		struct param_saneInv &Inv_param, struct param_sanePic &Pic_param) {

	default_param_common(dir);
	default_param_sanePos(Pos_param);
	default_param_saneProc(Proc_param);
	default_param_saneInv(Inv_param);
	default_param_sanePic(Pic_param);

}

void default_param_common(struct param_common &dir) {

	// param_common

	dir.data_dir = "./";
	dir.output_dir = "./";
	dir.tmp_dir = "./";
	dir.input_dir = "./";

	dir.fits_filelist = "file.list";

	dir.bolo          = "";
	dir.bolo_suffix   = ".bolo";

	dir.parallel_scheme = "frame" ;
}

void default_param_sanePos(struct param_sanePos &Pos_param) {
	//	Pos_param.flgdupl = 0; // map duplication factor
	//	Pos_param.flgdupl = 0; // map duplication factor
	//	Pos_param.maskfile = "";
	// param_sanePos
	Pos_param.maskfile = "";

	Pos_param.pixdeg   = 0.0;
	Pos_param.lon      = NAN;
	Pos_param.lat      = NAN;
	Pos_param.projcode = "TAN";
	Pos_param.axistype = "EQ";

	Pos_param.radesys  = "";  /* will be defaulted later to ICRS */
	Pos_param.equinox  = 0.0; /* will be defaulted later to J2000 */
	Pos_param.restwav  = 0.0;

	Pos_param.eq2gal   = false;
	Pos_param.gal2eq   = false;

	Pos_param.flgdupl  = true; // What to do with flagged data : (default : False -- map in a single pixel)
	Pos_param.projgaps = true; // What to do with gaps : (default : 0 -- no projection)
	Pos_param.fileFormat = 0;  // Default sanepic File Format

	// 0: sanepic format with reference position & offsets
	// 1: 'hipe' like format with LON/LAT for each time/bolo

}

void default_param_saneProc(struct param_saneProc &Proc_param) {

	//	Proc_param.napod  = 0; /*! number of samples to apodize, =0 -> no apodisation */
	//	Proc_param.NOFILLGAP = 0; /*! dont fill the gaps ? default is NO => the program fill */
	//	Proc_param.fsamp = 0.0;// sampling frequency
	//	Proc_param.remove_linear = 0; /*!  baseline is removed from the data, remove_linear = 1 else 0 */
	//	Proc_param.CORRon = 1; /*! correlation included in the analysis (=1), else 0, default 0*/
	//	Proc_param.remove_polynomia = 1; /*! remove a polynomia fitted to the data*/
	//	Proc_param.fhp = 0.0; // low pass filter frequency

	Proc_param.remove_linear = false;
	Proc_param.fill_gap = true;
	Proc_param.CORRon = true;
	Proc_param.remove_polynomia = true;
	Proc_param.highpass_filter = false;

	Proc_param.napod = 100;
	Proc_param.poly_order = 1;

	Proc_param.fsamp = -1.0;
	Proc_param.fhp = -1.0;
	Proc_param.fcut = -1.0;

	Proc_param.fsamp_file = "";
	Proc_param.fhp_file   = "";
	Proc_param.fcut_file  = "";

	Proc_param.wisdom = false;

}

void default_param_sanePS(struct param_sanePS &PS_param) {

	PS_param.ell_suffix = ".ell";
	PS_param.mix_suffix = ".mix";
	PS_param.ell = "";
	PS_param.mix = "";
	PS_param.signame = "";
	PS_param.ncomp = 1;
	PS_param.niter = 500;
	PS_param.save_data = 1;

	//TODO: Ugly turnaround until sanePS is released;

	PS_param.cov_matrix = "";
	PS_param.cov_matrix_suffix = "_psd.fits";

}

void default_param_saneInv(struct param_saneInv &Inv_param) {

	Inv_param.cov_matrix = "";
	Inv_param.cov_matrix_suffix = "_psd.fits";
	Inv_param.noise_dir = "./";
}

void default_param_sanePic(struct param_sanePic &Pic_param) {

	Pic_param.iterw = 0;
	Pic_param.itermax = 2000;
	Pic_param.tolerance =  1e-10;
	Pic_param.subtolerance =  1e-5;
	Pic_param.save_data = 0;
//	Pic_param.restore = 0; should not be defaulted here...
	Pic_param.map_prefix = "optimMap";
}

void fill_sanePS_struct(struct param_sanePS &PS_param,
		struct samples &samples_struct, struct param_common &dir) {

	vector<string> dummy_string;

	fillvect_strings(PS_param.mix, samples_struct.fitsvect, PS_param.mix_suffix, dir.input_dir, dummy_string);
	samples_struct.mix_names = dummy_string;
	dummy_string.clear();

	fillvect_strings(PS_param.ell, samples_struct.fitsvect, PS_param.ell_suffix, dir.input_dir, dummy_string);
	samples_struct.ell_names = dummy_string;
	dummy_string.clear();

}

#ifdef USE_MPI

int commit_dictionary(int rank, dictionary *dict) {

	int n;
	int size;

	char *key_buff=NULL;
	char *val_buff=NULL;

	if(rank==0) {
		n = dict->n;
		size = dict->size;

		MPI_Bcast(&n,    1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

		if((n>0)&&(dict!=NULL)) {
			for (int i=0; i<n; i++) {

				string key_s = dict->key[i];
				string val_s = (dict->val[i] ? dict->val[i] : "UNDEF"); // UNDEF

				int size_key = key_s.size();
				int size_val = val_s.size();

				MPI_Bcast(&size_key, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&size_val, 1, MPI_INT, 0, MPI_COMM_WORLD);

				if(size_key>0) {
					key_buff = (char *)key_s.c_str();
					MPI_Bcast(key_buff, size_key, MPI_CHAR, 0, MPI_COMM_WORLD);
				}

				if(size_val>0) {
					val_buff = (char*)val_s.c_str();
					MPI_Bcast(val_buff, size_val, MPI_CHAR, 0, MPI_COMM_WORLD);
				}
			}
		} else {
			cout << "n=0 or void !" << endl;
			return 1;
		}
	} else {
		int size_key=0;
		int size_val=0;

		MPI_Bcast(&n,    1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

		if((n>0)&&(dict!=NULL)) {
			for (int i=0; i<n; i++) {

				bool delete_val = 1;

				MPI_Bcast(&size_key, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&size_val, 1, MPI_INT, 0, MPI_COMM_WORLD);

				if(size_key>0) {
					key_buff = (char *)calloc(size_key+1, sizeof(char));
					fill(key_buff, key_buff + size_key+1, '\0');
					MPI_Bcast(key_buff, size_key, MPI_CHAR, 0, MPI_COMM_WORLD);
				}

				if(size_val>=0) {
					val_buff = (char *)calloc(size_val+1, sizeof(char));
					fill(val_buff, val_buff + size_val+1, '\0');
					MPI_Bcast(val_buff, size_val, MPI_CHAR, 0, MPI_COMM_WORLD);
				}

				if((size_val<0) || !strcmp(val_buff, (char*)"UNDEF")) {
					val_buff=NULL;
					delete_val=0;
				}

				if(dictionary_set(dict, key_buff, val_buff)) {
					cout << " dictionnary failure" << endl;
					return 1;
				}

				if(size_key>0)
					free(key_buff);
				if(delete_val)
					free(val_buff);

			}
		} else {
			cout << "dictionnary is NULL !" << endl;
			return 1;
		}
	}

	return 0;
}

#endif

//TODO: sanePS is not in the default distribution, so should not be in the master parser_function....
// Main driver
uint16_t parser_function(char * ini_name, std::string &output,
		struct param_common &dir, struct samples &samples_struct,
		struct param_sanePos &Pos_param, struct param_saneProc &Proc_param,
		struct param_sanePS &PS_param, struct param_saneInv &Inv_param,
		struct param_sanePic &Pic_param, int size, int rank) {

	dictionary * ini = NULL;
	string dummy;
	uint16_t parsed = 0;

#ifdef USE_MPI
	uint16_t mpi_parsed = 0;
#endif

	if (rank == 0) {
		// load dictionnary
		ini = iniparser_load(ini_name);

		if (ini == NULL) {
			output += "EE - cannot open file: " + StringOf(ini_name) + "\n";
			return INI_NOT_FOUND;
		}
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);

	if(size>1) {
		if(rank!=0)
			ini = dictionary_new(0);

		if(commit_dictionary(rank, ini)) {
			cout << "ERROR commit dictionary for rank : " << rank << ". EXITING..." << endl;
			return INI_NOT_FOUND;
		}
	}

#endif

	// default values :
	default_param(dir, Pos_param, Proc_param, Inv_param, Pic_param);

	//TODO sanePS should be out of here
	default_param_sanePS(PS_param);

	read_common(output, ini, dir);
	read_param_saneInv(output, ini, Inv_param);
	read_param_sanePos(output, ini, Pos_param);
	read_param_sanePic(output, ini, Pic_param);
	read_param_saneProc(output, ini, Proc_param);

	//TODO sanePS should be out of here (special case)
	read_param_sanePS(output, ini, PS_param);

	iniparser_freedict(ini);

	// Directory shall be available to all rank, so test for all..
	parsed |= check_common(output, dir, rank);

	// Now the ini file has been read, do the rest
	if (rank == 0) {
		parsed |= check_param_saneProc(output, dir, Proc_param);
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	parsed |=  MPI_Bcast(&( Proc_param.remove_polynomia), 1, MPI_INT,  0, MPI_COMM_WORLD);
	parsed |=  MPI_Bcast(&( Proc_param.highpass_filter), 1, MPI_INT,  0, MPI_COMM_WORLD);
#endif


	parsed |= check_param_sanePos(output, Pos_param);
	parsed |= check_param_saneInv(output, dir, Inv_param);

	if (rank==0){
		// Should probably be somewhere else...
		parsed |= check_param_sanePS(output, dir, PS_param);
	}

	parsed |= fill_samples_struct(output, samples_struct, dir, Inv_param, Proc_param, rank, size);

	// Test to know if all required files are present of not before doing the following... (based on the parsed value)
	parsed |= check_param_samples(output, samples_struct);

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&parsed, &mpi_parsed, 1, MPI_UNSIGNED_SHORT, MPI_BOR, MPI_COMM_WORLD);
	parsed = mpi_parsed;
#endif

	// Not needed at this stage, will be opened later...
	//	if (! parsed ){
	//		// Open the dirfile to read temporary files
	//		// Create it if it does not exist yet (case for saneFrameOrder and sanePre)
	//		for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++){
	//			string dirfile_basename = dir.tmp_dir + "dirfile/" + samples_struct.basevect[iframe];
	//
	//			if (! check_path(dummy, dirfile_basename , false, -1)) {
	//				samples_struct.dirfile_pointers[iframe] = gd_open((char *) dirfile_basename.c_str(),
	//						GD_RDWR | GD_CREAT | GD_VERBOSE | GD_UNENCODED);
	//				if (gd_error(samples_struct.dirfile_pointers[iframe]) != 0) {
	//					output += "error opening dirfile : " + dirfile_basename + "\n";
	//					return TMP_PATH_PROBLEM;
	//				}
	//			}
	//		}
	//	}

	return parsed;
}

// Structure print

void print_common(struct param_common dir) {

	cout << "Data Dir.        : " << dir.data_dir << endl;
	cout << "Input Dir.       : " << dir.input_dir << endl;
	cout << "Temp. Dir.       : " << dir.tmp_dir << endl;
	cout << "Output Dir.      : " << dir.output_dir << endl;
#ifdef USE_MPI
	cout << endl;
	cout << "Para. Mode       : " << dir.parallel_scheme << endl;
#endif
}

void print_param_sanePos(struct param_sanePos Pos_param) {

	cout << "Pixel Size       : " << prettyPrintAngle(Pos_param.pixdeg) << endl;

	if (Pos_param.flgdupl)
		cout << "Map Flags        : True" << endl;

	if (Pos_param.projgaps)
		cout << "Gap Filling      : Projected" << endl;
	else
		cout << "Gap Filling      : NOT projected" << endl;

	if (Pos_param.eq2gal)
		cout << "Converting to    : Galactic Coordinates" << endl;

	if (Pos_param.gal2eq)
		cout << "Converting to    : Equatorial Coordinates" << endl;

	cout << endl;
}

void print_param_saneProc(struct param_saneProc Proc_param) {

	if (Proc_param.fill_gap)
		cout << "Fill Gaps        : True\n";
	else
		cout << "Fill Gaps        : False\n";

	if (Proc_param.wisdom)
		cout << "Use Wisdom       : True\n";
	else
		cout << "Use Wisdom       : False\n";

	if (Proc_param.remove_linear)
		cout << "Simple Baseline  : will be removed (default)\n";
	else
		cout << "Simple Baseline  : will not be removed\n";

	if (Proc_param.CORRon)
		cout << "Correlations     : INCLUDED in the analysis" << endl;
	else
		cout << "Correlations     : NOT INCLUDED in the analysis" << endl;

	if (Proc_param.remove_polynomia)
		cout << "Poly. Order      : " << Proc_param.poly_order << endl;
	else
		cout << "Poly. Order      : None\n";

	if (Proc_param.napod > 0)
		cout << "# for Apodize    : " << Proc_param.napod << endl;

	cout << "HPF Freq.        : ";
	if (Proc_param.highpass_filter) {
		if (Proc_param.fhp_file == "")
			cout << Proc_param.fhp << " Hz" << endl;
		else
			cout << "[from file]" << endl;
	} else
		cout << "None" << endl;

	cout << "Sampling Freq.   : ";
	if ( Proc_param.fsamp_file == "" )
		cout << Proc_param.fsamp << " Hz" << endl;
	else
		cout << "[from file]" << endl;

}

void print_param_sanePic(struct param_sanePic Pic_param) {

	if (Pic_param.save_data)
		cout << "Write Iter. Maps : " << Pic_param.iterw << endl;
	else
		cout << "Write Iter. Maps : OFF \n";

	cout << "Max Iter.        : " << Pic_param.itermax << endl;

	cout << "Maps prefix      : " << Pic_param.map_prefix << endl;

	cout << endl;
}

void print_param_sanePS(struct param_sanePS PS_param) {

	if (PS_param.save_data)
		cout << "Save data.       : ON\n";
	else
		cout << "Save data.       : OFF\n";

	if (PS_param.signame != "")
		cout << "Removed map.     : " << PS_param.signame << endl;

	cout << "Noise comp.      : " << PS_param.ncomp << endl;

	cout << endl;
}

void print_param_saneInv(struct param_saneInv Inv_param) {

	cout << "Noise Dir.       : " << Inv_param.noise_dir << endl;

	cout << endl;
}

void parser_printOut(char * prog_name, struct param_common dir,
		struct samples samples_struct, struct param_sanePos Pos_param,
		struct param_saneProc Proc_param, struct param_sanePS PS_param,
		struct param_sanePic Pic_param,
		struct param_saneInv Inv_param) {

	string basename(prog_name);
	basename = FitsBasename(basename);
	int i;

	print_common(dir);

	i = basename.find("sanePos");
	if ((i >= 0) && (i < (int) basename.size())) {
		print_param_sanePos(Pos_param);
	}

	i = basename.find("sanePS");
	if ((i >= 0) && (i < (int) basename.size())) {
		print_param_sanePos(Pos_param);
		print_param_saneProc(Proc_param);
		print_param_sanePS(PS_param);
	}

	//	i=basename.find("saneProc");
	//	if((i>=0) && (i<(int)basename.size())){
	//		print_param_positions(Pos_param);
	//		print_param_process(Proc_param);
	//	}

	i = basename.find("sanePic");
	if ((i >= 0) && (i < (int) basename.size())) {
		print_param_sanePos(Pos_param);
		print_param_saneProc(Proc_param);
		print_param_sanePic(Pic_param);
	}

	i = basename.find("saneInv");
	if ((i >= 0) && (i < (int) basename.size())) {
		print_param_saneInv(Inv_param);
	}

	cout << "# of Scans       : " << StringOf(samples_struct.ntotscan) << endl;;

}


void export_param_sanePos(struct param_sanePos Pos_param, std::vector<string> &key, std::vector<string> &value, std::vector<string> &comment) {

	key.push_back("pixsize");
	value.push_back(StringOf(Pos_param.pixdeg));
	comment.push_back("size of pixels (degrees)");

	key.push_back("map_flagged_data");
	value.push_back(StringOf( Pos_param.flgdupl ?"True" :"False"));
	comment.push_back("flagged data put in a separated map (default is False)");

	key.push_back("file_format");
	value.push_back(StringOf(Pos_param.fileFormat));
	comment.push_back("SANEPIC 0, HIPE 1");

	key.push_back("project_gaps");
	value.push_back(StringOf(Pos_param.projgaps ? "True" :"False"));
	comment.push_back("gaps are projected to a pixel in the map, if so gap filling of noise only is performed iteratively. Default is False");

	key.push_back("mask_file");
	value.push_back(Pos_param.maskfile);
	comment.push_back("mask (fits file) : use to remove cross constraint due to strong sources");

	key.push_back("lon");
	value.push_back(StringOf(Pos_param.lon));
	comment.push_back("LON nominal of the input scan data");

	key.push_back("lat");
	value.push_back(StringOf(Pos_param.lat));
	comment.push_back("LAT nominal of the input scan data");

	key.push_back("proj_code");
	value.push_back(StringOf(Pos_param.projcode));
	comment.push_back("projection type, see wcslib projection code (default: TAN)");

	key.push_back("axis_type");
	value.push_back(StringOf(Pos_param.axistype));
	comment.push_back("axis type, EQ|GAL, (default EQ) ");

	key.push_back("eq2gal");
	value.push_back(StringOf(Pos_param.eq2gal ? "True" : "False"));
	comment.push_back("Convert J2000.0 EQ to GAL");

	key.push_back("gal2eq");
	value.push_back(StringOf(Pos_param.gal2eq ? "True" : "False"));
	comment.push_back("Convert GAL to J2000.0");

	key.push_back("restwav");
	value.push_back(StringOf(Pos_param.restwav));
	comment.push_back("rest wavelength in vacuo [m]");

	key.push_back("equinox");
	value.push_back(StringOf(Pos_param.equinox));
	comment.push_back("Equinox of celestial coordinate system");

	key.push_back("radesys");
	value.push_back(StringOf(Pos_param.radesys));
	comment.push_back("Coordinate reference frame for the RA and DEC");


}

void export_param_common(struct param_common dir, std::vector<string> &key, std::vector<string> &value, std::vector<string> &comment) {

	key.push_back("data_directory");
	value.push_back(dir.data_dir);
	comment.push_back("source data directory");

	key.push_back("input_directory");
	value.push_back(dir.input_dir);
	comment.push_back("input directory with all the configurations files");

	key.push_back("output_dir");
	value.push_back(dir.output_dir);
	comment.push_back("output directory");

	key.push_back("temp_dir");
	value.push_back(dir.tmp_dir);
	comment.push_back("temporary directory");

	key.push_back("fits_filelist");
	value.push_back(dir.fits_filelist);
	comment.push_back("file containing fits file names");

	key.push_back(".bolo");
	value.push_back(dir.bolo);
	comment.push_back("every scans have the same detector list");

	key.push_back("bolo_suffix");
	value.push_back(dir.bolo_suffix);
	comment.push_back("bolometers filelist suffix");

	key.push_back("para");
	value.push_back(dir.parallel_scheme);
	comment.push_back("parallel scheme");

}

void export_param_saneProc(struct param_saneProc Proc_param, std::vector<string> &key, std::vector<string> &value, std::vector<string> &comment) {

	key.push_back("fsamp");
	value.push_back(StringOf(Proc_param.fsamp));
	comment.push_back("detectors sampling frequency [Hz]");

	key.push_back("fsamp_file");
	value.push_back(Proc_param.fsamp_file);
	comment.push_back("filename containing the detectors sampling frequencies [Hz]");

	key.push_back("fhp");
	value.push_back(StringOf(Proc_param.fhp));
	comment.push_back("frequency of the high pass filter applied to the data [Hz]");

	key.push_back("fhp_file");
	value.push_back(Proc_param.fhp_file);
	comment.push_back("filename containing the frequencies of the high pass filter applied to the data [Hz]");

	key.push_back("apodize_Nsamples");
	value.push_back(StringOf(Proc_param.napod));
	comment.push_back("number of samples to apodize");

	key.push_back("fcut");
	value.push_back(StringOf(Proc_param.fcut));
	comment.push_back("frequency at which noise power spectra are thresholded");

	key.push_back("fcut_file");
	value.push_back(Proc_param.fcut_file);
	comment.push_back("filename containing the frequencies at which noise power spectra are thresholded");

	key.push_back("poly_order");
	value.push_back(StringOf(Proc_param.poly_order));
	comment.push_back("baseline polynomia order (default : 0 no baseline)");

	key.push_back("linear_baseline");
	value.push_back(StringOf(Proc_param.remove_linear ? "True" : "False"));
	comment.push_back("simple linear baseline removed from the data (False : default)");

	key.push_back("correlation");
	value.push_back(StringOf(Proc_param.CORRon ? "True" : "False"));
	comment.push_back("Set this keyword to False if correlations between detectors are not included in the analysis (True : default)");

	key.push_back("fill_gap");
	value.push_back(StringOf(Proc_param.fill_gap ? "True" : "False"));
	comment.push_back("Do we fill the gaps ?  (True : default)");

	key.push_back("wisdom");
	value.push_back(StringOf(Proc_param.wisdom ? "True" : "False"));
	comment.push_back("Do we use FFTW wisdom ?  (False: default)");


}

void export_param_saneInv(struct param_saneInv Inv_param, std::vector<string> &key, std::vector<string> &value, std::vector<string> &comment) {

	key.push_back("noise_dir");
	value.push_back(Inv_param.noise_dir);
	comment.push_back("cov matrix directory");

	key.push_back("cov_matrix");
	value.push_back(Inv_param.cov_matrix);
	comment.push_back("this file contains the matrix you want to invert");

	key.push_back("cov_matrix_suffix");
	value.push_back(Inv_param.cov_matrix_suffix);
	comment.push_back("this file contains the matrix you want to invert");

}

void export_param_sanePS(struct param_sanePS PS_param, std::vector<string> &key, std::vector<string> &value, std::vector<string> &comment) {

	key.push_back("MixingMatrix_suffix");
	value.push_back(PS_param.mix_suffix);
	comment.push_back("Mixing matrix files suffix");

	key.push_back("MixingMatrix");
	value.push_back(PS_param.mix);
	comment.push_back("the MixingMatrix file common to all scans");

	key.push_back("ncomp");
	value.push_back(StringOf(PS_param.ncomp));
	comment.push_back("number of component(s) to estimate (default : 1)");

	key.push_back("niter");
	value.push_back(StringOf(PS_param.niter));
	comment.push_back("number of iteration in the expectation minimization loop (default : 500"
			")");

	key.push_back("map_file");
	value.push_back(PS_param.signame);
	comment.push_back("map substracted from the data (sanePos needed)");

	key.push_back("ell");
	value.push_back(PS_param.ell);
	comment.push_back("the ell file common to all scans");

	key.push_back("ell_suffix");
	value.push_back(PS_param.ell_suffix);
	comment.push_back("ell files suffix");

	key.push_back("save_data");
	value.push_back(StringOf(PS_param.save_data));
	comment.push_back("set save_data to 1 if you wish to save the processing session after each sanePS step");

}

void export_param_sanePic(struct param_sanePic Pic_param, std::vector<string> &key, std::vector<string> &value, std::vector<string> &comment) {

	key.push_back("iterW");
	value.push_back(StringOf(Pic_param.iterw));
	comment.push_back("Write temporary map files on disk every iterW number of loop");

	key.push_back("tolerance");
	value.push_back(StringOf(Pic_param.tolerance));
	comment.push_back("Tolerance to reach the end of the loop");

	key.push_back("subtolerance");
	value.push_back(StringOf(Pic_param.subtolerance));
	comment.push_back("Tolerance to reach the first map");

	key.push_back("iterMAX");
	value.push_back(StringOf(Pic_param.itermax));
	comment.push_back("Maximum number of conjugate gradient loops (default: 2000)");

	key.push_back("map_prefix");
	value.push_back(Pic_param.map_prefix);
	comment.push_back("prefix for the fits file (default : optimMap");

}

void export_param_saneCheck(struct param_saneCheck Check_param, std::vector<string> &key, std::vector<string> &value, std::vector<string> &comment) {

	key.push_back("check_NAN");
	value.push_back(StringOf(Check_param.checkNAN ? "True" : "False"));
	comment.push_back("check whether there are NANs in every tables and flag them if needed");

	key.push_back("check_time_gaps");
	value.push_back(StringOf(Check_param.checktime ? "True" : "False"));
	comment.push_back("check whether there are gaps in time table and fill those gaps with flagged data to ensure continuity");

	key.push_back("check_flag");
	value.push_back(StringOf(Check_param.checkflag ? "True" : "False"));
	comment.push_back("check whether there are detectors that are fully or more than 80% flagged");

	key.push_back("check_bolo_gain");
	value.push_back(StringOf(Check_param.checkGain ? "True" : "False"));
	comment.push_back("Compute detector gain and print to screen");

}

string rebuild_ini(struct param_common dir, struct param_saneProc proc_param, struct param_sanePos pos_param,
		struct samples samples_struct, struct param_sanePS PS_param,
		struct param_sanePic Pic_param, struct param_saneInv Inv_param) {

	string text;

	// Generates ini file model
	text = "# Sanepic ini file :\n";
	text += "# Edit the lines to modify the options !\n";
	text += "# Some options are required : Check the comments !\n";
	text += "\n";

	std::vector<string> key_vect;
	std::vector<string> val_vect;
	std::vector<string> com_vect;

	text += "\n[common]\n";
	export_param_common(dir, key_vect, val_vect, com_vect);
	for (unsigned long ii=0; ii < key_vect.size(); ii++)
		text += key_vect[ii] + " = " + val_vect[ii] + " ; " + com_vect[ii] + "\n";
	key_vect.clear(); val_vect.clear(); com_vect.clear();

	text += "\n[sanePos]\n";
	export_param_sanePos(pos_param,   key_vect, val_vect, com_vect);
	for (unsigned long ii=0; ii < key_vect.size(); ii++)
		text += key_vect[ii] + " = " + val_vect[ii] + " ; " + com_vect[ii] + "\n";
	key_vect.clear(); val_vect.clear(); com_vect.clear();

	text += "\n[saneProc]\n";
	export_param_saneProc(proc_param, key_vect, val_vect, com_vect);
	for (unsigned long ii=0; ii < key_vect.size(); ii++)
		text += key_vect[ii] + " = " + val_vect[ii] + " ; " + com_vect[ii] + "\n";
	key_vect.clear(); val_vect.clear(); com_vect.clear();

	text += "\n[saneInv]\n";
	export_param_saneInv(Inv_param,   key_vect, val_vect, com_vect);
	for (unsigned long ii=0; ii < key_vect.size(); ii++)
		text += key_vect[ii] + " = " + val_vect[ii] + " ; " + com_vect[ii] + "\n";
	key_vect.clear(); val_vect.clear(); com_vect.clear();

	text += "\n[sanePic]\n";
	export_param_sanePic(Pic_param,   key_vect, val_vect, com_vect);
	for (unsigned long ii=0; ii < key_vect.size(); ii++)
		text += key_vect[ii] + " = " + val_vect[ii] + " ; " + com_vect[ii] + "\n";
	key_vect.clear(); val_vect.clear(); com_vect.clear();

	text += "\n[sanePS]\n";
	export_param_sanePS(PS_param,     key_vect, val_vect, com_vect);
	for (unsigned long ii=0; ii < key_vect.size(); ii++)
		text += key_vect[ii] + " = " + val_vect[ii] + " ; " + com_vect[ii] + "\n";
	key_vect.clear(); val_vect.clear(); com_vect.clear();

	return text;

}
