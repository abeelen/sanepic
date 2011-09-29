#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <typeinfo>
#include <sys/types.h>  // For stat()
#include <sys/stat.h>   // For stat()
#include <typeinfo>

#include <wordexp.h>

extern "C" {
#include "iniparser.h"
#include "dictionary.h"
#include "getdata.h"
}

#include "utilities.h"
#include "inputFileIO.h"
#include "parser_functions.h"
#include "mpi_architecture_builder.h"
#include "crc.h"

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

	s = iniparser_getstring(ini, "common:bolo_global_file", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "common:bolos_global_file : default value [" + StringOf(
				common.bolo_global_filename) + "]\n";
	else
		common.bolo_global_filename = StringOf(s);

#ifdef DEBUG
	output += output2;
#endif

}

void read_param_sanePos(string &output, dictionary *ini,
		struct param_sanePos &Pos_param) {

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


#ifdef DEBUG
	output += output2;
#endif

}

void read_param_saneProc(string &output, dictionary *ini,
		struct param_saneProc &Proc_param) {

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

	i = iniparser_getboolean(ini, "saneProc:fill_gap", -1);
	if (i == -1)
		output2 += "saneProc:fill_gap : default value [" + StringOf(
				Proc_param.fill_gap) + "]\n";
	else
		Proc_param.fill_gap = (bool) i;

	d = iniparser_getdouble(ini, (char*) "saneProc:sampling_frequency", -1.0);
	if (d == -1.0)
		output2 += "saneProc:sampling_frequency: default value [" + StringOf(
				Proc_param.fsamp) + "]\n";
	else
		Proc_param.fsamp = d;

	d = iniparser_getdouble(ini, (char*) "saneProc:filter_frequency", -1.0);
	if (d == -1.0)
		output2 += "saneProc:filter_frequency: default value [" + StringOf(
				Proc_param.f_lp) + "]\n";
	else
		Proc_param.f_lp = d;

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

	i = iniparser_getint(ini, "saneProc:poly_order", -1);
	if (i == -1)
		output2 += "saneProc:poly_order: default value [" + StringOf(
				Proc_param.poly_order) + "]\n";
	else
		Proc_param.poly_order = i;

	if (Proc_param.poly_order >= 0)
		Proc_param.remove_polynomia = 1;
	else
		Proc_param.remove_polynomia = 0;

	if (Proc_param.f_lp > 0)
		Proc_param.highpass_filter = 1;
	else
		Proc_param.highpass_filter = 0;

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

void read_param_saneInv(std::string &output, dictionary *ini,
		struct param_saneInv &Inv_param) {

	char *s;
	string output2 = "";

	s = iniparser_getstring(ini, "saneInv:cov_matrix_file", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "saneInv:cov_matrix_file : default value [" + StringOf(
				Inv_param.cov_matrix_file) + "]\n";
	else
		Inv_param.cov_matrix_file = StringOf(s);

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

void read_param_sanePS(std::string &output, dictionary *ini,
		struct param_sanePS &PS_param) {

	int i;
	char *s;
	string output2 = "";

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

	s = iniparser_getstring(ini, "sanePS:ell_global_file", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePS:ell_global_file : default value [" + StringOf(
				PS_param.ell_global_file) + "]\n";
	else
		PS_param.ell_global_file = StringOf(s);

	s = iniparser_getstring(ini, "sanePS:MixingMatrix_global_file",
			(char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePS:MixingMatrix_global_file : default value ["
				+ StringOf(PS_param.mix_global_file) + "]\n";
	else
		PS_param.mix_global_file = StringOf(s);

	i = iniparser_getboolean(ini, "sanePS:save_data", -1);
	if (i == -1)
		output2 += "sanePS:save_data: default value [" + StringOf(
				PS_param.save_data) + "]\n";
	else
		PS_param.save_data = (bool) i;

	//TODO: Ugly turnaround until sanePS is released;
	s = iniparser_getstring(ini, "saneInv:cov_matrix_file", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "saneInv:cov_matrix_file : default value [" + StringOf(
				PS_param.cov_matrix_file) + "]\n";
	else
		PS_param.cov_matrix_file = StringOf(s);

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
	char *s;
	string output2 = "";

	i = iniparser_getint(ini, "sanePic:iterW", -1);
	if (i == -1)
		output2 += "sanePic:iterW : default value [" + StringOf(
				Pic_param.iterw) + "]\n";
	else
		Pic_param.iterw = i;

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

int check_path(string &output, string strPath, bool create) {

	if (access(strPath.c_str(), 0) == 0) {
		struct stat status;
		stat(strPath.c_str(), &status);

		if (status.st_mode & S_IFDIR) {
			return 0;
		} else {
			output += "EE - " + strPath + " is a file.\n";
			return 1;
		}
	} else {
		if (! create) {
			output += "EE - " + strPath + " does not exist\n";
			return 1;
		}

		int status;
		status = mkdir(strPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (status == 0)
			output += "WW - " + strPath + " created\n";
		else {
			output += "EE - " + strPath + " failed to create\n";
			return 1;
		}
	}
	return 0;
}

int compute_dirfile_format_file(std::string tmp_dir,
		struct samples samples_param, int format) {

	string filedir = tmp_dir + "dirfile";

	DIRFILE *D, *H, *F, *I, *I2, *J, *K, *S, *R, *R2;

	// create folders
	D = gd_open((char *) filedir.c_str(),
			GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED);

	for (long iframe = 0; iframe < samples_param.ntotscan; iframe++) {

		string scan_name = samples_param.basevect[iframe];
		string scan_folder = filedir + "/" + scan_name;
		string fdata = filedir + "/" + scan_name + "/fData";
		string index_path = filedir + "/" + scan_name + "/Indexes";
		string data = filedir + "/" + scan_name + "/data";
		string flag_dir = filedir + "/" + scan_name + "/flag";
		string LON = filedir + "/" + scan_name + "/LON";
		string LAT = filedir + "/" + scan_name + "/LAT";
		string noise_path = filedir + "/" + scan_name + "/Noise_data";
		string ell_path = filedir + "/" + scan_name + "/Noise_data/ell";

		// create folders
		S = gd_open((char *) scan_folder.c_str(),
				GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED);
		H = gd_open((char *) index_path.c_str(),
				GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED);
		F = gd_open((char *) fdata.c_str(),
				GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED);
		J = gd_open((char *) data.c_str(),
				GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED);
		K = gd_open((char *) flag_dir.c_str(),
				GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED);
		I = gd_open((char *) noise_path.c_str(),
				GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED);
		I2 = gd_open((char *) ell_path.c_str(),
				GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED);
		if (format == 1) {

			// create folders
			R = gd_open((char *) LON.c_str(),
					GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED);
			R2 = gd_open((char *) LAT.c_str(),
					GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED);

			// close subdirfiles
			gd_close(R);
			gd_close(R2);
		}

		// close subdirfiles
		gd_close(H);
		gd_close(F);
		gd_close(J);
		gd_close(K);
		gd_close(S);
		gd_close(I);
		gd_close(I2);

		// include subdir and create format files
		gd_include(D, (char *) (scan_name + "/format").c_str(), 0,
				GD_RDWR | GD_CREAT | GD_UNENCODED);

		S = gd_open((char *) scan_folder.c_str(),
				GD_RDWR | GD_VERBOSE | GD_UNENCODED);

		gd_include(S, (char *) ("Indexes/format"), 0,
				GD_RDWR | GD_CREAT | GD_UNENCODED);
		gd_include(S, (char *) ("fData/format"), 0,
				GD_RDWR | GD_CREAT | GD_UNENCODED);
		gd_include(S, (char *) ("flag/format"), 0,
				GD_RDWR | GD_CREAT | GD_UNENCODED);
		gd_include(S, (char *) ("data/format"), 0,
				GD_RDWR | GD_CREAT | GD_UNENCODED);
		gd_include(S, (char *) ("Noise_data/format"), 0,
				GD_RDWR | GD_CREAT | GD_UNENCODED);
		gd_include(S, (char *) ("Noise_data/ell/format"), 0,
				GD_RDWR | GD_CREAT | GD_UNENCODED);
		if (format == 1) {
			gd_include(S, (char *) ("LON/format"), 0,
					GD_RDWR | GD_CREAT | GD_UNENCODED);
			gd_include(S, (char *) ("LAT/format"), 0,
					GD_RDWR | GD_CREAT | GD_UNENCODED);
		}

		gd_flush(S, NULL);

		gd_close(S);
	}

	// close dirfile
	gd_close(D);

	return 0;
}

int cleanup_dirfile_sanePos(std::string tmp_dir, struct samples samples_param,
		std::vector<std::vector<std::string> > bolo_vect) {

	std::vector<string> det_vect;

	for (long iframe = 0; iframe < samples_param.ntotscan; iframe++) {

		det_vect = bolo_vect[iframe];

		string scan_name = samples_param.basevect[iframe];
		string index_path = tmp_dir + "dirfile/" + scan_name + "/Indexes";

		DIRFILE *S = gd_open((char *) index_path.c_str(),
				GD_RDWR | GD_TRUNC | GD_VERBOSE | GD_UNENCODED);

		// then generate binaries and fill format file
		for (long idet = 0; idet < (long) det_vect.size(); idet++) {
			string outfile = scan_name + "_" + det_vect[idet];
			//configure dirfile field
			gd_entry_t E;
			E.field = (char*) outfile.c_str();
			E.field_type = GD_RAW_ENTRY;
			E.fragment_index = 0;
			E.spf = 1;
			E.data_type = GD_INT64;
			E.scalar[0] = NULL;

			// add to the dirfile
			gd_add(S, &E);
			gd_flush(S, E.field);
		}

		if (gd_close(S))
			cout << "error closing " << index_path << "-> memory leaks ..."
			<< endl;

	}
	return 0;
}

int cleanup_dirfile_saneInv(std::string tmp_dir, struct samples samples_param,
		long nframe, string noise_suffix,
		std::vector<std::vector<std::string> > bolo_vect) {

	std::vector<string> det_vect;

	for (long ii = 0; ii < nframe; ii++) {

		det_vect = bolo_vect[ii];

		string base_name = samples_param.basevect[ii];
		string noise_path = tmp_dir + "dirfile/" + base_name + "/Noise_data";
		string ell_path = noise_path + "/ell";

		DIRFILE *S = gd_open((char *) noise_path.c_str(),
				GD_RDWR | GD_TRUNC | GD_UNENCODED);
		DIRFILE *D = gd_open((char *) ell_path.c_str(),
				GD_RDWR | GD_TRUNC | GD_UNENCODED);

		string suffix = base_name + noise_suffix; // base_name instead of noisevect[ii]

		for (int idet = 0; idet < (long) det_vect.size(); idet++) {

			// ell binary filename
			string outfile = det_vect[idet] + "_" + suffix + "_ell";

			// configure dirfile field for ell
			gd_entry_t E;
			E.field = (char*) outfile.c_str();
			E.field_type = GD_RAW_ENTRY;
			E.fragment_index = 0;
			E.spf = 1;
			E.data_type = GD_DOUBLE;
			E.scalar[0] = NULL;

			// add to the dirfile
			gd_add(D, &E);

			// spectra filename
			outfile = det_vect[idet] + "_" + suffix;

			// set field information for spectra
			E.field = (char*) outfile.c_str();

			// add to the dirfile
			gd_add(S, &E);
			gd_flush(S, NULL);
		}

		if (gd_close(S))
			cout << "error closing " << noise_path << "-> memory leaks ..."
			<< endl;
		if (gd_close(D))
			cout << "error closing " << ell_path << "-> memory leaks ..."
			<< endl;

	}

	return 0;
}

int cleanup_dirfile_fdata(std::string tmp_dir, struct samples samples_param,
		std::vector<std::vector<std::string> > bolo_vect) {

	for (long iframe = 0; iframe < samples_param.ntotscan; iframe++) {

		std::vector<string> det_vect = bolo_vect[iframe];

		//get fourier transform dirfile names !
		string scan_name = samples_param.basevect[iframe];
		string fdata_path = tmp_dir + "dirfile/" + scan_name + "/fData";

		// clean up the dirfiles with TRUNC option
		DIRFILE *S = gd_open((char *) fdata_path.c_str(),
				GD_RDWR | GD_TRUNC | GD_VERBOSE | GD_UNENCODED);

		// then generate binaries and fill format file
		string prefixe[2] = { "fdata_", "fPs_" };
		for (long ip = 0; ip < 2; ip++)
			for (long idet = 0; idet < (long) det_vect.size(); idet++) {
				string outfile = prefixe[ip] + scan_name + "_" + det_vect[idet];
				//configure dirfile field
				gd_entry_t E;
				E.field = (char*) outfile.c_str();
				E.field_type = GD_RAW_ENTRY;
				E.fragment_index = 0;
				E.spf = 1;
				E.data_type = GD_COMPLEX128;
				E.scalar[0] = NULL;

				// add to the dirfile
				gd_add(S, &E);
				gd_flush(S, E.field);
			}

		gd_close(S);

		// check sizes in Indexes, data and flag format
		string indexes_path = tmp_dir + "dirfile/" + scan_name + "/Indexes";
		string data_path = tmp_dir + "dirfile/" + scan_name + "/data";
		DIRFILE *I = gd_open((char *) indexes_path.c_str(),
				GD_RDWR | GD_VERBOSE | GD_UNENCODED);
		DIRFILE *D = gd_open((char *) data_path.c_str(),
				GD_RDWR | GD_VERBOSE | GD_UNENCODED);

		long nframeI = gd_nframes(I);
		long nframeD = gd_nframes(I);

		long ns = samples_param.nsamples[iframe];
		gd_close(I);
		gd_close(D);

		if ((nframeI != ns) || (nframeD != ns)) {
			cout << "Error... Dirfile data or Indexes has incorrect size !!\n";
			cout << indexes_path << " : " << nframeI << " (vs ns= " << ns <<")" << endl;
			cout << data_path    << " : " << nframeD << " (vs ns= " << ns <<")" << endl;
			return 1;
		}

	}

	return 0;
}

uint16_t check_common(string &output, struct param_common dir) {

	if ((dir.bolo_global_filename == "") && (dir.bolo_suffix == "")) {
		output += "EE - You must mention one of those parameters :\n";
		output += "     param_common:bolo_suffix or param_common:bolo_blobal_file\n";
		return BOLOFILE_NOT_FOUND;
	}
	if (check_path(output, dir.data_dir, false))
		return DATA_INPUT_PATHS_PROBLEM;
	if (check_path(output, dir.input_dir, false))
		return DATA_INPUT_PATHS_PROBLEM;
	if (check_path(output, dir.output_dir, true))
		return OUPUT_PATH_PROBLEM;
	if (check_path(output, dir.tmp_dir, true))
		return TMP_PATH_PROBLEM;

	return 0;
}

uint16_t check_param_sanePos(string &output, struct param_sanePos Pos_param) {

	if (Pos_param.pixdeg <= 0 && Pos_param.maskfile == "") {
		output += "EE - Pixsize cannot be negative ! \n";
		return PIXDEG_WRONG_VALUE;
	}
	if ((Pos_param.fileFormat != 0) && (Pos_param.fileFormat != 1)) {
		output += "EE - Fileformat must be 0 (SANEPIC) or 1 (HIPE) \n";
		return FILEFORMAT_NOT_FOUND;
	}

	// Force axis type to GAL if converting from EQ to GAL
	if (Pos_param.eq2gal)
		Pos_param.axistype = "GAL";

	// Force axis type to EQ if converting from GAL to EQ
	if (Pos_param.gal2eq)
		Pos_param.axistype = "EQ";


	return 0;
}

uint16_t check_param_saneProc(string &output, struct param_saneProc Proc_param) {

	if (Proc_param.napod < 0) {
		output
		+= "EE - You must choose a positive number of samples to apodize\n";
		return NAPOD_WRONG_VALUE;
	}
	if (Proc_param.fsamp <= 0.0) {
		output += "EE - Sampling_frequency cannot be negative or 0 ! \n";
		return FSAMP_WRONG_VALUE;
	}

	return 0;
}

uint16_t check_param_sanePS(string &output, struct param_sanePS PS_param) {

	if (PS_param.ncomp <= 0) {
		output += "EE - Number of component ncomp cannot be negative or zero ! \n";
		return NCOMP_WRONG_VALUE;
	}
	if ((PS_param.ell_global_file == "") && (PS_param.ell_suffix == "")) {
		output += "EE - You must mention one of those parameters :\n";
		output += "     sanePS:ell_global_file or sanePS:ell_suffix\n";
		return ELL_FILE_NOT_FOUND;
	}
	if ((PS_param.mix_global_file == "") && (PS_param.mix_suffix == "")) {
		output += "EE - You must mention one of those parameters :\n";
		output += "     sanePS:mix_global_file or sanePS:mix_suffix\n";
		return MIX_FILE_NOT_FOUND;
	}

	return 0;
}

uint16_t check_param_saneInv(string &output,
		struct param_saneInv Inv_param) {

	if (check_path(output, Inv_param.noise_dir, false))
		return SANEINV_INPUT_ERROR;

	if ((Inv_param.cov_matrix_file == "")
			&& (Inv_param.cov_matrix_suffix == "")) {
		output += "EE - You must mention one of those parameters :\n";
		output += "     saneInv:cov_matrix_suffix or saneInv:cov_matrix_global_file\n";
		return SANEINV_INPUT_ERROR;
	}

	return 0;
}

void default_param(struct param_common &dir, struct samples &samples_param,
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
	dir.bolo_global_filename = "";
	dir.bolo_suffix = ".bolo";

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
	//	Proc_param.f_lp = 0.0; // low pass filter frequency

	Proc_param.remove_linear = false;
	Proc_param.fill_gap = true;
	Proc_param.CORRon = true;
	Proc_param.remove_polynomia = true;
	Proc_param.highpass_filter = false;

	Proc_param.fcut_file = "";
	Proc_param.napod = 100;
	Proc_param.poly_order = 1;
	Proc_param.fsamp = 0.0;
	Proc_param.f_lp = 0.0;

}

void default_param_sanePS(struct param_sanePS &PS_param) {

	PS_param.ell_suffix = ".ell";
	PS_param.mix_suffix = ".mix";
	PS_param.ell_global_file = "";
	PS_param.mix_global_file = "";
	PS_param.signame = "";
	PS_param.ncomp = 1;
	PS_param.save_data = 1;

	//TODO: Ugly turnaround until sanePS is released;

	PS_param.cov_matrix_file = "";
	PS_param.cov_matrix_suffix = "_ps.fits";

}

void default_param_saneInv(struct param_saneInv &Inv_param) {

	Inv_param.cov_matrix_file = "";
	Inv_param.cov_matrix_suffix = "_psd.fits";
	Inv_param.noise_dir = "./";
}

void default_param_sanePic(struct param_sanePic &Pic_param) {

	Pic_param.iterw = 0;
	Pic_param.itermax = 2000;
	Pic_param.save_data = 0;
	Pic_param.map_prefix = "optimMap";
}

void fill_sanePS_struct(struct param_sanePS &PS_param,
		struct samples &samples_param, struct param_common &dir) {

	for (long iframe = 0; iframe < samples_param.ntotscan; iframe++) {
		if (PS_param.mix_global_file != "")
			samples_param.mix_names.push_back(
					dir.input_dir + PS_param.mix_global_file);
		else
			samples_param.mix_names.push_back(
					dir.input_dir + FitsBasename(
							samples_param.fitsvect[iframe])
							+ PS_param.mix_suffix);

		if (PS_param.ell_global_file != "")
			samples_param.ell_names.push_back(
					dir.input_dir + PS_param.ell_global_file);
		else
			samples_param.ell_names.push_back(
					dir.input_dir + FitsBasename(
							samples_param.fitsvect[iframe])
							+ PS_param.ell_suffix);
	}

}

uint16_t fill_samples_param(string &output, struct samples &samples_param,
		struct param_common &dir, struct param_saneInv &Inv_param,
		string fcut_file) {

	string filename;
	filename = dir.input_dir + dir.fits_filelist;
	if (read_fits_list(output, filename, samples_param) != 0)
		return 0x4000;

	samples_param.ntotscan = (samples_param.fitsvect).size();

	// Fill basevect
	for (long iframe = 0; iframe < samples_param.ntotscan; iframe++) {
		samples_param.basevect.push_back(
				dirfile_Basename(samples_param.fitsvect[iframe]));
	}
	// Fill bolovect
	for (long iframe = 0; iframe < samples_param.ntotscan; iframe++) {
		if (dir.bolo_global_filename != "")
			samples_param.bolovect.push_back(dir.bolo_global_filename);
		else
			samples_param.bolovect.push_back(
					FitsBasename(samples_param.fitsvect[iframe])
					+ dir.bolo_suffix);

	}

	// Fill noisevect
	for (long iframe = 0; iframe < samples_param.ntotscan; iframe++) {
		if ((Inv_param.cov_matrix_file != ""))
			samples_param.noisevect.push_back(Inv_param.cov_matrix_file);
		else
			samples_param.noisevect.push_back(
					FitsBasename(samples_param.fitsvect[iframe])
					+ Inv_param.cov_matrix_suffix);
	}

	if (FitsBasename(fcut_file).size() == 0) {
		output
		+= "Warning ! fcut_file filename is missing in the ini file ...\n";
		return 0x8000;
	}

	// Read the fcut file
	std::vector<string> dummy2;
	if (read_strings(fcut_file, dummy2))
		return 0x8000;

	if (((int) dummy2.size()) == 0 || ((int) dummy2.size() != 1
			&& ((int) dummy2.size() != samples_param.ntotscan))) {
		output
		+= "You must provide at least one number of noise cut frequency (or one per scan) in fcut_file !\n";
		return 0x8000;
	}

	for (int iframe_dummy = 0; iframe_dummy < (int) dummy2.size(); iframe_dummy++)
		samples_param.fcut.push_back(atof(dummy2[iframe_dummy].c_str()));

	// If There is only one value, use it everywhere
	if ((dummy2.size() == 1) && (samples_param.ntotscan > 1)) {
		// if only one fcut, extend to all scans
		samples_param.fcut.resize(samples_param.ntotscan,
				samples_param.fcut[0]);
	}

	return 0;

}

int get_noise_bin_sizes(std::string tmp_dir, struct samples &samples_param, int rank) {

	long nframe_long;

	for (long ii = 0; ii < samples_param.ntotscan; ii++) {

		if (rank == 0) {
			string scan_name = samples_param.basevect[ii];
			// dirfile path
			string filedir = tmp_dir + "dirfile/" + scan_name
					+ "/Noise_data/ell/";

			// open dirfile
			DIRFILE* H = gd_open((char *) filedir.c_str(),
					GD_RDWR | GD_VERBOSE | GD_UNENCODED);
			unsigned int nframe = gd_nframes(H);

			// close dirfile
			if (gd_close(H)) {
				cout << "Dirfile gd_close error in get_noise_bin_sizes for : "
						<< filedir << endl;
				return 1;
			}

			nframe_long = (long) (nframe - 1);
		}

#ifdef USE_MPI
		MPI_Bcast(&nframe_long,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
#endif

		// get nbins value
		samples_param.nbins.push_back(nframe_long);

		if (rank == 0) {

			string scan_name = samples_param.basevect[ii];

			// get ndet value
			string filedir = tmp_dir + "dirfile/" + scan_name + "/Noise_data/";
			DIRFILE* H = gd_open((char *) filedir.c_str(),
					GD_RDWR | GD_VERBOSE | GD_UNENCODED);
			unsigned int nframe = gd_nframes(H);

			// close dirfile
			if (gd_close(H)) {
				cout << "Dirfile gd_close error in get_noise_bin_sizes for : "
						<< filedir << endl;
				return 1;
			}
			nframe_long = (long) nframe;
		}

#ifdef USE_MPI
		MPI_Bcast(&nframe_long,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
#endif

		// compute ndet considering entry size and nbins
		samples_param.ndet.push_back(nframe_long / samples_param.nbins[ii]);

	}

	return 0;
}

int channel_list_to_vect_list(struct samples samples_param,
		std::vector<std::vector<std::string> > &bolo_vect, int rank) {

#ifdef USE_MPI
	long size_max;
#endif

	long ndet;
	std::vector<string> det_vect;

	for (long iframe = 0; iframe < samples_param.ntotscan; iframe++) {

		if (rank == 0) {

			string output_read = "";
			if (read_channel_list(output_read, samples_param.bolovect[iframe],
					det_vect)) {
				cout << output_read << endl;
#ifdef USE_MPI
				MPI_Abort(MPI_COMM_WORLD, 1);
#endif
				return 1;
			}
#ifdef USE_MPI
			compute_bololist_size(det_vect, size_max);
#endif
			ndet = (long) det_vect.size();
		}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD); // other procs wait untill rank 0 has created dirfile architecture.
		MPI_Bcast(&ndet, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&size_max, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

		char* temp=new char[size_max+1];

		for(long ii=0; ii< ndet; ii++) {
			fill(temp,temp+size_max+1,'\0');

			if(rank==0)
				strcpy (temp, det_vect[ii].c_str());

			MPI_Barrier(MPI_COMM_WORLD); // other procs wait untill rank 0 has created dirfile architecture.
			MPI_Bcast(temp, size_max+1, MPI_CHAR, 0, MPI_COMM_WORLD);

			if(rank!=0)
				det_vect.push_back(temp);
		}

		delete [] temp;
#endif

		bolo_vect.push_back(det_vect);
		det_vect.clear();
	}

	return 0;
}

#ifdef USE_MPI

long compute_bololist_size(std::vector<std::string> str_vect, long &size_max)
{
	long size_buff=0;
	size_max=0;

	for(long idet=0; idet < (long)str_vect.size(); idet++) {
		string temp = str_vect[idet];
		size_buff +=temp.size();

		if((long)temp.size()>size_max)
			size_max=temp.size();
	}

	return size_buff;
}

int commit_dictionary(int rank, dictionary *dict) {

	int n;
	int size;

	char *key_buff=NULL;
	char *val_buff=NULL;

	if(rank==0) {
		n = dict->n;
		size = dict->size;

		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
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

		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
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
uint16_t parser_function(char * ini_name, std::string &output,
		struct param_common &dir, struct samples &samples_param,
		struct param_sanePos &Pos_param, struct param_saneProc &Proc_param,
		struct param_sanePS &PS_param, struct param_saneInv &Inv_param,
		struct param_sanePic &Pic_param, int size, int rank) {

	dictionary * ini = NULL;
	string filename;
	uint16_t parsed = 0x0000;

	if (rank == 0) {
		// load dictionnary
		ini = iniparser_load(ini_name);

		if (ini == NULL) {
			fprintf(stderr, "cannot parse file: %s\n", ini_name);
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
	default_param(dir, samples_param, Pos_param, Proc_param, Inv_param,
			Pic_param);

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

	// Fill fitsvec, noisevect, scans_index with values read from the 'str' filename
	filename = dir.input_dir + Proc_param.fcut_file;
	parsed += fill_samples_param(output, samples_param, dir, Inv_param,
			filename);

	for (int iframe = 0; iframe < (int) ((samples_param.fitsvect).size()); iframe++) {
		samples_param.bolovect[iframe] = dir.input_dir
				+ samples_param.bolovect[iframe]; // better for bolovect cause you dont need to handle path in every function call !
	}
	// TODO: this should be done elsewhere as the input files might be missing...
	// Store scan sizes so that we dont need to read it again and again in the loops !
	readFrames(dir.data_dir, samples_param.fitsvect, samples_param.nsamples);


	// Now the ini file has been read, do the rest
	if (rank == 0)
		parsed += check_common(output, dir);

#ifdef PARA_FRAME
	MPI_Barrier(MPI_COMM_WORLD);
#endif


	parsed += check_param_sanePos(output, Pos_param);

	parsed += check_param_saneProc(output, Proc_param);

	parsed += check_param_saneInv(output, Inv_param);

	parsed += check_param_sanePS(output, PS_param);

	return parsed;
}

void print_common(struct param_common dir) {

	cout << "Data Dir.        : " << dir.data_dir << endl;
	cout << "Input Dir.       : " << dir.input_dir << endl;
	cout << "Temp. Dir.       : " << dir.tmp_dir << endl;
	cout << "Output Dir.      : " << dir.output_dir << endl;
}

void print_param_sanePos(struct param_sanePos Pos_param) {

	cout << "Pixel Size       : " << setprecision(14) << Pos_param.pixdeg
			<< " deg\n";

	if (Pos_param.flgdupl)
		cout << "Map Flags        : True" << endl;

	if (Pos_param.projgaps)
		cout << "Gap Filling      : PROJECTED" << endl;
	else
		cout << "Gap Filling      : NOT projected (default)" << endl;

	if (Pos_param.eq2gal)
		cout << "Converting to    : Galactic Coordinates" << endl;

	if (Pos_param.gal2eq)
		cout << "Converting to    : Equatorial Coordinates" << endl;

	cout << endl;
}

void print_param_process(struct param_saneProc Proc_param) {

	if (Proc_param.fill_gap)
		cout << "Fill Gaps        : True\n";
	else
		cout << "Fill Gaps        : False\n";

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

	if (Proc_param.highpass_filter)
		cout << "HPF Freq.        : " << Proc_param.f_lp << " Hz" << endl;
	else
		cout << "HPF Freq.        : None" << endl;

	cout << "Sampling Freq.   : " << Proc_param.fsamp << " Hz\n";

	cout << endl;
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
		struct samples samples_param, struct param_sanePos Pos_param,
		struct param_saneProc Proc_param, struct param_sanePS PS_param,
		struct param_sanePic Pic_param,
		struct param_saneInv Inv_param) {

	string basename(prog_name);
	basename = FitsBasename(basename);
	int i;

	print_common(dir);
	cout << endl;

	i = basename.find("sanePos");
	if ((i >= 0) && (i < (int) basename.size())) {
		print_param_sanePos(Pos_param);
	}

	i = basename.find("sanePS");
	if ((i >= 0) && (i < (int) basename.size())) {
		print_param_sanePos(Pos_param);
		print_param_process(Proc_param);
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
		print_param_process(Proc_param);
		print_param_sanePic(Pic_param);
	}

	i = basename.find("saneInv");
	if ((i >= 0) && (i < (int) basename.size())) {
		print_param_saneInv(Inv_param);
	}

	printf("# of Scans       : %ld\n", samples_param.ntotscan);

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

	key.push_back("bolo_global_file");
	value.push_back(dir.bolo_global_filename);
	comment.push_back("every scans have the same detector list");

	key.push_back("bolo_suffix");
	value.push_back(dir.bolo_suffix);
	comment.push_back("bolometers filelist suffix");

}

void export_param_saneProc(struct param_saneProc Proc_param, std::vector<string> &key, std::vector<string> &value, std::vector<string> &comment) {

	key.push_back("sampling_frequency");
	value.push_back(StringOf(Proc_param.fsamp));
	comment.push_back("detectors sampling frequency [Hz]");

	key.push_back("filter_frequency");
	value.push_back(StringOf(Proc_param.f_lp));
	comment.push_back("frequency of the high pass filter applied to the data [Hz]");

	key.push_back("apodize_Nsamples");
	value.push_back(StringOf(Proc_param.napod));
	comment.push_back("number of samples to apodize");

	key.push_back("fcut_file");
	value.push_back(Proc_param.fcut_file);
	comment.push_back("filename containing the frequency at which noise power spectra are thresholded");

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

}

void export_param_saneInv(struct param_saneInv Inv_param, std::vector<string> &key, std::vector<string> &value, std::vector<string> &comment) {

	key.push_back("noise_dir");
	value.push_back(Inv_param.noise_dir);
	comment.push_back("cov matrix directory");

	key.push_back("cov_matrix_file");
	value.push_back(Inv_param.cov_matrix_file);
	comment.push_back("this file contains the matrix you want to invert");

	key.push_back("cov_matrix_suffix");
	value.push_back(Inv_param.cov_matrix_suffix);
	comment.push_back("this file contains the matrix you want to invert");

}

void export_param_sanePS(struct param_sanePS PS_param, std::vector<string> &key, std::vector<string> &value, std::vector<string> &comment) {

	key.push_back("MixingMatrix_suffix");
	value.push_back(PS_param.mix_suffix);
	comment.push_back("Mixing matrix files suffix");

	key.push_back("MixingMatrix_global_file");
	value.push_back(PS_param.mix_global_file);
	comment.push_back("the MixingMatrix file common to all scans");

	key.push_back("ncomp");
	value.push_back(StringOf(PS_param.ncomp));
	comment.push_back("number of component(s) to estimate (default : 1)");

	key.push_back("map_file");
	value.push_back(PS_param.signame);
	comment.push_back("map substracted from the data (sanePos needed)");

	key.push_back("ell_global_file");
	value.push_back(PS_param.ell_global_file);
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
		struct samples samples_param, struct param_sanePS PS_param,
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
