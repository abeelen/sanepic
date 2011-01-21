#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <typeinfo>
#include <sys/types.h>  // For stat()
#include <sys/stat.h>   // For stat()
#include <typeinfo>

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
#include "getdata.h"
}

#include "inputFileIO.h"
#include "parser_functions.h"
#include "mpi_architecture_builder.h"
#include "crc.h"

using namespace std;


string checkTrailingDir(string str){
	if (str[str.length()-1] != '/')
		str = str + '/';

	return str;
}

void read_common(string &output, dictionary	*ini, struct param_common &common){

	char *s;

	s = iniparser_getstring(ini,"commons:data_directory", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "commons:data_directory : default value [" + StringOf(common.dirfile) +"]\n";
	else
		common.dirfile = checkTrailingDir(StringOf(s));

	s = iniparser_getstring(ini,"commons:input_directory", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "commons:input_directory : default value [" + StringOf(common.input_dir) +"]\n";
	else
		common.input_dir = checkTrailingDir(StringOf(s));

	s = iniparser_getstring(ini,"commons:output_dir", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "commons:output_dir : default value [" + StringOf(common.output_dir) +"]\n";
	else
		common.output_dir = checkTrailingDir(StringOf(s));

	s = iniparser_getstring(ini,"commons:temp_dir", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "commons:temp_dir : default value [" + StringOf(common.tmp_dir) +"]\n";
	else
		common.tmp_dir = checkTrailingDir(StringOf(s));

	//TODO: move noise_dir to [common] or struct_saneInv ?
	s = iniparser_getstring(ini,"saneInv:noise_dir", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "saneInv:noise_dir : default value [" + StringOf(common.noise_dir) +"]\n";
	else
		common.noise_dir = checkTrailingDir(StringOf(s));

	s = iniparser_getstring(ini,"commons:fits_filelist", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "commons:fits_filelist : default value [" + StringOf(common.fits_filelist) +"]\n";
	else
		common.fits_filelist = StringOf(s);

	s = iniparser_getstring(ini,"commons:bolo_suffix", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "commons:bolos_suffix : default value [" + StringOf(common.bolo_suffix) +"]\n";
	else
		common.bolo_suffix = StringOf(s);

	s = iniparser_getstring(ini,"commons:bolo_global_file", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "commons:bolos_global_file : default value [" + StringOf(common.bolo_global_filename) +"]\n";
	else
		common.bolo_global_filename = StringOf(s);

}

void read_param_sanePos(string &output, dictionary *ini, struct param_sanePos &pos_param){

	char *s;
	double d;
	int i;

	d = iniparser_getdouble(ini,(char*)"sanePos:ra_nom", -1.0);
	if (d == -1.0)
		output += "sanePos:ra_nom : default value [" + StringOf(pos_param.ra_nom) +"]\n";
	else
		pos_param.ra_nom = d;

	d = iniparser_getdouble(ini,(char*)"sanePos:dec_nom", -1.0);
	if (d == -1.0)
		output += "sanePos:dec_nom : default value [" + StringOf(pos_param.dec_nom) +"]\n";
	else
		pos_param.dec_nom = d;

	d = iniparser_getdouble(ini,(char*)"sanePos:pixsize", -1.0);
	if (d == -1.0)
		output += "sanePos:pixsize : default value [" + StringOf(pos_param.pixdeg) +"]\n";
	else
		pos_param.pixdeg = d;

	s = iniparser_getstring(ini,(char*)"sanePos:proj_type", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output += "sanePos:proj_type : default value [" + StringOf(pos_param.projtype) +"]\n";
	else
		pos_param.projtype = StringOf(s);

	s = iniparser_getstring(ini,"sanePos:mask_file", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output += "sanePos:mask_file : default value [" + string(pos_param.maskfile) +"]\n";
	else
		pos_param.maskfile=StringOf(s);

	i = iniparser_getint(ini,(char*)"sanePos:file_format", -1);
	if (i == -1)
		output += "sanePos:file_format : default value [" + StringOf(pos_param.fileFormat) +"]\n";
	else
		pos_param.fileFormat = i;


	i = iniparser_getboolean(ini, "sanePos:map_flagged_data", -1);
	if (i == -1)
		output += "sanePos:flgdupl : default value [" + StringOf(pos_param.flgdupl) +"]\n";
	else
		pos_param.flgdupl = (bool)i;

	i = iniparser_getboolean(ini, "sanePos:project_gaps", -1);
	if (i == -1)
		output += "sanePos:project_gaps : default value [" + StringOf(pos_param.projgaps) +"]\n";
	else
		pos_param.projgaps = (bool)i;

}

void read_param_sanePre(string &output, dictionary *ini, struct param_sanePre &proc_param){

	int i;
	double d;
	char *s;

	i = iniparser_getint(ini,(char*)"sanePre:apodize_Nsamples", -1);
	if ( i == -1)
		output += "sanePre:apodize_Nsamples : default value [" + StringOf(proc_param.napod) +"]\n";
	else
		proc_param.napod=i;

	i = iniparser_getboolean(ini, "sanePre:nofill_gap", -1);
	if ( i == -1)
		output += "sanePre:nofill_gap : default value [" + StringOf(proc_param.NOFILLGAP) +"]\n";
	else
		proc_param.NOFILLGAP = (bool)i;

	d = iniparser_getdouble(ini,(char*)"sanePre:sampling_frequency", -1.0);
	if ( d == -1.0)
		output += "sanePre:sampling_frequency: default value [" + StringOf(proc_param.fsamp)+ "]\n";
	else
		proc_param.fsamp=d;

	d = iniparser_getdouble(ini,(char*)"sanePre:filter_frequency", -1.0);
	if (d == -1.0)
		output += "sanePre:filter_frequency: default value [" + StringOf(proc_param.f_lp)+ "]\n";
	else
		proc_param.f_lp=d;

	i = iniparser_getboolean(ini, "sanePre:no_baseline", -1);
	if (i == -1)
		output += "sanePre:no_baseline: default value [" + StringOf(proc_param.NORMLIN)+ "]\n";
	else
		proc_param.NORMLIN = (bool)i;

	i = iniparser_getboolean(ini, "sanePre:correlation", -1);
	if (i == -1)
		output += "sanePre:correlation: default value [" + StringOf(proc_param.CORRon)+ "]\n";
	else
		proc_param.CORRon = (bool)i;

	i = iniparser_getint(ini, "sanePre:poly_order", -1);
	if (i == -1)
		output += "sanePre:poly_order: default value [" + StringOf(proc_param.poly_order)+ "]\n";
	else
		proc_param.poly_order = i;

	if(proc_param.poly_order>=0)
		proc_param.remove_polynomia=1;
	else
		proc_param.remove_polynomia=0;

	s = iniparser_getstring(ini,"sanePre:fcut_file", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "sanePre:fcut_file : default value [" + StringOf(proc_param.fcut_file) +"]\n";
	else
		proc_param.fcut_file = StringOf(s);

}

void read_param_saneInv(std::string &output, dictionary *ini, struct param_saneInv &saneInv_struct){

	char *s;

	s = iniparser_getstring(ini,"saneInv:cov_matrix_file", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "saneInv:cov_matrix_file : default value [" + StringOf(saneInv_struct.cov_matrix_file) +"]\n";
	else
		saneInv_struct.cov_matrix_file = StringOf(s);

	s = iniparser_getstring(ini,"saneInv:cov_matrix_suffix", (char *) NULL);
	if(s == (char *) NULL || strlen(s) == 0)
		output += "saneInv:cov_matrix_suffix : default value [" + StringOf(saneInv_struct.cov_matrix_suffix) +"]\n";
	else
		saneInv_struct.cov_matrix_suffix = StringOf(s);

}

void read_param_sanePS(std::string &output, dictionary *ini, struct param_sanePS &sanePS_struct){

	int i;
	char *s;

	i = iniparser_getint(ini, "sanePS:ncomp", -1);
	if(i == -1)
		output += "sanePS:ncomp : default value [" + StringOf(sanePS_struct.ncomp) +"]\n";
	else
		sanePS_struct.ncomp = i;

	s = iniparser_getstring(ini,"sanePS:map_file", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "sanePS:map_file : default value [" + StringOf(sanePS_struct.signame) +"]\n";
	else
		sanePS_struct.signame = StringOf(s);


	s = iniparser_getstring(ini,"sanePS:MixingMatrix_Suffix", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "sanePS:MixingMatrix_Suffix : default value [" + StringOf(sanePS_struct.mix_suffix) +"]\n";
	else
		sanePS_struct.mix_suffix = StringOf(s);

	s = iniparser_getstring(ini,"sanePS:ell_suffix", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "sanePS:ell_suffix : default value [" + StringOf(sanePS_struct.ell_suffix) +"]\n";
	else
		sanePS_struct.ell_suffix = StringOf(s);

	s = iniparser_getstring(ini,"sanePS:ell_global_file", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "sanePS:ell_global_file : default value [" + StringOf(sanePS_struct.ell_global_file) +"]\n";
	else
		sanePS_struct.ell_global_file = StringOf(s);

	s = iniparser_getstring(ini,"sanePS:MixingMatrix_global_file", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "sanePS:MixingMatrix_global_file : default value [" + StringOf(sanePS_struct.mix_global_file) +"]\n";
	else
		sanePS_struct.mix_global_file = StringOf(s);


	//TODO: Ugly turnaround until sanePS is released;
	s = iniparser_getstring(ini,"saneInv:cov_matrix_file", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output += "saneInv:cov_matrix_file : default value [" + StringOf(sanePS_struct.cov_matrix_file) +"]\n";
	else
		sanePS_struct.cov_matrix_file = StringOf(s);

	s = iniparser_getstring(ini,"saneInv:cov_matrix_suffix", (char *) NULL);
	if(s == (char *) NULL || strlen(s) == 0)
		output += "saneInv:cov_matrix_suffix : default value [" + StringOf(sanePS_struct.cov_matrix_suffix) +"]\n";
	else
		sanePS_struct.cov_matrix_suffix = StringOf(s);


}

void read_param_sanePic(std::string &output, dictionary *ini, struct param_sanePic &sanePic_struct){

	int i;

	i = iniparser_getint(ini, "sanePic:iterW", -1);
	if ( i == -1)
		output += "sanePic:iterW : default value [" + StringOf(sanePic_struct.iterw) +"]\n";
	else
		sanePic_struct.iterw = i;

	if(sanePic_struct.iterw==0){
		sanePic_struct.save_data=0;
		sanePic_struct.iterw=10;
	}else{
		if(sanePic_struct.iterw<0)
			sanePic_struct.save_data=0;
		else
			sanePic_struct.save_data=1;
	}

}



int check_path(string &output, string strPath, string path_type){

	if ( access( strPath.c_str(), 0 ) == 0 )
	{
		struct stat status;
		stat( strPath.c_str(), &status );

		if ( status.st_mode & S_IFDIR )
		{
			//			cout << "The directory " << path_type << " : " << strPath << " exists." << endl;
			return 0;
		}
		else
		{
			output += "Warning : The path " + path_type + " : " + strPath + " is a file.\n";
			return 1;
		}
	}
	else
	{
		output += "Warning : Path " + path_type + " : " + strPath + " doesn't exist.\n";
		string make_it = "mkdir " + strPath;
		if(system((char*)make_it.c_str())==0)
			output += "Path : " + strPath + " created\n"; // TODO : if para_frame : only rank 0 should create folders !
		else
			return 1;
	}
	return 0;
}

int check_common(string &output, struct param_common dir){

	if((dir.bolo_global_filename=="") && (dir.bolo_suffix=="")){
		output += "You must mention one of those parameters :\nparam_common:suffix or param_common:suffix\n";
		return 1;
	}
	if(check_path(output, dir.dirfile, "Data directory"))
		return 1;
	if(check_path(output, dir.input_dir, "Input directory"))
		return 1;
	if(check_path(output, dir.output_dir, "Output directory"))
		return 1;
	if(check_path(output, dir.noise_dir, "Covariance Matrix directory"))
		return 1;
	if(check_path(output, dir.tmp_dir, "Temporary directory"))
		return 1;
	//	if(check_dirfile_paths(output, dir.tmp_dir))
	//		return 1;

	return 0;
}

int compute_dirfile_format_file(std::string tmp_dir, struct samples samples_struct){

	string filedir = tmp_dir + "dirfile";

	DIRFILE *D, *H, *F, *I, *I2, *J, *K, *S;

	// create folders
	D = gd_open((char *)filedir.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);


	for(long iframe=0; iframe<samples_struct.ntotscan; iframe ++){

		string scan_name = FitsBasename(samples_struct.fitsvect[iframe]);
		string scan_folder = filedir + "/" + scan_name;
		string fdata = filedir + "/" + scan_name + "/Fourier_transform";
		string index_path = filedir + "/" + scan_name + "/Indexes";
		string data = filedir + "/" + scan_name + "/data";
		string flag_dir = filedir + "/" + scan_name + "/flag";
		string noise_path = filedir + "/" + scan_name + "/Noise_data";
		string ell_path = filedir + "/" + scan_name + "/Noise_data/ell";

		// create folders and
		S = gd_open((char *)scan_folder.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		H = gd_open((char *)index_path.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		F = gd_open((char *)fdata.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		J = gd_open((char *)data.c_str(), GD_RDWR | GD_CREAT  | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		K = gd_open((char *)flag_dir.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);

		I = gd_open((char *)noise_path.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		I2 = gd_open((char *)ell_path.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);




		// close subdirfiles
		gd_close(H);
		gd_close(F);
		gd_close(J);
		gd_close(K);
		gd_close(S);
		gd_close(I);
		gd_close(I2);

		// include subdir and create format files
		gd_include(D, (char *)(scan_name + "/format").c_str(), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);

		S = gd_open((char *)scan_folder.c_str(), GD_RDWR | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);

		gd_include(S, (char *)("Indexes/format"), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);
		gd_include(S, (char *)("Fourier_transform/format"), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);
		gd_include(S, (char *)("flag/format"), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);
		gd_include(S, (char *)("data/format"), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);
		gd_include(S, (char *)("Noise_data/format"), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);
		gd_include(S, (char *)("Noise_data/ell/format"), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);

		gd_close(S);
	}
	//	gd_alter_protection(D, GD_PROTECT_NONE, dirf.fragment_index);
	//	gd_alter_protection(D, GD_PROTECT_NONE, dirf.fragment_index);
	//	gd_alter_protection(D, GD_PROTECT_NONE, dirf.fragment_fdata);
	//	gd_alter_protection(D, GD_PROTECT_NONE, dirf.fragment_noise);
	//	gd_alter_protection(D, GD_PROTECT_NONE, dirf.fragment_flag);
	//	gd_alter_protection(D, GD_PROTECT_NONE, dirf.fragment_data);

	// close dirfile
	gd_close(D);

	return 0;
}

int cleanup_dirfile_sanePos(string tmp_dir, struct samples samples_struct)
{

	for(long iframe=0; iframe<samples_struct.ntotscan; iframe ++){
		string scan_name = FitsBasename(samples_struct.fitsvect[iframe]);
		string index_path = tmp_dir + "dirfile/" + scan_name + "/Indexes";

		DIRFILE *S = gd_open((char *)index_path.c_str(), GD_RDWR | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);

		string output_read = "";
		std::vector<string> det_vect;
		if(read_channel_list(output_read, samples_struct.bolovect[iframe], det_vect)){
			cout << output_read << endl;
			return 1;
		}

		// then generate binaries and fill format file
		for(long idet=0; idet<(long)det_vect.size(); idet ++){
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

		}

		gd_close(S);

	}
	return 0;
}

int cleanup_dirfile_saneInv(string tmp_dir, struct samples samples_struct, long n_iter, string noise_suffix)
{

	for (long ii=0; ii< n_iter; ii++){

		string base_name = FitsBasename(samples_struct.fitsvect[ii]);
		string noise_path = tmp_dir + "dirfile/" + base_name + "/Noise_data";
		string ell_path =  noise_path + "/ell";

		DIRFILE *S = gd_open((char *)noise_path.c_str(), GD_RDWR | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		DIRFILE *D = gd_open((char *)ell_path.c_str(), GD_RDWR | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);

		string suffix = base_name + noise_suffix; // base_name instead of noisevect[ii] FitsBasename
		std::vector<string> bolos;

		string output_read = "";
		if(read_channel_list(output_read, samples_struct.bolovect[ii], bolos)){
			cout << output_read << endl;
			return 1;
		}
		long ndet = (long)bolos.size();

		for (int idet = 0; idet < ndet; idet++) {

			// ell binary filename
			string outfile = bolos[idet] + "_" + suffix + "_ell";

			// configure dirfile field for ell
			gd_entry_t E;
			E.field = (char*)outfile.c_str();
			E.field_type = GD_RAW_ENTRY;
			E.fragment_index = 0;
			E.spf = 1;
			E.data_type = GD_DOUBLE;
			E.scalar[0] = NULL;

			// add to the dirfile
			gd_add(D, &E);


			// spectra filename
			outfile = bolos[idet] + "_" + suffix;

			// set field information for spectra
			E.field = (char*)outfile.c_str();

			// add to the dirfile
			gd_add(S, &E);

		}

		gd_close(S);
		gd_close(D);

	}

	return 0;
}

int cleanup_dirfile_fdata(std::string tmp_dir, struct samples samples_struct){

	for(long iframe=0; iframe<samples_struct.ntotscan; iframe ++){

		//get fourier transform dirfile names !
		string scan_name = FitsBasename(samples_struct.fitsvect[iframe]);
		string fdata_path = tmp_dir + "dirfile/" + scan_name + "/Fourier_transform";

		// clean up the dirfiles with TRUNC option
		DIRFILE *S = gd_open((char *)fdata_path.c_str(), GD_RDWR | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);

		string output_read = "";
		std::vector<string> det_vect;
		if(read_channel_list(output_read, samples_struct.bolovect[iframe], det_vect)){
			cout << output_read << endl;
			return 1;
		}

		// then generate binaries and fill format file
		string prefixe[2] = {"fdata_","fPs_"};
		for(long ip=0;ip<2;ip++)
			for(long idet=0; idet<(long)det_vect.size(); idet ++){
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

			}

		gd_close(S);

		// check sizes in Indexes, data and flag format
		string indexes_path = tmp_dir + "dirfile/" + scan_name + "/Indexes";
		string data_path = tmp_dir + "dirfile/" + scan_name + "/data";
		DIRFILE *I = gd_open((char *)indexes_path.c_str(), GD_RDWR | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		DIRFILE *D = gd_open((char *)data_path.c_str(), GD_RDWR | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);

		long nframeI = gd_nframes(I);
		long nframeD = gd_nframes(I);

		long ns = samples_struct.nsamples[iframe];
		gd_close(I);
		gd_close(D);

		if((nframeI != ns) ||  (nframeD != ns)){
			cout << "Error... Dirfile data or Indexes has incorrect size !!\n";
			return 1;
		}else{
			cout << "Checked and size ok : " << iframe << endl; // TODO : remove after final changes
		}

	}


	return 0;
}

int check_param_positions(string &output, struct param_sanePos pos_param){

	if(pos_param.pixdeg < 0){
		output += "Pixsize cannot be negative ! \n";
		return 1 ;
	}
	if((pos_param.fileFormat!=0) && (pos_param.fileFormat!=1)){
		output += "Fileformat must be 0 (SANEPIC) or 1 (HIPE) \n";
		return 1;
	}

	return 0;
}

int check_param_process(string &output, struct param_sanePre proc_param){

	if(proc_param.napod<0){
		output += "You must choose a positive number of samples to apodize\n";
		return 1 ;
	}
	if (proc_param.fsamp<=0.0){
		output += "sampling_frequency cannot be negative or 0 ! \n";
		return 1;
	}
	if (proc_param.f_lp<0.0){
		output += "filter_frequency cannot be negative ! \n";
		return 1;
	}

	return 0;
}

int check_param_sanePS(string &output, struct param_sanePS structPS){

	if (structPS.ncomp<=0){
		output += "number of component ncomp cannot be negative or zero ! \n";
		return 1;
	}
	if((structPS.ell_global_file=="") && (structPS.ell_suffix=="")){
		output += "You must mention one of those parameters :\nsanePS:ell_global_file or sanePS:ell_suffix\n";
		return 1;
	}
	if((structPS.mix_global_file=="") && (structPS.mix_suffix=="")){
		output += "You must mention one of those parameters :\nsanePS:mix_global_file or sanePS:mix_suffix\n";
		return 1;
	}

	return 0;
}

int check_param_saneInv(string &output, struct param_saneInv saneInv_struct){

	if((saneInv_struct.cov_matrix_file=="") && (saneInv_struct.cov_matrix_suffix=="")){
		output += "You must mention one of those parameters :\nsaneInv:cov_matrix_suffix or saneInv:cov_matrix_global_file\n";
		return 1;
	}

	return 0;
}


void default_param(struct param_common &dir, struct samples &samples_struct, struct param_sanePos &pos_param, struct param_sanePre &proc_param,
		struct param_saneInv &inv_param, struct param_sanePic &pic_param){

	default_param_common(dir);
	default_param_sanePos(pos_param);
	default_param_sanePre(proc_param);
	default_param_saneInv(inv_param);
	default_param_sanePic(pic_param);

}

void default_param_common(struct param_common &dir){

	// param_common

	dir.dirfile    = "./";
	dir.output_dir = "./";
	dir.tmp_dir    = "./";
	dir.noise_dir  = "./";
	dir.input_dir  = "./";

	dir.fits_filelist        = "";
	dir.bolo_global_filename = "";
	dir.bolo_suffix          = ".bolo";


}

void default_param_sanePos(struct param_sanePos &pos_param){
	//	pos_param.flgdupl = 0; // map duplication factor
	//	pos_param.flgdupl = 0; // map duplication factor
	//	pos_param.maskfile = "";
	// param_sanePos
	pos_param.maskfile   = "" ;

	pos_param.pixdeg     = 0.0 ;
	pos_param.ra_nom     = 0.0 ;
	pos_param.dec_nom    = 0.0 ;
	pos_param.projtype   = "EQ";

	pos_param.flgdupl    = false ; // What to do with flagged data : (default : False -- map in a single pixel)
	pos_param.projgaps   = false ; // What to do with gaps : (default : 0 -- no projection)
	pos_param.fileFormat = 0 ;     // Default sanepic File Format

	// 0: sanepic format with reference position & offsets
	// 1: 'hipe' like format with RA/DEC for each time/bolo

}

void default_param_sanePre(struct param_sanePre &proc_param){

	//	proc_param.napod  = 0; /*! number of samples to apodize, =0 -> no apodisation */
	//	proc_param.NOFILLGAP = 0; /*! dont fill the gaps ? default is NO => the program fill */
	//	proc_param.fsamp = 0.0;// sampling frequency
	//	proc_param.NORMLIN = 0; /*!  baseline is removed from the data, NORMLIN = 1 else 0 */
	//	proc_param.CORRon = 1; /*! correlation included in the analysis (=1), else 0, default 0*/
	//	proc_param.remove_polynomia = 1; /*! remove a polynomia fitted to the data*/
	//	proc_param.f_lp = 0.0; // low pass filter frequency

	proc_param.NORMLIN          = false;
	proc_param.NOFILLGAP        = false;
	proc_param.CORRon           = true;
	proc_param.remove_polynomia = true;
	proc_param.napod            = 100;
	proc_param.poly_order       = 1;
	proc_param.fsamp            = 0.0;
	proc_param.f_lp             = 0.0;
	proc_param.fcut_file        = "";

}

void default_param_sanePS(struct param_sanePS &ps_param){

	ps_param.ell_suffix      = ".ell";
	ps_param.mix_suffix      = ".mix";
	ps_param.ell_global_file = "";
	ps_param.mix_global_file = "";
	ps_param.signame         = "";
	ps_param.ncomp           = 1;

	//TODO: Ugly turnaround until sanePS is released;

	ps_param.cov_matrix_file   = "";
	ps_param.cov_matrix_suffix = "_ps.fits";

}

void default_param_saneInv(struct param_saneInv &inv_param){

	inv_param.cov_matrix_file   = "";
	inv_param.cov_matrix_suffix = "_psd.fits";
}

void default_param_sanePic(struct param_sanePic &pic_param){

	//	sanePic_struct.iterw=10;
	//	sanePic_struct.save_data=0;
	pic_param.iterw     = 0;
	pic_param.save_data = 0;
}


void fill_sanePS_struct(struct param_sanePS &structPS, struct samples &samples_struct, struct param_common &dir){


	for(long iframe=0;iframe<samples_struct.ntotscan;iframe++){
		if(structPS.mix_global_file != "")
			samples_struct.mix_names.push_back(dir.input_dir+structPS.mix_global_file);
		else
			samples_struct.mix_names.push_back(dir.input_dir+FitsBasename(samples_struct.fitsvect[iframe]) + structPS.mix_suffix);

		if(structPS.ell_global_file != "")
			samples_struct.ell_names.push_back(dir.input_dir+structPS.ell_global_file);
		else
			samples_struct.ell_names.push_back(dir.input_dir+FitsBasename(samples_struct.fitsvect[iframe]) + structPS.ell_suffix);
	}

}

int fill_samples_struct(string &output, struct samples &samples_struct, struct param_common &dir, struct param_saneInv &inv_param, string fcut_file){

	string filename;
	filename = dir.input_dir+dir.fits_filelist;
	if(read_fits_list(output, filename, samples_struct)!=0)
		return 1;

	samples_struct.ntotscan = (samples_struct.fitsvect).size();


	// Fill bolovect
	for(long iframe=0;iframe<samples_struct.ntotscan;iframe++){
		if(dir.bolo_global_filename!="")
			samples_struct.bolovect.push_back(dir.bolo_global_filename);
		else
			samples_struct.bolovect.push_back(FitsBasename(samples_struct.fitsvect[iframe]) + dir.bolo_suffix);

	}

	// Fill noisevect
	for(long iframe = 0; iframe < samples_struct.ntotscan ; iframe ++){
		if((inv_param.cov_matrix_file!=""))
			samples_struct.noisevect.push_back(inv_param.cov_matrix_file);
		else
			samples_struct.noisevect.push_back(FitsBasename(samples_struct.fitsvect[iframe]) + inv_param.cov_matrix_suffix);
	}


	// Read the fcut file
	std::vector<string> dummy2;
	if(read_strings(fcut_file,dummy2))
		return 1;

	if(((int)dummy2.size())==0 || ((int) dummy2.size() != 1 && ((int) dummy2.size() != samples_struct.ntotscan) ) ){
		output += "You must provide at least one number of noise cut frequency (or one per scan) in fcut_file !\n";
		return 1;
	}

	for(int iframe_dummy=0; iframe_dummy<(int)dummy2.size(); iframe_dummy++)
		samples_struct.fcut.push_back(atof(dummy2[iframe_dummy].c_str()));

	// If There is only one value, use it everywhere
	if((dummy2.size() == 1) && (samples_struct.ntotscan > 1)){
		// if only one fcut, extend to all scans
		samples_struct.fcut.resize(samples_struct.ntotscan, samples_struct.fcut[0]);
	}

	return 0;


}


int get_noise_bin_sizes(std::string tmp_dir, struct samples &samples_struct)
{

	for(long ii=0; ii< samples_struct.ntotscan; ii++){

		string scan_name= FitsBasename(samples_struct.fitsvect[ii]);
		// dirfile path
		string filedir = tmp_dir + "dirfile/" + scan_name + "/Noise_data/ell/";

		// open dirfile
		DIRFILE* H = gd_open((char *)filedir.c_str(), GD_RDWR | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		unsigned int nframe = gd_nframes(H);

		// close dirfile
		if(gd_close(H)){
			cout << "Dirfile gd_close error in get_noise_bin_sizes for : " << filedir << endl;
			return 1;
		}

		// get nbins value
		samples_struct.nbins.push_back((long)(nframe-1));

		// get ndet value
		filedir = tmp_dir + "dirfile/" + scan_name + "/Noise_data/";
		H = gd_open((char *)filedir.c_str(), GD_RDWR | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		nframe = gd_nframes(H);

		// close dirfile
		if(gd_close(H)){
			cout << "Dirfile gd_close error in get_noise_bin_sizes for : " << filedir << endl;
			return 1;
		}

		// compute ndet considering entry size and nbins
		samples_struct.ndet.push_back(nframe / samples_struct.nbins[ii]);

	}

	return 0;
}


//TODO: sanePS is not in the default distribution, so should not be in the master parser_function....
int parser_function(char * ini_name, std::string &output, struct param_common &dir,
		struct samples &samples_struct,
		struct param_sanePos &pos_param, struct param_sanePre &proc_param,
		struct param_sanePS &structPS, struct param_saneInv &saneInv_struct, struct param_sanePic &sanePic_struct, int size, int rank){

	dictionary	*	ini ;
	string filename;


	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return 2;
	}

	// default values :
	default_param( dir, samples_struct, pos_param, proc_param, saneInv_struct, sanePic_struct);

	//TODO sanePS should be out of here
	default_param_sanePS(structPS);

	read_common(output, ini, dir);
	read_param_saneInv(output, ini, saneInv_struct);
	read_param_sanePos(output, ini, pos_param);
	read_param_sanePic(output, ini, sanePic_struct);
	read_param_sanePre(output, ini, proc_param);
	//TODO sanePS should be out of here (special case)
	read_param_sanePS(output, ini, structPS);

	iniparser_freedict(ini);

	// Now the ini file has been read, do the rest

	// Fill fitsvec, noisevect, scans_index with values read from the 'str' filename
	filename=dir.input_dir+proc_param.fcut_file;
	if (fill_samples_struct(output, samples_struct, dir, saneInv_struct, filename) !=0 )
		return 1;


	// TODO: Fundamental reason for that here ? Why not keep it simple ?
	for(int iframe=0;iframe<(int)((samples_struct.fitsvect).size());iframe++){
		samples_struct.fitsvect[iframe] = dir.dirfile + samples_struct.fitsvect[iframe];
		samples_struct.bolovect[iframe] = dir.input_dir + samples_struct.bolovect[iframe];
	}


	// TODO: Why is that needed in sample_struct, why not just for saneFrameOrder ?
	readFrames(samples_struct.fitsvect, samples_struct.nsamples);

	if(check_common(output, dir) && (rank==0)) // TODO : indiquer aux autres rank qu'il y a eu erreur et qu'ils doivent sortir !
		return 1;
	if(check_param_positions(output, pos_param))
		return 1;
	if(check_param_process(output, proc_param))
		return 1;
	if(check_param_saneInv(output, saneInv_struct))
		return 1;
	if(check_param_sanePS(output, structPS))
		return 1;

	return 0;
}

void print_common(struct param_common dir){

	cout << "Data directory      : " << dir.dirfile << endl;
	cout << "Input directory     : " << dir.input_dir << endl;
	cout << "Temporary directory : " << dir.tmp_dir << endl;
	cout << "Output directory    : " << dir.output_dir << endl;
	cout << "Noise directory     : " << dir.noise_dir << endl;
	cout << endl;
}

void print_param_positions(struct param_sanePos pos_param) {

	cout << "Pixel size : " << setprecision(14) << pos_param.pixdeg << " deg\n";

	if(pos_param.flgdupl)
		cout << "Separate map : " << setw(10) << "map_flagged_data = True\n";

	if(pos_param.projgaps)
		cout << "GAP FILLING : " << setw(10) << "PROJECTED\n";
	else
		cout << "GAP FILLING : " << setw(10) << "NOT projected (default)\n";

	cout << endl;
}

void print_param_process(struct param_sanePre proc_param){

	if(proc_param.NOFILLGAP)
		cout << "NOFILLGAPS : " << setw(27) << "the gaps in data timeline WILL NOT be filled\n";
	else
		cout << "NOFILLGAPS : " << setw(27) << "the gaps in data timeline WILL be filled\n";


	if(proc_param.NORMLIN)
		cout << "Baseline : " << setw(10) << "NOT removed from the data\n";
	else
		cout << "Baseline : " << setw(10) << "WILL BE removed from the data (default)\n";


	if(proc_param.CORRon)
		cout << "Correlations : " << setw(10) << "INCLUDED in the analysis\n";
	else
		cout << "Correlations : " << setw(10) << "NOT INCLUDED in the analysis\n";

	if(proc_param.remove_polynomia)
		cout << "Polynomia order : " << setw(18) << proc_param.poly_order << endl;
	else
		cout << "Polynomia : " << setw(10) << "No polynomia will be used\n";

	if(proc_param.napod>0)
		cout << "Number of samples to apodize : " << setw(7) << proc_param.napod << "\n";

	cout << "HP Filter Frequency : " << setw(18) << proc_param.f_lp << " Hz\n";

	cout << "Sampling frequency : " << setw(16) <<proc_param.fsamp << " Hz\n";

	cout << endl;
}

void print_param_sanePic(struct param_sanePic sanepic_struct)
{

	if(sanepic_struct.save_data)
		cout << "Saving temporary map every iterW : " << sanepic_struct.iterw << " iterations\n";
	else
		cout << "Save_data is OFF : You will not be able to restore the iterations done in case of crashes \n";

	cout << endl;
}

void print_param_sanePS(struct param_sanePS structPS)
{

	if(structPS.signame!="")
		cout << "A map will be removed from the data signal before estimation of the noise : " << structPS.signame << endl;

	cout << "Number of noise component to estimate : " << structPS.ncomp << endl;

	cout << endl;
}

void parser_printOut(char * prog_name, struct param_common dir, struct samples samples_struct,
		struct param_sanePos pos_param, struct param_sanePre proc_param,
		struct param_sanePS structPS, struct param_sanePic sanePic_struct){

	string basename (prog_name);
	basename=FitsBasename(basename);
	int i;

	cout << "\nYou have specified the following options : \n\n";

	print_common(dir);
	cout << endl;

	i=basename.find("sanePos");
	if((i>=0) && (i<(int)basename.size())){
		//		cout << "sanePos detected !!!\n";
		print_param_positions(pos_param);
	}

	i=basename.find("sanePS");
	if((i>=0) && (i<(int)basename.size())){
		//		cout << "sanePS detected !!!\n";
		print_param_positions(pos_param);
		print_param_process(proc_param);
		print_param_sanePS(structPS);
	}

	i=basename.find("sanePre");
	if((i>=0) && (i<(int)basename.size())){
		//		cout << "sanePre detected !!!\n";
		print_param_positions(pos_param);
		print_param_process(proc_param);
	}

	i=basename.find("sanePic");
	if((i>=0) && (i<(int)basename.size())){
		//		cout << "sanePic detected !!!\n";
		print_param_positions(pos_param);
		print_param_process(proc_param);
		print_param_sanePic(sanePic_struct);
	}

	printf("Number of scans      : %ld\n",samples_struct.ntotscan);

}
