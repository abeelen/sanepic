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
	string output2="";

	s = iniparser_getstring(ini,"commons:data_directory", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "commons:data_directory : default value [" + StringOf(common.dirfile) +"]\n";
	else
		common.dirfile = checkTrailingDir(StringOf(s));

	s = iniparser_getstring(ini,"commons:input_directory", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "commons:input_directory : default value [" + StringOf(common.input_dir) +"]\n";
	else
		common.input_dir = checkTrailingDir(StringOf(s));

	s = iniparser_getstring(ini,"commons:output_dir", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "commons:output_dir : default value [" + StringOf(common.output_dir) +"]\n";
	else
		common.output_dir = checkTrailingDir(StringOf(s));

	s = iniparser_getstring(ini,"commons:temp_dir", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "commons:temp_dir : default value [" + StringOf(common.tmp_dir) +"]\n";
	else
		common.tmp_dir = checkTrailingDir(StringOf(s));

	s = iniparser_getstring(ini,"commons:fits_filelist", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "commons:fits_filelist : default value [" + StringOf(common.fits_filelist) +"]\n";
	else
		common.fits_filelist = StringOf(s);

	s = iniparser_getstring(ini,"commons:bolo_suffix", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "commons:bolos_suffix : default value [" + StringOf(common.bolo_suffix) +"]\n";
	else
		common.bolo_suffix = StringOf(s);

	s = iniparser_getstring(ini,"commons:bolo_global_file", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "commons:bolos_global_file : default value [" + StringOf(common.bolo_global_filename) +"]\n";
	else
		common.bolo_global_filename = StringOf(s);

#ifdef DEBUG
	output += output2;
#endif

}

void read_param_sanePos(string &output, dictionary *ini, struct param_sanePos &pos_param){

	char *s;
	double d;
	int i;
	string output2="";

	d = iniparser_getdouble(ini,(char*)"sanePos:ra_nom", -1.0);
	if (d == -1.0)
		output2 += "sanePos:ra_nom : default value [" + StringOf(pos_param.ra_nom) +"]\n";
	else
		pos_param.ra_nom = d;

	d = iniparser_getdouble(ini,(char*)"sanePos:dec_nom", -1.0);
	if (d == -1.0)
		output2 += "sanePos:dec_nom : default value [" + StringOf(pos_param.dec_nom) +"]\n";
	else
		pos_param.dec_nom = d;

	d = iniparser_getdouble(ini,(char*)"sanePos:pixsize", -1.0);
	if (d == -1.0)
		output2 += "sanePos:pixsize : default value [" + StringOf(pos_param.pixdeg) +"]\n";
	else
		pos_param.pixdeg = d;

	s = iniparser_getstring(ini,(char*)"sanePos:proj_type", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePos:proj_type : default value [" + StringOf(pos_param.projtype) +"]\n";
	else
		pos_param.projtype = StringOf(s);

	s = iniparser_getstring(ini,"sanePos:mask_file", (char *) NULL);
	if (s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePos:mask_file : default value [" + string(pos_param.maskfile) +"]\n";
	else
		pos_param.maskfile=StringOf(s);

	i = iniparser_getint(ini,(char*)"sanePos:file_format", -1);
	if (i == -1)
		output2 += "sanePos:file_format : default value [" + StringOf(pos_param.fileFormat) +"]\n";
	else
		pos_param.fileFormat = i;


	i = iniparser_getboolean(ini, "sanePos:map_flagged_data", -1);
	if (i == -1)
		output2 += "sanePos:flgdupl : default value [" + StringOf(pos_param.flgdupl) +"]\n";
	else
		pos_param.flgdupl = (bool)i;

	i = iniparser_getboolean(ini, "sanePos:project_gaps", -1);
	if (i == -1)
		output2 += "sanePos:project_gaps : default value [" + StringOf(pos_param.projgaps) +"]\n";
	else
		pos_param.projgaps = (bool)i;


#ifdef DEBUG
	output += output2;
#endif

}

void read_param_sanePre(string &output, dictionary *ini, struct param_sanePre &proc_param){

	int i;
	double d;
	char *s;
	string output2="";


	i = iniparser_getint(ini,(char*)"sanePre:apodize_Nsamples", -1);
	if ( i == -1)
		output2 += "sanePre:apodize_Nsamples : default value [" + StringOf(proc_param.napod) +"]\n";
	else
		proc_param.napod=i;

	i = iniparser_getboolean(ini, "sanePre:nofill_gap", -1);
	if ( i == -1)
		output2 += "sanePre:nofill_gap : default value [" + StringOf(proc_param.NOFILLGAP) +"]\n";
	else
		proc_param.NOFILLGAP = (bool)i;

	d = iniparser_getdouble(ini,(char*)"sanePre:sampling_frequency", -1.0);
	if ( d == -1.0)
		output2 += "sanePre:sampling_frequency: default value [" + StringOf(proc_param.fsamp)+ "]\n";
	else
		proc_param.fsamp=d;

	d = iniparser_getdouble(ini,(char*)"sanePre:filter_frequency", -1.0);
	if (d == -1.0)
		output2 += "sanePre:filter_frequency: default value [" + StringOf(proc_param.f_lp)+ "]\n";
	else
		proc_param.f_lp=d;

	i = iniparser_getboolean(ini, "sanePre:no_baseline", -1);
	if (i == -1)
		output2 += "sanePre:no_baseline: default value [" + StringOf(proc_param.NORMLIN)+ "]\n";
	else
		proc_param.NORMLIN = (bool)i;

	i = iniparser_getboolean(ini, "sanePre:correlation", -1);
	if (i == -1)
		output2 += "sanePre:correlation: default value [" + StringOf(proc_param.CORRon)+ "]\n";
	else
		proc_param.CORRon = (bool)i;

	i = iniparser_getint(ini, "sanePre:poly_order", -1);
	if (i == -1)
		output2 += "sanePre:poly_order: default value [" + StringOf(proc_param.poly_order)+ "]\n";
	else
		proc_param.poly_order = i;

	if(proc_param.poly_order>=0)
		proc_param.remove_polynomia=1;
	else
		proc_param.remove_polynomia=0;

	s = iniparser_getstring(ini,"sanePre:fcut_file", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePre:fcut_file : default value [" + StringOf(proc_param.fcut_file) +"]\n";
	else
		proc_param.fcut_file = StringOf(s);


#ifdef DEBUG
	output += output2;
#endif
}

void read_param_saneInv(std::string &output, dictionary *ini, struct param_saneInv &saneInv_struct){

	char *s;
	string output2="";

	s = iniparser_getstring(ini,"saneInv:cov_matrix_file", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "saneInv:cov_matrix_file : default value [" + StringOf(saneInv_struct.cov_matrix_file) +"]\n";
	else
		saneInv_struct.cov_matrix_file = StringOf(s);

	s = iniparser_getstring(ini,"saneInv:cov_matrix_suffix", (char *) NULL);
	if(s == (char *) NULL || strlen(s) == 0)
		output2 += "saneInv:cov_matrix_suffix : default value [" + StringOf(saneInv_struct.cov_matrix_suffix) +"]\n";
	else
		saneInv_struct.cov_matrix_suffix = StringOf(s);

	s = iniparser_getstring(ini,"saneInv:noise_dir", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "saneInv:noise_dir : default value [" + StringOf(saneInv_struct.noise_dir) +"]\n";
	else
		saneInv_struct.noise_dir = checkTrailingDir(StringOf(s));

#ifdef DEBUG
	output += output2;
#endif

}

void read_param_sanePS(std::string &output, dictionary *ini, struct param_sanePS &sanePS_struct){

	int i;
	char *s;
	string output2="";

	i = iniparser_getint(ini, "sanePS:ncomp", -1);
	if(i == -1)
		output2 += "sanePS:ncomp : default value [" + StringOf(sanePS_struct.ncomp) +"]\n";
	else
		sanePS_struct.ncomp = i;

	s = iniparser_getstring(ini,"sanePS:map_file", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePS:map_file : default value [" + StringOf(sanePS_struct.signame) +"]\n";
	else
		sanePS_struct.signame = StringOf(s);


	s = iniparser_getstring(ini,"sanePS:MixingMatrix_Suffix", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePS:MixingMatrix_Suffix : default value [" + StringOf(sanePS_struct.mix_suffix) +"]\n";
	else
		sanePS_struct.mix_suffix = StringOf(s);

	s = iniparser_getstring(ini,"sanePS:ell_suffix", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePS:ell_suffix : default value [" + StringOf(sanePS_struct.ell_suffix) +"]\n";
	else
		sanePS_struct.ell_suffix = StringOf(s);

	s = iniparser_getstring(ini,"sanePS:ell_global_file", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePS:ell_global_file : default value [" + StringOf(sanePS_struct.ell_global_file) +"]\n";
	else
		sanePS_struct.ell_global_file = StringOf(s);

	s = iniparser_getstring(ini,"sanePS:MixingMatrix_global_file", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePS:MixingMatrix_global_file : default value [" + StringOf(sanePS_struct.mix_global_file) +"]\n";
	else
		sanePS_struct.mix_global_file = StringOf(s);

	i = iniparser_getboolean(ini, "sanePS:save_data", -1);
	if (i == -1)
		output2 += "sanePS:save_data: default value [" + StringOf(sanePS_struct.save_data)+ "]\n";
	else
		sanePS_struct.save_data = (bool)i;

	//TODO: Ugly turnaround until sanePS is released;
	s = iniparser_getstring(ini,"saneInv:cov_matrix_file", (char *) NULL);
	if( s == (char *) NULL || strlen(s) == 0)
		output2 += "saneInv:cov_matrix_file : default value [" + StringOf(sanePS_struct.cov_matrix_file) +"]\n";
	else
		sanePS_struct.cov_matrix_file = StringOf(s);

	s = iniparser_getstring(ini,"saneInv:cov_matrix_suffix", (char *) NULL);
	if(s == (char *) NULL || strlen(s) == 0)
		output2 += "saneInv:cov_matrix_suffix : default value [" + StringOf(sanePS_struct.cov_matrix_suffix) +"]\n";
	else
		sanePS_struct.cov_matrix_suffix = StringOf(s);

#ifdef DEBUG
	output += output2;
#endif

}

void read_param_sanePic(std::string &output, dictionary *ini, struct param_sanePic &sanePic_struct){

	int i;
	char *s;
	string output2="";

	i = iniparser_getint(ini, "sanePic:iterW", -1);
	if ( i == -1)
		output2 += "sanePic:iterW : default value [" + StringOf(sanePic_struct.iterw) +"]\n";
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


	i = iniparser_getint(ini, "sanePic:iterMAX", -1);
	if ( i <= 0 )
		output2 += "sanePic:iterMAX : default value [" + StringOf(sanePic_struct.itermax) +"]\n";
	else
		sanePic_struct.itermax = i;

	s = iniparser_getstring(ini,"sanePic:map_prefix", (char *) NULL);
	if(s == (char *) NULL || strlen(s) == 0)
		output2 += "sanePic:map_prefix : default value [" + StringOf(sanePic_struct.map_prefix) +"]\n";
	else
		sanePic_struct.map_prefix = StringOf(s);

#ifdef DEBUG
	output += output2;
#endif

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
			output += "Path : " + strPath + " created\n";
		else
			return 1;
	}
	return 0;
}

int compute_dirfile_format_file(std::string tmp_dir, struct samples samples_struct, int format){

	string filedir = tmp_dir + "dirfile";

	DIRFILE *D, *H, *F, *I, *I2, *J, *K, *S, *R, *R2;

	// create folders
	D = gd_open((char *)filedir.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);


	for(long iframe=0; iframe<samples_struct.ntotscan; iframe ++){

		string scan_name = samples_struct.basevect[iframe];
		string scan_folder = filedir + "/" + scan_name;
		string fdata = filedir + "/" + scan_name + "/fData";
		string index_path = filedir + "/" + scan_name + "/Indexes";
		string data = filedir + "/" + scan_name + "/data";
		string flag_dir = filedir + "/" + scan_name + "/flag";
		string RA = filedir + "/" + scan_name + "/RA";
		string DEC = filedir + "/" + scan_name + "/DEC";
		string noise_path = filedir + "/" + scan_name + "/Noise_data";
		string ell_path = filedir + "/" + scan_name + "/Noise_data/ell";

		// create folders
		S = gd_open((char *)scan_folder.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		H = gd_open((char *)index_path.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		F = gd_open((char *)fdata.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		J = gd_open((char *)data.c_str(), GD_RDWR | GD_CREAT  | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		K = gd_open((char *)flag_dir.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		I = gd_open((char *)noise_path.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		I2 = gd_open((char *)ell_path.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
		if(format==1){

			// create folders
			R = gd_open((char *)RA.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
			R2 = gd_open((char *)DEC.c_str(), GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);

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
		gd_include(D, (char *)(scan_name + "/format").c_str(), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);

		S = gd_open((char *)scan_folder.c_str(), GD_RDWR | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);

		gd_include(S, (char *)("Indexes/format"), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);
		gd_include(S, (char *)("fData/format"), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);
		gd_include(S, (char *)("flag/format"), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);
		gd_include(S, (char *)("data/format"), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);
		gd_include(S, (char *)("Noise_data/format"), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);
		gd_include(S, (char *)("Noise_data/ell/format"), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);
		if(format==1){
			gd_include(S, (char *)("RA/format"), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);
			gd_include(S, (char *)("DEC/format"), 0, GD_RDWR | GD_CREAT | GD_UNENCODED | GD_BIG_ENDIAN);
		}

		gd_flush(S,NULL);

		gd_close(S);
	}

	// close dirfile
	gd_close(D);

	return 0;
}

int cleanup_dirfile_sanePos(std::string tmp_dir, struct samples samples_struct, std::vector<std::vector<std::string> > bolo_vect)
{

	std::vector<string> det_vect;

	for(long iframe=0; iframe<samples_struct.ntotscan; iframe ++){

		det_vect=bolo_vect[iframe];

		string scan_name = samples_struct.basevect[iframe];
		string index_path = tmp_dir + "dirfile/" + scan_name + "/Indexes";

		DIRFILE *S = gd_open((char *)index_path.c_str(), GD_RDWR | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);

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
			gd_flush(S,NULL);
		}

		if(gd_close(S))
			cout << "error closing " << index_path << "-> memory leaks ..." << endl;

	}
	return 0;
}

int cleanup_dirfile_saneInv(std::string tmp_dir, struct samples samples_struct, long nframe, string noise_suffix, std::vector<std::vector<std::string> > bolo_vect)
{

	std::vector<string> det_vect;

	for (long ii=0; ii< nframe; ii++){

		det_vect=bolo_vect[ii];

		string base_name = samples_struct.basevect[ii];
		string noise_path = tmp_dir + "dirfile/" + base_name + "/Noise_data";
		string ell_path =  noise_path + "/ell";

		DIRFILE *S = gd_open((char *)noise_path.c_str(), GD_RDWR | GD_TRUNC | GD_UNENCODED | GD_BIG_ENDIAN);
		DIRFILE *D = gd_open((char *)ell_path.c_str(), GD_RDWR | GD_TRUNC | GD_UNENCODED | GD_BIG_ENDIAN);

		string suffix = base_name + noise_suffix; // base_name instead of noisevect[ii]

		for (int idet = 0; idet < (long)det_vect.size(); idet++) {

			// ell binary filename
			string outfile = det_vect[idet] + "_" + suffix + "_ell";

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
			outfile = det_vect[idet] + "_" + suffix;

			// set field information for spectra
			E.field = (char*)outfile.c_str();

			// add to the dirfile
			gd_add(S, &E);
			gd_flush(S,NULL);
		}

		if(gd_close(S))
			cout << "error closing " << noise_path << "-> memory leaks ..." << endl;
		if(gd_close(D))
			cout << "error closing " << ell_path << "-> memory leaks ..." << endl;

	}

	return 0;
}

int cleanup_dirfile_fdata(std::string tmp_dir, struct samples samples_struct, std::vector<std::vector<std::string> > bolo_vect){

	for(long iframe=0; iframe<samples_struct.ntotscan; iframe ++){

		std::vector<string> det_vect=bolo_vect[iframe];

		//get fourier transform dirfile names !
		string scan_name = samples_struct.basevect[iframe];
		string fdata_path = tmp_dir + "dirfile/" + scan_name + "/fData";

		// clean up the dirfiles with TRUNC option
		DIRFILE *S = gd_open((char *)fdata_path.c_str(), GD_RDWR | GD_TRUNC | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);

		// then generate binaries and fill format file
		string prefixe[2] = {"fdata_","fPs_"};
		for(long ip=0;ip<2;ip++)
			for(long idet=0; idet< (long)det_vect.size(); idet ++){
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
				gd_flush(S,NULL);
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
		}

	}

	return 0;
}

uint16_t check_common(string &output, struct param_common dir){

	if((dir.bolo_global_filename=="") && (dir.bolo_suffix=="")){
		output += "You must mention one of those parameters :\nparam_common:suffix or param_common:suffix\n";
		return BOLOFILE_NOT_FOUND;
	}
	if(check_path(output, dir.dirfile, "Data directory"))
		return DATA_INPUT_PATHS_PROBLEM;
	if(check_path(output, dir.input_dir, "Input directory"))
		return DATA_INPUT_PATHS_PROBLEM;
	if(check_path(output, dir.output_dir, "Output directory"))
		return OUPUT_PATH_PROBLEM;
	if(check_path(output, dir.tmp_dir, "Temporary directory"))
		return TMP_PATH_PROBLEM;

	return 0;
}

uint16_t check_param_positions(string &output, struct param_sanePos pos_param){

	if(pos_param.pixdeg < 0){
		output += "Pixsize cannot be negative ! \n";
		return  PIXDEG_WRONG_VALUE;
	}
	if((pos_param.fileFormat!=0) && (pos_param.fileFormat!=1)){
		output += "Fileformat must be 0 (SANEPIC) or 1 (HIPE) \n";
		return FILEFORMAT_NOT_FOUND;
	}

	return 0;
}

uint16_t check_param_process(string &output, struct param_sanePre proc_param){

	if(proc_param.napod<0){
		output += "You must choose a positive number of samples to apodize\n";
		return NAPOD_WRONG_VALUE;
	}
	if (proc_param.fsamp<=0.0){
		output += "sampling_frequency cannot be negative or 0 ! \n";
		return FSAMP_WRONG_VALUE;
	}
	if (proc_param.f_lp<0.0){
		output += "filter_frequency cannot be negative ! \n";
		return F_LP_WRONG_VALUE;
	}

	return 0;
}

uint16_t check_param_sanePS(string &output, struct param_sanePS structPS){

	if (structPS.ncomp<=0){
		output += "number of component ncomp cannot be negative or zero ! \n";
		return NCOMP_WRONG_VALUE;
	}
	if((structPS.ell_global_file=="") && (structPS.ell_suffix=="")){
		output += "You must mention one of those parameters :\nsanePS:ell_global_file or sanePS:ell_suffix\n";
		return ELL_FILE_NOT_FOUND;
	}
	if((structPS.mix_global_file=="") && (structPS.mix_suffix=="")){
		output += "You must mention one of those parameters :\nsanePS:mix_global_file or sanePS:mix_suffix\n";
		return MIX_FILE_NOT_FOUND;
	}

	return 0;
}

uint16_t check_param_saneInv(string &output, struct param_saneInv saneInv_struct){

	if(check_path(output, saneInv_struct.noise_dir, "Covariance Matrix directory"))
		return SANEINV_INPUT_ERROR;

	if((saneInv_struct.cov_matrix_file=="") && (saneInv_struct.cov_matrix_suffix=="")){
		output += "You must mention one of those parameters :\nsaneInv:cov_matrix_suffix or saneInv:cov_matrix_global_file\n";
		return SANEINV_INPUT_ERROR;
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
	ps_param.save_data       = 1;

	//TODO: Ugly turnaround until sanePS is released;

	ps_param.cov_matrix_file   = "";
	ps_param.cov_matrix_suffix = "_ps.fits";

}

void default_param_saneInv(struct param_saneInv &inv_param){

	inv_param.cov_matrix_file   = "";
	inv_param.cov_matrix_suffix = "_psd.fits";
	inv_param.noise_dir= "./";
}

void default_param_sanePic(struct param_sanePic &pic_param){

	pic_param.iterw      = 0;
	pic_param.itermax    = 2000;
	pic_param.save_data  = 0;
	pic_param.map_prefix = "optimMap";
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

uint16_t fill_samples_struct(string &output, struct samples &samples_struct, struct param_common &dir, struct param_saneInv &inv_param, string fcut_file){

	string filename;
	filename = dir.input_dir+dir.fits_filelist;
	if(read_fits_list(output, filename, samples_struct)!=0)
		return 0x4000;

	samples_struct.ntotscan = (samples_struct.fitsvect).size();


	// Fill basevect
	for(long iframe=0;iframe<samples_struct.ntotscan;iframe++){
		samples_struct.basevect.push_back(dirfile_Basename(samples_struct.fitsvect[iframe]));
	}
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


	if(FitsBasename(fcut_file).size()==0){
		output += "Warning ! fcut_file filename is missing in the ini file ...\n";
		return 0x8000;
	}

	// Read the fcut file
	std::vector<string> dummy2;
	if(read_strings(fcut_file,dummy2))
		return 0x8000;

	if(((int)dummy2.size())==0 || ((int) dummy2.size() != 1 && ((int) dummy2.size() != samples_struct.ntotscan) ) ){
		output += "You must provide at least one number of noise cut frequency (or one per scan) in fcut_file !\n";
		return 0x8000;
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

int get_noise_bin_sizes(std::string tmp_dir, struct samples &samples_struct, int rank)
{

	long nframe_long;


	for(long ii=0; ii< samples_struct.ntotscan; ii++){

		if(rank==0){
			string scan_name= samples_struct.basevect[ii];
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

			nframe_long=(long)(nframe-1);
		}

#ifdef USE_MPI
		MPI_Bcast(&nframe_long,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
#endif

		// get nbins value
		samples_struct.nbins.push_back(nframe_long);

		if(rank==0){

			string scan_name= samples_struct.basevect[ii];

			// get ndet value
			string filedir = tmp_dir + "dirfile/" + scan_name + "/Noise_data/";
			DIRFILE* H = gd_open((char *)filedir.c_str(), GD_RDWR | GD_VERBOSE | GD_UNENCODED | GD_BIG_ENDIAN);
			unsigned int nframe = gd_nframes(H);

			// close dirfile
			if(gd_close(H)){
				cout << "Dirfile gd_close error in get_noise_bin_sizes for : " << filedir << endl;
				return 1;
			}
			nframe_long = (long)nframe;
		}

#ifdef USE_MPI
		MPI_Bcast(&nframe_long,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
#endif

		// compute ndet considering entry size and nbins
		samples_struct.ndet.push_back(nframe_long / samples_struct.nbins[ii]);

	}

	return 0;
}

int channel_list_to_vect_list(struct samples samples_struct, std::vector<std::vector<std::string> > &bolo_vect, int rank)
{

#ifdef USE_MPI
	long size_max;
#endif

	long ndet;
	std::vector<string> det_vect;

	for (long iframe=0;iframe<samples_struct.ntotscan;iframe++){

		if(rank==0){

			string output_read = "";
			if(read_channel_list(output_read, samples_struct.bolovect[iframe], det_vect)){
				cout << output_read << endl;
#ifdef USE_MPI
				MPI_Abort(MPI_COMM_WORLD, 1);
#endif
				return 1;
			}
#ifdef USE_MPI
			compute_bololist_size(det_vect, size_max);
#endif
			ndet = (long)det_vect.size();
		}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD); // other procs wait untill rank 0 has created dirfile architecture.
		MPI_Bcast(&ndet, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(&size_max, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

		char* temp=new char[size_max+1];

		for(long ii=0; ii< ndet; ii++){
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

	for(long idet=0; idet < (long)str_vect.size(); idet++){
		string temp = str_vect[idet];
		size_buff +=temp.size();

		if((long)temp.size()>size_max)
			size_max=temp.size();
	}

	return size_buff;
}

void fill_var_sizes_struct(struct param_common dir, struct param_sanePos pos_param, struct param_sanePre proc_param,
		struct param_saneInv inv_param, struct param_sanePS ps_param, struct samples samples_struct, struct ini_var_strings &ini_v)
{

	// common
	ini_v.dirfile = (int)dir.dirfile.size();
	ini_v.sizemax = ini_v.dirfile;

	ini_v.output_dir = (int)dir.output_dir.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.output_dir ? ini_v.sizemax : ini_v.output_dir;

	ini_v.input_dir = (int)dir.input_dir.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.input_dir ? ini_v.sizemax : ini_v.input_dir;

	ini_v.tmp_dir = (int)dir.tmp_dir.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.tmp_dir ? ini_v.sizemax : ini_v.tmp_dir;

	ini_v.fits_filelist = (int)dir.fits_filelist.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.fits_filelist ? ini_v.sizemax : ini_v.fits_filelist;

	ini_v.bolo_global_filename = (int)dir.bolo_global_filename.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.bolo_global_filename ? ini_v.sizemax : ini_v.bolo_global_filename;

	ini_v.bolo_suffix = (int)dir.bolo_suffix.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.bolo_suffix ? ini_v.sizemax : ini_v.bolo_suffix;

	// sanePos
	ini_v.maskfile = (int)pos_param.maskfile.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.maskfile ? ini_v.sizemax : ini_v.maskfile;

	ini_v.projtype = (int)pos_param.projtype.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.projtype ? ini_v.sizemax : ini_v.projtype;

	// sanePre
	ini_v.fcut_file = (int)proc_param.fcut_file.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.fcut_file ? ini_v.sizemax : ini_v.fcut_file;

	// saneInv
	ini_v.noise_dir = (int)inv_param.noise_dir.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.noise_dir ? ini_v.sizemax : ini_v.noise_dir;

	ini_v.cov_matrix_file = (int)inv_param.cov_matrix_file.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.cov_matrix_file ? ini_v.sizemax : ini_v.cov_matrix_file;

	ini_v.cov_matrix_suffix = (int)inv_param.cov_matrix_suffix.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.cov_matrix_suffix ? ini_v.sizemax : ini_v.cov_matrix_suffix;

	// sanePS
	ini_v.ell_suffix = (int)ps_param.ell_suffix.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.ell_suffix ? ini_v.sizemax : ini_v.ell_suffix;

	ini_v.ell_global_file = (int)ps_param.ell_global_file.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.ell_global_file ? ini_v.sizemax : ini_v.ell_global_file;

	ini_v.signame = (int)ps_param.signame.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.signame ? ini_v.sizemax : ini_v.signame;

	ini_v.mix_global_file = (int)ps_param.mix_global_file.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.mix_global_file ? ini_v.sizemax : ini_v.mix_global_file;

	ini_v.mix_suffix = (int)ps_param.mix_suffix.size();
	ini_v.sizemax = ini_v.sizemax > ini_v.mix_suffix ? ini_v.sizemax : ini_v.mix_suffix;

	ini_v.ntotscan = (int)samples_struct.ntotscan;

	ini_v.fitsvect = new int[samples_struct.ntotscan];
	ini_v.noisevect = new int[samples_struct.ntotscan];
	ini_v.bolovect = new int[samples_struct.ntotscan];


	for(long ii = 0; ii < samples_struct.ntotscan; ii++){
		ini_v.fitsvect[ii] = (int)(samples_struct.fitsvect[ii]).size();
		ini_v.sizemax = ini_v.sizemax > ini_v.fitsvect[ii] ? ini_v.sizemax : ini_v.fitsvect[ii];
		ini_v.noisevect[ii] = (int)(samples_struct.noisevect[ii]).size();
		ini_v.sizemax = ini_v.sizemax > ini_v.noisevect[ii] ? ini_v.sizemax : ini_v.noisevect[ii];
		ini_v.bolovect[ii] = (int)(samples_struct.bolovect[ii]).size();
		ini_v.sizemax = ini_v.sizemax > ini_v.bolovect[ii] ? ini_v.sizemax : ini_v.bolovect[ii];

	}

	ini_v.sizemax = ini_v.sizemax +1;
}

void Build_derived_type_ini_var (struct ini_var_strings *ini_v,
		MPI_Datatype* message_type_ptr){

	int block_lengths[24];
	MPI_Aint displacements[24];
	MPI_Aint addresses[25];
	MPI_Datatype typelist[24];
	/* Specification des types*/
	for(int ii=0;ii<24;ii++)
		typelist[ii] = MPI_INT;

	/* Specification du nombre d’elements de chaque type */
	for(int ii=0;ii<20;ii++)
		block_lengths[ii] = 1;

	block_lengths[20] = ini_v->ntotscan;
	block_lengths[21] = ini_v->ntotscan;
	block_lengths[22] = ini_v->ntotscan;
	block_lengths[23] = ini_v->ntotscan;

	/* Calcul du deplacement de chacun des membres
	 * relativement a indata */
	MPI_Address(ini_v, &addresses[0]);

	// common
	MPI_Address(&(ini_v->dirfile), &addresses[1]);
	MPI_Address(&(ini_v->output_dir), &addresses[2]);
	MPI_Address(&(ini_v->input_dir), &addresses[3]);
	MPI_Address(&(ini_v->tmp_dir), &addresses[4]);
	MPI_Address(&(ini_v->fits_filelist), &addresses[5]);
	MPI_Address(&(ini_v->bolo_global_filename), &addresses[6]);
	MPI_Address(&(ini_v->bolo_suffix), &addresses[7]);

	// sanePos
	MPI_Address(&(ini_v->maskfile), &addresses[8]);
	MPI_Address(&(ini_v->projtype), &addresses[9]);

	// sanePre
	MPI_Address(&(ini_v->fcut_file), &addresses[10]);

	// saneInv
	MPI_Address(&(ini_v->noise_dir), &addresses[11]);
	MPI_Address(&(ini_v->cov_matrix_file), &addresses[12]);
	MPI_Address(&(ini_v->cov_matrix_suffix), &addresses[13]);

	// sanePS
	MPI_Address(&(ini_v->ell_suffix), &addresses[14]);
	MPI_Address(&(ini_v->ell_global_file), &addresses[15]);
	MPI_Address(&(ini_v->signame), &addresses[16]);
	MPI_Address(&(ini_v->mix_global_file), &addresses[17]);
	MPI_Address(&(ini_v->mix_suffix), &addresses[18]);

	// all projects
	MPI_Address(&(ini_v->sizemax), &addresses[19]);

	// samples
	MPI_Address(&(ini_v->ntotscan), &addresses[20]);
	MPI_Address(&(ini_v->fitsvect[0]), &addresses[21]);
	MPI_Address(&(ini_v->noisevect[0]), &addresses[22]);
	MPI_Address(&(ini_v->bolovect[0]), &addresses[23]);
	MPI_Address(&(ini_v->basevect[0]), &addresses[24]);

	for(int ii=0;ii<24;ii++)
		displacements[ii] = addresses[ii+1] - addresses[0];

	/* Creation du type derive */
	MPI_Type_struct(24, block_lengths, displacements, typelist,
			message_type_ptr);
	/* Remise du type pour qu’il puisse etre utilise */
	MPI_Type_commit(message_type_ptr);

}

int commit_struct_from_root(struct param_common &dir, struct param_sanePos &pos_param, struct param_sanePre &proc_param,
		struct param_saneInv &inv_param, struct param_sanePic &pic_param, struct param_sanePS &ps_param, struct samples &samples_struct, struct ini_var_strings ini_v, int rank)
{

	if(commit_param_common(dir, ini_v,  rank))
		return 1;

	commit_param_sanePos(pos_param, ini_v, rank);
	commit_param_sanePre(proc_param, ini_v, rank);
	commit_param_saneInv(inv_param, ini_v, rank);
	commit_param_sanePic(pic_param, ini_v, rank);
	commit_param_sanePS(ps_param, ini_v, rank);
	commit_samples_struct(samples_struct, ini_v, rank);

	return 0;
}


int commit_param_common(struct param_common &dir, struct ini_var_strings ini_v, int rank)
{
	int position=0 ;
	int size_buff = ini_v.dirfile + ini_v.input_dir + ini_v.output_dir + ini_v.tmp_dir + ini_v.fits_filelist +
			ini_v.bolo_global_filename + ini_v.bolo_suffix;

	char buffer[size_buff];

	if (rank == 0){

		MPI_Pack((char*)dir.dirfile.c_str(), ini_v.dirfile, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)dir.input_dir.c_str(), ini_v.input_dir, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)dir.output_dir.c_str(), ini_v.output_dir, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)dir.tmp_dir.c_str(), ini_v.tmp_dir, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)dir.fits_filelist.c_str(), ini_v.fits_filelist, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)dir.bolo_global_filename.c_str(), ini_v.bolo_global_filename, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)dir.bolo_suffix.c_str(), ini_v.bolo_suffix, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		/* Diffusion du contenu du buffer */
		MPI_Bcast(buffer, size_buff, MPI_PACKED, 0, MPI_COMM_WORLD);


	} else {

		MPI_Bcast(buffer, size_buff, MPI_PACKED, 0, MPI_COMM_WORLD);

		char *temp_char;
		temp_char = new char[ini_v.sizemax];

		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.dirfile,
				MPI_CHAR, MPI_COMM_WORLD);

		dir.dirfile = (string)temp_char;

		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.input_dir,
				MPI_CHAR, MPI_COMM_WORLD);

		dir.input_dir = (string)temp_char;

		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.output_dir,
				MPI_CHAR, MPI_COMM_WORLD);

		dir.output_dir = (string)temp_char;

		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.tmp_dir,
				MPI_CHAR, MPI_COMM_WORLD);

		dir.tmp_dir = (string)temp_char;

		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.fits_filelist,
				MPI_CHAR, MPI_COMM_WORLD);

		dir.fits_filelist = (string)temp_char;

		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.bolo_global_filename,
				MPI_CHAR, MPI_COMM_WORLD);

		dir.bolo_global_filename = (string)temp_char;

		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.bolo_suffix,
				MPI_CHAR, MPI_COMM_WORLD);

		dir.bolo_suffix = (string)temp_char;

		delete [] temp_char;

	}

	return 0;
}

int commit_param_sanePos(struct param_sanePos &pos_param, struct ini_var_strings ini_v, int rank)
{

	int position=0 ;
	int size_buff = sizeof(double)+ 2*sizeof(bool)+sizeof(int) + ini_v.maskfile + ini_v.projtype + 2*sizeof(double);

	char buffer[size_buff];

	if (rank == 0){

		MPI_Pack(&(pos_param.pixdeg), 1, MPI_DOUBLE, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(pos_param.ra_nom), 1, MPI_DOUBLE, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(pos_param.dec_nom), 1, MPI_DOUBLE, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)pos_param.maskfile.c_str(), ini_v.maskfile, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)pos_param.projtype.c_str(), ini_v.projtype, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(pos_param.flgdupl), 1, MPI_C_BOOL, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(pos_param.projgaps), 1, MPI_C_BOOL, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(pos_param.fileFormat), 1, MPI_INT, buffer, size_buff, &position, MPI_COMM_WORLD);

		/* Diffusion du contenu du buffer */
		MPI_Bcast(buffer, size_buff, MPI_PACKED, 0, MPI_COMM_WORLD);


	} else {

		MPI_Bcast(buffer, size_buff, MPI_PACKED, 0, MPI_COMM_WORLD);

		char *temp_char;
		temp_char = new char[ini_v.sizemax];

		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		/* Unpack des donnees depuis le buffer */

		MPI_Unpack(buffer, size_buff, &position, &(pos_param.pixdeg), 1,
				MPI_DOUBLE, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(pos_param.ra_nom), 1,
				MPI_DOUBLE, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(pos_param.dec_nom), 1,
				MPI_DOUBLE, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.maskfile,
				MPI_CHAR, MPI_COMM_WORLD);

		pos_param.maskfile = (string)temp_char;

		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.projtype,
				MPI_CHAR, MPI_COMM_WORLD);

		pos_param.projtype = (string)temp_char;

		MPI_Unpack(buffer, size_buff, &position, &(pos_param.flgdupl), 1,
				MPI_C_BOOL, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(pos_param.projgaps), 1,
				MPI_C_BOOL, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(pos_param.fileFormat), 1,
				MPI_INT, MPI_COMM_WORLD);


		delete [] temp_char;
	}

	return 0;
}

int commit_param_sanePre(struct param_sanePre &proc_param, struct ini_var_strings ini_v, int rank)
{

	int position=0 ;
	int size_buff = 2*sizeof(double) + 4*sizeof(bool) + sizeof(int) + sizeof(long) + ini_v.fcut_file;

	char buffer[size_buff];

	if (rank == 0){

		MPI_Pack(&(proc_param.NORMLIN), 1, MPI_C_BOOL, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(proc_param.NOFILLGAP), 1, MPI_C_BOOL, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(proc_param.CORRon), 1, MPI_C_BOOL, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(proc_param.remove_polynomia), 1, MPI_C_BOOL, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(proc_param.napod), 1, MPI_LONG, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(proc_param.poly_order), 1, MPI_INT, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(proc_param.fsamp), 1, MPI_DOUBLE, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(proc_param.f_lp), 1, MPI_DOUBLE, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)proc_param.fcut_file.c_str(), ini_v.fcut_file, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		/* Diffusion du contenu du buffer */
		MPI_Bcast(buffer, size_buff, MPI_PACKED, 0, MPI_COMM_WORLD);


	} else {

		MPI_Bcast(buffer, size_buff, MPI_PACKED, 0, MPI_COMM_WORLD);

		char *temp_char;
		temp_char = new char[ini_v.sizemax];

		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		/* Unpack des donnees depuis le buffer */

		MPI_Unpack(buffer, size_buff, &position, &(proc_param.NORMLIN), 1,
				MPI_C_BOOL, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(proc_param.NOFILLGAP), 1,
				MPI_C_BOOL, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(proc_param.CORRon), 1,
				MPI_C_BOOL, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(proc_param.remove_polynomia), 1,
				MPI_C_BOOL, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(proc_param.napod), 1,
				MPI_LONG, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(proc_param.poly_order), 1,
				MPI_INT, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(proc_param.fsamp), 1,
				MPI_DOUBLE, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(proc_param.f_lp), 1,
				MPI_DOUBLE, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.fcut_file,
				MPI_CHAR, MPI_COMM_WORLD);

		proc_param.fcut_file = (string)temp_char;

		delete [] temp_char;
	}

	return 0;
}

int commit_param_saneInv(struct param_saneInv &inv_param, struct ini_var_strings ini_v, int rank)
{

	int position=0 ;
	int size_buff = ini_v.noise_dir + ini_v.cov_matrix_file + ini_v.cov_matrix_suffix;

	char buffer[size_buff];

	if (rank == 0){

		MPI_Pack((char*)(inv_param.noise_dir).c_str(), ini_v.noise_dir, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)(inv_param.cov_matrix_file).c_str(), ini_v.cov_matrix_file, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)(inv_param.cov_matrix_suffix).c_str(), ini_v.cov_matrix_suffix, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		/* Diffusion du contenu du buffer */
		MPI_Bcast(buffer, size_buff, MPI_PACKED, 0, MPI_COMM_WORLD);


	} else {

		MPI_Bcast(buffer, size_buff, MPI_PACKED, 0, MPI_COMM_WORLD);

		char *temp_char;
		temp_char = new char[ini_v.sizemax];

		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.noise_dir,
				MPI_CHAR, MPI_COMM_WORLD);

		inv_param.noise_dir = (string)temp_char;
		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.cov_matrix_file,
				MPI_CHAR, MPI_COMM_WORLD);

		inv_param.cov_matrix_file = (string)temp_char;
		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.cov_matrix_suffix,
				MPI_CHAR, MPI_COMM_WORLD);

		inv_param.cov_matrix_suffix = (string)temp_char;
		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		delete [] temp_char;
	}

	return 0;
}

int commit_param_sanePic(struct param_sanePic &pic_param, struct ini_var_strings ini_v, int rank)
{

	int position=0 ;
	int size_buff = 4*sizeof(int);

	char buffer[size_buff];

	if (rank == 0){

		MPI_Pack(&(pic_param.iterw), 1, MPI_INT, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(pic_param.itermax), 1, MPI_INT, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(pic_param.save_data), 1, MPI_INT, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(pic_param.restore), 1, MPI_INT, buffer, size_buff, &position, MPI_COMM_WORLD);

		/* Diffusion du contenu du buffer */
		MPI_Bcast(buffer, size_buff, MPI_PACKED, 0, MPI_COMM_WORLD);


	} else {

		MPI_Bcast(buffer, size_buff, MPI_PACKED, 0, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(pic_param.iterw), 1,
				MPI_INT, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(pic_param.itermax), 1,
				MPI_INT, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(pic_param.save_data), 1,
				MPI_INT, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(pic_param.restore), 1,
				MPI_INT, MPI_COMM_WORLD);


	}

	return 0;
}

int commit_param_sanePS(struct param_sanePS &ps_param, struct ini_var_strings ini_v, int rank)
{

	int position=0 ;
	int size_buff = ini_v.ell_suffix + ini_v.mix_suffix + ini_v.ell_global_file + ini_v.mix_global_file + ini_v.cov_matrix_file
			+ ini_v.cov_matrix_suffix + ini_v.signame + sizeof(int) + 2*sizeof(bool);

	char buffer[size_buff];

	if (rank == 0){

		MPI_Pack((char*)(ps_param.ell_suffix).c_str(), ini_v.ell_suffix, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)(ps_param.mix_suffix).c_str(), ini_v.mix_suffix, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)(ps_param.ell_global_file).c_str(), ini_v.ell_global_file, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)(ps_param.mix_global_file).c_str(), ini_v.mix_global_file, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)(ps_param.cov_matrix_file).c_str(), ini_v.cov_matrix_file, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)(ps_param.cov_matrix_suffix).c_str(), ini_v.cov_matrix_suffix, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack((char*)(ps_param.signame).c_str(), ini_v.signame, MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(ps_param.ncomp), 1, MPI_INT, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(ps_param.save_data), 1, MPI_C_BOOL, buffer, size_buff, &position, MPI_COMM_WORLD);

		MPI_Pack(&(ps_param.restore), 1, MPI_C_BOOL, buffer, size_buff, &position, MPI_COMM_WORLD);

		/* Diffusion du contenu du buffer */
		MPI_Bcast(buffer, size_buff, MPI_PACKED, 0, MPI_COMM_WORLD);


	} else {

		MPI_Bcast(buffer, size_buff, MPI_PACKED, 0, MPI_COMM_WORLD);

		char *temp_char;
		temp_char = new char[ini_v.sizemax];

		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		/* Unpack des donnees depuis le buffer */

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.ell_suffix,
				MPI_CHAR, MPI_COMM_WORLD);

		ps_param.ell_suffix = (string)temp_char;
		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.mix_suffix,
				MPI_CHAR, MPI_COMM_WORLD);

		ps_param.mix_suffix = (string)temp_char;
		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.ell_global_file,
				MPI_CHAR, MPI_COMM_WORLD);

		ps_param.ell_global_file = (string)temp_char;
		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.mix_global_file,
				MPI_CHAR, MPI_COMM_WORLD);

		ps_param.mix_global_file = (string)temp_char;
		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.cov_matrix_file,
				MPI_CHAR, MPI_COMM_WORLD);

		ps_param.cov_matrix_file = (string)temp_char;
		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.cov_matrix_suffix,
				MPI_CHAR, MPI_COMM_WORLD);

		ps_param.cov_matrix_suffix = (string)temp_char;
		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.signame,
				MPI_CHAR, MPI_COMM_WORLD);

		ps_param.signame = (string)temp_char;
		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		MPI_Unpack(buffer, size_buff, &position, &(ps_param.ncomp), 1,
				MPI_INT, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(ps_param.save_data), 1,
				MPI_C_BOOL, MPI_COMM_WORLD);

		MPI_Unpack(buffer, size_buff, &position, &(ps_param.restore), 1,
				MPI_C_BOOL, MPI_COMM_WORLD);

		delete [] temp_char;
	}

	return 0;
}


int commit_samples_struct(struct samples &samples_struct, struct ini_var_strings ini_v, int rank){


	int position=0 ;
	int size_buff = sizeof(bool) + // framegiven
			ini_v.ntotscan*sizeof(long) + // nsamples
			ini_v.ntotscan*sizeof(double) + // fcut
			ini_v.ntotscan*sizeof(int);  // scans_index

	for(long ii=0; ii< ini_v.ntotscan; ii++){
		size_buff += ini_v.fitsvect[ii] + ini_v.noisevect[ii] + ini_v.bolovect[ii] +  ini_v.basevect[ii];
	}

	char buffer[size_buff];

	if (rank == 0){

		MPI_Pack(&(samples_struct.framegiven), 1, MPI_C_BOOL, buffer, size_buff, &position, MPI_COMM_WORLD);

		for(long ii=0; ii< ini_v.ntotscan; ii++){

			MPI_Pack(&(samples_struct.nsamples[ii]), 1, MPI_LONG, buffer, size_buff, &position, MPI_COMM_WORLD);

			MPI_Pack(&(samples_struct.fcut[ii]), 1, MPI_DOUBLE, buffer, size_buff, &position, MPI_COMM_WORLD);

			MPI_Pack((char*)(samples_struct.fitsvect[ii]).c_str(), ini_v.fitsvect[ii], MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

			MPI_Pack((char*)(samples_struct.noisevect[ii]).c_str(), ini_v.noisevect[ii], MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

			MPI_Pack((char*)(samples_struct.bolovect[ii]).c_str(), ini_v.bolovect[ii], MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

			MPI_Pack((char*)(samples_struct.basevect[ii]).c_str(), ini_v.basevect[ii], MPI_CHAR, buffer, size_buff, &position, MPI_COMM_WORLD);

			if(samples_struct.framegiven)
				MPI_Pack(&(samples_struct.scans_index[ii]), 1, MPI_INT, buffer, size_buff, &position, MPI_COMM_WORLD);
		}

		/* Diffusion du contenu du buffer */
		MPI_Bcast(buffer, size_buff, MPI_PACKED, 0, MPI_COMM_WORLD);


	} else {

		MPI_Bcast(buffer, size_buff, MPI_PACKED, 0, MPI_COMM_WORLD);

		samples_struct.ntotscan = (long)ini_v.ntotscan;
		long nsamp_tmp = 0;
		double fcut_tmp = 0.0;
		int scan_index_tmp = 0;
		char *temp_char;
		temp_char = new char[ini_v.sizemax];
		fill(temp_char, temp_char+ini_v.sizemax,'\0');

		/* Unpack des donnees depuis le buffer */

		MPI_Unpack(buffer, size_buff, &position, &(samples_struct.framegiven), 1,
				MPI_C_BOOL, MPI_COMM_WORLD);

		for(long ii=0; ii< samples_struct.ntotscan; ii++){

			MPI_Unpack(buffer, size_buff, &position, &nsamp_tmp, 1,
					MPI_LONG, MPI_COMM_WORLD);

			samples_struct.nsamples.push_back(nsamp_tmp);


			MPI_Unpack(buffer, size_buff, &position, &fcut_tmp, 1,
					MPI_DOUBLE, MPI_COMM_WORLD);

			samples_struct.fcut.push_back(fcut_tmp);


			MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.fitsvect[ii],
					MPI_CHAR, MPI_COMM_WORLD);

			samples_struct.fitsvect.push_back(temp_char);
			fill(temp_char, temp_char+ini_v.sizemax,'\0');

			MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.noisevect[ii],
					MPI_CHAR, MPI_COMM_WORLD);

			samples_struct.noisevect.push_back(temp_char);
			fill(temp_char, temp_char+ini_v.sizemax,'\0');

			MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.bolovect[ii],
					MPI_CHAR, MPI_COMM_WORLD);

			samples_struct.bolovect.push_back(temp_char);
			fill(temp_char, temp_char+ini_v.sizemax,'\0');

			MPI_Unpack(buffer, size_buff, &position, temp_char, ini_v.basevect[ii],
					MPI_CHAR, MPI_COMM_WORLD);

			samples_struct.basevect.push_back(temp_char);
			fill(temp_char, temp_char+ini_v.sizemax,'\0');

			if(samples_struct.framegiven){
				MPI_Unpack(buffer, size_buff, &position, &scan_index_tmp, 1,
						MPI_INT, MPI_COMM_WORLD);
				samples_struct.scans_index.push_back(scan_index_tmp);
			}
		}

		delete [] temp_char;
	}



	return 0;
}

#endif



//TODO: sanePS is not in the default distribution, so should not be in the master parser_function....
uint16_t parser_function(char * ini_name, std::string &output, struct param_common &dir,
		struct samples &samples_struct,
		struct param_sanePos &pos_param, struct param_sanePre &proc_param,
		struct param_sanePS &structPS, struct param_saneInv &saneInv_struct, struct param_sanePic &sanePic_struct, int size, int rank){

	dictionary	*	ini ;
	string filename;
	uint16_t parsed=0x0000;


	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return INI_NOT_FOUND;
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
	parsed+=fill_samples_struct(output, samples_struct, dir, saneInv_struct, filename);


	for(int iframe=0;iframe<(int)((samples_struct.fitsvect).size());iframe++){
		//		samples_struct.fitsvect[iframe] = dir.dirfile + samples_struct.fitsvect[iframe];
		samples_struct.bolovect[iframe] = dir.input_dir + samples_struct.bolovect[iframe]; // better for bolovect cause you dont need to handle path in every function call !
	}


	// Store scan sizes so that we dont need to read it again and again in the loops !
	readFrames(dir.dirfile, samples_struct.fitsvect, samples_struct.nsamples);

	if(rank==0)
		parsed+=check_common(output, dir);

	parsed+=check_param_positions(output, pos_param);

	parsed+=check_param_process(output, proc_param);

	parsed+=check_param_saneInv(output, saneInv_struct);

	parsed+=check_param_sanePS(output, structPS);

	return parsed;
}

void print_common(struct param_common dir){

	cout << "Data Dir.        : " << dir.dirfile << endl;
	cout << "Input Dir.       : " << dir.input_dir << endl;
	cout << "Temp. Dir.       : " << dir.tmp_dir << endl;
	cout << "Output Dir.      : " << dir.output_dir << endl;
}

void print_param_positions(struct param_sanePos pos_param) {

	cout << "Pixel Size       : " << setprecision(14) << pos_param.pixdeg << " deg\n";

	if(pos_param.flgdupl)
		cout << "Map Flags        : True\n";

	if(pos_param.projgaps)
		cout << "Gap Filling      : PROJECTED\n";
	else
		cout << "Gap Filling      : NOT projected (default)\n";

	cout << endl;
}

void print_param_process(struct param_sanePre proc_param){

	if(proc_param.NOFILLGAP)
		cout << "NOFILLGAPS : the gaps in data timeline WILL NOT be filled\n";
	else
		cout << "NOFILLGAPS : the gaps in data timeline WILL be filled\n";


	if(proc_param.NORMLIN)
		cout << "Simple Baseline. : NOT removed from the data\n";
	else
		cout << "Simple Baseline. : WILL BE removed from the data (default)\n";


	if(proc_param.CORRon)
		cout << "Correlations     : INCLUDED in the analysis\n";
	else
		cout << "Correlations     : NOT INCLUDED in the analysis\n";

	if(proc_param.remove_polynomia)
		cout << "Poly. Order      : " << proc_param.poly_order << endl;
	else
		cout << "Poly. Order      : No polynomia will be used\n";

	if(proc_param.napod>0)
		cout << "# for Apodize    : " << proc_param.napod << "\n";

	cout << "HP Freq.         : " << proc_param.f_lp << " Hz\n";

	cout << "Sampling Freq.   : " << proc_param.fsamp << " Hz\n";

	cout << endl;
}

void print_param_sanePic(struct param_sanePic sanepic_struct)
{

	if(sanepic_struct.save_data)
		cout << "Write Iter. Maps : " << sanepic_struct.iterw << endl;
	else
		cout << "Write Iter. Maps : OFF \n";


	cout << "Max Iter.        : " << sanepic_struct.itermax << endl;

	cout << "Maps suffix      : " << sanepic_struct.map_prefix << endl;

	cout << endl;
}

void print_param_sanePS(struct param_sanePS structPS)
{

	if(structPS.save_data)
		cout << "Save data.       : ON\n";
	else
		cout << "Save data.       : OFF\n";

	if(structPS.signame!="")
		cout << "Removed map.     : " << structPS.signame << endl;

	cout << "Noise comp.      : " << structPS.ncomp << endl;

	cout << endl;
}

void print_param_saneInv(struct param_saneInv saneInv_struct)
{

	cout << "Noise Dir.       : " << saneInv_struct.noise_dir << endl;

	cout << endl;
}

void parser_printOut(char * prog_name, struct param_common dir, struct samples samples_struct,
		struct param_sanePos pos_param, struct param_sanePre proc_param,
		struct param_sanePS structPS, struct param_sanePic sanePic_struct, struct param_saneInv saneInv_struct){

	string basename (prog_name);
	basename=FitsBasename(basename);
	int i;

	print_common(dir);
	cout << endl;

	i=basename.find("sanePos");
	if((i>=0) && (i<(int)basename.size())){
		print_param_positions(pos_param);
	}

	i=basename.find("sanePS");
	if((i>=0) && (i<(int)basename.size())){
		print_param_positions(pos_param);
		print_param_process(proc_param);
		print_param_sanePS(structPS);
	}

	i=basename.find("sanePre");
	if((i>=0) && (i<(int)basename.size())){
		print_param_positions(pos_param);
		print_param_process(proc_param);
	}

	i=basename.find("sanePic");
	if((i>=0) && (i<(int)basename.size())){
		print_param_positions(pos_param);
		print_param_process(proc_param);
		print_param_sanePic(sanePic_struct);
	}

	i=basename.find("saneInv");
	if((i>=0) && (i<(int)basename.size())){
		print_param_saneInv(saneInv_struct);
	}


	printf("# of Scans       : %ld\n",samples_struct.ntotscan);

}
