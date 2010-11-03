#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <unistd.h>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/types.h>  // For stat()
#include <sys/stat.h>   // For stat()


extern "C"{
#include "iniparser.h"
#include "dictionary.h"
#include "nrutil.h"
}

#include "inputFileIO.h"
#include "parser_functions.h"
#include "struct_definition.h"
#include "mpi_architecture_builder.h"
#include "covMatrix_IO.h"
#include "crc.h"


using namespace std;


template <class T>
bool from_string(T &Value,const std::string &str,std::ios_base & (*f)(std::ios_base&))
// converts string to streamable value, and returns true on success and false otherwise.
{
	std::istringstream stream(str);
	stream>>f>>Value;
	return (!stream.fail()) && stream.get()==std::istringstream::traits_type::eof();
}

template <class T>
T convert_from_string(string str)
{
	T value;
	string type = (typeid (T).name());
	char type_c = (char)type[0];
	//	cout << "type : " << type_c << endl;
	switch(type_c){
	case 'i' :
		value = atoi((char*)str.c_str());
		break;
	case 'l' :
		value = atol((char*)str.c_str());
		break;
	case 'd' :
		value = atof((char*)str.c_str());
		break;
	case 'b' :
		if(str.length()>1){
			std::transform(str.begin(), str.end(),
					str.begin(), ::tolower);
			if(!from_string(value,str,std::boolalpha))
				cout << "from_string failed !!!" << endl;
		}else
			if(!from_string(value,str,std::noboolalpha))
				cout << "from_string failed !!!" << endl;
		break;
	}

	return value;
}

template <class T>
int read_item(string &output, dictionary *ini, T &value, T default_value, string field, bool mandatory)
{

	string str;

	switch(read_parser_string(ini, field, str)){
	case 2:
		if(mandatory)
			return 2;
		else{
			output += "WARNING ! You must add a line in ini file specifying : " + field + "\n";
			output +=  "Using default : " + StringOf(default_value) + "\n";
			value=default_value;
		}
		break;
	case 1:
		if(mandatory)
			return 1;
		else{
			output += "Key is empty : You must specify : " + field + "\n";
			output += "Using default : " + StringOf(default_value) + "\n";
			value=default_value;
		}
		break;
	case 0:
		value=convert_from_string<T>(str);
		break;
	}

	return 0;
}

int read_dir(string &output, dictionary	*ini, struct param_common &dir, string dirtype){

	string str;

	switch(read_parser_string(ini, dirtype, str)){
	case 2:
		output += "WARNING ! You must add a line in ini file specifying : " + dirtype + "\n";
		output +=  "Using default directory : ./ \n";
		str="./";
		break;
	case 1:
		output += "Key is empty : You must specify : commons:data_directory\n";
		output += "Using default directory : ./ \n";
		str="./";
		break;
	case 0:
		if (str[str.length()-1] != '/')
			str = str + '/';
		break;
	}

	switch((unsigned int)checksum((char*)dirtype.c_str(), (size_t)dirtype.size(), 0)){
	case (unsigned int)2308: // commons:data_directory
	dir.dirfile=str;
	break;

	case (unsigned int)1674: // commons:temp_dir
	dir.tmp_dir=str;
	break;

	case (unsigned int)1925: // commons:output_dir
	dir.output_dir=str;
	break;

	case (unsigned int)1738: 	// saneInv:noise_dir
	dir.noise_dir=str;
	break;

	case (unsigned int)2458: 	// commons:input_dir
	dir.input_dir=str;
	break;

	}


	return 0;


}

int read_fits_file_list(string &output, dictionary	*ini, struct param_common &dir, struct samples &samples_str){

	string str;


	switch(read_parser_string(ini, "commons:fits_filelist", str)){
	case 2:
		output += "You must add a line in ini file specifying : commons:fits_filelist\n";
		return 1;
	case 1:
		output += "Key is empty : You must specify : commons:fits_filelist\n";
		return 1;
	case 0:
		samples_str.filename=dir.input_dir + str;

		// Fill fitsvec, noisevect, scans_index with values read from the 'str' filename
		if(read_fits_list(output, samples_str.filename, \
				samples_str.fitsvect, samples_str.scans_index, \
				samples_str.framegiven)!=0)
			return 1;


		for(int ii=0;ii<(int)((samples_str.fitsvect).size());ii++){
#ifdef DEBUG_PRINT
			cout << dir.dirfile + samples_str.fitsvect[ii] << endl;
#endif
			samples_str.fitsvect[ii] = dir.dirfile + samples_str.fitsvect[ii];
		}

		// Populate the nsamples vector
		readFrames(samples_str.fitsvect, samples_str.nsamples);
	}
	return 0;
}

int read_apodize_samples(string &output, dictionary	*ini, struct param_sanePre &proc_param){

	if(read_item<long>(output, ini, proc_param.napod, 100, "sanePre:apodize_Nsamples", 1))
		return 2;

	if(proc_param.napod<0){
		output += "You must choose a positive number of samples to apodize\n";
		return 1 ;
	}

	return 0;

}

int read_nofillgap(std::string &output, dictionary	*ini, struct param_sanePre &proc_param){

	if(read_item<bool>(output, ini, proc_param.NOFILLGAP, 0, "sanePre:nofill_gap", 0))
		return 2;

	return 0;

}

int read_sampling_frequency(string &output, dictionary	*ini, struct param_sanePre &proc_param){

	if(read_item<double>(output, ini, proc_param.fsamp, -1.0, "sanePre:sampling_frequency", 1))
		return 2;

	if (proc_param.fsamp<=0.0){
		output += "sampling_frequency cannot be negative or 0 ! \n";
		return 1;
	}

	return 0;
}


int read_filter_frequency(string &output, dictionary *ini, struct param_sanePre &proc_param){

	if(read_item<double>(output, ini, proc_param.f_lp, -1.0, "sanePre:filter_frequency", 1))
		return 2;

	if (proc_param.f_lp<0.0){
		output += "filter_frequency cannot be negative ! \n";
		return 1;
	}

	return 0;
}


int read_noise_cut_freq(string &output, dictionary	*ini, struct param_common dir, struct param_sanePre &proc_param, std::vector<double> &fcut){

	string str;


	switch(read_parser_string(ini, "sanePre:fcut_file", str)){
	case 2:
		output += "Warning ! You must add a line in ini file specifying : sanePre:fcut_file\n";
		output += "Using min covariance frequency value as default\n";
		break;
	case 1:
		output += "Key is empty : You must specify : sanePre:fcut_file\n";
		output += "Using min covariance frequency value as default\n";
		break;
	case 0:
		proc_param.fcut_file = dir.input_dir + str;

		std::vector<string> dummy2;
		if(read_strings(proc_param.fcut_file,dummy2))
			return 1;

		if(((int)dummy2.size())==0){
			output += "You must provide at least one number of noise cut frequency (or one per scan) in fcut_file !\n";
			return 1;}
		for(int ii=0; ii<(int)dummy2.size(); ii++)
			fcut.push_back(atof(dummy2[ii].c_str()));
	}

	return 0;

}


int read_baseline(std::string &output, dictionary	*ini, struct param_sanePre &proc_param){

	if(read_item<bool>(output, ini, proc_param.NORMLIN, 0, "sanePre:no_baseline", 0))
		return 2;

	return 0;

}

int read_correlation(std::string &output, dictionary	*ini, struct param_sanePre &proc_param){

	if(read_item<bool>(output, ini, proc_param.CORRon, 1, "sanePre:correlation", 0))
		return 2;

	return 0;
}

int read_remove_poly(std::string &output, dictionary	*ini, struct param_sanePre &proc_param){

	if(read_item<int>(output, ini, proc_param.poly_order, 1, "sanePre:poly_order", 0))
		return 2;

	if(proc_param.poly_order>=0)
		proc_param.remove_polynomia=1;
	else
		proc_param.remove_polynomia=0;

	return 0;
}


int read_iter(std::string &output, dictionary	*ini, int &iterw){

	if(read_item<int>(output, ini, iterw, 10, "sanePic:iterW", 0))
		return 2;

	return 0;
}

int read_ell_suffix(string &output, dictionary	*ini, string &ell_suffix){

	string str;

	switch(read_parser_string(ini, "sanePS:ell_suffix", str)){
	case 2:
		output += "You must add a line in ini file specifying : sanePS:ell_suffix\n";
		ell_suffix="";
		break;
	case 1:
		ell_suffix="";
		break;
	case 0:
		ell_suffix=str;
		break;
	}

	return 0;
}

int read_ell_global_file(string &output, dictionary	*ini, string &ell_global_file){

	string str;

	switch(read_parser_string(ini, "sanePS:ell_global_file", str)){
	case 2:
		output += "You must add a line in ini file specifying : sanePS:ell_global_file\n";
		ell_global_file="";
		break;
	case 1:
		ell_global_file="";
		break;
	case 0:
		ell_global_file=str;
		break;
	}

	return 0;
}

int read_map_file(dictionary	*ini, string &signame){

	string str;

	if (read_parser_string(ini,"sanePS:map_file", str)){
		signame="NOSIGFILE";
	}else{
		signame=str;
	}
	return 0;
}

int read_cov_matrix_file(string &output, dictionary	*ini, string &fname){

	string str;

	switch(read_parser_string(ini, "saneInv:cov_matrix_file", str)){
	case 2:
		output += "You must add a line in ini file specifying : saneInv:cov_matrix_file\n";
		break;
	case 0:
		fname=str;
		break;
	}
	return 0;

}

int read_cov_matrix_suffix(string &output, dictionary	*ini, string &fname){

	string str;

	switch(read_parser_string(ini, "saneInv:cov_matrix_suffix", str)){
	case 2:
		output += "You must add a line in ini file specifying : saneInv:cov_matrix_suffix\n";
		break;
	case 0:
		fname=str;
		break;
	}
	return 0;

}

int read_mixmatfile_suffix(dictionary	*ini, string &MixMat_suffix){

	string str;

	if (read_parser_string(ini,"sanePS:MixingMatrix_Suffix", str)){
		MixMat_suffix = "";
	}else{
		MixMat_suffix = str;
	}

	return 0;
}

int read_mixmat_global_file(dictionary	*ini, string &MixMat_global){

	string str;

	if (read_parser_string(ini,"sanePS:MixingMatrix_global_file", str)){
		MixMat_global = "";
	}else{
		MixMat_global = str;
	}

	return 0;
}

int read_ncomp(string &output, dictionary	*ini, long &ncomp){

	if(read_item<long>(output, ini, ncomp, 1, "sanePS:ncomp", 1))
		return 2;

	if (ncomp<=0){
		output += "number of component ncomp cannot be negative or zero ! \n";
		return 1;
	}

	return 0;
}


int read_fcut(string &output, dictionary	*ini, double &fcut){

	if(read_item<double>(output, ini, fcut, 12.0, "sanePS:fcut", 0))
		return 2;

	if (fcut<0.0){
		output += "noise cut frequency cannot be negative ! \n";
		return 1;
	}

	return 0;
}

int read_parser_string(dictionary	*ini, string line, string & str){
	char *s;

	s = iniparser_getstring(ini, line.c_str(), (char*)NULL);

	// Key is not present :
	if(s==(char*)NULL)
		return 2;

	// Key is empty...
	if (s[0] == '\0')
		return 1;

	str=(string)s;
	return 0;
}

int read_common(string &output, dictionary	*ini, struct param_common &dir){

	read_dir(output, ini, dir, "commons:data_directory") || \
			read_dir(output, ini, dir, "commons:input_directory") || \
			read_dir(output, ini, dir, "commons:output_dir") || \
			read_dir(output, ini, dir, "commons:temp_dir") || \
			read_dir(output, ini, dir, "saneInv:noise_dir");

	read_bolo_suffix(ini, dir.suffix);
	read_bolo_global_file(ini, dir.bolo_global_filename);


	if((dir.bolo_global_filename=="") && (dir.suffix=="")){
		output += "You must mention one of those parameters :\nparam_common:suffix or param_common:suffix\n";
		return 2;
	}

	return 0;

}

int read_param_process(string &output, dictionary *ini, struct param_common dir, struct param_sanePre &proc_param, std::vector<double> &fcut){

	int returned=0;

	returned = read_apodize_samples(output, ini, proc_param) + \
			read_nofillgap(output, ini, proc_param)              + \
			read_sampling_frequency(output, ini, proc_param)     + \
			read_filter_frequency(output, ini, proc_param)       + \
			read_baseline(output, ini, proc_param)               + \
			read_correlation(output, ini,proc_param)              + \
			read_remove_poly(output, ini, proc_param);

	read_noise_cut_freq(output, ini, dir, proc_param, fcut);

	returned > 0 ? returned = 1 : returned = 0;
	return returned;

}

int read_param_positions(string &output, dictionary *ini, struct param_sanePos &pos_param){

	string str;

	if(read_item<double>(output, ini, pos_param.pixdeg, 6/3600, "sanePos:pixsize", 1))
		return 2;

	if(pos_param.pixdeg < 0){
		output += "Pixsize cannot be negative ! \n";
		return 1 ;
	}

	// Read the mask_file if present
	if (read_parser_string(ini,"sanePos:mask_file",str)==0)
		pos_param.maskfile=str;

	// Read the file format if present
	// default to SANEPIC file format (with reference position & offsets)
	// 0: sanepic format with reference position & offsets
	// 1: 'hipe' like format with RA/DEC for each time/bolo
	if(read_item<int>(output, ini, pos_param.fileFormat, 0, "sanePos:file_format", 0))
		return 2;

	if((pos_param.fileFormat!=0) && (pos_param.fileFormat!=1)){
		output += "Fileformat must be 0 (SANEPIC) or 1 (HIPE) \n";
		return 1;
	}

	// Read what to do with flagged data : (default : 0 -- map in a single pixel)
	if(read_item<bool>(output, ini, pos_param.flgdupl, 0, "sanePos:map_flagged_data", 0))
		return 2;


	// Read what to do with gaps : (default : 0 -- no projection)
	if(read_item<bool>(output, ini, pos_param.projgaps, 0, "sanePos:project_gaps", 0))
		return 2;

	return 0;
}


int read_param_saneInv(std::string &output, dictionary *ini, struct param_saneInv &saneInv_struct){

	read_cov_matrix_file(output, ini, saneInv_struct.cov_matrix_file);
	read_cov_matrix_suffix(output, ini, saneInv_struct.cov_matrix_suffix);

	if((saneInv_struct.cov_matrix_file=="") && (saneInv_struct.cov_matrix_suffix=="")){
		output += "You must mention one of those parameters :\nsaneInv:cov_matrix_suffix or saneInv:cov_matrix_global_file\n";
		return 2;
	}

	return 0;
}

int read_param_sanePS(std::string &output, dictionary *ini, struct param_sanePS &sanePS_struct){

	if(read_map_file(ini, sanePS_struct.signame) ||
			read_mixmatfile_suffix(ini, sanePS_struct.mix_suffix) ||
			read_ell_suffix(output, ini, sanePS_struct.ell_suffix) ||
			read_ell_global_file(output, ini, sanePS_struct.ell_global_file) ||
			read_fcut(output, ini, sanePS_struct.fcutPS) ||
			read_ncomp(output, ini, sanePS_struct.ncomp) ||
			read_mixmat_global_file(ini, sanePS_struct.mix_global_file))
		return 2;

	if((sanePS_struct.ell_global_file=="") && (sanePS_struct.ell_suffix=="")){
		output += "You must mention one of those parameters :\nsanePS:ell_global_file or sanePS:ell_suffix\n";
		return 2;
	}
	if((sanePS_struct.mix_global_file=="") && (sanePS_struct.mix_suffix=="")){
		output += "You must mention one of those parameters :\nsanePS:mix_global_file or sanePS:mix_suffix\n";
		return 2;
	}

	return 0;
}

int read_param_sanePic(std::string &output, dictionary *ini, struct param_sanePic &sanePic_struct){

	read_iter(output, ini, sanePic_struct.iterw);

	if(sanePic_struct.iterw==0){
		sanePic_struct.save_data=0;
		sanePic_struct.iterw=10;
	}else{
		if(sanePic_struct.iterw<0)
			sanePic_struct.save_data=0;
		else
			sanePic_struct.save_data=1;
	}

	return 0;
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

void print_common(struct param_common dir){

	cout << "Data directory : " << dir.dirfile << "\n";
	cout << "Input directory : " << dir.input_dir << "\n";
	cout << "Temporary directory : " << dir.tmp_dir << "\n";
	cout << "Output directory : " << dir.output_dir << "\n";
	cout << "Noise directory : " << dir.noise_dir << endl;

	cout << endl;

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
		system((char*)make_it.c_str());
		output += "Path : " + strPath + " created\n";
		return 0;
	}

}

int check_dirfile_paths(string &output,string strPath){

	return check_path(output, strPath + "Fourier_data/","Fourier data binaries") || \
			check_path(output, strPath + "Noise_data/","Noise data binaries") || \
			check_path(output, strPath + "Indexes/","Indexes");

}

int read_bolo_suffix(dictionary	*ini, string &suffix){


	string str;

	if (read_parser_string(ini,"commons:bolo_suffix", str)){
		suffix="";
	}else{
		suffix=str;
	}
	return 0;
}

int read_bolo_global_file(dictionary *ini, string &bolo_global_filename){

	string str;

	if (read_parser_string(ini,"commons:bolo_global_file", str)){
		bolo_global_filename="";
	}else{
		bolo_global_filename=str;
	}
	return 0;

}


int read_bolo_gain_global_file(string &output, dictionary *ini, string dir, string &bolo_global_filename){

	string str;

	if (read_parser_string(ini,"saneCheck:bolo_gain_global_file", str)){
		output += "Warning! saneCheck:bolo_global_file field is void or absent in the ini file. Assuming 1 bolo_gain file per scan :\n";
		bolo_global_filename="";
	}else{
		bolo_global_filename=dir + str;
	}
	return 0;


}


void fill_sanePS_struct(struct param_sanePS &structPS, struct samples samples_struct){


	for(long ii=0;ii<samples_struct.ntotscan;ii++){
		if(structPS.mix_global_file!="")
			structPS.mix_names.push_back(structPS.mix_global_file);
		else
			structPS.mix_names.push_back(FitsBasename(samples_struct.fitsvect[ii]) + structPS.mix_suffix);

		if(structPS.ell_global_file!="")
			structPS.ell_names.push_back(structPS.ell_global_file);
		else
			structPS.ell_names.push_back(FitsBasename(samples_struct.fitsvect[ii]) + structPS.ell_suffix);
	}

}

void fill_noisevect_fcut(struct param_common dir, struct samples &samples_str, struct param_saneInv &saneInv_struct, std::vector<double> &fcut){

	if((saneInv_struct.cov_matrix_file!="")){
		samples_str.noisevect.push_back(saneInv_struct.cov_matrix_file);

	}else{
		for(long iframe = 0; iframe < samples_str.ntotscan ; iframe ++)
			samples_str.noisevect.push_back(FitsBasename(samples_str.fits_table[iframe]) + saneInv_struct.cov_matrix_suffix);
	}


	if((int)fcut.size()==0){

		for(long iframe = 0; iframe < (long)samples_str.noisevect.size() ; iframe ++){

			std::vector<string> bolos;
			long nbins;
			double *ell;
			double **Rellth;

			read_CovMatrix(dir.noise_dir + samples_str.noisevect[iframe], bolos, nbins, ell, Rellth);
			fcut.push_back(ell[0]);
			delete [] ell;
			long nBolos=bolos.size();
			free_dmatrix(Rellth, 0, nBolos * nBolos - 1, 0, nbins - 1);

		}
	}

	if((samples_str.noisevect.size() == 1) && (samples_str.ntotscan > 1)){
		// same noise file for all the scans
		(samples_str.noisevect).resize(samples_str.fitsvect.size(),samples_str.noisevect[0]);

		// if only one fcut, extend to all scans
		fcut.resize(samples_str.ntotscan, fcut[0]);
	}


}

int parser_function(char * ini_name, std::string &output, struct param_common &dir,
		struct samples &samples_struct,
		struct param_sanePos &pos_param, struct param_sanePre &proc_param, std::vector<double> &fcut,
		struct param_sanePS &structPS, struct param_saneInv &saneInv_struct, struct param_sanePic &sanePic_struct, int size){

	dictionary	*	ini ;
	string filename;
	struct detectors det;

	// default values :
	proc_param.napod  = 0; /*! number of samples to apodize, =0 -> no apodisation */
	proc_param.NOFILLGAP = 0; /*! dont fill the gaps ? default is NO => the program fill */
	samples_struct.ntotscan=0; /*! total number of scans */
	pos_param.flgdupl = 0; // map duplication factor
	proc_param.fsamp = 0.0;// 25.0; /*! sampling frequency : BLAST Specific*/
	proc_param.NORMLIN = 0; /*!  baseline is removed from the data, NORMLIN = 1 else 0 */
	proc_param.CORRon = 1; /*! correlation included in the analysis (=1), else 0, default 0*/
	proc_param.remove_polynomia = 1; /*! remove a polynomia fitted to the data*/
	proc_param.f_lp = 0.0; // low pass filter frequency
	pos_param.flgdupl = 0; // map duplication factor
	pos_param.maskfile = "";
	structPS.ncomp=1;
	sanePic_struct.iterw=10;
	sanePic_struct.save_data=0;


	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return 2;
	}

	if(read_common(output, ini, dir)==1)
		return 2;


	check_path(output, dir.dirfile, "Data directory");
	check_path(output, dir.input_dir, "Input directory");
	check_path(output, dir.output_dir, "Output directory");
	check_path(output, dir.noise_dir, "Covariance Matrix directory");
	check_path(output, dir.tmp_dir, "Temporary directory");
	check_dirfile_paths(output, dir.tmp_dir);


	if(read_fits_file_list(output, ini, dir,samples_struct)==1)
		return 2;

	samples_struct.ntotscan = (samples_struct.fitsvect).size();

	if (samples_struct.ntotscan == 0) {
		output += "Must provide at least one scan.\n\n";
		return 2;
	}

	if(read_param_positions(output, ini, pos_param))
		return 2;
	if(read_param_process(output, ini, dir, proc_param, fcut))
		return 2;
	if(read_param_sanePS(output, ini, structPS))
		return 2;
	if(read_param_saneInv(output, ini, saneInv_struct))
		return 2;
	if(read_param_sanePic(output, ini, sanePic_struct))
		return 2;

	if(pos_param.maskfile!="")
		pos_param.maskfile = dir.input_dir + pos_param.maskfile;

	iniparser_freedict(ini);

	return 0;
}

void parser_printOut(struct param_common dir, struct samples samples_struct, std::vector<detectors> detector_tab,
		struct param_sanePos pos_param, struct param_sanePre proc_param,
		struct param_sanePS structPS, struct param_sanePic sanePic_struct){

	cout << "\nYou have specified the following options : \n\n";

	print_common(dir);
	cout << endl;
	print_param_process(proc_param);
	print_param_positions(pos_param);

	printf("Number of scans      : %ld\n",samples_struct.ntotscan);
	printf("Number of bolometers : \n");
	for(long iframe=0;iframe<samples_struct.ntotscan;iframe++){
		if(samples_struct.fits_table!=NULL)
			printf("Scan number %ld : %s %ld\n", iframe,(char*)(FitsBasename(samples_struct.fits_table[iframe]).c_str()), detector_tab[iframe].ndet);
		else
			printf("Scan number %ld : %s %ld\n", iframe,(char*)(FitsBasename(samples_struct.fitsvect[iframe]).c_str()), detector_tab[iframe].ndet);
	}


}
