#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
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



int read_dir(dictionary	*ini, struct common &dir, string dirtype ,int rank){

	string str;

	switch(read_parser_string(ini, dirtype, str)){
	case 2:
		if(rank==0){
			cout << "WARNING ! You must add a line in ini file specifying : " << dirtype << endl;
			cout << "Using default directory : ./" << endl;
		}
		str="./";
		break;
		//			return 1;
	case 1:
		if(rank==0){
			cout <<"Key is empty : You must specify : commons:data_directory" << endl;
			cout << "Using default directory : ./" << endl;
		}
		//			return 1;
		str="./";
		break;
	case 0:
		if (str[str.length()-1] != '/')
			str = str + '/';
		break;
	}


	//	unsigned int chkk=(checksum((char*)dirtype.c_str(), (size_t)dirtype.size(), 0));
	//	cout << "chkk : " << chkk << endl;

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

	case (unsigned int)1738: 	//saneInv:noise_dir
	dir.noise_dir=str;
	break;

	case (unsigned int)2458: 	//commons:input_dir
	dir.input_dir=str;
	break;

	}


	return 0;


}

//int read_dirfile(dictionary	*ini, struct common &dir, int rank){
//
//	string str;
//
//	switch(read_parser_string(ini, "commons:data_directory", str)){
//	case 2:
//		if(rank==0)
//			cout <<"You must add a line in ini file specifying : commons:data_directory" << endl;
//		return 1;
//	case 1:
//		if(rank==0)
//			cout <<"Key is empty : You must specify : commons:data_directory" << endl;
//		return 1;
//	case 0:
//		if (str[str.length()-1] != '/')
//			str = str + '/';
//		dir.dirfile = str;
//
//	}
//	return 0;
//}
//
//
//int read_tmpdir(dictionary	*ini, struct common &dir, int rank){
//
//	char *pPath;
//	string str;
//
//	pPath = getenv ("TMPBATCH");
//	if (pPath!=NULL){
//		dir.tmp_dir=pPath;
//	}else{
//		switch(read_parser_string(ini, "commons:temp_dir", str)){
//		case 2:
//			if(rank==0)
//				cout <<"You must add a line in ini file specifying : commons:temp_dir" << endl;
//			return 1;
//		case 1:
//			if(rank==0)
//				cout <<"Key is empty : You must specify : commons:temp_dir" << endl;
//			return 1;
//		case 0:
//			if (str[str.length()-1] != '/')
//				str = str + '/';
//			dir.tmp_dir=str;
//
//		}
//	}
//	return 0;
//
//}
//
//
//int read_outdir(dictionary	*ini, struct common &dir, int rank){
//
//	string str;
//
//
//	switch(read_parser_string(ini, "commons:output_dir", str)){
//	case 2:
//		if(rank==0)
//			cout <<"You must add a line in ini file specifying : commons:output_dir" << endl;
//		return 1;
//	case 1:
//		if(rank==0)
//			cout <<"Key is empty : You must specify : commons:output_dir" << endl;
//		return 1;
//	case 0:
//		if (str[str.length()-1] != '/')
//			str = str + '/';
//		dir.output_dir=str;
//	}
//	return 0;
//
//}
//
//int read_noisedir(dictionary	*ini, struct common &dir, int rank){
//
//	string str;
//	switch(read_parser_string(ini, "saneInv:noise_dir", str)){
//	case 2:
//		if(rank==0)
//			cout <<"You must add a line in ini file specifying : saneInv:noise_dir" << endl;
//		return 1;
//	case 1:
//		if(rank==0)
//			cout <<"Key is empty : You must specify : saneInv:noise_dir" << endl;
//		return 1;
//	case 0:
//		if (str[str.length()-1] != '/')
//			str = str + '/';
//		dir.noise_dir=str;
//	}
//	return 0;
//
//}

//int read_channel_list(dictionary	*ini, struct common &dir, std::vector<string> &bolonames, int rank){
//
//	string str;
//
//
//	switch(read_parser_string(ini, "commons:channel", str)){
//	case 2:
//		if(rank==0)
//			cout <<"You must add a line in ini file specifying : commons:channel" << endl;
//		return 1;
//	case 1:
//		if(rank==0)
//			cout <<"Key is empty : You must specify : commons:channel" << endl;
//		return 1;
//	case 0:
//		dir.channel=str;
//		if(read_strings(str, bolonames))
//			return 1;
//	}
//
//
//	return 0;
//
//}


int read_channel_list(std::string fname, std::vector<string> &bolonames, int rank){

	if(read_strings(fname,bolonames)){
		if(rank==0)
			cout <<"You must create file specifying bolometer list named " << fname << endl;
		return 1;
	}
	return 0;
}


int read_fits_file_list(dictionary	*ini, struct common &dir, struct samples &samples_str, int rank){

	string str;


	switch(read_parser_string(ini, "commons:fits_filelist", str)){
	case 2:
		if(rank==0)
			cout <<"You must add a line in ini file specifying : commons:fits_filelist" << endl;
		return 1;
	case 1:
		if(rank==0)
			cout <<"Key is empty : You must specify : commons:fits_filelist" << endl;
		return 1;
	case 0:
		samples_str.filename=dir.input_dir + str;

		// Fill fitsvec, noisevect, scans_index with values read from the 'str' filename
		if(read_fits_list(samples_str.filename, \
				samples_str.fitsvect, samples_str.noisevect, samples_str.scans_index, \
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


		// read the possible noise file from the ini file

		// if no noise file is given in fits_filelist
		if((int)(samples_str.noisevect).size()==0 ){

			// read the noise file name in the ini file

			switch(read_parser_string(ini, "saneInv:cov_matrix_file", str)){
			case 2:
				if(rank==0)
					cout <<"You must add a line in ini file specifying : saneInv:cov_matrix_file" << endl;
				return 1;
			case 1:
				if(rank==0)
					cout <<"Key is empty : You must specify : saneInv:cov_matrix_file" << endl;
				return 1;
			case 0:
				samples_str.cov_matrix_file = str;

				samples_str.noisevect.push_back(str);
			}

			if((int)((samples_str.fitsvect).size())>1)
				// meme fichier de bruit pour tous les scans
				(samples_str.noisevect).resize(samples_str.fitsvect.size(),samples_str.noisevect[0]);

		}
	}
	return 0;
}

int read_fits_list(string fname, std::vector<string> &fitsfiles, std::vector<string> &noisefiles, std::vector<int> &frameorder, bool &framegiven) {



	ifstream file;
	file.open(fname.c_str(), ios::in);
	if(!file.is_open()){
		cerr << "File [" << fname << "] Invalid." << endl;
		exit(-1);
	}

	framegiven=0;

	string s, p, line, temp;
	int d;
	char *pch;
	int nb_elem = 0;

	// count number of elements on the first line !
	getline(file, line);
	line.erase(0, line.find_first_not_of(" \t")); // remove leading white space
	pch = strtok ((char*) line.c_str()," ,;\t");

	while (pch != NULL) {
		pch = strtok (NULL, " ,;\t");
		nb_elem++; 	}

	// set pointer back to the beginning of file in order to parse the first line too
	file.seekg (0, ios::beg);

	switch(nb_elem) {
	case 3:
		framegiven=1;
		while(file >> s >> p >> d) {
			size_t found;
			s.erase(0, s.find_first_not_of(" \t")); // remove leading white space in the first name
			found = s.find_first_of("!#;"); 		// Check for comment character at the beginning of the filename
			if (found == 0) continue;

			fitsfiles.push_back(s);
			noisefiles.push_back(p);
			frameorder.push_back(d);
		}
		break;

	case 2:
		while(file >> s >> p){
			size_t found;
			s.erase(0, s.find_first_not_of(" \t")); // remove leading white space in the first name
			found = s.find_first_of("!#;"); 		// Check for comment character at the beginning of the filename

			if (found == 0) continue;

			//			cout << "2 : " << s << " " << p << endl;
			fitsfiles.push_back(s);
			noisefiles.push_back(p);}
		break;

	case 1:
		while(file >> s){
			size_t found;
			s.erase(0, s.find_first_not_of(" \t")); // remove leading white space in the first name
			found = s.find_first_of("!#;"); 		// Check for comment character at the beginning of the filename
			if (found == 0) continue;

			//			cout << "1 : " << s << endl;
			fitsfiles.push_back(s);}
		break;

	default:
		cerr << "File [" << fname << "] must have at least one row and 2 colums. Exiting\n";
		return 1;
		break;
	}

	if(fitsfiles.size()==0){
		cerr << "File [" << fname << "] must have at least one row. Exiting\n";
		return 1;
	}

	if (file>>s){
		cerr << "File [" << fname << "]. Each line must have the same number of rows. Exiting\n";
		return 1;
	}


#ifdef DEBUG_PRINT
	cout << "read fits list ok !!!\n";
#endif

	file.close();

	return 0;
}


int read_apodize_samples(dictionary	*ini, struct param_process &proc_param, int rank){

	int i;

	i = iniparser_getint(ini, "sanePre:apodize_Nsamples", -1);
	proc_param.napod=i;

	if( i<0 && rank==0 ){
		printf("You must choose a positive number of samples to apodize\n");
		return 1 ;
	}

	return 0;

}

int read_nofillgap(dictionary	*ini, struct param_process &proc_param, int rank){

	proc_param.NOFILLGAP = iniparser_getboolean(ini, "sanePre:nofill_gap", 0);

	return 0;

}

int read_sampling_frequency(dictionary	*ini, struct param_process &proc_param, int rank){

	double d;

	d = iniparser_getdouble(ini,(char*)"sanePre:sampling_frequency", -1.0);
	if (d<=0.0){
		if(rank==0)
			printf("sampling_frequency cannot be negative or 0 ! Or maybe you forgot to mention sampling frequency \n");
		return 1;
	}else
		proc_param.fsamp=d;

	return 0;
}


int read_filter_frequency(dictionary *ini, struct param_process &proc_param, int rank){

	double d;

	d = iniparser_getdouble(ini,(char*)"sanePre:filter_frequency", -1.0);
	if (d<0.0){
		if(rank==0)
			printf("filter_frequency cannot be negative ! or maybe you have to mention filter frequency \n");
		return 1;
	}else{
		proc_param.f_lp=d;
	}
	return 0;
}


int read_noise_cut_freq(dictionary	*ini, struct common dir, struct param_process &proc_param, std::vector<double> &fcut, int rank){

	string str;


	switch(read_parser_string(ini, "sanePre:fcut_file", str)){
	case 2:
		if(rank==0)
			cout <<"Warning ! You must add a line in ini file specifying : sanePre:fcut_file" << endl;
		cout << "Using min covariance frequency as default\n";
		//return 1;
		break;
	case 1:
		if(rank==0)
			cout <<"Key is empty : You must specify : sanePre:fcut_file" << endl;
		return 1;
	case 0:
		proc_param.fcut_file = dir.input_dir + str;

		std::vector<string> dummy2;
		if(read_strings(proc_param.fcut_file,dummy2))
			return 1;

		if(((int)dummy2.size())==0){
			if(rank==0)
				printf("You must provide at least one number of noise cut frequency (or one per scan) in fcut_file !\n");
			return 1;}
		for(int ii=0; ii<(int)dummy2.size(); ii++)
			fcut.push_back(atof(dummy2[ii].c_str()));
	}
	return 0;

}


int read_baseline(dictionary	*ini, struct param_process &proc_param, int rank){

	bool b;

	b = iniparser_getboolean(ini, "sanePre:no_baseline", 0);
	proc_param.NORMLIN=b;

	return 0;

}

int read_correlation(dictionary	*ini, struct param_process &proc_param, int rank){

	bool b;

	b = iniparser_getboolean(ini, "sanePre:correlation", 1);
	proc_param.CORRon=b;

	return 0;
}

int read_remove_poly(dictionary	*ini, struct param_process &proc_param, int rank){

	int	i = -1;
	i = iniparser_getint(ini, "sanePre:poly_order", -1);
	if(i>=0){
		proc_param.remove_polynomia=1;
		proc_param.poly_order=i;
	}else{
		proc_param.remove_polynomia=0;
	}

	return 0;

}


int read_iter(dictionary	*ini, int &iterw, int rank){

	int i;

	i = iniparser_getint(ini, "sanePic:iterW", 0);
	iterw=i;
	//iterw =  ;

	return 0;
}

int read_ell_suffix(dictionary	*ini, string &ell_suffix, int rank){


	string str;


	switch(read_parser_string(ini, "sanePS:ell_suffix", str)){
	case 2:
		if(rank==0)
			cout <<"You must add a line in ini file specifying : sanePS:ell_suffix" << endl;
		ell_suffix="";
		//		return 1;
	case 1:
		if(rank==0)
			cout <<"Key is empty : You must specify : sanePS:ell_suffix" << endl;
		ell_suffix="";
		//		return 1;
	case 0:
		ell_suffix=str;
	}

	return 0;
}

int read_ell_global_file(dictionary	*ini, string &ell_global_file, int rank){


	string str;


	switch(read_parser_string(ini, "sanePS:ell_global_file", str)){
	case 2:
		if(rank==0)
			cout <<"You must add a line in ini file specifying : sanePS:ell_global_file" << endl;
		ell_global_file="";
		//		return 1;
	case 1:
		if(rank==0)
			cout <<"Key is empty : You must specify : sanePS:ell_global_file" << endl;
		ell_global_file="";
		//		return 1;
	case 0:
		ell_global_file=str;
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

int read_cov_matrix_file(dictionary	*ini, string &fname, int rank){

	string str;

	switch(read_parser_string(ini, "saneInv:cov_matrix_file", str)){
	case 2:
		if(rank==0)
			cout <<"You must add a line in ini file specifying : saneInv:cov_matrix_file" << endl;
		return 1;
	case 1:
		if(rank==0)
			cout <<"Key is empty : You must specify : saneInv:cov_matrix_file" << endl;
		return 1;
	case 0:
		fname=str;
	}
	return 0;

}

int read_mixmatfile_suffix(dictionary	*ini, string &MixMat_suffix, int rank){

	string str;

	if (read_parser_string(ini,"sanePS:MixingMatrix_Suffix", str)){
		MixMat_suffix = "";
	}else{
		MixMat_suffix = str;
	}

	return 0;
}

int read_mixmat_global_file(dictionary	*ini, string &MixMat_global, int rank){

	string str;

	if (read_parser_string(ini,"sanePS:MixingMatrix_global_file", str)){
		MixMat_global = "";
	}else{
		MixMat_global = str;
	}

	return 0;
}

int read_ncomp(dictionary	*ini, long &ncomp, int rank){

	double d;

	d = iniparser_getdouble(ini,(char*)"sanePS:ncomp", -1.0);
	if (d<0.0){
		if(rank==0)
			printf("number of component cannot be negative ! or maybe you have to mention it in the ini file \n");
		return 1;
	}else{

		ncomp=(long)d;
	}

	return 0;
}


int read_fcut(dictionary	*ini, double &fcut, int rank){

	double d;

	d = iniparser_getdouble(ini,(char*)"sanePS:fcut", 12.0);
	if (d<0.0){
		if(rank==0)
			printf("noise cut frequency cannot be negative ! or maybe you have to mention it in the ini file \n");
		return 1;
	}else{
		fcut=d;
	} // default = 12

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

int read_common(dictionary	*ini, struct common &dir, int rank){



	return read_dir(ini, dir, "commons:data_directory" , rank) || \
			read_dir(ini, dir, "commons:input_directory" , rank) || \
			read_dir(ini, dir, "commons:output_dir" , rank) || \
			read_dir(ini, dir, "commons:temp_dir" , rank) || \
			read_dir(ini, dir, "saneInv:noise_dir" , rank);


}

// when a function is wrong, it crashes the reading procedure and the other parameters are not read so I changed || to +
int read_param_process(dictionary *ini,struct param_process &proc_param, int rank){

	int returned=0;

	returned = read_apodize_samples(ini, proc_param, rank) + \
			read_nofillgap(ini, proc_param, rank)              + \
			read_sampling_frequency(ini, proc_param, rank)     + \
			read_filter_frequency(ini, proc_param, rank)       + \
			read_baseline(ini, proc_param, rank)               + \
			read_correlation(ini,proc_param,rank)              +\
			read_remove_poly(ini, proc_param, rank);

	returned > 0 ? returned = 1 : returned = 0;
	return returned;

}

int read_param_positions(dictionary *ini, struct param_positions &pos_param, int rank){

	string str;
	bool b;

	// read the pixelsize
	if (read_parser_string(ini, "sanePos:pixsize", str)==0)
		pos_param.pixdeg=atof(str.c_str());
	else
		return 1;

	// Read the mask_file if present
	if (read_parser_string(ini,"sanePos:mask_file",str)==0)
		pos_param.maskfile=str;

	// Read the file format if present
	// default to SANEPIC file format (with reference position & offsets)
	// 0: sanepic format with reference position & offsets
	// 1: 'hipe' like format with RA/DEC for each time/bolo
	pos_param.fileFormat = iniparser_getint(ini, "sanePos:file_format", 0);

	if(pos_param.pixdeg < 0){
		if(rank==0)
			printf("Pixsize cannot be negative ! or you forgot to mention pixel size\n");
		return 1 ;
	}

	// Read what to do with flagged data : (default : 0 -- map in a single pixel)
	b = (bool)iniparser_getboolean(ini, "sanePos:map_flagged_data", (bool)0);
	pos_param.flgdupl=b;


	// Read what to do with gaps : (default : 0 -- no projection)
	b = (bool)iniparser_getboolean(ini, "sanePos:project_gaps", (bool)0);
	pos_param.projgaps=b;


	return 0;
}

void print_param_positions(struct param_positions pos_param) {

	cout << "Pixel size : " << setprecision(14) << pos_param.pixdeg << " deg\n";

	if(pos_param.flgdupl)
		cout << "Separate map : " << setw(10) << "map_flagged_data = True\n";


	if(pos_param.projgaps)
		cout << "GAP FILLING : " << setw(10) << "PROJECTED\n";
	else
		cout << "GAP FILLING : " << setw(10) << "NOT projected (default)\n";

	cout << endl;

}


void print_param_process(struct param_process proc_param){




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

void print_common(struct common dir){

	cout << "Data directory : " << dir.dirfile << "\n";
	cout << "Input directory : " << dir.input_dir << "\n";
	cout << "Temporary directory : " << dir.tmp_dir << "\n";
	cout << "Output directory : " << dir.output_dir << "\n";
	cout << "Noise directory : " << dir.noise_dir << endl;

	cout << endl;

}

int check_path(string strPath, string path_type){



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
			cout << "Warning : The path " << path_type << " : " << strPath << " is a file." << endl;
			return 1;
		}
	}
	else
	{
		cout << "Warning : Path " << path_type << " : " << strPath << " doesn't exist." << endl;
		string make_it = "mkdir " + strPath;
		system((char*)make_it.c_str());
		cout << "Path : " << strPath << " created" << endl;
		return 0;
	}



}

int check_dirfile_paths(string strPath){

	return check_path(strPath + "Fourier_data/","Fourier data binaries") || \
			check_path(strPath + "Noise_data/","Noise data binaries") || \
			check_path(strPath + "Indexes/","Indexes");

}

//int read_save_data(dictionary *ini, int &save_data, int rank){
//
//	string str;
//
//	if (read_parser_string(ini, "sanePic:save_data", str)==0){
//		save_data=atoi(str.c_str());
//	}else{
//		cout << "Please mention sanePic:save_data !\n";
//		return 1;
//	}
//
//	return 0;
//
//}

int read_restore(dictionary *ini, int &restore, int rank){

	string str;

	if (read_parser_string(ini, "sanePic:restore", str)==0){
		restore=atoi(str.c_str());
	}else{
		cout << "Please mention sanePic:restore !\n";
		return 1;
	}

	return 0;

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

void fill_sanePS_struct(std::string dir, struct PS &structPS, struct samples samples_struct){


	for(long ii=0;ii<samples_struct.ntotscan;ii++){
		if(structPS.mix_global_file!="")
			structPS.mix_names.push_back(dir + structPS.mix_global_file);
		else
			structPS.mix_names.push_back(dir + FitsBasename(samples_struct.fitsvect[ii]) + structPS.mix_suffix);

		if(structPS.ell_global_file!="")
			structPS.ell_names.push_back(dir + structPS.ell_global_file);
		else
			structPS.ell_names.push_back(dir + FitsBasename(samples_struct.fitsvect[ii]) + structPS.ell_suffix);
	}

}

int parser_function(char * ini_name, struct common &dir,
		std::vector<detectors> &detector_tab,struct samples &samples_struct,
		struct param_positions &pos_param, struct param_process &proc_param, std::vector<double> &fcut,
		struct PS &structPS, struct sanePic &sanePic_struct, int rank, int size){
	//		double &fcut_sanePS, string &MixMatSuffix, string &ell_suffix, string &signame, long &ncomp, int &iterw,
	//		int &save_data, int &restore, int rank, int size){

	dictionary	*	ini ;
	string filename;
	struct detectors det;
	string suffix, bolo_global_filename;

	// default values :
	proc_param.napod  = 0; /*! number of samples to apodize, =0 -> no apodisation */
	proc_param.NOFILLGAP = 0; /*! dont fill the gaps ? default is NO => the program fill */
	samples_struct.ntotscan=0; /*! total number of scans */
	//	det.ndet=0; /*! number of channels used*/
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
	sanePic_struct.restore=0;


	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return 2;
	}

	if(read_common(ini, dir, 0)==1)
		return 2;

	//	if(read_ell_dir(ini, dir.ell_path, rank))
	//		return 2;

	if(rank==0){
		check_path(dir.dirfile, "Data directory");
		check_path(dir.input_dir, "Input directory");
		check_path(dir.output_dir, "Output directory");
		check_path(dir.noise_dir, "Covariance Matrix directory");
		check_path(dir.tmp_dir, "Temporary directory");
		check_dirfile_paths(dir.tmp_dir);
	}

	if(read_fits_file_list(ini, dir,samples_struct, rank)==1)
		return 2;

	samples_struct.ntotscan = (samples_struct.fitsvect).size();

	if (samples_struct.ntotscan == 0) {
		if(rank==0)
			cerr << "Must provide at least one scan.\n\n";
		return 2;
	}

	read_bolo_suffix(ini, suffix);

	read_bolo_global_file(ini, bolo_global_filename);

	for(long oo=0;oo<samples_struct.ntotscan;oo++){

		if(bolo_global_filename!="")
			filename=dir.input_dir + bolo_global_filename;
		else
			filename= dir.dirfile + FitsBasename(samples_struct.fitsvect[oo]) + suffix ; //  + ".bolo"


		//		cout << filename << endl;
		if(read_channel_list(filename, det.boloname, rank)==1)
			return 2;
		det.ndet = (long)((det.boloname).size());
		if (det.ndet == 0) {
			if(rank==0)
				cerr << "Must provide at least one channel.\n\n";
			return 2;
		}
		if(size>det.ndet){
			cerr << "You are using too many processors : " << size << " processors for only " << detector_tab[oo].ndet << " detectors!" << endl;
			return 3;
		}
		detector_tab.push_back(det);
		det.ndet=0;
		det.boloname.clear();

	}
	//	cout << detector_tab[0].ndet << endl;


	if(	read_param_positions(ini, pos_param, rank) ||
			read_param_process(ini, proc_param, rank) ||
			read_map_file(ini, structPS.signame) ||
			read_mixmatfile_suffix(ini, structPS.mix_suffix, rank) ||
			read_ell_suffix(ini, structPS.ell_suffix, rank) ||
			read_ell_global_file(ini, structPS.ell_global_file, rank) ||
			read_fcut(ini, structPS.fcutPS, rank) ||
			read_ncomp(ini, structPS.ncomp, rank) ||
			read_mixmat_global_file(ini, structPS.mix_global_file, rank) ||
			read_restore(ini, sanePic_struct.restore, rank))
		return 2;

	//	cout << "parser save data : " << save_data << endl;
	//	cout << "parser load data : " << restore << endl;

	read_iter(ini, sanePic_struct.iterw, rank);

	if(sanePic_struct.iterw==0){
		sanePic_struct.save_data=0;
		sanePic_struct.iterw=10;
	}else{
		if(sanePic_struct.iterw<0)
			sanePic_struct.save_data=0;
		else
			sanePic_struct.save_data=1;
	}


	read_noise_cut_freq(ini, dir, proc_param, fcut,rank);

	if((int)fcut.size()==0){
		string fname;
		std::vector<string> bolos;
		long nbins;
		double *ell;
		double **Rellth;

		read_cov_matrix_file(ini, fname, rank);
		fname = dir.noise_dir + fname;
		read_CovMatrix(fname, bolos, nbins, ell, Rellth);
		fcut.push_back(ell[0]);
		delete [] ell;
		long nBolos=bolos.size();
		free_dmatrix(Rellth, 0, nBolos * nBolos - 1, 0, nbins - 1);
	}


	iniparser_freedict(ini);

	// if only one fcut, extend to all scans
	if((int)fcut.size()==1)
		fcut.resize(samples_struct.ntotscan, fcut[0]);


	if(rank==0){
		cout << "\nYou have specified the following options : \n\n";

		print_common(dir);
		cout << endl;
		print_param_process(proc_param);
		print_param_positions(pos_param);

		cout << "sanePic save data : " << sanePic_struct.save_data << endl;
		cout << "sanePic restore : " << sanePic_struct.restore << endl;

		printf("Number of scans      : %ld\n",samples_struct.ntotscan);
		printf("Number of bolometers : \n");
		for(long iframe=0;iframe<samples_struct.ntotscan;iframe++)
			printf("Scan number %ld : %ld\n", iframe, detector_tab[iframe].ndet);
	}


	return 0;
}


