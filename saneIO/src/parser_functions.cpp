/*
 * parser_functions.cpp
 *
 *  Created on: 20 oct. 2009
 *      Author: matthieu
 */


#include <iostream>
#include <iomanip>

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#include "inputFileIO.h"
#include "parser_functions.h"
#include "struct_definition.h"
//TODO : Shall we move read_fits_file here ?
#include "mpi_architecture_builder.h"

using namespace std;

int read_dirfile(dictionary	*ini, struct directories &dir, int rank){

	string str;

	if(read_parser_string(ini, "commons:data_directory", rank, str))
		return 1;

	cout << str << endl;

	if (str[str.length()-1] != '/')
		str = str + '/';
	dir.dirfile = str;

	return 0;
}


int read_tmpdir(dictionary	*ini, struct directories &dir, int rank){

	char *pPath;
	string str;

	pPath = getenv ("TMPBATCH");
	if (pPath!=NULL){
		dir.tmp_dir=pPath;
		if(rank==0)
			printf ("The current path is: %s\n",pPath);
	}else{
		if (read_parser_string(ini, "commons:temp_dir",rank,str))
			return 1;
		// TODO: Do we want this ???
		//		if(s==NULL){
		//			if(rank==0){
		//				printf("Warning : The line corresponding to temporary directory in the ini file has been erased : commons:output_dir\n");
		//				cout << "Using default output directory : " << dir.dirfile << endl;
		//			}
		//			dir.tmp_dir=dir.dirfile;
		//		}else{
		if (str[str.length()-1] != '/')
			str = str + '/';
		dir.tmp_dir=str;
	}

	return 0;
}


int read_outdir(dictionary	*ini, struct directories &dir, int rank){

	string str;

	if(read_parser_string(ini, "commons:data_directory", rank, str))
		return 1;

	// TODO: Do we really want this ??
	//	if(s==NULL){
	//		if(rank==0){
	//			printf("Warning : The line corresponding to output directory in the ini file has been erased : commons:output_dir\n");
	//			cout << "Using default output directory : " << dir.dirfile << endl;
	//		}
	//		dir.outdir=dir.dirfile;
	//	}else{
	if (str[str.length()-1] != '/')
		str = str + '/';
	dir.outdir=str;

	return 0;
}


int read_channel_list(dictionary	*ini, std::vector<string> &bolonames, int rank){

	string str;
	//	char *s;

	if(read_parser_string(ini, "commons:channel", rank, str))
		return 1;

	//TODO : Why only rank 0 ?
	if(rank==0)
		read_strings(str, bolonames);

	return 0;

}


int read_fits_file_list(dictionary	*ini, struct directories &dir, struct samples &samples_str, int rank){

	string str;

	if(read_parser_string(ini, "commons:fits_filelist", rank, str))
		return 1;

	samples_str.filename=str;

	// TODO: read_fits_list should return an error in case...
	// Fill fitsvec, noisevect, scans_index with values read from the 'str' filename
	read_fits_list(samples_str.filename, \
			samples_str.fitsvect, samples_str.noisevect, samples_str.scans_index, \
			samples_str.framegiven);


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
		if(read_parser_string(ini, "commons:noise_prefixe",rank,str))
			return 1;
		samples_str.noisevect.push_back(str);
		// meme fichier de bruit pour tous les scans
		(samples_str.noisevect).resize(samples_str.fitsvect.size(),samples_str.noisevect[0]);
	}

	return 0;
}

int read_apodize_samples(dictionary	*ini, struct param_process &proc_param, int rank){

	int i;

	i = iniparser_getint(ini, "commons:apodize_Nsamples", -1);
	proc_param.napod=i;

	if( i<0 && rank==0 ){
		printf("You must choose a positive number of samples to apodize\n");
		return 1 ;
	}

	return 0;

}

int read_nofillgap(dictionary	*ini, struct param_process &proc_param, int rank){

	proc_param.NOFILLGAP = iniparser_getboolean(ini, "commons:nofill_gap", 0);

	return 0;

}

//
//int read_mask_file(dictionary	*ini, std::string maskfile, int rank){
//
//	if (read_parser_string(ini, "positions:mask_file", rank, maskfile))
//		return 1;
//
//	return 0;
//}



int read_sampling_frequency(dictionary	*ini, struct param_process &proc_param, int rank){

	double d;

	d = iniparser_getdouble(ini,(char*)"sanepic_preprocess:sampling_frequency", -1.0);
	if (d<=0.0){
		if(rank==0)
			printf("sampling_frequency cannot be negative or 0 ! Or maybe you forgot to mention sampling frequency \n");
		return 1;
	}else
		proc_param.fsamp=d;

	return 0;
}


int read_filter_frequency(dictionary	*ini, struct param_process &proc_param, int rank){

	double d;

	d = iniparser_getdouble(ini,(char*)"sanepic_preprocess:filter_frequency", -1.0);
	if (d<0.0){
		if(rank==0)
			printf("filter_frequency cannot be negative ! or maybe you have to mention filter frequency \n");
		return 1;
	}else{
		//printf("filter_frequency  :   [%g]\n", d);
		proc_param.f_lp=d;
	}//filter_frequency = 0.005 ;

	return 0;
}


int read_noise_cut_freq(dictionary	*ini, std::vector<double> &fcut, int rank){

	string str;

	if (read_parser_string(ini, "sanepic_preprocess:fcut_file", rank, str))
		return 1;

	std::vector<string> dummy2;
	read_strings(str,dummy2);

	if(((int)dummy2.size())==0){
		if(rank==0)
			printf("You must provide at least one number of noise cut frequency (or one per scan) in fcut_file !\n");
		return 1;}
	for(int ii=0; ii<(int)dummy2.size(); ii++)
		fcut.push_back(atof(dummy2[ii].c_str()));

	return 0;

}

/*
int read_noise_file_list(dictionary	*ini, std::vector<string> &extentnoiseSP){


	string str;
	char *s;

	s = iniparser_getstring(ini, "commons:noise_prefixe_file",NULL);
	if(s==NULL){
		printf("You must add a line corresponding to noise_prefixe file in the ini file : commons:noise_prefixe_file\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		//printf("noise_prefixe_file : [%s]\n",s);
		//extentnoiseSP.push_back(s);

		std::vector<string> dummy3;
		read_strings(str,dummy3);
		if(((int)dummy3.size())==0){
			printf("You must provide at least one noise prefixe (or one per scan) in noise_prefixe_file !\n");
			return -1;}
		for(int ii=0; ii<(int)dummy3.size(); ii++)
			extentnoiseSP.push_back(dummy3[ii].c_str());
	}else{
		printf("you must specify noise_prefixe_file\n");
		return -1;
	}//noise_prefixe = NoisePS ;

	return 0;

}*/

int read_baseline(dictionary	*ini, struct param_process &proc_param, int rank){

	bool b;

	b = iniparser_getboolean(ini, "sanepic_preprocess:no_baseline", 0);
	//if(b!=0){
	//printf("no_baseline:    [%d]\n", b);
	proc_param.NORMLIN=b;
	//}
	//NORMLIN = False ;

	return 0;

}

int read_correlation(dictionary	*ini, struct param_process &proc_param, int rank){

	bool b;

	b = iniparser_getboolean(ini, "sanepic_preprocess:correlation", 1);
	//if(b!=1){
	//printf("correlation:    [%d]\n", b);
	proc_param.CORRon=b;
	//}//CORRon = True

	return 0;
}

int read_remove_poly(dictionary	*ini, struct param_process &proc_param, int rank){

	//	bool b = 1;
	int	i = -1;
	//	b = iniparser_getboolean(ini, "sanepic_preprocess:remove_poly", 1);
	i = iniparser_getint(ini, "sanepic_preprocess:poly_order", -1);
	//if(b!=1){
	//printf("remove_poly:    [%d]\n", b);
	if(i>=0){
		proc_param.remove_polynomia=1;
		proc_param.poly_order=i;
		cout << "poly_order : " << proc_param.poly_order << endl;
	}else{
		proc_param.remove_polynomia=0;
		//		proc_param.poly_order=4;
	}
	//}//remove_poly = True

	return 0;

}


int read_iter(dictionary	*ini, int &iterw, int rank){

	int i;

	i = iniparser_getint(ini, "sanepic_conjugate_gradient:iterW", 0);
	//if(isnan((double)i)){
	if(i>0){
		//printf("iterw :      [%d]\n", i);
		iterw=i;
	}//iterw =  ;

	return 0;
}

int read_ell_file(dictionary	*ini, string &ellFile, int rank){


	string str;
	//char *s;

	if (read_parser_string(ini, "sanepic_estim_PS:ell_file", rank, str))
		return 1;

	ellFile=str;

	return 0;
}


int read_map_file(dictionary	*ini, string &signame, int rank){


	string str;

	if (read_parser_string(ini,"sanepic_estim_PS:map_file", rank,str)){
		signame="NOSIGFILE";
	}else{
		signame=str;
	}
	return 0;
}

int read_cov_matrix_file(dictionary	*ini, string &fname, int rank){

	string str;
	//char*s;

	if (read_parser_string(ini, "sanepic_inv_matrix:cov_matrix_file",rank,str))
		return 1;

	fname=str;
	return 0;

}

int read_mixmatfile(dictionary	*ini, string &MixMatfile, int rank){

	string str;
	if (read_parser_string(ini,"sanepic_estim_PS:noise_estim", rank,str))
		return 1;

	MixMatfile=str;
	return 0;
}

int read_ncomp(dictionary	*ini, long &ncomp, int rank){

	double d;

	d = iniparser_getdouble(ini,(char*)"sanepic_estim_PS:ncomp", -1.0);
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

	d = iniparser_getdouble(ini,(char*)"sanepic_estim_PS:fcut", 12.0);
		if (d<0.0){
			if(rank==0)
				printf("noise cut frequency cannot be negative ! or maybe you have to mention it in the ini file \n");
			return 1;
		}else{
			fcut=d;
		} // default = 12

		return 0;
}

int read_parser_string(dictionary	*ini, string line, int rank, string &str){
	char *s;

	s = iniparser_getstring(ini, line.c_str(), NULL);
	if(s==NULL){
		if(rank==0)
			cout <<"You must add a line in ini file specifying : " << line << endl;
		return 1;
	}
	str=(string)s;
	return 0;
}

int read_directories(dictionary	*ini, struct directories &dir, int rank){


	return read_dirfile(ini, dir, rank) || \
	read_tmpdir(ini, dir, rank)  ||	\
	read_outdir(ini, dir, rank);

}

int read_param_process(dictionary *ini,struct param_process &proc_param, int rank){

	return read_apodize_samples(ini, proc_param, rank) || \
	read_nofillgap(ini, proc_param, rank)              || \
	read_sampling_frequency(ini, proc_param, rank)     || \
	read_filter_frequency(ini, proc_param, rank)       || \
	read_baseline(ini, proc_param, rank)               || \
	read_remove_poly(ini, proc_param, rank);


}

int read_param_positions(dictionary *ini, struct param_positions & pos_param, int rank){

	string str;
	bool b;

	// read the pixelsize
	if (read_parser_string(ini, "sanepic_compute_positions:pixsize", rank,str))
		return 1;
	pos_param.pixdeg=atof(str.c_str());
	if(pos_param.pixdeg < 0 && rank==0){
		printf("Pixsize cannot be negative ! or you forgot to mention pixel size\n");
		return -1 ;
	}

	// Read what to do with flagged data : (default : 0 -- map in a single pixel)
	b = iniparser_getboolean(ini, "commons:map_flagged_data", 0);
	pos_param.flgdupl=b;

	// Read what to do with gaps : (default : 0 -- nogaps)
	b = iniparser_getboolean(ini, "sanepic_conjugate_gradient:project_gaps", 0);
	pos_param.projgaps=b;

	return 0;
}

void print_param_positions(struct param_positions &pos_param) {

	cout << "You have specified a pixel size : [" << setprecision(14) << pos_param.pixdeg << " deg]\n";

	if(pos_param.flgdupl)
		cout << "Flagged data are put in a separate map : map_flagged_data = True\n";

	if(pos_param.projgaps)
		cout << "Gaps are projected to a pixel in the map, gap filling of noise only is performed iteratively\n";
	else
		cout << "Gaps are NOT projected to a pixel in the map (default)\n";

}


void print_param_process(struct param_process & proc_param){


	if(proc_param.napod>0)
		cout << "You have specified a number of samples to apodize : [" << proc_param.napod << "]\n";

	if(proc_param.NOFILLGAP)
		cout << "You have set nofill_gap to True : the gaps in data timeline will NOT be filled\n";
	else
		cout << "You have set nofill_gap to False (default) : the gaps in data timeline will be filled\n";


	cout << "You have chosen a sampling frequency equal to : [" << proc_param.fsamp << " Hz]\n";

	if(proc_param.NORMLIN)
		cout << "No baseline will be removed from the data\n";
	else
		cout << "A baseline will be removed from the data (default)\n";

	if(proc_param.remove_polynomia)
		cout << "Polynomia order :  " << proc_param.poly_order << endl;
	else
		cout << "No polynomia will be used\n";

	if(proc_param.CORRon)
		cout << "Correlations between detectors are included in the analysis\n";
	else
		cout << "Correlations between detectors are NOT included in the analysis\n";

	cout << "Frequency of the high pass filter applied to the data : [" << proc_param.f_lp << " Hz]\n";


}

void print_directories(struct directories dir){

	cout << "Data directory : [" << dir.dirfile << "]\n";
	cout << "Temporary directory : [" << dir.tmp_dir << "]\n";
	cout << "Output directory : [" << dir.outdir << "]\n";

}

