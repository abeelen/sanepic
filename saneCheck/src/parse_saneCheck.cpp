

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <unistd.h>
#include <string>
#include <vector>
#include <algorithm>

#include "dataIO.h"
#include "struct_definition.h"
#include "mpi_architecture_builder.h"
#include "parser_functions.h"
#include "tools.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#include "inputFileIO.h"
#include "parse_saneCheck.h"

using namespace std;


template <class T>
std::string StringOf(const T& object){
	std::ostringstream os;
	os << object;
	return os.str();
}

int parse_saneCheck_ini_file(char * ini_name, struct common &dir,
		std::vector<detectors> &detector_tab,struct samples &samples_struct, double &fsamp, int rank)
{


	dictionary	*	ini ;

	struct param_positions pos_param;
	struct param_process proc_param;
	std::vector<double> fcut;
	double fcut_double;
	string MixMatfile, ellFile, signame;
	long ncomp=1;
	int iterw=10;


	string text;
	string str;
	string filename;
	struct detectors det;

	string s;
	ofstream file;


	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	if(read_common(ini, dir, rank)==1)
		return -1;

	check_path(dir.dirfile, "Input directory");
	check_path(dir.output_dir, "Output directory");
	check_path(dir.noise_dir, "Covariance Matrix directory");
	check_path(dir.tmp_dir, "Temporary directory");
	check_dirfile_paths(dir.tmp_dir);


	if(read_fits_file_list(ini, dir,samples_struct, rank)==1)
		return -1;

	samples_struct.ntotscan = (samples_struct.fitsvect).size();


	for(long oo=0;oo<samples_struct.ntotscan;oo++){
		filename= dir.dirfile + FitsBasename(samples_struct.fitsvect[oo]) + ".bolo";
		//		cout << filename << endl;
		if(read_channel_list(filename, det.boloname, rank)==1)
			return -1;
		det.ndet = (long)((det.boloname).size());
		if (det.ndet == 0) {
			if(rank==0)
				cerr << "Must provide at least one channel.\n\n";
			return -1;
		}
		detector_tab.push_back(det);
		det.ndet=0;
		det.boloname.clear();

	}

	if(read_iter(ini, iterw, rank)==-1)
		return -1;

	read_param_positions(ini, pos_param, rank);

	read_param_process(ini, proc_param, rank);


	read_noise_cut_freq(ini, proc_param, fcut,rank);


	read_ell_file(ini, ellFile, rank);
	read_map_file(ini, signame);
	read_mixmatfile(ini, MixMatfile, rank);
	read_ncomp(ini, ncomp, rank);
	read_fcut(ini, fcut_double, rank);

	fsamp=proc_param.fsamp;



	if(rank==0){
		cout << "\nYou have specified the following options : \n\n";


		print_common(dir);
		if(check_path(dir.output_dir, "output directory"))
			return -1;
		check_path( dir.tmp_dir, "Temporary directory");
		check_path( dir.noise_dir, "Covariance Matrix directory");
		if(check_path( dir.dirfile, "Input directory"))
			return -1;
		cout << endl;
		print_param_process(proc_param);
		print_param_positions(pos_param);

		printf("Number of scans      : %ld\n",samples_struct.ntotscan);
		printf("Number of bolometers : \n");
		for(long iframe=0;iframe<samples_struct.ntotscan;iframe++)
			printf("Scan number %ld : %ld\n", iframe, detector_tab[iframe].ndet);




		// Generates ini file model
		text = "# Sanepic ini file :\n";
		text += "# Edit the lines to modify the options !\n";
		text += "# Some options are required : Check the comments !\n";

		text += "\n\n";

		text += "[commons]\n\n";
		text += "data_directory = " + dir.dirfile + " ; source data directory\n";
		text += "channel = " + dir.channel + " ; file listing bolometers name\n";
		text += "output_dir = " + dir.output_dir + " ; output directory\n";
		text += "temp_dir = " + dir.tmp_dir +" ; temporary directory\n";
		text += "fits_filelist = " + samples_struct.filename + " ; file containing fits file names, [corresponding noise file, [processors indexes]]\n";

		text += "\n\n";

		text += "[sanePos]\n\n";
		text += "pixsize = " + StringOf(pos_param.pixdeg) + " ; size of pixels (degrees)\n";
		text += "map_flagged_data = " + StringOf( pos_param.flgdupl ? "True" : "False" ) + " ; flagged data put in a separated map (default is False)\n";
		text += "file_format = " + StringOf(pos_param.fileFormat) + " ; HIPE = 1, SANEPIC = 2\n";
		text += "project_gaps = " + StringOf(pos_param.projgaps ? "True" : "False" ) + "; Keyword specifying if gaps are projected to a pixel in the map, if so gap filling of noise only is performed iteratively. Default is False\n";
		text += "mask_file = " + pos_param.maskfile + "; mask (fits file) : use to remove cross constraint due to strong sources\n";

		text += "\n\n";

		text += "[sanePre]\n\n";
		text += "sampling_frequency = " + StringOf(proc_param.fsamp) + " ; detectors sampling frequency (Hz)\n";
		text += "filter_frequency = " + StringOf(proc_param.f_lp) + " ; frequency of the high pass filter applied to the data\n";
		text += "apodize_Nsamples = " + StringOf(proc_param.napod) + " ; number of samples to apodize\n";
		text += "fcut_file = " + proc_param.fcut_file + " ; noise power spectra are thresholded. Default is the frequency cut of the high pass filter applied to the data\n";
		text += "poly_order = " + StringOf(proc_param.poly_order) + " ; baseline polynomia order (default = 0; no baseline)\n";
		text += "no_baseline = " + StringOf(proc_param.NORMLIN ? "True" : "False") + " ; no simple baseline removed from the data (False : default)\n";
		text += "correlation = " + StringOf(proc_param.CORRon ? "True" : "False") + " ; Set this keyword to False if correlations between detectors are not included in the analysis (default = True)\n";
		text += "nofill_gap = " + StringOf(proc_param.NOFILLGAP ? "True" : "False") +" ; Do we fill the gaps ? F = fill (default), T = don't fill\n";

		text += "\n\n";

		text += "[saneInv]\n\n";
		text += "noise_dir = " + dir.noise_dir + " ; cov matrix directory\n";
		text += "cov_matrix_file = " + samples_struct.cov_matrix_file + " ; this file contains the matrix you want to invert\n";

		text += "\n\n";

		text += "[sanePS]\n\n";
		text += "noise_estim = " + MixMatfile + " ; Enter filename containing the mixing matrix of noise components.\n";
		text += "ell_file = " + ellFile + " ; file containing the bin for the noise spectrum\n";
		text += "ncomp = " + StringOf(ncomp) + " ; number of component(s) to estimate\n";
		text += "fcut = " + StringOf(fcut_double) + " ; freq above which value of the noise will not be estimated\n";
		text += "map_file = " + signame + " ; fits file containing the map that should be substracted to the data for a second noise estimation step\n\n";

		text += "\n\n";

		text += "[sanePic]\n\n";

		if(read_iter(ini, iterw, 0)==-1)
			text += "iterW =  ; Write temporary map files on disk every iterW number of loop\n";
		else{
			text += "iterW = " + StringOf(iterw) + " ; Write temporary map files on disk every iterW number of loop\n";
		}

		string outfile = dir.output_dir + "sanepic_ini_model.txt";
		file.open(outfile.c_str(), ios::out);
		if(!file.is_open()){
			cerr << "File [" << outfile << "] Invalid." << endl;
			return -1;
		}

		file << text;
		cout << "Writing model ini file in : " << outfile << endl;

		file.close();
	}

	iniparser_freedict(ini);
	return 0;
}
