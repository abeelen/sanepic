

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
#include "mpi_architecture_builder.h"
#include "parser_functions.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#include "inputFileIO.h"
#include "parse_saneCheck.h"

using namespace std;


int parse_saneCheck_ini_file(char * ini_name, struct directories &dir,
		struct detectors &det,struct samples &samples_struct, int rank)
{


	dictionary	*	ini ;

	struct param_positions pos_param;
	struct param_process proc_param;
	std::vector<double> fcut;
	double fcut_double;
	string MixMatfile, ellFile, signame;
	long ncomp;
	int iterw;


	string text;
	string str;

	string s;
	ostringstream out;

	ofstream file;


	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	if(read_directories(ini, dir, rank)==1)
		return -1;

	if(read_channel_list(ini,det.boloname, rank)==1)
		return -1;

	if(read_fits_file_list(ini, dir,samples_struct, rank)==1)
		return -1;

	samples_struct.ntotscan = (samples_struct.fitsvect).size();
	det.ndet = (long)((det.boloname).size());

	if (det.ndet == 0) {
		cerr << "Must provide at least one channel.\n\n";
		return -1 ;
	}

	if(read_param_positions(ini, pos_param, rank)==1)
		return -1;

	if(read_param_process(ini, proc_param, rank)==1)
		return -1;

	if(read_noise_cut_freq(ini, fcut,rank)==1)
		return -1;

	if(read_ell_file(ini, ellFile, rank) ||
			read_map_file(ini, signame, rank) ||
			read_mixmatfile(ini, MixMatfile, rank)||
			read_ncomp(ini, ncomp, rank) ||
			read_fcut(ini, fcut_double, rank))
		return -1;



	if(rank==0){

		//		printf("\nsaneCheck parser operations completed :\n");
		cout << "You have specified the following options : \n\n";

		print_directories(dir);
		print_param_process(proc_param);
		print_param_positions(pos_param);


		printf("Number of scans      : %ld\n",samples_struct.ntotscan);
		printf("Number of bolometers : %ld\n",det.ndet);
	}






	text = "# Sanepic ini file :\n";
	text = text + "# Edit the lines to modify the options !\n";
	text = text + "# Some options are required : Check the comments !\n";
	text = text + "\n\n";

	text = text + "[commons]\n\n" + "# Mandatory fields\n\n";
	text = text + "data_directory =   " + dir.dirfile + " ; source data directory\n\n";

	str = "";
	read_parser_string(ini, "commons:channel", rank, str);
	text = text + "channel =   " + str + " ; a .txt file that includes bolometers name\n\n";
	text = text + "output_dir =   " + dir.outdir + " ; output directory\n\n";
	text = text + "temp_dir =   " + dir.tmp_dir +" ; temp directory\n\n";

	str = "";
	read_parser_string(ini, "commons:fits_filelist", rank, str);
	text = text + "fits_filelist =   " + str + " ; a .txt file containing the fits file names, (the noise filenames and processors indexes)\n\n";
	out << proc_param.napod;
	s=out.str();
	text = text + "apodize_Nsamples =   " + s + " ; number of samples to apodize\n\n";

	text = text + "# Optional fields\n\n";
	if(pos_param.flgdupl)
		text = text + "map_flagged_data =   " + "True" + " ; Keyword specifying if flagged data are put in a separated map, (default is False)\n\n";
	else
		text = text + "map_flagged_data =   " + "False" + " ; Keyword specifying if flagged data are put in a separated map, (default is False)\n\n";

	if(proc_param.NOFILLGAP)
		text = text + "nofill_gap =   " + "True" +" ; Do we fill the gaps with white noise + baseline ? F = fill (default), T = don't fill\n\n";
	else
		text = text + "nofill_gap =   " + "False" +" ; Do we fill the gaps with white noise + baseline ? F = fill (default), T = don't fill\n";

	text = text + "\n\n\n";

	text = text + "[sanePos]\n\n";
	text = text + "# Mandatory fields\n\n";
	out.str("");
	out << pos_param.pixdeg;
	s=out.str();
	text = text + "pixsize = " + s + " ; size of pixels (degrees)\n";

	text = text + "\n\n\n";

	text = text + "[sanePre]\n\n";
	text = text + "# Mandatory fields\n\n";
	out.str("");
	out << proc_param.fsamp;
	s=out.str();
	text = text + "sampling_frequency = " + s + " ; detectors sampling frequency (Hz)\n\n";
	out.str("");
	out << proc_param.f_lp;
	s=out.str();
	text = text + "filter_frequency = " + s + " ; frequency of the high pass filter applied to the data\n\n";

	read_parser_string(ini, "sanepic_preprocess:fcut_file", rank, str);
	text = text + "fcut_file = " + str + " ; frequency under which noise power spectra are thresholded. Default is the frequency cut of the high pass filter applied to the data\n\n";

	text = text + "# Optional fields\n\n";
	out.str("");
	out << proc_param.poly_order;
	s=out.str();
	text = text + "poly_order = " + s + " ; Remove a polynomia fitted to the data to reduce fluctuations on timescales larger than the length of the considered segment (default = 0)\n\n";
	if (proc_param.NORMLIN)
		text = text + "no_baseline = " + "True" + " ; Keyword specifiying if a baseline is removed from the data or not (0 for YES (default), 1 for NO)\n\n";
	else
		text = text + "no_baseline = " + "False" + " ; Keyword specifiying if a baseline is removed from the data or not (0 for YES (default), 1 for NO)\n\n";

	if(proc_param.CORRon)
		text = text + "correlation = " + "True" + " ; Set this keyword to False if correlations between detectors are not included in the analysis (default = True)\n\n";
	else
		text = text + "correlation = " + "False" + " ; Set this keyword to False if correlations between detectors are not included in the analysis (default = True)\n\n";

	text = text + "\n\n\n";

	text = text + "[saneInv]\n\n";
	text = text + "# Mandatory fields\n\n";
	text = text + "noise_dir = " + dir.noise_dir + " ; cov matrix directory\n\n";

	text = text + "# Optional fields\n\n";
	str = "";
	read_parser_string(ini, "sanepic_inv_matrix:cov_matrix_file",rank,str);
	text = text + "cov_matrix_file = " + str + " ; this file contains the matrix you want to invert\n\n";

	text = text + "\n\n\n";

	text = text + "[sanePS]\n\n";
	text = text + "# Mandatory fields\n\n";

	text = text + "noise_estim = " + MixMatfile + " ; Enter filename containing the mixing matrix of noise components.\n\n";
	text = text + "ell_file = " + ellFile + " ; file containing the bin for the noise spectrum\n\n";
	out.str("");
	out << ncomp;
	s=out.str();
	text = text + "ncomp = " + s + " ; number of component(s) to estimate\n\n";
	out.str("");
	out << fcut_double;
	s=out.str();
	text = text + "fcut = " + s + " ; freq above which value of the noise will not be estimated\n\n";

	text = text + "# Optional fields\n\n";
	text = text + "map_file = " + signame + " ; fits file containing the map that should be substracted to the data for a second noise estimation step\n\n";

	text = text + "\n\n\n";

	text = text + "[sanePic]\n\n";
	text = text + "# Optional fields\n\n";
	if(pos_param.projgaps)
		text = text + "project_gaps = " + "True" + " ; Keyword specifying if gaps are projected to a pixel in the map, if so gap filling of noise only is performed iteratively. Default is False\n\n";
	else
		text = text + "project_gaps = " + "False" + " ; Keyword specifying if gaps are projected to a pixel in the map, if so gap filling of noise only is performed iteratively. Default is False\n\n";
	if(read_iter(ini, iterw, 0)==-1)
		text = text + "iterW = " + " ; Write temporary map files on disk every iterW number of loop\n\n";
	else{
		out.str("");
		out << iterw;
		s=out.str();
		text = text + "iterW = " + s + " ; Write temporary map files on disk every iterW number of loop\n\n";
	}


	string outfile = dir.outdir + "sanepic_ini_model.txt";
	file.open(outfile.c_str(), ios::out);
	if(!file.is_open()){
		cerr << "File [" << outfile << "] Invalid." << endl;
		return -1;
	}

	file << text;
	cout << "Writing model ini file in : " << outfile << endl;

	file.close();

	iniparser_freedict(ini);
	return 0 ;
}
