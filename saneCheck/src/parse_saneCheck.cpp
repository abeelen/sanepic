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
#include "utilities.h"
#include "tools.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#include "inputFileIO.h"
#include "parse_saneCheck.h"

using namespace std;

uint16_t parse_saneCheck_ini_file(char * ini_name, string &output, struct param_common &dir,
		struct samples &samples_struct, struct param_sanePos &pos_param, struct param_sanePre &proc_param,
		struct param_sanePS &structPS, struct param_saneInv &saneInv_struct, struct param_sanePic &sanePic_struct, struct saneCheck &check_struct, int rank, int size)
{


	dictionary	*	ini ;

	string bolo_gain_file="";
	string text;
	string filename;
	ofstream file;
	uint16_t parsed = 0x0000;

	parsed+=parser_function(ini_name, output, dir, samples_struct, pos_param, proc_param,
			structPS, saneInv_struct, sanePic_struct, size, rank);


	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return parsed;
	}

	read_saneCheck_ini(ini , check_struct);
	//	read_bolo_gain_global_file(output, ini, dir.input_dir, bolo_gain_file, rank);

	iniparser_freedict(ini);



	if(rank==0){

		// Generates ini file model
		text = "# Sanepic ini file :\n";
		text += "# Edit the lines to modify the options !\n";
		text += "# Some options are required : Check the comments !\n";

		text += "\n\n";

		text += "[commons]\n\n";
		text += "data_directory = " + dir.data_dir + " ; source data directory (only the data themselves)\n";
		text += "input_directory = " + dir.input_dir + " ; input directory with all the configurations files\n";
		text += "output_dir = " + dir.output_dir + " ; output directory\n";
		text += "temp_dir = " + dir.tmp_dir +" ; temporary directory\n";
		text += "fits_filelist = " + dir.fits_filelist + " ; file containing fits file names, [corresponding noise file, [processors indexes]]\n";
		text += "bolo_global_file = " + bolo_gain_file + " ; every scans have the same detector list which name is filled in this field\n";
		text += "bolo_suffix = " + dir.bolo_suffix + " ; bolometers filelist suffix : can be void\n";


		text += "\n\n";

		text += "[sanePos]\n\n";
		text += "pixsize = " + StringOf(pos_param.pixdeg) + " ; size of pixels (degrees)\n";
		text += "map_flagged_data = " + StringOf( pos_param.flgdupl ? "True" : "False" ) + " ; flagged data put in a separated map (default is False)\n";
		text += "file_format = " + StringOf(pos_param.fileFormat) + " ; SANEPIC = 0, HIPE = 1 \n";
		text += "project_gaps = " + StringOf(pos_param.projgaps ? "True" : "False" ) + "; Keyword specifying if gaps are projected to a pixel in the map, if so gap filling of noise only is performed iteratively. Default is False\n";
		text += "mask_file = " + pos_param.maskfile + "; mask (fits file) : use to remove cross constraint due to strong sources\n";
		text += "lon = " + StringOf(pos_param.lon) + "; LON nominal of the input scan data\n";
		text += "lat = " + StringOf(pos_param.lat) + "; LAT nominal of the input scan data\n";
		text += "proj_code = " + StringOf(pos_param.projcode) + "; projection type, see wcslib projection names for the projection you wish to use\n";
		text += "axis_type" + StringOf(pos_param.axistype) + "; axis type, EQ/GAL\n";
		text += "eq2gal" + StringOf(pos_param.eq2gal) + "; Convert J2000.0 EQ to GAL\n";
		text += "gal2eq" + StringOf(pos_param.gal2eq) + "; Convert GAL to J2000.0 \n";
		text += "\n\n";

		text += "[sanePre]\n\n";
		text += "sampling_frequency = " + StringOf(proc_param.fsamp) + " ; detectors sampling frequency (Hz)\n";
		text += "filter_frequency = " + StringOf(proc_param.f_lp) + " ; frequency of the high pass filter applied to the data\n";
		text += "apodize_Nsamples = " + StringOf(proc_param.napod) + " ; number of samples to apodize\n";
		text += "fcut_file = " + proc_param.fcut_file + " ; noise power spectra are thresholded. Default is the frequency cut of the high pass filter applied to the data\n";
		text += "poly_order = " + StringOf(proc_param.poly_order) + " ; baseline polynomia order (default = 0; no baseline)\n";
		text += "linear_baseline = " + StringOf(proc_param.remove_linear ? "True" : "False") + " ; simple linear baseline removed from the data (False : default)\n";
		text += "correlation = " + StringOf(proc_param.CORRon ? "True" : "False") + " ; Set this keyword to False if correlations between detectors are not included in the analysis (default = True)\n";
		text += "nofill_gap = " + StringOf(proc_param.fill_gap ? "True" : "False") +" ; Do we fill the gaps ? F = fill (default), T = don't fill\n";


		text += "\n\n";

		text += "[saneInv]\n\n";
		text += "noise_dir = " + saneInv_struct.noise_dir + " ; cov matrix directory\n";
		text += "cov_matrix_file = " + saneInv_struct.cov_matrix_file + " ; this file contains the matrix you want to invert\n";
		text += "cov_matrix_suffix = " + saneInv_struct.cov_matrix_suffix + " ; this file contains the matrix you want to invert\n";

		text += "\n\n";

		text += "[sanePS]\n\n";
		text += "MixingMatrix_suffix = " + structPS.mix_suffix + " ; Mixing matrix files suffix : the mixmat files are not the same : each scan has a mixmat file named : basename(scan_filename) + MixingMatrix_suffix\n";
		text += "MixingMatrix_global_file = " + structPS.mix_global_file + " ;  the MixingMatrix file  (fill this field if the MixingMatrix file is the same for all the scans !)\n";
		text += "ncomp = " + StringOf(structPS.ncomp) + " ; number of component(s) to estimate\n";
		text += "map_file = " + structPS.signame + " ; fits file containing the map that should be substracted to the data for a second noise estimation step\n";
		text += "ell_global_file = " + structPS.ell_global_file + "; the ell file  (fill this field if the ell file is the same for all the scans !)\n";
		text += "ell_suffix = " + structPS.ell_suffix + " ; ell files suffix : the ell files are not the same : each scan has an ell file named : basename(scan_filename) + ell_suffix\n";
		text += "save_data = " + StringOf(structPS.save_data) + " ; Set save_data to 1 if you wish to save the processing session after each sanePS step \n";

		text += "\n\n";

		text += "[sanePic]\n\n";

		text += "iterW = " + StringOf(sanePic_struct.iterw) + " ; Write temporary map files on disk every iterW number of loop\n";
		text += "iterMAX = " + StringOf(sanePic_struct.itermax) + " ; Maximum number of conjugate gradient loops \n";

		text += "\n\n";


		text += "[saneCheck]\n\n";
		text += "check_NAN = " + StringOf(check_struct.checkNAN ? "True" : "False") + "; check whether there are NANs in every tables and flag them if needed\n";
		text += "check_time_gaps = " + StringOf(check_struct.checktime ? "True" : "False") + "; check whether there are gaps in time table and fill those gaps with flagged data to ensure continuity\n";
		text += "check_flag = " + StringOf(check_struct.checkflag ? "True" : "False") + "; check whether there are detectors that are fully or more than 80% flagged\n";
		text += "check_bolo_gain = " + StringOf(check_struct.checkGain ? "True" : "False") + "; Compute detector gain and print to screen \n";

		text += "\n\n";



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


	return parsed;
}


void read_saneCheck_ini(dictionary	*ini , struct saneCheck &check_struct){


	//	check_struct.bolo_gain_check="";

	check_struct.checkNAN = iniparser_getboolean(ini, "saneCheck:check_NAN", 1);

	check_struct.checktime = iniparser_getboolean(ini, "saneCheck:check_time_gaps", 1);

	check_struct.checkGain = iniparser_getboolean(ini, "saneCheck:check_bolo_gain", 1);

	check_struct.checkflag = iniparser_getboolean(ini, "saneCheck:check_flag", 1);

}

void print_saneCheck_ini(struct saneCheck check_struct){


	cout << endl << "Checks : ..." << endl;

	//	if(check_struct.bolo_gain_check!="")
	//	check_struct.bolo_gain_check=""; // to be printed
	cout <<  "NaNs in data     : ";
	if(!check_struct.checkNAN)
		cout <<  "no" << endl;
	else
		cout <<  "yes" << endl;

	cout << "time gaps        : ";
	if(!check_struct.checktime)
		cout <<  "no" << endl;
	else
		cout <<  "yes" << endl;

	cout << "Gains            : ";
	if(!check_struct.checkGain)
		cout <<  "no" << endl;
	else
		cout <<  "yes" << endl;

	cout << "flags            : ";
	if(!check_struct.checkflag)
		cout <<  "no" << endl;
	else
		cout <<  "yes" << endl;

	cout << endl;

}
