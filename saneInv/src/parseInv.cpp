
#include "parseInv.h"
#include "inputFileIO.h"

#include "parser_functions.h"
#include "struct_definition.h"
#include <iostream>
#include <string>

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

using namespace std;


int parse_saneInv_ini_file(char * ini_name, struct samples &samples_struct,struct common &dir,std::vector<detectors> &detector_tab, int rank)
{
	dictionary	*	ini ;

	/* Some temporary variables to hold query results */
	//	char		*	s ;
	string filename;
	struct detectors det;

	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) { // if dictionnary was not found, exit + error msg
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	//	s = iniparser_getstring(ini, "commons:channel",NULL);
	//	if(s!=NULL){
	//		boloname=s;
	//	}else{
	//		if(rank==0)
	//			printf("You must specify a bolometer file : commons:channel\n");
	//		return(-1);
	//	}


	if(read_common(ini, dir, 0)==-1)
		return -1;

	if(rank==0)
		if(check_path(dir.dirfile, "Data directory") ||
				check_path(dir.input_dir, "Input directory") ||
				check_path(dir.noise_dir, "Noise directory") ||
				check_path(dir.output_dir, "Output directory") ||
				check_path(dir.tmp_dir, "Temporary directory") ||
				check_dirfile_paths(dir.tmp_dir))
			return(-1);


	if(read_fits_file_list(ini, dir,samples_struct, 0)==-1)
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


	if(rank==0){
		cout << "\nYou have specified the following options : \n";
		print_common(dir);
	}

	// cleaning up
	iniparser_freedict(ini);

	return 0;
}


