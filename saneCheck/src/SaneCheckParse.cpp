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

#include "DataIO.h"
#include "StructDefinition.h"
#include "MPIConfiguration.h"
#include "ParserFunctions.h"
#include "Utilities.h"
#include "SaneCheckTools.h"
#include "InputFileIO.h"
#include "SaneCheckParse.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

using namespace std;

uint32_t parse_saneCheck_ini_file(char * ini_name, string &output, struct param_common &dir,
		struct samples &samples_struct, struct param_sanePos &pos_param, struct param_saneProc &proc_param,
		struct param_sanePS &PS_param, struct param_saneInv &Inv_param, struct param_sanePic &Pic_param, struct param_saneCheck &Check_param,
		int size, int rank)
{

	dictionary	*	ini ;
	uint32_t parsed = 0x0000;

	parsed+=parser_function(ini_name, output, dir, samples_struct, pos_param, proc_param,
			PS_param, Inv_param, Pic_param, size, rank);

	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return parsed;
	}
	default_param_saneCheck(Check_param);
	read_param_saneCheck(output, ini , Check_param);
	// string bolo_gain_file="";
	//	read_bolo_gain(output, ini, dir.input_dir, bolo_gain_file, rank);

	struct param_saneFix Fix_param;
	default_param_saneFix(Fix_param);
	read_param_saneFix(output, ini, Fix_param);

	iniparser_freedict(ini);



	if(rank==0){

		string text;
		std::vector<string> key_vect;
		std::vector<string> val_vect;
		std::vector<string> com_vect;

		text = rebuild_ini(dir, proc_param,  pos_param, samples_struct, PS_param, Pic_param, Inv_param);
		text += rebuild_ini_saneCheck(Check_param);
		text += rebuild_ini_saneFix(Fix_param);


		string outfile = dir.output_dir + "sanepic_model.ini";
		ofstream file;
		file.open(outfile.c_str(), ios::out);
		if(!file.is_open()){
			cerr << "File [" << outfile << "] Invalid." << endl;
			return -1;
		}

		file << text;
		cout << endl << "Writing model ini file in : " << outfile << endl;

		file.close();
	}


	return parsed;
}

