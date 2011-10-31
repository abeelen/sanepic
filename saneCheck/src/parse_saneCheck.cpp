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
		struct samples &samples_struct, struct param_sanePos &pos_param, struct param_saneProc &proc_param,
		struct param_sanePS &PS_param, struct param_saneInv &Inv_param, struct param_sanePic &Pic_param, struct param_saneCheck &Check_param, int rank, int size)
{


	dictionary	*	ini ;

	string bolo_gain_file="";
	string filename;
	ofstream file;
	uint16_t parsed = 0x0000;

	parsed+=parser_function(ini_name, output, dir, samples_struct, pos_param, proc_param,
			PS_param, Inv_param, Pic_param, size, rank);


	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return parsed;
	}

	read_saneCheck_ini(ini , Check_param);
	//	read_bolo_gain(output, ini, dir.input_dir, bolo_gain_file, rank);

	iniparser_freedict(ini);



	if(rank==0){

		string text;
		std::vector<string> key_vect;
		std::vector<string> val_vect;
		std::vector<string> com_vect;


		text = rebuild_ini(dir, proc_param,  pos_param, samples_struct, PS_param, Pic_param, Inv_param);

		//TODO: Shall we include that in rebuild_ini ?
		text += "\n[saneCheck]\n";
		export_param_saneCheck(Check_param,   key_vect, val_vect, com_vect);
		for (unsigned long ii=0; ii < key_vect.size(); ii++)
			text += key_vect[ii] + " = " + val_vect[ii] + " ; " + com_vect[ii] + "\n";
		key_vect.clear(); val_vect.clear(); com_vect.clear();


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


void read_saneCheck_ini(dictionary	*ini , struct param_saneCheck &Check_param){


	//	Check_param.bolo_gain_check="";

	Check_param.checkNAN = iniparser_getboolean(ini, "saneCheck:check_NAN", 1);

	Check_param.checktime = iniparser_getboolean(ini, "saneCheck:check_time_gaps", 1);

	Check_param.checkGain = iniparser_getboolean(ini, "saneCheck:check_bolo_gain", 1);

	Check_param.checkflag = iniparser_getboolean(ini, "saneCheck:check_flag", 1);

}

void print_saneCheck_ini(struct param_saneCheck Check_param){


	cout << endl << "Checks : ..." << endl;

	//	if(Check_param.bolo_gain_check!="")
	//	Check_param.bolo_gain_check=""; // to be printed
	cout <<  "NaNs in data     : ";
	if(!Check_param.checkNAN)
		cout <<  "no" << endl;
	else
		cout <<  "yes" << endl;

	cout << "time gaps        : ";
	if(!Check_param.checktime)
		cout <<  "no" << endl;
	else
		cout <<  "yes" << endl;

	cout << "Gains            : ";
	if(!Check_param.checkGain)
		cout <<  "no" << endl;
	else
		cout <<  "yes" << endl;

	cout << "flags            : ";
	if(!Check_param.checkflag)
		cout <<  "no" << endl;
	else
		cout <<  "yes" << endl;

	cout << endl;

}
