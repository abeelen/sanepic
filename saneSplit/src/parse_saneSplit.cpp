#include <iostream>

#include "parser_functions.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#include "parse_saneSplit.h"

using namespace std;

int parse_saneSplit_ini_file(char * ini_name, struct param_common &dir)
{


	dictionary	*	ini ;


	// load dictionnary
	ini = iniparser_load(ini_name);


	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return 2 ;
	}

	if(read_common(ini, dir, 0)==1) /* read directories infos */
		return 2;

	cout << "You have specified the following options : \n\n";

	print_common(dir); /* print dir locations on console */

	return 0;
}





