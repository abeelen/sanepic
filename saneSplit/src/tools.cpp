

#include "covMatrixIO.h"
#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parser_functions.h"
#include "tools.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()


extern "C" {
#include "nrutil.h"
#include "nrcode.h"
#include <fitsio.h>
}



using namespace std;



void read_Split_file(string fname, std::vector< long > &cut_sample, struct samples sample_struct){

	std::ifstream file;
	string s, line;
	int nb_elem=0;

	file.open(fname.c_str(), ios::in);
	if(!file.is_open()){
		cerr << "File [" << fname << "] Invalid." << endl;
		exit(-1);
	}

	nb_elem=0;


	cout << " while \n";
	while(file >> s){
		size_t found;
		s.erase(0, s.find_first_not_of(" \t")); // remove leading white space in the first name
		found = s.find_first_of("!#;"); 		// Check for comment character at the beginning of the filename
		if (found == 0) continue;
		cut_sample.push_back(atol(s.c_str()));
		nb_elem++;
		cout << s << endl;
	}


	if(cut_sample[0]<0){
		cout << "Warning : in " << fname << " you must provide Positives cut limits ! Exiting\n";
		exit(EXIT_FAILURE);
	}

	for(int ii=1;ii<nb_elem;ii++)
		if((cut_sample[ii]<0)||((cut_sample[ii]-cut_sample[ii-1])<0)){
			cout << "Warning : in " << fname << " you must provide a crescent order of Positives cut limits !\nExiting\n";
			exit(EXIT_FAILURE);
		}

	file.close();
}

