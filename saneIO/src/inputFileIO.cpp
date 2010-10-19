/*
 * inputFileIO.cpp
 *
 *  Created on: 25 sept. 2009
 *      Author: abeelen
 */


#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>

#include "inputFileIO.h"

#define EscapeChar "!#;"

using namespace std;

/*!
 * Reads a detector list in a .txt file
 * Returns a vector of string containing the name of the considered channels
 */
int read_strings(string fname, std::vector<string> &bolos) {
	string line;
	size_t found;

	ifstream inputFile(fname.c_str());
	if (!inputFile.is_open()) {
		cerr << "Error opening file '" << fname << "'. Exiting.\n";
		return 1;
	}

	while (!inputFile.eof()) {
		getline(inputFile, line);
		line.erase(0, line.find_first_not_of(" \t"));	// remove leading white space
		found = line.find_first_of(EscapeChar);          // Check for comment character at the beginning of the line

		if (line.empty() || found == 0 ) continue; 		// skip if empty or commented

		line = line.substr(0, line.find_first_of(" \t")); // pick out first word
		bolos.push_back(line);
	}

	inputFile.close();

	return 0;
}

int read_double(string fname, double *& array, long & size){
	string line;
	vector<double> temp;
	size_t found;

	ifstream inputFile(fname.c_str(), ios::in);
	if (!inputFile.is_open()) {
		cerr << "Error opening file '" << fname << "'. Exiting.\n";
		return 1;
	}

	// Count the number of lines ;
	while(! inputFile.eof()){

		getline(inputFile,line);
		line.erase(0, line.find_first_not_of(" \t"));	// remove leading white space
		found = line.find_first_of(EscapeChar);				// Check for comment character at the beginning of the line
		if (line.empty() || found == 0 ) continue; 		// skip if empty or commented

		temp.push_back(atof(line.c_str()));
	}
	// Last element is an empty line
	temp.pop_back();
	inputFile.close();

	size = temp.size();
	// Memory allocation
	array = new double[size];
	for (long ii=0; ii< size; ii++)
		array[ii] = temp[ii];

	return 0;
}


std::string FitsBasename(std::string path)
{

	size_t found;
	string filename;

	// Strip the path and get the filename
	// Find the last " directory separator
	found = path.find_last_of("/\\");

	if (found != string::npos)
		filename = 	path.substr(found+1);
	else
		filename = path;

	// Strip the file extension (whatever is after the last ".fits"

	found = filename.find_last_of(".fits");
	if (found != string::npos)
			filename = 	filename.substr(0,found-4);

	return filename;
}

