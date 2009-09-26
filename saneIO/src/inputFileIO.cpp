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

#include "inputFileIO.h"

using namespace std;

/*!
 * Reads a detector list in a .txt file
 * Returns a vector of string containing the name of the considered channels
 */
void read_strings(string fname, std::vector<string> &bolos) {
	string line;

	ifstream inputFile(fname.c_str());
	if (!inputFile.is_open()) {
		cerr << "Error opening bolometer file '" << fname << "'. Exiting.\n";
		exit(1);
	}

	while (!inputFile.eof()) {
		size_t found;
		getline(inputFile, line);
		line.erase(0, line.find_first_not_of(" \t"));	// remove leading white space
		found = line.find_first_of("!#;");				// Check for comment character at the beginning of the line

		if (line.empty() || found == 0 ) continue; 		// skip if empty or commented

		line = line.substr(0, line.find_first_of(" \t")); // pick out first word
		bolos.push_back(line);
	}

	inputFile.close();
}
