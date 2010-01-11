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

using namespace std;

/*!
 * Reads a detector list in a .txt file
 * Returns a vector of string containing the name of the considered channels
 */
void read_strings(string fname, std::vector<string> &bolos) {
	string line;

	ifstream inputFile(fname.c_str());
	if (!inputFile.is_open()) {
		cerr << "Error opening file '" << fname << "'. Exiting.\n";
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

void read_double(string fname, double *& array, long & size){
	string line;
	vector<double> temp;

	ifstream inputFile(fname.c_str(), ios::in);
	if (!inputFile.is_open()) {
		cerr << "Error opening file '" << fname << "'. Exiting.\n";
		exit(1);
	}

	// Count the number of lines ;
	while(! inputFile.eof()){
		getline(inputFile,line);
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
}


std::string		Basename(std::string path)
//std::string		Basename(char* path)
{
	//	std::string		Result;

	cout << "starts with :" << path << ":" << endl;


	char * pch=NULL;
	char * temp=NULL;
	//	char * tmp=NULL;
	string basename;
	//	strcpy(path.c_str(),tmp);
	//	string temp_str=path;
	//	printf ("Splitting string \"%s\" into tokens:\n",path.c_str());
	pch = strtok ((char *) path.c_str(),"/.");
	while (pch != NULL)
	{
		if (strcmp (pch,(char*)"fits") == 0)
			break;
		temp=pch;
		//		printf ("%s\n",temp);
		pch = strtok (NULL, "/.");

	}

	//	cout << "result : " << temp << endl;
	basename = (string)temp;
	cout << "result :" << basename << ":" << endl;
	//		getchar();
	//	Result = fname;
	//	Result += ext;
	return basename;

}

std::string FitsBasename(std::string path)
{
//	cout << "FitsBasename starts with :" << path << ":" << endl;

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

//	cout << "FitsBasename ends with :" << filename << ":" << endl;
	return filename;
}

