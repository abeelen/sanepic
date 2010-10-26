#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>

#include "inputFileIO.h"

#define EscapeChar "!#;"

using namespace std;


template <typename T>
std::string StringOf(const T& object){
	std::ostringstream os;
	os << object;
	return os.str();
}


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


int read_bolo_for_all_scans(std::vector<detectors> &detector_tab, struct param_common dir, struct samples samples_struct, int rank, int size){

	string output = "";
	string filename;
	struct detectors det;
	for(long oo=0;oo<samples_struct.ntotscan;oo++){

		if(dir.bolo_global_filename!="")
			filename=dir.input_dir + dir.bolo_global_filename;
		else
			filename=dir.input_dir + FitsBasename(samples_struct.fitsvect[oo]) + dir.suffix ; //  + ".bolo"

		if(read_channel_list(output, filename, det.boloname, rank)==1){
			if(rank==0)
				cout << endl << output << endl;
			return 1;
		}
		det.ndet = (long)((det.boloname).size());
		if (det.ndet == 0) {
			if(rank==0)
				cout << "Must provide at least one channel.\n\n";
			return 1;
		}
		if(size>det.ndet){
			if(rank==0)
				cout << "You are using too many processors : " + StringOf(size) + " processors for only " + StringOf(detector_tab[oo].ndet) + " detectors!";
			return 1;
		}
		detector_tab.push_back(det);
		det.ndet=0;
		det.boloname.clear();
		if(dir.bolo_global_filename!="") {
			detector_tab.resize(samples_struct.ntotscan, detector_tab[0]);
			break;
		}
	}

	return 0;
}


int read_channel_list(std::string &output, std::string fname, std::vector<string> &bolonames, int rank){

	if(read_strings(fname,bolonames)){
		if(rank==0)
			output += "You must create file specifying bolometer list named " + fname + "\n";
		return 1;
	}
	return 0;
}

