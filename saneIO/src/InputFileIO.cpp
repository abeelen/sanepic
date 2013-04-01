#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>

#include <stdint.h>

#include "ImageIO.h"
#include "InputFileIO.h"
#include "CovMatrixIO.h"
#include "ErrorCode.h"

extern "C"{
#include "nrutil.h"
}

#define EscapeChar "!#;"
//#define DelimiterChar " ,;\t"

using namespace std;

std::string replace_all(std::string str, std::string tobe_replace, std::string with_this)
{
	int len = tobe_replace.size(), pos;
	while((pos=str.find(tobe_replace)) != (int)string::npos)
	{
		str.replace(pos, len, with_this);
	}

	return str;
}

std::string Basename(std::string path)
{
	string filename="";
	int i;

	// Strip the path and get the filename
	// Find the last " directory separator
	i=path.rfind("/");

	if((i<0) || (i>(int)path.size()))
		filename = path;
	else
		filename.assign(path.begin()+i+1,path.end());

	return filename;
}

std::string FitsBasename(std::string path)
{
	string filename="";
	int i;

	filename = Basename(path);

	// Strip the file extension (whatever is after the last ".fits"
	i=filename.rfind(".fits");

	if((i>0) && (i<(int)path.size()))
		filename.erase(filename.begin()+i,filename.end());


	return filename;
}

std::string dirfile_Basename(std::string path)
{
	//	size_t found;
	string filename="";

	// remove .fits extension
	filename = FitsBasename(path);

	// From dirfile standard
	// http://getdata.sourceforge.net/dirfile.html
	// change every reserved character by a "_"
	filename=replace_all(filename, "/", "_");
	filename=replace_all(filename, "&", "_");
	filename=replace_all(filename, ";", "_");
	filename=replace_all(filename, "<", "_");
	filename=replace_all(filename, ">", "_");
	filename=replace_all(filename, "|", "_");
	filename=replace_all(filename, ".", "_");


	return filename;
}

uint32_t read_file_line(std::string &output, std::string fname, std::vector<std::string> & content ) {
	/*
	 * Read an ASCII file and return its content in vector of line stripped from comments
	 */

	ifstream file;
	file.open(fname.c_str(), ios::in);
	if(!file.is_open()){
		output += "EE - Could not open file [" + fname + "].\n";
		return FILE_PROBLEM;
	}

	std::string line;
	size_t found;
	content.clear();

	while(!file.eof()){
		getline(file,line);
		line.erase(0, line.find_first_not_of(" \t")); // remove leading white space in the first name
		found = line.find_first_of(EscapeChar); 	  // Check for comment character at the beginning of the filename
		if (found == 0) continue;                     // Skip commented line...
		if (found != std::string::npos) line.erase(found); // remove comments in line if any
		if (line.empty()) continue;                   // Skip empty line...

		content.push_back(line);
	}
	file.close();

	return 0;
}

long readFitsLength(string filename){

	fitsfile *fptr;
	int status = 0;
	long ns;
	char comment[80];

	//	Open the fits file
	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// Go to the signal Extension ...
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ... and check for the NAXIS1 keyword...
	if (fits_read_key(fptr,TLONG, (char *) "NAXIS1", &ns, (char *) &comment, &status)){
		fits_report_error(stderr, status);
		cout << "naxis\n";
		return 1;
	}

	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

	return ns;
}

void readFramesFromFits(struct samples &samples_struct) {

	long nScan  = samples_struct.fitsvect.size();
	samples_struct.nsamples.clear();
	samples_struct.nsamples.assign(nScan,-1); // default value = -1

	// read the size of this particular frame ...

	for (long iframe=0; iframe< samples_struct.ntotscan; iframe++)
		samples_struct.nsamples[iframe] = readFitsLength(samples_struct.fitsvect[iframe]);

}

int readChannelList(std::string &output, std::string fname, std::vector<string> &bolonames){

	bolonames.clear();

	if(read_file(output,fname,bolonames)){
		output += "EE - You must create file specifying bolometer list named " + fname + "\n";
		return FILE_PROBLEM;
	}
	return 0;
}

uint32_t readFitsList(std::string &output, string fname, struct samples &samples_struct ) {

	std::vector<string> &fitsvect = samples_struct.fitsvect;
	std::vector<long> &scans_index = samples_struct.scans_index;
	int &framegiven = samples_struct.framegiven;

	string i_string;
	int    i_int;

	// read file ...
	std::vector<std::string> file_content;
	if ( read_file_line(output, fname, file_content) )
		return FILE_PROBLEM;

	// ... check for number of element in the first line  ...
	// (as reference)
	size_t nb_elem = min(word_count(file_content[0]), (size_t) 2);

	if (nb_elem < 1){
		output += "EE - "+fname+" must contain at least one column\n";
		return FILE_PROBLEM;
	}

	if (nb_elem > 1)
		framegiven = true;
	else
		framegiven= false;

	fitsvect.clear();
	scans_index.clear();

	// fill output vectors

	for (std::vector<string>::iterator it = file_content.begin(); it!=file_content.end(); ++it) {

		// Check for the required number of column  ...
		if ( min(word_count(*it), (size_t) 2) < nb_elem ){
			output += "EE - "+fname+" must contain the same number of columns on each line\n";
			return FILE_PROBLEM;
		}
		// Warn for additionnal (unused column)
		if ( min(word_count(*it), (size_t) 2) > nb_elem ){
			output += "WW - "+fname+" contains additionnal unused column\n";
		}


		// stream the columns...
		istringstream iline(*it);

		switch(nb_elem) {
		case 2:
			iline >> i_string >> i_int;
			fitsvect.push_back(i_string);
			scans_index.push_back(i_int);
			break;
		case 1:
			iline >> i_string;
			fitsvect.push_back(i_string);
			scans_index.push_back(0);
			break;
		default: // for compleness, should never happens
			output += "EE - File [" + fname + "] must have at least one row and at most 2 colums. Exiting\n";
			return FILE_SIZE_PROBLEM;
			break;
		}
	}

	if(fitsvect.size()==0){
		output += "EE - File [" + fname + "] must have at least one row with the correct type : \"string\" \"int\" . Exiting\n";
		return FILE_SIZE_PROBLEM;
	}

	samples_struct.ntotscan = fitsvect.size();

	return 0;
}

