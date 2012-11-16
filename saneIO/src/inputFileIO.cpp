#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>

#include "imageIO.h"
#include "inputFileIO.h"
#include "covMatrix_IO.h"
#include "error_code.h"

extern "C"{
#include "nrutil.h"
}

#define EscapeChar "!#;"
#define DelimiterChar " ,;\t"

using namespace std;

/*!
 * Reads a detector list in a .txt file
 * Returns a vector of string containing the name of the considered channels
 */
int read_strings(string fname, std::vector<string> &bolos) {
	string line;
	size_t found;

	ifstream inputFile(fname.c_str(), ios::in);
	if (!inputFile.is_open()) {
		cerr << "Error opening file '" << fname << "'. Exiting.\n";
		return FILE_PROBLEM;
	}

	while (!inputFile.eof()) {
		getline(inputFile, line);
		line.erase(0, line.find_first_not_of(" \t"));	// remove leading white space
		found = line.find_first_of(EscapeChar);          // Check for comment character at the beginning of the line

		if (line.empty() || found == 0 ) continue; 		// skip if empty or commented

		line = line.substr(0, line.find_first_of(EscapeChar)); // pick out line until Escape Character

		line = line.substr(0, line.find_first_of(" \t")); // pick out first word

		bolos.push_back(line);
	}

	inputFile.close();

	return 0;
}

int read_double(string fname, std::vector<double> &array ){
	string line;
	size_t found;

	ifstream inputFile(fname.c_str(), ios::in);
	if (!inputFile.is_open()) {
		cerr << "Error opening file '" << fname << "'. Exiting.\n";
		return FILE_PROBLEM;
	}

	// Count the number of lines ;
	while(! inputFile.eof()){

		getline(inputFile,line);
		line.erase(0, line.find_first_not_of(" \t"));	// remove leading white space
		found = line.find_first_of(EscapeChar);			// Check for comment character at the beginning of the line

		if (line.empty() || found == 0 ) continue; 		// skip if empty or commented

		line = line.substr(0, line.find_first_of(EscapeChar)); // pick out line until Escape Character

		array.push_back(atof(line.c_str()));
	}

	inputFile.close();

	return 0;
}

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

void readFramesFromFits(struct samples &samples_struct, int rank){

	long * nsamples;
#ifdef PARA_FRAME
	long * nsamples_tot ;
#endif

	long nScan  = samples_struct.fitsvect.size();

	// Create a simple array which will be filled by each rank
	nsamples = new long[nScan];
	fill(nsamples,nsamples+nScan,0);


//TODO This should be made by subrank 0
	// read the size of this particular frame ...
	for (long iframe=0; iframe< samples_struct.ntotscan; iframe++){
	  nsamples[iframe] = readFitsLength(samples_struct.fitsvect[iframe]);
	}

// //TODO : This will not work, as saneIO is compiled with PARA_FRAME or PARA_BOLO for the moment... i.e. USE_MPI in both cases.
// #ifdef PARA_FRAME

// 	// read the size of this particular frame ...
// 	for (long iframe=samples_struct.iframe_min; iframe<samples_struct.iframe_max; iframe++){
// 		nsamples[iframe] = readFitsLength(samples_struct.fitsvect[iframe]);
// 	}

// 	// ... retrieve all the data...
// 	if (rank == 0) {
// 		nsamples_tot = new long[nScan];
// 		fill(nsamples_tot,nsamples_tot+nScan,0);
// 	}
// 	MPI_Reduce(nsamples, nsamples_tot, nScan, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
// 	// .. copy it back
// 	if (rank == 0){
// 		for (long iframe=0; iframe < nScan; iframe++)
// 			nsamples[iframe] = nsamples_tot[iframe];
// 		delete [] nsamples_tot;
// 	}
// 	// .. BCast it...
// 	MPI_Bcast(nsamples, nScan, MPI_LONG, 0, MPI_COMM_WORLD);
// #endif

	// ... and put it into vector format ..
	samples_struct.nsamples.resize(nScan);
	for (long iframe = 0; iframe < nScan; iframe++)
		samples_struct.nsamples[iframe] = nsamples[iframe];

	delete [] nsamples;
}

int read_channel_list(std::string &output, std::string fname, std::vector<string> &bolonames){

	bolonames.clear();

	if(read_strings(fname,bolonames)){
		output += "You must create file specifying bolometer list named " + fname + "\n";
		return 1;
	}
	return 0;
}

void skip_comment(ifstream &file, string &line){

	size_t found;

	// Skip any commented line at the beginning
	while(! file.eof()){
		getline(file,line);
		line.erase(0, line.find_first_not_of(" \t"));	// remove leading white space
		found = line.find_first_of(EscapeChar);			// Check for comment character at the beginning of the line

		if ( (! line.empty()) && (found != 0) ) break; 		// skip if empty or commented
	}
}


uint16_t read_fits_list(std::string &output, string fname, struct samples &samples_struct ) {

	std::vector<string> &fitsvect = samples_struct.fitsvect;
	std::vector<int> &scans_index = samples_struct.scans_index;
	bool &framegiven = samples_struct.framegiven;

	ifstream file;
	file.open(fname.c_str(), ios::in);
	if(!file.is_open()){
		output += "File [" + fname + "] Invalid.\n";
		return FILE_PROBLEM;
	}

	framegiven=0;

	string s, line, temp;
	int d;
	char *pch;
	int nb_elem = 0;

	size_t found;

	// Skip any comments line at the beginning of the file
	skip_comment(file, line);


	// count number of elements on the first valid line !
	line.erase(0, line.find_first_not_of(" \t")); // remove leading white space
	pch = strtok ((char*) line.c_str(), DelimiterChar);

	while (pch != NULL) {
		pch = strtok (NULL, DelimiterChar);
		nb_elem++;
	}

	// set pointer back to the beginning of file in order to parse the first line too and...
	file.seekg (0, ios::beg);

#ifdef DEBUG
	cout << "case :  " << nb_elem << endl;
#endif

	switch(nb_elem) {
	case 2:
		framegiven=1;
		while(!file.eof()){
			getline(file,line);
			line.erase(0, line.find_first_not_of(" \t")); // remove leading white space in the first name
			found = line.find_first_of(EscapeChar); 		// Check for comment character at the beginning of the filename

			if ( line.empty() || found == 0) continue;
#ifdef DEBUG
			cout << "frame_read : " << s << " " << d << endl;
#endif
			istringstream iline(line);
			iline >> s >> d;
			fitsvect.push_back(s);
			scans_index.push_back(d);
		}
		break;

	case 1:
		while(!file.eof()){
			getline(file,line);
			line.erase(0, line.find_first_not_of(" \t")); // remove leading white space in the first name
			found = line.find_first_of(EscapeChar); 		// Check for comment character at the beginning of the filename

			if ( line.empty() || found == 0) continue;

#ifdef DEBUG
			cout << "frame_read : " << s << endl;
#endif
			fitsvect.push_back(line);
			scans_index.push_back(0);
		}
		break;

	default:
		output += "File [" + fname + "] must have at least one row and at most 2 colums. Exiting\n";
		return FILE_SIZE_PROBLEM;
		break;
	}

	if(fitsvect.size()==0){
		output += "File [" + fname + "] must have at least one row with the correct type : \"string\" \"int\" . Exiting\n";
		return FILE_SIZE_PROBLEM;
	}

	if (file>>s){
		output += "File [" + fname + "]. Each line must have the same number of rows. Exiting\n";
		return FILE_SIZE_PROBLEM;
	}

	file.close();

	samples_struct.ntotscan = fitsvect.size();
	return 0;
}

int read_mixmat_txt(string MixMatfile, long ndet, long ncomp, double **&mixmat) {
	FILE *fp;
	long ncomp2;
	double dummy1; // used to read mixing matrix

	if ((fp = fopen(MixMatfile.c_str(), "r")) == NULL) {
		cerr << "ERROR: Can't find Mixing Matrix file. Exiting. \n";
		cout
		<< "Advice : verify the file is in your noise directory and that his name is : "
		<< MixMatfile << endl;
		return 1;
	}
	if (1 != fscanf(fp, "%ld", &ncomp2))
			cerr << "EE - failed reading ncomp2 from file" << MixMatfile << endl;


	mixmat = dmatrix(0, ndet - 1, 0, ncomp - 1);

	for (long ii = 0; ii < ndet; ii++) {
		for (long jj = 0; jj < ncomp2; jj++) {
			if (1 != fscanf(fp, "%lf", &dummy1))
				cerr << "EE - failed reading element from file " << MixMatfile << endl;
			if (jj < ncomp)
				mixmat[ii][jj] = dummy1;
		}
	}
	fclose(fp);
	return 0;

}


