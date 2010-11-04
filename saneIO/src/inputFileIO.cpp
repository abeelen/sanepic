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

extern "C"{
#include "nrutil.h"
}

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
	std::vector<double> temp;
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

std::string remplace_all(std::string str, std::string tobe_replace, std::string with_this)
{
	int len = tobe_replace.size(), pos;
	while((pos=str.find(tobe_replace)) != (int)string::npos)
	{
		str.replace(pos, len, with_this);
	}

	return str;
}

std::string FitsBasename(std::string path)
{

	//	size_t found;
	string filename;
	int i;


//	cout << path << " deb fonction" << endl;

	path=remplace_all(path, "\\", "/");

	// Strip the path and get the filename
	// Find the last " directory separator
	//	found = path.find_last_of("/\\");
	i=path.rfind("/");

//	cout << " i : " << i << endl;

	if((i<0) || (i>(int)path.size()))
		filename = path;
	else
		filename.assign(path.begin()+i+1,path.end());

//	cout << filename << " apres path" << endl;

	//	if (found != string::npos)
	//		filename = 	path.substr(found+1);
	//	else
	//		filename = path;

	// Strip the file extension (whatever is after the last ".fits"

	//	found = filename.find_last_of(".fits");
	i=filename.rfind(".fits");
	//	cout << " i : " << i << endl;

	//	if (found != string::npos)
	//		filename = 	filename.substr(0,found-4);

	if((i>0) && (i<(int)path.size()))
		filename.erase(filename.begin()+i,filename.end());

//	cout << filename << " fin fonction" << endl;

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
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", NULL, &status)){
		fits_report_error(stderr, status);
		exit(0);
	}

	// ... and check for the NAXIS1 keyword...
	if (fits_read_key(fptr,TLONG, (char *) "NAXIS1", &ns, (char *) &comment, &status)){
		fits_report_error(stderr, status);
		cout << "naxis\n";
		exit(0);
	}

	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

	return ns;
}

void readFrames(std::vector<string> &inputList, long *& nsamples){

	long nScan  = inputList.size();
	nsamples = new long[nScan];
	for (long i=0; i<nScan; i++){
		nsamples[i] = readFitsLength(inputList[i]);
	}

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

		if(read_channel_list(output, filename, det.boloname)==1){
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
				cout << "You are using too many processors : " << size << " processors for only " << detector_tab[oo].ndet << " detectors!";
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


int read_channel_list(std::string &output, std::string fname, std::vector<string> &bolonames){

	if(read_strings(fname,bolonames)){
		output += "You must create file specifying bolometer list named " + fname + "\n";
		return 1;
	}
	return 0;
}


int read_fits_list(string &output, string fname, std::vector<string> &fitsfiles, std::vector<int> &frameorder, bool &framegiven) {



	ifstream file;
	file.open(fname.c_str(), ios::in);
	if(!file.is_open()){
		output += "File [" + fname + "] Invalid.\n";
		return 1;
	}

	framegiven=0;

	string s, line, temp;
	int d;
	char *pch;
	int nb_elem = 0;

	// count number of elements on the first line !
	getline(file, line);
	line.erase(0, line.find_first_not_of(" \t")); // remove leading white space
	pch = strtok ((char*) line.c_str()," ,;\t");

	while (pch != NULL) {
		pch = strtok (NULL, " ,;\t");
		nb_elem++; 	}

	// set pointer back to the beginning of file in order to parse the first line too
	file.seekg (0, ios::beg);

	switch(nb_elem) {
	case 2:
		framegiven=1;
		while(file >> s >> d){
			size_t found;
			s.erase(0, s.find_first_not_of(" \t")); // remove leading white space in the first name
			found = s.find_first_of("!#;"); 		// Check for comment character at the beginning of the filename

			if (found == 0) continue;

			fitsfiles.push_back(s);
			frameorder.push_back(d);
		}
		break;

	case 1:
		while(file >> s){
			size_t found;
			s.erase(0, s.find_first_not_of(" \t")); // remove leading white space in the first name
			found = s.find_first_of("!#;"); 		// Check for comment character at the beginning of the filename
			if (found == 0) continue;

			fitsfiles.push_back(s);
		}
		break;

	default:
		output += "File [" + fname + "] must have at least one row and 2 colums. Exiting\n";
		return 1;
		break;
	}

	if(fitsfiles.size()==0){
		output += "File [" + fname + "] must have at least one row with the correct type : \"string\" \"int\" . Exiting\n";
		return 1;
	}

	if (file>>s){
		output += "File [" + fname + "]. Each line must have the same number of rows. Exiting\n";
		return 1;
	}


#ifdef DEBUG_PRINT
	cout << "read fits list ok !!!\n";
#endif

	file.close();

	return 0;
}


void fill_sanePS_struct(struct param_sanePS &structPS, struct samples samples_struct){


	for(long ii=0;ii<samples_struct.ntotscan;ii++){
		if(structPS.mix_global_file!="")
			structPS.mix_names.push_back(structPS.mix_global_file);
		else
			structPS.mix_names.push_back(FitsBasename(samples_struct.fitsvect[ii]) + structPS.mix_suffix);

		if(structPS.ell_global_file!="")
			structPS.ell_names.push_back(structPS.ell_global_file);
		else
			structPS.ell_names.push_back(FitsBasename(samples_struct.fitsvect[ii]) + structPS.ell_suffix);
	}

}

void fill_noisevect_fcut(struct param_common dir, struct samples &samples_str, struct param_saneInv &saneInv_struct, std::vector<double> &fcut){

	if((saneInv_struct.cov_matrix_file!="")){
		samples_str.noisevect.push_back(saneInv_struct.cov_matrix_file);

	}else{
		for(long iframe = 0; iframe < samples_str.ntotscan ; iframe ++)
			samples_str.noisevect.push_back(FitsBasename(samples_str.fits_table[iframe]) + saneInv_struct.cov_matrix_suffix);
	}


	if((int)fcut.size()==0){

		for(long iframe = 0; iframe < (long)samples_str.noisevect.size() ; iframe ++){

			std::vector<string> bolos;
			long nbins;
			double *ell;
			double **Rellth;

			read_CovMatrix(dir.noise_dir + samples_str.noisevect[iframe], bolos, nbins, ell, Rellth);
			fcut.push_back(ell[0]);
			delete [] ell;
			long nBolos=bolos.size();
			free_dmatrix(Rellth, 0, nBolos * nBolos - 1, 0, nbins - 1);

		}
	}

	if((samples_str.noisevect.size() == 1) && (samples_str.ntotscan > 1)){
		// same noise file for all the scans
		(samples_str.noisevect).resize(samples_str.fitsvect.size(),samples_str.noisevect[0]);

		// if only one fcut, extend to all scans
		fcut.resize(samples_str.ntotscan, fcut[0]);
	}


}
