

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

#include "boloIO.h"

using namespace std;

/*!
 * Reads a detector list in a .txt file
 * Returns a vector of string containing the name of the considered channels
 */
void read_bolofile(string fname, std::vector<string> &bolos) {
	string line;

	ifstream inputFile(fname.c_str());
	if (!inputFile.is_open()) {
		cerr << "Error opening bolometer file '" << fname << "'. Exiting.\n";
		exit(1);
	}

	while (!inputFile.eof()) {
		getline(inputFile, line);

		line.erase(0, line.find_first_not_of(" \t")); // remove leading white space
		if (line.empty() || line[0] == '#')
			continue; // skip if empty or commented
		line = line.substr(0, line.find_first_of(" \t")); // pick out first word

		bolos.push_back(line);
	}

	inputFile.close();
}

/*!
 * Reads the detectors offsets in a .txt file
 * Returns an array containing the considered channel offsets + the source offsets
 */
void read_bolo_offsets(string field, string file_BoloOffsets, float *scoffsets, double *offsets){

  double lel, xel;
  long temp1, temp2, temp3;
  int nobolo = 1;

  char boloname[100];
  FILE *fp;


  if ((fp = fopen(file_BoloOffsets.c_str(),"r")) == NULL){
    cerr << "ERROR: Can't find offset file. Exiting. \n";
    exit(1);
  }
  while (fscanf(fp, "%s%ld%ld%ld%lf%lf\n", boloname, &temp1, &temp2, &temp3, &lel, &xel) != EOF) {
    if (field == boloname) {
      nobolo = 0;
      if (temp3 == 250){
	offsets[0] = xel/60.0/60.0 - scoffsets[1];
	offsets[1] = lel/60.0/60.0 + scoffsets[0];
      }
      if (temp3 == 350){
	offsets[0] = xel/60.0/60.0 - scoffsets[3];
	offsets[1] = lel/60.0/60.0 + scoffsets[2];
      }
      if (temp3 == 500){
	offsets[0] = xel/60.0/60.0 - scoffsets[5];
	offsets[1] = lel/60.0/60.0 + scoffsets[4];
      }
    }
  }
  fclose (fp);


  if (nobolo){
    cerr << "Bolometer name not found in offset list" << endl;
    exit(1);
  }


}
