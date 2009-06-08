#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include "blastSpecific.h"

using namespace std;

foffset* read_mapoffsets(string fname, float *scoffsets, int *nfoff)
{
  char buffer[256];
  float p, y;
  string line, word;
  int f, s0, s1, i, fcount;
  foffset *foffsets;

  ifstream FILE (fname.c_str());
  if (! FILE.is_open()) {
    cerr << "Error opening bolometer offset file '" << fname << "'.\n";
    exit(1);
  }

  // get overall offsets
  while (! FILE.eof()) {
    FILE.getline(buffer, 255);
    line = buffer;

    s0 = 20;
    if (line.substr(0,s0) != "StarCamera2BoreSight") continue;

    // extract 6 fields from line
    i = 0;
    while (i < 6) {
      // find beginning of word
      s0 = line.find_first_not_of(" \t", s0);

      // find end of word
      s1 = line.find_first_of(" \t", s0);

      // get and storeword
      word = line.substr(s0, s1-s0);
      scoffsets[i++] = atof(word.c_str());

      // shift placeholder
      s0 = s1;
    }

    break;
  }

  // find "Begin" tag
  while (! FILE.eof()) {
    FILE.getline(buffer, 255);
    line = buffer;

    if (line.substr(0,6) == "Begin:") break;
  }

  // count frame lines
  fcount = 0;
  while (! FILE.eof()) {
    FILE.getline(buffer, 255);
    line = buffer;

    if (line[0] == '#') continue;
    if (line.substr(0,4) == "End:") break;

    fcount++;
  }

  // allocate memory
  foffsets = new foffset [fcount];
  *nfoff = fcount;

  // reset pointer
  FILE.seekg(0);

  // find "Begin" tag
  while (! FILE.eof()) {
    FILE.getline(buffer, 255);
    line = buffer;

    if (line.substr(0,6) == "Begin:") break;
  }

  // store data
  fcount = 0;
  while (! FILE.eof()) {
    FILE.getline(buffer, 255);
    line = buffer;

    if (line[0] == '#') continue;
    if (line.substr(0,4) == "End:") break;

    sscanf(line.c_str(), "%d%f%f",  &f, &p, &y);

    (foffsets[fcount]).frame = f;
    (foffsets[fcount]).pitch = p;
    (foffsets[fcount]).yaw   = y;

    fcount++;
  }

  return(foffsets);
}




void correctFrameOffsets(int nfoff, long ff, double *offsets, foffset *foffsets, double *froffsets){

	int find ;

	for (find=0; find<nfoff-1; find++) {
		if (ff < (foffsets[find+1]).frame) break;
	}

	froffsets[0] = offsets[0] - (foffsets[find]).yaw;
	froffsets[1] = offsets[1] + (foffsets[find]).pitch;

}




