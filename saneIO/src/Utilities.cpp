#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <unistd.h>

#include <stdint.h>

#include "Utilities.h"
#include "StructDefinition.h"

#ifdef USE_MPI
#include "mpi.h"
#endif


using namespace std;


char* tableFormat(std::vector<string> strings) {
	// Return the table Format for the given vector of string
	stringstream stream_tform;
	string string_tform;
	char* p_tform;

	stream_tform << maxStringLength(strings) << "A";
	string_tform = stream_tform.str();

	p_tform = new char[string_tform.size() + 1];
	strcpy(p_tform, string_tform.c_str());
	return p_tform;
}

long maxStringLength(std::vector<string> strings) {
	// Return the longest string length of a string vector
	unsigned long maxSize = 0;
	std::vector<string>::iterator itString;

	for (itString = strings.begin(); itString != strings.end(); itString++) {
		string iString = *(itString);
		if (iString.size() > maxSize)
			maxSize = iString.size();

	}
	return maxSize;
}

char** vString2carray(std::vector<string> strings) {
	// Transform a vector of string into a array of char

	int stringLength = maxStringLength(strings);
	int nBolos = strings.size();

	char **data;

	data = new char*[nBolos];

	for (int i = 0; i < nBolos; i++) {
		data[i] = new char[stringLength+1]; // TODO : +\0 ?
		strcpy(data[i], strings[i].c_str());
	}

	return data;

}

void vDouble2carray(std::vector<double> doubles, double ** data, long *size){
	// Transform a vector of double into a array of double, safe way

	*size = doubles.size();
	*data = new double[(*size)];
	for (long ii=0; ii<(*size); ii++){
		(*data)[ii] = doubles[ii];
	}

}

unsigned long max_size_str_vect(std::vector<std::string> str_vect)
{
  unsigned long size_max=0;

  for(unsigned long ii=0; ii < str_vect.size(); ii++)
    if( str_vect[ii].size()>size_max)
      size_max=str_vect[ii].size();

  return size_max;
}

long getTotalSystemMemory()
{
	long pages = sysconf(_SC_PHYS_PAGES);
	long page_size = sysconf(_SC_PAGE_SIZE);
	return pages * page_size;
}

long getAvailableSystemMemory()
{
	// This does not work as expected as Linux will use most of the RAM for IO caching, thus leaving very little available memory
	long pages = sysconf(_SC_AVPHYS_PAGES);
	long page_size = sysconf(_SC_PAGE_SIZE);
	return pages * page_size;
}


void computeMemoryRequirement_sanePic(struct samples samples_struct, long npix, long indpix_size){
  for (long iframe  = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++){
    samples_struct.memory[iframe] = indpix_size + npix * 7 + samples_struct.nsamples[iframe] * 9;
    // poly_order = 5 to be safe
    samples_struct.memory[iframe] += max( npix + (samples_struct.ndet[iframe]+2)*samples_struct.nbins[iframe], samples_struct.nsamples[iframe] * (1+5) );
  }
}



std::string prettyPrintSize(double size){
  /*
   * Return a pettry string from a size in byte
   */
  ostringstream returnString;
  int power;
  const char * SI_prefix[] = { "", "k", "M", "G", "T"};

  power = min((int) floor(log(size)/log(1024.)),4);
  returnString << fixed << setprecision(2) << size/(pow(1024.,power)) << " ";
  returnString << SI_prefix[power] << "B";

  return returnString.str();

}


std::string prettyPrintAngle(float angle){
  /*
   * Return a pettry string from an angle in degree
   */
  ostringstream returnString;
  int power;

  const char * angle_suffix[] = { "deg", "arcmin", "arcsec"};

  power = max((int) floor(log(angle)/log(60)),-2);
  returnString << fixed << setprecision(2) << angle/(pow((double) 60.,power)) << " ";
  returnString << angle_suffix[(int) abs(power)] ;

  return returnString.str();

}




#ifdef USE_MPI

uint32_t MPI_Bcast_vector_string(std::vector<std::string> & strvect, int root, MPI_Comm Comm){

	int irank;
	uint32_t returnCode = 0;

	MPI_Comm_rank(Comm,&irank);


	long strmaxsize;
	long vectlength;

	if (irank == root){
		vectlength = strvect.size();
		strmaxsize = max_size_str_vect(strvect);
	}

	MPI_Barrier(Comm);
	returnCode |= MPI_Bcast(&strmaxsize, 1, MPI_LONG, root, Comm);
	returnCode |= MPI_Bcast(&vectlength, 1, MPI_LONG, root, Comm);

	char* temp = new char[vectlength*(strmaxsize+1)];

	if(irank==root){
	  fill(temp,temp+vectlength*(strmaxsize+1),'\0');
	  for(long ii=0; ii< vectlength; ii++) {
	    strcpy (temp+ii*(strmaxsize+1), strvect[ii].c_str());
	  }
	}

	MPI_Barrier(Comm);
	returnCode |= MPI_Bcast(temp, vectlength*(strmaxsize+1), MPI_CHAR, root, Comm);

	if (irank != root){
	  strvect.assign(vectlength,"");
	  for(long ii=0; ii< vectlength; ii++) {
	    strvect[ii] = temp+ii*(strmaxsize+1);
	  }
	}
	delete [] temp;
	return returnCode;
}

uint32_t MPI_Bcast_vector_int(std::vector<int> & intvect, int root, MPI_Comm Comm){

	int irank;
	uint32_t returnCode = 0;

	MPI_Comm_rank(Comm,&irank);

	long vectlength;
	int* dummy;

	if (irank == root){
		vectlength = intvect.size();
	}

	MPI_Barrier(Comm);
	returnCode |= MPI_Bcast(&vectlength, 1, MPI_LONG, root, Comm);

	if (irank == root)
		dummy = &(intvect)[0];
	else
		dummy = new int[vectlength];

	MPI_Barrier(Comm);
	returnCode |= MPI_Bcast(dummy, vectlength, MPI_INT, root, Comm);

	if (irank != root) {
		intvect.clear();
		for(long ii=0; ii<vectlength; ii++)
			intvect.push_back(dummy[ii]);
		delete [] dummy;
	}
	return returnCode;
}

uint32_t MPI_Bcast_vector_long(std::vector<long> & longvect, int root, MPI_Comm Comm){

	int irank;
	uint32_t returnCode = 0;

	MPI_Comm_rank(Comm,&irank);

	long vectlength;
	long * dummy;

	if (irank == root){
		vectlength = longvect.size();
	}

	MPI_Barrier(Comm);
	returnCode |= MPI_Bcast(&vectlength, 1, MPI_LONG, root, Comm);

	if (irank == root)
		dummy = &(longvect)[0];
	else
		dummy = new long[vectlength];

	MPI_Barrier(Comm);
	returnCode |= MPI_Bcast(dummy, vectlength, MPI_LONG, root, Comm);

	if (irank != root) {
		longvect.clear();
		for(long ii=0; ii<vectlength; ii++)
			longvect.push_back(dummy[ii]);
		delete [] dummy;
	}
	return returnCode;
}

uint32_t MPI_Bcast_dmatrix(double ** buffer, long nrow, long ncol, int root, MPI_Comm Comm){

  uint32_t returnCode = 0;

  for (long iRow=0; iRow < nrow; iRow++)
    returnCode |= MPI_Bcast(buffer[iRow], (int) ncol, MPI_DOUBLE, root, Comm);

  return returnCode;

}


uint32_t MPI_Reduce_dmatrix(double **sendrecvbuf, long nrow, long ncol, MPI_Op op, int root, MPI_Comm Comm){

  uint32_t returnCode = 0;

  int irank;
  MPI_Comm_rank(Comm,&irank);

  for (long iRow=0; iRow < nrow; iRow++){
    if (irank == root)
      returnCode |= MPI_Reduce(MPI_IN_PLACE, sendrecvbuf[iRow], (int) ncol, MPI_DOUBLE, op, root, Comm);
    else
      returnCode |= MPI_Reduce(sendrecvbuf[iRow], NULL, (int) ncol, MPI_DOUBLE, op, root, Comm);
  }

  return returnCode;

}

uint32_t MPI_Allreduce_dmatrix(double **sendrecvbuf, long nrow, long ncol, MPI_Op op, MPI_Comm Comm){

  uint32_t returnCode = 0;

  for (long iRow=0; iRow < nrow; iRow++)
      returnCode |= MPI_Allreduce(MPI_IN_PLACE, sendrecvbuf[iRow], (int) ncol, MPI_DOUBLE, op, Comm);

  return returnCode;

}
#endif
