/*
 * mpi_architecture_builder.h
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */

#ifndef MPI_ARCHITECTURE_BUILDER_H_
#define MPI_ARCHITECTURE_BUILDER_H_

#include <vector>
#include <cstdlib>
#include <string>

#include "stdio.h"

struct corner{
	double x;
	double y;
};

struct box {
	struct corner blc;
	struct corner trc;
};

struct sortclass_int {
	bool operator() (int i,int j) { return (i<j);}
};

struct sortclass_string {
	bool operator() (std::string i,std::string j) { return (i<j);}
};

struct directories {
	std::string dirfile;
	std::string outdir;
	std::string tmp_dir;
};

struct samples {
	std::vector<std::string> fitsvect;
	std::vector<std::string> noisevect;
	std::vector<long> scans_index;
	bool framegiven;
	std::string *fits_table;
	std::string *noise_table;
	long *index_table;
	long *nsamples;
	long ntotscan;
};

struct detectors {
	long ndet;
	std::vector<std::string> boloname;
};

struct input_commons {
	long napod;
	bool NOFILLGAP;
	bool flgdupl;
	double pixdeg;
};

struct user_options {
	double fsamp;
	bool NORMLIN;
	bool remove_polynomia;
	bool CORRon;
	double f_lp;
	//std::string noiseSppreffile;
	bool projgaps;
};


using namespace std;

// scans are distributed over processors
void find_best_order_frames(long *position, long *frnum, long *ns, long ntotscan, int size);
//double* randg(long nombre, int seedpass); // on garde

int compare_array_double (const void *array_1, const void *array_2);
double randg_archi(long nombre, int seedpass);

int write_ParallelizationScheme(string fname, long  *position, long  *frnum, long  *ns,  long  ntotscan, int  size,
		std::vector<string> fitsvect, std::vector<string> noisevect, std::vector<long> &scans_index);
//void read_ParallelizationScheme(string fname,  long **position, long **frnum, long **ns,  long *ntotscan, int *size);
//int check_ParallelizationScheme(string fname, long *ns, long ntotscan, int size, long **position, long **frnum);
//int define_parallelization_scheme(int rank,string fname,long **frnum,long ntotscan,int size,long *nsamples,long *fframes);

///////////////
int check_ParallelizationScheme(string fname, string dirfile,long ntotscan, int size, long *&nsamples, std::vector<string> fitsfiles, std::vector<string> noisefiles,string *&fits_table, string *&noise_table, long *&index_table);
int define_parallelization_scheme(int rank,string fname,string dirfile,long ntotscan,int size, long *&nsamples, std::vector<string> fitsfiles, std::vector<string> noisefiles, string *&fits_table, string *&noise_table, long *&index_table);
///////////////////////////

long readFitsLength(string filename);
void readFrames(std::vector<string> &inputFiles, long *& nsamples);
void read_fits_list(string fname, std::vector<string> &fitsfiles, std::vector<string> &noisefiles, std::vector<long> &frameorder, bool &framegiven);

void readBoxFile(string filename, std::vector<struct box> & boxList);

#define parallel_scheme_filename  "parallel_scheme.txt";

#endif /* MPI_ARCHITECTURE_BUILDER_H_ */
