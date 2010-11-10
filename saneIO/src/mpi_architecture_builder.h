#ifndef MPI_ARCHITECTURE_BUILDER_H_
#define MPI_ARCHITECTURE_BUILDER_H_

#include <vector>
#include <string>
#include "struct_definition.h"

#define parallel_scheme_filename  "parallel_scheme.txt";


struct sortclass_int {
	bool operator() (int i,int j) { return (i<j);}
};

struct sortclass_double {
	bool operator() (double i,double j) { return (i<j);}
};

struct sortclass_long {
	bool operator() (long i,long j) { return (i<j);}
};

struct sortclass_string {
	bool operator() (std::string i,std::string j) { return (i<j);}
};

template<class T> void vector2array(std::vector<T> vect, T* a)
{
	// copy list of type T to array of type T
	typename std::vector<T>::iterator iter;
	int i;

	for (iter=vect.begin(), i=0; iter != vect.end(); iter++, i++) {
		a[i] = *iter;
	}
}

// scans are distributed over processors
void find_best_order_frames(long *position, long *frnum, long *ns, long ntotscan, int size);

int compare_array_double (const void *array_1, const void *array_2);
double randg_archi(long nombre, int seedpass);
double* randg(long nombre, int seedpass);

int write_ParallelizationScheme(std::string fname, long *position, long *frnum, int size, struct samples samples_struct);

int check_ParallelizationScheme(std::string fname, std::string dirfile,struct samples &samples_struct, int size);
int define_parallelization_scheme(int rank,std::string fname, std::string dirfile, string data_dir, struct samples &samples_struct,int size, long &iframe_min, long &iframe_max);
int verify_parallelization_scheme(int rank, std::string outdir,struct samples samples_struct, int size, long &iframe_min, long &iframe_max);

/*! this function determines which processor has to treat the given fits file referenced by his number in the input list */
int who_do_it(int size, int rank, int ii);

#endif /* MPI_ARCHITECTURE_BUILDER_H_ */
