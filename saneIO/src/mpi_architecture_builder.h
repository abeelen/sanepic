#ifndef MPI_ARCHITECTURE_BUILDER_H_
#define MPI_ARCHITECTURE_BUILDER_H_

#include <vector>
#include <string>
#include "struct_definition.h"

#define parallel_scheme_filename  "parallel_scheme.txt";

//! A structure that is used to sort integers vectors with standard routine "sort"
struct sortclass_int {
	bool operator() (int i,int j) { return (i<j);}
};

//! A structure that is used to sort integers vectors with standard routine "sort"
struct sortclass_string {
	bool operator() (std::string i,std::string j) { return (i<j);}
};

//! Finds an optimized way to distribute scans over processors
/*!
  \param position The processor indice array. This array will be filled by find_best_order_frames
  \param frnum Frame number (proc 0 has to do scans which frame's index goes from frnum(0) to frnum(1), proc 1 from frnum(1) to frnum(2), etc...)
  \param ns A vector containing the number of samples for each input fits filenames
  \param ntotscan The total number of scans
  \param size The number of processors to be used
 */
void find_best_order_frames(long *position, long *frnum, std::vector<long> ns, long ntotscan, int size);


int compare_array_double (const void *array_1, const void *array_2);


double randg_archi(long nombre, int seedpass);


double* randg(long nombre, int seedpass);


int write_ParallelizationScheme(std::string fname, long *position, long *frnum, int size, struct samples samples_struct);

// Mpi samples_struct reorganisation functions
int verify_parallelization_scheme(int rank, struct samples &samples_struct, int size);


int configure_PARA_FRAME_samples_struct(std::string outdir, struct samples &samples_struct, int rank, int size, long &iframe_min, long &iframe_max);


int check_filelist_validity(struct samples samples_str, struct samples samples_str_para);


void reorder_samples_struct(int rank, struct samples &samples_struct, int size, long &iframe_min, long &iframe_max);


/*! this function determines which processor has to treat the given fits file referenced by his number in the input list */
int who_do_it(int size, int rank, int ii);

#endif /* MPI_ARCHITECTURE_BUILDER_H_ */
