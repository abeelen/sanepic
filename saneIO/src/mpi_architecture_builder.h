#ifndef MPI_ARCHITECTURE_BUILDER_H_
#define MPI_ARCHITECTURE_BUILDER_H_

#include <vector>
#include <string>
#include "struct_definition.h"

#define parallel_scheme_filename  "parallel_scheme";

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

//! A routine that is used with standard sort function, to get a sorted global vector (dat_compare) indexes
/*!
  \param array_1 A pointer to dat_compare vector value
  \param array_2 An other pointer to dat_compare vector value
  \return An integer which is >0 if dat_compare[*array_1] > dat_compare[*array_2], <0 otherwise
 */
int compare_array_double (const void *array_1, const void *array_2);


//! Generates a random number between o and 1
/*!
 * If seedpass = 0 : time is used as a seed
 * If seedpass = -1 : the seed used rand() is 0
 * If seedpass !=0 and !=-1, the seed = seedpass value
  \param nombre A pointer to dat_compare vector value
  \param seedpass The seed that is used for standard rand() routines is determined by seedpass
  \return A random double between 0/1
 */
double* randg(long nombre, int seedpass);

//! Writes the parallel scheme file to disk
/*!
  \param outdir Output directory : outdir/parallel_scheme
  \param position The processor indice array. This array will be filled by find_best_order_frames
  \param frnum Frame number (proc 0 has to do scans which frame's index goes from frnum(0) to frnum(1), proc 1 from frnum(1) to frnum(2), etc...)
  \param size The number of processors to be used
  \param samples_struct A samples structure that contains everything about frames, noise files and frame processing order
  \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int write_ParallelizationScheme(std::string outdir, long *position, long *frnum, int size, struct samples samples_struct);

//! Verify that samples_struct is correct in regards to the number of processors launched with mpirun
/*!
  \param rank The processor rank given by MPI_Comm_rank, in case paraframe or parabolo is defined
  \param size The number of processors to be used
  \param samples_struct A samples structure that contains everything about frames, noise files and frame processing order
  \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int verify_parallelization_scheme(struct samples &samples_struct, int rank, int size);


//! Test and fill and reorganize samples_struct structure
/*!
 * Reads parallel_scheme.txt first, then calls check_filelist_validity(...), verify_parallelization_scheme(...) and reorder_samples_struct(...)
 * The values stored in iframe_min and iframe_max are different according to which "rank" is calling reorder_samples_struct
  \param outdir output directory
  \param rank The processor rank given by MPI_Comm_rank, in case paraframe or parabolo is defined
  \param size The number of processors to be used
  \param samples_struct A samples structure that contains everything about frames, noise files and frame processing order
  \param iframe_max Ending frame index for processor "rank"
 */
int configure_PARA_FRAME_samples_struct(std::string outdir, struct samples &samples_struct, int rank, int size);

//! Reorder whole samples_struct vectors (fitsvect, bolovect, ...) using processor orders and attributes, to each rank, its frame begin and end indexes
/*!
 * The values stored in iframe_min and iframe_max are different according to which "rank" is calling reorder_samples_struct
  \param rank The processor rank given by MPI_Comm_rank, in case paraframe or parabolo is defined
  \param size The number of processors to be used
  \param samples_struct A samples structure that contains everything about frames, noise files and frame processing order
  \param iframe_max Ending frame index for processor "rank"
 */
void reorder_samples_struct(struct samples & samples_struct, int rank, int size);

#endif /* MPI_ARCHITECTURE_BUILDER_H_ */
