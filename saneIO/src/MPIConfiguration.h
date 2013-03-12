#ifndef MPI_ARCHITECTURE_BUILDER_H_
#define MPI_ARCHITECTURE_BUILDER_H_

#include <utility>
#include <vector>
#include <cstring>
#include <string>
#include <map>

#include "StructDefinition.h"
#include "Utilities.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

#define parallelScheme_filename  "parallel_scheme";
#define processorName_filename   "parallel_names";

using namespace std;

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

//! Verify that samples_struct is correct in regards to the number of processors launched with mpirun
/*!
  \param rank The processor rank given by MPI_Comm_rank, in case paraframe or parabolo is defined
  \param size The number of processors to be used
  \param samples_struct A samples structure that contains everything about frames, noise files and frame processing order
  \return An integer >0 if there were a problem, or 0 if everything went OK
 */
int verifyParallelScheme(struct samples &samples_struct, int rank, int size);

//! Reorder whole samples_struct vectors (fitsvect, bolovect, ...) using processor orders and attributes, to each rank, its frame begin and end indexes
/*!
 * The values stored in iframe_min and iframe_max are different according to which "rank" is calling reorder_samples_struct
  \param rank The processor rank given by MPI_Comm_rank, in case paraframe or parabolo is defined
  \param size The number of processors to be used
  \param samples_struct A samples structure that contains everything about frames, noise files and frame processing order
  \param iframe_max Ending frame index for processor "rank"
 */
void reorderSamplesStruct(struct samples & samples_struct, int rank, int size);

uint16_t printNodeUsage(std::vector<std::string> nodeName, std::vector<long> nodeSizes, std::vector<long> order);


#ifdef USE_MPI

void AssignNodeByProcname(char * proc_names, int size, int *nodeID, std::vector<std::string> & nodeName, std::vector<long> & nodeSizes);

void gatherProcName(int rank, int size, char *proc_names);
uint16_t checkProcName(int rank, int size, std::string outdir);

//! Test and fill and reorganize samples_struct structure
/*! Reads parallel_scheme.txt first, then calls check_filelist_validity(...), verify_parallelization_scheme(...) and reorder_samples_struct(...) */
uint16_t configureMPI(string outdir, struct samples & samples_struct, int rank, int size,
		int  &bolo_rank, int  &bolo_size, int &node_rank, int &node_size,
		MPI_Comm & MPI_COMM_NODE, MPI_Comm & MPI_COMM_MASTER_NODE);

#endif

#endif /* MPI_ARCHITECTURE_BUILDER_H_ */
