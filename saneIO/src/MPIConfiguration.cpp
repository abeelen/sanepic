#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <list>
#include <string>
#include <cstring>

#include <algorithm>
#include <cstring>
#include <map>
#include <numeric>

#include <assert.h>
#include <sys/stat.h>
#include <stdint.h>

#include "MPIConfiguration.h"
#include "ParserFunctions.h"
#include "Utilities.h"
#include "InputFileIO.h"
#include "ErrorCode.h"
#include "Crc.h"

extern "C" {
#include <fitsio.h>
}

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;

double *dat_compare;


double* randg(long nombre, int seedpass) {

	double* nombre_hasard;
	time_t temps;
	temps = time(NULL);

	unsigned int seed = 0;

	if (seedpass == 0) seed = (unsigned int) temps;
	if (seedpass != 0 && seedpass != -1) seed = (unsigned int) seedpass;
	if (seedpass != -1) srandom(seed);

	nombre_hasard= new double[nombre];

	for (long i=0;i<nombre/2;i++) {
		double t1 = (double(rand())/RAND_MAX);
		double t2 = (double(rand())/RAND_MAX);
		nombre_hasard[2*i]=sqrt(-2*log(t1))*cos(2*M_PI*t2);
		nombre_hasard[2*i+1]=sqrt(-2*log(t1))*sin(2*M_PI*t2);
	}

	if (nombre/2!=nombre/2.) {
		double t1 = (double(rand())/RAND_MAX);
		double t2 = (double(rand())/RAND_MAX);
		nombre_hasard[nombre-1]=sqrt(-2*log(t1))*cos(2*M_PI*t2);
	}


	return nombre_hasard;
}

int compare_array_double (const void *array_1, const void *array_2)
{

	const long *casted_array_1 = (const long *) array_1;
	const long *casted_array_2 = (const long *) array_2;

	return (dat_compare[*casted_array_1] > dat_compare[*casted_array_2]) - (dat_compare[*casted_array_1] < dat_compare[*casted_array_2]);
}

void find_best_order_frames(long *position, long *frnum, std::vector<long> ns, long ntotscan, int size){

	long count, a, b;
	double maxproctmp, valmin, stdmin, stdtmp, valtmp;

	long *ns_order;
	double *std, *sizeperproc, *maxproc;
	double *valtemp;

	int nessai = 10000;

	long ntot = 0; // total number of frames (*20 to get sample)

	for (long ii=0;ii<ntotscan;ii++)
		ntot += ns[ii];



	maxproc = new double[nessai];
	std = new double[nessai];
	ns_order = new long[ntotscan];
	sizeperproc = new double[ntotscan];
	dat_compare = new double[ntotscan];

	for(long ii=0;ii<ntotscan;ii++){
		sizeperproc[ii]=0.0;
		dat_compare[ii]=0.0;
		ns_order[ii]=0;
	}

	for(int ii=0;ii<nessai;ii++){
		maxproc[ii]=0.0;
		std[ii]=0.0;
	}

	//init random generator
	valtemp = randg(1,0); // valtemp is a random double between 0/1. with 0 the seed of random function is fixed

	//init arrays
	for (long ii=0;ii<size+1;ii++)
		frnum[ii] = 0;


	for (long jj=0;jj<nessai;jj++){

		for (long kk=0;kk<ntotscan;kk++){
			valtemp = randg(1,-1); // return a random value between 0/1
			dat_compare[kk] = valtemp[0];
			delete [] valtemp;
		}

		for (long ii=0;ii<ntotscan;ii++)
			ns_order[ii] = ns[ii];
		for (long ii=0;ii<ntotscan;ii++)
			position[ii] = ii;

		qsort(position,ntotscan,sizeof(long),compare_array_double);

		for (long ii=0;ii<ntotscan;ii++)
			ns_order[ii] = ns[position[ii]];

		count = 0;
		a = ns_order[0];
		for (long ii=1;ii<=ntotscan;ii++){
			b = 0;
			for(long kk=0;kk<ii;kk++)
				b += ns_order[kk];
			if (abs(b-ntot*double(count+1)/double(size)) <= abs(a-ntot*double(count+1)/double(size))){
				a = b;
			} else {
				count += 1;
				frnum[count] = ii-1;
				sizeperproc[count-1] = 0.0;
				for (long kk=frnum[count-1];kk<ii-1;kk++)
					sizeperproc[count-1] += double(ns_order[kk]);
				a = b;
			}
		}
		frnum[size] = ntotscan;
		sizeperproc[count] = 0.0;
		for (long kk=frnum[count];kk<ntotscan;kk++)
			sizeperproc[count] += double(ns_order[kk]);

		//*********** check values
		maxproctmp = *max_element(sizeperproc, sizeperproc+ntotscan);

		maxproc[jj] = maxproctmp;
		std[jj] = 0.0;
		for(long kk=0;kk<ntotscan;kk++)
			if (sizeperproc[kk] > 0.5)
				std[jj] += (sizeperproc[kk]-double(ntot)/size)*(sizeperproc[kk]-double(ntot)/size)/size;


	}// end of first ntotscan loop

	valmin = *min_element(maxproc, maxproc+nessai);

	stdmin=(double)ntot; // was working on 64 bits but not on 32 ...
	stdmin=stdmin*stdmin;
	for (int ii=0;ii<nessai;ii++)
		if (long(valmin) == long(maxproc[ii]))
			if (std[ii] < stdmin)
				stdmin = std[ii];



	valtmp = 2.0*valmin;
	stdtmp = 2.0*stdmin;
#ifdef DEBUG
	printf("max range min = %lf, std range = %lf\n",valtmp,sqrt(stdtmp));
#endif

	while ((stdtmp > stdmin) || ((long)valtmp > (long)valmin)){


		for (long kk=0;kk<ntotscan;kk++){
			valtemp = randg(1,-1);
			dat_compare[kk] = valtemp[0];
			delete [] valtemp;
		}

		for (long ii=0;ii<ntotscan;ii++)
			ns_order[ii] = ns[ii];
		for (long ii=0;ii<ntotscan;ii++)
			position[ii] = ii;

		qsort(position,ntotscan,sizeof(long),compare_array_double);

		for (long ii=0;ii<ntotscan;ii++)
			ns_order[ii] = ns[position[ii]];

		count = 0;
		a = ns_order[0];
		for (long ii=1;ii<=ntotscan;ii++){
			b = 0;
			for(long kk=0;kk<ii;kk++)
				b += ns_order[kk];
			if (abs(b-ntot*double(count+1)/double(size)) <= abs(a-ntot*double(count+1)/double(size))){
				a = b;
			} else {
				count += 1;
				frnum[count] = ii-1;
				sizeperproc[count-1] = 0.0;
				for (long kk=frnum[count-1];kk<ii-1;kk++)
					sizeperproc[count-1] += double(ns_order[kk]);
				a = b;
			}
		}
		frnum[size] = ntotscan;
		sizeperproc[count] = 0.0;
		for (long kk=frnum[count];kk<ntotscan;kk++)
			sizeperproc[count] += double(ns_order[kk]);


		//*********** check values
		maxproctmp = *max_element(sizeperproc, sizeperproc+ntotscan);
		valtmp = maxproctmp;
		stdtmp = 0.0;
		for(long kk=0;kk<ntotscan;kk++)
			if (sizeperproc[kk] > 0.5)
				stdtmp += (sizeperproc[kk]-double(ntot)/size)*(sizeperproc[kk]-double(ntot)/size)/size;
#ifdef DEBUG
		printf("max range = %lf, std range = %lf\n",valtmp,sqrt(stdtmp));
#endif
	}

	delete [] maxproc;
	delete [] std;
	delete [] ns_order;
	delete [] sizeperproc;
	delete [] dat_compare;

#ifdef DEBUG
	printf("max range = %lf, std range = %lf\n",valtmp,sqrt(stdtmp));
#endif

}

int verifyParallelScheme(struct samples &samples_struct, int rank, int size){


	string origin_file;
	if(samples_struct.framegiven)
		origin_file = "fits_filelist";
	else
		origin_file = "parallel_scheme";

	long size_tmp = 0;
	int num_frame = 0;
	//	char c;

	// Retrieve the MPI size from the scans_indexes...

	std::vector<long> index_copy(samples_struct.scans_index);
	struct sortclass_int sortobject;
	sort(index_copy.begin(), index_copy.end(), sortobject);
	std::vector<long>::iterator it = unique(index_copy.begin(), index_copy.end());
	size_tmp = it - index_copy.begin();

	//	std::list<index> uniq_index(samples_struct.scans_index.begin(),samples_struct.scans_index.end());
	//	uniq_index.sort();
	//	uniq_index.uniq();
	//	size_tmp = uniq_index.size();

#ifdef DEBUG
	//	if(rank==0){
	cout << " my rank : " << rank << endl;
	cout << "size unique : " << size_tmp << endl;
	//	}
#endif

	// Too much rank defined compared to available processors
	if((size_tmp)>size){
		if (rank == 0){
			cerr << "EE - You must use at least " << size_tmp << " processors as defined in " << origin_file <<  endl;
			cerr << "EE - Exiting" << endl;
		}
		return 1;
	}

	// Too few rank defined compared to what available processors
	if((size_tmp)<size){
		if(rank==0){
			cout << "WW - The number of processors defined in " << origin_file << " is lower than the number of processor used by MPI !\n";
			//			cout << "     Do you wish to continue ? (y/[n]) ";
			//			c=getchar();
			//			switch (c){
			//			case 'y':
			cout << "WW - Continuing with only " << (size_tmp) << "/" << size << " processors !" << endl;
			//				break;
			//			default:
			//				cout << "EE - Please modify " << origin_file << " to use the correct number of processors" << endl;
			//				cout << "EE - Exiting ! " << endl;
			//				return 1;
			//				break;
			//			}

			// Check that we are using rank 0
			index_copy.resize( size_tmp );
			for(long ii=0;ii<size_tmp;ii++)
				if(index_copy[ii]==0)
					num_frame++;

			if(num_frame==0){
				cerr << "EE - Please modify " << origin_file << " to use at least processor 0 (master rank)" << endl;
				cerr << "EE - Exiting !" << endl;
				return 1;
			}
		}
	}

	return 0;

}


void reorderSamplesStruct( struct samples & samples_struct, int rank, int size){

	long nFrames = samples_struct.nsamples.size();

	// Reordering all vectors ...
	vector<size_t> index(nFrames);
	// ... first sort scans_index ...
	paired_sort(samples_struct.scans_index, index );

	reorder(samples_struct.scans_index,      index);
	reorder(samples_struct.fitsvect,         index);
	reorder(samples_struct.noisevect,        index);
	reorder(samples_struct.bolovect,         index);
	reorder(samples_struct.basevect,         index);
	reorder(samples_struct.bolo_list,        index);
	reorder(samples_struct.dirfile_pointers, index);

	reorder(samples_struct.fcut,             index);
	reorder(samples_struct.fsamp,            index);
	reorder(samples_struct.fhp,              index);

	if (samples_struct.nsamples.size() != 0)
		reorder(samples_struct.nsamples ,    index);

	if (samples_struct.ell_names.size() != 0)
		reorder(samples_struct.ell_names,    index);
	if (samples_struct.mix_names.size() != 0)

		reorder(samples_struct.mix_names,    index);
	if (samples_struct.nbins.size() != 0)

		reorder(samples_struct.nbins,        index);
	if (samples_struct.ndet.size() != 0)
		reorder(samples_struct.ndet,         index);

}

uint16_t printNodeUsage(vector<string> nodeName, vector<long> nodeSizes, vector<long> order){

	assert(nodeName.size() == nodeSizes.size());

	size_t nframes = order.size();
	size_t nNodes  = nodeName.size();

	// used_nodeSizes will contains the given nodeSize or 0 if one node is not used...
	vector<long> used_nodeSizes;
	used_nodeSizes.assign(nodeSizes.size(), 0);
	for (size_t ii=0; ii < nframes; ii++)
		used_nodeSizes[order[ii]] = nodeSizes[order[ii]];


	vector<string> uniq_nodeName(nodeName);
	vector<string>::iterator it;
	sort(uniq_nodeName.begin(), uniq_nodeName.end());
	it = unique(uniq_nodeName.begin(), uniq_nodeName.end());
	uniq_nodeName.resize(it-uniq_nodeName.begin());

	vector<long> usedProc;
	vector<long> totProc;
	usedProc.assign(uniq_nodeName.size(), 0);
	totProc.assign(uniq_nodeName.size(), 0);

	for (it=uniq_nodeName.begin(); it!= uniq_nodeName.end(); ++it){
		long id = distance(uniq_nodeName.begin(), it);
		for (size_t ii = 0; ii < nNodes ; ii++){
			if (nodeName[ii].compare(*it) == 0){
				usedProc[id] += used_nodeSizes[ii];
				totProc[id]  += nodeSizes[ii];
			}
		}
	}
	cout << "MPI using : " << endl;
	for (size_t ii=0; ii < uniq_nodeName.size(); ii++)
		cout << "\t" << uniq_nodeName[ii] << "\twith " << usedProc[ii] << "/" << totProc[ii] << " process(es)"<< endl;

	if ( accumulate(usedProc.begin(), usedProc.end(), 0) != accumulate(totProc.begin(), totProc.end(), 0) ){
		cerr << endl << "EE - Some processor are not used, please check" << endl;
		return PARA_PROBLEM;
	}

	return EXIT_SUCCESS;
}

#ifdef USE_MPI
void AssignNodeByProcname(char * proc_names, int size, int *nodeID, vector<string> & nodeName, vector<long> & nodeSizes){
	/**
	 *   proc_names : a size*MPI_MAX_PROCESSOR_NAME char array with all proc_names
	 *   nodeID     : a size int array with group index
	 *  nodeSizes   : a long vector of size the number of sub_groups, with item being number of proc in each sub_groups
	 */

	uint32_t  * proc_ids;

	// ... transforming the processor name to checksums ...
	proc_ids = new uint32_t[size];
	for (int ii=0; ii < size ; ii++)
		proc_ids[ii] = checksum(proc_names+(ii*MPI_MAX_PROCESSOR_NAME), strlen(proc_names+(ii*MPI_MAX_PROCESSOR_NAME)), 0);

	// ... finding unique machines name ...
	list<uint32_t> myList (proc_ids,proc_ids+size);
	myList.sort();
	myList.unique();

	nodeSizes.assign(myList.size(), 0);
	nodeName.assign(myList.size(), string(""));

	// ... assign group to rank ...
	// ... find size and name (not optimal...) of each group ...
	for (list<uint32_t>::iterator it=myList.begin(); it!= myList.end(); ++it){
		int id = distance(myList.begin(), it);
		uint32_t checksum  = *it;
		for (int ii=0; ii< size; ii++){
			if (proc_ids[ii] == checksum ){
				nodeID[ii] = id;
				nodeSizes[id] += 1;
				nodeName[id] = proc_names+(ii*MPI_MAX_PROCESSOR_NAME);
			}
		}
	}
}


void gatherProcName(int rank, int size, char *proc_names){
	/** Retrieve all proc name
	 *  \param rank, size : MPI rank & size
	 *  \param proc_names : returned char array must be declared by rank 0 as new char[size*MPI_MAX_PROCESSOR_NAME]
	 */

	char *proc_name;
	int procname_length;

	// Retrieve all the processor_name
	proc_name = new char[MPI_MAX_PROCESSOR_NAME];

	MPI_Get_processor_name(proc_name, &procname_length);
	MPI_Gather(proc_name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, proc_names, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, MPI_COMM_WORLD);

	delete [] proc_name;

}

uint16_t checkProcName(int rank, int size, std::string outdir){
	/** check or store the processor names
	 * \param rank, size : MPI rank & size
	 * \param outdir : output dir for check file
	 */

	char *proc_names = NULL;
	if (rank == 0){
		proc_names = new char[size*MPI_MAX_PROCESSOR_NAME];
	}

	gatherProcName(rank, size, proc_names);

	MPI_Barrier(MPI_COMM_WORLD);

	// Store processor_name or check if file already there
	if (rank == 0){

		string filename=outdir+processorName_filename;
		struct stat buf;

		if (stat(filename.c_str(), &buf) != -1){
			// check consistancy
			string line;
			vector<string> buffer;
			ifstream checkFile(filename.c_str());
			if (checkFile.is_open()) {
				while ( checkFile >> line ) buffer.push_back(line);
				checkFile.close();
			}
			if (size != (int) buffer.size() ) {
				cerr << "EE - Saved file differ in number of processor used (" << buffer.size() << " vs " <<size << " used)" << endl;
				cerr << "EE - Please check " << filename << endl;
				return EXIT_FAILURE;
			}
			bool check = true;
			for (int ii=0; ii< size; ii++){
				check &= (buffer[ii].compare(proc_names+(ii*MPI_MAX_PROCESSOR_NAME)) == 0);
			}
			if (! check){
				cerr << "EE - Saved file differ in processor order or processor name" << endl;
				cerr << "EE - Please check " << filename << endl;
				return EXIT_FAILURE;
			}

		} else {
			ofstream checkFile(filename.c_str());
			for (int ii=0; ii < size; ii++)
				checkFile << proc_names+(ii*MPI_MAX_PROCESSOR_NAME) << endl;
			checkFile.close();
		}

		delete [] proc_names;

	}

	return EXIT_SUCCESS;
}

uint16_t configureMPI(string outdir, struct samples & samples_struct, int rank, int size,
		int &bolo_rank, int &bolo_size, int &node_rank, int &node_size,
		MPI_Comm & MPI_COMM_NODE, MPI_Comm & MPI_COMM_MASTER_NODE){

	if (rank==0)
		cout << endl << "MPI initialization... " << endl;

	uint16_t returnCode = 0;
	char * proc_names = NULL;
	int * nodeID;
	vector<long> nodeSizes;
	vector<string> nodeName;

	uint16_t return_value;

	if(! samples_struct.framegiven){

		if (rank == 0){

			struct samples samples_struct_scheme;

			// get scans order from parallel_scheme
			string scheme_file = outdir + parallelScheme_filename;
			string output="";

			// read parallel scheme file
			if(readFitsList(output , scheme_file, samples_struct_scheme)){
				cerr << endl << output;
				cerr << "EE - Please run saneFrameOrder first" << endl;
				return_value = FILE_PROBLEM;
			} else {

				if(!samples_struct_scheme.framegiven){
					cerr << "EE - node indexes are absent in " << scheme_file << endl;
					cerr << "EE - You must run saneFrameorder first" << endl;
					return_value = FILE_PROBLEM;
				}

				// check validity between para_scheme and fits_filelist file
				if(samples_struct.ntotscan != samples_struct_scheme.ntotscan){
					cerr << endl << "EE - The number of scans is different between the input file list and the saved parallel scheme" << endl;
					cerr << "EE - Please check or rerun saneFrameorder" << endl;
					return_value = FILE_PROBLEM;
				}

				for(long ii=0;ii<samples_struct.ntotscan;ii++){
					if( Basename(samples_struct.fitsvect[ii]) != Basename(samples_struct_scheme.fitsvect[ii]) ){
						cerr << endl << "EE - Your input file list and saved parallel scheme do not have the same sample files" << endl;
						return_value = FILE_PROBLEM;
					}
				}
			}

			// IF validity was ok : copy the index in samples_struct.scans_index
			samples_struct.scans_index = samples_struct_scheme.scans_index;
		}

		// all rank exit if problem
		MPI_Bcast(&return_value, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (return_value)
			return return_value;

		MPI_Barrier(MPI_COMM_WORLD); // other procs wait untill rank 0 has read the fits_list

		returnCode |=  MPI_Bcast_vector_long(samples_struct.scans_index,  0, MPI_COMM_WORLD);

	}


	// reorder samples_struct
	reorderSamplesStruct(samples_struct, rank, size);


	MPI_Barrier(MPI_COMM_WORLD); // other procs wait untill rank 0 has read the fits_list

	if ( checkProcName( rank, size, outdir) )
		return PARA_PROBLEM;

	if (rank == 0)
		proc_names = new char[size*MPI_MAX_PROCESSOR_NAME];

	gatherProcName(rank, size, proc_names);

	MPI_Barrier(MPI_COMM_WORLD);

	nodeID = new int[size];
	if (rank ==0 )
		AssignNodeByProcname(proc_names, size, nodeID, nodeName, nodeSizes);
	MPI_Bcast(nodeID, size, MPI_INT, 0, MPI_COMM_WORLD);

	// Node Groups splitting...
	MPI_Comm_split(MPI_COMM_WORLD, nodeID[rank], rank, &MPI_COMM_NODE);
	MPI_Comm_size(MPI_COMM_NODE,&node_size); // get mpi number of processors in COMM_NODE
	MPI_Comm_rank(MPI_COMM_NODE,&node_rank); // each proc is given a rank within the node
	MPI_Comm_split(MPI_COMM_WORLD, node_rank, 0, &MPI_COMM_MASTER_NODE); // Communicator for rank0 of each node

	switch(samples_struct.parallel_scheme){

	case 0: {
		// ... in the mixed case, each computer become a node with a nodeSizes cpus ...
		bolo_size = node_size;
		bolo_rank = node_rank;

		break;
	}
	case 1: {
		//  .. in the para case, each proc become a node of size 1 (1 proc) ...

		// Reset nodeName for further print out...
		nodeName.clear();
		nodeName.resize(size);
		for (int ii=0; ii< size; ii++){
			nodeID[ii] = ii;

			if(rank == 0)
				nodeName[ii] = proc_names+(ii*MPI_MAX_PROCESSOR_NAME);
		}

		nodeSizes.clear();
		nodeSizes.assign(size, 1);

		bolo_size = 1;
		bolo_rank = 0;

		break;
	}
	}

	if (rank == 0 && (bolo_rank != 0  || node_rank !=0) ){
		cerr << "EE - rank 0 is not sub_rank 0 or node_rank 0" << endl;
		cerr << "     Unpredictable problem may occur" << endl;
		return PARA_PROBLEM;
	}

	// ... and determine iframe_min/iframe_max

	samples_struct.iframe_min = -1;
	samples_struct.iframe_max = -1;

	long iFrame = 0;
	while( samples_struct.scans_index[iFrame] != nodeID[rank] && iFrame < samples_struct.ntotscan ) iFrame++;
	samples_struct.iframe_min = iFrame;
	while( samples_struct.scans_index[iFrame] == nodeID[rank] && iFrame < samples_struct.ntotscan ) iFrame++;
	samples_struct.iframe_max = iFrame;


	if (rank == 0)
		delete [] proc_names;

	delete [] nodeID;

	if (rank == 0)
		if ( printNodeUsage(nodeName, nodeSizes, samples_struct.scans_index) )
			return PARA_PROBLEM;


	return EXIT_SUCCESS;
}

#endif
