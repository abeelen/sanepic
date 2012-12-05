#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
#include <map>



#include "mpi_architecture_builder.h"
#include "parser_functions.h"
#include "inputFileIO.h"


extern "C" {
#include <fitsio.h>
}

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;

#ifdef USE_MPI
extern MPI_Comm MPI_COMM_SUB;
#endif

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

int write_ParallelizationScheme(string outdir, long *position, long *frnum, int size, struct samples samples_struct)
// Write the Parallelization Scheme for further use.
{

	ofstream file;
	string fname = outdir + parallel_scheme_filename;

	cout << "parallelization scheme written in file : " << fname << endl;


	long *proc_index, *proc_index_sorted;
	proc_index        = new long [samples_struct.ntotscan];
	proc_index_sorted = new long[samples_struct.ntotscan];

	for (long ii=0; ii<size; ii++){
		for (long jj=frnum[ii]; jj< frnum[ii+1]; jj++)
			proc_index_sorted[jj] = ii;
	}

	proc_index = new long[samples_struct.ntotscan];
	for(long hh=0; hh<samples_struct.ntotscan;hh++){
		proc_index[position[hh]] = proc_index_sorted[hh];
	}

	file.open(fname.c_str(), ios::out);
	if(!file.is_open()){
		cerr << "File [" << fname << "] Invalid." << endl;
		return 1;
	}

	for (long ii=0;ii<samples_struct.ntotscan;ii++){
		file << Basename(samples_struct.fitsvect[ii]) << " " << proc_index[ii] << endl;
	}

	file.close();
	delete [] proc_index;
	delete [] proc_index_sorted;

	return 0;
}

int verify_parallelization_scheme(struct samples &samples_struct, int rank, int size){


	string origin_file;
	if(samples_struct.framegiven)
		origin_file = "fits_filelist";
	else
		origin_file = "parallel_scheme";

	long size_tmp = 0;
	int num_frame = 0;
	//	char c;

	// Retrieve the MPI size from the scans_indexes...

	std::vector<int> index_copy(samples_struct.scans_index);
	struct sortclass_int sortobject;
	sort(index_copy.begin(), index_copy.end(), sortobject);
	std::vector<int>::iterator it = unique(index_copy.begin(), index_copy.end());
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

int configure_PARA_FRAME_samples_struct(string outdir, struct samples & samples_struct, int rank, int size){

	struct samples samples_str_para;
	if(!samples_struct.framegiven){
		// get scans order from parallel_scheme
		string para_file = outdir + parallel_scheme_filename;
		string output="";

		// read parallel scheme file
		if(readFitsList(output , para_file, samples_str_para)){
			cout << output << endl;
			return 1;
		}

		// add path to data to insure validity with read_fits_list
		if(!samples_str_para.framegiven){
			cout << "ERROR in " << para_file << ". You must run saneFrameorder first because indexes are absent !\n\n Exiting ...\n";
			return 1;
		}

		// check validity between para_scheme and fits_filelist file
		if(samples_struct.ntotscan!=samples_str_para.ntotscan){
			cerr << "number of scans are different between your fits file_list and the parallel_scheme" << endl;
			return 1;
		}

		for(int ii=0;ii<samples_struct.ntotscan;ii++){
			if( Basename(samples_struct.fitsvect[ii]) != Basename(samples_str_para.fitsvect[ii]) ){
				cerr << "Your fits_filelist file and parallel_scheme do not have the same sample files. Exiting\n";
				return 1;
			}
		}

		// IF validity was ok : copy the index in samples_struct.scans_index
		samples_struct.scans_index = samples_str_para.scans_index;

	}

	// check validity between indexes and mpi #
	if(verify_parallelization_scheme( samples_struct, rank, size))
		return 1;

	// reorder samples_struct
	reorder_samples_struct(samples_struct, rank, size);

	if (samples_struct.iframe_max==samples_struct.iframe_min){ // ifram_min=iframe_max => This processor will not do anything
		cout << "WW - rank " << rank << " not used... (run saneFrameorder to fix)" << endl 	;
	}

	if (samples_struct.iframe_min < 0 || samples_struct.iframe_min > samples_struct.iframe_max || samples_struct.iframe_max	> samples_struct.ntotscan) {
		cerr << "EE - Error distributing frame ranges. Check iframe_min and iframe_max. " << endl;
		cerr << "EE - Exiting" 	<< endl;
		return 1;
	}

	return 0;
}

void reorder_samples_struct( struct samples & samples_struct, int rank, int size){
	// TODO : This routine does not do what it is supposed to do...

	// copy the whole vectors
	std::vector<int>                     scans_index_copy(samples_struct.scans_index);
	std::vector<long>                       nsamples_copy(samples_struct.nsamples);

	std::vector<DIRFILE *>          dirfile_pointers_copy(samples_struct.dirfile_pointers);
	std::vector<string>                     fitsvect_copy(samples_struct.fitsvect);
	std::vector<std::string>               noisevect_copy(samples_struct.noisevect);
	std::vector<std::string>                bolovect_copy(samples_struct.bolovect);
	std::vector<std::string>                basevect_copy(samples_struct.basevect);
	std::vector<std::vector<std::string> > bolo_list_copy(samples_struct.bolo_list);

	std::vector<double>                         fcut_copy(samples_struct.fcut);
	std::vector<double>                        fsamp_copy(samples_struct.fsamp);
	std::vector<double>                          fhp_copy(samples_struct.fhp);


	std::vector<std::string>               ell_names_copy(samples_struct.ell_names);
	std::vector<std::string>              mix_names_copy(samples_struct.mix_names);
	std::vector<long>                         nbins_copy(samples_struct.nbins);
	std::vector<long>                          ndet_copy(samples_struct.ndet);

	long frame_index=0;

	// reorganize them and define each processor iframe_min and _max !
	for(long ii = 0; ii<size; ii++){

		if(rank==ii)
			samples_struct.iframe_min=frame_index;
		for(long jj = 0; jj<samples_struct.ntotscan; jj++){
			if(scans_index_copy[jj]==ii){

				samples_struct.scans_index[frame_index]      = scans_index_copy[jj];

				samples_struct.fitsvect[frame_index]         = fitsvect_copy[jj];
				samples_struct.noisevect[frame_index]        = noisevect_copy[jj];
				samples_struct.bolovect[frame_index]         = bolovect_copy[jj];
				samples_struct.basevect[frame_index]         = basevect_copy[jj];
				samples_struct.bolo_list[frame_index]        = bolo_list_copy[jj];
				samples_struct.dirfile_pointers[frame_index] = dirfile_pointers_copy[jj];

				samples_struct.fcut[frame_index]             = fcut_copy[jj];
				samples_struct.fsamp[frame_index]            = fsamp_copy[jj];
				samples_struct.fhp[frame_index]              = fhp_copy[jj];


				// nsamples is usually filled AFTER reordering..
				if (samples_struct.nsamples.size() != 0)
					samples_struct.nsamples[frame_index]    = nsamples_copy[jj];


				// Only present for sanePS
				if (samples_struct.ell_names.size() != 0)
					samples_struct.ell_names[frame_index]   = ell_names_copy[jj];
				if (samples_struct.mix_names.size() != 0)
					samples_struct.mix_names[frame_index]   = mix_names_copy[jj];
				if (samples_struct.nbins.size() != 0)
					samples_struct.nbins[frame_index]       = nbins_copy[jj];
				if (samples_struct.ndet.size() != 0)
					samples_struct.ndet[frame_index]        = ndet_copy[jj];

				frame_index++;
			}
		}
		if(rank==ii)
			samples_struct.iframe_max=frame_index;
	}
}

