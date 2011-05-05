#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>


#include "mpi_architecture_builder.h"
#include "parser_functions.h"
#include "inputFileIO.h"


extern "C" {
#include <fitsio.h>
}


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
	for (long ii=0;ii<ntotscan+1;ii++)
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

int write_ParallelizationScheme(string fname, long *position, long *frnum, int size, struct samples samples_struct)
// Write the Parallelization Scheme for further use.
{

	ofstream file;
	file.open(fname.c_str(), ios::out);
	if(!file.is_open()){
		cerr << "File [" << fname << "] Invalid." << endl;
		return 1;
	}

	cout << "parallelization scheme written in file : \n" << fname << endl;


	string *fitsvect_temp;
	int *scans_index_temp;
	string temp;
	size_t found;

	fitsvect_temp = new string [samples_struct.ntotscan];
	scans_index_temp = new int [samples_struct.ntotscan];


	for (long ii=0;ii<samples_struct.ntotscan;ii++){
		temp = samples_struct.fitsvect[position[ii]];
		found=temp.find_last_of('/');


		fitsvect_temp[ii] = temp.substr(found+1);
	}

	int val_proc = 0;

	for (long ii=1;ii<samples_struct.ntotscan+1;ii++){
		if(frnum[ii]==0)
			break;

		for(int jj=frnum[ii-1];jj<frnum[ii];jj++)
			scans_index_temp[jj]=val_proc;
		val_proc++;
	}
#ifdef DEBUG
	cout << "nb proc : " << val_proc << endl;
#endif


	if(val_proc>size){
		cerr << "Error in frame order repartition, number of processor are not equal to mpi size\n";
		return -1;
	}

	for (long ii=0;ii<samples_struct.ntotscan;ii++){
		file << fitsvect_temp[ii] << " " << scans_index_temp[ii] << endl;
	}

	file.close();

	delete [] fitsvect_temp;

	return 0;
}

int verify_parallelization_scheme(int rank, struct samples &samples_struct, int size){


	string origin_file;
	if(samples_struct.framegiven)
		origin_file = "fits_filelist";
	else
		origin_file = "parallel_scheme";

	long size_tmp = 0;
	int num_frame = 0;
	char c;
	std::vector<int> index_copy(samples_struct.scans_index);

	struct sortclass_int sortobject;
	sort(index_copy.begin(), index_copy.end(), sortobject);

	std::vector<int>::iterator it;

	// using default comparison:
	it = unique(index_copy.begin(), index_copy.end());
	size_tmp = it - index_copy.begin();

#ifdef DEBUG
	//	if(rank==0){
	cout << " my rank : " << rank << endl;
	cout << "size unique : " << size_tmp << endl;
	cout << size << " vs size : " <<  size_tmp << endl;
	//	}
#endif

	if((size_tmp)>size){
		cerr << "Number of processors are different between MPI and " << origin_file <<  ". Exiting\n";
		return 1;
	}else{

		if((size_tmp)<size){
			index_copy.resize( size_tmp );
			if(rank==0){
				cout << "Warning. The number of processors used in " << origin_file << " is < to the number of processor used by MPI !\n";
				cout << "Do you wish to continue ? (y/n)\n";
				c=getchar();
				switch (c){
				case 'y':
					cout << "Let's continue with only " << (size_tmp) << " processor(s) !\n";
					break;
				default:
					cout << "Exiting ! Please modify " << origin_file << " to use the correct number of processors\n";
					return 1;
					break;
				}
			}
			for(long ii=0;ii<size_tmp;ii++)
				if(index_copy[ii]==0)
					num_frame++;

			if(num_frame==0){
				cout << "Exiting ! Please modify " << origin_file << " to use at least processor 0 \n";
				return 1;
			}


		}else{
			for(int ii=0;ii<size_tmp;ii++)
				if(index_copy[ii]!=ii){
					cerr << "There is a problem in " << origin_file << " : you have forgot a processor to use or the processor numbers are not continuous. Exiting" << endl;
					return 1;
				}
		}
	}

	return 0;

}


int configure_PARA_FRAME_samples_struct(string outdir, struct samples &samples_struct, int rank, int size, long &iframe_min, long &iframe_max){


	struct samples samples_str_para;
	if(!samples_struct.framegiven){

		// get scans order from parallel_scheme
		string para_file = outdir + parallel_scheme_filename;
		string output="";

		// read parallel scheme file
		if(read_fits_list(output , para_file, samples_str_para)){
			cout << output << endl;
			return 1;
		}

		samples_str_para.ntotscan = samples_str_para.fitsvect.size();

		if(!samples_str_para.framegiven){
			cout << "ERROR in " << para_file << ". You must run saneFrameorder because indexes are absent ! Exiting ...\n";
			return 1;
		}

		// check validity between para_scheme and fits_filelist file
		if(check_filelist_validity(samples_struct, samples_str_para))
			return 1;

		// IF validity was ok : copy the index in samples_struct.scans_index
		samples_struct.scans_index.insert(samples_struct.scans_index.begin(),samples_str_para.scans_index.begin(),samples_str_para.scans_index.end());
		//		for(int ii=0;ii<samples_struct.ntotscan;ii++)
		//			samples_struct.scans_index.push_back(samples_str_para.scans_index[ii]);
	}

	// check validity between indexes and mpi #
	if(verify_parallelization_scheme(rank, samples_struct, size))
		return 1;

	// reorder samples_struct
	reorder_samples_struct(rank, samples_struct, size, iframe_min, iframe_max);

	return 0;
}


int check_filelist_validity(struct samples samples_str, struct samples samples_str_para){



#ifdef DEBUG
	cout << "ntotscan" << endl;
	cout << samples_str.ntotscan << " vs " << samples_str_para.ntotscan << endl;
#endif

	if(samples_str.ntotscan!=samples_str_para.ntotscan){
		cerr << "number of scans are different between your fits file_list and the parallel_scheme" << endl;
		return 1;
	}


	struct sortclass_string sortobject;
	sort((samples_str_para.fitsvect).begin(), (samples_str_para.fitsvect).end(), sortobject);
	sort((samples_str.fitsvect).begin(), (samples_str.fitsvect).end(), sortobject);


	for(int ii=0;ii<samples_str.ntotscan;ii++)
		if(samples_str_para.fitsvect[ii]!=(FitsBasename(samples_str.fitsvect[ii])+ ".fits")){

#ifdef DEBUG
			cout << "comparaison triée : " << endl;
			cout << samples_str_para.fitsvect[ii] << " vs " << FitsBasename(samples_str.fitsvect[ii]) + ".fits" << endl;
#endif

			cerr << "Your fits_filelist file and parallel_scheme do not have the same sample files. Exiting\n";
			return 1;
		}

	return 0;
}

void reorder_samples_struct(int rank, struct samples &samples_struct,  int size, long &iframe_min, long &iframe_max){

	long num_frame=0;
	iframe_min=0;
	iframe_max=0;

	// copy the whole vectors
	std::vector<string> fits_copy(samples_struct.fitsvect);
	std::vector<long> nsamples_copy(samples_struct.nsamples);
	std::vector<std::string> bolovect_copy(samples_struct.bolovect);
	std::vector<double> fcut_copy(samples_struct.fcut);


	// reorganize them and define each processor iframe_min and _max !
	for(long ii = 0; ii<size; ii++){
		if(rank==ii)
			iframe_min=num_frame;
		for(long jj = 0; jj<samples_struct.ntotscan; jj++){

			if(samples_struct.scans_index[jj]==ii){
				samples_struct.fitsvect[num_frame]=fits_copy[jj];
				samples_struct.nsamples[num_frame]=nsamples_copy[jj];
				samples_struct.bolovect[num_frame]=bolovect_copy[jj];
				samples_struct.fcut[num_frame]=fcut_copy[jj];

				num_frame++;
			}
		}
		if(rank==ii)
			iframe_max=num_frame;
	}
}

int who_do_it(int size, int rank, int ii)
{

	if(size==1) // if there is only 1 proc, he has to do the job
		return 0;

	if(size>=ii) // if the loop number is smaller than the number of MPI processors
		return ii;

	if(size<ii){ // if the loop number is larger than the number of MPI processors
		while(ii>size)
			ii=ii-size; // find the processor that will do the job by substracting iteratively the number of MPI procs
		return ii;
	}

	return -1;
}
