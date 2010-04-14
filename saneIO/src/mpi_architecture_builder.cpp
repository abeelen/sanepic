#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>


//#include "time.h"
//#include "dataIO.h"
#include "mpi_architecture_builder.h"
#include "struct_definition.h"
#include "parser_functions.h"

//#ifdef USE_MPI
//#include "mpi.h"
//#endif



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



double randg_archi(long nombre, int seedpass) {

	double nombre_hasard=0.5;
	time_t temps;
	temps = time(NULL);

	unsigned int seed = 0;

	if (seedpass == 0) seed = (unsigned int) temps;
	if (seedpass != 0 && seedpass != -1) seed = (unsigned int) seedpass;
	if (seedpass != -1) srandom(seed);

	//	nombre_hasard= new double[nombre];

	for (long i=0;i<nombre/2;i++) {
		cout << "hmm problem" << endl;
		exit(0);
		//double t1 = (double(rand())/RAND_MAX);
		//double t2 = (double(rand())/RAND_MAX);
		//	nombre_hasard[2*i]=sqrt(-2*log(t1))*cos(2*M_PI*t2);
		///	nombre_hasard[2*i+1]=sqrt(-2*log(t1))*sin(2*M_PI*t2);
	}

	if (nombre/2!=nombre/2.) {
		double t1 = (double(rand())/RAND_MAX);
		double t2 = (double(rand())/RAND_MAX);
		nombre_hasard=sqrt(-2*log(t1))*cos(2*M_PI*t2);
	}//


	return nombre_hasard;
}


int compare_array_double (const void *array_1, const void *array_2)
{

	const long *casted_array_1 = (const long *) array_1;
	const long *casted_array_2 = (const long *) array_2;

	return (dat_compare[*casted_array_1] > dat_compare[*casted_array_2]) - (dat_compare[*casted_array_1] < dat_compare[*casted_array_2]);
}


void find_best_order_frames(long *position, long *frnum, long *ns, long ntotscan, int size){

	long count, a, b;
	double maxproctmp, valmin, stdmin, stdtmp, valtmp;

	long *ns_order;
	double *std, *sizeperproc, *maxproc;
	double *valtemp;

	int nessai = 10000;

	long ntot = 0; // total number of frames (*20 to get sample)

	for (long ii=0;ii<ntotscan;ii++)
		ntot += ns[ii];


	//cout << ntot << endl;

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
	//cout << valtemp << endl;

	//init arrays
	for (long ii=0;ii<ntotscan+1;ii++)
		frnum[ii] = 0;


	for (long jj=0;jj<nessai;jj++){

		for (long kk=0;kk<ntotscan;kk++){
			valtemp = randg(1,-1); // return a random value between 0/1
			//cout << valtemp[0] << endl;
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

		//	for(long kk=0;kk<ntotscan;kk++)
		//  cout << kk << " " <<  sizeperproc[kk] << endl;

		//exit(0);
		//*********** check values
		maxproctmp = *max_element(sizeperproc, sizeperproc+ntotscan);
		//		minmax(sizeperproc,ntotscan,&temp,&maxproctmp,&tmpposmin,&tmpposmax);
		//cout << "maxproc " << jj << " " << maxproctmp << endl;

		maxproc[jj] = maxproctmp;
		//seeds(jj+1,*) = seed;
		std[jj] = 0.0;
		for(long kk=0;kk<ntotscan;kk++)
			if (sizeperproc[kk] > 0.5)
				std[jj] += (sizeperproc[kk]-double(ntot)/size)*(sizeperproc[kk]-double(ntot)/size)/size;


	}// end of first ntotscan loop

	valmin = *min_element(maxproc, maxproc+nessai);
	//	minmax(maxproc,nessai,&valmin,&temp,&tmpposmin,&tmpposmax);

	//stdmin = double(ntot*ntot);
	stdmin=(double)ntot; // was working on 64 bits but not on 32 ...
	stdmin=stdmin*stdmin;
	for (int ii=0;ii<nessai;ii++)
		if (long(valmin) == long(maxproc[ii]))
			if (std[ii] < stdmin)
				stdmin = std[ii];



	valtmp = 2.0*valmin;
	stdtmp = 2.0*stdmin;
#ifdef DEBUG_PRINT
	printf("max range min = %lf, std range = %lf\n",valtmp,sqrt(stdtmp));
#endif
	//	getchar();

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
		//		minmax(sizeperproc,ntotscan,&temp,&maxproctmp,&tmpposmin,&tmpposmax);
		valtmp = maxproctmp;
		stdtmp = 0.0;
		for(long kk=0;kk<ntotscan;kk++)
			if (sizeperproc[kk] > 0.5)
				stdtmp += (sizeperproc[kk]-double(ntot)/size)*(sizeperproc[kk]-double(ntot)/size)/size;
#ifdef DEBUG_PRINT
		printf("max range = %lf, std range = %lf\n",valtmp,sqrt(stdtmp));
#endif
		//		getchar();
	}

	delete [] maxproc;
	delete [] std;
	delete [] ns_order;
	delete [] sizeperproc;
	delete [] dat_compare;

#ifdef DEBUG_PRINT
	printf("max range = %lf, std range = %lf\n",valtmp,sqrt(stdtmp));
#endif

}

//int write_ParallelizationScheme(string fname, long *position, long *frnum, long *ns, long ntotscan, int size,
//		std::vector<string> fitsvect, std::vector<string> noisevect, std::vector<int> &scans_index)
int write_ParallelizationScheme(string fname, long *position, long *frnum, int size, struct samples samples_struct)
// Write the Parrallelization Scheme for further use.
{
	//FILE *fp;

	/*	if ((fp = fopen(fname.c_str(),"w")) == NULL){
	cerr << "Error : couldn't open file to write parallelization Scheme : " << fname << endl << "Exiting..." << endl;
	return -1;
	}*/

	ofstream file;
	file.open(fname.c_str(), ios::out);
	if(!file.is_open()){
		cerr << "File [" << fname << "] Invalid." << endl;
		return -1;
	}

	cout << "parallelization scheme file " << fname << endl;


	string *fitsvect_temp, *noisevect_temp;
	int *scans_index_temp;
	string temp;
	size_t found;

	fitsvect_temp = new string [samples_struct.ntotscan];
	noisevect_temp = new string [samples_struct.ntotscan];
	scans_index_temp = new int [samples_struct.ntotscan];

//	cout << "write" << endl;
//	cout << samples_struct.fitsvect[0] << " "  << samples_struct.fitsvect[1] << endl;
//	cout << samples_struct.noisevect[0] << endl; //" "  << samples_struct.noisevect[1] << endl;

	for (long ii=0;ii<samples_struct.ntotscan;ii++){
		//cout << position[ii] << endl;
		temp = samples_struct.fitsvect[position[ii]];
		found=temp.find_last_of('/');
		// cout << " file: " << str.substr(found+1) << endl;


		fitsvect_temp[ii] = temp.substr(found+1);
		noisevect_temp[ii] = samples_struct.noisevect[position[ii]];
		//cout << fitsvect_temp[ii] << " " << noisevect_temp[ii] << endl;
		// scans_index_temp[ii] = scans_index[position[ii]];
	}

	int val_proc = 0;
	//scans_index[0] = 0;

	for (long ii=1;ii<samples_struct.ntotscan+1;ii++){
		if(frnum[ii]==0)
			break;

		for(int jj=frnum[ii-1];jj<frnum[ii];jj++)
			scans_index_temp[jj]=val_proc;
		val_proc++;
	}
#ifdef DEBUG_PRINT
	cout << "nb proc : " << val_proc << endl;
#endif

	//	//DEBUG
	//	cout << fitsvect[0] << " "  << fitsvect[1] << endl;
	//	cout << noisevect[0] << " "  << noisevect[1] << endl;
	//	cout << scans_index[0] << " "  << scans_index[1] << endl;

	if(val_proc>size){
		cerr << "Error in frame order repartition, number of processor are not equal to mpi size\n";
		return -1;
	}

	for (long ii=0;ii<samples_struct.ntotscan;ii++){
		file << fitsvect_temp[ii] << " " << noisevect_temp[ii] << " " << scans_index_temp[ii] << endl;
		//cout << fitsvect_temp[ii] << " " << noisevect_temp[ii] << " " << endl;
	}

	file.close();

	delete [] fitsvect_temp;
	delete [] noisevect_temp;

//	cout << "fin fonction\n";

	return 0;

}
/*
int read_ParallelizationScheme(string fname,string dirfile, long ntotscan, int size, long *&nsamples,string *&fits_table, long* &index_table)
// read the Parrallelization Scheme for further use.

{
	std::vector<string> fits_dummy;
	std::vector<string> noise_dummy;
	std::vector<long> index_dummy;
	long ntotscan_dummy;
	bool framegiven;
	long *fframes, *nsamples_dummy;

	read_fits_list(fname, fits_dummy, noise_dummy, index_dummy, framegiven);

	if((framegiven==0)||((int)fits_dummy.size()==0))
		return -1;


	for(int ii=0;ii<(int)fits_dummy.size();ii++){
		cout << dirfile + fits_dummy[ii] << endl;
		fitsvect[ii] = dirfile + fits_dummy[ii];}

	readFrames( &ntotscan_dummy , fits_dummy, fframes, nsamples_dummy);


}*/

//int check_ParallelizationScheme(string fname, string dirfile,long ntotscan, int size, long *&nsamples, std::vector<string> fitsfiles, std::vector<string> noisefiles, string *&fits_table, string *&noise_table, long *&index_table)
int check_ParallelizationScheme(string fname, string dirfile,struct samples &samples_struct, int size)
// read and check that the saved Parallelization Scheme corresponds to the actual data
{


	std::vector<string> fits_dummy;
	std::vector<string> noise_dummy;
	std::vector<int> index_dummy;
	long ntotscan_dummy;
	long size_tmp;
	//size_t found;

	bool framegiven;
	long *nsamples_dummy;
	string temp;

	if(read_fits_list(fname, fits_dummy, noise_dummy, index_dummy, framegiven))
		return 1;

#ifdef DEBUG_PRINT
	cout <<" readed list : " << endl;
	for(int ii = 0; ii< (int)fits_dummy.size();ii++)
		cout << fits_dummy[ii] << " " << noise_dummy[ii] << " " << index_dummy[ii] << endl;


	//	cout << fits_dummy[0] << " " << fits_dummy[1] << " " << fits_dummy[2] << " " << fits_dummy[3] << endl;
	//	cout <<  noise_dummy[0] << " " <<  noise_dummy[1] << " " <<  noise_dummy[2] << " " <<  noise_dummy[3] << endl;
	//	cout <<   index_dummy[0] << " " <<  index_dummy[1] << " " <<   index_dummy[2] << " " <<  index_dummy[3] << endl;
	for(int ii = 0; ii< (int)fits_dummy.size();ii++)
		cout << fits_dummy[ii] << " " << noise_dummy[ii] << " " << index_dummy[ii] << endl;
#endif

	//	cout << fits_dummy[0] << " " << fits_dummy[1] << " " << fits_dummy[2] << " " << fits_dummy[3] << endl;
	//	cout <<  noise_dummy[0] << " " <<  noise_dummy[1] << " " <<  noise_dummy[2] << " " <<  noise_dummy[3] << endl;
	//	cout <<   index_dummy[0] << " " <<  index_dummy[1] << " " <<   index_dummy[2] << " " <<  index_dummy[3] << endl;
#ifdef DEBUG_PRINT
	cout << "framegiven : " << framegiven << endl;
#endif

	if((framegiven==0)||((int)fits_dummy.size()==0)){
		cerr << "The file " << fname <<  " is empty\n.Exiting\n";
		return -1;
	}

	ntotscan_dummy=(long)fits_dummy.size();
	if(samples_struct.ntotscan!=ntotscan_dummy){
		cerr << "number of scans are different between your fits file list and the mpi scheme" << endl;
		return -1;
	}

	vector2array(fits_dummy, samples_struct.fits_table);
	vector2array(noise_dummy, samples_struct.noise_table);
	if(index_dummy.size()>0)
		vector2array(index_dummy,  samples_struct.index_table);

	for(int ii=0;ii<(int)fits_dummy.size();ii++){
		//	cout << dirfile + fits_dummy[ii] << endl;
		fits_dummy[ii] = dirfile + fits_dummy[ii];
	}

	readFrames( fits_dummy, nsamples_dummy);

#ifdef DEBUG_PRINT
	cout << "ntotscan" << endl;
	cout << samples_struct.ntotscan << " vs " << ntotscan_dummy << endl;
#endif

	struct sortclass_string sortobject;
	sort(fits_dummy.begin(), fits_dummy.end(), sortobject);
	sort((samples_struct.fitsvect).begin(), (samples_struct.fitsvect).end(), sortobject);

	//	cout << "comparaison triée : " << endl;




	if ((int)((samples_struct.noisevect).size())>0){

		sort (noise_dummy.begin(), noise_dummy.end(), sortobject);
		sort ((samples_struct.noisevect).begin(), (samples_struct.noisevect).end(), sortobject);

		for(int ii=0;ii<samples_struct.ntotscan;ii++)
			if((noise_dummy[ii]!=samples_struct.noisevect[ii])||(fits_dummy[ii]!=samples_struct.fitsvect[ii])){
#ifdef DEBUG_PRINT
				cout << "comparaison triée : " << endl;
				cout << noise_dummy[ii] << endl;
				cout << samples_struct.noisevect[ii] << endl;
				cout << fits_dummy[ii] << endl;
				cout << samples_struct.fitsvect[ii] << endl;
#endif
				cerr << "Parallelscheme file and " << fname << " do not have the same sample files or noise files. Exiting\n";
				return -1;
			}
	}else{
		for(int ii=0;ii<samples_struct.ntotscan;ii++)
			if(fits_dummy[ii]!=samples_struct.fitsvect[ii]){
#ifdef DEBUG_PRINT
				cout << "comparaison triée : " << endl;
				cout << fits_dummy[ii] << endl;
				cout << samples_struct.fitsvect[ii] << endl;
#endif
				cerr << "Parallelscheme file and " << fname << " do not have the same sample files. Exiting\n";
				return -1;
			}
	}

	if(index_dummy.size()>0){
		size_tmp = *max_element(samples_struct.index_table, samples_struct.index_table+samples_struct.ntotscan);

#ifdef DEBUG_PRINT
		cout << size << " vs size : " <<  size_tmp+1 << endl;
#endif

		if((size_tmp+1)!=size){
			cerr << "Number of processors are different between MPI and parallel scheme. Exiting\n";
			return -1;
		}
	}


	for(int ii=0;ii<samples_struct.ntotscan;ii++)
		samples_struct.nsamples[ii]=nsamples_dummy[ii];

	delete [] nsamples_dummy;
	return 0;
}



//int define_parallelization_scheme(int rank,string fname,string dirfile,long ntotscan,int size, long *&nsamples, std::vector<string> fitsfiles, std::vector<string> noisefiles, string *&fits_table, string *&noise_table, long *&index_table, long iframe_min, long iframe_max){
int define_parallelization_scheme(int rank,string fname,string dirfile,struct samples &samples_struct,int size, long &iframe_min, long &iframe_max){

	//cout << "rank" << rank << endl;
	int test=0;

	test=check_ParallelizationScheme(fname,dirfile,samples_struct,size);
	if (test==-1)
		return test;


#ifdef DEBUG_PRINT
	cout << "Et ca donne ca !" << endl;

	cout << samples_struct.fits_table[0] << " " << samples_struct.fits_table[1] << " " << samples_struct.fits_table[2] << " " << samples_struct.fits_table[3] << endl;
	cout << samples_struct.noise_table[0] << " " << samples_struct.noise_table[1] << " " << samples_struct.noise_table[2] << " " << samples_struct.noise_table[3] << endl;
	cout << samples_struct.index_table[0] << " " << samples_struct.index_table[1] << " " << samples_struct.index_table[2] << " " << samples_struct.index_table[3] << endl;
	cout << samples_struct.nsamples[0] << " " << samples_struct.nsamples[1] << " " << samples_struct.nsamples[2] << " " << samples_struct.nsamples[3] << endl;


	cout << "mon rank : " << rank << endl;

#endif

	iframe_min = -1;

	for(long ii=0;ii<samples_struct.ntotscan;ii++){
		if((samples_struct.index_table[ii]==rank)&&(iframe_min == -1)){
			iframe_min=ii;
			break;
		}
	}

	iframe_max=iframe_min;
	for(iframe_max=iframe_min;iframe_max<samples_struct.ntotscan-1;iframe_max++)
		if(samples_struct.index_table[iframe_max]!=rank){
			iframe_max--;
			break;
		}

	iframe_max++;

#ifdef DEBUG_PRINT
	cout << rank << " iframe_min : " << iframe_min << endl;
	cout << rank << " iframe_max : " << iframe_max << endl;
#endif

	for(long ii=0;ii<samples_struct.ntotscan;ii++)
		samples_struct.fits_table[ii] = dirfile + samples_struct.fits_table[ii];


	return 0;

}


int verify_parallelization_scheme(int rank, string outdir,struct samples samples_struct, int size, long &iframe_min, long &iframe_max){


	ofstream file;


	long size_tmp = 0;
	int return_error = 0;
	int num_frame = 0;
	char c;
	vector2array(samples_struct.scans_index,  samples_struct.index_table);
	//	if(rank==0){

	struct sortclass_int sortobject;
	sort(samples_struct.scans_index.begin(), samples_struct.scans_index.end(), sortobject);

	std::vector<int>::iterator it;
	//			int size_tmp=0;

	// using default comparison:
	it = unique(samples_struct.scans_index.begin(), samples_struct.scans_index.end());
	size_tmp = it - samples_struct.scans_index.begin();

#ifdef DEBUG_PRINT
	cout << "size unique : " << size_tmp << endl;

	cout << size << " vs size : " <<  size_tmp << endl;
#endif

	if((size_tmp)>size){
		cerr << "Number of processors are different between MPI and parallel scheme. Exiting\n";
		return return_error =1;
	}else{

		samples_struct.scans_index.resize( size_tmp );

		cout << "trié + unique : " << samples_struct.scans_index[0] <<  " " << samples_struct.scans_index[1] << endl;


		if((size_tmp)<size){
			if(rank==0){
				cout << "Warning. The number of processors used in fits_filelist is < to the number of processor used by MPI !\n";
				cout << "Do you wish to continue ? (y/n)\n";
				c=getchar();
				switch (c){
				case('y') :
					cout << "Let's continue with only " << (size_tmp) << " processor(s) !\n";
				break;
				default:
					cout << "Exiting ! Please modify fits filelist to use the correct number of processors\n";
					return return_error =1;
					break;
				}
			}
			for(long ii=0;ii<size_tmp;ii++)
				if(samples_struct.scans_index[ii]==0)
					num_frame++;

			if(num_frame==0){
				cout << "Exiting ! Please modify fits filelist to use at least processor 0 \n";
				return return_error =1;
			}


		}else{


			for(long ii=0;ii<size_tmp;ii++)
				if(samples_struct.scans_index[ii]!=ii){
					cerr << "There is a problem in the fits filelist : you have forgot a processor to use. Exiting" << endl;
					return return_error =1;
				}
		}
	}
	//	}



	if(rank==0){

		//		string outfile = outdir + samples_struct.filename + ".txt";
		string outfile = outdir + parallel_scheme_filename;
		//		cout << "outfile : " << outfile;

		file.open(outfile.c_str(), ios::out);
		if(!file.is_open()){
			cerr << "File [" << outfile << "] Invalid : Unable to create it !" << endl;
			return return_error = 1;
		}
	}


	//	MPI_Barrier(MPI_COMM_WORLD);
	//	MPI_Bcast(&return_error,1,MPI_INT,0,MPI_COMM_WORLD);

	//	if(return_error>0){
	//		MPI_Finalize();
	//		exit(0);
	//
	//	}

	string temp;
	size_t found;

	num_frame=0;
	iframe_min=0;
	iframe_max=0;

	long * nsamples_temp;
	nsamples_temp = new long[samples_struct.ntotscan];

	for(long jj = 0; jj<samples_struct.ntotscan; jj++)
		nsamples_temp[jj]= samples_struct.nsamples[jj];


	for(long ii = 0; ii<size; ii++){
		if(rank==ii)
			iframe_min=num_frame;
		for(long jj = 0; jj<samples_struct.ntotscan; jj++){
			if(samples_struct.index_table[jj]==ii){

				samples_struct.fits_table[num_frame]=samples_struct.fitsvect[jj];
				samples_struct.noise_table[num_frame]=samples_struct.noisevect[jj];
				samples_struct.nsamples[num_frame]=nsamples_temp[jj];
				if(rank==0){
					temp = samples_struct.fits_table[num_frame];
					found=temp.find_last_of('/');
					file << temp.substr(found+1) << " " << samples_struct.noise_table[num_frame] << " " << ii << endl;

				}
				num_frame++;
			}
		}
		if(rank==ii)
			iframe_max=num_frame;
	}

	delete [] nsamples_temp;

	if(rank==0)
		file.close();

	return 0;

}



long readFitsLength(string filename){

	fitsfile *fptr;
	int status = 0;
	long ns;
	char comment[80];

	//	Open the fits file
	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// Go to the signal Extension ...
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", NULL, &status)){
		fits_report_error(stderr, status);
		cout << "signal\n";
		exit(0);
	}

	// ... and check for the NAXIS1 keyword...
	if (fits_read_key(fptr,TLONG, (char *) "NAXIS1", &ns, (char *) &comment, &status)){
		fits_report_error(stderr, status);
		cout << "naxis\n";
		exit(0);
	}

	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

	return ns;
}

void readFrames(std::vector<string> &inputList, long *& nsamples){

	//read_strings(filename, inputList);

	long nScan  = inputList.size();
	nsamples = new long[nScan];
	for (long i=0; i<nScan; i++){
		nsamples[i] = readFitsLength(inputList[i]);
	}

}





void readBoxFile(string filename, std::vector<struct box> & boxList){
	// Read a file with 4 number on a line, describing the boxes
	// 2 numbers for the bottom_left_corner (blc)
	// 2 numbers for the top_right_corner (trc)

	ifstream file;
	file.open(filename.c_str(), ios::in);
	if(!file.is_open()){
		cerr << "File [" << filename << "] Invalid." << endl;
		exit(-1);
	}

	double x_min, x_max, y_min, y_max;

	while(file >> x_min >> y_min >> x_max >> y_max) {
		struct box ibox;
		struct corner icorn;

		icorn.x = x_min;
		icorn.y = y_min;
		ibox.blc = icorn;

		icorn.x = x_max;
		icorn.y = y_max;
		ibox.trc = icorn;

		boxList.push_back(ibox);
	}

	file.close();

}
