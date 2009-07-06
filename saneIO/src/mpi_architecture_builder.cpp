#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include "time.h"

#include "mpi_architecture_builder.h"

using namespace std;


template<class T> void vector2array(vector<T> l, T* a)
{
	// copy list of type T to array of type T
	typename vector<T>::iterator iter;
	int i;

	for (iter=l.begin(), i=0; iter != l.end(); iter++, i++) {
		a[i] = *iter;
	}
}


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


int compare_array_double (const void *a, const void *b)
{

	const long *da = (const long *) a;
	const long *db = (const long *) b;

	return (dat_compare[*da] > dat_compare[*db]) - (dat_compare[*da] < dat_compare[*db]);
}


void find_best_order_frames(long *pos, long *frnum, long *ns, long ntotscan, int size){

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

	//init random generator
	valtemp = randg(1,0); // valtemp is a random double between 0/1. with 0 the seed of random function is fixed

	//init arrays
	for (long ii=0;ii<=ntotscan;ii++)
		frnum[ii] = 0;


	for (long jj=0;jj<nessai;jj++){

		for (long kk=0;kk<ntotscan;kk++){
			valtemp = randg(1,-1); // return a random value between 0/1
			dat_compare[kk] = valtemp[0];
		}

		for (long ii=0;ii<ntotscan;ii++)
			ns_order[ii] = ns[ii];
		for (long ii=0;ii<ntotscan;ii++)
			pos[ii] = ii;

		qsort(pos,ntotscan,sizeof(long),compare_array_double);

		for (long ii=0;ii<ntotscan;ii++)
			ns_order[ii] = ns[pos[ii]];

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
		maxproc[jj] = maxproctmp;
		//seeds(jj+1,*) = seed;
		std[jj] = 0.0;
		for(long kk=0;kk<ntotscan;kk++)
			if (sizeperproc[kk] > 0.5)
				std[jj] += (sizeperproc[kk]-double(ntot)/size)*(sizeperproc[kk]-double(ntot)/size)/size;


	}// end of first ntotscan loop

	valmin = *min_element(maxproc, maxproc+nessai);
	//	minmax(maxproc,nessai,&valmin,&temp,&tmpposmin,&tmpposmax);

	stdmin = double(ntot*ntot);
	for (int ii=0;ii<nessai;ii++)
		if (long(valmin) == long(maxproc[ii]))
			if (std[ii] < stdmin)
				stdmin = std[ii];



	valtmp = 2.0*valmin;
	stdtmp = 2.0*stdmin;
	while ((stdtmp > stdmin) || ((long)valtmp > (long)valmin)){


		for (long kk=0;kk<ntotscan;kk++){
			valtemp = randg(1,-1);
			dat_compare[kk] = valtemp[0];
		}

		for (long ii=0;ii<ntotscan;ii++)
			ns_order[ii] = ns[ii];
		for (long ii=0;ii<ntotscan;ii++)
			pos[ii] = ii;

		qsort(pos,ntotscan,sizeof(long),compare_array_double);

		for (long ii=0;ii<ntotscan;ii++)
			ns_order[ii] = ns[pos[ii]];

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

	}

	printf("max range = %lf, std range = %lf\n",valtmp,sqrt(stdtmp));


}

void write_ParallelizationScheme(string fname, long *pos, long *frnum, long *ns,  long ntotscan, int size)
// Write the Parrallelization Scheme for further use.
// working , tested 18/06
{
	FILE *fp;

	if ((fp = fopen(fname.c_str(),"w")) == NULL){
		cerr << "Error : couldn't open file to write parallelization Scheme." << endl << "Exiting..." << endl;
		exit(1);
	}
	// test
	/*int ii;
	fprintf(fp,"%d\n",size);
	fprintf(fp,"%ld\n",ntotscan);
	for (ii=0;ii<ntotscan; ii++)
		fprintf(fp,"%ld ",ns[ii]);
	fprintf(fp,"\n");
	for (ii=0;ii<ntotscan; ii++)
		fprintf(fp,"%ld ",pos[ii]);
	fprintf(fp,"\n");
	for (ii=0;ii<ntotscan+1; ii++)
		fprintf(fp,"%ld ",frnum[ii]);
	fprintf(fp,"\n");*/

	//
	fwrite(&size,     sizeof(int),  1, fp);
	fwrite(&ntotscan, sizeof(long), 1, fp);

	fwrite(ns,        sizeof(long), ntotscan,fp);
	fwrite(pos,       sizeof(long), ntotscan,fp);
	fwrite(frnum,     sizeof(long), ntotscan+1,fp);
	fclose(fp);

}

void read_ParallelizationScheme(string fname,  long **pos, long **frnum, long **ns,  long *ntotscan, int *size)
// read the Parrallelization Scheme for further use.
//TODO: test
{
	FILE *fp;

	if ((fp = fopen(fname.c_str(),"r")) == NULL){
		cerr << "Error : couldn't open file to read parallelization Scheme." << endl << "Exiting..." << endl;
		exit(1);
	}

	fread(size,      sizeof(int),  1, fp);
	fread(ntotscan,  sizeof(long), 1, fp);

	*ns    = new long[*ntotscan];
	*pos   = new long[*ntotscan];
	*frnum = new long[(*ntotscan)+1];

	fread(*ns,       sizeof(long), *ntotscan,fp);
	fread(*pos,      sizeof(long), *ntotscan,fp);
	fwrite(*frnum,   sizeof(long), (*ntotscan)+1,fp);
	fclose(fp);

}

void check_ParallelizationScheme(string fname, long *ns, long ntotscan, int size, long **pos, long **frnum)
// read and check that the saved Parallelization Scheme correspond to the actual data
{

	long /**dummy_frnum, *dummy_pos,*/ *dummy_nsamples, dummy_ntotscan;
	int dummy_size;
	read_ParallelizationScheme(fname, pos, frnum, &dummy_nsamples, &dummy_ntotscan, &dummy_size);



	//compare array function or
	if(ntotscan==(dummy_ntotscan)){
		for(int ii=0;ii<ntotscan;ii++)
			if(ns[ii]!=dummy_nsamples[ii]){
				cerr << "You must use a parallel scheme that fits the actual mpi configuration (number of samples are different)\n" << endl;
				exit(1);
			}
	}else{
		cerr << "You must use a parallel scheme that fits the actual mpi configuration (number of scans are different)\n" << endl;
		exit(1);
	}
	if(size!=(dummy_size)){
		cerr << "You must use a parallel scheme that fits the actual mpi configuration (size is different)\n" << endl;
		exit(1);
	}



}

