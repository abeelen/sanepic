#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include "time.h"

#include "mpi_architecture_builder.h"


using namespace std;


template<class T> void vector2array(std::vector<T> vect, T* a)
{
	// copy list of type T to array of type T
	typename std::vector<T>::iterator iter;
	int i;

	for (iter=vect.begin(), i=0; iter != vect.end(); iter++, i++) {
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

	maxproc = new double[nessai];
	std = new double[nessai];
	ns_order = new long[ntotscan];
	sizeperproc = new double[ntotscan];
	dat_compare = new double[ntotscan];

	//init random generator
	valtemp = randg(1,0); // valtemp is a random double between 0/1. with 0 the seed of random function is fixed
	cout << valtemp << endl;

	//init arrays
	for (long ii=0;ii<=ntotscan;ii++)
		frnum[ii] = 0;


	for (long jj=0;jj<nessai;jj++){

		for (long kk=0;kk<ntotscan;kk++){
			valtemp = randg(1,-1); // return a random value between 0/1
			cout << valtemp << endl;
			dat_compare[kk] = valtemp[0];
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

	}

	printf("max range = %lf, std range = %lf\n",valtmp,sqrt(stdtmp));


}

int write_ParallelizationScheme(string fname, long *position, long *frnum, long *ns,  long ntotscan, int size)
// Write the Parrallelization Scheme for further use.
// working , tested 18/06
{
	FILE *fp;

	if ((fp = fopen(fname.c_str(),"w")) == NULL){
		cerr << "Error : couldn't open file to write parallelization Scheme." << endl << "Exiting..." << endl;
		return -1;
	}
	// test
	/*int ii;
	fprintf(fp,"%d\n",size);
	fprintf(fp,"%ld\n",ntotscan);
	for (ii=0;ii<ntotscan; ii++)
		fprintf(fp,"%ld ",ns[ii]);
	fprintf(fp,"\n");
	for (ii=0;ii<ntotscan; ii++)
		fprintf(fp,"%ld ",position[ii]);
	fprintf(fp,"\n");
	for (ii=0;ii<ntotscan+1; ii++)
		fprintf(fp,"%ld ",frnum[ii]);
	fprintf(fp,"\n");*/

	//
	fwrite(&size,     sizeof(int),  1, fp);
	fwrite(&ntotscan, sizeof(long), 1, fp);

	fwrite(ns,        sizeof(long), ntotscan,fp);
	fwrite(position,       sizeof(long), ntotscan,fp);
	fwrite(frnum,     sizeof(long), ntotscan+1,fp);
	fclose(fp);

	return 0;

}

void read_ParallelizationScheme(string fname,  long **position, long **frnum, long **ns,  long *ntotscan, int *size)
// read the Parrallelization Scheme for further use.

{
	FILE *fp;
	size_t result;

	if ((fp = fopen(fname.c_str(),"r")) == NULL){
		cerr << "Error : couldn't open file to read parallelization Scheme." << endl << "Exiting..." << endl;
		exit(1);
	}

	result = fread(size,      sizeof(int),  1, fp);
	result = fread(ntotscan,  sizeof(long), 1, fp);

	*ns    = new long[*ntotscan];
	*position   = new long[*ntotscan];
	*frnum = new long[(*ntotscan)+1];

	result = fread(*ns,       sizeof(long), *ntotscan,fp);
	result = fread(*position,      sizeof(long), *ntotscan,fp);
	result = fread(*frnum,   sizeof(long), (*ntotscan)+1,fp);
	fclose(fp);

}

void check_ParallelizationScheme(string fname, long *ns, long ntotscan, int size, long **position, long **frnum)
// read and check that the saved Parallelization Scheme corresponds to the actual data
{

	long /**dummy_frnum, *dummy_pos,*/ *dummy_nsamples, dummy_ntotscan;
	int dummy_size;
	read_ParallelizationScheme(fname, position, frnum, &dummy_nsamples, &dummy_ntotscan, &dummy_size);



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


void define_parallelization_scheme(int rank,string fname,long **frnum,long ntotscan,int size,long *nsamples,long *fframes){



	if (rank == 0){

		long *ruleorder ;
		long *fframesorder ;
		long *nsamplesorder ;
		//string *extentnoiseSp_allorder;

		check_ParallelizationScheme(fname,nsamples,ntotscan,size, &ruleorder, frnum);
		// reorder nsamples
		//find_best_order_frames(ruleorder,frnum,nsamples,ntotscan,size);
		//cout << "ruleorder : " << ruleorder[0] << " " << ruleorder[1] << " " << ruleorder[2] << " \n";


		fframesorder  = new long[ntotscan];
		//extentnoiseSp_allorder = new string[ntotscan];
		nsamplesorder = new long[ntotscan];

		for (long ii=0;ii<ntotscan;ii++){
			nsamplesorder[ii] = nsamples[ruleorder[ii]];
			fframesorder[ii] = fframes[ruleorder[ii]];
			//extentnoiseSp_allorder[ii] = extentnoiseSp_all[ruleorder[ii]];
		}
		for (long ii=0;ii<ntotscan;ii++){
			nsamples[ii] = nsamplesorder[ii];
			fframes[ii] = fframesorder[ii];
			//extentnoiseSp_all[ii] = extentnoiseSp_allorder[ii];
			//printf("frnum[%d] = %d\n",ii,frnum[ii]);
		}

		delete [] fframesorder;
		delete [] nsamplesorder;
		//delete [] extentnoiseSp_allorder;

		delete [] ruleorder;


	}else{
		*frnum = new long[ntotscan+1];
	}


}
