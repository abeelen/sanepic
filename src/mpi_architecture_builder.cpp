/*
 * mpi_architecture_builder.cpp
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */


#include "mpi_architecture_builder.h"

using namespace std;

double *dat_compare;

int compare_array_double (const void *a, const void *b)
{

	const long *da = (const long *) a;
	const long *db = (const long *) b;

	return (dat_compare[*da] > dat_compare[*db]) - (dat_compare[*da] < dat_compare[*db]);
}

void find_best_order_frames(long *pos, long *frnum, long *ns, int ntotscan, int size){

	long ii, jj, kk, count, a, b;
	double temp, maxproctmp, valmin, stdmin, stdtmp, valtmp;
	int tmpposmin, tmpposmax;

	long *ns_order;
	double *std, *sizeperproc, *maxproc;
	double *valtemp;

	int nessai = 10000;

	long ntot = 0; // total number of frames (*20 to get sample)

	for (ii=0;ii<ntotscan;ii++)
		ntot += ns[ii];

	maxproc = new double[nessai];
	std = new double[nessai];
	ns_order = new long[ntotscan];
	sizeperproc = new double[ntotscan];
	dat_compare = new double[ntotscan];

	//init random generator
	valtemp = randg(1,0); // valtemp is a random double between 0/1. with 0 the seed of random function is fixed

	//init arrays
	for (ii=0;ii<=ntotscan;ii++)
		frnum[ii] = 0;


	for (jj=0;jj<nessai;jj++){

		for (kk=0;kk<ntotscan;kk++){
			valtemp = randg(1,-1); // return a random value between 0/1
			dat_compare[kk] = valtemp[0];
		}

		for (ii=0;ii<ntotscan;ii++)
			ns_order[ii] = ns[ii];
		for (ii=0;ii<ntotscan;ii++)
			pos[ii] = ii;

		qsort(pos,ntotscan,sizeof(long),compare_array_double);

		for (ii=0;ii<ntotscan;ii++)
			ns_order[ii] = ns[pos[ii]];

		count = 0;
		a = ns_order[0];
		for (ii=1;ii<=ntotscan;ii++){
			b = 0;
			for(kk=0;kk<ii;kk++)
				b += ns_order[kk];
			if (abs(b-ntot*double(count+1)/double(size)) <= abs(a-ntot*double(count+1)/double(size))){
				a = b;
			} else {
				count += 1;
				frnum[count] = ii-1;
				sizeperproc[count-1] = 0.0;
				for (kk=frnum[count-1];kk<ii-1;kk++)
					sizeperproc[count-1] += double(ns_order[kk]);
				a = b;
			}
		}
		frnum[size] = ntotscan;
		sizeperproc[count] = 0.0;
		for (kk=frnum[count];kk<ntotscan;kk++)
			sizeperproc[count] += double(ns_order[kk]);


		//*********** check values
		minmax(sizeperproc,ntotscan,&temp,&maxproctmp,&tmpposmin,&tmpposmax);
		maxproc[jj] = maxproctmp;
		//seeds(jj+1,*) = seed;
		std[jj] = 0.0;
		for(kk=0;kk<ntotscan;kk++)
			if (sizeperproc[kk] > 0.5)
				std[jj] += (sizeperproc[kk]-double(ntot)/size)*(sizeperproc[kk]-double(ntot)/size)/size;


	}// end of first ntotscan loop


	minmax(maxproc,nessai,&valmin,&temp,&tmpposmin,&tmpposmax);

	stdmin = double(ntot*ntot);
	for (ii=0;ii<nessai;ii++)
		if (long(valmin) == long(maxproc[ii]))
			if (std[ii] < stdmin)
				stdmin = std[ii];



	valtmp = 2.0*valmin;
	stdtmp = 2.0*stdmin;
	while ((stdtmp > stdmin) || ((long)valtmp > (long)valmin)){


		for (kk=0;kk<ntotscan;kk++){
			valtemp = randg(1,-1);
			dat_compare[kk] = valtemp[0];
		}

		for (ii=0;ii<ntotscan;ii++)
			ns_order[ii] = ns[ii];
		for (ii=0;ii<ntotscan;ii++)
			pos[ii] = ii;

		qsort(pos,ntotscan,sizeof(long),compare_array_double);

		for (ii=0;ii<ntotscan;ii++)
			ns_order[ii] = ns[pos[ii]];

		count = 0;
		a = ns_order[0];
		for (ii=1;ii<=ntotscan;ii++){
			b = 0;
			for(kk=0;kk<ii;kk++)
				b += ns_order[kk];
			if (abs(b-ntot*double(count+1)/double(size)) <= abs(a-ntot*double(count+1)/double(size))){
				a = b;
			} else {
				count += 1;
				frnum[count] = ii-1;
				sizeperproc[count-1] = 0.0;
				for (kk=frnum[count-1];kk<ii-1;kk++)
					sizeperproc[count-1] += double(ns_order[kk]);
				a = b;
			}
		}
		frnum[size] = ntotscan;
		sizeperproc[count] = 0.0;
		for (kk=frnum[count];kk<ntotscan;kk++)
			sizeperproc[count] += double(ns_order[kk]);


		//*********** check values
		minmax(sizeperproc,ntotscan,&temp,&maxproctmp,&tmpposmin,&tmpposmax);
		valtmp = maxproctmp;
		stdtmp = 0.0;
		for(kk=0;kk<ntotscan;kk++)
			if (sizeperproc[kk] > 0.5)
				stdtmp += (sizeperproc[kk]-double(ntot)/size)*(sizeperproc[kk]-double(ntot)/size)/size;

	}

	printf("max range = %lf, std range = %lf\n",valtmp,sqrt(stdtmp));


}
