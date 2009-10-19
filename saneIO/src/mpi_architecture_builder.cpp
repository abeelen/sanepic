#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <string.h>

#include "time.h"
#include "positionsIO.h"
#include "mpi_architecture_builder.h"

extern "C" {
#include <fitsio.h>
}


//using namespace std;


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


	cout << ntot << endl;

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
	stdmin=double(ntot);
	stdmin=stdmin*stdmin;
	for (int ii=0;ii<nessai;ii++)
		if (long(valmin) == long(maxproc[ii]))
			if (std[ii] < stdmin)
				stdmin = std[ii];



	valtmp = 2.0*valmin;
	stdtmp = 2.0*stdmin;
	printf("max range = %lf, std range = %lf\n",valtmp,sqrt(stdtmp));
	getchar();

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

		printf("max range = %lf, std range = %lf\n",valtmp,sqrt(stdtmp));
		getchar();
	}

	delete [] maxproc;
	delete [] std;
	delete [] ns_order;
	delete [] sizeperproc;
	delete [] dat_compare;

	printf("max range = %lf, std range = %lf\n",valtmp,sqrt(stdtmp));

}

int write_ParallelizationScheme(string fname, long *position, long *frnum, long *ns, long ntotscan, int size,
		std::vector<string> fitsvect, std::vector<string> noisevect, std::vector<long> &scans_index)
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

	cout << "file " << fname << endl;

	//std::vector<string> fitsvect_temp;
	//std::vector<string> noisevect_temp;
	//std::vector<long> scans_index_temp;
	string *fitsvect_temp, *noisevect_temp;
	long *scans_index_temp;
	string temp;
	size_t found;

	fitsvect_temp = new string [ntotscan];
	noisevect_temp = new string [ntotscan];
	scans_index_temp = new long [ntotscan];

	//cout << "write" << endl;
	//cout << fitsvect[0] << " "  << fitsvect[1] << " "  << fitsvect[2] << " "  << fitsvect[3] << endl;
	//cout << noisevect[0] << " "  << noisevect[1] << " "  << noisevect[2] << " "  << noisevect[3] << endl;

	for (long ii=0;ii<ntotscan;ii++){
		cout << "loop " << ii << " : " << endl;
		cout << position[ii] << endl;
		temp = fitsvect[position[ii]];
		found=temp.find_last_of("/");
		// cout << " file: " << str.substr(found+1) << endl;


		fitsvect_temp[ii] = temp.substr(found+1);
		noisevect_temp[ii] = noisevect[position[ii]];
		//cout << fitsvect_temp[ii] << " " << noisevect_temp[ii] << endl;
		// scans_index_temp[ii] = scans_index[position[ii]];
	}

	long val_proc = 0;
	//scans_index[0] = 0;

	for (long ii=1;ii<ntotscan+1;ii++){
		if(frnum[ii]==0)
			break;

		for(int jj=frnum[ii-1];jj<frnum[ii];jj++)
			scans_index_temp[jj]=val_proc;
		val_proc++;
	}

	cout << "nb proc : " << val_proc << endl;

	if(val_proc>size){
		cerr << "Error in frame order repartition, number of processor are not equal to mpi size\n";
		return -1;
	}

	for (long ii=0;ii<ntotscan;ii++){
		// fitsvect[ii] = fitsvect_temp[ii];
		//noisevect[ii] = noisevect_temp[ii];
		//scans_index[ii] = scans_index_temp[ii];
		file << fitsvect_temp[ii] << " " << noisevect_temp[ii] << " " << scans_index_temp[ii] << endl;
		//cout << fitsvect_temp[ii] << " " << noisevect_temp[ii] << " " << endl;
	}

	file.close();

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
	//fwrite(&size,     sizeof(int),  1, fp);
	//fwrite(&ntotscan, sizeof(long), 1, fp);

	//fwrite(ns,        sizeof(long), ntotscan,fp);
	//fwrite(position,       sizeof(long), ntotscan,fp);
	//fwrite(frnum,     sizeof(long), ntotscan+1,fp);




	//fclose(fp);

	delete [] fitsvect_temp;
	delete [] noisevect_temp;

	cout << "fin fonction\n";

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

int check_ParallelizationScheme(string fname, string dirfile,long ntotscan, int size, long *&nsamples, std::vector<string> fitsfiles, std::vector<string> noisefiles, string *&fits_table, string *&noise_table, long *&index_table)
// read and check that the saved Parallelization Scheme corresponds to the actual data
{




	std::vector<string> fits_dummy;
	std::vector<string> noise_dummy;
	std::vector<long> index_dummy;
	long ntotscan_dummy;
	long size_tmp;
	//size_t found;

	bool framegiven;
	long *fframes, *nsamples_dummy;
	string temp;

	//string *fits_table, *noise_table;
	//long *index_table;



	read_fits_list(fname, fits_dummy, noise_dummy, index_dummy, framegiven);

	cout <<" readed list : " << endl;
	cout << fits_dummy[0] << " " << fits_dummy[1] << " " << fits_dummy[2] << " " << fits_dummy[3] << endl;
	cout <<  noise_dummy[0] << " " <<  noise_dummy[1] << " " <<  noise_dummy[2] << " " <<  noise_dummy[3] << endl;
	cout <<   index_dummy[0] << " " <<  index_dummy[1] << " " <<   index_dummy[2] << " " <<  index_dummy[3] << endl;
	cout << "framegiven : " << framegiven << endl;

	if((framegiven==0)||((int)fits_dummy.size()==0))
		return -1;

	ntotscan_dummy=(long)fits_dummy.size();
	if(ntotscan!=ntotscan_dummy){
		cerr << "number of scans are different between your fits file list and the mpi scheme" << endl;
		return -1;
	}

	vector2array(fits_dummy, fits_table);
	vector2array(noise_dummy, noise_table);
	vector2array(index_dummy,  index_table);

	for(int ii=0;ii<(int)fits_dummy.size();ii++){
	  //	cout << dirfile + fits_dummy[ii] << endl;
		fits_dummy[ii] = dirfile + fits_dummy[ii];
	}

	readFrames( &ntotscan_dummy , fits_dummy, fframes, nsamples_dummy);
	delete [] fframes;

	cout << "ntotscan" << endl;
	cout << ntotscan << " vs " << ntotscan_dummy << endl;

	/*if(ntotscan!=ntotscan_dummy){
		cerr << "number of scans are different between your fits file list and the mpi scheme" << endl;
		return -1;
		}*/

	//fits_table = new string[ntotscan];
	//index_table= new long[ntotscan];
	//noise_table = new string[ntotscan];



	/*	for(int ii=0;ii<ntotscan;){
	  temp = fits_table[ii];
	  cout << "temp " << temp << endl;
	  found=temp.find_last_of("/");
	  fits_table[ii] = temp.substr(found+1);
	  cout << "fits_table " << fits_table[ii] << endl;
	  }*/

	struct sortclass_string sortobject;
	sort(fits_dummy.begin(), fits_dummy.end(), sortobject);
	sort(fitsfiles.begin(), fitsfiles.end(), sortobject);

	cout << "comparaison triée : " << endl;

	for(int ii=0;ii<ntotscan;ii++)
	  if(fits_dummy[ii]!=fitsfiles[ii]){
	    cout << fits_dummy[ii] << endl;
	    cout << fitsfiles[ii] << endl;
	    return -1;
	  }


	if ((int)noisefiles.size()>0){

	cout << "comparaison triée bruit : " << endl;

	sort (noise_dummy.begin(), noise_dummy.end(), sortobject);
	sort (noisefiles.begin(), noisefiles.end(), sortobject);

	for(int ii=0;ii<ntotscan;ii++)
	  if(noise_dummy[ii]!=noisefiles[ii]){
	    cout << noise_dummy[ii] << endl;
	    cout << noisefiles[ii] << endl;
	    return -1;
	  }
	}

	//struct sortclass_int sortobject;
	//sort(index_dummy.begin(), index_dummy.end(), sortobject);
	size_tmp = *max_element(index_table, index_table+ntotscan);

	cout << size << " vs size : " <<  size_tmp+1 << endl;

	if((size_tmp+1)!=size){
	  cerr << "Number of processors are different between MPI and parallel scheme. Exiting\n";
	  return -1;
	}


	for(int ii=0;ii<ntotscan;ii++)
	  nsamples[ii]=nsamples_dummy[ii];



	//exit(0);


	return 0;
}



int define_parallelization_scheme(int rank,string fname,string dirfile,long ntotscan,int size,long *&nsamples, std::vector<string> fitsfiles, std::vector<string> noisefiles, string *&fits_table, string *&noise_table, long *&index_table){


  // cout << "avant check" << endl;
  cout << "rank" << rank << endl;
  // if (rank == 0){
    //cout << "avant check" << endl;
    int test=0;
    //long *ruleorder ;
    //long *fframesorder ;
    //long *nsamplesorder ;
    //string *extentnoiseSp_allorder;

    test=check_ParallelizationScheme(fname,dirfile,ntotscan,size,nsamples,fitsfiles,noisefiles,fits_table,noise_table,index_table);
    if (test==-1)
      return test;
    // reorder nsamples
    //find_best_order_frames(ruleorder,frnum,nsamples,ntotscan,size);
	  //cout << "ruleorder : " << ruleorder[0] << " " << ruleorder[1] << " " << ruleorder[2] << " \n";


    // }else{
    //	*frnum = new long[ntotscan+1];
    // }

  return 0;

}

/*int define_parallelization_scheme(int rank,string fname,long **frnum,long ntotscan,int size,long *nsamples,long *fframes){



	if (rank == 0){

		int test=0;
		long *ruleorder ;
		long *fframesorder ;
		long *nsamplesorder ;
		//string *extentnoiseSp_allorder;

		test=check_ParallelizationScheme(fname,nsamples,ntotscan,size, &ruleorder, frnum);
		if (test==-1)
			return test;
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

	return 0;

}*/

long readFitsLength(string filename){

	fitsfile *fptr;
	int status = 0;
	int hdu_type;
	long ns;

	//	Open the fits file
	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// Go to the first Extension (should be an image named "Primary")
	if (fits_movabs_hdu(fptr, 1, &hdu_type,  &status))
		fits_report_error(stderr, status);
	if (hdu_type != IMAGE_HDU)
		fits_report_error(stderr, BAD_HDU_NUM);

	if(fits_read_key(fptr, TLONG, (char *) "NSAMP", &ns, NULL, &status))
		fits_report_error(stderr ,status);

	return ns;

}

void readFrames(long * nScan, std::vector<string> &inputList, long *& fframes, long *& nsamples){

	//read_strings(filename, inputList);

	*nScan    = inputList.size();
	fframes  = new long[*nScan];
	nsamples = new long[*nScan];
	for (long i=0; i<*nScan; i++){
		fframes[i]  = i;
		nsamples[i] = readFitsLength(inputList[i]);
	}

}



void read_fits_list(string fname, std::vector<string> &fitsfiles, std::vector<string> &noisefiles, std::vector<long> &frameorder, bool &framegiven) {



	ifstream file;
	file.open(fname.c_str(), ios::in);
	if(!file.is_open()){
		cerr << "File [" << fname << "] Invalid." << endl;
		exit(-1);
	}


	framegiven=0;

	string s,p, line, temp;
	long d;
	char *pch;
	int nb_elem = 0;
	//int num=0; // framecounter

	// count number of elements on the first line !
	getline(file, line);
	line.erase(0, line.find_first_not_of(" \t")); // remove leading white space
	pch = strtok ((char*) line.c_str()," ,");

	while (pch != NULL) {
		pch = strtok (NULL, " ,-");
		nb_elem++; 	}

	// set pointer back to the beginning of file in order to parse the first line too
	file.seekg (0, ios::beg);

	switch(nb_elem) {
	case 3:
		framegiven=1;
		while(file >> s >> p >> d) {
			size_t found;
			s.erase(0, s.find_first_not_of(" \t")); // remove leading white space in the first name
			found = s.find_first_of("!#;"); 		// Check for comment character at the beginning of the filename
			if (found == 0) continue;

			cout << "3 : " << s << " " << p << " " << d << endl;
			fitsfiles.push_back(s);
			noisefiles.push_back(p);
			frameorder.push_back(d);
		}
		break;

	case 2:
		while(file >> s >> p){
			size_t found;
			s.erase(0, s.find_first_not_of(" \t")); // remove leading white space in the first name
			found = s.find_first_of("!#;"); 		// Check for comment character at the beginning of the filename

			if (found == 0) continue;

			cout << "2 : " << s << " " << p << endl;
			fitsfiles.push_back(s);
			noisefiles.push_back(p);}
		//frameorder.push_back(num++); }
		break;

	case 1:
		//noisefiles.push_back("read_file_in_the_ini");

		while(file >> s){
			size_t found;
			s.erase(0, s.find_first_not_of(" \t")); // remove leading white space in the first name
			found = s.find_first_of("!#;"); 		// Check for comment character at the beginning of the filename
			if (found == 0) continue;

			cout << "1 : " << s << endl;
			fitsfiles.push_back(s);}
		//frameorder.push_back(num++); }
		break;

	default:
		cerr << "File [" << fname << "] must have at least one row and 2 colums. Exiting\n";
		exit(0);
		break;
	}

	if(fitsfiles.size()==0){
		cerr << "File [" << fname << "] must have at least one row. Exiting\n";
		exit(0);
	}

	if (file>>s){
		cerr << "File [" << fname << "]. Each line must have the same number of rows. Exiting\n";
		exit(0);
	}

	cout << "read fits list ok !!!\n";
	file.close();

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

}
