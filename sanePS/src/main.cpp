/*
 * main.cpp
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>






#include "imageIO.h"
#include "inline_IO2.h"
#include "mpi_architecture_builder.h"
#include "parsePS.h"
#include "estimPS.h"


extern "C" {
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
}


#ifdef USE_MPI
#include "mpi.h"
#include <algorithm>
#include <fstream>
#endif

using namespace std;



int main(int argc, char *argv[])
{


	// read framesorder

	int size;//,size_det;
	int rank;//,rank_det;
#ifdef USE_MPI
	// int tag = 10;
	//MPI_Status status;

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	cout << size << endl;
	cout << rank << endl;
	if(rank==0)
		printf("\nSanepic Noise Estimation Procedure:\n");

#else
	size = 1;
	rank = 0;
	printf("\nSanepic Noise Estimation Procedure:\n\n");
	cout << "Mpi will not be used for the main loop" << endl;
#endif


	struct param_process proc_param;
	struct samples samples_struct;
	struct param_positions pos_param;
	struct directories dir;
	struct detectors det;

	int flagon;
	long long ind_size;
	long long *indpix;

	// map making parameters


	long NAXIS1, NAXIS2;
	long long npix; // nn = side of the map, npix = number of filled pixels

	//internal data params
	long ns, ff; // number of samples for this scan, first frame number of this scan


	string field; // actual boloname in the bolo loop
	string ellFile; // file containing the ells
	string extentnoiseSp; // noise file
	string prefixe; // prefix used for temporary name file creation


	string MixMatfile = "NOFILE";
	string signame;


	long iframe_min=0, iframe_max=0;

	// main loop variables
	double *S;


	long ncomp;
	double fcut;

	string fname;

	struct wcsprm * wcs;


	int parsed=0;
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {

		parsed=parse_sanePS_ini_file(argv[1], proc_param, dir, samples_struct,det,
				MixMatfile, ellFile, signame, rank, ncomp, fcut);


		if (parsed==-1){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			exit(1);
		}

	}

	long *frames_index;

	frames_index = new long [samples_struct.ntotscan];

	samples_struct.fits_table = new string[samples_struct.ntotscan];
	samples_struct.index_table= new int[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];

	//First time run S=0, after sanepic, S = Pure signal
	if(signame != "NOSIGFILE"){
		// if second launch of estimPS, read S and nn in the previously generated fits map
		cout << "Reading maps ??" << endl;
		S = new double[npix];

		// read pixel indexes
		read_indpix(ind_size, npix, indpix, dir.tmp_dir, flagon);

		read_fits_signal(signame, S, indpix, NAXIS1, NAXIS2, npix);
		FILE * fp;
		fp = fopen("test_signal.txt","w");
		for (int i =0;i<npix;i++)
			fprintf(fp,"%lf\n",S[i]);

		fclose(fp);

		if(rank==0)
			cout << "Map size :" << NAXIS1 << "x" << NAXIS2 << endl;

	}


#ifdef USE_MPI

	ofstream file;

	if(samples_struct.scans_index.size()==0){

		int test=0;
		fname = dir.outdir + parallel_scheme_filename;
		cout << fname << endl;

		test = define_parallelization_scheme(rank,fname,dir.dirfile,samples_struct,size, iframe_min, iframe_max);

		if(test==-1){
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			exit(1);
		}
	}else{
		long size_tmp = 0;
		int return_error = 0;
		int num_frame = 0;
		char c;
		vector2array(samples_struct.scans_index,  samples_struct.index_table);

		if(rank==0){
			//check the processor order given is correct
			//			size_tmp = *max_element(samples_struct.index_table, samples_struct.index_table+samples_struct.ntotscan);

			struct sortclass_int sortobject;
			sort(samples_struct.scans_index.begin(), samples_struct.scans_index.end(), sortobject);

			std::vector<int>::iterator it;
			//			int size_tmp=0;

			// using default comparison:
			it = unique(samples_struct.scans_index.begin(), samples_struct.scans_index.end());
			size_tmp = it - samples_struct.scans_index.begin();

			cout << "size unique : " << size_tmp << endl;

			cout << size << " vs size : " <<  size_tmp << endl;

			if((size_tmp)>size){
				cerr << "Number of processors are different between MPI and parallel scheme. Exiting\n";
				return_error =1;
			}else{

				samples_struct.scans_index.resize( size_tmp );

				cout << "trié + unique : " << samples_struct.scans_index[0] <<  " " << samples_struct.scans_index[1] << endl;


				if((size_tmp)<size){
					cout << "Warning. The number of processors used in fits_filelist is < to the number of processor used by MPI !\n";
					cout << "Do you wish to continue ? (y/n)\n";
					c=getchar();
					switch (c){
					case('y') :
						cout << "Let's continue with only " << (size_tmp) << " processor(s) !\n";
					break;
					default:
						cout << "Exiting ! Please modify fits filelist to use the correct number of processors\n";
						return_error =1;
						break;
					}

					for(long ii=0;ii<size_tmp;ii++)
						if(samples_struct.scans_index[ii]==0)
							num_frame++;

					if(num_frame==0){
						cout << "Exiting ! Please modify fits filelist to use at least processor 0 \n";
						return_error =1;
					}


				}else{


					for(long ii=0;ii<size_tmp;ii++)
						if(samples_struct.scans_index[ii]!=ii){
							cerr << "There is a problem in the fits filelist : you have forgot a processor to use. Exiting" << endl;
							return_error =1;
						}
				}
			}
		}





		if(rank==0){

			string outfile = dir.outdir + samples_struct.filename + "_sanepre.txt";
			cout << "outfile : " << outfile;
			file.open(outfile.c_str(), ios::out);
			if(!file.is_open()){
				cerr << "File [" << fname << "] Invalid." << endl;
				return_error = 1;
			}
		}


		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&return_error,1,MPI_INT,0,MPI_COMM_WORLD);

		if(return_error>0){
			MPI_Finalize();
			exit(0);

		}

		string temp;
		size_t found;

		num_frame=0;
		iframe_min=0;
		iframe_max=0;

		long * nsamples_temp;
		nsamples_temp = new long[samples_struct.ntotscan];

		for(long jj = 0; jj<samples_struct.ntotscan; jj++)
			nsamples_temp[jj]= samples_struct.nsamples[jj];


		for(int ii = 0; ii<size; ii++){
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
	}

	if(rank==0){
		file.close();
		cout << "on aura : \n";
		cout << samples_struct.fits_table[0] << " " << samples_struct.fits_table[1] << " " << samples_struct.fits_table[2] << " " << samples_struct.fits_table[3] << endl;
		cout << samples_struct.noise_table[0] << " " << samples_struct.noise_table[1] << " " << samples_struct.noise_table[2] << " " << samples_struct.noise_table[3] << endl;
		cout << samples_struct.nsamples[0] << " " << samples_struct.nsamples[1] << " " << samples_struct.nsamples[2] << " " << samples_struct.nsamples[3] << endl;
		//cout << samples_struct.filename << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (iframe_max==iframe_min){ // test
		cout << "Warning. Rank " << rank << " will not do anything ! please run saneFrameorder\n";
		//		MPI_Finalize();
		//		exit(0);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	for(long ii=0;ii<size;ii++){
		if(rank==ii)
			cout << "[ " << rank << " ]. iframe min : " << iframe_min << " iframemax : " << iframe_max << endl;
		else
			MPI_Barrier(MPI_COMM_WORLD);
	}

	//////// temp
	//	MPI_Barrier(MPI_COMM_WORLD);
	//	MPI_Finalize();
	//	exit(0);

#else
	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;

	//convert vector to standard C array to speed up memory accesses
	vector2array(samples_struct.noisevect,  samples_struct.noise_table);
	vector2array(samples_struct.fitsvect, samples_struct.fits_table);
	vector2array(samples_struct.scans_index,  samples_struct.index_table);
//	cout << iframe_min << " " << iframe_max << endl;
#endif

	string fits_filename;

	//	if (MixMatfile != "NOFILE"){
	for (long iframe=iframe_min;iframe<iframe_max;iframe++){
		ns = samples_struct.nsamples[iframe];
		ff = iframe;
		extentnoiseSp = samples_struct.noise_table[iframe];
		fits_filename=samples_struct.fits_table[iframe];

		EstimPowerSpectra(proc_param,det,dir, pos_param, ns, ff, NAXIS1,NAXIS2, npix,
				iframe, indpix,	S, MixMatfile, ellFile,
				extentnoiseSp, fits_filename, ncomp, fcut);
		// fsamp = bolometers sampling freq
		// ns = number of samples in the "iframe" scan
		// ff = first sample number
		// ndet = total number of detectors
		// nn = side of the map
		// npix = total number of filled pixels
		// napod = number of border pixels used to apodize data
		// iframe == scan number
		// flgdupl = flagged data map duplication indicator
		// factdupl = duplication factor (1 or 2)
		// indpix = pixels index
		// S = Pnd
		// MixMatfile = this file contains the number of components that interviene in the common-mode component of the noise
		// and the value of alpha, the amplitude factor which depends on detectors but not on time (see formulae (3) in "Sanepic:[...], Patanchon et al.")
		// bolonames = detectors names
		// dirfile = data directory
		// bextension = -B option : "_data" for example
		// fextension = "NOFLAG" or -G option ("_flag" for example)
		// cextension = "NOCALP" or -R option ("_calp" for example)
		// shift_data_to_point (default 0), for subtracting a time offset to the data to match the pointing
		// poutdir = outpout dir or current path (default)
		// termin = output file suffix
		// NORMLIN = baseline is remove from the data, default =0, option -L
		// NOFILLGAP = fill the gap ? default yes => 0
		// tmp_dir = noise power spectrum file suffix = path
		// extentnoiseSp = noise file
		// outdir = output directory
	}
	//	}




#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

	//clean up
	delete [] samples_struct.nsamples;
	delete [] S;

	delete [] samples_struct.fits_table;
	delete [] samples_struct.index_table;
	delete [] samples_struct.noise_table;

	delete [] frames_index;

	delete [] indpix;
	if(signame == "NOSIGFILE"){
		int nwcs = 1;
		wcsvfree(&nwcs, &wcs);
	}
	printf("\nEnd of sanePS\n");
	return 0;
}
