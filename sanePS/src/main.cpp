#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>






#include "imageIO.h"
#include "inline_IO2.h"
#include "mpi_architecture_builder.h"
#include "parser_functions.h"
#include "estimPS.h"
#include "struct_definition.h"
#include "inputFileIO.h"

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
	//	cout << size << endl;
	//	cout << rank << endl;
	if(rank==0)
		printf("\nSanepic Noise Estimation Procedure:\n");

#else
	size = 1;
	rank = 0;
	printf("\nSanepic Noise Estimation Procedure:\n\n");
	cout << "Mpi will not be used for the main loop" << endl;
#endif


	struct param_process proc_param; /*! A structure that contains user options about preprocessing properties */
	struct samples samples_struct;  /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_positions pos_param; /*! A structure that contains user options about map projection and properties */
	struct common dir; /*! structure that contains output input temp directories */
	//	struct detectors det; /*! A structure that contains everything about the detectors names and number */
	std::vector<detectors> detector_tab;

	int flagon; /*!  if one sample is rejected, flagon=1 */
	long long ind_size; // indpix size
	long long *indpix; // map index

	// map making parameters


	long NAXIS1, NAXIS2; // map size
	long long npix; // npix = number of filled pixels

	//internal data params
	long ns; // number of samples for this scan


	string field; // actual boloname in the bolo loop
	string prefixe; // prefix used for temporary name file creation


	string MixMatfile = "NOFILE"; // mixing matrix file
	string signame; // map filename


	long iframe_min=0, iframe_max=0;

	// main loop variables
	double *S = NULL; // signal


	long ncomp; // number of noise component to estimate
	double fcut; // cut-off freq : dont focus on freq larger than fcut for the estimation !

	int parsed=0; // parser error code
	if (argc<2) { // too few arguments
		printf("Please run %s using a *.ini file\n",argv[0]);
		parsed=1;
	} else {

		// those variables will not be used by sanePre but they are read in ini file (to check his conformity)
		int iterw=10;
		int save_data, load_data;
		std::vector<double> fcut_vector;

		/* parse ini file and fill structures */
		parsed=parser_function(argv[1], dir, detector_tab, samples_struct, pos_param, proc_param, fcut_vector,
				fcut, MixMatfile, signame, ncomp, iterw, save_data, load_data, rank, size);
	}

	if (parsed>0){ // error during parser phase
		if (rank==0)
			switch (parsed){

			case 1: printf("Please run %s using a *.ini file\n",argv[0]);
			break;

			case 2 : printf("Wrong program options or argument. Exiting !\n");
			break;

			case 3 : printf("Exiting...\n");
			break;

			default :;
			}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return 1;
	}

	samples_struct.fits_table = new string[samples_struct.ntotscan];
	samples_struct.index_table= new int[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];

	//First time run S=0, after sanepic, S = Pure signal
	if(signame != "NOSIGFILE"){

		//TODO : Add some check for the map size/ind_size/npix

		if(read_indpix(ind_size, npix, indpix, dir.tmp_dir, flagon)){ // read map indexes
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			return(EXIT_FAILURE);
		}

		// if second launch of estimPS, read S and NAXIS1/2 in the previously generated fits map
		if(rank==0)
			cout << "Reading model map : " << signame << endl;
		S = new double[npix]; // pure signal

		// read pure signal
		if(read_fits_signal(signame, S, indpix, NAXIS1, NAXIS2)){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			return(EXIT_FAILURE);
		}

#ifdef DEBUG
		FILE * fp;
		fp = fopen("reconstructed_1dsignal.txt","w");
		for (int i =0;i<npix;i++)
			fprintf(fp,"%lf\n",S[i]);
		fclose(fp);
#endif

		if(rank==0)
			cout << "         map size : " << NAXIS1 << "x" << NAXIS2 << endl;

	}


#ifdef USE_MPI

	ofstream file;

	if(samples_struct.scans_index.size()==0){

		int test=0;
		string fname = dir.output_dir + parallel_scheme_filename;
		//		cout << fname << endl;

		test = define_parallelization_scheme(rank,fname,dir.dirfile,samples_struct,size, iframe_min, iframe_max);

		if(test==-1){
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			exit(1);
		}
	}else{
		int test=0;
		test = verify_parallelization_scheme(rank,dir.output_dir,samples_struct, size, iframe_min, iframe_max);


		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&test,1,MPI_INT,0,MPI_COMM_WORLD);

		if(test>0){
			MPI_Finalize();
			exit(0);

		}

	}


	MPI_Barrier(MPI_COMM_WORLD);

	if (iframe_max==iframe_min){ // test
		cout << "Warning. Rank " << rank << " will not do anything ! please run saneFrameorder\n";

	}

	MPI_Barrier(MPI_COMM_WORLD);

	//	for(long ii=0;ii<size;ii++){
	//		if(rank==ii)
	//			cout << "[ " << rank << " ]. iframemin : " << iframe_min << " iframemax : " << iframe_max << endl;
	//		else
	//			MPI_Barrier(MPI_COMM_WORLD);
	//	}
#else
	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;

	//convert vector to standard C array to speed up memory accesses
	vector2array(samples_struct.noisevect,  samples_struct.noise_table);
	vector2array(samples_struct.fitsvect, samples_struct.fits_table);
	vector2array(samples_struct.scans_index,  samples_struct.index_table);
#endif

	string fits_filename; // input scan filename (fits file)
	std::ostringstream oss; // we need to store the string in a stringstream before using basename

	for (long iframe=iframe_min;iframe<iframe_max;iframe++){ // proceed scan by scan
		ns = samples_struct.nsamples[iframe];
		fits_filename=samples_struct.fits_table[iframe];
		cout << "[ " << rank << " ] " << fits_filename << endl;
		oss << fits_filename;
		string filename = oss.str();
		MixMatfile=Basename(filename);
		oss.str("");

		struct detectors det = detector_tab[iframe];

		EstimPowerSpectra(proc_param,det,dir, pos_param, ns, NAXIS1,NAXIS2, npix,
				iframe, indpix,	S, MixMatfile, MixMatfile,
				fits_filename, ncomp, fcut, rank);
		// ns = number of samples in the "iframe" scan
		// npix = total number of filled pixels
		// iframe = scan number
		// indpix = pixels index
		// MixMatfile = this file contains the number of components in the common-mode component of the noise
		// and the value of alpha, the amplitude factor which depends on detectors but not on time (see formulae (3) in "Sanepic:[...], Patanchon et al.")
	}




#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

	//clean up
	delete [] samples_struct.nsamples;
	delete [] samples_struct.fits_table;
	delete [] samples_struct.index_table;
	delete [] samples_struct.noise_table;



	if(signame == "NOSIGFILE"){
		//		int nwcs = 1;
		//		wcsvfree(&nwcs, &wcs);
	}else{
		delete [] S;
		delete [] indpix;
	}
	printf("\nEnd of sanePS\n");
	return 0;
}
