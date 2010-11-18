#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <list>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <sysexits.h>

#include "mpi_architecture_builder.h"
#include "struct_definition.h"

#include "dataIO.h"
#include "imageIO.h"
#include "temporary_IO.h"
#include "inputFileIO.h"
#include "parser_functions.h"
#include "map_mak_tools.h"

extern "C" {
#include "nrutil.h"
}

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;


void print_struct(struct param_sanePre proc_param, struct samples samples_struct, struct param_sanePos pos_param, struct param_common dir,
		struct param_saneInv saneInv_struct, struct param_sanePic struct_sanePic, struct param_sanePS structPS);

int main(int argc, char *argv[])
{

	int size;
	int rank;
#ifdef USE_MPI
	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#else
	size = 1;
	rank = 0;
#endif


	struct param_sanePre proc_param; /*! A structure that contains user options about preprocessing properties */
	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_sanePos pos_param; /*! A structure that contains user options about map projection and properties */
	struct param_common dir; /*! structure that contains output input temp directories */

	string field; /*! actual boloname in the bolo loop */
	struct param_sanePic struct_sanePic;
	std::vector<double> fcut; /*! noise cutting frequency vector */

	// those variables will not be used by sanePic but they are read in ini file (to check his conformity)
	struct param_sanePS structPS;
	struct param_saneInv saneInv_struct;
	string parser_output = "";


	// parser variables
	int parsed = 0;

	if (argc < 2) // no enough arguments
		parsed = 1;
	else {

		/* parse ini file and fill structures */
		parsed = parser_function(argv[1], parser_output, dir,
				samples_struct, pos_param, proc_param, fcut, structPS, saneInv_struct,
				struct_sanePic, size);

		if(rank==0)
			// print parser warning and/or errors
			cout << endl << parser_output << endl;


	}

	if (parsed > 0) { // error during parser phase
		if (rank == 0)
			switch (parsed) {

			case 1:
				printf(
						"Please run %s using the following options : sanepic_ini.ini \n",
						argv[0]);
				break;

			case 2:
				printf("Wrong program options or argument. Exiting !\n");
				break;

			case 3:
				printf("Exiting...\n");
				break;

			default:
				;
			}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_CONFIG;
	}

	// REORDER WITH MPI_ARCHI

#ifdef USE_MPI

	ofstream file;
	long iframe_min, iframe_max;

	// User has not given a processor order
	if(samples_struct.scans_index.size()==0){

		int test=0;
		string fname = dir.output_dir + parallel_scheme_filename;
		if(rank==0)
			cout << "Getting configuration and frame order from file : " << fname << endl;
		test = define_parallelization_scheme(rank,fname,dir.input_dir,samples_struct,size, iframe_min, iframe_max);

		if(test==-1){
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
			return EX_CONFIG;
		}

	}else{ // user has given a processor order
		int test=0;
		test = verify_parallelization_scheme(rank,dir.output_dir,samples_struct, size, iframe_min, iframe_max);


		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&test,1,MPI_INT,0,MPI_COMM_WORLD);

		if(test>0){
			MPI_Finalize();
			return EX_CONFIG;

		}

	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (iframe_max==iframe_min){
		cout << "Warning. Rank " << rank << " will not do anything ! please run saneFrameorder\n";

	}

	MPI_Barrier(MPI_COMM_WORLD);

#else

	//	vector2array(samples_struct.fitsvect, samples_struct.fits_table);
	//	vector2array(samples_struct.scans_index,  samples_struct.index_table);


#endif





	if(read_bolo_for_all_scans(dir, samples_struct, rank, size)){
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return(EX_IOERR);
	}


	fill_noisevect_fcut(dir, samples_struct, saneInv_struct, fcut);
	//	vector2array(samples_struct.noisevect,  samples_struct.noise_table);

	fill_sanePS_struct(structPS, samples_struct);

	if(rank==0)
		// parser print screen function
		print_struct(proc_param, samples_struct, pos_param, dir, saneInv_struct,
				struct_sanePic,structPS);

	for(long iframe=0; iframe< samples_struct.ntotscan; iframe ++){
		long ns = samples_struct.nsamples[iframe];
		double f_lppix = proc_param.f_lp*double(ns)/proc_param.fsamp;string output_read = "";
		std::vector<string> det_vect;
		if(read_channel_list(output_read, samples_struct.bolovect[iframe], det_vect)){
			cout << output_read << endl;
#ifdef USE_MPI
			//			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			return EX_CONFIG;
		}
		long ndet = (long)det_vect.size();
		if(test_map_mak(samples_struct,  ns, proc_param, dir.output_dir, det_vect, ndet, f_lppix, iframe, rank, size)){
			cout << "ERROR IN test_map_mak... EXITING ... \n";
			return EXIT_FAILURE;
		}
	}

	cout << "\n\nEXIT SUCCESS  !!!\n";
	return EXIT_SUCCESS;
}



void print_struct(struct param_sanePre proc_param, struct samples samples_struct, struct param_sanePos pos_param, struct param_common dir,
		struct param_saneInv saneInv_struct, struct param_sanePic struct_sanePic, struct param_sanePS structPS)
{

	cout << "\ndir : struct common\n";
	print_common(dir);

	cout << "\npos_param : struct param_sanePos\n";
	cout << "maskfile = "  << pos_param.maskfile << endl;
	cout << "pixdeg = " << pos_param.pixdeg << endl;
	cout << "flgdupl = " << pos_param.flgdupl << endl;
	cout << "projgaps = " << pos_param.projgaps << endl;
	cout << "fileFormat = " << pos_param.fileFormat << endl;

	cout << "\nproc_param : struct param_sanePre\n";
	cout << "NORMLIN = " << proc_param.NORMLIN << endl;
	cout << "NOFILLGAP = " << proc_param.NOFILLGAP << endl;
	cout << "CORRon = " << proc_param.CORRon << endl;
	cout << "remove_polynomia = " << proc_param.remove_polynomia << endl;
	cout << "napod = " << proc_param.napod << endl;
	cout << "poly_order = " << proc_param.poly_order << endl;
	cout << "fsamp = " << proc_param.fsamp << endl;
	cout << "f_lp = " << proc_param.f_lp << endl;
	cout << "fcut_file = " << proc_param.fcut_file << endl;

	cout << "\nsaneInv_struct : struct param_saneInv\n";
	cout << "cov_matrix_file = " << saneInv_struct.cov_matrix_file << endl;
	cout << "cov_matrix_suffix = " << saneInv_struct.cov_matrix_suffix << endl;

	cout << "\nstruct_sanePic : struct param_sanePic\n";
	cout << "iterW = " << struct_sanePic.iterw << endl;

	cout << "\nstructPS : struct param_sanePS\n";
	cout << "fcutPS = " << structPS.fcutPS << endl;
	cout << "ell_suffix = " << structPS.ell_suffix << endl;
	cout << "mix_suffix = " << structPS.mix_suffix << endl;
	cout << "ell_global_file = " << structPS.ell_global_file << endl;
	cout << "mix_global_file = " << structPS.mix_global_file << endl;
	cout << "signame = " << structPS.signame << endl;
	cout << "ncomp = " << structPS.ncomp << endl;

	cout << "\nsamples_struct : struct samples\n";
	cout << "filename = " << samples_struct.filename << endl;
	cout << "ntotscan = " << samples_struct.ntotscan << endl;
	cout << "fitsvect = " << endl;
	for(long ii=0; ii< (long)samples_struct.fitsvect.size(); ii++)
		cout << samples_struct.fitsvect[ii] << endl;
#ifdef USE_MPI
	cout << "index_table = " << endl;
	for(long ii=0; ii< samples_struct.ntotscan; ii++)
		cout << samples_struct.index_table[ii] << endl;
#endif
	cout << "noisevect = " << endl;
	for(long ii=0; ii< (long)samples_struct.noisevect.size(); ii++)
		cout << samples_struct.noisevect[ii] << endl;
	cout << "framegiven = " << samples_struct.framegiven << endl;
#ifdef USE_MPI
	cout << "fits_table = " << endl;
	for(long ii=0; ii< samples_struct.ntotscan; ii++)
		cout << samples_struct.fits_table[ii] << endl;
	cout << "noise_table = " << endl;
	for(long ii=0; ii< samples_struct.ntotscan; ii++)
		cout << samples_struct.noise_table[ii] << endl;
#endif

	printf("Number of scans      : %ld\n",samples_struct.ntotscan);
	//	printf("Number of bolometers : \n");
	//	for(long iframe=0;iframe<samples_struct.ntotscan;iframe++){
	//		if(samples_struct.fits_table!=NULL)
	//			printf("Scan number %ld : %s %ld\n", iframe,(char*)(FitsBasename(samples_struct.fits_table[iframe]).c_str()), detector_tab[iframe].ndet);
	//		else
	//			printf("Scan number %ld : %s %ld\n", iframe,(char*)(FitsBasename(samples_struct.fitsvect[iframe]).c_str()), detector_tab[iframe].ndet);
	//	}


}
