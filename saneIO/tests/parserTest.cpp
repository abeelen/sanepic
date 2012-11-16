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

#include <limits>

#include "mpi_architecture_builder.h"
#include "struct_definition.h"

#include "dataIO.h"
#include "imageIO.h"
#include "temporary_IO.h"
#include "inputFileIO.h"
#include "parser_functions.h"

extern "C" {
#include "nrutil.h"
#include "getdata.h"
}

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;


void print_struct(struct param_saneProc proc_param, struct samples samples_struct, struct param_sanePos pos_param, struct param_common dir,
		struct param_saneInv inv_param, struct param_sanePic pic_param, struct param_sanePS ps_param);

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


	struct param_saneProc proc_param; /*! A structure that contains user options about preprocessing properties */
	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_sanePos pos_param; /*! A structure that contains user options about map projection and properties */
	struct param_common dir; /*! structure that contains output input temp directories */

	string field; /*! actual boloname in the bolo loop */
	struct param_sanePic pic_param;

	// those variables will not be used by sanePic but they are read in ini file (to check his conformity)
	struct param_sanePS ps_param;
	struct param_saneInv inv_param;
	string parser_output = "";


	// parser variables
	int parsed = 0;

	if (argc < 2) // no enough arguments
		parsed = 1;
	else {

		/* parse ini file and fill structures */
		parsed = parser_function(argv[1], parser_output, dir,
				samples_struct, pos_param, proc_param, ps_param, inv_param,
				pic_param, size, rank);


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


	if(rank==12){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				ps_param, pic_param, inv_param);

	}


	// this should be done by subrank 0 only
//	init_dirfile(dir.tmp_dir, samples_struct, pos_param.fileFormat, rank);
	//	fill_sanePS_struct(ps_param, samples_struct, dir);

	cout << rank << " parser done." << endl;


#ifdef PARA_FRAME
	//	MPI_Barrier(MPI_COMM_WORLD);

	if(configure_PARA_FRAME_samples_struct(dir.tmp_dir, samples_struct, rank, size)){
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return EX_IOERR;
	}

	MPI_Barrier(MPI_COMM_WORLD);

#endif

	readFramesFromFits(samples_struct, rank);

//	get_noise_bin_sizes(dir.tmp_dir, samples_struct, rank);

	if (rank==0){
		//		print_common(dir);
		//		print_param_sanePos(pos_param);
		//		print_param_saneProc(proc_param);
		//		print_param_saneInv(inv_param);
		//		print_param_sanePic(pic_param);
		//		print_param_sanePS(ps_param);

		cout << " ---" << endl;
		print_struct(proc_param, samples_struct, pos_param, dir, inv_param, pic_param,ps_param);

	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif
	cout << rank << " EXIT SUCCESS  !!!\n";
	return 0;
}




void print_struct(struct param_saneProc proc_param, struct samples samples_struct, struct param_sanePos pos_param, struct param_common dir,
		struct param_saneInv inv_param, struct param_sanePic pic_param, struct param_sanePS ps_param)
{

	cout << "\ndir : struct common\n";
	print_common(dir);

	cout.precision(std::numeric_limits< double >::digits10);

	cout << "\npos_param : struct param_sanePos\n";
	cout << "maskfile   = "  << pos_param.maskfile << endl;
	cout << "pixdeg     = " << pos_param.pixdeg << endl;
	cout << "flgdupl    = " << pos_param.flgdupl << endl;
	cout << "projgaps   = " << pos_param.projgaps << endl;
	cout << "fileFormat = " << pos_param.fileFormat << endl;
	cout << "lon        = " << pos_param.lon << endl;
	cout << "lat        = " << pos_param.lat << endl;

	cout << "\nproc_param : struct param_saneProc\n";
	cout << "CORRon           = " << proc_param.CORRon << endl;
	cout << "remove_polynomia = " << proc_param.remove_polynomia << endl;
	cout << "poly_order       = " << proc_param.poly_order << endl;
	cout << "napod            = " << proc_param.napod << endl;
	cout << "fsamp            = " << proc_param.fsamp << endl;
	cout << "fsamp_file       = " << proc_param.fsamp_file << endl;
	cout << "fhp              = " << proc_param.fhp << endl;
	cout << "fhp_file         = " << proc_param.fhp_file << endl;
	cout << "highpass_filter  = " << proc_param.highpass_filter << endl;
	cout << "fcut             = " << proc_param.fcut << endl;
	cout << "fcut_file        = " << proc_param.fcut_file << endl;

	cout << "\ninv_param : struct param_saneInv\n";
	cout << "cov_matrix        = " << inv_param.cov_matrix << endl;
	cout << "cov_matrix_suffix = " << inv_param.cov_matrix_suffix << endl;

	cout << "\npic_param : struct param_sanePic\n";
	cout << "iterW = " << pic_param.iterw << endl;

	cout << "\nps_param : struct param_sanePS\n";
	cout << "ell        = " << ps_param.ell << endl;
	cout << "ell_suffix = " << ps_param.ell_suffix << endl;
	cout << "mix        = " << ps_param.mix << endl;
	cout << "mix_suffix = " << ps_param.mix_suffix << endl;
	cout << "signame    = " << ps_param.signame << endl;
	cout << "ncomp      = " << ps_param.ncomp << endl;

	cout << "\nsamples_struct : struct samples\n";
	cout << "ntotscan = " << samples_struct.ntotscan << endl;


	cout << endl << "samples_struct :" << endl;
	cout << "framegiven = " << samples_struct.framegiven << endl;
	cout << endl;
	cout << "iframe_min = " << samples_struct.iframe_min << endl;
	cout << "iframe_max = " << samples_struct.iframe_max << endl;
	cout << endl;

	cout << endl << "scans_index :" << endl;
	for(long ii=0; ii< (long)samples_struct.scans_index.size(); ii++) {
		cout << ii << " : " ;
		cout << samples_struct.scans_index[ii]  << endl;
	}

	cout << endl << "fitsvect :" << endl;
	for(long ii=0; ii< (long)samples_struct.fitsvect.size(); ii++) {
		cout << ii << " : " ;
		cout << samples_struct.fitsvect[ii]  << endl;
	}

	cout << endl << "basevect :" << endl;
	for(long ii=0; ii< (long)samples_struct.basevect.size(); ii++) {
		cout << ii << " : " ;
		cout << samples_struct.basevect[ii]  << endl;
	}

	cout << endl << "fcut : fsamp : fhp : nsamples :" << endl;
	for(long ii=0; ii< (long)samples_struct.fsamp.size(); ii++) {
		cout << ii << " : " ;
		cout << samples_struct.fcut[ii] << " : ";
		cout << samples_struct.fsamp[ii] << " : ";
		cout << samples_struct.fhp[ii] << " : ";
		cout << samples_struct.nsamples[ii]  << endl;
	}

	cout << endl << "ell_names :" << endl;
	for(long ii=0; ii< (long)samples_struct.ell_names.size(); ii++) {
		cout << ii << " : ";
		cout << samples_struct.ell_names[ii] << endl;
	}

	cout << endl << "mix_names :" << endl;
	for(long ii=0; ii< (long)samples_struct.mix_names.size(); ii++) {
		cout << ii << " : ";
		cout << samples_struct.mix_names[ii] << endl;
	}

	cout << endl << "nbins :" << endl;
	for(long ii=0; ii< (long)samples_struct.nbins.size(); ii++) {
		cout << ii << " : ";
		cout << samples_struct.nbins[ii] << endl;
	}

	cout << endl << "ndet :" << endl;
	for(long ii=0; ii< (long)samples_struct.ndet.size(); ii++) {
		cout << ii << " : ";
		cout << samples_struct.ndet[ii] << endl;
	}
	cout << endl << "noisevect :" << endl;
	for(long ii=0; ii< (long)samples_struct.noisevect.size(); ii++) {
		cout << ii << " : ";
		cout << samples_struct.noisevect[ii] << " : ";
		cout << samples_struct.fcut[ii] << endl;
	}
	cout << endl << "bolovect :" << endl;
	for(long ii=0; ii< (long)samples_struct.bolovect.size(); ii++) {
		cout << ii << " : ";
		cout << samples_struct.bolovect[ii]  << endl;
	}

	cout << endl << "bolo_list :" << endl;
	for (long ii=0; ii< (long) samples_struct.bolo_list.size(); ii++){
		cout << ii << " [" << samples_struct.bolo_list[ii].size() << "] : ";
		for (long jj=0; jj<4; jj++)
			cout << samples_struct.bolo_list[ii][jj] << ", " ;
		cout << "... " << endl;
	}



}
