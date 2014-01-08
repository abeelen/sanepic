#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <list>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <ctime>

#include <fcntl.h>
#include <unistd.h>
#include <sysexits.h>

#include <limits>

#include "MPIConfiguration.h"
#include "StructDefinition.h"

#include "DataIO.h"
#include "ImageIO.h"
#include "TemporaryIO.h"
#include "InputFileIO.h"
#include "ParserFunctions.h"

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

	int      rank,      size; /* MPI processor rank and MPI total number of used processors */
	int  bolo_rank,  bolo_size; /* As for parallel scheme */
	int node_rank, node_size; /* On a node basis, same as *sub* but for frame scheme */

#ifdef USE_MPI

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	MPI_Comm MPI_COMM_NODE, MPI_COMM_MASTER_NODE;
#else
	size = 1;
	rank = 0;
	bolo_size  = 1;
	bolo_rank  = 0;
	node_size = 1;
	node_rank = 0;
#endif


	struct param_saneProc proc_param; /*! A structure that contains user options about preprocessing properties */
	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_sanePos Pos_param; /*! A structure that contains user options about map projection and properties */
	struct param_common dir; /*! structure that contains output input temp directories */

	string field; /*! actual boloname in the bolo loop */
	struct param_sanePic Pic_param;

	// those variables will not be used by sanePic but they are read in ini file (to check his conformity)
	struct param_sanePS PS_param;
	struct param_saneInv inv_param;
	string parser_output = "";

	int flagon = 0; /*  if at least one sample is rejected, flagon=1 */
	int factdupl = 1; /* map duplication factor */
	long long addnpix = 0; /* number of pix to add to compute the final maps in case of duplication + box constraint */
	long long npixsrc = 0; /* number of pix in box constraint */

	// map making parameters
	long long indpix_size; /* indpix read size */
	long long indpsrc_size; /* indpsrc read size */
	long long npix; /* nn = side of the map, npix = number of filled pixels */

	int nwcs=1;             // We will only deal with one wcs....
	struct wcsprm * wcs;    // wcs structure of the image
	long NAXIS1, NAXIS2;  // size of the image

	char * subheader;       // Additionnal header keywords
	int nsubkeys=0;           //

	long long *indpix, *indpsrc; /* pixels indices, mask pixels indices */




	// parser variables
	int parsed = 0;

	if (argc < 2) // no enough arguments
		parsed = 1;
	else {

		/* parse ini file and fill structures */
		parsed = parser_function(argv[1], parser_output, dir,
				samples_struct, Pos_param, proc_param, PS_param, inv_param,
				Pic_param, size, rank);


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


	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, Pos_param,  proc_param,
				PS_param, Pic_param, inv_param);

	}


	// this should be done by subrank 0 only
	//	init_dirfile(dir.tmp_dir, samples_struct, Pos_param.fileFormat, rank);
	//	fill_sanePS_struct(PS_param, samples_struct, dir);

	//	cout << rank << " parser done." << endl;


#ifdef USE_MPI

	MPI_Barrier(MPI_COMM_WORLD);

	if(configureMPI(dir.output_dir, samples_struct, rank, size,
			bolo_rank,  bolo_size, node_rank, node_size,
			MPI_COMM_NODE, MPI_COMM_MASTER_NODE)){
		if (rank==0)
			cerr << endl << endl << "Exiting... (not)" << endl;

	}
	MPI_Barrier(MPI_COMM_WORLD);

#endif


	if (bolo_rank == 0) {
		if ( readNoiseBinSizeFromDirfile(dir.tmp_dir, samples_struct) ) {
			cerr << "EE - Error in reading Noise sizes - Did you run saneInv ?" << endl;
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, EX_CONFIG);
#endif
			return EX_CONFIG;
		}
	}

#ifdef USE_MPI
	std::vector<long> ndet_tot, nbins_tot;
	ndet_tot.assign(samples_struct.ntotscan,0);
	nbins_tot.assign(samples_struct.ntotscan,0);

	MPI_Barrier(MPI_COMM_NODE);

	MPI_Allreduce(&(samples_struct.ndet[0]),   &ndet_tot[0], samples_struct.ntotscan, MPI_LONG, MPI_SUM, MPI_COMM_NODE);
	MPI_Allreduce(&(samples_struct.nbins[0]), &nbins_tot[0], samples_struct.ntotscan, MPI_LONG, MPI_SUM, MPI_COMM_NODE);

	samples_struct.ndet  =  ndet_tot;
	samples_struct.nbins = nbins_tot;
	ndet_tot.clear();
	nbins_tot.clear();

	//      MPI_Bcast_vector_long(samples_struct.ndet, 0, MPI_COMM_SUB);
	//      MPI_Bcast_vector_long(samples_struct.nbins, 0, MPI_COMM_SUB);
	MPI_Barrier(MPI_COMM_WORLD);

#endif

	// 	if (rank == 0)
	// 		readFramesFromFits(samples_struct);
	// #ifdef USE_MPI
	// 	MPI_Bcast_vector_long(samples_struct.nsamples, 0, MPI_COMM_WORLD);
	// #endif
#ifdef USE_MPI
	// Broadcast all position related values
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// Read file sizes once for all
	if (bolo_rank == 0) {
		if ( readFramesFromDirfile(dir.tmp_dir, samples_struct)) {
			cerr << "EE - Error in frame size - Did you run sanePre ?" << endl;
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, EX_CONFIG);
#endif
			return EX_CONFIG;
		}
	}
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_NODE);
	std::vector<long> nsamples_tot;
	nsamples_tot.assign(samples_struct.ntotscan, 0);

	MPI_Allreduce(&(samples_struct.nsamples[0]), &nsamples_tot[0], samples_struct.ntotscan, MPI_LONG, MPI_SUM, MPI_COMM_NODE);

	samples_struct.nsamples = nsamples_tot;
	nsamples_tot.clear();

	//      MPI_Bcast_vector_long(samples_struct.nsamples, 0, MPI_COMM_SUB);
	MPI_Barrier(MPI_COMM_WORLD);

#endif

#ifdef USE_MPI
	// Broadcast all position related values
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if (Pos_param.flgdupl)
		factdupl = 2; // default 0 : if flagged data are put in a duplicated map

	// Read indpsrc and checks...

	//      read pointing header
	if(read_keyrec(dir.tmp_dir, wcs, &NAXIS1, &NAXIS2, &subheader, &nsubkeys, rank)){ // read keyrec file
		cerr << "EE - Error reading saved keyrec -- Did you run sanePre ?" << endl;
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_IOERR;
	}

#ifdef USE_MPI
	// Broadcast all position related values
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if (rank == 0){

		if (readIndexCCR(indpsrc_size, npixsrc, indpsrc, dir.tmp_dir)) { // read mask index
			cerr << "EE - Please run sanePos" << endl;

#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return (EX_IOERR);
		}

		if (indpsrc_size != NAXIS1 * NAXIS2) { // check size compatibility
			if (rank == 0)
				cerr << "EE - indpsrc size is not the right size : Check indpsrc.bin file or run sanePos"
				<< endl;
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return (EX_IOERR);
		}

		addnpix = samples_struct.ntotscan * npixsrc;

		// read indpix
		if (read_indpix(indpix_size, npix, indpix, dir.tmp_dir, flagon)) { // read map index
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return (EX_IOERR);
		}

		if (indpix_size != (factdupl * NAXIS1 * NAXIS2 + 2 + addnpix + 1)) { // check size compatibility
			if (rank == 0)
				cout
				<< "indpix size is not the right size : Check Indpix_*.bi file or run sanePos"
				<< endl;
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return (EX_IOERR);
		}

	} // rank ==0


#ifdef USE_MPI
	// Broadcast all position related values
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&npix,         1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&npixsrc,      1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&addnpix,      1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&indpix_size,  1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&indpsrc_size, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank!=0){
		indpix  = new long long[indpix_size];
		indpsrc = new long long[indpsrc_size];
	}

	MPI_Bcast(indpix,  indpix_size,  MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(indpsrc, indpsrc_size, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
#endif


#ifdef USE_MPI
	// Broadcast all position related values
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if (rank == 0) {
		computeMemoryRequirement_sanePic(samples_struct, npix, indpix_size);

		for (int ii=0; ii< samples_struct.ntotscan; ii++)
			cout << ii << " " << samples_struct.memory[ii]*8. << endl;
	//	cout << ii << " " << prettyPrintSize(samples_struct.memory[ii]*8.) << endl;
	}

	//	get_noise_bin_sizes(dir.tmp_dir, samples_struct, rank);

#ifdef USE_MPI
	// Broadcast all position related values
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if (rank==0){
		//		print_common(dir);
		//		print_param_sanePos(Pos_param);
		//		print_param_saneProc(proc_param);
		//		print_param_saneInv(inv_param);
		//		print_param_sanePic(Pic_param);
		//		print_param_sanePS(PS_param);

		cout << " ---" << endl;
		print_struct(proc_param, samples_struct, Pos_param, dir, inv_param, Pic_param,PS_param);

	}

#ifdef USE_MPI
	// Broadcast all position related values
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if (rank == 0)
		cout << endl << "iframe_min/max per node :" << endl;

	if (bolo_rank == 0){
		ostringstream line;
		line << rank << " " << bolo_rank << " : " << "iframe from " << samples_struct.iframe_min << " to " << samples_struct.iframe_max;
		cout << line.str() << endl;
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif
	// cout << rank << " EXIT SUCCESS  !!!\n";
	return 0;
}




void print_struct(struct param_saneProc proc_param, struct samples samples_struct, struct param_sanePos Pos_param, struct param_common dir,
		struct param_saneInv inv_param, struct param_sanePic Pic_param, struct param_sanePS PS_param)
{

	cout << "\ndir : struct common\n";
	print_common(dir);

	cout.precision(std::numeric_limits< double >::digits10);

	cout << "\nPos_param : struct param_sanePos\n";
	cout << "maskfile   = "  << Pos_param.maskfile << endl;
	cout << "pixdeg     = " << Pos_param.pixdeg << endl;
	cout << "flgdupl    = " << Pos_param.flgdupl << endl;
	cout << "projgaps   = " << Pos_param.projgaps << endl;
	cout << "fileFormat = " << Pos_param.fileFormat << endl;
	cout << "lon        = " << Pos_param.lon << endl;
	cout << "lat        = " << Pos_param.lat << endl;

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

	cout << "\nPic_param : struct param_sanePic\n";
	cout << "iterW = " << Pic_param.iterw << endl;

	cout << "\nPS_param : struct param_sanePS\n";
	cout << "ell        = " << PS_param.ell << endl;
	cout << "ell_suffix = " << PS_param.ell_suffix << endl;
	cout << "mix        = " << PS_param.mix << endl;
	cout << "mix_suffix = " << PS_param.mix_suffix << endl;
	cout << "signame    = " << PS_param.signame << endl;
	cout << "ncomp      = " << PS_param.ncomp << endl;

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
