
#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parse_saneCheck.h"
#include "tools.h"
#include "struct_definition.h"

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()


extern "C" {
#include "nrutil.h"
#include "nrcode.h"
}

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;


int main(int argc, char *argv[]) {


	int rank, size;

#ifdef USE_MPI

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#else
	size = 1;
	rank = 0;
	cout << "Mpi is not used for this step" << endl;
#endif

	/*
	 *  This should be organized as :
	 *
	 *  - parse the input ini file (OK)
	 *  - check for existence of directory/files pointed from the ini file (OK)
	 *  - for each file, check for:
	 *      - presence of the different hdus (only the position change between the two format) (OK)
	 *      - consistent sizes of the different hdu (OK)
	 *      - presence of non flagged NaN values (OK)
	 *
	 *      - bad channels [latter, noisier channels] (OK)
	 *
	 *      - check the fsamp from the ini file, should be the most frequent one (OK)
	 *      - check for time gaps > 1.9 * the most frequent time gap (OK)
	 */


	int parsed=1;

	struct samples samples_struct;
	struct common dir;
	struct detectors det;
	double fsamp;

	string outname;

	if(rank==0)
		printf("\nBeginning of saneCheck:\n\n");


	parsed=parse_saneCheck_ini_file(argv[1],dir,
			det, samples_struct, fsamp, rank);

	if(parsed==-1){
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		exit(1);
	}

	readFrames(samples_struct.fitsvect, samples_struct.nsamples);

	struct detectors bolo_fits_0;
	read_bolo_list(samples_struct.fitsvect[0],bolo_fits_0);

	long *bolo_bad_tot;
	long *bolo_bad_80_tot;

	if(rank==0){ // only for MPI_reduce
		bolo_bad_tot= new long [bolo_fits_0.ndet];
		bolo_bad_80_tot= new long [bolo_fits_0.ndet];
		fill(bolo_bad_tot, bolo_bad_tot + bolo_fits_0.ndet ,0);
		fill(bolo_bad_80_tot, bolo_bad_80_tot + bolo_fits_0.ndet ,0);
	}


	for(int ii=0;ii<samples_struct.ntotscan;ii++){

		int do_it=who_do_it(size, rank, ii);

		if(rank==do_it){

			struct detectors bolo_fits;
			long *bolo_bad;
			long *bolo_bad_80;

			int format_fits=0;

			cout << endl << endl << "[" << rank <<  "] Checking : " << samples_struct.fitsvect[ii] << endl << endl;

			format_fits=test_format(samples_struct.fitsvect[ii]); // format = 1 => HIPE, else Sanepic

			read_bolo_list(samples_struct.fitsvect[ii],bolo_fits);
			if(check_bolos(bolo_fits.boloname, bolo_fits_0.boloname)){
				cout << "Skipping file : " << samples_struct.fitsvect[ii] << ". Please run only together scans that correspond to the same field\n";
				continue;
			}

			check_detector_is_in_fits(det,bolo_fits,samples_struct.fitsvect[ii]);

			bolo_bad = new long[bolo_fits.ndet];
			bolo_bad_80 = new long[bolo_fits.ndet];
			fill(bolo_bad,bolo_bad+bolo_fits.ndet,0);
			fill(bolo_bad_80,bolo_bad_80+bolo_fits.ndet,0);

			cout << "\n[" << rank <<  "] Checking presence of common HDU and position HDU\n";
			check_commonHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits);
			check_positionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits, format_fits);

			if(format_fits==1){
				cout << "[" << rank <<  "] HIPE format found, Checking Alt position HDU presence\n";
				check_altpositionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits);
			}

			cout << "\n[" << rank <<  "] Checking NANs in common HDU and position HDU\n";
			check_NAN_commonHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits);
			check_NAN_positionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits);

			if(format_fits==1){
				cout << "[" << rank <<  "] HIPE format found, Checking NANs in Alt position HDU\n";
				check_NAN_altpositionHDU(samples_struct.fitsvect[ii],samples_struct.nsamples[ii],bolo_fits);
			}

			cout << "\n[" << rank <<  "] Checking time gaps in time table\n";
			check_time_gaps(samples_struct.fitsvect[ii],samples_struct.nsamples[ii], fsamp, dir);

			check_flag(samples_struct.fitsvect[ii],bolo_fits, samples_struct.nsamples[ii],outname, bolo_bad,bolo_bad_80);

			// reduce
#ifdef USE_MPI
			MPI_Reduce(bolo_bad,bolo_bad_tot,bolo_fits_0.ndet,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Reduce(bolo_bad_80,bolo_bad_80_tot,bolo_fits_0.ndet,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
#else
			for(long kk = 0; kk< bolo_fits_0.ndet; kk++){
				bolo_bad_tot[kk]+=bolo_bad[kk];
				bolo_bad_80_tot[kk]+=bolo_bad_80[kk];
			}
#endif
			delete [] bolo_bad;
			delete [] bolo_bad_80;
		}
	}

	if(rank==0){
		cout << endl;

		// generating log files :
		outname = dir.output_dir + "bolo_totally_flagged.txt";
		cout << "Writing informations in :\n" << outname << endl << endl;
		log_gen(bolo_bad_tot,outname, bolo_fits_0);


		outname = dir.output_dir + "bolo_80_percent_flagged.txt";
		cout << "Writing informations in :\n" << outname << endl;
		log_gen(bolo_bad_80_tot, outname, bolo_fits_0);


		cout << "\nEnd of saneCheck\n";
	}

	//cleaning
	delete [] samples_struct.nsamples;


#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

}
