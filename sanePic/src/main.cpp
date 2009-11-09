/*
 * sanepic_Main_Loop.cpp
 *
 *  Created on: 29 mai 2009
 *      Author: matthieu
 */




#include "parseSanepic.h"
#include "conjugate_gradient.h"
#include "mpi_architecture_builder.h"

#include <iostream>
#include <vector>

#include "imageIO.h"
#include "inline_IO2.h"





#ifdef USE_MPI
#include "mpi.h"
#include <algorithm>
#include <fstream>
#endif


using namespace std;



//**********************************************************************************//
//**********************************************************************************//
//***************** Beginning of conjugate gradient program ************************//
//**********************************************************************************//
//**********************************************************************************//




int main(int argc, char *argv[])
{

	//

	int size;//, size_det;
	int rank;//, rank_det;
#ifdef USE_MPI
	// int tag = 10;
	//MPI_Status status;

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#else
	size = 1;
	rank = 0;
	cout << "Mpi will not be used for the main loop" << endl;
#endif



	//************************************************************************//
	//************************************************************************//
	//main conjugate gradient loop
	//************************************************************************//
	//************************************************************************//

	struct user_options u_opt;
	struct samples samples_struct;
	struct input_commons com;
	struct directories dir;
	struct detectors det;



	//DEFAULT PARAMETERS
	com.napod = 0; /*!  number of samples to apodize */
	u_opt.fsamp = 0.0; //25.0; /*!  sampling frequency : BLAST Specific */
	u_opt.projgaps = 0; /*!1: data flagged are put in a single pixel  (assume no signal in this pixel),
	0: data flagged are not reprojected */

	long iframe_min, iframe_max; /*! For mpi usage : defines min/max number of frame for each processor */
	int flagon = 0; /*!  if one sample is rejected, flagon=1 */
	int iterw = 10; /*!  period in iterations to which the data are written to disk, 0 = no intermediate map to be written*/
	u_opt.NORMLIN = 0; /*!  baseline is removed from the data, NORMLIN = 1 else 0 */
	com.NOFILLGAP = 0; /*!  fill the gap ? default is YES (debug parameter) */
	com.flgdupl = 0; /*!  1 if flagged data are put in a separate map */
	u_opt.remove_polynomia = 1; /*! Remove a fitted polynomia from the data ? */
	u_opt.CORRon = 1; /*!  correlation included in the analysis (=1), else 0, default 0 */
	int factdupl = 1; /*! map duplication factor */
	long long addnpix=0; /*! number of pix to add to compute the final maps in case of duplication + box constraint */
	long long npixsrc = 0; /*! number of pix in box constraint */



	samples_struct.ntotscan=0; /*! total number of scans */
	det.ndet=0; /*! number of channels */


	// map making parameters
	com.pixdeg=-1.0; /*! size of pixels (degree) */

	long long npix2; /*! used to check PNd reading was correct */
	long long ind_size; /*! indpix read size */
	long NAXIS1, NAXIS2;
	long long npix; /*! nn = side of the map, npix = number of filled pixels */


	//internal data params
	u_opt.f_lp=0.0; /*! frequencies : filter knee freq, noise PS threshold freq ; frequencies converted in a number of samples */




	double *PNdtot; /*! to deal with mpi parallelization : Projected noised data */
	long long *indpix, *indpsrc; /*! pixels indices, mask pixels indices */



	string field; /*! actual boloname in the bolo loop */
	string prefixe; /*! prefix used for temporary name file creation */

	std::vector<double> fcut; /*! noise cutting frequency vector */
	//std::vector<string> extentnoiseSP; /*! noise filenames vector of string */
	std::vector<struct box> boxFile;



	// main loop variables
	double *S; /*! Pure signal */

	// parallel scheme file
	string fname; /*! parallel scheme filename */



	// Parse ini file
	if (argc<2) {
		printf("Please run %s using a *.ini file\n",argv[0]);
		exit(0);
	} else {

		int parsed=1;

		parsed=parse_sanePic_ini_file(argv[1],u_opt,iterw, dir, samples_struct,com,
				det,boxFile, fcut, rank);
		if (parsed==-1){
#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
#endif
			exit(1);
		}
		//exit(0);

	}
	////////////////////////////////////////////////////////////////


	//frames_index = new long [samples_struct.ntotscan];
	//extentnoiseSp_all = new string[samples_struct.ntotscan];


	samples_struct.fits_table = new string[samples_struct.ntotscan];
	samples_struct.index_table= new long[samples_struct.ntotscan];
	samples_struct.noise_table = new string[samples_struct.ntotscan];


	//vector2array(extentnoiseSP,  extentnoiseSp_all);
	cout << samples_struct.nsamples[0] << endl;


	/********************* Define parallelization scheme   *******/

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
		vector2array(samples_struct.scans_index,  samples_struct.index_table); // TODO : passer index_table en int plutot que long

		if(rank==0){
			//check the processor order given is correct
			//			size_tmp = *max_element(samples_struct.index_table, samples_struct.index_table+samples_struct.ntotscan);

			struct sortclass_long sortobject;
			sort(samples_struct.scans_index.begin(), samples_struct.scans_index.end(), sortobject);

			std::vector<long>::iterator it;
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


		for(long ii = 0; ii<size; ii++){
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

#endif

	//
	//	// allocate memory
	//	tancoord = new double[2];
	//	tanpix = new double[2];

	//	read_info_pointing(NAXIS1, NAXIS2, u_opt.outdir, tanpix, tancoord);
	struct wcsprm * wcs;
	read_MapHeader(dir.tmp_dir,wcs, &NAXIS1, &NAXIS2);

	// read nn, coordsyst, tanpix, tancoord
	//	read_info_pointing(NAXIS1, NAXIS2, u_opt.tmp_dir, tanpix, tancoord);
	//cout << tanpix[0] << " " << tanpix[1] << endl;
	//cout << tancoord[0] << " " << tancoord[1] << endl;

	cout << "Map size :" << NAXIS1 << "x" << NAXIS2 << endl;


	// TODO : The mask should not be remade everytime...
	//******************************** some preprocess again  ****************/

	unsigned short *mask;
	mask    = new unsigned short[NAXIS1*NAXIS2];
	indpsrc = new long long[NAXIS1*NAXIS2];

	// Initialize the masks
	addnpix=0;
	npixsrc=0;
	for (long long ii=0; ii<NAXIS1*NAXIS2; ii++){
		mask[ii]    =  1;
		indpsrc[ii] = -1;
	}

	// TODO : untested....
	// if a box for crossing constraint removal is given in ini file
	// TODO : save mask in fits file
	// TODO : being able to read a mask in fits file format
	for (unsigned long iBox = 0; iBox < boxFile.size(); iBox++){
		for (long ii=boxFile[iBox].blc.x; ii<boxFile[iBox].trc.x ; ii++)
			for (long jj=boxFile[iBox].blc.y; jj<boxFile[iBox].trc.y; jj++){
				mask[jj*NAXIS1 + ii] = 0;
				indpsrc[jj*NAXIS1 + ii] = npixsrc++;
			}
	}

	// each frame contains npixsrc pixels with index indsprc[] for which
	// crossing constraint are removed
	// thus
	// addnpix = number of pix to add in pixon
	//         = number of scans * number of pix in box crossing constraint removal
	addnpix = samples_struct.ntotscan*npixsrc;

	// compute indpsrc and addnpix
	//sanepic_preprocess(NAXIS1,NAXIS2, xxi, xxf, yyi, yyf, indpsrc, npixsrc, ntotscan, addnpix);


	if (com.flgdupl) factdupl = 2; // -M =1, default 0 : if flagged data are put in a duplicated map

	// read npix, PNdtot from file
	read_PNd(PNdtot, npix,  dir.tmp_dir);
//	for (int ii=0;ii<20;ii++)
//			cout << PNdtot[ii] << " ";
//		cout << endl << "avant read indpix\n";
//		getchar();



	// read indpix
	read_indpix(ind_size, npix2, indpix,  dir.tmp_dir, flagon);

	if(ind_size!=(factdupl*NAXIS1*NAXIS2+2 + addnpix)){
		cout << "indpix size is not the right size : Check Indpix_*.bi file or run sanePos" << endl;
		exit(0);
	}
	if (npix!=npix2){
		cout << "Warning ! Indpix_for_conj_grad.bi and PNdCorr_*.bi are not compatible, npix!=npix2" << endl;
		exit(0);
	}


	/*************************************************************/

	printf("[%2.2i] iframe_min %ld\tiframe_max %ld \n",rank,iframe_min,iframe_max);

	if (iframe_min < 0 || iframe_min >= iframe_max || iframe_max > samples_struct.ntotscan){
		cerr << "Error distributing frame ranges. Check iframe_min and iframe_max. Exiting" << endl;
		exit(1);
	}

	/* END PARAMETER PROCESSING */




	printf("[%2.2i] Main Conjugate gradient loop\n",rank);

	//MALLOC
	S = new double[npix];
	fill(S,S+npix,0.0);

	// conjugate GRADIENT LOOP
	//	sanepic_conjugate_gradient(com.flgdupl, npix, S, iframe_min, iframe_max,
	//			samples_struct.nsamples, fcut,u_opt.f_lp, u_opt.fsamp,
	//			indpix,
	//			wcs, NAXIS1, NAXIS2,
	//			factdupl, dir.tmp_dir, det.ndet,
	//			samples_struct.noise_table,dir.tmp_dir, det.boloname, iterw,
	//			indpsrc, npixsrc,flagon, u_opt.projgaps, rank, u_opt.CORRon,
	//			dir.dirfile, PNdtot, samples_struct.ntotscan,addnpix,u_opt.NORMLIN,com.NOFILLGAP,
	//			com.napod, u_opt.remove_polynomia, dir.outdir,samples_struct.fits_table);

	sanepic_conjugate_gradient(samples_struct,com,det,dir,u_opt, npix, S, iframe_min, iframe_max,
			fcut,indpix,wcs, NAXIS1, NAXIS2, iterw,
			indpsrc, npixsrc,flagon, rank, PNdtot, addnpix);



	//******************************************************************//
	//******************************************************************//
	//*********************** End of program ***************************//
	//******************************************************************//
	//******************************************************************//




	// TODO : This will be rewrite differently
	//	if (rank == 0){
	//		//write infos for second part
	//		write_info_for_second_part(u_opt.outdir, NAXIS1, NAXIS2, npix,u_opt.pixdeg, tancoord, tanpix, coordsyst, flagon, indpix);
	//	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif


	// clean up
	delete [] S;

	delete [] samples_struct.nsamples;
	delete [] samples_struct.fits_table;
	delete [] samples_struct.index_table;
	delete [] samples_struct.noise_table;


	delete [] indpsrc;
	delete [] indpix;
	delete [] PNdtot;

	//delete [] frames_index;





#ifdef USE_MPI
	MPI_Finalize();
#endif



	return 0;
}



