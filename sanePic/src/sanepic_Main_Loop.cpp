/*
 * sanepic_Main_Loop.cpp
 *
 *  Created on: 29 mai 2009
 *      Author: matthieu
 */





// liste des variables a donner :


// relancer le mpi et faire un best frame order

/*includes*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include "todprocess.h"
#include "map_making.h"
#include "sane_io.h"
#include "estimPS.h"
#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"
#include "mpi_architecture_builder.h"
#include <time.h>
#include <fftw3.h>
#include <fcntl.h>
#include <unistd.h>
#include <list>
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;


void usage(char *name)
{
	cerr << "USAGE: " << name << ": simulate time stream data" << endl << endl;
	cerr << name << " [-F <data dir>] [-f <first frame>]" << endl;
	cerr << "\t[-l <last frame> | -n <num frames>]" << endl;
	cerr << "\t[-y <channel>] [-C <channel file>] [-p <pixsize>]" << endl;
	cerr << "-H <filter freq>       frequency of the high pass filter applied to the data" << endl;
	cerr << "-J <noise freq cut>    frequency under which noise power spectra are thresholded. Default is the frequency cut of the high pass filter applied to the data (-H option)" << endl; // ok
	cerr << "\t[-O <offset file>] [-B <bolo ext>] [-G <flag ext>] [-P <pointing ext>]" << endl;
	cerr << "\t[-o <output dir>] [-e <output file>]" << endl << endl;
	cerr << "-A <apodize nsamp>     number of samples to apodize (supercedes -a)\n"; // ok
	cerr << "-B <bolo ext>          bolo field extension" << endl;
	cerr << "-C <channel file>      include bolometers named in file" << endl;
	cerr << "-e <output file str>   id string in output file" << endl;
	cerr << "-F <data dir>          source data directory" << endl;
	cerr << "-f <first frame>       first frame (multiple allowed)" << endl;
	cerr << "-G <flag ext>          flag field extension" << endl;
	cerr << "-H <filter freq>       filter frequency (multiple allowed)" << endl;
	cerr << "-R <sampling frequency> detectors sampling frequency. Required" << endl;
	cerr << "-K <noise prefixe>     noise power spectrum file prefixe " << endl;
	cerr << "-k <noise suffixe>     noise power spectrum file suffixe " << endl;
	cerr << "-l <last frame>        last frame (multiple allowed, exclusive from -n)" << endl;
	cerr << "-n <num frames>        number of frames (multiple allowed, exclusive from -l)" << endl;
	cerr << "-x <box coord x1>      if set, defines the x coordinate of the left side of the box(es) in which crossing constraints are removed" << endl;
	cerr << "-X <box coord x2>      if set, defines the x coordinate of the right side of the box(es) in which crossings constraints are removed" << endl;
	cerr << "-z <box coord y1>      if set, defines the y coordinate of the top side of the box(es) in which crossing constraints are removed" << endl;
	cerr << "-Z <box coord y2>      if set, defines the y coordinate of the bottom side of the box(es) in which crossings constraints are removed" << endl;
	//cerr << "-O <offset file>       file containing bolometer offsets" << endl;
	cerr << "-o <output dir>        output directory" << endl;
	cerr << "-P <pointing ext>      pointing field extension" << endl;
	cerr << "-p <pixsize>           size of pixels (degrees)" << endl;
	cerr << "-y <channel>           include bolometer (multiple allowed)" << endl;
	//cerr << "-S <Scan offset file>  file containing offsets for different frame ranges" << endl;
	cerr << "-t <RA source>         RA (or Gal l if keyword c is 2) of the tangent point and of the source for telescope coordinate maps" << endl;
	cerr << "-T <DEC source>        DEC (or Gal b if keyword c is 2) of the tangent point and of the source for telescope coordinate maps" << endl;
	cerr << "-u <RA min>            fixed minimum RA (or Gal l) of the map, if none of -N, -t, and -T are defined" << endl;
	cerr << "-U <RA max>            fixed maximum RA (or Gal l) of the map, if none of -N, -t, and -T are defined" << endl;
	cerr << "-v <DEC min>           fixed minimum DEC (or Gal b) of the map, if none of -N, -t, and -T are defined" << endl;
	cerr << "-V <DEC max>           fixed maximum DEC (or Gal b) of the map, if none of -N, -t, and -T are defined" << endl;
	cerr << "-c <coord system>      Coordinate system, 1: RA/DEC, 2: L/B (Gal), 3: Telescope coordinates, Default is RA/DEC" << endl;
	//cerr << "-N <map radius>        fixed radius (half a side) of the map in degrees (in case of data outside, they are constrained to have constant values)" << endl;
	cerr << "-L <no baseline>       Keyword specifiying if a baseline is removed from the data or not (0 for YES, 1 for NO)" << endl;
	cerr << "-r <precompute PNd>    Keyword specifying if PNd is precomputed and read from disk" << endl;
	cerr << "-g <project gaps>      Keyword specifying if gaps are projected to a pixel in the map, if so gap filling of noise only is performed iteratively. Default is 0" << endl;
	cerr << "-M <map flagged data>   Keyword specifying if flagged data are put in a separate map, default is 0" << endl;
	cerr << "-s <time offset>       Keyword for subtracting a time offset to the data to match the pointing "<< endl;
	cerr << "-E <no corr>           Set this keyword to 0 if correlations are not included in the analysis" << endl;
	//cerr << "-I <parallel scheme>   Set this keyword to 1 in order to distribute different field visit (or frame ranges) to different processors, instead of different detector data chuncks as it is by default" << endl;
	cerr << "-j <write unconverged maps>  The value specified by this keyword indicates the period in iterations to which the data are written to disk. Default is 10. Set this keyword to zero if you don't want intermediate maps to be written." << endl;
	cerr << "-a <noise estim>       Optional. Enter filename containing the mixing matrix of noise components. If set, the noise power spectra files for each field are computed after map-making. Default is no re-estimation" << endl; // ok
	cerr << "-i <Noise PS etimation> Optional. Computes power spectra estimation from the data and final map and saves it in a file. Default is no estimation." << endl;
	exit(1);
}



template<class T> void list2array(list<T> l, T* a)
{
	// copy list of type T to array of type T
	typename list<T>::iterator iter;
	int i;

	for (iter=l.begin(), i=0; iter != l.end(); iter++, i++) {
		a[i] = *iter;
	}
}





//**********************************************************************************//
//**********************************************************************************//
//***************** Beginning of conjugate gradient program ************************//
//**********************************************************************************//
//**********************************************************************************//




int main(int argc, char *argv[])
{

	// read framesorder

	int size, size_det;
	int rank, rank_det;
#ifdef USE_MPI
	// int tag = 10;
	MPI_Status status;

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	cout << size << endl;
	cout << rank << endl;

#else
	size = 1;
	rank = 0;
	cout << "Mpi will not be used for the main loop" << endl;
#endif

	char * pPath;
	pPath = getenv ("TMPBATCH");
	if (pPath!=NULL)
		printf ("The current path is: %s",pPath);





	/*Variables*/
	time_t t1;

	//************************************************************************//
	//************************************************************************//
	//main conjugate gradient loop
	//************************************************************************//
	//************************************************************************//


	bool projgaps = 0; //1: data flagged are put in a single pixel
	//   (assume no signal in this pixel),
	//0: data flagged are not reprojected

	//default value of the data to pointing shift
	int shift_data_to_point = 0;

	int samples_per_frames = 20;

	//DEFAULT PARAMETERS
	long napod = 0; // number of samples to apodize
	double fsamp = 0.0; //25.0; // sampling frequency : BLAST Specific
	double errarcsec = 15.0; // rejection criteria : scerr[ii] > errarcsec, sample is rejected
	// source error

	long ii, jj, ll, ib, iframe; // loop indices
	long iframe_min, iframe_max;
	int flagon = 0; // if rejectsample [ii]==3, flagon=1
	int iterw = 10; // period in iterations to which the data are written to disk, 0 = no intermediate map to be written
	bool NORMLIN = 0; // baseline is removed from the data, NORMLIN = 1 else 0
	bool NOFILLGAP = 0; // fill the gap ? default is YES (debug parameter)
	//bool PND_ready = 0; // PNd precomputed ? read on disk if =1
	bool flgdupl = 0; // 1 if flagged data are put in a separate map
	bool CORRon = 1; // correlation included in the analysis (=1), else 0, default 0
	//bool parallel_frames = 0; // a parallel scheme is used : mpi has been launched
	int factdupl = 1;
	long addnpix=0;
	long npixsrc = 0;
	//bool bfixc = 0; // indicates that 4 corners are given for the cross corelation removal box

	//set coordinate system
	//double *srccoord, *coordscorner;
	//srccoord = new double[2]; // RA/DEC tangent point/source
	//coordscorner = new double[4]; // map min/max RA/DEC coords (-N,-t,-T absents)
	//srccoord[0] = -1000; // RA tangent point/source
	//srccoord[1] = -1000; // DEC tangent point/source
	//double radius = -1.0; // map radius (half a side) in degrees


	// data parameters
	long *fframes  ; // first frames table ff_in list -> fframes
	long *fframesorder ;
	long *nsamples ; // number of samples table nf_in list -> nsamples
	long *nsamplesorder ;
	long *ruleorder ;
	long *frnum ;

	// box for crossing constraints removal
	long *xxi_boxes, *xxf_boxes; // left x, right x
	long *yyi_boxes, *yyf_boxes; // top y, bottom y

	long ntotscan; // total number of scans
	long ndet; // number of channels
	int nnf; // extentnoiseSp_list number of elements


	// map making parameters
	double pixdeg; // size of pixels (degree)


	int nn, npix; // nn = side of the map, npix = number of filled pixels
	double *tancoord; // tangent point coordinates
	double *tanpix; // tangent pixel

	//internal data params
	long ns, ff; // number of samples for this scan, first frame number of this scan
	double f_lp, f_lp_Nk, f_lppix, f_lppix_Nk; // frequencies : filter knee freq, noise PS threshold freq ; frequencies converted in a number of samples


	FILE *fp;

	char testfile[100];



	//char type='d'; // returned type of read_data functions, d=64bit double
	unsigned char *mask; // samples flags, pointing flags, rejected samples list
	// mask = box for crossing constraint removal mask in the map
	double *PNd, *PNdtot; //
	long *indpix, *indpsrc; // pixels indices, mask pixels indices


	//long *pixon; // this array is used to store the rules for pixels : they are seen or not
	//long *samptopix; // sample to pixel conversion array




	string field; // actual boloname in the bolo loop
	string *bolonames; // channel list -> bolonames array : considered bolometers names
	string *extentnoiseSp_all; // ((list -> string*))
	//string *extentnoiseSp_allorder;
	string bolofield; // bolofield = boloname + bextension
	//string calfield; // calfield  = field+cextension;
	//string flagfield; // flagfield = field+fextension;
	string dirfile; // data directory
	string outdir; // output directory
	string poutdir; // current path (pPath) or output dir (outdir)
	string bextension; // bolometer field extension
	//string cextension = "NOCALP"; // needed for calfield calculation !
	string fextension = "NOFLAG"; // flag field extension
	string pextension; // pointing extension
	//string file_offsets; // bolometer offsets file
	//string file_frame_offsets = "NOOFFS"; // offset file
	string termin; // output file suffix
	string noiseSppreffile; // noise file suffix
	string extentnoiseSp; // noise file
	string prefixe; // prefix used for temporary name file creation

	string MixMatfile = "NOFILE";
	int doInitPS = 0;

	/* DEFAULT PARAMETERS */
	int coordsyst = 1; /// Default is RA/DEC
	int coordsyst2 = -1;



	/* COMMAND-LINE PARAMETER PROCESSING */


	int retval; // parser variable
	int tmpcount = 0; // parser variable (check if parsing was performed well)
	int tmpcount2 = 0; // parser variable (check if parsing was performed well)
	int temp; // parser variable (check if parsing was performed well)

	list<string> channel; // bolometer list
	list<long> ff_in, nf_in, xxi_in, xxf_in, yyi_in, yyf_in; // first frame list, number of frames per sample, box for crossing constraints removal coordinates lists (left x, right x, top y, bottom y)
	list<double> fcut_in;
	list<string> extentnoiseSp_list; // noise file prefixe


	//time t2, t3, t4, t5, dt;



	f_lp = 0.0; // low pass filter frequency
	f_lp_Nk = 0.0; // noise PS frequency threshold
	pixdeg = -1.0; // "Size of pixels (deg)"

	// main loop variables
	int npixeff, iter, idupl;
	double var0, var_n, delta0, delta_n, delta_o, rtq, alpha, beta;

	double *S;
	double *PtNPmatS,  *PtNPmatStot, *r, *q, *qtot, *d, *Mp, *Mptot, *s;
	long *hits, *hitstot;




	long mi;
	double *map1d;



	string fname;
	char iterchar[30];
	char iframechar[30];
	string iterstr, iframestr;


	// Parse command line options
	while ( (retval = getopt(argc, argv, "F:f:l:n:y:C:H:J:o:O:B:R:G:P:S:e:p:A:k:K:t:T:u:U:v:V:c:N:L:g:r:M:x:X:z:Z:s:E:I:j:a:D:i:")) != -1) {
		switch (retval) {
		case 'F':
			dirfile = optarg;
			break;
		case 'f':
			if (ff_in.size() != nf_in.size()) {
				cerr << "'-f' must be followed by '-l' or '-n'. Exiting.\n";
				exit(1);
			}
			ff_in.push_back(atoi(optarg));
			break;
		case 'l':
		case 'n':
			if (nf_in.size() != ff_in.size()-1) {
				cerr << "'-l' or '-n' must follow '-f'. Exiting.\n";
				exit(1);
			}
			temp = atoi(optarg);
			if (retval == 'l') temp -= ff_in.back() - 1;
			nf_in.push_back(temp);
			break;
		case 'x':
			if (xxi_in.size() != yyf_in.size()){
				cerr << "'-x' must be followed by '-X -z -Z'. Exiting.\n";
				exit(1);
			}
			xxi_in.push_back(atoi(optarg));
			break;
		case 'X':
			if (xxf_in.size() != xxi_in.size()-1){
				cerr << "'-X' must follow '-x' and be followed by '-z -Z'. Exiting. \n";
				exit(1);
			}
			xxf_in.push_back(atoi(optarg));
			break;
		case 'z':
			if (yyi_in.size() != xxi_in.size()-1){
				cerr << "'-z' must follow by '-x -X' and be followed by 'Z'. Exiting.\n";
				exit(1);
			}
			yyi_in.push_back(atoi(optarg));
			break;
		case 'Z':
			if (yyf_in.size() != xxi_in.size()-1){
				cerr << "'-Z' must follow '-x -X -z'. Exiting. \n";
				exit(1);
			}
			yyf_in.push_back(atoi(optarg));
			break;
		case 'p':
			pixdeg = atof(optarg);
			break;
		case 'H':
			f_lp = atof(optarg);
			break;
		case 'J':
			f_lp_Nk = atof(optarg);
			break;
		case 'A':
			napod = atoi(optarg);
			break;
		case 'y':
			channel.push_back(optarg);
			break;
		case 'C':
			read_bolofile(optarg, channel);
			//      cerr << "num ch: "<< channel.size() << endl;
			break;
		case 'o':
			outdir = optarg;
			break;
		case 'B':
			bextension = optarg;
			break;
		case 'R':
			fsamp = atof(optarg);
			//cextension = optarg;
			break;
		case 'G':
			fextension = optarg;
			break;
		case 'P':
			pextension = optarg;
			break;
		case 'O':// plus utile
			//file_offsets = optarg;
			break;
		case 'e':
			termin = optarg;
			break;
		case 'k':
			noiseSppreffile = optarg;
			break;
		case 'K':
			extentnoiseSp_list.push_back(optarg);
			break;
		case 'S':// plus utile
			//file_frame_offsets = optarg;
			break;
		case 't':// plus utile
			//srccoord[0] = atof(optarg);
			//tmpcount2 += 1;
			break;
		case 'T':// plus utile
			//srccoord[1] = atof(optarg);
			//tmpcount2 += 1;
			break;
		case 'u':// plus utile
			//coordscorner[0] = atof(optarg);
			//tmpcount += 1;
			break;
		case 'U':// plus utile
			//coordscorner[1] = atof(optarg);
			//tmpcount += 1;
			break;
		case 'v':// plus utile
			//coordscorner[2] = atof(optarg);
			//tmpcount += 1;
			break;
		case 'V':// plus utile
			//coordscorner[3] = atof(optarg);
			//tmpcount += 1;
			break;
		case 'c':
			coordsyst = atoi(optarg);
			break;
		case 'N':// plus utile
			//radius = atof(optarg);
			//tmpcount2 += 1;
			break;
		case 'L':
			NORMLIN = 1;
			break;
		case 'g':
			projgaps = 1;
			break;
		case 'r': // plus utile
			//PND_ready = atoi(optarg);
			break;
		case 'M':
			flgdupl = 1;
			break;
		case 's':
			shift_data_to_point = atoi(optarg);
			break;
		case 'E':
			CORRon = 1;
			break;
		case 'I': //useless
			//parallel_frames = 1;
			break;
		case 'j':
			iterw = atoi(optarg);
			break;
		case 'a':
			MixMatfile = optarg;
			break;
		case 'i':
			doInitPS = 1;
			break;
		default:
			cerr << "Option '" << (char)retval << "' not valid. Exiting.\n\n";
			usage(argv[0]);
		}
	}

	////////////////////////////////////////////////////////////////

	if (CORRon) printf("[%2.2i] CORRELATIONS BETWEEN DETECTORS INCLUDED\n", rank);
	if (!CORRon) printf("[%2.2i] NO CORRELATIONS BETWEEN DETECTORS INCLUDED\n", rank);


	// Set default parameter values
	if (ff_in.size() == 0) {
		ff_in.push_back(0);
		nf_in.push_back(-1);
	}
	if (ff_in.size() == 1 && nf_in.size() == 0)
		nf_in.push_back(-1);

	// Check improper usage
	if (dirfile == "") usage(argv[0]);
	if (channel.size() == 0) {
		cerr << "Must provide at least one channel.\n\n";
		usage(argv[0]);
	}
	if (ff_in.size() != nf_in.size()) {
		cerr << "'-l' or '-n' must follow '-f'. Exiting.\n";
		exit(1);
	}
	if (xxi_in.size() != xxf_in.size() || xxi_in.size() != yyi_in.size() || xxi_in.size() != yyf_in.size()) {
		cerr << "'-x' must be followed by '-X -z -Z'. Exiting.\n";
		exit(1);
	}
	if (tmpcount == 1 || tmpcount == 2 || tmpcount == 3){
		cerr << "ERROR: None or all the following keywords must be set: -u -U -v -V. Exiting. \n";
		exit(1);
	}
	if (tmpcount2 == 1 || tmpcount2 == 2){
		cerr << "ERROR: None or all the following keywords must be set: -t -T -N. Exiting. \n";
		exit(1);
	}


	if (f_lp_Nk == 0.0)
		f_lp_Nk = f_lp;

	if (napod){
		printf("[%2.2i] Data are apodized\n", rank);
	} else {
		printf("[%2.2i] Data are not apodized\n", rank);
	}


	if (pixdeg < 0){
		cerr << "ERROR: enter pixel size -p keyword\n";
		exit(1);
	}

	if (fsamp<=0.0){
		cerr << "ERROR: enter a correct sampling frequency -R keyword\n";
		exit(1);
	}

	ntotscan = ff_in.size();
	ndet = channel.size();

	nnf = extentnoiseSp_list.size();
	if (nnf != 1 && nnf != ntotscan){
		cerr << "ERROR: There should be one noise power spectrum file per scan, or a single one for all the scans. Check -K options" << endl;
		exit(1);
	}
	//  printf("%d\n",nnf);


	// convert lists to regular arrays
	fframes  = new long[ntotscan];
	fframesorder = new long[ntotscan];
	nsamples = new long[ntotscan];
	nsamplesorder = new long[ntotscan];
	ruleorder = new long[ntotscan];
	frnum = new long[ntotscan+1];
	bolonames = new string [ndet];
	extentnoiseSp_all = new string[ntotscan];
	//extentnoiseSp_allorder = new string[ntotscan];


	list2array(ff_in, fframes);
	list2array(nf_in, nsamples);
	list2array(channel, bolonames);
	list2array(extentnoiseSp_list,extentnoiseSp_all);
	if (nnf == 1 && ntotscan > 1)
		for (ii=1;ii<ntotscan;ii++)
			extentnoiseSp_all[ii] = extentnoiseSp_all[0];


	//  printf("xxi_in.size() = %d\n",xxi_in.size());


	if (xxi_in.size() != 0){
		xxi_boxes = new long[xxi_in.size()];
		xxf_boxes = new long[xxi_in.size()];
		yyi_boxes = new long[xxi_in.size()];
		yyf_boxes = new long[xxi_in.size()];
		list2array(xxi_in, xxi_boxes);
		list2array(xxf_in, xxf_boxes);
		list2array(yyi_in, yyi_boxes);
		list2array(yyf_in, yyf_boxes);
	}


	if(samples_per_frames>1){
		for (int ii=0; ii<ntotscan; ii++) {
			nsamples[ii] *= samples_per_frames;      // convert nframes to nsamples
		}
	}

	// utilisé lors de la lecture des coord de la map en pixel (dans la f° read_data)
	string ra_field;
	string dec_field;
	string phi_field;
	string scerr_field = "ERR"+pextension;
	string flpoint_field = "FLPOINTING";

	if (coordsyst == 2){
		ra_field = "L"+pextension;
		dec_field = "B"+pextension;
		phi_field = "PHIG"+pextension;
		printf("[%2.2i] Coordinate system: Galactic\n",rank );
	}else{
		ra_field = "RA"+pextension;
		dec_field = "DEC"+pextension;
		phi_field = "PHI"+pextension;
		if (coordsyst == 3){
			printf("[%2.2i] Map in Telescope coordinates. Reference coordinate system is RA/DEC (J2000)\n", rank);
		} else {
			printf("[%2.2i] Coordinate system: RA/DEC (J2000)\n", rank);
		}
	}


	if (NORMLIN)
		printf("NO BASELINE REMOVED\n");


	if (projgaps)
		printf("Flaged data are binned. iterative solution to fill gaps with noise only.\n");


	// path in which data are written
	if (pPath != NULL){
		poutdir = pPath;
	} else {
		poutdir = outdir;
	}

	///////////////////////////////////////////////////////////////////

	// allocate memory
	tancoord = new double[2];
	tanpix = new double[2];

	/********************* Define parallelization scheme   *******/

#ifdef USE_MPI
cout << "parallel_frames : 1 " << "   size : " << size  << endl;

	if (rank == 0){

		//mem alloc
		ruleorder = new long [ntotscan];
		frnum = new long [ntotscan];

		// read parallel schema in a file
		sprintf(testfile,"%s%s%s%s",outdir.c_str(),"parallel_for_Sanepic_",termin.c_str(),".txt");
		if ((fp = fopen(testfile,"w"))!=NULL){
			//fprintf(fp,"%d\n",size);
			fread(&size,sizeof(int), 1, fp);
			fread(ruleorder,sizeof(long),ntotscan,fp);
			fread(frnum,sizeof(long),ntotscan,fp);
			fclose(fp);
		}else{
			cerr << "Error : couldn't open file to write parallel options. Exiting" << endl;
			exit(1);
		}

		// reorder nsamples
		//find_best_order_frames(ruleorder,frnum,nsamples,ntotscan,size);


		for (ii=0;ii<ntotscan;ii++){
			nsamplesorder[ii] = nsamples[ruleorder[ii]];
			fframesorder[ii] = fframes[ruleorder[ii]];
		}
		for (ii=0;ii<ntotscan;ii++){
			nsamples[ii] = nsamplesorder[ii];
			fframes[ii] = fframesorder[ii];
		}

	}

#endif

	// read nn, coordsyst, tanpix, tancoord
	sprintf(testfile,"%s%s%s%s%d%s",outdir.c_str(),"InfoPointing_for_Sanepic_",termin.c_str(),"_",rank,".txt");
	if ((fp = fopen(testfile,"r")) == NULL){
		cerr << "File InfoPointing_for_sanepic :" << testfile << "not found. Exiting" << endl;
		exit(1);
	}

	fscanf(fp,"%d\n",&nn);
	fscanf(fp,"%d\n",&coordsyst2);
	fscanf(fp,"%lf\n",tanpix);
	fscanf(fp,"%lf\n",tanpix+1);
	fscanf(fp,"%lf\n",tancoord);
	fscanf(fp,"%lf\n",tancoord+1);
	fclose(fp);

	cout << "Map size :" << nn << "x" << nn << endl;
	//cout << tanpix[0] << " " << tanpix[1] << endl;
	//cout << tancoord[0] << " " << tancoord[1] << endl;

	if (coordsyst!=coordsyst2){
		cerr << "Error : coordinates systems must be the same for preprocessing and mapmaking" << endl;
		exit(0);
	}

	//******************************** some preprocess again // compute indpsrc and addnpix ****************/

	//************************************* Deal with masking the point sources
	// define the mask
	mask = new unsigned char[nn*nn];
	for (ii=0;ii<nn*nn;ii++)
		mask[ii] = 1;


	if (xxi_in.size() != 0){
		for (ib = 0;ib < (long)xxi_in.size(); ib++){ // to avoid warning, mat-27/05
			// for each box crossing constraint removal
			for (ii=xxi_boxes[ib];ii<xxf_boxes[ib];ii++)
				for (ll=yyi_boxes[ib];ll<yyf_boxes[ib];ll++)
					mask[ll*nn + ii] = 0;  // mask is initialised to 0
		}
	}



	//long npixsrc = 0;
	indpsrc = new long[nn*nn];
	for (ii=0;ii<nn*nn;ii++){
		if (mask[ii] == 0){
			indpsrc[ii] = npixsrc;
			npixsrc += 1;
		} else {
			indpsrc[ii] = -1;
		}
	}
	addnpix = ntotscan*npixsrc; // addnpix = number of pix to add in pixon = number of scans * number of pix in box crossing constraint removal

	cout  << "addnpix : " << addnpix << endl;
	//*************************************************************************************************//

	if (flgdupl) factdupl = 2; // -M =1, default 0 : if flagged data are put in a duplicated map

	// read npix, PNdtot from file
	sprintf(testfile,"%s%s%s%s",outdir.c_str(),"PNdCorr_",termin.c_str(),".bi");
	if ((fp = fopen(testfile,"r"))!=NULL){
		fread(&npix,sizeof(int),1,fp);
		PNdtot= new double[npix]; // mat 04/06
		PNd= new double[npix]; // mat 04/06
		fread(PNdtot,sizeof(double),npix,fp);
		fclose(fp);
	}else{
		cerr << "erreur : impossible de trouver le fichier " << testfile << endl;
		exit(0);
	}
	/*for (ii=0;ii<20;ii++)
		cout << PNdtot[ii] << " ";
	cout << endl << "avant read indpix\n";
	exit(0);*/


	indpix=new long[factdupl*nn*nn+2 + addnpix];
	// read indpix
	sprintf(testfile,"%s%s%s%s",outdir.c_str(),"Indpix_for_conj_grad_",termin.c_str(),".bi");
	if ((fp = fopen(testfile,"r"))!=NULL){
		fread(&flagon,sizeof(int),1,fp); // mat 04/06
		fread(indpix,sizeof(long), factdupl*nn*nn+2 + addnpix, fp);
		fclose(fp);
	}else{
		cerr << "Error : cannot find Indpix file " << testfile << endl;
		exit(0);
	}

	/*cout << "flagon : " << flagon << endl;
	for (ii=0;ii<nn*nn;ii=ii+nn+1)
		cout << indpix[ii] << " ";
	cout << endl;
	exit(1);*/
	/*
	if (rank == 0){
		sprintf(testfile,"%s%s%s%s",outdir.c_str(),"testfile_",termin.c_str(),".txt");
		fp = fopen(testfile,"a");
		for (ii=0;ii<=ntotscan;ii++) fprintf(fp,"frnum[%ld] = %ld \n",ii,frnum[ii]);
		fclose(fp);
	}*/


	// read npix et nn
	S = new double[npix];
	r = new double[npix];
	q = new double[npix];
	qtot = new double[npix];
	d = new double[npix];
	Mp = new double[npix];
	Mptot = new double[npix];
	s = new double[npix];
	PtNPmatS = new double[npix];
	PtNPmatStot = new double[npix];
	hits = new long[npix];
	hitstot = new long[npix];

	map1d = new double[nn*nn];




#ifdef USE_MPI

		MPI_Bcast(nsamples,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
		MPI_Bcast(fframes,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
		MPI_Bcast(frnum,ntotscan+1,MPI_LONG,0,MPI_COMM_WORLD);

		iframe_min = frnum[rank];
		iframe_max = frnum[rank+1];
		rank_det = 0;
		size_det = 1;

#else
		iframe_min = 0;
		iframe_max = ntotscan;
		rank_det = rank;
		size_det = size;
#endif
	/*************************************************************/

	printf("[%2.2i] iframe_min %ld\tiframe_max %ld \n",rank,iframe_min,iframe_max);

	if (iframe_min < 0 || iframe_min >= iframe_max || iframe_max > ntotscan){
		cerr << "Error distributing frame ranges. Check iframe_min and iframe_max. Exiting" << endl;
		exit(1);
	}



	/* END PARAMETER PROCESSING */




	printf("[%2.2i] Main Conjugate gradient loop\n",rank);



	for (idupl = 0;idupl<=flgdupl;idupl++){

		for (ii=0;ii<npix;ii++) S[ii] = 0.0;//PNd[ii];

		//Conjugate gradien Inversion
		if (projgaps || !flagon){
			npixeff = npix;
		} else {
			npixeff = npix-1;
		}



		printf("[%2.2i] npix = %d, npixeff = %d\n", rank, npix, npixeff);


		t1 = time(0);


		init1D_double(PtNPmatS,0,npix,0.0);
		init1D_double(PtNPmatStot,0,npix,0.0);
		init1D_double(Mp,0,npix,0.0);
		init1D_double(Mptot,0,npix,0.0);
		init1D_long(hits,0,npix,0);
		init1D_long(hitstot,0,npix,0);

		prefixe = "fPs";

		for (iframe=iframe_min;iframe<iframe_max;iframe++){
			ns = nsamples[iframe];
			ff = fframes[iframe];
			f_lppix_Nk = f_lp_Nk*double(ns)/fsamp;
			//prefixe = "fPs"; now outside of the loop

			//    cout << "[" << rank << "] " << iframe << "/" << iframe_max << endl;

			if (CORRon){
				write_tfAS(S,indpix,nn,npix,flgdupl,factdupl, poutdir,termin,ff,ns,ndet,iframe);
				// read pointing + deproject + fourier transform

				do_PtNd(PtNPmatS,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,
						f_lppix_Nk,fsamp,ff,ns,ndet,size_det,rank_det,indpix,nn,npix,iframe,Mp,hits);
				// return Pnd = At N-1 d
			} else {

				do_PtNPS_nocorr(S,extentnoiseSp_all,noiseSppreffile,poutdir,termin,dirfile,bolonames,
						f_lppix_Nk,fsamp,flgdupl,factdupl,ff,ns,ndet,size_det,rank_det,indpix,
						nn,npix,iframe,PtNPmatS,Mp,hits);
			}

		} // end of iframe loop




#ifdef USE_MPI
		MPI_Reduce(PtNPmatS,PtNPmatStot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(hits,hitstot,npix,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(Mp,Mptot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
		PtNPmatStot=PtNPmatS; // ajout Mat 02/06
		hitstot=hits;
		Mptot=Mp;
#endif


		if (rank == 0) {

			//     t2 = time(0);
			//     printf("temps de calcul: %ld\n",t2-t1);

			for (ii=0;ii<npixeff;ii++)
				if (Mptot[ii] == 0)
					printf("ERROR: Mp[%ld] has elements = 0\n",ii);


			for (ii=0;ii<npixeff;ii++)
				Mptot[ii] = 1.0/Mptot[ii];


			for (ii=0;ii<npixeff;ii++)
				r[ii] = PNdtot[ii] - PtNPmatStot[ii];

			for (ii=0;ii<npixeff;ii++)
				d[ii] =  Mptot[ii] * r[ii];


			delta_n = 0.0;
			for (ii=0;ii<npixeff;ii++)
				delta_n += r[ii]*d[ii];

			var_n = 0.0;
			for (ii=0;ii<npixeff;ii++)
				var_n += r[ii]*r[ii];


			delta0 = delta_n;
			var0 = var_n;
			printf("[%2.2i] var0 = %lf\n",rank, var0);

		}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&var0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&var_n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(d,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

		printf("[%2.2i] Main Conjugate gradient loop started\n",rank);


		//start loop
		iter = 0;
		while(((iter < 2000) && (var_n/var0 > 1e-10) && (idupl || !flgdupl)) || (!idupl && var_n/var0 > 1e-4)){
			// added brackets in order to avoid warning, mat-27/05
			init1D_double(q,0,npixeff,0.0);
			init1D_double(qtot,0,npixeff,0.0);
cout << "dans while\n";
exit(0);
			//t1 = time(0);
			//sprintf(testfile,"%s%s%s%s",outdir.c_str(),"Timing_",termin.c_str(),".txt");
			//fp = fopen(testfile,"a");
			//fprintf(fp,"starting while loop: t = %ld\n",t1);
			//fclose(fp);

			prefixe = "fPs";

			for (iframe=iframe_min;iframe<iframe_max;iframe++){
				ns = nsamples[iframe];
				ff = fframes[iframe];
				f_lppix_Nk = f_lp_Nk*double(ns)/fsamp;
				//prefixe = "fPs";

				if (CORRon){
					write_tfAS(d,indpix,nn,npix,flgdupl,factdupl, poutdir,termin,ff,ns,ndet,iframe);
					// read pointing + deproject + fourier transform

					do_PtNd(q,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,f_lppix_Nk,
							fsamp,ff,ns,ndet,size_det,rank_det,indpix,nn,npix,iframe,NULL,NULL);
					// return Pnd = At N-1 d
				} else {

					do_PtNPS_nocorr(d,extentnoiseSp_all,noiseSppreffile,poutdir,termin,dirfile,bolonames,
							f_lppix_Nk,fsamp,flgdupl,factdupl,ff,ns,ndet,size_det,rank_det,indpix,
							nn,npix,iframe,q,NULL,NULL);
				}
			} // end of iframe loop




#ifdef USE_MPI
			MPI_Reduce(q,qtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
			qtot=q; // ajout mat 02/06
#endif



			if (rank == 0){
				rtq= 0.0;
				for (ii=0;ii<npixeff;ii++)
					rtq += qtot[ii] * d[ii];

				alpha = delta_n/rtq;


				for (ii=0;ii<npixeff;ii++)
					S[ii] += alpha*d[ii];
			}

#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(S ,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif



			if ((iter % 10) == 0){
				init1D_double(PtNPmatS,0,npixeff,0.0);
				init1D_double(PtNPmatStot,0,npixeff,0.0);

				prefixe = "fPs";

				for (iframe=iframe_min;iframe<iframe_max;iframe++){
					ns = nsamples[iframe];
					ff = fframes[iframe];
					f_lppix_Nk = f_lp_Nk*double(ns)/fsamp;
					//prefixe = "fPs";

					if (CORRon){
						write_tfAS(S,indpix,nn,npix,flgdupl,factdupl, poutdir,termin,ff,ns,ndet,iframe);
						// read pointing + deproject + fourier transform

						do_PtNd(PtNPmatS,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,
								f_lppix_Nk,fsamp,ff,ns,ndet,size_det,rank_det,indpix,nn,npix,iframe,
								NULL,NULL);
						// return Pnd = At N-1 d
					} else {

						do_PtNPS_nocorr(S,extentnoiseSp_all,noiseSppreffile,poutdir,termin,dirfile,bolonames,
								f_lppix_Nk,fsamp,flgdupl,factdupl,ff,ns,ndet,size_det,rank_det,
								indpix,nn,npix,iframe,PtNPmatS,NULL,NULL);
					}
				} // end of iframe loop



#ifdef USE_MPI
				MPI_Reduce(PtNPmatS,PtNPmatStot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
				PtNPmatStot=PtNPmatS;
#endif


				if (rank == 0){
					for (ii=0;ii<npixeff;ii++)
						r[ii] = PNdtot[ii] - PtNPmatStot[ii];
				}


			} else {

				if (rank == 0){
					for (ii=0;ii<npixeff;ii++)
						r[ii] -= alpha*qtot[ii];
				}
			}





			if (rank == 0){

				for (ii=0;ii<npixeff;ii++)
					s[ii] = Mptot[ii]*r[ii];


				delta_o = delta_n;

				delta_n = 0.0;
				for (ii=0;ii<npixeff;ii++)
					delta_n += r[ii]*s[ii];

				var_n = 0.0;
				for (ii=0;ii<npixeff;ii++)
					var_n += r[ii]*r[ii];



				beta = delta_n/delta_o;
				for (ii=0;ii<npixeff;ii++)
					d[ii] = s[ii] + beta*d[ii];


				cout << "iter = " << iter;
				cout << ", crit  = " << setiosflags(ios::scientific) << setiosflags(ios::floatfield) << var_n/var0;
				cout << ", crit2 = " << setiosflags(ios::scientific) << setiosflags(ios::floatfield) << delta_n/delta0;
				cout << "\r " << flush;

				//      printf("[%2.2i] iter = %d, crit = %10.15g, crit2 = %10.15g     \n",rank, iter,var_n/var0,delta_n/delta0);



				if (iter == 0){
					for (ii=0; ii<nn; ii++) {
						for (jj=0; jj<nn; jj++) {
							mi = jj*nn + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = Mptot[indpix[mi]];
							} else {
								map1d[mi] = 0.0;
							}
						}
					}

					fname = '!' + outdir + "optimMap_" + termin + "_noisevar.fits";
					write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


					for (ii=0; ii<nn ; ii++){
						for (jj=0; jj<nn; jj++){
							mi = jj*nn + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = Mptot[indpix[mi]] * PNdtot[indpix[mi]];
							} else {
								map1d[mi] = 0.0;
							}
						}
					}
					fname = '!' + outdir + "binMap_" + termin + "_flux.fits";
					write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);



					for (ii=0; ii<nn; ii++) {
						for (jj=0; jj<nn; jj++) {
							mi = jj*nn + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = hitstot[indpix[mi]];
							} else {
								map1d[mi] = 0.0;
							}
						}
					}

					if (addnpix){
						for (iframe = 0;iframe<ntotscan;iframe++){
							for (ii=0; ii<nn; ii++) {
								for (jj=0; jj<nn; jj++) {
									mi = jj*nn + ii;
									if ((indpsrc[mi] != -1) && (indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc] != -1))
										map1d[mi] += hitstot[indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc]];
								}
							}
						}
					}

					fname = '!' + outdir + "optimMap_" + termin + "_hits.fits";
					write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


					for (ii=0; ii<nn; ii++) {
						for (jj=0; jj<nn; jj++) {
							mi = jj*nn + ii;
							map1d[mi] = 0.0;
						}
					}

					if (addnpix){
						for (iframe = 0;iframe<ntotscan;iframe++){
							for (ii=0; ii<nn; ii++) {
								for (jj=0; jj<nn; jj++) {
									mi = jj*nn + ii;
									if ((indpsrc[mi] != -1) && (indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc] != -1))
										map1d[mi] += 1.0/Mptot[indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc]];
								}
							}
						}

						fname = '!' + outdir + "optimMap_" + termin + "_invnoisevaruncpix.fits";
						write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
					}

				}


				if (iterw && (iter % iterw) == 0){

					// make the map
					for (ii=0; ii<nn; ii++) {
						for (jj=0; jj<nn; jj++) {
							mi = jj*nn + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = -S[indpix[mi]];
							} else {
								map1d[mi] = 0.0;
							}
						}
					}

					sprintf(iterchar,"%d",iter);
					iterstr = iterchar;
					fname = '!' + outdir + "optimMap_" + termin + "_flux" + iterstr + "b.fits";
					write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


					if (flgdupl){
						for (ii=0; ii<nn; ii++) {
							for (jj=0; jj<nn; jj++) {
								mi = jj*nn + ii;
								if (indpix[mi] >= 0){
									map1d[mi] = -S[indpix[mi+nn*nn]]; //-finalmap[ii][jj];
								} else {
									map1d[mi] = 0.0;
								}
							}
						}


						sprintf(iterchar,"%d",iter);
						iterstr = iterchar;
						fname = '!' + outdir + "optimMap_" + termin + "_fluxflags" + iterstr + "b.fits";
						write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
					}
					if (addnpix){
						for (ii=0; ii<nn ; ii++){
							for (jj=0; jj<nn ; jj++){
								mi = jj*nn + ii;
								map1d[mi] = 0.0;
							}
						}
						for (iframe = 0;iframe<ntotscan;iframe++){
							for (ii=0; ii<nn; ii++) {
								for (jj=0; jj<nn; jj++) {
									mi = jj*nn + ii;
									if ((indpsrc[mi] != -1)  && (indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc] != -1)){
										map1d[mi] += -S[indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc]]/Mptot[indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc]];
									}
								}
							}
						}


						sprintf(iterchar,"%d",iter);
						iterstr = iterchar;
						fname = '!' + outdir + "optimMap_" + termin + "_fluxuncpix_" + iterstr + "b.fits";
						write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);

					}
				}



				sprintf(testfile,"%s%s%s%s",outdir.c_str(),"ConvFile_",termin.c_str(),".txt");
				fp = fopen(testfile,"a");
				fprintf(fp,"iter = %d, crit = %10.15g, crit2 = %10.15g\n",iter,var_n/var0, delta_n/delta0);
				fclose(fp);

			}

#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&var_n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(d ,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

			iter++;

		} // end of while loop
		printf("\n");




		sprintf(testfile,"%s%s%s%s",outdir.c_str(),"testfile_",termin.c_str(),".txt");
		fp = fopen(testfile,"a");
		fprintf(fp,"test en sortie de la boucle while \n");
		fclose(fp);





		if  ((projgaps || (flgdupl)) && !idupl){

			init1D_double(PNd,0,npix,0.0);
			init1D_double(PNdtot,0,npix,0.0);

			prefixe = "fdata";

			for (iframe=iframe_min;iframe<iframe_max;iframe++){

				ns = nsamples[iframe];
				ff = fframes[iframe];
				f_lppix = f_lp*double(ns)/fsamp;
				f_lppix_Nk = f_lp_Nk*double(ns)/fsamp;
				//prefixe = "fdata"; now out of loop

				if (CORRon){

					write_ftrProcesdata(S,indpix,indpsrc,nn,npix,npixsrc,ntotscan,addnpix,flgdupl,factdupl,2,
							poutdir,termin,errarcsec,dirfile,scerr_field,flpoint_field,bolonames,
							bextension,fextension,/*cextension,*/shift_data_to_point,f_lppix,ff,ns,
							napod,ndet,NORMLIN,NOFILLGAP,iframe);// fillgaps + butterworth filter + fourier transform
					// "fdata_" files generation (fourier transform of the data)

					do_PtNd(PNd,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,f_lppix_Nk,
							fsamp,ff,ns,ndet,size_det,rank_det,indpix,nn,npix,iframe,NULL,NULL);
					// return Pnd = At N-1 d
				} else {

					do_PtNd_nocorr(PNd,extentnoiseSp_all,noiseSppreffile,poutdir,termin,errarcsec,dirfile,
							scerr_field,flpoint_field,bolonames,bextension,fextension,
							/*cextension,*/shift_data_to_point,f_lppix,f_lppix_Nk,fsamp,ntotscan,addnpix,
							flgdupl,factdupl,2,ff,ns,napod,ndet,size_det,rank_det,indpix,indpsrc,
							nn,npix,npixsrc,NORMLIN,NOFILLGAP,iframe,S);
				}
			} // end of iframe loop




#ifdef USE_MPI
			MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
			PNdtot=PNd; // ajout Mat 02/06
#endif
		}



	}// end of idupl loop



#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	sprintf(testfile,"%s%s%s%s",outdir.c_str(),"testfile_",termin.c_str(),".txt");
	fp = fopen(testfile,"a");
	fprintf(fp,"test avant ecriture \n");
	fclose(fp);





	//******************************  write final map in file ********************************



	if (rank == 0){

		printf(" after CC INVERSION %d\n",npix*(npix+1)/2);



		bool fru;

		for (ii=0; ii<nn; ii++) {
			for (jj=0; jj<nn; jj++) {
				mi = jj*nn + ii;
				if (indpix[mi] >= 0){
					map1d[mi] = -S[indpix[mi]];
				} else {
					map1d[mi] = 0.0;
				}
			}
		}


		fname = '!' + outdir + "optimMap_" + termin + "_flux.fits";
		write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


		for (ii=0; ii<nn; ii++) {
			for (jj=0; jj<nn; jj++) {
				mi = jj*nn + ii;
				if (indpix[mi] >= 0){
					map1d[mi] = Mptot[indpix[mi]];
				} else {
					map1d[mi] = 0.0;
				}
			}
		}


		fname = '!' + outdir + "optimMap_" + termin + "_noisevar.fits";
		write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);

		if (addnpix){
			for (iframe = 0;iframe<ntotscan;iframe++){
				fru = 0;
				for (ii=0; ii<nn; ii++) {
					for (jj=0; jj<nn; jj++) {
						mi = jj*nn + ii;
						if ((indpsrc[mi] != -1)  && (indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc] != -1)){
							fru = 1;
							map1d[mi] = -S[indpix[indpsrc[mi]+factdupl*nn*nn+npixsrc*iframe]];
						} else {
							map1d[mi] = 0.0;
						}
					}
				}


				if (fru){
					sprintf(iframechar,"%ld",iframe);
					iframestr = iframechar;
					fname = '!' + outdir + "optimMap_" + termin + "_flux_fr" + iframestr + ".fits";
					write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


					for (ii=0; ii<nn; ii++) {
						for (jj=0; jj<nn; jj++) {
							mi = jj*nn + ii;
							if ((indpsrc[mi] != -1)  && (indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc] != -1)){
								map1d[mi] = Mptot[indpix[indpsrc[mi]+factdupl*nn*nn+npixsrc*iframe]];
							} else {
								map1d[mi] = 0.0;
							}
						}
					}

					fname = '!' + outdir + "optimMap_" + termin + "_noisevar_fr" + iframestr + ".fits";
					write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
				}
			}
		}
	}



	//*******************************************************************//
	//******************  Update noise power spectra  *******************//

if (doInitPS){
	printf("%s\n",MixMatfile.c_str());

	if (MixMatfile != "NOFILE"){
		for (iframe=iframe_min;iframe<iframe_max;iframe++){
			ns = nsamples[iframe];
			ff = fframes[iframe];
			extentnoiseSp = extentnoiseSp_all[iframe];

			EstimPowerSpectra(fsamp,ns,ff,ndet,nn,npix,napod,iframe,flgdupl,factdupl,indpix,S,
					MixMatfile,bolonames,dirfile,bextension,fextension/*,cextension*/,shift_data_to_point,
					poutdir,termin,NORMLIN,NOFILLGAP,noiseSppreffile,extentnoiseSp,outdir);

		}
	}
}



	//******************************************************************//
	//******************************************************************//
	//*********************** End of program *************************//
	//******************************************************************//
	//******************************************************************//





	if (rank == 0){


		//write infos for second part
		sprintf(testfile,"%s%s%s%s",outdir.c_str(),"InfoFor2ndStep_",termin.c_str(),".txt");
		if((fp = fopen(testfile,"w"))==NULL){
			cerr << "Cannot open file :" << testfile << "\tExiting." << endl;
			exit(1);
		}
		fprintf(fp,"%d\n",nn);
		fprintf(fp,"%d\n",npix);
		fprintf(fp,"%lf\n",pixdeg);
		fprintf(fp,"%lf\n",tancoord[0]);
		fprintf(fp,"%lf\n",tancoord[1]);
		fprintf(fp,"%lf\n",tanpix[0]);
		fprintf(fp,"%lf\n",tanpix[1]);
		fprintf(fp,"%d\n",coordsyst);
		fprintf(fp,"%i\n",flagon);
		fprintf(fp,"\n");
		for (ii=0;ii<nn*nn;ii++)
			fprintf(fp,"%ld\n",indpix[ii]);
		fclose(fp);



		//write command line in a file
		sprintf(testfile,"%s%s%s%s",outdir.c_str(),"CommandLine_",termin.c_str(),".txt");
		fp = fopen(testfile,"w");
		fprintf(fp,"%d\n",argc);
		for (ii=1;ii<argc;ii++)
			fprintf(fp,"%s \t",argv[ii]);
		fprintf(fp,"\n");
		fclose(fp);

	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif


	// clean up





#ifdef USE_MPI
	MPI_Finalize();
#endif



	return 0;
}



