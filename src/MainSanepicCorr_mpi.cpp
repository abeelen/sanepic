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
	cerr << "-J <noise freq cut>    frequency under which noise power spectra are thresholded. Default is the frequency cut of the high pass filter applied to the data (-H option)" << endl;
	cerr << "\t[-O <offset file>] [-B <bolo ext>] [-G <flag ext>] [-P <pointing ext>]" << endl;
	cerr << "\t[-o <output dir>] [-e <output file>]" << endl << endl;
	cerr << "-A <apodize nsamp>     number of samples to apodize (supercedes -a)\n";
	//cerr << "-m <padding interval>  number of samples extrapolated before (and after) data" << endl;
	cerr << "-B <bolo ext>          bolo field extension" << endl;
	cerr << "-C <channel file>      include bolometers named in file" << endl;
	cerr << "-e <output file str>   id string in output file" << endl;
	cerr << "-F <data dir>          source data directory" << endl;
	cerr << "-f <first frame>       first frame (multiple allowed)" << endl;
	cerr << "-G <flag ext>          flag field extension" << endl;
	cerr << "-H <filter freq>       filter frequency (multiple allowed)" << endl;
	cerr << "-K <noise prefixe>     noise power spectrum file prefixe " << endl;
	cerr << "-k <noise suffixe>     noise power spectrum file suffixe " << endl;
	cerr << "-l <last frame>        last frame (multiple allowed, exclusive from -n)" << endl;
	cerr << "-n <num frames>        number of frames (multiple allowed, exclusive from -l)" << endl;
	cerr << "-x <box coord x1>      if set, defines the x coordinate of the left side of the box(es) in which crossing constraints are removed" << endl;
	cerr << "-X <box coord x2>      if set, defines the x coordinate of the right side of the box(es) in which crossings constraints are removed" << endl;
	cerr << "-z <box coord y1>      if set, defines the y coordinate of the top side of the box(es) in which crossing constraints are removed" << endl;
	cerr << "-Z <box coord y2>      if set, defines the y coordinate of the bottom side of the box(es) in which crossings constraints are removed" << endl;
	cerr << "-O <offset file>       file containing bolometer offsets" << endl;
	cerr << "-o <output dir>        output directory" << endl;
	cerr << "-P <pointing ext>      pointing field extension" << endl;
	cerr << "-p <pixsize>           size of pixels (degrees)" << endl;
	cerr << "-y <channel>           include bolometer (multiple allowed)" << endl;
	cerr << "-S <Scan offset file>  file containing offsets for different frame ranges" << endl;
	cerr << "-t <RA source>         RA (or Gal l if keyword c is 2) of the tangent point and of the source for telescope coordinate maps" << endl;
	cerr << "-T <DEC source>        DEC (or Gal b if keyword c is 2) of the tangent point and of the source for telescope coordinate maps" << endl;
	cerr << "-u <RA min>            fixed minimum RA (or Gal l) of the map, if none of -N, -t, and -T are defined" << endl;
	cerr << "-U <RA max>            fixed maximum RA (or Gal l) of the map, if none of -N, -t, and -T are defined" << endl;
	cerr << "-v <DEC min>           fixed minimum DEC (or Gal b) of the map, if none of -N, -t, and -T are defined" << endl;
	cerr << "-V <DEC max>           fixed maximum DEC (or Gal b) of the map, if none of -N, -t, and -T are defined" << endl;
	cerr << "-c <coord system>      Coordinate system, 1: RA/DEC, 2: L/B (Gal), 3: Telescope coordinates, Default is RA/DEC" << endl;
	cerr << "-N <map radius>        fixed radius (half a side) of the map in degrees (in case of data outside, they are constrained to have constant values)" << endl;
	cerr << "-L <no baseline>       Keyword specifiying if a baseline is removed from the data or not (0 for YES, 1 for NO)" << endl;
	cerr << "-r <precompute PNd>    Keyword specifying if PNd is precomputed and read from disk" << endl;
	cerr << "-g <project gaps>      Keyword specifying if gaps are projected to a pixel in the map, if so gap filling of noise only is performed iteratively. Default is 0" << endl;
	cerr << "-M <map flagged data>   Keyword specifying if flagged data are put in a separate map, default is 0" << endl;
	cerr << "-s <time offset>       Keyword for subtracting a time offset to the data to match the pointing "<< endl;
	cerr << "-E <no corr>           Set this keyword to 0 if correlations are not included in the analysis" << endl;
	cerr << "-I <parallel scheme>   Set this keyword to 1 in order to distribute different field visit (or frame ranges) to different processors, instead of different detector data chuncks as it is by default" << endl;
	cerr << "-j <write unconverged maps>  The value specified by this keyword indicates the period in iterations to which the data are written to disk. Default is 10. Set this keyword to zero if you don't want intermediate maps to be written." << endl;
	cerr << "-a <noise estim>       Optional. Enter filename containing the mixing matrix of noise components. If set, the noise power spectra files for each field are computed after map-making. Default is no re-estimation" << endl;
	cerr << "-i <Noise PS etimation> Optional. Computes power spectra estimation from the data and saves it in a file. Default is no estimation. You must run sanepic again without -i option to compute mapmaking" << endl;
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
//*************************** Beginning of main program ****************************//
//**********************************************************************************//
//**********************************************************************************//




int main(int argc, char *argv[])
{



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
	size = 0;
	rank = 0;
	cout << "pas de mpi" << endl;
#endif

	char * pPath;
	pPath = getenv ("TMPBATCH");
	if (pPath!=NULL)
		printf ("The current path is: %s",pPath);


	bool projgaps = 0; //1: data flagged are put in a single pixel
	//   (assume no signal in this pixel),
	//0: data flagged are not reprojected

	//default value of the data to pointing shift
	int shift_data_to_point = 0;


	//DEFAULT PARAMETERS
	long napod = 0; // number of samples to apodize
	double fsamp = 25.0; // sampling frequency : BLAST Specific
	double errarcsec = 15.0; // rejection criteria : scerr[ii] > errarcsec, sample is rejected
	// source error

	long ii, jj, ll, iframe, idet, ib; // loop indices
	long iframe_min, iframe_max;
	int flagon = 0; // if rejectsample [ii]==3, flagon=1
	int iterw = 10; // period in iterations to which the data are written to disk, 0 = no intermediate map to be written
	bool bfixc = 0; // indicates that 4 corners are given for the cross corelation removal box
	bool pixout = 0; // indicates that at least one pixel has been flagged and is out
	bool NORMLIN = 0; // baseline is removed from the data, NORMLIN = 1 else 0
	bool NOFILLGAP = 0; // fill the gap ? default is YES
	bool PND_ready = 0; // PNd precomputed ? read on disk if =1
	bool flgdupl = 0; // 1 if flagged data are put in a separate map
	bool CORRon = 1; // correlation included in the analysis (=1), else 0, default 0
	bool parallel_frames = 0; // a parallel scheme is used : mpi has been launched

	//set coordinate system
	double *srccoord, *coordscorner;
	srccoord = new double[2]; // RA/DEC tangent point/source
	coordscorner = new double[4]; // map min/max RA/DEC coords (-N,-t,-T absents)
	srccoord[0] = -1000; // RA tangent point/source
	srccoord[1] = -1000; // DEC tangent point/source
	double radius = -1.0; // map radius (half a side) in degrees


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
	int samples_per_frames=1; // default = 1, BLAST = 20
	samples_per_frames=20; // Blast specific, a supprimer par la suite


	// map making parameters
	double pixdeg; // size of pixels (degree)


	int nn, npix; // nn = side of the map, npix = number of filled pixels
	double ra_min, ra_max, dec_min, dec_max; // coord ra/dec of the map
	double *offsets, *froffsets, *offmap; // offsets par rapport au bolo de ref, froffsets = ?, offmap = ?
	double *tancoord; // tangent point coordinates
	double *tanpix; // tangent pixel
	double gra_min, gra_max, gdec_min, gdec_max; // global ra/dec min and max (to get the min and max of all ra/dec min/max computed by different processors)

	//internal data params
	long ns, ff; // number of samples for this scan, first frame number of this scan
	double f_lp, f_lp_Nk, f_lppix, f_lppix_Nk; // frequencies : filter knee freq, noise PS threshold freq ; frequencies converted in a number of samples


	FILE *fp;

	char testfile[100];



	char type='d'; // returned type of read_data functions, d=64bit double
	double *ra, *dec, *phi, *scerr; // RA/DEC, phi (angle) coordinates of the bolo, source errors
	unsigned char *flag, *flpoint, *rejectsamp, *mask; // samples flags, pointing flags, rejected samples list
	// mask = box for crossing constraint removal mask in the map
	double *PNd, *PNdtot; //
	long *indpix, *indpsrc; // pixels indices, mask pixels indices

	int *xx, *yy; // data coordinates in the map
	long *pixon; // this array is used to store the rules for pixels : they are seen or not
	long *samptopix; // sample to pixel conversion array




	string field; // actual boloname in the bolo loop
	string *bolonames; // channel list -> bolonames array : considered bolometers names
	string *extentnoiseSp_all; // ((list -> string*))
	string *extentnoiseSp_allorder;
	string bolofield; // bolofield = boloname + bextension
	string calfield; // calfield  = field+cextension;
	string flagfield; // flagfield = field+fextension;
	string dirfile; // data directory
	string outdir; // output directory
	string poutdir; // current path (pPath) or output dir (outdir)
	string bextension; // bolometer field extension
	string cextension = "NOCALP"; // needed for calfield calculation !
	string fextension = "NOFLAG"; // flag field extension
	string pextension; // pointing extension
	string file_offsets; // bolometer offsets file
	string file_frame_offsets = "NOOFFS"; // offset file
	string termin; // output file suffix
	string noiseSppreffile; // noise file suffix
	string extentnoiseSp; // noise file
	string prefixe; // prefix used for temporary name file creation

	string MixMatfile = "NOFILE";
	int doInitPS = 0;

	/* DEFAULT PARAMETERS */
	int coordsyst = 1; /// Default is RA/DEC
	int coordsyst2 = 1;



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
		case 'D':
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
			cextension = optarg;
			break;
		case 'G':
			fextension = optarg;
			break;
		case 'P':
			pextension = optarg;
			break;
		case 'O':
			file_offsets = optarg;
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
		case 'S':
			file_frame_offsets = optarg;
			break;
		case 't':
			srccoord[0] = atof(optarg);
			tmpcount2 += 1;
			break;
		case 'T':
			srccoord[1] = atof(optarg);
			tmpcount2 += 1;
			break;
		case 'u':
			coordscorner[0] = atof(optarg);
			tmpcount += 1;
			break;
		case 'U':
			coordscorner[1] = atof(optarg);
			tmpcount += 1;
			break;
		case 'v':
			coordscorner[2] = atof(optarg);
			tmpcount += 1;
			break;
		case 'V':
			coordscorner[3] = atof(optarg);
			tmpcount += 1;
			break;
		case 'c':
			coordsyst = atoi(optarg);
			break;
		case 'N':
			radius = atof(optarg);
			tmpcount2 += 1;
			break;
		case 'L':
			NORMLIN = atoi(optarg);
			break;
		case 'g':
			projgaps = atoi(optarg);
			break;
		case 'r':
			PND_ready = atoi(optarg);
			break;
		case 'M':
			flgdupl = atoi(optarg);
			break;
		case 's':
			shift_data_to_point = atoi(optarg);
			break;
		case 'E':
			CORRon = atoi(optarg);
			break;
		case 'I':
			parallel_frames = atoi(optarg);
			break;
		case 'j':
			iterw = atoi(optarg);
			break;
		case 'a':
			MixMatfile = optarg;
			break;
		case 'i':
			doInitPS = atoi(optarg);
			break;
		default:
			cerr << "Option '" << (char)retval << "' not valid. Exiting.\n\n";
			usage(argv[0]);
		}
	}



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
	extentnoiseSp_allorder = new string[ntotscan];


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
	if (tmpcount == 4)
		bfixc = 1;

	if (tmpcount2 == 3){
		if (tmpcount == 4){
			cerr << "ERROR: Conflicting input parameter: -u -U -v -V keywords are not compatible with -N -t -T keywords . Exiting. \n";
			exit(1);
		}
		bfixc = 1;
		coordscorner[0] = srccoord[0];
		coordscorner[1] = srccoord[0];
		coordscorner[2] = srccoord[1];
		coordscorner[3] = srccoord[1];
	}

	if (coordsyst != 3){
		srccoord[0] = -1000;
		srccoord[1] = -1000;
	}
	if ((coordsyst == 3) && (tmpcount2 != 3)){
		cerr << "ERROR: You must provide coordinates of the source in RA/DEC for telescope coordinates, use -t -T -N\n";
		exit(1);
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


	// map offsets
	float scoffsets[6]; // offsets depending on wavelength
	int nfoff; // number of offsets
	foffset *foffsets; // tableau d'offsets
	if (file_frame_offsets != "NOOFFS") // option -S file containing offsets for different frame ranges
		foffsets = read_mapoffsets(file_frame_offsets, scoffsets, &nfoff); //-> return foffsets, scoffsets, nfoff
	else {
		printf("[%2.2i] No offsets between visits\n", rank);
		nfoff = 2;
		ff = fframes[0];
		for (ii=0;ii<ntotscan;ii++) if (fframes[ii] > ff) ff = fframes[ii];
		foffsets = new foffset [2];
		(foffsets[0]).frame = 0;
		(foffsets[0]).pitch = 0.0;
		(foffsets[0]).yaw   = 0.0;
		(foffsets[1]).frame = ff+1; // ff = last frame number here
		(foffsets[1]).pitch = 0.0;
		(foffsets[1]).yaw   = 0.0;
		for (ii=0;ii<6;ii++)
			scoffsets[ii] = 0.0;
	}


	// path in which data are written
	if (pPath != NULL){
		poutdir = pPath;
	} else {
		poutdir = outdir;
	}
	printf("[%2.2i] Data written in %s\n",rank, poutdir.c_str());





	/********************* Define parallelization scheme   *******/

	if (parallel_frames && (rank == 0)){

		// reorder nsamples
		find_best_order_frames(ruleorder,frnum,nsamples,ntotscan,size);

		for (ii=0;ii<ntotscan;ii++){
			nsamplesorder[ii] = nsamples[ruleorder[ii]];
			fframesorder[ii] = fframes[ruleorder[ii]];
			//extentnoiseSp_allorder[ii] = extentnoiseSp_all[ii];
		}
		for (ii=0;ii<ntotscan;ii++){
			nsamples[ii] = nsamplesorder[ii];
			fframes[ii] = fframesorder[ii];
			//extentnoiseSp_all[ii] = extentnoiseSp_allorder[ii];
			//printf("frnum[%d] = %d\n",ii,frnum[ii]);
		}

	}




	if (rank == 0){
		sprintf(testfile,"%s%s%s%s",outdir.c_str(),"testfile_",termin.c_str(),".txt");
		fp = fopen(testfile,"a");
		for (ii=0;ii<=ntotscan;ii++) fprintf(fp,"frnum[%ld] = %ld \n",ii,frnum[ii]);
		fclose(fp);


	// write parallel schema in a file
		sprintf(testfile,"%s%s%s%s",outdir.c_str(),"parallel_for_Sanepic_",termin.c_str(),".txt");
		if ((fp = fopen(testfile,"w"))!=NULL){
		//fprintf(fp,"%d\n",size);
		fwrite(&size,sizeof(int), 1, fp);
		fwrite(ruleorder,sizeof(long),ntotscan,fp);
		fwrite(frnum,sizeof(long),ntotscan,fp);
		fclose(fp);
		}else{
			cerr << "Error : couldn't open file to write parallel options. Exiting" << endl;
			exit(1);
		}
	}

	if (parallel_frames){
#ifdef USE_MPI

		MPI_Bcast(nsamples,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
		MPI_Bcast(fframes,ntotscan,MPI_LONG,0,MPI_COMM_WORLD);
		MPI_Bcast(frnum,ntotscan+1,MPI_LONG,0,MPI_COMM_WORLD);
#endif
		iframe_min = frnum[rank];
		iframe_max = frnum[rank+1];
		rank_det = 0;
		size_det = 1;
	} else {
		iframe_min = 0;
		iframe_max = ntotscan;
		rank_det = rank;
		size_det = size;
	}
	/*************************************************************/

	printf("[%2.2i] iframe_min %ld\tiframe_max %ld \n",rank,iframe_min,iframe_max);

	if (iframe_min < 0 || iframe_min >= iframe_max || iframe_max > ntotscan){
		cerr << "Error distributing frame ranges. Check iframe_min and iframe_max. Exiting" << endl;
		exit(1);
	}



	/* END PARAMETER PROCESSING */




	/********** Alocate memory ***********/
	printf("[%2.2i] Allocating Memory\n",rank);

	ns = nsamples[0];
	// seek maximum number of samples
	for (ii=0;ii<ntotscan;ii++) if (nsamples[ii] > ns) ns = nsamples[ii];

	ra = new double[2*ns]; // RA bolo de ref
	dec = new double[2*ns]; // DEc du bolo de ref
	phi = new double[2*ns]; // (du bolo de ref) angle de la matrice de detecteur par rapport a RA/dec
	scerr = new double[2*ns]; // BLAST SPECIFIC : mesure l'erreur de pointage, si trop grande on flag la donnée
	xx = new int[2*ns]; // sample column coordinates in the map
	yy = new int[2*ns]; // sample row coordinates in the map
	samptopix = new long[2*ns]; // sample to pixel conversion index
	flag = new unsigned char[2*ns]; // flag data => =1
	rejectsamp = new unsigned char[2*ns]; // rejected samples after flag conditions
	flpoint = new unsigned char[2*ns]; // flpoint est un flag du pointage/time. Savoir au temps t, si tu prends ces données là, ou non.
	tancoord = new double[2]; // coordinates in ra/dec of the tangent point
	tanpix = new double[2]; // coordinates in the map of the tangent point
	offsets = new double[2]; //
	froffsets = new double[2]; //


	offmap = new double[2]; // map offsets


	// init some mapmaking variables
	ra_min  = 1000.0;
	ra_max  = -1000.0;
	dec_min = 1000.0;
	dec_max = -1000.0;

	offmap[0] = 0.0;
	offmap[1] = 0.0;






	//********************************************************************************
	//*************  find coordinates of pixels in the map
	//********************************************************************************

	printf("[%2.2i] Finding coordinates of pixels in the map\n",rank);

	coordsyst2 = coordsyst;
	if (coordsyst2 != 4){ // coordsyst never = 4 => debug mode => delete
		// ndet = number of channels
		for (idet=0;idet<ndet;idet++){

			// field = actual boloname
			field = bolonames[idet];
			// bolofield = boloname + bextension
			bolofield = field+bextension;

			//printf("%s\n",bolofield.c_str());

			if (cextension != "NOCALP")
				calfield  = field+cextension;
			if (fextension != "NOFLAG")
				flagfield = field+fextension;

			//read bolometer offsets
			read_bolo_offsets(field,file_offsets,scoffsets,offsets);

			// for each scan
			for (iframe=iframe_min;iframe<iframe_max;iframe++){

				//	cout << "[" << rank << "] " << idet << "/" << ndet << " " << iframe << "/" << ntotscan << endl;

				// read pointing files
				ns = nsamples[iframe];
				ff = fframes[iframe];

				// lis les données RA/DEC, Phi, et scerr (sky coordinates err, BLASPEC) pour ff et ns
				read_data_std(dirfile, ff, 0, ns, ra,   ra_field,  type); // type = 'd' 64-bit double
				read_data_std(dirfile, ff, 0, ns, dec,  dec_field, type);
				read_data_std(dirfile, ff, 0, ns, phi,  phi_field, type);
				read_data_std(dirfile, ff, 0, ns, scerr, scerr_field, type);
				read_data_std(dirfile, ff, 0, ns, flpoint, flpoint_field, 'c'); // flpoint = donnee vs time, take the data or not
				for (ii=0;ii<ns;ii++)
					if (isnan(ra[ii]) || isnan(dec[ii]) || isnan(phi[ii]))
						flpoint[ii] = 1; // sample is flagged, dont take this sample



				// find offset based on frame range
				correctFrameOffsets(nfoff,ff,offsets,foffsets,froffsets);


				// spheric coords to squared coords (tangent plane)
				sph_coord_to_sqrmap(pixdeg, ra, dec, phi, froffsets, ns, xx, yy, &nn, coordscorner,
						tancoord, tanpix, bfixc, radius, offmap, srccoord);

				// store ra/dec min and max (of the map) and for this processor
				if (coordscorner[0] < ra_min) ra_min = coordscorner[0];
				if (coordscorner[1] > ra_max) ra_max = coordscorner[1];
				if (coordscorner[2] < dec_min) dec_min = coordscorner[2];
				if (coordscorner[3] > dec_max) dec_max = coordscorner[3];

			}

		} //// end of idet loop

#ifdef USE_MPI
		MPI_Reduce(&ra_min,&gra_min,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
		MPI_Reduce(&ra_max,&gra_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
		MPI_Reduce(&dec_min,&gdec_min,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
		MPI_Reduce(&dec_max,&gdec_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&gra_min,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&gra_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&gdec_min,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&gdec_max,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

#else
		gra_min=ra_min;
		gra_max=ra_max;
		gdec_min=dec_min;
		gdec_max=dec_max;
#endif

		//set coordinates
		coordscorner[0] = gra_min; // store ra/dec min/max of the final map
		coordscorner[1] = gra_max;
		coordscorner[2] = gdec_min;
		coordscorner[3] = gdec_max;

		if (rank == 0) {
			printf("[%2.2i] ra  = [ %7.3f, %7.3f ] \n",rank, gra_min, gra_max );
			printf("[%2.2i] dec = [ %7.3f, %7.3f ] \n",rank, gdec_min, gdec_max);
		}

		/// just to set nn in order to compute map-making matrices and vectors
		sph_coord_to_sqrmap(pixdeg, ra, dec, phi, froffsets, ns, xx, yy, &nn, coordscorner,
				tancoord, tanpix, 1, radius, offmap, srccoord);



		sprintf(testfile,"%s%s%s%s",outdir.c_str(),"InfoPointing_for_Sanepic_",termin.c_str(),".txt");
		fp = fopen(testfile,"w");
		fprintf(fp,"%d\n",nn);
		fprintf(fp,"%d\n",coordsyst);
		fprintf(fp,"%lf\n",tanpix[0]);
		fprintf(fp,"%lf\n",tanpix[1]);
		fprintf(fp,"%lf\n",tancoord[0]);
		fprintf(fp,"%lf\n",tancoord[1]);
		fclose(fp);


	} else {
		// read those parameters from a file : -c = 4 option
		sprintf(testfile,"%s%s%s%s%d%s",outdir.c_str(),"InfoPointing_for_Sanepic_",termin.c_str(),"_",rank,".txt");
		if ((fp = fopen(testfile,"r")) == NULL){
			cerr << "File InfoPointing_for_sanepic... not found. Exiting" << endl;
			exit(1);
		}
		fscanf(fp,"%d\n",&nn);
		fscanf(fp,"%d\n",&coordsyst);
		fscanf(fp,"%lf\n",tanpix);
		fscanf(fp,"%lf\n",tanpix+1);
		fscanf(fp,"%lf\n",tancoord);
		fscanf(fp,"%lf\n",tancoord+1);
		fclose(fp);
	}







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



	long npixsrc = 0;
	indpsrc = new long[nn*nn];
	for (ii=0;ii<nn*nn;ii++){
		if (mask[ii] == 0){
			indpsrc[ii] = npixsrc;
			npixsrc += 1;
		} else {
			indpsrc[ii] = -1;
		}
	}
	long addnpix = ntotscan*npixsrc; // addnpix = number of pix to add in pixon = number of scans * number of pix in box crossing constraint removal

	//******************************************






	// map duplication factor
	int factdupl = 1;
	if (flgdupl) factdupl = 2; // -M =1, default 0 : if flagged data are put in a duplicated map


	//pixon indicates pixels that are seen
	pixon = new long[factdupl*nn*nn+2 + addnpix];   // last pixel is for flagged samples
	init1D_long(pixon,0,factdupl*nn*nn+2 + addnpix,0); // pixon[0->end] = 0, initialisation




	//**********************************************************************************
	//loop to get coordinates of pixels that are seen
	//**********************************************************************************

	/// loop again on detectors
	//for (idet=rank*ndet/size;idet<(rank+1)*ndet/size;idet++){
	for (idet=0;idet<ndet;idet++){


		field = bolonames[idet];
		//printf("%s\n",field.c_str());
		bolofield = field+bextension;
		calfield  = field+cextension;
		flagfield = field+fextension;


		// read bolometer offsets
		read_bolo_offsets(field,file_offsets,scoffsets,offsets);



		//loop to get coordinates of pixels that are seen
		for (iframe=0;iframe<ntotscan;iframe++){
			ns = nsamples[iframe];
			ff = fframes[iframe];

			read_data_std(dirfile, ff, 0, ns, ra,   ra_field,  type);
			read_data_std(dirfile, ff, 0, ns, dec,  dec_field, type);
			read_data_std(dirfile, ff, 0, ns, phi,  phi_field, type);
			read_data_std(dirfile, ff, 0, ns, scerr, scerr_field, type);
			read_data_std(dirfile, ff, 0, ns, flpoint, flpoint_field, 'c');
			for (ii=0;ii<ns;ii++)
				if (isnan(ra[ii]) || isnan(dec[ii]))
					flpoint[ii] = 1;
			if (fextension != "NOFLAG"){
				read_data_std(dirfile, ff, shift_data_to_point, ns, flag, flagfield,  'c');
			} else {
				for (ii=0;ii<ns;ii++)
					flag[ii] = 0;
			}



			// find offset based on frame range
			correctFrameOffsets(nfoff,ff,offsets,foffsets,froffsets);


			// return coordscorner, nn ,xx, yy, tancoord, tanpix
			sph_coord_to_sqrmap(pixdeg, ra, dec, phi, froffsets, ns, xx, yy, &nn, coordscorner,
					tancoord, tanpix, 1, radius, offmap, srccoord);




			// create flag field -----> conditions to flag data,
					// flag
					// scerr = pointing error : BLASPEC
					// flpoint = do we consider those data at time t
					// napod = number of samples to apodize
					// errarcsec critere de flag, default = 15.0, BLAST specific : pointing error threshold
					// NOFILLGAP = fill the gap ? default = yes
					// rejectsamples: rejected samples array, rejectsample value =0,1,2 or 3
			flag_conditions(flag,scerr,flpoint,ns,napod,xx,yy,nn,errarcsec,NOFILLGAP,rejectsamp);
			//flag_conditions(flag,scerr,flpoint,ns,napod,xx,yy,nn,errarcsec,NOFILLGAP,rejectsamp);

			// returns rejectsamp



			for (ii=0;ii<ns;ii++){
				//sample is rejected
				if (rejectsamp[ii] == 2){
					pixout = 1; // pixel is out of map
					pixon[factdupl*nn*nn + addnpix] += 1; // add one to a pixel outside the map
					samptopix[ii] = factdupl*nn*nn + addnpix; // the ii sample corresponds to a pixel outside the map


					printf("[%2.2i] PIXEL OUT, ii = %ld, xx = %i, yy = %i\n",rank, ii,xx[ii],yy[ii]);

				}
				// sample is not rejected
				if (rejectsamp[ii] == 0) {
					// if not in the  box constraint removal mask (defined by user)
					if (mask[yy[ii]*nn + xx[ii]] == 1){
						ll = yy[ii]*nn + xx[ii]; // ll=map indice
					} else {//if pixel is in the mask
						ll = indpsrc[yy[ii]*nn + xx[ii]]+factdupl*nn*nn+iframe*npixsrc;
					}
					pixon[ll] += 1;
					samptopix[ii] = ll; // the ii sample ii coresponds to pixel ll

					//printf("ll=%ld\n",ll);

				}
				// sample is flagged or scerr > scerr_threshold or flpoint =0
				if (rejectsamp[ii] == 1){
					if (flgdupl){ // if flagged pixels are in a duplicated map
						ll = yy[ii]*nn + xx[ii];
						pixon[nn*nn+ll] += 1;
						samptopix[ii] = nn*nn+ll;
					} else { // else every flagged sample is projected to the same pixel (outside the map)
						pixon[nn*nn+1 + addnpix] += 1;
						samptopix[ii] = nn*nn+1 + addnpix;
					}
				}
				// pixel dans la zone d'apodisation //a suppr
				if (rejectsamp[ii] == 3){
					pixon[factdupl*nn*nn+1 + addnpix] += 1;
					flagon = 1;
					samptopix[ii] = factdupl*nn*nn+1 + addnpix;
				}
			}



			//printf("pixon[nn*nn] = %d/n",pixon[nn*nn]);


			//if (rank == 0){

			sprintf(testfile,"%s%s%ld%s%ld%s%s%s",poutdir.c_str(),"samptopix_",iframe,"_",idet,"_",termin.c_str(),".bi");
			fp = fopen(testfile,"w");
			fwrite(samptopix,sizeof(long), ns, fp);
			fclose(fp);
			//}

		} // end of iframe loop

	}//end of idet loop






	//***********************************************************************************



	//************** init mapmaking variables *************//

	printf("[%2.2i] Init map making variables\n",rank);

	// useless
	//ns = nsamples[0];
	//for (ii=0;ii<ntotscan;ii++) if (nsamples[ii] > ns) ns = nsamples[ii];
	// a supprimer

	// pixel indices
	indpix = new long[factdupl*nn*nn+2 + addnpix];
	init1D_long(indpix,0,factdupl*nn*nn+2 + addnpix,-1);



	ll=0;
	for (ii=0;ii<factdupl*nn*nn+2 + addnpix;ii++){
		if (pixon[ii] != 0){
			indpix[ii] = ll; // pixel indice : 0 -> number of seen pixels
			ll++;
		}
	}
	npix = ll;  // npix = number of filled pixels


	delete [] pixon;

	// write in a file for conjugate gradient step // ajout Mat 02/06
	sprintf(testfile,"%s%s%s%s",outdir.c_str(),"Indpix_for_conj_grad_",termin.c_str(),".txt");
			fp = fopen(testfile,"w");
		//	long indpix_size = factdupl*nn*nn+2 + addnpix;
		//	fwrite(&indpix_size,sizeof(long),1,fp);
			fwrite(indpix,sizeof(long), factdupl*nn*nn+2 + addnpix, fp);
			fclose(fp);

	//  printf("[%2.2i] indpix[nn*nn] = %d\n",rank, indpix[nn*nn]);



	PNd = new double[npix];
	PNdtot = new double[npix];
	init1D_double(PNd,0,npix,0.0);
	init1D_double(PNdtot,0,npix,0.0);



	if (pixout)
		printf("THERE ARE SAMPLES OUTSIDE OF MAP LIMITS: ASSUMING CONSTANT SKY EMISSION FOR THOSE SAMPLES, THEY ARE PUT IN A SINGLE PIXEL\n");
	printf("[%2.2i] Total number of detectors : %d\t Total number of Scans : %d \n",rank,(int)ndet, (int) ntotscan);
	printf("[%2.2i] Size of the map : %d x %d\t Total Number of filled pixels : %d\n",rank, nn,nn, npix);


#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	//************************************************************************//
	//************************************************************************//
	//Pre-processing of the data
	//************************************************************************//
	//************************************************************************//



	printf("[%2.2i] Pre-processing of the data\n",rank);



	prefixe = "fdata";

	// loop over the scans
	for (iframe=iframe_min;iframe<iframe_max;iframe++){

		ns = nsamples[iframe]; // number of samples for this scan
		ff = fframes[iframe]; //first frame of this scan
		f_lppix = f_lp*double(ns)/fsamp; // knee freq of the filter in terms of samples in order to compute fft
		f_lppix_Nk = f_lp_Nk*double(ns)/fsamp; // noise PS threshold freq, in terms of samples
		//prefixe = "fdata"; => now outside the loop


		cout << "[" << rank << "] " << iframe << "/" << iframe_max ;
		cout << " ( -f " << ff << " -l " << ff+ns/20 << " )" << endl;

		// if there is correlation between detectors
		if (CORRon){
			// var double *S =  NULL
			    // indpix = pixel indice
			    // indpsrc = box crossing constraint removal pixel indice
			    // nn = size of map
			    // npix = number of pixels that are seen
			    // npixsrc = number of pixels in box CCRemoval
			    // ntotscan = total number of scans
			    // addnpix = number of added pixels in the map
			    // flgdupl = flaggued pixels are in a duplicate map : 1/0
			    // factdupl = duplication de la map : 2/1
			    // fillg =2 ????
			    // poutdir = outpout dir or current path (default)
			    // termin = output file suffix
			    // errarcsec = pointing error threshold
			    // dirfile = data directory
			    // scerr_field = "ERR" + pextension (_def for example)
			    // flpoint_field = "FLPOINTING"
			    // bolonames = bolo names array
			    // bextension = -B option : "_data" for example
			    // fextension = "NOFLAG" or -G option ("_flag" for example)
			    // cextension = "NOCALP" or -R option ("_calp" for example)
			    // shift_data_to_point (default 0), for subtracting a time offset to the data to match the pointing
			    // f_lppix = filter freq in term of sample
			    // ff = first frame number of this scan
			    // ns = number of sample for this scan
			    // napod = number of samples to apodize -A option
			    // ndet = bolo total number
			    // NORMLIN = baseline is remove from the data, default =0, option -L
			    // NOFILLGAP = fill the gap ? default yes => 0
			    // iframe = scan number : 0=> ntotscan
			write_ftrProcesdata(NULL,indpix,indpsrc,nn,npix,npixsrc,ntotscan,addnpix,flgdupl,factdupl,2,
					poutdir,termin,errarcsec,dirfile,scerr_field,flpoint_field,bolonames,
					bextension,fextension,cextension,shift_data_to_point,f_lppix,ff,ns,
					napod,ndet,NORMLIN,NOFILLGAP,iframe);// fillgaps + butterworth filter + fourier transform
			// "fdata_" files generation (fourier transform of the data)


			do_PtNd(PNd,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,f_lppix_Nk,
								fsamp,ff,ns,ndet,size_det,rank_det,indpix,nn,npix,iframe,NULL,NULL);

			// PNd = npix dimension, initialised to 0.0
				// extentnoiseSp_all = list of power spectrum file names (for each scan or same PS for all the scans)
				// noiseSppreffile = noise power spectrum file suffix = path
				// poutdir = outpout dir or current path (default)
				// prefixe = "fdata"; => prefixe de lecture/sauvegarde des données
				// termin = output file suffix
				// bolonames = bolo names array
				// f_lppix_Nk = freq threshold noise en terme de sample
				// fsamp = freq echantillonage des data
				// ff = n° premier sample du scan
				// ns = nombre de sample ds le scan
				// ndet = nombre de bolo
				// size = 1 // cf mpi
				// rank = 0 // cf mpi
				// indpix = pixel indice double[nn*nn]
				// nn = taille de la carte (1 coté)
				// npix = total number of filled pixels (pixel dont on a les data correspondantes)
				// iframe = indice du scan
				// *Mp = Null : la map ???
				// *Hits = Null

				// return Pnd = At N-1 d
		} else {

			do_PtNd_nocorr(PNd,extentnoiseSp_all,noiseSppreffile,poutdir,termin,errarcsec,dirfile,
					scerr_field,flpoint_field,bolonames,bextension,fextension,
					cextension,shift_data_to_point,f_lppix,f_lppix_Nk,fsamp,ntotscan,addnpix,
					flgdupl,factdupl,2,ff,ns,napod,ndet,size_det,rank_det,indpix,indpsrc,
					nn,npix,npixsrc,NORMLIN,NOFILLGAP,iframe,NULL);

		}


	} // end of iframe loop

	//  cout << "[" << rank << "] end of iframe loop" << endl;

#ifdef USE_MPI
	MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
	PNdtot=PNd;
#endif


	if (rank == 0){
		sprintf(testfile,"%s%s%s%s",outdir.c_str(),"PNdCorr_",termin.c_str(),".bi");
		fp = fopen(testfile,"w");
		fwrite(&npix,sizeof(int),1,fp); // ajout Mat 02/06
		fwrite(PNdtot,sizeof(double),npix,fp);
		fclose(fp);
	}

/* ----------------------------------------------------------------------------------------------*/
	// Noise Power spectra	estimation loop
	double *S;
	S = new double[npix];


	for (ii=0;ii<npix;ii++) S[ii] = 0.0;//PNd[ii];

	if (doInitPS == 1){

		if (MixMatfile != "NOFILE"){

			for (iframe=iframe_min;iframe<iframe_max;iframe++){
				ns = nsamples[iframe];
				ff = fframes[iframe];
				extentnoiseSp = extentnoiseSp_all[iframe];

				// estimate noise power spectra from data
				EstimPowerSpectra(fsamp,ns,ff,ndet,nn,npix,napod,iframe,flgdupl,factdupl,indpix,
						S,MixMatfile,bolonames,dirfile,bextension,fextension,cextension,
						shift_data_to_point,poutdir,termin,NORMLIN,NOFILLGAP,noiseSppreffile,
						extentnoiseSp,outdir);

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
				// noiseSppreffile = noise power spectrum file suffix = path
				// extentnoiseSp = noise file
				// outdir = output directory

			}
		}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		cout << "exit after first EstimPowerSpectra" << endl;
		exit(0);

	}

/* ---------------------------------------------------------------------------------------------*/

	/*
// Close MPI process


#ifdef USE_MPI
  MPI_Finalize();
#endif

  // clean up
    delete [] ra;
    delete [] dec;
    delete [] phi;
    delete [] scerr; //scerr_field needed
    delete [] xx;
    delete [] yy;
    delete [] flag;
    delete [] rejectsamp;
    delete [] samptopix;
    delete [] flpoint; // flpoint_field needed but not flpoint
    delete [] mask;
	delete [] offsets;
	delete [] froffests;
	delete [] offmap;
	delete [] scoffsets;
	delete [] foffsets;
	delete [] srccoord;
	delete [] coordscorner;

	delete [] extentnoiseSp_allorder;
	delete [] xxi_boxes;
	delete [] xxf_boxes;
	delete [] yyi_boxes;
	delete [] yyf_boxes;

	free(testfile);

	delete [] fframesorder; //will be needed
	delete [] nsamplesorder; //will be needed
	delete [] ruleorder; //will be needed
	delete [] frnum; //will be needed

	delete [] PNd; //needed
    delete [] PNdtot; //needed
    delete [] fframes; // needed
    delete [] nsamples; //needed
	delete [] tancoord; //needed
	delete [] tanpix; //needed
	delete [] indpix; //needed
	delete [] indpsrc; //needed
	delete [] bolonames; // needed
	delete [] extentnoiseSp_all; // needed

    return 0;
  }*/


	//  cout << "[" << rank << "] End of Init Loop" << endl;


	//******************************************************************//
	//******************************************************************//
	//**********************  End of init loop *************************//
	//******************************************************************//
	//******************************************************************//





	// liste des variables a donner :
	// npix, nn, indpix, projgaps, flagon, iframemax,

	// relancer le mpi et faire un best frame order

	/*includes*/
	/*
	//#include <iostream>
	//#include <iomanip>
	//#include <fstream>
	#include "todprocess.h"
	#include "map_making.h"
	#include "sane_io.h"
	#include <time.h>
	//#include <fftw3.h>
	#include <fcntl.h>
	#include <unistd.h>
	//#include <list>
	#include <stdio.h>
	#include <stdlib.h>

	 */

	/*Variables*/
	time_t t1;

	//************************************************************************//
	//************************************************************************//
	//main conjugate gradient loop
	//************************************************************************//
	//************************************************************************//

	printf("[%2.2i] Main Conjugate gradient loop\n",rank);

	int npixeff, iter, idupl;
	double var0, var_n, delta0, delta_n, delta_o, rtq, alpha, beta;

	//double *S
	double *PtNPmatS,  *PtNPmatStot, *r, *q, *qtot, *d, *Mp, *Mptot, *s;
	long *hits, *hitstot;


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


	long mi;
	double *map1d;
	map1d = new double[nn*nn];


	string fname;
	char iterchar[30];
	char iframechar[30];
	string iterstr, iframestr;



	for (idupl = 0;idupl<=flgdupl;idupl++){
		/*

		for (ii=0;ii<npix;ii++) S[ii] = 0.0;//PNd[ii];

		if (doInitPS == 1){

			if (MixMatfile != "NOFILE"){

				for (iframe=iframe_min;iframe<iframe_max;iframe++){
					ns = nsamples[iframe];
					ff = fframes[iframe];
					extentnoiseSp = extentnoiseSp_all[iframe];

					EstimPowerSpectra(fsamp,ns,ff,ndet,nn,npix,napod,iframe,flgdupl,factdupl,indpix,
							S,MixMatfile,bolonames,dirfile,bextension,fextension,cextension,
							shift_data_to_point,poutdir,termin,NORMLIN,NOFILLGAP,noiseSppreffile,
							extentnoiseSp,outdir);
				}
			}

#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif
			cout << "exit after first EstimPowerSpectra" << endl;
			exit(0);

		}
		 */


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

				do_PtNd(PtNPmatS,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,
						f_lppix_Nk,fsamp,ff,ns,ndet,size_det,rank_det,indpix,nn,npix,iframe,Mp,hits);
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

					do_PtNd(q,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,f_lppix_Nk,
							fsamp,ff,ns,ndet,size_det,rank_det,indpix,nn,npix,iframe,NULL,NULL);
				} else {

					do_PtNPS_nocorr(d,extentnoiseSp_all,noiseSppreffile,poutdir,termin,dirfile,bolonames,
							f_lppix_Nk,fsamp,flgdupl,factdupl,ff,ns,ndet,size_det,rank_det,indpix,
							nn,npix,iframe,q,NULL,NULL);
				}
			} // end of iframe loop




#ifdef USE_MPI
			MPI_Reduce(q,qtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
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

						do_PtNd(PtNPmatS,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,
								f_lppix_Nk,fsamp,ff,ns,ndet,size_det,rank_det,indpix,nn,npix,iframe,
								NULL,NULL);
					} else {

						do_PtNPS_nocorr(S,extentnoiseSp_all,noiseSppreffile,poutdir,termin,dirfile,bolonames,
								f_lppix_Nk,fsamp,flgdupl,factdupl,ff,ns,ndet,size_det,rank_det,
								indpix,nn,npix,iframe,PtNPmatS,NULL,NULL);
					}
				} // end of iframe loop



#ifdef USE_MPI
				MPI_Reduce(PtNPmatS,PtNPmatStot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
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
							bextension,fextension,cextension,shift_data_to_point,f_lppix,ff,ns,
							napod,ndet,NORMLIN,NOFILLGAP,iframe);

					do_PtNd(PNd,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,f_lppix_Nk,
												fsamp,ff,ns,ndet,size_det,rank_det,indpix,nn,npix,iframe,NULL,NULL);
				} else {

					do_PtNd_nocorr(PNd,extentnoiseSp_all,noiseSppreffile,poutdir,termin,errarcsec,dirfile,
												scerr_field,flpoint_field,bolonames,bextension,fextension,
												cextension,shift_data_to_point,f_lppix,f_lppix_Nk,fsamp,ntotscan,addnpix,
												flgdupl,factdupl,2,ff,ns,napod,ndet,size_det,rank_det,indpix,indpsrc,
												nn,npix,npixsrc,NORMLIN,NOFILLGAP,iframe,S);
				}
			} // end of iframe loop




#ifdef USE_MPI
			MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
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


	printf("%s\n",MixMatfile.c_str());

	if (MixMatfile != "NOFILE"){
		for (iframe=iframe_min;iframe<iframe_max;iframe++){
			ns = nsamples[iframe];
			ff = fframes[iframe];
			extentnoiseSp = extentnoiseSp_all[iframe];

			EstimPowerSpectra(fsamp,ns,ff,ndet,nn,npix,napod,iframe,flgdupl,factdupl,indpix,S,
					MixMatfile,bolonames,dirfile,bextension,fextension,cextension,shift_data_to_point,
					poutdir,termin,NORMLIN,NOFILLGAP,noiseSppreffile,extentnoiseSp,outdir);

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
		fp = fopen(testfile,"w");
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
	//delete [] ra;
	//delete [] dec;
	//delete [] phi;
	//delete [] scerr;
	//delete [] xx;
	//delete [] yy;
	//delete [] flag;
	//delete [] rejectsamp;
	//delete [] samptopix;
	//delete [] flpoint;
	//delete [] mask;




#ifdef USE_MPI
	MPI_Finalize();
#endif



	return 0;
}









