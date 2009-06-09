#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <cstdlib>

#include <fcntl.h>
#include <unistd.h>
#include <list>
#include <vector>
#include <stdio.h>

#include "todprocess.h"
#include "map_making.h"
// #include "sane_io.h"
//#include "estimPS.h"
#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"
#include "mpi_architecture_builder.h"

#include "boloIO.h"

extern "C" {
#include <fftw3.h>
}


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
	cerr << "-R <sampling frequency> detectors sampling frequency. Required" << endl;
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
	cerr << "-r <precompute PNd>    Keyword specifying if PNd is precomputed and read from disk" << endl; //add this option to work
	cerr << "-g <project gaps>      Keyword specifying if gaps are projected to a pixel in the map, if so gap filling of noise only is performed iteratively. Default is 0" << endl;
	cerr << "-M <map flagged data>  Keyword specifying if flagged data are put in a separate map, default is 0" << endl;
	cerr << "-s <time offset>       Keyword for subtracting a time offset to the data to match the pointing "<< endl;
	cerr << "-E <no corr>           Set this keyword to 0 if correlations are not included in the analysis" << endl;
	//cerr << "-I <parallel scheme>   Set this keyword to 1 in order to distribute different field visit (or frame ranges) to different processors, instead of different detector data chuncks as it is by default" << endl;
	//cerr << "-j <write unconverged maps>  The value specified by this keyword indicates the period in iterations to which the data are written to disk. Default is 10. Set this keyword to zero if you don't want intermediate maps to be written." << endl;
	cerr << "-a <noise estim>       Optional. Enter filename containing the mixing matrix of noise components. If set, the noise power spectra files for each field are computed after map-making. Default is no re-estimation" << endl;
	cerr << "-i <Noise PS etimation> Optional. Computes power spectra estimation from the data and saves it in a file. Default is no estimation." << endl;
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
	size = 1;
	rank = 0;
	cout << "Mpi is not used for this step" << endl;
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
	double fsamp = 0.0;// 25.0; // sampling frequency : BLAST Specific
	double errarcsec = 15.0; // rejection criteria : scerr[ii] > errarcsec, sample is rejected
	// source error

	long ii, ll, iframe, idet, ib; // loop indices
	long iframe_min, iframe_max;
	int flagon = 0; // if rejectsample [ii]==3, flagon=1
	//int iterw = 10; // period in iterations to which the data are written to disk, 0 = no intermediate map to be written
	bool bfixc = 0; // indicates that 4 corners are given for the cross corelation removal box
	bool pixout = 0; // indicates that at least one pixel has been flagged and is out
	bool NORMLIN = 0; // baseline is removed from the data, NORMLIN = 1 else 0
	bool NOFILLGAP = 0; // fill the gap ? default is YES
	bool PND_ready = 0; // PNd precomputed ? read on disk if =1
	bool flgdupl = 0; // 1 if flagged data are put in a separate map
	bool CORRon = 1; // correlation included in the analysis (=1), else 0, default 0
	//bool parallel_frames = 0; // a parallel scheme is used : mpi has been launched

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
	// moved array of strings to vector of strings...
	//string *bolonames; // channel list -> bolonames array : considered bolometers names
	string *extentnoiseSp_all; // ((list -> string*))
	string *extentnoiseSp_allorder;
	string bolofield; // bolofield = boloname + bextension
	//string calfield; // calfield  = field+cextension;
	string flagfield; // flagfield = field+fextension;
	string dirfile; // data directory
	string outdir; // output directory
	string poutdir; // current path (pPath) or output dir (outdir)
	string bextension; // bolometer field extension
	//string cextension = "NOCALP"; // needed for calfield calculation !
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

	std::vector<string> bolonames; // bolometer list

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
		case 'J':
			f_lp_Nk = atof(optarg);
			break;
		case 'A':
			napod = atoi(optarg);
			break;
		case 'y':
			bolonames.push_back(optarg);
			break;
		case 'C':
			read_bolofile(optarg, bolonames);
			//      cerr << "num ch: "<< channel.size() << endl;
			break;
		case 'o':
			outdir = optarg;
			break;
		case 'B':
			bextension = optarg;
			break;
		case 'R':
			//cextension = optarg;
			fsamp = atof(optarg);
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
			NORMLIN = 1;
			break;
		case 'g':
			projgaps = 1;
			break;
		case 'r':
			PND_ready = 1;
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
		case 'I': // no more used
			//parallel_frames = 1;
			break;
		case 'j':
		//	iterw = atoi(optarg);
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
	if (bolonames.size() == 0) {
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

	cout << "Sampling frequency : " << fsamp << endl;

	ntotscan = ff_in.size();
	ndet = bolonames.size();

	nnf = extentnoiseSp_list.size();
	if (nnf != 1 && nnf != ntotscan){
		cerr << "ERROR: There should be one noise power spectrum file per scan, or a single one for all the scans. Check -K options" << endl;
		exit(1);
	}
	//  printf("%d\n",nnf);


//	TODO: read needed data from position/parallelization stuff....


	PNd = new double[npix];
	PNdtot = new double[npix];
	init1D_double(PNd,0,npix,0.0);
	init1D_double(PNdtot,0,npix,0.0);


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

		printf("[%2.2i] iframe : %ld/%ld",rank,iframe,iframe_max);
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
					bextension,fextension/*,cextension*/,shift_data_to_point,f_lppix,ff,ns,
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
					/*cextension,*/shift_data_to_point,f_lppix,f_lppix_Nk,fsamp,ntotscan,addnpix,
					flgdupl,factdupl,2,ff,ns,napod,ndet,size_det,rank_det,indpix,indpsrc,
					nn,npix,npixsrc,NORMLIN,NOFILLGAP,iframe,NULL);

		}


	} // end of iframe loop

	//  cout << "[" << rank << "] end of iframe loop" << endl;

	printf("[%2.2i] End of Pre-Processing\n",rank);

#ifdef USE_MPI
	MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
	PNdtot=PNd;
#endif
/*
for (ii=0;ii<20;ii++)
	cout << PNdtot[ii] << " ";
cout << endl;*/

	if (rank == 0){
		sprintf(testfile,"%s%s%s%s",outdir.c_str(),"PNdCorr_",termin.c_str(),".bi");
		if((fp = fopen(testfile,"w"))!=NULL){
		fwrite(&npix,sizeof(int),1,fp); // ajout Mat 02/06
		fwrite(PNdtot,sizeof(double),npix,fp);
		fclose(fp);
		}else{
			cerr << "Error : cannot open " << testfile << endl;
			exit(1);
		}
	}


	/*
	 *
	 * Noise Power spectra estimation has been moved to another program
	 *
	 *

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
						S,MixMatfile,bolonames,dirfile,bextension,fextension, //cextension,
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
		cout << "Exit after first EstimPowerSpectra" << endl;
		exit(0);

	}

*/


/* ---------------------------------------------------------------------------------------------*/


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
	delete [] froffsets;
	delete [] offmap;
	//delete [] scoffsets;
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
//	delete [] bolonames; // needed
	delete [] extentnoiseSp_all; // needed

    return 0;
  }


	//  cout << "[" << rank << "] End of Init Loop" << endl;


	//******************************************************************//
	//******************************************************************//
	//**********************  End of init loop *************************//
	//******************************************************************//
	//******************************************************************//



