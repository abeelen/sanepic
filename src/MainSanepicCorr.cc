#include <iostream>
#include <fstream>
#include "todprocess.h"
#include "map_making.h"
#include <time.h>
#include <fftw3.h>
#include <getdata.h>
//#include "/global/software/pgi-6.0/linux86/6.0/include/mpi.h"
#include <fcntl.h>
#include <unistd.h>
#include <list>
#include <stdio.h>
#include <stdlib.h>



extern "C" {
#include <fitsio.h>
}



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
  cerr << "-m <padding interval>  number of samples extrapolated before (and after) data" << endl;
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
  cerr << "-M <map flagged data>   Keyword specififying if flagged data are put in a separate map, default is 0" << endl;
  cerr << "-s <time offset>       Keyword for subtracting a time offset to the data to match the pointing "<< endl;
  cerr << "-E <no corr>           Set this keyword to 0 if correlations are not included in the analysis" << endl;
  cerr << "-j <write unconverged maps>  The value specified by this keyword indicates the period in iterations to which the data are written to disk. Default is 10. Set this keyword to zero if you don't want intermediate maps to be written." << endl;
  cerr << "-a <noise estim>       Optional. Enter filename containing the mixing matrix of noise components. If set, the noise power spectra files for each field are computed after map-making. Default is no re-estimation" << endl;
  exit(1);
}




void print_fits_error(int status){
  if(status){
    fits_report_error(stderr, status); /* print error report */
    exit(status);    /* terminate the program, returning error status */
  }
  return;
}



void write_fits(string fname, double pixsize, long nx, long ny,
		double *tancoord, double *tanpix, int coordsyst, char dtype, void *data)
{
  // all angles in degrees
  // coordcenter is a 2-element array containing RA/DEC (or l/b) of the central pixel

  fitsfile *fp;
  int fits_status = 0;

  long naxis = 2;           // number of dimensions
  long naxes[] = {nx, ny};  // size of dimensions
  long fpixel[] = {1, 1};   // index for write_pix
  long ndata = nx * ny;     // number of data points

  double dtmp;
  char *strx, *stry;


  // create fits file
  if ( fits_create_file(&fp, fname.c_str(), &fits_status) )
    print_fits_error(fits_status);

  // create fits image (switch on data type)
  switch (dtype) {
  case 'd':    // double
    if ( fits_create_img(fp, DOUBLE_IMG, naxis, naxes, &fits_status) )
      print_fits_error(fits_status);
    break;
  case 'l':    // long
    if ( fits_create_img(fp, LONG_IMG, naxis, naxes, &fits_status) )
      print_fits_error(fits_status);
    break;
  default:
    cerr << "write_fits: data type '" << dtype << "' not supported. Exiting.\n";
    exit(1);
  }

  // write date to file
  if ( fits_write_date(fp, &fits_status) )
    print_fits_error(fits_status);

  // write map parameters (keywords)
  if ( fits_write_key(fp, TLONG, "NROW", &nx, "Number of rows", &fits_status) )
    print_fits_error(fits_status);

  if ( fits_write_key(fp, TLONG, "NCOL", &ny, "Number of columns", &fits_status) )
    print_fits_error(fits_status);

  if ( fits_write_key(fp, TDOUBLE, "PIXSIZE", &pixsize, "Size of pixels (deg)", &fits_status) )
    print_fits_error(fits_status);

  if ( fits_write_comment(fp, "Galactic coordinates",  &fits_status) )
    print_fits_error(fits_status);

  dtmp = (tanpix[0]); // 0-based index to 1
  if ( fits_write_key(fp, TDOUBLE, "CRPIX1", &dtmp, "X PIXEL OF TANGENT POINT", &fits_status) )
    print_fits_error(fits_status);

  dtmp = (tanpix[1]); // 0-based index to 1
  if ( fits_write_key(fp, TDOUBLE, "CRPIX2", &dtmp, "Y PIXEL OF TANGENT POINT", &fits_status) )
    print_fits_error(fits_status);

  dtmp = -pixsize;
  if ( fits_write_key(fp, TDOUBLE, "CDELT1", &dtmp, "COORD VALUE INCR DEG/PIXEL AT ORIGIN ON LINE AXIS",
		      &fits_status) )
    print_fits_error(fits_status);

  if ( fits_write_key(fp, TDOUBLE, "CDELT2", &pixsize, "COORD VALUE INCR DEG/PIXEL AT ORIGIN ON LINE AXIS",
		      &fits_status) )
    print_fits_error(fits_status);

  if (coordsyst == 2){
    if ( fits_write_key(fp, TDOUBLE, "CRVAL1", tancoord, "GLON AT TANGENT POINT (DEG)", &fits_status) )
      print_fits_error(fits_status);

    if ( fits_write_key(fp, TDOUBLE, "CRVAL2", tancoord+1, "GLAT AT TANGENT POINT (DEG)", &fits_status) )
      print_fits_error(fits_status);

    strx = "GLON-TAN";
    stry = "GLAT-TAN";

  } else {
    if ( fits_write_key(fp, TDOUBLE, "CRVAL1", tancoord, "RA AT TANGENT POINT (DEG)", &fits_status) )
      print_fits_error(fits_status);

    if ( fits_write_key(fp, TDOUBLE, "CRVAL2", tancoord+1, "DEC AT TANGENT POINT (DEG)", &fits_status) )
      print_fits_error(fits_status);

    strx = "RA---TAN";
    stry = "DEC--TAN";

  }

  if ( fits_write_key(fp, TSTRING, "CTYPE1", strx, "TANGENT PLANE PROJECTION", &fits_status) )
    print_fits_error(fits_status);

  if ( fits_write_key(fp, TSTRING, "CTYPE2", stry, "TANGENT PLANE PROJECTION", &fits_status) )
    print_fits_error(fits_status);


  // write map data
  switch (dtype) {
  case 'd':    // double
    if ( fits_write_pix(fp, TDOUBLE, fpixel, ndata, (double*) data, &fits_status) )
      print_fits_error(fits_status);
    break;
  case 'l':    // long
    if ( fits_write_pix(fp, TLONG, fpixel, ndata, (long*) data, &fits_status) )
      print_fits_error(fits_status);
    break;
  }

  // close file
  if(fits_close_file(fp, &fits_status))
    print_fits_error(fits_status);

}



void write_vector(char *filename, void *data, int typesize, long nn) {

  FILE *fp;

  fp = fopen(filename,"w");
  fwrite(data,typesize, nn, fp);
  fclose(fp);

}


void read_vector(char *filename, void *data, int typesize, long nn) {

  FILE *fp;

  fp = fopen(filename,"r");
  fread(data,typesize, nn, fp);
  fclose(fp);

}



void read_bolofile(string fname, list<string>& bolos) {
  char buff[256];
  string line;

  ifstream BOLO (fname.c_str());
  if (! BOLO.is_open()) {
    cerr << "Error opening bolometer file '" << fname << "'. Exiting.\n";
    exit(1);
  }

  while (! BOLO.eof()) {
    BOLO.getline(buff,255);
    line = buff;

    line.erase(0, line.find_first_not_of(" \t"));       // remove leading white space
    if (line.empty() || line[0] == '#') continue;       // skip if empty or commented
    line = line.substr(0, line.find_first_of(" \t"));   // pick out first word

    bolos.push_back(line);
  }

  BOLO.close();
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




void read_bolo_offsets(string field, string file_BoloOffsets, float *scoffsets, double *offsets){


  double lel, xel;
  long temp1, temp2, temp3;
  int nobolo = 1;

  char boloname[100];
  FILE *fp;


  if ((fp = fopen(file_BoloOffsets.c_str(),"r")) == NULL){
    cerr << "ERROR: Can't find offset file. Exiting. \n";
    exit(1);
  }
  while (fscanf(fp, "%s%ld%ld%ld%lf%lf\n", boloname, &temp1, &temp2, &temp3, &lel, &xel) != EOF) {
    if (field == boloname) {
      nobolo = 0;
      if (temp3 == 250){
	offsets[0] = xel/60.0/60.0 - scoffsets[1];
	offsets[1] = lel/60.0/60.0 + scoffsets[0];
      }
      if (temp3 == 350){
	offsets[0] = xel/60.0/60.0 - scoffsets[3];
	offsets[1] = lel/60.0/60.0 + scoffsets[2];
      }
      if (temp3 == 500){
	offsets[0] = xel/60.0/60.0 - scoffsets[5];
	offsets[1] = lel/60.0/60.0 + scoffsets[4];
      }
    }
  }
  fclose (fp);


  if (nobolo){
    cerr << "Bolometer name not found in offset list" << endl;
    exit(1);
  }


}







void do_PtNd(double *PNd, string *extentnoiseSp_all, string noiseSppreffile,
	     string dir, string prefixe, string termin, string *bolonames,
	     double f_lppix, double fsamp, long ff, long ns, long marge, long ndet, int size,
	     int rank, long *indpix, long nn, long npix, long iframe, double *Mp, long *hits);




void write_ftrProcesdata(double *S, long *indpix, long *indpsrc, int nn, long npix,
			 long npixsrc, long ntotscan, long addnpix, bool flgdupl, int factdupl,
			 int fillg, string dir, string termin, double errarcsec, string dirfile,
			 string scerr_field, string flpoint_field, string *bolonames,
			 string bextension, string fextension, string cextension,
			 int shift_data_to_point, double f_lppix, long ff, long ns,
			 long marge, long napod, long ndet, bool NORMLIN, bool NOFILLGAP,
			 long iframe);



void write_tfAS(double *S, long *indpix, int nn, long npix, bool flgdupl, int factdupl,
		string dir, string termin, long ff, long ns, long marge, long ndet, long iframe);




void do_PtNd_nocorr(double *PNd, string *extentnoiseSp_all, string noiseSppreffile,
		    string dir, string termin, double errarcsec, string dirfile,
		    string scerr_field, string flpoint_field, string *bolonames,
		    string bextension, string fextension, string cextension,
		    int shift_data_to_point, double f_lppix, double f_lppix_Nk,
		    double fsamp, long ntotscan, long addnpix, bool flgdupl, int factdupl,
		    int fillg, long ff, long ns, long marge, long napod, long ndet,
		    int size, int rank, long *indpix, long *indpsrc, long nn, long npix,
		    long npixsrc, bool NORMLIN, bool NOFILLGAP, long iframe, double *S);




void do_PtNPS_nocorr(double *S, string *extentnoiseSp_all, string noiseSppreffile, string dir,
		     string termin, string dirfile, string *bolonames, double f_lppix,
		     double fsamp, bool flgdupl, int factdupl, long ff, long ns, long marge,
		     long ndet, int size, int rank, long *indpix, long nn, long npix,
		     long iframe, double *PtNPmatS, double *Mp, long *hits);




void EstimPowerSpectra(double fsamp, long ns, long ff, long ndet, int nn, long npix, long napod,
		       long marge, long iframe, bool flgdupl, int factdupl, long *indpix,
		       double *S, string MixMatfile, string *bolonames, string dirfile, string bextension,
		       string fextension, string cextension, int shift_data_to_point, string dir,
		       string termin, bool NORMLIN, bool NOFILLGAP, string noiseSppreffile,
		       string extentnoiseSp, string outdirSpN);


double fdsf(double **Rellexp, double *w, double **A, double **P, double **N, long ndet, long ncomp, long nbins);

void rescaleAP(double **A, double **P, long ndet, long ncomp, long nbins);


//**********************************************************************************//
//**********************************************************************************//
//*************************** Beginning of main program ****************************//
//**********************************************************************************//
//**********************************************************************************//




int main(int argc, char *argv[])
{



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
  long napod = 0;
  long marge = 0;
  double fsamp = 100.0;
  double errarcsec = 15.0;


  long ii, jj, ll, iframe, idet, ib;
  int flagon = 0;
  int iterw = 10;
  bool bfixc = 0;
  bool pixout = 0;
  bool NORMLIN = 0;
  bool NOFILLGAP = 0;
  bool PND_ready = 0;
  bool flgdupl = 0;
  bool CORRon = 1;

  //set coordinate system
  double *srccoord, *coordscorner;
  srccoord = new double[2];
  coordscorner = new double[4];
  srccoord[0] = -1000;
  srccoord[1] = -1000;
  double radius = -1.0;


  // data parameters
  long *fframes  ;
  long *nsamples ;

  long *xxi_boxes, *xxf_boxes;
  long *yyi_boxes, *yyf_boxes;

  long ntotscan;
  long ndet;
  int nnf;

  // map making parameters
  double pixdeg;


  int nn, npix;
  double ra_min, ra_max, dec_min, dec_max;
  double *offsets, *froffsets, *offmap;
  double *tancoord;
  double *tanpix;

  //internal data params
  long ns, ff;
  double f_lp, f_lp_Nk, f_lppix, f_lppix_Nk;


  FILE *fp;

  char testfile[100];



  char type='d';
  double *ra, *dec, *phi, *scerr;
  unsigned char *flag, *flpoint, *rejectsamp, *mask;
  double *PNd;
  long *indpix, *indpsrc;

  int *xx, *yy;
  long *pixon;
  long *samptopix;


  string field;
  string *bolonames;
  string *extentnoiseSp_all;
  string bolofield;
  string calfield;
  string flagfield;
  string dirfile;
  string outdir;
  string poutdir;
  string bextension;
  string cextension = "NOCALP";
  string fextension = "NOFLAG";
  string pextension;
  string file_offsets;
  string file_frame_offsets = "NOOFFS";
  string termin;
  string noiseSppreffile;
  string extentnoiseSp;
  string prefixe;

  string MixMatfile = "NOFILE";

  /* DEFAULT PARAMETERS */
  int coordsyst = 1; /// Default is RA/DEC
  int coordsyst2 = 1;


   /* COMMAND-LINE PARAMETER PROCESSING */

  int retval;
  int tmpcount = 0;
  int tmpcount2 = 0;
  int temp;

  list<string> channel;
  list<long> ff_in, nf_in, xxi_in, xxf_in, yyi_in, yyf_in;
  list<double> fcut_in;
  list<string> extentnoiseSp_list;


  time_t t1, t2, t3, t4, t5, dt;



  f_lp = 0.0;
  f_lp_Nk = 0.0;
  pixdeg = -1.0;


  // Parse command line options
  while ( (retval = getopt(argc, argv, "F:f:l:n:y:C:H:J:o:O:B:R:G:P:S:e:p:A:m:k:K:t:T:u:U:v:V:c:N:L:g:r:M:x:X:z:Z:s:E:j:a:")) != -1) {
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
    case 'm':
      marge = atoi(optarg);
      break;
    case 'y':
      channel.push_back(optarg);
      break;
    case 'C':
      read_bolofile(optarg, channel);
      cerr << "num ch: "<< channel.size() << endl;
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
    case 'j':
      iterw = atoi(optarg);
      break;
    case 'a':
      MixMatfile = optarg;
      break;
    default:
      cerr << "Option '" << (char)retval << "' not valid. Exiting.\n\n";
      usage(argv[0]);
    }
  }


  if (CORRon) printf("CORRELATIONS BETWEEN DETECTORS INCLUDED\n");
  if (!CORRon) printf("NO CORRELATIONS BETWEEN DETECTORS INCLUDED\n");


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
    printf("Data are apodized\n");
  } else {
    printf("Data are not apodized\n");
  }
  if (marge){
    printf("Data are extrapolated using a linear predictor\n");
  } else {
    printf("Data are not extrapolated outside the edges\n");
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
  printf("%d\n",nnf);


  // convert lists to regular arrays
  fframes  = new long[ntotscan];
  nsamples = new long[ntotscan];
  bolonames = new string [ndet];
  extentnoiseSp_all = new string[ntotscan];


  list2array(ff_in, fframes);
  list2array(nf_in, nsamples);
  list2array(channel, bolonames);
  list2array(extentnoiseSp_list,extentnoiseSp_all);
  if (nnf == 1 && ntotscan > 1)
    for (ii=1;ii<ntotscan;ii++)
      extentnoiseSp_all[ii] = extentnoiseSp_all[0];


  printf("xxi_in.size() = %d\n",xxi_in.size());


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



  for (int ii=0; ii<ntotscan; ii++) {
    nsamples[ii] *= 20;      // convert nframes to nsamples
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


  string ra_field;
  string dec_field;
  string phi_field;
  string scerr_field = "ERR"+pextension;
  string flpoint_field = "FLPOINTING";

  if (coordsyst == 2){
    ra_field = "L"+pextension;
    dec_field = "B"+pextension;
    phi_field = "PHIG"+pextension;
    printf("Coordinate system: Galactic\n");
  }else{
    ra_field = "RA"+pextension;
    dec_field = "DEC"+pextension;
    phi_field = "PHI"+pextension;
    if (coordsyst == 3){
      printf("Map in Telescope coordinates. Reference coordinate system is RA/DEC (J2000)\n");
    } else {
      printf("Coordinate system: RA/DEC (J2000)\n");
    }
  }


  if (NORMLIN)
    printf("NO BASELINE REMOVED\n");


  if (projgaps)
    printf("Flaged data are binned. iterative solution to fill gaps with noise only.\n");


  // map offsets
  float scoffsets[6];
  int nfoff;
  foffset *foffsets;
  if (file_frame_offsets != "NOOFFS")
    foffsets = read_mapoffsets(file_frame_offsets, scoffsets, &nfoff);
  else {
    printf("No offsets between visits\n");
    nfoff = 2;
    ff = fframes[0];
    for (ii=0;ii<ntotscan;ii++) if (fframes[ii] > ff) ff = fframes[ii];
    foffsets = new foffset [2];
    (foffsets[0]).frame = 0;
    (foffsets[0]).pitch = 0.0;
    (foffsets[0]).yaw   = 0.0;
    (foffsets[1]).frame = ff+1;
    (foffsets[1]).pitch = 0.0;
    (foffsets[1]).yaw   = 0.0;
    for (ii=0;ii<6;ii++)
      scoffsets[ii] = 0.0;
  }



  if (pPath != NULL){
    poutdir = pPath;
  } else {
    poutdir = outdir;
    //poutdir = "/scratch/";
  }
  printf("Data written in %s\n",poutdir.c_str());



  /* END PARAMETER PROCESSING */





  /********** Alocate memory ***********/
  ns = nsamples[0];
  for (ii=0;ii<ntotscan;ii++) if (nsamples[ii] > ns) ns = nsamples[ii];

  ra = new double[2*ns];
  dec = new double[2*ns];
  phi = new double[2*ns];
  scerr = new double[2*ns];
  xx = new int[2*ns];
  yy = new int[2*ns];
  samptopix = new long[2*ns];
  flag = new unsigned char[2*ns];
  rejectsamp = new unsigned char[2*ns];
  flpoint = new unsigned char[2*ns];
  tancoord = new double[2];
  tanpix = new double[2];
  offsets = new double[2];
  froffsets = new double[2];


  offmap = new double[2];


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

  coordsyst2 = coordsyst;
  if (coordsyst2 != 4){

    for (idet=0;idet<ndet;idet++){


      field = bolonames[idet];
      bolofield = field+bextension;

      //printf("%s\n",bolofield.c_str());

      if (cextension != "NOCALP")
	calfield  = field+cextension;
      if (fextension != "NOFLAG")
	flagfield = field+fextension;

      //read bolometer offsets
      read_bolo_offsets(field,file_offsets,scoffsets,offsets);


      for (iframe=0;iframe<ntotscan;iframe++){
	// read pointing files
	ns = nsamples[iframe];
	ff = fframes[iframe];

	read_data(dirfile, ff, 0, ns, ra,   ra_field,  type);
	read_data(dirfile, ff, 0, ns, dec,  dec_field, type);
	read_data(dirfile, ff, 0, ns, phi,  phi_field, type);
	read_data(dirfile, ff, 0, ns, scerr, scerr_field, type);
	read_data(dirfile, ff, 0, ns, flpoint, flpoint_field, 'c');
	for (ii=0;ii<ns;ii++)
	  if (isnan(ra[ii]) || isnan(dec[ii]) || isnan(phi[ii]))
	    flpoint[ii] = 1;



	// find offset based on frame range
	correctFrameOffsets(nfoff,ff,offsets,foffsets,froffsets);



	sph_coord_to_sqrmap(pixdeg, ra, dec, phi, froffsets, ns, xx, yy, &nn, coordscorner,
			    tancoord, tanpix, bfixc, radius, offmap, srccoord);
	if (coordscorner[0] < ra_min) ra_min = coordscorner[0];
	if (coordscorner[1] > ra_max) ra_max = coordscorner[1];
	if (coordscorner[2] < dec_min) dec_min = coordscorner[2];
	if (coordscorner[3] > dec_max) dec_max = coordscorner[3];

      }

    }//// end of idet loop


    //set coordinates
    coordscorner[0] = ra_min;
    coordscorner[1] = ra_max;
    coordscorner[2] = dec_min;
    coordscorner[3] = dec_max;

    /// just to set nn in order to compute map-making matrices and vectors
    sph_coord_to_sqrmap(pixdeg, ra, dec, phi, froffsets, ns, xx, yy, &nn, coordscorner,
			tancoord, tanpix, 1, radius, offmap, srccoord);





    printf("inside main: ra_min  = %lf\n",ra_min );
    printf("inside main: ra_max  = %lf\n",ra_max );
    printf("inside main: dec_min = %lf\n",dec_min);
    printf("inside main: dec_max = %lf\n",dec_max);



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

    sprintf(testfile,"%s%s%s%s",outdir.c_str(),"InfoPointing_for_Sanepic_",termin.c_str(),".txt");
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
    for (ib = 0;ib<long(xxi_in.size()); ib++){
      for (ii=xxi_boxes[ib];ii<xxf_boxes[ib];ii++)
	for (ll=yyi_boxes[ib];ll<yyf_boxes[ib];ll++)
	  mask[ll*nn + ii] = 0;
    }
  }

  //for (ii=361;ii<418;ii++)
  // for (jj=598;jj<661;jj++)
  //   mask[jj*nn + ii] = 0;
  //for (ii=328;ii<371;ii++)
  // for (jj=663;jj<718;jj++)
  //   mask[jj*nn + ii] = 0;


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
  long addnpix = ntotscan*npixsrc;

  //******************************************







    int factdupl = 1;
    if (flgdupl) factdupl = 2;


    //pixon indicates pixels that are seen
    pixon = new long[factdupl*nn*nn+2 + addnpix];   // last pixel is for flagged samples
    init1D_long(pixon,0,factdupl*nn*nn+2 + addnpix,0);




    //**********************************************************************************
    //loop to get coordinates of pixels that are seen
    //**********************************************************************************

    /// loop again on detectors
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

	read_data(dirfile, ff, 0, ns, ra,   ra_field,  type);
	read_data(dirfile, ff, 0, ns, dec,  dec_field, type);
	read_data(dirfile, ff, 0, ns, phi,  phi_field, type);
	read_data(dirfile, ff, 0, ns, scerr, scerr_field, type);
	read_data(dirfile, ff, 0, ns, flpoint, flpoint_field, 'c');
	for (ii=0;ii<ns;ii++)
	  if (isnan(ra[ii]) || isnan(dec[ii]))
	    flpoint[ii] = 1;
	if (fextension != "NOFLAG"){
	  read_data(dirfile, ff, shift_data_to_point, ns, flag, flagfield,  'c');
	} else {
	  for (ii=0;ii<ns;ii++)
	    flag[ii] = 0;
	}



	// find offset based on frame range
	correctFrameOffsets(nfoff,ff,offsets,foffsets,froffsets);



	sph_coord_to_sqrmap(pixdeg, ra, dec, phi, froffsets, ns, xx, yy, &nn, coordscorner,
			    tancoord, tanpix, 1, radius, offmap, srccoord);




	//create flag field -----> conditions to flag data
	flag_conditions(flag,scerr,flpoint,ns,napod,marge,xx,yy,nn,errarcsec,NOFILLGAP,rejectsamp);




	for (ii=0;ii<ns;ii++){
	  if (rejectsamp[ii] == 2){
	    pixout = 1;
	    pixon[factdupl*nn*nn + addnpix] += 1;
	    samptopix[ii] = factdupl*nn*nn + addnpix;


	    //printf("PIXEL OUT, ii = %d, xx = %d, yy = %ld\n",ii,xx[ii],yy[ii]);

	  }
	  if (rejectsamp[ii] == 0) {
	    if (mask[yy[ii]*nn + xx[ii]] == 1){
	      ll = yy[ii]*nn + xx[ii];
	    } else {
	      ll = indpsrc[yy[ii]*nn + xx[ii]]+factdupl*nn*nn+iframe*npixsrc;
	    }
	    pixon[ll] += 1;
	    samptopix[ii] = ll;
	  }
	  if (rejectsamp[ii] == 1){
	    if (flgdupl){
	      ll = yy[ii]*nn + xx[ii];
	      pixon[nn*nn+ll] += 1;
	      samptopix[ii] = nn*nn+ll;
	    } else {
	      pixon[nn*nn+1 + addnpix] += 1;
	      samptopix[ii] = nn*nn+1 + addnpix;
	    }
	  }
	  if (rejectsamp[ii] == 3){
	    pixon[factdupl*nn*nn+1 + addnpix] += 1;
	    flagon = 1;
	    samptopix[ii] = factdupl*nn*nn+1 + addnpix;
	  }
	}

	pixon[factdupl*nn*nn + addnpix] += 2*marge;


	//printf("pixon[nn*nn] = %d/n",pixon[nn*nn]);


	//if (rank == 0){

	  sprintf(testfile,"%s%s%ld%s%ld%s%s%s",poutdir.c_str(),"samptopix_",iframe,"_",idet,"_",termin.c_str(),".bi");
	  fp = fopen(testfile,"w");
	  fwrite(samptopix,sizeof(long), ns, fp);
	  fclose(fp);
	  //}

      }

    }//end of idet loop






  //***********************************************************************************



  //************** init mapmaking variables *************//

  ns = nsamples[0];
  for (ii=0;ii<ntotscan;ii++) if (nsamples[ii] > ns) ns = nsamples[ii];



  indpix = new long[factdupl*nn*nn+2 + addnpix];
  init1D_long(indpix,0,factdupl*nn*nn+2 + addnpix,-1);



  ll=0;
  for (ii=0;ii<factdupl*nn*nn+2 + addnpix;ii++){
      if (pixon[ii] != 0){
	indpix[ii] = ll;
	ll++;
      }
  }
  npix = ll;


  delete [] pixon;


  printf("indpix[nn*nn] = %ld\n",indpix[nn*nn]);



  PNd = new double[npix];
  init1D_double(PNd,0,npix,0.0);


  if (pixout)
    printf("THERE ARE SAMPLES OUTSIDE OF MAP LIMITS: ASSUMING CONSTANT SKY EMISSION FOR THOSE SAMPLES, THEY ARE PUT IN A SINGLE PIXEL\n");
  printf("TOTAL NUMBER OF DETECTORS: %d\n",(int)ndet);
  printf("TOTAL NUMBER OF SCANS: %d\n", (int)ntotscan);
  printf("SIZE OF THE MAP: %d * %d\n",nn,nn);
  printf("TOTAL NUMBER OF FILLED PIXELS: %d\n",npix);

















  double *S0;
  S0 = new double[npix];
  init1D_double(S0,0,npix,0.0);


  if (MixMatfile != "NOFILE"){
    for (iframe=0;iframe<ntotscan;iframe++){
      ns = nsamples[iframe];
      ff = fframes[iframe];
      extentnoiseSp = extentnoiseSp_all[iframe];

      EstimPowerSpectra(fsamp,ns,ff,ndet,nn,npix,napod,marge,iframe,flgdupl,factdupl,indpix,S0,
			MixMatfile,bolonames,dirfile,bextension,fextension,cextension,shift_data_to_point,
			poutdir,termin,NORMLIN,NOFILLGAP,noiseSppreffile,extentnoiseSp,outdir);

    }
  }
  delete [] S0;











  //************************************************************************//
  //************************************************************************//
  //Pre-processing of the data
  //************************************************************************//
  //************************************************************************//



  for (iframe=0;iframe<ntotscan;iframe++){

      ns = nsamples[iframe];
      ff = fframes[iframe];
      f_lppix = f_lp*double(ns+2*marge)/fsamp;
      f_lppix_Nk = f_lp_Nk*double(ns+2*marge)/fsamp;
      prefixe = "fdata";


      if (CORRon){
	write_ftrProcesdata(NULL,indpix,indpsrc,nn,npix,npixsrc,ntotscan,addnpix,flgdupl,
			    factdupl,2,poutdir,termin,errarcsec,dirfile,scerr_field,
			    flpoint_field,bolonames,bextension,fextension,cextension,
			    shift_data_to_point,f_lppix,ff,ns,marge,napod,ndet,NORMLIN,
			    NOFILLGAP,iframe);

	do_PtNd(PNd,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,
		f_lppix_Nk,fsamp,ff,ns,marge,ndet,1,0,indpix,nn,npix,iframe,NULL,NULL);

      } else {

	do_PtNd_nocorr(PNd,extentnoiseSp_all,noiseSppreffile,poutdir,termin,errarcsec,dirfile,
		       scerr_field,flpoint_field,bolonames,bextension,fextension,
		       cextension,shift_data_to_point,f_lppix,f_lppix_Nk,fsamp,ntotscan,addnpix,
		       flgdupl,factdupl,2,ff,ns,marge,napod,ndet,1,0,indpix,indpsrc,
		       nn,npix,npixsrc,NORMLIN,NOFILLGAP,iframe,NULL);

      }

  } // end of iframe loop






  //******************************************************************//
  //******************************************************************//
  //**********************  End of init loop *************************//
  //******************************************************************//
  //******************************************************************//





  //************************************************************************//
  //************************************************************************//
  //main loop here over the fields and detectors
  //************************************************************************//
  //************************************************************************//


  int npixeff, iter, idupl;
  double var0, var_n, delta0, delta_n, delta_o, rtq, alpha, beta;

  double *S, *PtNPmatS, *r, *q, *d, *Mp, *s;
  long *hits;


  S = new double[npix];
  r = new double[npix];
  q = new double[npix];
  d = new double[npix];
  Mp = new double[npix];
  s = new double[npix];
  PtNPmatS = new double[npix];
  hits = new long[npix];


  long mi;
  double *map1d;
  map1d = new double[nn*nn];


  string fname;
  char iterchar[30];
  char iframechar[30];
  string iterstr, iframestr;




  for (idupl = 0;idupl<=flgdupl;idupl++){


  for (ii=0;ii<npix;ii++) S[ii] = 0.0;//PNd[ii];




  //Conjugate gradien Inversion
  if (projgaps || !flagon){
    npixeff = npix;
  } else {
    npixeff = npix-1;
  }




  printf("npix = %d, npixeff = %d\n", npix, npixeff);


  t1 = time(0);


  init1D_double(PtNPmatS,0,npix,0.0);
  init1D_double(Mp,0,npix,0.0);
  init1D_long(hits,0,npix,0);



  for (iframe=0;iframe<ntotscan;iframe++){
    ns = nsamples[iframe];
    ff = fframes[iframe];
    f_lppix_Nk = f_lp_Nk*double(ns+2*marge)/fsamp;
    prefixe = "fPs";

    if (CORRon){
      write_tfAS(S,indpix,nn,npix,flgdupl,factdupl, poutdir,termin,ff,ns,marge,ndet,iframe);


      do_PtNd(PtNPmatS,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,
	      f_lppix_Nk,fsamp,ff,ns,marge,ndet,1,0,indpix,nn,npix,iframe,Mp,hits);
    } else {
      do_PtNPS_nocorr(S,extentnoiseSp_all,noiseSppreffile,poutdir,termin,dirfile,bolonames,
		      f_lppix_Nk,fsamp,flgdupl,factdupl,ff,ns,marge,ndet,1,0,indpix,nn,npix,
		      iframe,PtNPmatS,Mp,hits);
    }
  } // end of iframe loop




  t2 = time(0);
  printf("temps de calcul: %ld\n",t2-t1);

  for (ii=0;ii<npixeff;ii++)
    if (Mp[ii] == 0)
      printf("ERROR: Mp[%ld] has elements = 0\n",ii);


  for (ii=0;ii<npixeff;ii++)
    Mp[ii] = 1.0/Mp[ii];


  for (ii=0;ii<npixeff;ii++)
    r[ii] = PNd[ii] - PtNPmatS[ii];

  for (ii=0;ii<npixeff;ii++)
    d[ii] =  Mp[ii] * r[ii];


  delta_n = 0.0;
  for (ii=0;ii<npixeff;ii++)
    delta_n += r[ii]*d[ii];

  var_n = 0.0;
  for (ii=0;ii<npixeff;ii++)
    var_n += r[ii]*r[ii];


  delta0 = delta_n;
  var0 = var_n;
  printf("var0 = %lf\n",var0);





  //start loop
  iter = 0;
  while(iter < 2000 && var_n/var0 > 1e-5 && (idupl || !flgdupl) || !idupl && var_n/var0 > 1e-4){

    init1D_double(q,0,npixeff,0.0);



    for (iframe=0;iframe<ntotscan;iframe++){
      ns = nsamples[iframe];
      ff = fframes[iframe];
      f_lppix_Nk = f_lp_Nk*double(ns+2*marge)/fsamp;
      prefixe = "fPs";

      if (CORRon){
	write_tfAS(d,indpix,nn,npix,flgdupl,factdupl, poutdir,termin,ff,ns,marge,ndet,iframe);

	do_PtNd(q,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,f_lppix_Nk,
		fsamp,ff,ns,marge,ndet,1,0,indpix,nn,npix,iframe,NULL,NULL);
      } else {

	do_PtNPS_nocorr(d,extentnoiseSp_all,noiseSppreffile,poutdir,termin,dirfile,bolonames,
			f_lppix_Nk,fsamp,flgdupl,factdupl,ff,ns,marge,ndet,1,0,indpix,nn,npix,
			iframe,q,NULL,NULL);
      }
    } // end of iframe loop



    rtq= 0.0;
    for (ii=0;ii<npixeff;ii++)
      rtq += q[ii] * d[ii];

    alpha = delta_n/rtq;


    for (ii=0;ii<npixeff;ii++)
      S[ii] += alpha*d[ii];




    if ((iter % 10) == 0){
      init1D_double(PtNPmatS,0,npixeff,0.0);


      for (iframe=0;iframe<ntotscan;iframe++){
	ns = nsamples[iframe];
	ff = fframes[iframe];
	f_lppix_Nk = f_lp_Nk*double(ns+2*marge)/fsamp;
	prefixe = "fPs";

	if (CORRon){
	  write_tfAS(S,indpix,nn,npix,flgdupl,factdupl, poutdir,termin,ff,ns,marge,ndet,iframe);


	  do_PtNd(PtNPmatS,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,
		  f_lppix_Nk,fsamp,ff,ns,marge,ndet,1,0,indpix,nn,npix,iframe,NULL,NULL);
	} else {
	  do_PtNPS_nocorr(S,extentnoiseSp_all,noiseSppreffile,poutdir,termin,dirfile,bolonames,
			  f_lppix_Nk,fsamp,flgdupl,factdupl,ff,ns,marge,ndet,1,0,indpix,nn,npix,
			  iframe,PtNPmatS,NULL,NULL);
	}
      } // end of iframe loop



      for (ii=0;ii<npixeff;ii++)
	r[ii] = PNd[ii] - PtNPmatS[ii];



    } else {


	for (ii=0;ii<npixeff;ii++)
	  r[ii] -= alpha*q[ii];

    }




    for (ii=0;ii<npixeff;ii++)
      s[ii] = Mp[ii]*r[ii];


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


    printf("iter = %d, crit = %10.15g, crit2 = %10.15g     \n",iter,var_n/var0,delta_n/delta0);




    if (iter == 0){
      for (ii=0; ii<nn; ii++) {
	for (jj=0; jj<nn; jj++) {
	  mi = jj*nn + ii;
	  if (indpix[mi] >= 0){
	    map1d[mi] = Mp[indpix[mi]];
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
	    map1d[mi] = Mp[indpix[mi]] * PNd[indpix[mi]];
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
	    map1d[mi] = hits[indpix[mi]];
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
		map1d[mi] += hits[indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc]];
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
		map1d[mi] += 1.0/Mp[indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc]];
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
		map1d[mi] += -S[indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc]]/Mp[indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc]];
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


    iter++;

  }
  printf("\n");







  if  ((projgaps || (flgdupl)) && !idupl){

    init1D_double(PNd,0,npix,0.0);

    for (iframe=0;iframe<ntotscan;iframe++){

      ns = nsamples[iframe];
      ff = fframes[iframe];
      f_lppix = f_lp*double(ns+2*marge)/fsamp;
      f_lppix_Nk = f_lp_Nk*double(ns+2*marge)/fsamp;
      prefixe = "fdata";

      if (CORRon){

	write_ftrProcesdata(S,indpix,indpsrc,nn,npix,npixsrc,ntotscan,addnpix,flgdupl,factdupl,2,
			    poutdir,termin,errarcsec,dirfile,scerr_field,flpoint_field,bolonames,
			    bextension,fextension,cextension,shift_data_to_point,f_lppix,ff,ns,
			    marge,napod,ndet,NORMLIN,NOFILLGAP,iframe);

	do_PtNd(PNd,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,f_lppix_Nk,
		fsamp,ff,ns,marge,ndet,1,0,indpix,nn,npix,iframe,NULL,NULL);
      } else {

	do_PtNd_nocorr(PNd,extentnoiseSp_all,noiseSppreffile,poutdir,termin,errarcsec,dirfile,
		       scerr_field,flpoint_field,bolonames,bextension,fextension,
		       cextension,shift_data_to_point,f_lppix,f_lppix_Nk,fsamp,ntotscan,addnpix,
		       flgdupl,factdupl,2,ff,ns,marge,napod,ndet,1,0,indpix,indpsrc,
		       nn,npix,npixsrc,NORMLIN,NOFILLGAP,iframe,S);
      }
    } // end of iframe loop



  }

  }









  //******************************  write final map in file ********************************


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
	  map1d[mi] = Mp[indpix[mi]];
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
		map1d[mi] = Mp[indpix[indpsrc[mi]+factdupl*nn*nn+npixsrc*iframe]];
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




  //*******************************************************************//
  //******************  Update noise power spectra  *******************//


  if (MixMatfile != "NOFILE"){
    for (iframe=0;iframe<ntotscan;iframe++){
      ns = nsamples[iframe];
      ff = fframes[iframe];
      extentnoiseSp = extentnoiseSp_all[iframe];

      EstimPowerSpectra(fsamp,ns,ff,ndet,nn,npix,napod,marge,iframe,flgdupl,factdupl,indpix,S,
			MixMatfile,bolonames,dirfile,bextension,fextension,cextension,shift_data_to_point,
			poutdir,termin,NORMLIN,NOFILLGAP,noiseSppreffile,extentnoiseSp,outdir);

    }
  }





  //******************************************************************//
  //******************************************************************//
  //*********************** End of program *************************//
  //******************************************************************//
  //******************************************************************//







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
    fprintf(fp,"%s\t",argv[ii]);
  fprintf(fp,"\n");
  fclose(fp);




  // clean up
  delete [] ra;
  delete [] dec;
  delete [] phi;
  delete [] scerr;
  delete [] xx;
  delete [] yy;
  delete [] flag;
  delete [] rejectsamp;
  delete [] samptopix;
  delete [] flpoint;
  delete [] mask;



  return 0;
}





void write_tfAS(double *S, long *indpix, int nn, long npix, bool flgdupl, int factdupl, string dir, string termin, long ff, long ns, long marge, long ndet, long iframe){


  long idet1;
  long ndata = ns+2*marge;

  FILE *fp;
  char testfile[100];

  double *Ps;
  long *samptopix;

  fftw_plan fftplan;
  fftw_complex *fdata;

  samptopix = new long[ns];
  Ps = new double[ndata];
  fdata = new fftw_complex[ndata/2+1];


  //for (idet1=rank*ndet/size;idet1<(rank+1)*ndet/size;idet1++){
  for (idet1=0;idet1<ndet;idet1++){

    //Read pointing data
    sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet1,"_",termin.c_str(),".bi");
    fp = fopen(testfile,"r");
    fread(samptopix,sizeof(long),ns,fp);
    fclose(fp);

    deproject(S,indpix,samptopix,ndata,marge,nn,npix,Ps,flgdupl,factdupl);

    //Fourier transform of the data
    fftplan = fftw_plan_dft_r2c_1d(ndata, Ps, fdata, FFTW_ESTIMATE);
    fftw_execute(fftplan);
    fftw_destroy_plan(fftplan);

    sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"fPs_",iframe,"_",idet1,"_",termin.c_str(),".bi");
    fp = fopen(testfile,"w");
    fwrite(fdata,sizeof(double), (ndata/2+1)*2, fp);
    fclose(fp);

  }

  delete[] samptopix;
  delete[] Ps;
  delete[] fdata;

}







void write_ftrProcesdata(double *S, long *indpix, long *indpsrc, int nn, long npix,
			 long npixsrc, long ntotscan, long addnpix, bool flgdupl, int factdupl,
			 int fillg, string dir, string termin, double errarcsec, string dirfile,
			 string scerr_field, string flpoint_field, string *bolonames,
			 string bextension, string fextension, string cextension,
			 int shift_data_to_point, double f_lppix, long ff, long ns,
			 long marge, long napod, long ndet, bool NORMLIN, bool NOFILLGAP,
			 long iframe){


  long ii, idet1;
  long ndata = ns+2*marge;

  double *scerr, *data, *calp, *bfilter, *data_lp, *Ps;
  unsigned char *flpoint, *flag, *rejectsamp;
  long *samptopix;

  fftw_plan fftplan;
  fftw_complex *fdata;

  string field1;

  char testfile[100];

  FILE *fp;

  scerr = new double[ndata];
  data =  new double[ndata];
  data_lp = new double[ndata];
  calp =  new double[ndata];
  flag =  new unsigned char[ndata];
  flpoint = new unsigned char[ndata];
  rejectsamp = new unsigned char[ndata];

  samptopix = new long[ns];
  Ps = new double[ndata];
  bfilter = new double[ndata/2+1];
  fdata = new fftw_complex[ndata/2+1];



  for (idet1=0;idet1<ndet;idet1++){

    field1 = bolonames[idet1];

    if (S != NULL){
      read_data(dirfile, ff, 0, ns, scerr, scerr_field, 'd');
      read_data(dirfile, ff, 0, ns, flpoint, flpoint_field, 'c');
    }

    read_data(dirfile, ff, shift_data_to_point, ns, data, field1+bextension, 'd');

    if (fextension != "NOFLAG"){
      read_data(dirfile, ff, shift_data_to_point, ns, flag, field1+fextension,  'c');
    } else {
      printf("NOFLAG\n");
      for (ii=0;ii<ns;ii++)
	flag[ii] = 0;
    }

    if (cextension != "NOCALP"){
      read_data(dirfile, ff, 0, ns/20, calp, field1+cextension, 'd');
    } else {
      printf("NOCALP\n");
      for (ii=0;ii<ns/20;ii++)
	calp[ii] = 1.0;
    }



    if (S != NULL){
      //// Read pointing
      sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet1,"_",termin.c_str(),".bi");
      fp = fopen(testfile,"r");
      fread(samptopix,sizeof(long),ns,fp);
      fclose(fp);

      if (addnpix){
	deproject(S,indpix,samptopix,ns+2*marge,marge,nn,npix,Ps,fillg,factdupl,ntotscan,indpsrc,npixsrc);
      } else {
	deproject(S,indpix,samptopix,ns+2*marge,marge,nn,npix,Ps,fillg,factdupl);
      }

      for (ii=0;ii<ns;ii++) rejectsamp[ii] = 0;
      for (ii=0;ii<ns;ii++)
	if ((flag[ii] & 1) != 0 || (scerr[ii] > errarcsec) || (flpoint[ii] & 1) != 0)
	  rejectsamp[ii] = 1;
    }


    if (S != NULL){
      //********************  pre-processing of data ********************//
      MapMakPreProcessData(data,flag,calp,ns,marge,napod,4,f_lppix,data_lp,bfilter,
			   NORMLIN,NOFILLGAP,Ps);
    }
    else {
      MapMakPreProcessData(data,flag,calp,ns,marge,napod,4,f_lppix,data_lp,bfilter,
			   NORMLIN,NOFILLGAP);
    }


    //Fourier transform of the data
    fftplan = fftw_plan_dft_r2c_1d(ndata, data_lp, fdata, FFTW_ESTIMATE);
    fftw_execute(fftplan);
    fftw_destroy_plan(fftplan);


    sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"fdata_",iframe,"_",idet1,"_",termin.c_str(),".bi");
    fp = fopen(testfile,"w");
    fwrite(fdata,sizeof(double), (ndata/2+1)*2, fp);
    fclose(fp);

  }


  delete[] scerr;
  delete[] data;
  delete[] data_lp;
  delete[] calp;
  delete[] flag;
  delete[] flpoint;
  delete[] samptopix;
  delete[] Ps;
  delete[] bfilter;
  delete[] fdata;
  delete[] rejectsamp;


}





void do_PtNd(double *PNd, string *extentnoiseSp_all, string noiseSppreffile,
	     string dir, string prefixe, string termin, string *bolonames,
	     double f_lppix, double fsamp, long ff, long ns, long marge, long ndet, int size,
	     int rank, long *indpix, long nn, long npix, long iframe, double *Mp, long *hits){


  long ii, jj, idet1, idet2, nbins;
  double dnbins;
  long ndata = ns+2*marge;
  string field1, field2;
  string extentnoiseSp;

  char nameSpfile[100];
  char testfile[100];

  long *samptopix;
  double *ell, *SpN, *bfilter, *bfilter_, *Nk, *Nd;

  fftw_plan fftplan;
  fftw_complex *fdata, *Ndf;


  samptopix = new long[ns];
  Nd = new double[ndata];
  bfilter = new double[ndata/2+1];
  bfilter_ = new double[ndata/2+1];
  Nk = new double[ndata/2+1];
  fdata = new fftw_complex[ndata/2+1];
  Ndf = new fftw_complex[ndata/2+1];

  double **SpN_all;

  FILE *fp;




  for (idet1=rank*ndet/size;idet1<(rank+1)*ndet/size;idet1++){
    field1 = bolonames[idet1];

    //Read pointing data
    sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet1,"_",termin.c_str(),".bi");
    fp = fopen(testfile,"r");
    fread(samptopix,sizeof(long),ns,fp);
    fclose(fp);


    //**************************************** Noise power spectrum
    extentnoiseSp = extentnoiseSp_all[iframe];
    sprintf(nameSpfile,"%s%s%s%s",noiseSppreffile.c_str(),field1.c_str(),"-all",extentnoiseSp.c_str());
    if ((fp = fopen(nameSpfile,"r")) == NULL){
      cerr << "ERROR: Can't find noise power spectra file" << nameSpfile << " , check -k or -K in command line. Exiting. \n";
      exit(1);
    }
    fread(&dnbins,sizeof(double), 1, fp);
    nbins = (long)dnbins;
    SpN_all = dmatrix(0,ndet-1,0,nbins-1);
    ell = new double[nbins+1];
    SpN = new double[nbins];
    fread(ell,sizeof(double), nbins+1, fp);
    fread(*SpN_all,sizeof(double), nbins*ndet, fp);
    fclose(fp);
    //*****************************************

    for (ii=0;ii<ndata/2+1;ii++)
      bfilter[ii] = pow(double(ii)/f_lppix, 16) /(1.0+pow(double(ii)/f_lppix, 16));
    for (ii=0;ii<ndata/2+1;ii++)
      bfilter_[ii] = 1.0/(bfilter[ii]+0.000001);


    //Init N-1d
    for (ii=0;ii<ndata/2+1;ii++){
      Ndf[ii][0] = 0;
      Ndf[ii][1] = 0;
    }



    for (idet2=0;idet2<ndet;idet2++){
      field2 = bolonames[idet2];

      //read Fourier transform of the data
      sprintf(testfile,"%s%s%s%ld%s%ld%s%s%s",dir.c_str(),prefixe.c_str(),"_",iframe,"_",idet2,"_",termin.c_str(),".bi");
      fp = fopen(testfile,"r");
      fread(fdata,sizeof(double), (ndata/2+1)*2, fp);
      fclose(fp);


      //****************** Cross power spectrum of the noise  ***************//
      for (ii=0;ii<nbins;ii++)
	SpN[ii] = SpN_all[idet2][ii];


      // interpolate logarithmically the noise power spectrum
      InvbinnedSpectrum2log_interpol(ell,SpN,bfilter_,nbins,ndata,fsamp,Nk);


      for (jj=0;jj<ndata/2+1;jj++)
	if (isnan(Nk[jj]))
	  printf("Ca ne va pas fr %ld, det1 %ld, det2 %ld\n",iframe, idet1, idet2);


      //********************************* compute N^-1 d  ***********************//
      for (ii=0;ii<ndata/2+1;ii++){
	Ndf[ii][0] += fdata[ii][0]*Nk[ii];
	Ndf[ii][1] += fdata[ii][1]*Nk[ii];
      }



      //Compute weight map for preconditioner
      if ((Mp != NULL) && (idet2 == idet1))
	compute_diagPtNPCorr(Nk,samptopix,ndata,marge,nn,indpix,npix,f_lppix,Mp);


    }// end of idet2 loop


    fftplan = fftw_plan_dft_c2r_1d(ndata, Ndf, Nd, FFTW_ESTIMATE);
    fftw_execute(fftplan);
    fftw_destroy_plan(fftplan);



    for (ii=-marge;ii<ndata-marge;ii++){
      if ((ii < 0) || (ii >= ndata-2*marge)){
	PNd[npix-2] += Nd[ii+marge];
      } else {
	PNd[indpix[samptopix[ii]]] += Nd[ii+marge];
      }
    }

    //compute hit counts
    if (hits != NULL){
      for (ii=0;ii<ndata-2*marge;ii++){
	hits[indpix[samptopix[ii]]] += 1;
      }
    }


    delete[] ell;
    delete[] SpN;
    free_dmatrix(SpN_all,0,ndet-1,0,nbins-1);

  }// end of idet1 loop



  delete[] samptopix;
  delete[] Nd;
  delete[] bfilter;
  delete[] bfilter_;
  delete[] Nk;
  delete[] fdata;
  delete[] Ndf;


}



void do_PtNd_nocorr(double *PNd, string *extentnoiseSp_all, string noiseSppreffile,
		    string dir, string termin, double errarcsec, string dirfile,
		    string scerr_field, string flpoint_field, string *bolonames,
		    string bextension, string fextension, string cextension,
		    int shift_data_to_point, double f_lppix, double f_lppix_Nk,
		    double fsamp, long ntotscan, long addnpix, bool flgdupl, int factdupl,
		    int fillg, long ff, long ns, long marge, long napod, long ndet,
		    int size, int rank, long *indpix, long *indpsrc, long nn, long npix,
		    long npixsrc, bool NORMLIN, bool NOFILLGAP, long iframe, double *S){



  long ii, idet;
  long ndata = ns+2*marge;
  string field;
  string extentnoiseSp;

  char nameSpfile[100];
  char testfile[100];

  long *samptopix;
  double *bfilter, *Nk, *data, *data_lp, *scerr, *calp, *Ps;
  unsigned char *flag, *flpoint, *rejectsamp;


  samptopix = new long[ns];
  bfilter = new double[ndata/2+1];
  Nk = new double[ndata/2+1];

  scerr = new double[ndata];
  data =  new double[ndata];
  data_lp = new double[ndata];
  calp =  new double[ndata];
  flag =  new unsigned char[ndata];
  flpoint = new unsigned char[ndata];
  rejectsamp = new unsigned char[ndata];
  Ps = new double[ndata];


  FILE *fp;



  for (idet=rank*ndet/size;idet<(rank+1)*ndet/size;idet++){

    field = bolonames[idet];



    if (S != NULL){
      read_data(dirfile, ff, 0, ns, scerr, scerr_field, 'd');
      read_data(dirfile, ff, 0, ns, flpoint, flpoint_field, 'c');
    }

    read_data(dirfile, ff, shift_data_to_point, ns, data, field+bextension, 'd');

    if (fextension != "NOFLAG"){
      read_data(dirfile, ff, shift_data_to_point, ns, flag, field+fextension,  'c');
    } else {
      printf("NOFLAG\n");
      for (ii=0;ii<ns;ii++)
	flag[ii] = 0;
    }

    if (cextension != "NOCALP"){
      read_data(dirfile, ff, 0, ns/20, calp, field+cextension, 'd');
    } else {
      printf("NOCALP\n");
      for (ii=0;ii<ns/20;ii++)
	calp[ii] = 1.0;
    }


    //// Read pointing
    sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet,"_",termin.c_str(),".bi");
    fp = fopen(testfile,"r");
    fread(samptopix,sizeof(long),ns,fp);
    fclose(fp);


    if (S != NULL){

      if (addnpix){
	deproject(S,indpix,samptopix,ns+2*marge,marge,nn,npix,Ps,fillg,factdupl,ntotscan,indpsrc,npixsrc);
      } else {
	deproject(S,indpix,samptopix,ns+2*marge,marge,nn,npix,Ps,fillg,factdupl);
      }

      for (ii=0;ii<ns;ii++) rejectsamp[ii] = 0;
      for (ii=0;ii<ns;ii++)
	if ((flag[ii] & 1) != 0 || (scerr[ii] > errarcsec) || (flpoint[ii] & 1) != 0)
	  rejectsamp[ii] = 1;
    }


    if (S != NULL){
      //********************  pre-processing of data ********************//
      MapMakPreProcessData(data,flag,calp,ns,marge,napod,4,f_lppix,data_lp,bfilter,
			   NORMLIN,NOFILLGAP,Ps);
    }
    else {
      MapMakPreProcessData(data,flag,calp,ns,marge,napod,4,f_lppix,data_lp,bfilter,
			   NORMLIN,NOFILLGAP);
    }


    for (ii=0;ii<ndata/2+1;ii++)
      bfilter[ii] = pow(double(ii)/f_lppix_Nk, 16) /(1.0+pow(double(ii)/f_lppix_Nk, 16));



    //****************** Compute (or read) input power spectrum of the NOISE  ***************//
    extentnoiseSp = extentnoiseSp_all[iframe];
    sprintf(nameSpfile,"%s%s%s",noiseSppreffile.c_str(),field.c_str(),extentnoiseSp.c_str());
    readNSpectrum(nameSpfile,bfilter,ns,marge,fsamp,Nk);


    //********************** compute P^t N-1 d ************************//
    compute_PtNmd(data_lp,Nk,ns+2*marge,marge,nn,indpix,samptopix,npix,PNd);


  }// end of idet loop


  delete[] samptopix;
  delete[] bfilter;
  delete[] Nk;
  delete[] scerr;
  delete[] data;
  delete[] data_lp;
  delete[] calp;
  delete[] flag;
  delete[] flpoint;
  delete[] rejectsamp;
  delete[] Ps;


}





void do_PtNPS_nocorr(double *S, string *extentnoiseSp_all, string noiseSppreffile, string dir,
		string termin, string dirfile, string *bolonames, double f_lppix,
		double fsamp, bool flgdupl, int factdupl, long ff, long ns, long marge,
		long ndet, int size, int rank, long *indpix, long nn, long npix,
		long iframe, double *PtNPmatS, double *Mp, long *hits){



  long ii, idet;
  long ndata = ns+2*marge;
  string field;
  string extentnoiseSp;

  char nameSpfile[100];
  char testfile[100];

  long *samptopix;
  double *bfilter, *Nk, *Ps;


  samptopix = new long[ns];
  bfilter = new double[ndata/2+1];
  Nk = new double[ndata/2+1];
  Ps = new double[ndata];

  FILE *fp;



  for (idet=rank*ndet/size;idet<(rank+1)*ndet/size;idet++){

    field = bolonames[idet];


    sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet,"_",termin.c_str(),".bi");
    fp = fopen(testfile,"r");
    fread(samptopix,sizeof(long),ns,fp);
    fclose(fp);


    // AS
    deproject(S,indpix,samptopix,ndata,marge,nn,npix,Ps,flgdupl,factdupl);


    extentnoiseSp = extentnoiseSp_all[iframe];
    sprintf(nameSpfile,"%s%s%s",noiseSppreffile.c_str(),field.c_str(),extentnoiseSp.c_str());
    for (ii=0;ii<(ns+2*marge)/2+1;ii++)
      bfilter[ii] = pow(double(ii)/f_lppix, 16) /(1.0+pow(double(ii)/f_lppix, 16));
    readNSpectrum(nameSpfile,bfilter,ns,marge,fsamp,Nk);





    //AtN-1A AS (espensive part)
    compute_PtNmd(Ps,Nk,ns+2*marge,marge,nn,indpix,samptopix,npix,PtNPmatS);


    //Compute weight map for preconditioner
    if ((Mp != NULL))
      compute_diagPtNP(Nk,samptopix,ndata,marge,nn,indpix,npix,f_lppix,Mp);


    //compute hit counts
    if (hits != NULL){
      for (ii=0;ii<ndata-2*marge;ii++){
	hits[indpix[samptopix[ii]]] += 1;
      }
    }


  } // end of idet loop


  delete[] samptopix;
  delete[] bfilter;
  delete[] Nk;
  delete[] Ps;


}







void EstimPowerSpectra(double fsamp, long ns, long ff, long ndet, int nn, long npix, long napod,
		       long marge, long iframe, bool flgdupl, int factdupl, long *indpix,
		       double *S, string MixMatfile, string *bolonames, string dirfile, string bextension,
		       string fextension, string cextension, int shift_data_to_point, string dir,
		       string termin, bool NORMLIN, bool NOFILLGAP, string noiseSppreffile,
		       string extentnoiseSp, string outdirSpN)
{



  double fcut = 3.0;
  long ncomp = 3;

  long ii, jj, kk, ll, idet, idet1, idet2, ib;
  long ncomp2;
  int nbins = 500;
  int nbins2;
  double tmpsign, mm, sign0, factapod;

  double dnbins, dummy1;

  double *data, *data_lp, *Ps, *calp, *apodwind, *commontmp, *commonm_f, *bfilter, *SPref;
  long *samptopix;
  unsigned char *flag;
  double **commonm, **commonm2, **common_f, **vect;
  double **P, **N, **Rellth, **Rellexp;
  double **aa, **Cov, **iCov, **iCov2, **SpN_all;
  double *Nk;
  double *ell, *SpN, *Nell;
  double *sign;
  double *p, *uvec, *ivec;

  fftw_plan fftplan;
  fftw_complex *fdata1, *fdata2;

  FILE *fp;

  char testfile[100];
  char nameSpfile[200];

  string field;
  string tempstr1;
  string tempstr2;


  data = new double[ns];
  data_lp = new double[ns];
  Ps = new double[ns];
  commontmp = new double[ns];
  commonm_f = new double[ns];
  calp = new double[ns];
  flag = new unsigned char[ns];
  samptopix = new long[ns];
  Nk = new double[ns/2+1];
  bfilter = new double[ns/2+1];
  fdata1 = new fftw_complex[ns/2+1];
  fdata2 = new fftw_complex[ns/2+1];

  Nell = new double[nbins];
  SPref = new double[nbins];
  P = dmatrix(0,ncomp-1,0,nbins-1);
  N = dmatrix(0,ndet-1,0,nbins-1);
  Rellexp = dmatrix(0,ndet*ndet-1,0,nbins-1);
  Rellth = dmatrix(0,ndet*ndet-1,0,nbins-1);
  aa = dmatrix(0,ndet-1,0,20);

  sign = new double[ndet];

  Cov = dmatrix(0,ncomp-1,0,ncomp-1);
  iCov = dmatrix(0,ncomp-1,0,ncomp-1);
  iCov2 = dmatrix(0,ncomp-1,0,ncomp-1);
  p = new double[ncomp];
  uvec = new double[ncomp];
  ivec = new double[ncomp];

  vect = dmatrix(0,ncomp-1,0,ndet-1);

  commonm = dmatrix(0,ncomp,0,ns-1);
  commonm2 = dmatrix(0,ncomp,0,ns-1);
  common_f = dmatrix(0,ncomp,0,ns-1);

  init2D_double(commonm,0,0,ncomp,ns,0.0);
  init2D_double(commonm2,0,0,ncomp,ns,0.0);
  init2D_double(common_f,0,0,ncomp,ns,0.0);
  init1D_double(commontmp,0,ns,0.0);
  init1D_double(commonm_f,0,ns,0.0);
  init2D_double(Cov,0,0,ncomp,ncomp,0.0);
  init2D_double(iCov,0,0,ncomp,ncomp,0.0);
  init2D_double(iCov2,0,0,ncomp,ncomp,0.0);

  for (ii=0;ii<ns/2+1;ii++)
    bfilter[ii] = 1.0;




  //////////////////////////////////////
  // read mixing parameters
  if ((fp = fopen(MixMatfile.c_str(),"r")) == NULL){
    cerr << "ERROR: Can't find Mixing Matrix file. Exiting. \n";
    return;
  }
  fscanf(fp,"%d",&ncomp2);
  for (ii=0;ii<ndet;ii++){
    for (jj=0;jj<ncomp2;jj++){
      fscanf(fp,"%lf",&dummy1);
      aa[ii][jj] = dummy1;
    }
  }

  if (ncomp2 < ncomp) ncomp = ncomp2;




  //**************************************************************************************
  //**************************** Read data and compute components

  //apodization
  apodwind = apodwindow(ns,int(ns*0.04));


  // loop over detectors
  for (idet=0;idet<ndet;idet++){

    field = bolonames[idet];


    read_data(dirfile, ff, shift_data_to_point, ns, data, field+bextension, 'd');

    if (fextension != "NOFLAG"){
      read_data(dirfile, ff, shift_data_to_point, ns, flag, field+fextension,  'c');
    } else {
      printf("NOFLAG\n");
      for (ii=0;ii<ns;ii++)
	flag[ii] = 0;
    }

    if (cextension != "NOCALP"){
      read_data(dirfile, ff, 0, ns/20, calp, field+cextension, 'd');
    } else {
      printf("NOCALP\n");
      for (ii=0;ii<ns/20;ii++)
	calp[ii] = 1.0;
    }



    //******************************* subtract signal

    //Read pointing data
    sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet,"_",termin.c_str(),".bi");
    fp = fopen(testfile,"r");
    fread(samptopix,sizeof(long),ns,fp);
    fclose(fp);

    deproject(S,indpix,samptopix,ns+2*marge,marge,nn,npix,Ps,flgdupl,factdupl);


    data[ii] = data[ii] - Ps[marge+ii]/calp[ii/20];


    MapMakPreProcessData(data,flag,calp,ns,marge,napod,4,1.0,data_lp,bfilter,
			 NORMLIN,NOFILLGAP);


    for (ii=0;ii<ns;ii++)
      data[ii] = data_lp[ii]*apodwind[ii];



    //sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"data_fortest_",idet,"_",ff,"_",termin.c_str(),".bi");
    //fp = fopen(testfile,"w");
    //fwrite(data,sizeof(double), ns, fp);
    //fclose(fp);



    // compute fft and save data to disk for later
    fftplan = fftw_plan_dft_r2c_1d(ns, data, fdata1, FFTW_ESTIMATE);
    fftw_execute(fftplan);
    fftw_destroy_plan(fftplan);



    // save fft to disk for later
    sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"fdata_",idet,"_",ff,"_",termin.c_str(),".bi");
    fp = fopen(testfile,"w");
    fwrite(fdata1,sizeof(double), (ns/2+1)*2, fp);
    fclose(fp);



    /// compute sigma of the noise
    mm = 0.0;
    for (ii=ns/2;ii<ns/2+500;ii++) mm += data[ii];
    mm = mm/501.0;
    tmpsign = 0.0;
    for (ii=ns/2;ii<ns/2+500;ii++) tmpsign += (data[ii]-mm)*(data[ii]-mm);
    sign[idet] = sqrt(tmpsign/500.0);
    // normalize to the first detector
    if (idet == 0) sign0 = sign[0];
    sign[idet] = sign[idet]/sign0;



    for (jj=0;jj<ncomp;jj++)
      for (ii=0;ii<ns;ii++)
	commonm[jj][ii] += aa[idet][jj]/(sign[idet]*sign[idet])*data[ii];


  }

  //***************************************************************************


  /////////// AtN-1A
  for (jj=0;jj<ncomp;jj++)
    for (kk=0;kk<ncomp;kk++)
      for (ii=0;ii<ndet;ii++)
	Cov[jj][kk] += aa[ii][jj] * aa[ii][kk]/sign[ii]/sign[ii];


  // invert AtN-1A
  dcholdc(Cov,ncomp,p);
  for (ii=0;ii<ncomp;ii++){
    for (jj=0;jj<ncomp;jj++)
      uvec[jj] = 0.0;
    uvec[ii] = 1.0;
    dcholsl(Cov,ncomp,p,uvec,ivec);
    for (jj=0;jj<ncomp;jj++)
      iCov[ii][jj] = ivec[jj];
  }


  for (ii=0;ii<ns;ii++)
    for (jj=0;jj<ncomp;jj++)
      for (kk=0;kk<ncomp;kk++)
	commonm2[jj][ii] += iCov[jj][kk] * commonm[kk][ii];



  factapod = 0.0;
  for (ii=0;ii<ns;ii++)
    factapod += apodwind[ii]*apodwind[ii]/ns;



  /// filter common mode to fcut Hz
  for (ii=0;ii<ncomp;ii++){
    for (jj=0;jj<ns;jj++)
      commontmp[jj] = commonm2[ii][jj];
    //butterworth(commontmp,ns,fcut/fsamp*double(ns),8,commonm_f,bfilter,1,0,0);
    for (jj=0;jj<ns;jj++)
      common_f[ii][jj] = commontmp[jj];                         //// - commonm_f[jj];
  }



  //**************************************** Read pre-estimated power spectrum for reference
  sprintf(nameSpfile,"%s%s%s%s",noiseSppreffile.c_str(),field.c_str(),"-all",extentnoiseSp.c_str());
  if ((fp = fopen(nameSpfile,"r")) == NULL){
    cerr << "ERROR: Can't find noise power spectra file" << nameSpfile << " , check -k or -K in command line. Exiting. \n";
    exit(1);
  }
  fread(&dnbins,sizeof(double), 1, fp);
  nbins = (long)dnbins;
  nbins2 = nbins;
  SpN_all = dmatrix(0,ndet-1,0,nbins-1);
  ell = new double[nbins+1];
  SpN = new double[nbins];
  fread(ell,sizeof(double), nbins+1, fp);
  fread(*SpN_all,sizeof(double), nbins*ndet, fp);
  fclose(fp);
  delete [] SpN;
  //*****************************************




  //************************************************************************//
  // second part: -- data - common mode
  //              -- estimate noise power spectra


  /////////////////////////////////////// loop again over detectors
  for (idet=0;idet<ndet;idet++){

    field = bolonames[idet];


    read_data(dirfile, ff, shift_data_to_point, ns, data, field+bextension, 'd');

    if (fextension != "NOFLAG"){
      read_data(dirfile, ff, shift_data_to_point, ns, flag, field+fextension,  'c');
    } else {
      printf("NOFLAG\n");
      for (ii=0;ii<ns;ii++)
	flag[ii] = 0;
    }

    if (cextension != "NOCALP"){
      read_data(dirfile, ff, 0, ns/20, calp, field+cextension, 'd');
    } else {
      printf("NOCALP\n");
      for (ii=0;ii<ns/20;ii++)
	calp[ii] = 1.0;
    }



    //******************************* subtract signal

    //Read pointing data
    sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet,"_",termin.c_str(),".bi");
    fp = fopen(testfile,"r");
    fread(samptopix,sizeof(long),ns,fp);
    fclose(fp);

    deproject(S,indpix,samptopix,ns+2*marge,marge,nn,npix,Ps,flgdupl,factdupl);


    data[ii] = data[ii] - Ps[marge+ii]/calp[ii/20];


    MapMakPreProcessData(data,flag,calp,ns,marge,napod,4,1.0,data_lp,bfilter,
			 NORMLIN,NOFILLGAP);


    for (ii=0;ii<ns;ii++)
      data[ii] = data_lp[ii] * apodwind[ii];


    // Subtract components
    for (ii=0;ii<ns;ii++)
      for (jj=0;jj<ncomp;jj++)
	data[ii] -= aa[idet][jj]*common_f[jj][ii];


    //Noise power spectra
    ///// measure power spectrum of the uncorrelated part of the noise
    noisepectrum_estim(data,ns,ell,nbins,fsamp,NULL,Nell,Nk);

    for (ii=0;ii<nbins;ii++){
      Rellth[idet*ndet+idet][ii] += Nell[ii]/factapod;
      N[idet][ii] = Nell[ii]/factapod;
    }


  }




  ////*********************** Init component power spectra

  for (ii=0;ii<ncomp;ii++){
    for (jj=0;jj<ns;jj++)
      commontmp[jj]=common_f[ii][jj];
    noisepectrum_estim(commontmp,ns,ell,nbins,fsamp,NULL,Nell,Nk);
    for (jj=0;jj<nbins;jj++)
      P[ii][jj] = Nell[jj]/factapod;

    // subtract a factor to correct from noise
    //for (jj=0;jj<nbins;jj++)
    //  P[ii][jj] -= iCov[ii][ii]*sign0*sign0;
    //for (jj=0;jj<nbins;jj++)
    //  if (P[ii][jj] < 0)
    //	P[ii][jj] = 0.0;

  }
  printf("noise var det 0 =  %10.15g\n",sign0*sign0);




  for (ii=0;ii<ndet;ii++)
    for (kk=0;kk<ndet;kk++)
      for (ll=0;ll<ncomp;ll++)
	for (jj=0;jj<nbins;jj++)
	  Rellth[ii*ndet+kk][jj] += aa[ii][ll] * aa[kk][ll] * P[ll][jj];





  //************************************************************************//
  // third part: -- estimate the covariance matrix of the data R_exp
  //



  /////////////////////////////////////// loop again over detectors
  for (idet1=0;idet1<ndet;idet1++){

    // read data from disk
    sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"fdata_",idet1,"_",ff,"_",termin.c_str(),".bi");
    fp = fopen(testfile,"r");
    fread(fdata1,sizeof(double), (ns/2+1)*2, fp);
    fclose(fp);

    for (idet2=0;idet2<ndet;idet2++) {

      // read data from disk
      sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"fdata_",idet2,"_",ff,"_",termin.c_str(),".bi");
      fp = fopen(testfile,"r");
      fread(fdata2,sizeof(double), (ns/2+1)*2, fp);
      fclose(fp);

      noisecrosspectrum_estim(fdata1,fdata2,ns,ell,nbins,fsamp,NULL,Nell,Nk);

      for (ii=0;ii<nbins;ii++)
	Rellexp[idet1*ndet+idet2][ii] += Nell[ii]/factapod;

    }
    printf("Computing Rellexp, idet = %d\n",idet1);
  }



  ////// normalize to the first detector power spectrum in order to avoid numerical problems
  for (ii=0;ii<nbins;ii++)
    SPref[ii] = Rellexp[0][ii];
  for (idet1=0;idet1<ndet;idet1++)
    for (idet2=0;idet2<ndet;idet2++)
      for (ii=0;ii<nbins;ii++)
	Rellexp[idet1*ndet+idet2][ii] = Rellexp[idet1*ndet+idet2][ii]/SPref[ii];
  for (jj=0;jj<ncomp;jj++)
    for (ii=0;ii<nbins;ii++)
      P[jj][ii] = P[jj][ii]/SPref[ii];
  for (jj=0;jj<ndet;jj++)
    for (ii=0;ii<nbins;ii++)
      N[jj][ii] = N[jj][ii]/SPref[ii];







  /*
  //// write Rellexp to disk and also first guess of parameters
  sprintf(testfile,"%s%s%d%s%s",outdirSpN.c_str(),"Rellexp_",ff,termin.c_str(),".txt");
  fp = fopen(testfile,"w");
  for (jj=0;jj<nbins;jj++)
    for (ii=0;ii<ndet;ii++)
      for (kk=0;kk<ndet;kk++)
	fprintf(fp,"%10.15g \t",Rellexp[ii*ndet+kk][jj]);
  fprintf(fp,"\n");
  fclose(fp);

  sprintf(testfile,"%s%s%d%s%s",outdirSpN.c_str(),"Ninit_",ff,termin.c_str(),".txt");
  fp = fopen(testfile,"w");
  for (ii=0;ii<ndet;ii++)
    for (jj=0;jj<nbins;jj++)
      fprintf(fp,"%10.15g \t",N[ii][jj]);
  fprintf(fp,"\n");
   fclose(fp);



  sprintf(testfile,"%s%s%d%s%s",outdirSpN.c_str(),"Pinit_",ff,termin.c_str(),".txt");
  fp = fopen(testfile,"w");
  for (ii=0;ii<ncomp;ii++)
    for (jj=0;jj<nbins;jj++)
      fprintf(fp,"%10.15g \t",P[ii][jj]);
  fprintf(fp,"\n");
  fclose(fp);


  sprintf(testfile,"%s%s%d%s%s",outdirSpN.c_str(),"Ainit_",ff,termin.c_str(),".txt");
  fp = fopen(testfile,"w");
  for (ii=0;ii<ndet;ii++)
    for (jj=0;jj<ncomp;jj++)
      fprintf(fp,"%10.15g \t",aa[ii][jj]);
  fprintf(fp,"\n");
  fclose(fp);
  */








  //***** Fourth part
  //*********************** fit component and noise power spectra, and mixing matrix *************//
  //********* Using Expectation/Maximization algorithm
  long nbiter = 500;
  long iter;

  double tottest=0.0;

  double f;

  double *iN, *Pr, *w;
  double **Rxs, **Rxsq, **RnRxsb, **Rxx, **Rxxq, **Rss, **Rssq, **RnRssb;
  double **Pr2, **AiNA, **Mattmp, **ImDR, **ACq, **Cq, **Wq;


  ib = 0;
  while ((ell[ib] < fcut) && (ib < nbins)){
    nbins2 = ib+1;
    ib++;
  }

  /////// to be removed:
  //nbins2 = nbins;

  printf("nbins2 = %d\n",nbins2);



  iN = new double[ndet];
  Pr = new double[ncomp];
  w = new double[nbins2];

  Rxs =    dmatrix(0,ndet-1 ,0,ncomp-1);
  Rxsq =   dmatrix(0,ndet-1 ,0,ncomp-1);
  RnRxsb = dmatrix(0,ndet-1 ,0,ncomp-1);
  Rxx =    dmatrix(0,ndet-1 ,0,ndet-1) ;
  Rxxq =   dmatrix(0,ndet-1 ,0,ndet-1) ;
  Rss =    dmatrix(0,ncomp-1,0,ncomp-1);
  Rssq =   dmatrix(0,ncomp-1,0,ncomp-1);
  RnRssb = dmatrix(0,ncomp-1,0,ncomp*ndet-1);
  Pr2 =    dmatrix(0,ncomp-1,0,ncomp-1);
  AiNA =   dmatrix(0,ncomp-1,0,ncomp-1);
  Mattmp = dmatrix(0,ndet-1 ,0,ndet-1) ;
  ImDR =   dmatrix(0,ndet-1, 0,ndet-1) ;
  ACq =    dmatrix(0,ndet-1, 0,ncomp-1);
  Cq =     dmatrix(0,ncomp-1,0,ncomp-1);
  Wq =     dmatrix(0,ndet-1, 0,ncomp-1);


  //// Compute weights
  for (ib=0;ib<nbins2;ib++)
    w[ib] = (ell[ib+1] - ell[ib])*ns/fsamp;



  sprintf(testfile,"%s%s%d%s%s",outdirSpN.c_str(),"weights_",ff,termin.c_str(),".txt");
  fp = fopen(testfile,"w");
  for (jj=0;jj<nbins2;jj++){
    fprintf(fp,"%10.15g \t",w[jj]);
  }
  fprintf(fp,"\n");
  fclose(fp);



  f = fdsf(Rellexp,w,aa,P,N,ndet,ncomp,nbins2) ;
  printf("Pre em:   obj: %10.15g\n", f) ;


  for (iter=1;iter<=nbiter;iter++){

    init1D_double(iN,0,ndet,0.0);
    init1D_double(Pr,0,ncomp,0.0);
    init2D_double(Rxs,0,0,ndet,ncomp,0.0);
    init2D_double(Rxsq,0,0,ndet,ncomp,0.0);
    init2D_double(RnRxsb,0,0,ndet,ncomp,0.0);
    init2D_double(Rxx,0,0,ndet,ndet,0.0);
    init2D_double(Rxxq,0,0,ndet,ndet,0.0);
    init2D_double(Rss,0,0,ncomp,ncomp,0.0);
    init2D_double(Rssq,0,0,ncomp,ncomp,0.0);
    init2D_double(RnRssb,0,0,ncomp,ncomp*ndet,0.0);
    init2D_double(Mattmp,0,0,ndet,ndet,0.0);
    init2D_double(ImDR,0,0,ndet,ndet,0.0);
    init2D_double(ACq,0,0,ndet,ncomp,0.0);
    init2D_double(Cq,0,0,ncomp,ncomp,0.0);
    init2D_double(Wq,0,0,ndet,ncomp,0.0);

    for (ib=0;ib<nbins2;ib++){

      for (idet=0;idet<ndet;idet++)
	iN[idet] = 1.0/N[idet][ib];

      for (idet1=0;idet1<ndet;idet1++)
	for (idet2=0;idet2<ndet;idet2++)
	  Rxxq[idet1][idet2] = Rellexp[idet1*ndet + idet2][ib];

      // Robust wrt Pq=0
      for (ii=0;ii<ncomp;ii++)
	Pr[ii] = sqrt(P[ii][ib]);

      for (ii=0;ii<ncomp;ii++)
	for (jj=0;jj<ncomp;jj++)
	  Pr2[ii][jj] = Pr[ii]*Pr[jj];

      for (ii=0;ii<ncomp;ii++)
	for (jj=0;jj<ncomp;jj++){
	  AiNA[ii][jj] = 0.0;
	  for (idet=0;idet<ndet;idet++)
	    AiNA[ii][jj] += aa[idet][jj] * aa[idet][ii] * iN[idet];
	}

      for (ii=0;ii<ncomp;ii++){
	for (jj=0;jj<ncomp;jj++)
	  Cov[ii][jj] = Pr2[ii][jj] * AiNA[ii][jj];
	Cov[ii][ii] += 1.0;
      }


      // invert matrix
      dcholdc(Cov,ncomp,p);
      for (ii=0;ii<ncomp;ii++){
	for (jj=0;jj<ncomp;jj++)
	  uvec[jj] = 0.0;
	uvec[ii] = 1.0;
	dcholsl(Cov,ncomp,p,uvec,ivec);
	for (jj=0;jj<ncomp;jj++)
	  iCov[ii][jj] = ivec[jj];
      }


      for (ii=0;ii<ncomp;ii++)
	for (jj=0;jj<ncomp;jj++)
	  Cq[ii][jj] = Pr2[ii][jj] * iCov[ii][jj];

      for (idet=0;idet<ndet;idet++)
	for (ii=0;ii<ncomp;ii++){
	  Wq[idet][ii] = 0.0;
	  for (jj=0;jj<ncomp;jj++)
	    Wq[idet][ii] += aa[idet][jj] * iN[idet] *Cq[jj][ii] ;
	}
      for (idet=0;idet<ndet;idet++)
	for (ii=0;ii<ncomp;ii++){
	  Rxsq[idet][ii] = 0.0;
	  for (jj=0;jj<ndet;jj++)
	    Rxsq[idet][ii] += Rxxq[idet][jj] * Wq[jj][ii];
	}

      for (kk=0;kk<ncomp;kk++)
	for (ii=0;ii<ncomp;ii++){
	  Rssq[kk][ii] = Cq[kk][ii];
	  for (jj=0;jj<ndet;jj++)
	    Rssq[kk][ii] += Rxsq[jj][kk] * Wq[jj][ii];
	}

      for (kk=1;kk<ncomp;kk++)
	for (ii=0;ii<kk;ii++){
	  Rssq[ii][kk] = 0.5*(Rssq[ii][kk]+Rssq[kk][ii]);
	  Rssq[kk][ii] = Rssq[ii][kk];
	}


     // update power spectra
     for (ii=0;ii<ncomp;ii++)
       P[ii][ib] = abs(Rssq[ii][ii]);



     for (ii=0;ii<ndet;ii++)
       for (jj=0;jj<ncomp;jj++)
	 RnRxsb[ii][jj] += w[ib] * iN[ii]*Rxsq[ii][jj];



     for (kk=0;kk<ncomp;kk++)
       for (ii=0;ii<ndet;ii++)
	 for (jj=0;jj<ncomp;jj++)
	   RnRssb[kk][jj+ii*ncomp] = RnRssb[kk][jj+ii*ncomp] + w[ib] * iN[ii] * Rssq[kk][jj] ;

    }



    // update mixing matrix
    for (idet=0;idet<ndet;idet++){

      for (ii=0;ii<ncomp;ii++){
	uvec[ii] = RnRxsb[idet][ii];
	for (jj=0;jj<ncomp;jj++)
	  Cov[ii][jj] = RnRssb[ii][jj+idet*ncomp];
      }

      // solving the linear system
      dcholdc(Cov,ncomp,p);
      dcholsl(Cov,ncomp,p,uvec,ivec);
      for (ii=0;ii<ncomp;ii++)
	aa[idet][ii] = ivec[ii];
    }



    // EM Step with respect to N, with the new values of A and P
    for (ib=0;ib<nbins2;ib++){

      for (idet = 0;idet<ndet;idet++)
	iN[idet] = 1.0/N[idet][ib];

      for (idet1=0;idet1<ndet;idet1++)
	for (idet2=0;idet2<ndet;idet2++)
	  Rxxq[idet1][idet2] = Rellexp[idet1*ndet + idet2][ib];

      // Robust wrt Pq=0
      for (ii=0;ii<ncomp;ii++)
	Pr[ii] = sqrt(P[ii][ib]);

      for (ii=0;ii<ncomp;ii++)
	for (jj=0;jj<ncomp;jj++)
	  Pr2[ii][jj] = Pr[ii]*Pr[jj];

      for (ii=0;ii<ncomp;ii++)
	for (jj=0;jj<ncomp;jj++){
	  AiNA[ii][jj] = 0.0;
	  for (idet=0;idet<ndet;idet++)
	    AiNA[ii][jj] += aa[idet][jj] * aa[idet][ii] * iN[idet];
	}

      for (ii=0;ii<ncomp;ii++){
	for (jj=0;jj<ncomp;jj++)
	  Cov[ii][jj] = Pr2[ii][jj] * AiNA[ii][jj];
	Cov[ii][ii] += 1.0;
      }



      // invert matrix
      dcholdc(Cov,ncomp,p);
      for (ii=0;ii<ncomp;ii++){
	for (jj=0;jj<ncomp;jj++)
	  uvec[jj] = 0.0;
	uvec[ii] = 1.0;
	dcholsl(Cov,ncomp,p,uvec,ivec);
	for (jj=0;jj<ncomp;jj++)
	  iCov[ii][jj] = ivec[jj];
      }


      for (ii=0;ii<ncomp;ii++)
	for (jj=0;jj<ncomp;jj++)
	  Cq[ii][jj] = Pr2[ii][jj] * iCov[ii][jj];

      for (idet=0;idet<ndet;idet++)
	for (ii=0;ii<ncomp;ii++){
	  ACq[idet][ii] = 0.0;
	  for (jj=0;jj<ncomp;jj++)
	    ACq[idet][ii] += aa[idet][jj]*Cq[jj][ii];
	}

      for (idet=0;idet<ndet;idet++)
	for (jj=0;jj<ndet;jj++){
	  ImDR[jj][idet] = 0.0;
	  if (jj == idet)
	    ImDR[jj][idet] = 1.0;
	  for (ii=0;ii<ncomp;ii++)
	    ImDR[jj][idet] -= ACq[jj][ii]*iN[idet]*aa[idet][ii] ;
	}


      for (idet=0;idet<ndet;idet++)
	for (ii=0;ii<ndet;ii++){
	  Mattmp[idet][ii] = 0.0;
	  for (jj=0;jj<ndet;jj++)
	    Mattmp[idet][ii] += Rellexp[idet+jj*ndet][ib] * ImDR[ii][jj];
	}

      for (idet=0;idet<ndet;idet++){
	N[idet][ib] = 0.0;
	for (ii=0;ii<ndet;ii++)
	  N[idet][ib] += ImDR[idet][ii] * Mattmp[ii][idet];
      }


      for (idet=0;idet<ndet;idet++){
	for (ii=0;ii<ncomp;ii++){
	  N[idet][ib] += ACq[idet][ii]*aa[idet][ii];
	}
	N[idet][ib] = abs(N[idet][ib]);
      }


    }


    tottest = 0.0;
    for (ib=0;ib<nbins2;ib++){
      for (idet=0;idet<ndet;idet++)
	tottest += N[idet][ib]/ndet;
      for (idet=0;idet<ndet;idet++)
	if (N[idet][ib] < tottest*1e-8)
	  N[idet][ib] = tottest*1e-8;
    }




    //printf("A[1][2] =  %10.15g\n", aa[1][2]) ;
    //printf("N[2][3] =  %10.15g\n", N[2][3]) ;
    //printf("P[2][3] =  %10.15g\n", P[2][3]) ;


    // Fixing the indeterminacies.  Is it useful here?
    //rescaleAP(aa, P, ndet, ncomp, nbins2) ;



    ///// here is the problem

    f = fdsf(Rellexp,w,aa,P,N,ndet,ncomp,nbins2) ;
    printf("em->iter: %5ld   obj: %10.15g\n", iter, f) ;


  }


  // Fixing the indeterminacies.  Is it useful here?
  rescaleAP(aa, P, ndet, ncomp, nbins2) ;




  //****************************** Compute covariance matrix from the fitted model


  printf("EM step completed\n");


  for (jj=0;jj<nbins;jj++)
    for (ii=0;ii<ndet*ndet;ii++)
      Rellth[ii][jj] = 0.0;

  for (jj=0;jj<nbins2;jj++){
    for (idet=0;idet<ndet;idet++){
      Rellth[idet*ndet+idet][jj] += N[idet][jj]*SPref[jj];
      for (ii=0;ii<ndet;ii++)
	for (ll=0;ll<ncomp;ll++)
	  Rellth[idet*ndet+ii][jj] += aa[idet][ll] * aa[ii][ll] * P[ll][jj]*SPref[jj];
    }
  }


  if (nbins2 < nbins)
    for (jj=nbins2;jj<nbins;jj++)
      for (idet=0;idet<ndet;idet++)
	Rellth[idet*ndet+idet][jj] = Rellexp[idet*ndet+idet][jj]*SPref[jj];



  //*****************  Write power spectra to disk  ********************//

  sprintf(nameSpfile,"%s%s%s%s%d%s",outdirSpN.c_str(),"BoloPS",termin.c_str(),"_",ff,".psd");
  fp = fopen(nameSpfile,"w");
  for (idet1=0;idet1<ndet;idet1++){
    for (idet2=0;idet2<ndet;idet2++){

      ///// write power spectrum to disk
      tempstr1 = bolonames[idet1];
      tempstr2 = bolonames[idet2];
      //sprintf(nameSpfile,"%s%s%s%s%s%ld%s",outdirSpN.c_str(),tempstr1.c_str(),"-",tempstr2.c_str(),"_",ff,".psd");
      //fp = fopen(nameSpfile,"w");
      fprintf(fp,"%s%s%s\n",tempstr1.c_str(),"-",tempstr2.c_str());
      fprintf(fp,"%d\n",nbins);
      for (ii=0;ii<nbins;ii++){
	fprintf(fp,"%g\t",ell[ii]);
	fprintf(fp,"%10.15g\n",(Rellth[idet1*ndet+idet2][ii]+Rellth[idet2*ndet+idet1][ii])/2.0);
      }
      fprintf(fp,"%g\n",ell[nbins]);
      //fprintf(fp,"\n");
      //fclose(fp);
    }
  }
  fclose(fp);




  sprintf(nameSpfile,"%s%s%s%s%d%s",outdirSpN.c_str(),"BoloPS",termin.c_str(),"_",ff,"_exp.psd");
  fp = fopen(nameSpfile,"w");
  for (idet1=0;idet1<ndet;idet1++){
    for (idet2=0;idet2<ndet;idet2++){

      ///// write power spectrum to disk
      tempstr1 = bolonames[idet1];
      tempstr2 = bolonames[idet2];
      //sprintf(nameSpfile,"%s%s%s%s%s%ld%s",outdirSpN.c_str(),tempstr1.c_str(),"-",tempstr2.c_str(),"_",ff,"_exp.psd");
      //fp = fopen(nameSpfile,"w");
      fprintf(fp,"%s%s%s\n",tempstr1.c_str(),"-",tempstr2.c_str());
      fprintf(fp,"%d\n",nbins);
      for (ii=0;ii<nbins;ii++){
	fprintf(fp,"%g\t",ell[ii]);
	fprintf(fp,"%10.15g\n",(Rellexp[idet1*ndet+idet2][ii]+Rellexp[idet2*ndet+idet1][ii])/2.0*SPref[ii]);
      }
      fprintf(fp,"%g\n",ell[nbins]);
      //fprintf(fp,"\n");
      //fclose(fp);
    }
  }
  fclose(fp);




  sprintf(testfile,"%s%s%d%s%s",outdirSpN.c_str(),"Afinal_",ff,termin.c_str(),".txt");
  fp = fopen(testfile,"w");
  for (ii=0;ii<ndet;ii++)
    for (jj=0;jj<ncomp;jj++)
      fprintf(fp,"%10.15g \t",aa[ii][jj]);
  fprintf(fp,"\n");
  fclose(fp);



  //**************** Write component power spectra to disk
  for (idet1=0;idet1<ndet;idet1++){

    tempstr1 = bolonames[idet1];
    sprintf(nameSpfile,"%s%s%s%ld%s",outdirSpN.c_str(),tempstr1.c_str(),"_uncnoise",ff,".psd");
    fp = fopen(nameSpfile,"w");
    //fprintf(fp,"%d\n",nbins);
    for (ii=0;ii<nbins;ii++){
      //fprintf(fp,"%g\t",ell[ii]);
      fprintf(fp,"%10.15g\n",N[idet1][ii]*SPref[ii]);
    }
    //fprintf(fp,"%g\n",ell[nbins]);
    fprintf(fp,"\n");
    fclose(fp);
  }
  for (jj=0;jj<ncomp;jj++){

    sprintf(nameSpfile,"%s%s%ld%s%ld%s",outdirSpN.c_str(),"Comp_",jj,"_uncnoise",ff,".psd");
    fp = fopen(nameSpfile,"w");
    //fprintf(fp,"%d\n",nbins);
    for (ii=0;ii<nbins;ii++){
      //fprintf(fp,"%g\t",ell[ii]);
      fprintf(fp,"%10.15g\n",P[jj][ii]*SPref[ii]);
    }
    //fprintf(fp,"%g\n",ell[nbins]);
    fprintf(fp,"\n");
    fclose(fp);
  }





  delete [] data;
  delete [] data_lp;
  delete [] fdata1;
  delete [] fdata2;
  delete [] Ps;
  delete [] samptopix;
  delete [] commontmp;
  delete [] commonm_f;
  delete [] calp;
  delete [] flag;
  delete [] bfilter;
  delete [] apodwind;
  delete [] Nell;
  free_dmatrix(Rellexp,0,ndet*ndet-1,0,nbins-1);
  free_dmatrix(Rellth,0,ndet*ndet-1,0,nbins-1);
  free_dmatrix(aa,0,ndet-1,0,20);
  delete [] sign;
  free_dmatrix(Cov,0,ncomp-1,0,ncomp-1);
  free_dmatrix(iCov,0,ncomp-1,0,ncomp-1);
  free_dmatrix(iCov2,0,ncomp-1,0,ncomp-1);
  delete [] p;
  delete [] uvec;
  delete [] ivec;
  free_dmatrix(commonm,0,ncomp,0,ns-1);
  free_dmatrix(commonm2,0,ncomp,0,ns-1);
  free_dmatrix(common_f,0,ncomp,0,ns-1);
  free_dmatrix(vect,0,ncomp-1,0,ndet-1);
  delete [] Nk;
  delete [] ell;
  free_dmatrix(SpN_all,0,ndet-1,0,nbins-1);
  free_dmatrix(P,0,ncomp-1,0,nbins-1);
  free_dmatrix(N,0,ndet-1,0,nbins-1);
  delete [] iN;
  delete [] Pr;
  free_dmatrix(Rxs,0,ndet-1 ,0,ncomp-1);
  free_dmatrix(Rxsq,0,ndet-1 ,0,ncomp-1);
  free_dmatrix(RnRxsb,0,ndet-1 ,0,ncomp-1);
  free_dmatrix(Rxx,0,ndet-1 ,0,ndet-1);
  free_dmatrix(Rxxq,0,ndet-1 ,0,ndet-1);
  free_dmatrix(Rss,0,ncomp-1,0,ncomp-1);
  free_dmatrix(Rssq,0,ncomp-1,0,ncomp-1);
  free_dmatrix(RnRssb,0,ncomp-1,0,ncomp*ndet-1);
  free_dmatrix(Pr2,0,ncomp-1,0,ncomp-1);
  free_dmatrix(AiNA,0,ncomp-1,0,ncomp-1);
  free_dmatrix(Mattmp,0,ndet-1,0,ndet-1);
  free_dmatrix(ImDR,0,ndet-1,0,ndet-1);
  free_dmatrix(ACq,0,ndet-1,0,ncomp-1);
  free_dmatrix(Cq,0,ncomp-1,0,ncomp-1);
  free_dmatrix(Wq,0,ndet-1,0,ncomp-1);


}






double fdsf(double **Rellexp, double *w, double **A, double **P, double **N, long ndet, long ncomp, long nbins){

  double f;

  long ib, ii, jj, kk;
  double triRhR, detiR;

  double *p, *uvec, *ivec;
  double *Pl, *Pnl;
  double **R, **hR, **eR, **iR, **iRhR;

  p = new double[ndet];
  uvec = new double[ndet];
  ivec = new double[ndet];

  Pl   = new double[ncomp] ;
  Pnl  = new double[ndet] ;
  R  = dmatrix(0,ndet-1,0,ndet-1);
  hR = dmatrix(0,ndet-1,0,ndet-1);
  eR = dmatrix(0,ndet-1,0,ndet-1);
  iR = dmatrix(0,ndet-1,0,ndet-1);
  iRhR = dmatrix(0,ndet-1,0,ndet-1);


  init2D_double(R,0,0,ndet,ndet,0.0);
  init2D_double(hR,0,0,ndet,ndet,0.0);
  init2D_double(eR,0,0,ndet,ndet,0.0);
  init2D_double(iR,0,0,ndet,ndet,0.0);
  init1D_double(Pl,0,ncomp,0.0);
  init1D_double(Pnl,0,ndet,0.0);
  init2D_double(iRhR,0,0,ndet,ndet,0.0);


  // init
  f   = 0. ;

  for (ib=0;ib<nbins;ib++){

    //// reading Rexp
    for (ii=0;ii<ncomp;ii++)
      Pl[ii] = P[ii][ib];
    for (ii=0;ii<ndet;ii++){
      Pnl[ii] = N[ii][ib];
      for (jj=0;jj<ndet;jj++)
	hR[ii][jj] = Rellexp[ii+jj*ndet][ib];
    }

    /// building Rth from A, P, N
    for (ii=0;ii<ndet;ii++)
      for (jj=0;jj<ndet;jj++){
	R[ii][jj] = 0.0;
	if (ii == jj)
	  R[ii][jj] = Pnl[ii];
	for (kk=0;kk<ncomp;kk++)
	  R[ii][jj] += A[ii][kk] * Pl[kk] * A[jj][kk];
      }


    //printf("ib=%d\n",ib);
    /// inverting Rth
    for (ii=0;ii<ndet;ii++)
      for (jj=0;jj<ndet;jj++)
	eR[ii][jj] = R[ii][jj];
    dcholdc(eR,ndet,p);
    for (ii=0;ii<ndet;ii++){
      for (jj=0;jj<ndet;jj++)
	uvec[jj] = 0.0;
      uvec[ii] = 1.0;
      dcholsl(eR,ndet,p,uvec,ivec);
      for (jj=0;jj<ndet;jj++)
	iR[ii][jj] = ivec[jj];
    }
    // printf("ib=%d\n",ib);


    /// computing mismatch from Rexp and Rth
    for (ii=0;ii<ndet;ii++)
      for (jj=0;jj<ndet;jj++){
	iRhR[ii][jj] = 0.0;
	for (kk=0;kk<ndet;kk++)
	  iRhR[ii][jj] += iR[ii][kk]*hR[kk][jj] ;
      }

    triRhR = 0.0;
    detiR = 1.0;
    for (ii=0;ii<ndet;ii++){
      triRhR += iRhR[ii][ii];
      detiR *= 1./(p[ii]*p[ii]);
    }


  //f   +=  w[ib] * ( triRhR - log(det(iRhR)) - ndet ) ;
    f   +=  w[ib] * (triRhR - log(detiR) - ndet ) ; //pb when hR non inversible


  }


  delete [] p;
  delete [] uvec;
  delete [] ivec;
  delete [] Pl;
  delete [] Pnl;

  free_dmatrix(R,0,ndet-1,0,ndet-1);
  free_dmatrix(hR,0,ndet-1,0,ndet-1);
  free_dmatrix(eR,0,ndet-1,0,ndet-1);
  free_dmatrix(iR,0,ndet-1,0,ndet-1);
  free_dmatrix(iRhR,0,ndet-1,0,ndet-1);


  return f;

}







void rescaleAP(double **A, double **P, long ndet, long ncomp, long nbins){

  long ii, jj, ib;
  double *norm2ratio;

  norm2ratio = new double[ncomp];

  init1D_double(norm2ratio,0,ncomp,0.0);

  for (ii=0;ii<ncomp;ii++){
    for (jj=0;jj<ndet;jj++)
      norm2ratio[ii] += A[jj][ii] * A[jj][ii];
    norm2ratio[ii] = 1.0/norm2ratio[ii];
  }

  for (ii=0;ii<ncomp;ii++){
    for (jj=0;jj<ndet;jj++)
      A[jj][ii] = A[jj][ii] * sqrt(norm2ratio[ii]) ;
    for (ib=0;ib<nbins;ib++)
      P[ii][ib] = P[ii][ib] / norm2ratio[ii] ;
  }

  delete [] norm2ratio;

}
