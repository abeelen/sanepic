#include <iostream>
#include <iomanip>
#include <fstream>
#include "todprocess.h"
#include "map_making.h"
#include <time.h>
#include <fftw3.h>
// #include <getdata.h>
//#include "/global/software/pgi-6.0/linux86/6.0/include/mpi.h"
#include <fcntl.h>
#include <unistd.h>
#include <vector>
#include <list>


using namespace std;




void read_bolofile(string fname, vector<string>& bolos) {
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




void read_noisefile(string fname, string bolo1bolo2, double *ell, double *SPN, long *nbins){
  char buff[256];
  string line;

  long ii;
  double dummy1, dummy2;

  ifstream Spfile (fname.c_str());
  if (!Spfile.is_open()){
    cerr << "Error opening Noise power spectra file '" << fname << "'Exiting.\n";
    exit(1);
  }

  while(! Spfile.eof()) {
    Spfile.getline(buff,255);
    line = buff;

    line.erase(0, line.find_first_not_of(" \t"));       // remove leading white space
    if (line.empty() || line[0] == '#') continue;       // skip if empty or commented
    line = line.substr(0, line.find_first_of(" \t"));   // pick out first word

    if (line == bolo1bolo2){
      Spfile.getline(buff,255);
      line = buff;
      *nbins = atoi(line.c_str());
      for (ii=0;ii<*nbins;ii++){
	Spfile.getline(buff,255);
	line = buff;
	sscanf(line.c_str(), "%lf%lf", &dummy1, &dummy2);
	ell[ii] = dummy1;
	SPN[ii] = dummy2;
      }
      Spfile.getline(buff,255);
      line = buff;
      ell[*nbins] = atof(line.c_str());

    }

  }

  Spfile.close();

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


  long ii, jj;
  long idet1, idet2, ibin;

  // data parameters
  long ndet, ndetOrig, nbins;

  FILE *fpw;

  char nameSpfile[200];

  double *ell1, *p, *uvec, *ivec, *temparray;
  double **Mat_k, **iMat_k, **Rellth, **RellthOrig, **iRellth;


  string field;
  string field1;
  string field2;
  string *bolonamesIn;
  string bolofield;
  string bolofield1;
  string bolofield2;
  string file_offsets;
  string file_frame_offsets;
  string termin;
  string noiseSpfile;
  string noiseSp_dir_output;
  string extentnoiseSp;

  vector<string> channelIn;
  vector<string> channelOut;
  vector<int> indexIn;


  // Input argument c: hannelList ....
  read_bolofile(argv[1], channelIn);
  //  cerr << "num ch: "<< channel.size() << endl;
  if (channelIn.size() == 0) {
    cerr << "Must provide a valid input channel list.\n\n";
  }
  printf("TOTAL NUMBER OF DETECTORS IN PS file: %d\n",(int)channelIn.size());

  // ... Corresponding NoisePS file
  noiseSpfile = argv[2];

  // Read input PS file
  ifstream Spfile(argv[2], ios::binary);
  if (!Spfile.is_open()){
    cerr << "Could not open PS file" << endl;
    exit(1);
  }

  Spfile.read(reinterpret_cast<char *>(&ndetOrig),sizeof(long));
  Spfile.read(reinterpret_cast<char *>(&nbins),sizeof(long));

  nbins = (int) nbins;

  if (ndetOrig != channelIn.size()) {
    cerr << "Input Channel List does not correspond to PS file\n";
  }

  //  printf("%d bins\n",(int)nbins);

  ell1 = new double[nbins+1];
  RellthOrig  = dmatrix(0,ndetOrig*ndetOrig-1,0,nbins-1);

  for (ii=0; ii<nbins+1; ii++)
    Spfile.read(reinterpret_cast<char *>(&ell1[ii]), sizeof(double));

  for (idet1 = 0;   idet1 < ndetOrig; idet1++)
    for (idet2 = 0; idet2 < ndetOrig; idet2++)
      for ( ii=0; ii<nbins; ii++)
	Spfile.read(reinterpret_cast<char *>(&RellthOrig[idet1*ndetOrig+idet2][ii]),sizeof(double));

  Spfile.close();

  // Input argument for output : channelList
  read_bolofile(argv[3],channelOut);
  if (channelOut.size() == 0) {
    cerr << "Must provide a valid output channel list.\n\n";
  }
  ndet = channelOut.size();
  printf("TOTAL NUMBER OF DETECTORS TO OUTPUT : %d\n",(int)ndet);

  noiseSp_dir_output = argv[4];
  extentnoiseSp = argv[5];

  //************************************************************************//
  //************************************************************************//
  //program starts here
  //************************************************************************//
  //************************************************************************//


  // Take subsample of the RellthOrig
  Rellth = dmatrix(0,ndet-1,0,ndet*nbins-1);

  // find indexes of input bolo file corresponding to the output bolo file
  indexIn.resize(channelOut.size(),-1);
  for ( idet1 = 0; idet1 < ndet; idet1++) {
    for ( idet2 = 0; idet2 < ndetOrig; idet2++) {
      if (channelOut[idet1] == channelIn[idet2])
	indexIn[idet1] = idet2;
    }
  }

//   for (idet1=0; idet1<ndet; idet1++)
//     cout << idet1 << " " << channelOut[idet1] << " " << indexIn[idet1] << endl;

  for ( idet1 = 0; idet1 < ndet; idet1++) {
    for ( idet2 = 0; idet2 < ndet; idet2++) {
      for ( ii = 0; ii < nbins; ii++) {
	//	cout << idet1 << " " << idet2 << " " << ii << endl;
	Rellth[idet1][idet2*nbins+ii] = (RellthOrig[indexIn[idet1]*ndetOrig+indexIn[idet2]][ii]+RellthOrig[indexIn[idet2]*ndetOrig+indexIn[idet1]][ii])/2;
      }
    }
  }


  Mat_k = dmatrix(0,ndet-1,0,ndet-1);
  iMat_k = dmatrix(0,ndet-1,0,ndet-1);
  iRellth = dmatrix(0,ndet-1,0,ndet*nbins-1);
  p = new double[ndet];
  uvec = new double[ndet];
  ivec = new double[ndet];
  temparray = new double[(ndet+1)*nbins + 2];
  //  SPN = new double[nbins];


//   // Read covariance matrix from disk
//   for (idet1=0;idet1<ndet;idet1++){
//     field1 = bolonames[idet1];

//     for (idet2=0;idet2<ndet;idet2++){

//       field2 = bolonames[idet2];
//       field = field1+"-"+field2;

//       read_noisefile(noiseSpfile, field, ell1, SPN, &nbins);

//       for (ibin=0;ibin<nbins;ibin++) Rellth[idet1][idet2*nbins + ibin] = SPN[ibin];


//       //printf("%s\n",field.c_str());
//       //for (ibin=0;ibin<nbins;ibin++) printf("%d, %10.15g\n", ibin, SPN[ibin]);
//       //printf("\n");

//     }
//     cout << "Progress : "  << setiosflags(ios::fixed) << setprecision(4) << idet1*100./ndet << " %\r" << flush ;

//   }








  for (ibin=0;ibin<nbins;ibin++){

    cout << "Progress : " << ibin*100./nbins << " %\r" << flush;

    for (idet1=0;idet1<ndet;idet1++){
      for (idet2=0;idet2<ndet;idet2++){
	Mat_k[idet1][idet2] = Rellth[idet1][idet2*nbins + ibin];
	iMat_k[idet1][idet2] = 0.0;
      }
    }


    for (idet1=0;idet1<ndet;idet1++){
      for (idet2=0;idet2<ndet;idet2++){
	if (Mat_k[idet1][idet2] > sqrt(Mat_k[idet1][idet1]*Mat_k[idet2][idet2]))
	  printf("ERROR %10.15g\t",Mat_k[idet1][idet2]);
      }
    }



    // invert the covariance matrix per mode
    dcholdc(Mat_k,ndet,p);
    for (idet1=0;idet1<ndet;idet1++){
      for (idet2=0;idet2<ndet;idet2++)
    	uvec[idet2] = 0.0;
      uvec[idet1] = 1.0;
      dcholsl(Mat_k,ndet,p,uvec,ivec);
      for (idet2=0;idet2<ndet;idet2++)
    	iMat_k[idet1][idet2] = ivec[idet2];
    }

    for (idet1=0;idet1<ndet;idet1++)
      for (idet2=0;idet2<ndet;idet2++)
	iRellth[idet1][idet2*nbins + ibin] = iMat_k[idet1][idet2];

  }


  /// write inverse covariance matrix to disk


  temparray[0] = double(nbins);
  for (ii=0;ii<=nbins;ii++)
    temparray[ii+1] = ell1[ii];


  for (idet1=0;idet1<ndet;idet1++){
    field1 = channelOut[idet1];
    sprintf(nameSpfile,"%s%s%s%s",noiseSp_dir_output.c_str(),field1.c_str(),"-all_Inv",extentnoiseSp.c_str());
    for (ii=0;ii<nbins*ndet;ii++)
      temparray[ii+nbins+2] = iRellth[idet1][ii];
    fpw = fopen(nameSpfile,"w");
    fwrite(temparray,sizeof(double), ndet*nbins+nbins+2, fpw);
    fclose(fpw);
  }



  //******************************************************************//
  //******************************************************************//
  //************************  End of loop ****************************//
  //******************************************************************//
  //******************************************************************//



  return 0;
}