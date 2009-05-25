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


using namespace std;




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


  long ii, idet1, idet2, ibin;

  // data parameters
  long ndet;

  //internal constants
  long nbins;

  FILE *fpw;

  char nameSpfile[200];

  double *ell1, *p, *uvec, *ivec, *SPN, *temparray;
  double **Mat_k, **iMat_k, **Rellth, **iRellth;


  string field;
  string field1;
  string field2;
  string *bolonames;
  string bolofield;
  string bolofield1;
  string bolofield2;
  string file_offsets;
  string file_frame_offsets;
  string termin;
  string noiseSpfile;
  string noiseSp_dir_output;
  string extentnoiseSp;

  list<string> channel;


  read_bolofile(argv[1], channel);
  cerr << "num ch: "<< channel.size() << endl;
  if (channel.size() == 0) {
    cerr << "Must provide at least one channel.\n\n";
  }
  ndet = channel.size();
  printf("TOTAL NUMBER OF DETECTORS: %d\n",(int)ndet);

  bolonames = new string [ndet];
  list2array(channel, bolonames);


  noiseSpfile = argv[2];
  noiseSp_dir_output = argv[3];
  extentnoiseSp = argv[4];



  //************************************************************************//
  //************************************************************************//
  //program starts here
  //************************************************************************//
  //************************************************************************//


  nbins = 500;

  
  
  
  Mat_k = dma(0,ndet-1,0,ndet-1);
  iMat_k = dma(0,ndet-1,0,ndet-1);
  Rellth  = dma(0,ndet-1,0,ndet*nbins-1);
  iRellth = dma(0,ndet-1,0,ndet*nbins-1);
  p = new double[ndet];
  uvec = new double[ndet];
  ivec = new double[ndet];
  temparray = new double[(ndet+1)*nbins + 2];
  ell1 = new double[nbins+1];
  SPN = new double[nbins];



  // Read covariance matrix from disk
  for (idet1=0;idet1<ndet;idet1++){
    field1 = bolonames[idet1];
    
    for (idet2=0;idet2<ndet;idet2++){
      
      field2 = bolonames[idet2];
      field = field1+"-"+field2;

      read_noisefile(noiseSpfile, field, ell1, SPN, &nbins);

      for (ibin=0;ibin<nbins;ibin++) Rellth[idet1][idet2*nbins + ibin] = SPN[ibin];


      //printf("%s\n",field.c_str());
      //for (ibin=0;ibin<nbins;ibin++) printf("%d, %10.15g\n", ibin, SPN[ibin]);
      //printf("\n");

    }
    printf("%d\n",idet1);

  }	
  


  



  for (ibin=0;ibin<nbins;ibin++){
    for (idet1=0;idet1<ndet;idet1++){
      for (idet2=0;idet2<ndet;idet2++){
	Mat_k[idet1][idet2] = Rellth[idet1][idet2*nbins + ibin];
	iMat_k[idet1][idet2] = 0.0;
      }
    }


    printf("%ld\n",ibin);

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
    field1 = bolonames[idet1];
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
