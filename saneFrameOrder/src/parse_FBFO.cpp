/*
 * parse_FBFO.cpp
 *
 *  Created on: 15 juin 2009
 *      Author: matthieu
 */

/*
 * testPos.cpp
 *
 *  Created on: 11 juin 2009
 *      Author: matthieu
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include "dataIO.h"
#include "mpi_architecture_builder.h"


extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#include "parse_FBFO.h"

using namespace std;

int parse_FBFO(char * ini_name, string &tmp_dir, long &ntotscan, long *&nsamples,
		std::vector<string> &fitsvect, std::vector<string> &noisevect, std::vector<int> &scans_index)
{
	dictionary	*	ini ;

	/* Some temporary variables to hold query results */
	//int				b ;
	//int				i ;
	//double			d ;
	char		*	s ;
	//long l;
	string str, dirfile;

	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	// printf dictionnary to stderr for debugging
	//iniparser_dump(ini, stderr);

	//const char *temp;


	s = iniparser_getstring(ini, "commons:data_directory", NULL);
	if(s==NULL){
		printf("You must add a line in ini file specifying a data directory : commons:data_directory\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("data_directory : [%s]\n",s);
		dirfile = s;
	}else{
		printf("You must specify a data directory : commons:data_directory\n");
		return -1 ;
	}//./RCW_120_M/

	s = iniparser_getstring(ini, "commons:output_dir",NULL);
	if(s==NULL){
		printf("Warning : The line corresponding to output directory in the ini file has been erased : commons:output_dir\n");
		cout << "Using default output directory : " << dirfile << endl;
		tmp_dir=dirfile;
	}else{
		str=(string)s;
		if(str.size()!=0){
			printf("output_dir : [%s]\n",s);
			tmp_dir=s;
		}else{
			cout << "Using default output directory : " << dirfile << endl;
			tmp_dir=dirfile;
		}//output_dir = ./RCW_120_M/ ;
	}


	s = iniparser_getstring(ini, "commons:fits_filelist",NULL);
	if(s==NULL){
		printf("You must add a line in the ini file corresponding to a frame file : commons:fits_filelist\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("fits filelist : [%s]\n",s);
		//read frame file function
		//std::vector<string> fitsvect;
		//std::vector<string> noisevect;
		//std::vector<long> scans_index;
		bool framegiven;
		read_fits_list(str, fitsvect, noisevect, scans_index, framegiven);
		//cout << "fitsvect " << fitsvect[0] << " " << fitsvect[1] << " " << fitsvect[2] << " " << fitsvect[3] << endl;
		//cout << "noisevect " << noisevect[0] << " " << noisevect[1] << " " << noisevect[2] << " " << noisevect[3] << endl;
		//cout << "scans_index " << scans_index[0] << " " << scans_index[1] << " " << scans_index[2] << " " << scans_index[3] << endl;

		for(int ii=0;ii<(int)fitsvect.size();ii++){
			cout << dirfile + fitsvect[ii] << endl;
			fitsvect[ii] = dirfile + fitsvect[ii];}

		readFrames(fitsvect, nsamples);
		ntotscan = fitsvect.size();

		//getchar();
	}else{
		printf("You must specify a fits filelist : commons:fits_filelist\n");
		return -1 ;
	}//frame_file =./RCW_120_M/fits_files.txt ;

	/*
	s = iniparser_getstring(ini, "commons:frame_file",NULL);
	if(s!=NULL){
		printf("frame file : [%s]\n",s);
		//read frame file function
		std::vector<string> dummy;
		read_strings((string)s,dummy);
		if(((int)dummy.size())==0){
			printf("You must provide one number of samples per scan !");
			return(-1);}
		vector<string>::iterator it, it2;
		//it=dummy.begin();
		int ind=1;
		it2=dummy.end();

		for(it=dummy.begin();it<it2;it++){
			if(ind%2==1){
				//cout << "ind : "<< ind/2 << endl;
				//cout << "it : " << (*it) << endl;
				//cout << "frames" << endl;
				temp=(*it).c_str();
				fframes_vec.push_back(atol(temp));
			}else{
				//cout << "ind : "<< ind/2 << endl;
				//cout << "it : " << (*it) << endl;
				//cout << "nsamples" << endl;
				temp=(*it).c_str();
				l=atol(temp);
				l -= fframes_vec.back() - 1;
				nsamples_vec.push_back(l);
			}

			ind++;
		}


	}else{
		printf("You must specify a frame file : commons:frame_file\n");
		return(-1);
	}//frame_file =./RCW_120_M/frame_file.txt ;
	 */

	/*char * pPath;
	pPath = getenv ("TMPBATCH");
	if (pPath!=NULL){
	  printf ("The current path is: %s\n",pPath);
	  tmp_dir=pPath;
	}else{
	  s = iniparser_getstring(ini, "commons:output_dir",NULL);
	  if(s!=NULL){
	    printf("temp_dir : [%s]\n",s);
	    tmp_dir=s;
	  }else{
	    printf("You must provide an output directory to write parallel_scheme file : commons:output_dir\n");
	    return(-1);
	  }
	  }*/

	//fname = tmp_dir + parallel_scheme_filename;

	// Check improper usage
	/*if (fframes_vec.size() != nsamples_vec.size()) {
		cerr << "Must give at least one first frame number. Exiting.\n";
		return(-1);
	}*/

	ntotscan = fitsvect.size();

	iniparser_freedict(ini);

	return 0 ;
}


