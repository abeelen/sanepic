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


#include "parse_FBFO.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

int parse_FBFO(char * ini_name, string &fname, long &ntotscan, std::vector<long> &fframes_vec, std::vector<long> &nsamples_vec)
{
	dictionary	*	ini ;

	/* Some temporary variables to hold query results */
	//int				b ;
	//int				i ;
	//double			d ;
	char		*	s ;
	long l;
	string str;

	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	// printf dictionnary to stderr for debugging
	//iniparser_dump(ini, stderr);

	const char *temp;

	s = iniparser_getstring(ini, "commons:frame_file",NULL);
	if(s!=NULL){
		printf("frame file : [%s]\n",s);
		//read frame file function
		std::vector<string> dummy;
		read_bolofile((string)s,dummy);
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


	s = iniparser_getstring(ini, "sanepic_parallel_scheme:fname",NULL);
	if(s!=NULL){
		printf("file_name : [%s]\n",s);
		fname=s;
	}else{
		printf("You must provide a file name for the parallel scheme data file : sanepic_parallel_scheme:file_name\n");
		return(-1);
	}

	// Check improper usage
	if (fframes_vec.size() != nsamples_vec.size()) {
		cerr << "Must give at least one first frame number. Exiting.\n";
		return(-1);
	}

	ntotscan = fframes_vec.size();

	iniparser_freedict(ini);

	return 0 ;
}


