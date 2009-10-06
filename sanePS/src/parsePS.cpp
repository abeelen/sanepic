/*
 * parsePS.cpp
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */

#include "inputFileIO.h"
#include "parsePS.h"
#include "mpi_architecture_builder.h"

int parse_sanePS_ini_file(char * ini_name, struct user_options &u_opt,
		long &ntotscan, long &ndet,
		std::vector<string> &bolonames,long *&fframes, long *&nsamples, std::vector<string> &extentnoiseSP, string &MixMatfile,string &signame,
		std::vector<string> &fitsvect,std::vector<string> &noisevect, std::vector<long> &scans_index)
{
	dictionary	*	ini ;

	/* Some temporary variables to hold query results */
	int				b ;
	int				i ;
	double			d ;
	char		*	s ;
	//long l;
	string str;
	//const char *temp;

	// for poutdir default value
	/*char * pPath;
	pPath = getenv ("TMPBATCH");
	if (pPath!=NULL)
		printf ("The current path is: %s",pPath);*/



	int nnf; // extentnoiseSp_list number of elements




	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	// printf dictionnary to stderr for debugging
	//iniparser_dump(ini, stderr);


	/*#ifdef USE_MPI
	s = iniparser_getstring(ini, "sanepic_parallel_scheme:fname",NULL);
	if(s==NULL){
		printf("You must add a line in ini file specifying Find_best_frame_order result : sanepic_parallel_scheme:fname\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("fname : [%s]\n",s);
		fname=s;
	}else{
		printf("You need to run Find_best_frame_order first and specify the generated file path and name !\n");
		return -1 ;
	}
#endif*/


	i = iniparser_getint(ini, "commons:time_offset", 0);
	//if(isnan((double)i)){
	if(i!=0){
		printf("time_offset :      [%d]\n", i);
		u_opt.shift_data_to_point=i;
	}//time_offset =  ;*/

	i = iniparser_getint(ini, "commons:apodize_Nsamples", -1);
	if(i>0){
		printf("apodize_Nsamples :      [%d]\n", i);
		u_opt.napod=i;
	}else{
		printf("You must choose a number of samples to apodize commons:apodize_Nsamples\n");
		return -1 ;
	}//apodize_Nsamples = 100 ;

	d = iniparser_getdouble(ini,(char*)"sanepic_preprocess:sampling_frequency", -1.0);
	if (d<=0.0){
		printf("sampling_frequency cannot be negative or 0 ! Or maybe you forgot to mention sampling frequency \n");
		return -1 ;
	}else{
		printf("sampling_frequency  :   [%g]\n", d);
		u_opt.fsamp=d;
	}//sampling_frequency = 25.0 ;


	b = iniparser_getboolean(ini, "commons:nofill_gap", -1);
	if(b!=-1){
		printf("nofill_gap:    [%d]\n", b);
		u_opt.NOFILLGAP=b;
	}	//NOFILLGAP = 0 ;

	b = iniparser_getboolean(ini, "sanepic_preprocess:no_baseline", -1);
	if(b!=-1){
		printf("no_baseline:    [%d]\n", b);
		u_opt.NORMLIN=b;
	}	//NORMLIN = False ;

	b = iniparser_getboolean(ini, "sanepic_preprocess:remove_poly", -1);
	if(b!=-1){
		printf("remove_poly:    [%d]\n", b);
		u_opt.remove_polynomia=b;
	} //remove_polynomia= True;


	b = iniparser_getboolean(ini, "commons:map_flagged_data", -1);
	if(b!=-1){
		printf("map_flagged_data:    [%d]\n", b);
		u_opt.flgdupl=b;
	}
	//flgdupl = False ;

	s = iniparser_getstring(ini, "commons:data_directory", NULL);
	if(s==NULL){
		printf("You must add a line in ini file specifying a data directory : commons:data_directory\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("data_directory : [%s]\n",s);
		u_opt.dirfile = s;
	}else{
		printf("You must specify a data directory : commons:data_directory\n");
		return -1 ;
	}//./RCW_120_M/

	s = iniparser_getstring(ini, "commons:output_dir",NULL);
	if(s==NULL){
		printf("Warning : The line corresponding to output directory in the ini file has been erased : commons:output_dir\n");
		cout << "Using default output directory : " << u_opt.dirfile << endl;
		u_opt.outdir=u_opt.dirfile;
	}else{
		str=(string)s;
		if(str.size()!=0){
			printf("output_dir : [%s]\n",s);
			u_opt.outdir=s;
		}else{
			cout << "Using default output directory : " << u_opt.dirfile << endl;
			u_opt.outdir=u_opt.dirfile;
		}//output_dir = ./RCW_120_M/ ;
	}



	// for poutdir default value
	char * pPath;
	pPath = getenv ("TMPBATCH");
	if (pPath!=NULL){
		u_opt.tmp_dir=pPath;
		printf ("The current path is: %s\n",pPath);
	}else{
		s = iniparser_getstring(ini, "commons:temp_dir",NULL);
		if(s==NULL){
			printf("Warning : The line corresponding to temporary directory in the ini file has been erased : commons:output_dir\n");
			cout << "Using default temporary directory : " << u_opt.dirfile << endl;
			u_opt.tmp_dir=u_opt.dirfile;
		}else{
			str=(string)s;
			if(str.size()!=0){
				printf("temp_dir : [%s]\n",s);
				u_opt.tmp_dir=s;
			}else{
				cout << "Using default temporary directory : " << u_opt.dirfile << endl;
				u_opt.tmp_dir=u_opt.dirfile;
			}//tmp_dir = ./internal_sanepic/ ;
		}

	}


/*	s = iniparser_getstring(ini, "commons:bolofield_extension",NULL);
	if(s==NULL){
		printf("You must choose add a line in the ini_file corresponding to the bolofield extension : commons:bolofield_extension\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("bolofield extension : [%s]\n",s);
		bextension=s;
	}else{
		printf("You must specify a bolo extension : commons:bolofield_extension\n");
		return -1 ;
	}//_data

	s = iniparser_getstring(ini, "commons:flag_field_extension",NULL);
	if(s==NULL){
		printf("You must add a line corresponding to flag_field_extension in the parser file : commons:flag_field_extension\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("flag_field_extension : [%s]\n",s);
		fextension=s;
	}else{
		printf("you must specify flag_field_extension\n");
		return -1 ;
	}//flag_field_extension = _flag ;
*/


	/*s = iniparser_getstring(ini, "commons:tmp_dir",NULL);
	if(s==NULL){
		printf("You must add a line corresponding to noise_suffixe in the ini file : commons:noise_suffixe\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("noise_suffixe : [%s]\n",s);
		noiseSppreffile=s;
	}else{
		printf("you must specify a noise_suffixe\n");
		return -1 ;
	}//noise_suffixe = ./RCW_120_M/ ;*/
	u_opt.noiseSppreffile=u_opt.tmp_dir;

	s = iniparser_getstring(ini, "commons:channel",NULL);
	if(s==NULL){
		printf("You must add a line in ini file corresponding to a bolometer file : commons:channel\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("channel file : [%s]\n",s);
		read_strings((string)s, bolonames);
	}else{
		printf("You must specify a bolometer file : commons:channel\n");
		return -1 ;
	}//	channel =./RCW_120_M/bolos_commons.txt ;

/*
	s = iniparser_getstring(ini, "commons:frame_file",NULL);
	if(s==NULL){
		printf("You must add a line in the ini file corresponding to a frame file : commons:frame_file\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("frame file : [%s]\n",s);
		//read frame file function
		std::vector<string> dummy;
		read_strings((string)s,dummy);
		if(((int)dummy.size())==0){
			printf("You must provide one number of samples per scan !");
			exit(0);}
		vector<string>::iterator it, it2;
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
		return -1 ;
	}//frame_file =./RCW_120_M/frame_file.txt ;

*/


	s = iniparser_getstring(ini, "commons:fits_filelist",NULL);
	if(s==NULL){
		printf("You must add a line in the ini file corresponding to a frame file : commons:fits_filelist\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("fits filelist : [%s]\n",s);
		//read frame file function
		//std::vector<string> fitsvect_path;
		std::vector<string> noisevect;
		//std::vector<long> scans_index;
		bool framegiven;

		read_fits_list(str, fitsvect, noisevect, scans_index, framegiven);
		//cout << "fitsvect " << fitsvect[0] << " " << fitsvect[1] << " " << fitsvect[2] << " " << fitsvect[3] << endl;
		//cout << "noisevect " << noisevect[0] << " " << noisevect[1] << " " << noisevect[2] << " " << noisevect[3] << endl;
		//cout << "scans_index " << scans_index[0] << " " << scans_index[1] << " " << scans_index[2] << " " << scans_index[3] << endl;
		for(int ii=0;ii<(int)fitsvect.size();ii++){
			cout << u_opt.dirfile + fitsvect[ii] << endl;
			fitsvect[ii] = u_opt.dirfile + fitsvect[ii];}

		readFrames( &ntotscan , fitsvect, fframes, nsamples);

		//getchar();
	}else{
		printf("You must specify a fits filelist : commons:fits_filelist\n");
		return -1 ;
	}//frame_file =./RCW_120_M/fits_files.txt ;




	s = iniparser_getstring(ini, "commons:noise_prefixe_file",NULL);
	if(s==NULL){
		printf("You must add a line corresponding to noise_prefixe file in the ini file : commons:noise_prefixe\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("noise_prefixe_file : [%s]\n",s);
		//extentnoiseSP.push_back(s);
		std::vector<string> dummy3;
		read_strings(str,dummy3);
		if(((int)dummy3.size())==0){
			printf("You must provide at least one noise prefixe (or one per scan) in noise_prefixe_file !\n");
			return -1;}
		for(int ii=0; ii<(int)dummy3.size(); ii++)
			extentnoiseSP.push_back(dummy3[ii].c_str());
	}else{
		printf("you must specify noise_prefixe_file\n");
		return -1;
	}//noise_prefixe = NoisePS ;


	s = iniparser_getstring(ini, (char*)"sanepic_estim_PS:noise_estim", NULL);
	if(s==NULL){
		printf("You must add a line corresponding to the mixing matrix of noise components in the ini file : sanepic_estim_PS:noise_estim\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("noise_estim :      [%s]\n", s);
		MixMatfile=s;

	}else{
		printf("You must give filename containing the mixing matrix of noise components : noise_estim\n");
		return(-1);
	}// MixMatfile = Mixlaboca


	s = iniparser_getstring(ini, (char*)"sanepic_estim_PS:map_file", NULL);
	if(s==NULL){
		printf("You must add a line corresponding to the fits map file (created by sanePic) in the ini file : sanepic_estim_PS:map_file\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("noise_estim :      [%s]\n", s);
		signame=s;

	}else{
		signame="NOSIGFILE";
	}//



	// path in which data are written
	/*if (pPath != NULL){
		poutdir = pPath;
	} else {
		poutdir = outdir;
	}*/



	// Set default parameter values
	/*if (fframes_vec.size() == 0) {
		fframes_vec.push_back(0);
		nsamples_vec.push_back(-1);
	}
	if (fframes_vec.size() == 1 && nsamples_vec.size() == 0)
		nsamples_vec.push_back(-1);*/

	// Check improper usage
	if (bolonames.size() == 0) {
		cerr << "Must provide at least one channel.\n\n";
		exit(1);
		//usage(argv[0]);
	}
/*	if (fframes_vec.size() != nsamples_vec.size()) {
		cerr << "Must give at least one first frame number. Exiting.\n";
		exit(1);
	}*/


	//ntotscan = fframes_vec.size();
	ndet = bolonames.size();


	nnf = (int)extentnoiseSP.size();
	//nnf=1; // Temporarily
	if (nnf != 1 && nnf != ntotscan){
		cerr << "ERROR: There should be one noise power spectrum file per scan, or a single one for all the scans. Check -K options" << endl;
		exit(1);
	}
	if (nnf == 1 && ntotscan > 1)
	extentnoiseSP.resize(ntotscan, extentnoiseSP[0]);

	//  printf("%d\n",nnf);



	iniparser_freedict(ini);
	return 0 ;
}
