/*
 * parseSanepic.cpp
 *
 *  Created on: 18 juin 2009
 *      Author: matthieu
 */

#include "parseSanepic.h"
#include "mpi_architecture_builder.h"
#include "inputFileIO.h"

int parse_sanePic_ini_file(char * ini_name, double &pixdeg, int  &shift_data_to_point, long  &napod,double &fsamp, bool &NOFILLGAP,bool &NORMLIN,bool &projgaps,bool &remove_polynomia, bool &flgdupl,
		bool &CORRon, int &iterw, long &ntotscan, long &ndet, double &f_lp, double &f_lp_Nk, string &dirfile, string &outdir, string &tmp_dir,
		string &termin,	int &coordsyst, string &MixMatfile, std::vector<string> &bolonames,long *&fframes,long *&nsamples, string &fname,
		std::vector<long> &xxi, std::vector<long> &xxf, std::vector<long> &yyi, std::vector<long> &yyf, std::vector<double> &fcut, std::vector<string> &extentnoiseSP,std::vector<string> &fitsvect,std::vector<string> &noisevect, std::vector<long> &scans_index)
{

	dictionary	*	ini ;

	/* Some temporary variables to hold query results */
	int				b ;
	int				i ;
	double			d ;
	char		*	s ;
	//long l;
	string str;
	const char *temp;


	/*char * pPath;
	pPath = getenv ("TMPBATCH");
	if (pPath!=NULL)
		printf ("The current path is: %s",pPath);*/






	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	// printf dictionnary to stderr for debugging
	//iniparser_dump(ini, stderr);

	printf("\nsanepic_conjugate_gradient:\n");

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
		return -1;
	}
#endif*/

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
			cout << dirfile + fitsvect[ii] << endl;
			fitsvect[ii] = dirfile + fitsvect[ii];}

		readFrames( &ntotscan , fitsvect, fframes, nsamples);

		//getchar();
	}else{
		printf("You must specify a fits filelist : commons:fits_filelist\n");
		return -1 ;
	}//frame_file =./RCW_120_M/fits_files.txt ;

	/*
	s = iniparser_getstring(ini, "commons:frame_file",NULL);
	if(s==NULL){
		printf("You must add a line in ini file corresponding to a bolometer file : commons:channel\n");
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
		return -1 ;
	}//frame_file =./RCW_120_M/frame_file.txt ;
	 */

	i = iniparser_getint(ini, "commons:coord_syst", -1);
	if((i==1)||(i==2)||(i==3)){
		printf("Coordinate system :      [%d]\n", i);
		coordsyst=i;
	}else{
		printf("Choose a coordinate system between 1 and 3 : commons:coord_syst\n");
		return -1 ;
	}//coord_syst = 1 ;

	/*
	s = iniparser_getstring(ini, "commons:bolofield_extension",NULL);
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
	 */

	s = iniparser_getstring(ini, "sanepic_compute_positions:pixsize",NULL);
	if(s==NULL){
		printf("You must add a line in the ini file corresponding to pixel size : sanepic_compute_positions:pixsize\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		temp=str.c_str();
		pixdeg=atof(temp);
		cout << "Pixsize : " << setprecision(str.size()) << pixdeg << endl;
		//printf("Pixsize :   [%lf]\n", pixdeg);

	}else{
		printf("Pixsize cannot be negative ! or you forgot to mention pixel size\n");
		return -1 ;
	}
	/*d = iniparser_getdouble(ini,(char*)"sanepic_compute_positions:pixsize", -1.0); // pb de precision
	if (d<=0){
		printf("Pixsize cannot be negative ! or you forgot to mention pixel size\n");
		return -1 ;
	}else{
		printf("Pixsize :   [%f]\n", d);
		pixdeg=d;
	}//pixsize =  0.00168725828819 ;*/

	i = iniparser_getint(ini, "commons:apodize_Nsamples", -1);
	if(i>0){
		printf("apodize_Nsamples :      [%d]\n", i);
		napod=i;
	}else{
		printf("You must choose a number of samples to apodize commons:apodize_Nsamples\n");
		return -1 ;
	}//apodize_Nsamples = 100 ;


	b = iniparser_getboolean(ini, "commons:nofill_gap", -1);
	if(b!=-1){
		printf("nofill_gap:    [%d]\n", b);
		NOFILLGAP=b;
	}
	//NOFILLGAP = 0 ;

	d = iniparser_getdouble(ini,(char*)"sanepic_preprocess:sampling_frequency", -1.0);
	if (d<=0.0){
		printf("sampling_frequency cannot be negative or 0 ! Or maybe you forgot to mention sampling frequency \n");
		return -1 ;
	}else{
		printf("sampling_frequency  :   [%g]\n", d);
		fsamp=d;
	}//sampling_frequency = 25.0 ;

	d = iniparser_getdouble(ini,(char*)"sanepic_preprocess:filter_frequency", -1.0);
	if (d<=0.0){
		printf("filter_frequency cannot be negative ! or maybe you have to mention filter frequency \n");
		return -1 ;
	}else{
		printf("filter_frequency  :   [%g]\n", d);
		f_lp=d;
	}//filter_frequency = 0.005 ;

	// a changer
	/*d = iniparser_getdouble(ini,(char*)"sanepic_preprocess:noise_cut_frequency", 0.0);
	if (d<0.0){
		printf("noise_cut_frequency cannot be negative !\n");
		return -1 ;
	}else{
		printf("noise_cut_frequency  :   [%g]\n", d);
		if(d==0.0){
			f_lp_Nk=f_lp;
			//fcut.push_back(f_lp);
			//cout << "merde\n";
		}else{
			//int ll=0;
			//while(d[ll]!='\0')
			//fcut.push_back(d);
			f_lp_Nk=d;
		}//noise_cut_frequency = 0.01 ;
	}*/

	s = iniparser_getstring(ini, "sanepic_preprocess:fcut_file",NULL);
	if(s==NULL){
		printf("You must add a line corresponding to noise cut frequency file in the parser file : sanepic_preprocess:fcut_file\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("fcut file : [%s]\n",s);
		//read frame file function
		std::vector<string> dummy2;
		read_strings(str,dummy2);
		if(((int)dummy2.size())==0){
			printf("You must provide at least one number of noise cut frequency (or one per scan) in fcut_file !\n");
			return -1;}
		for(int ii=0; ii<(int)dummy2.size(); ii++)
			fcut.push_back(atof(dummy2[ii].c_str()));

	}else{
		printf("You must specify a noise cut frequency file : sanepic_preprocess:noise_cut_frequency\n");
		return -1;
	}//frame_file =./RCW_120_M/fcut_file.txt ;


	/*
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

	s = iniparser_getstring(ini, "commons:pointing_field_extension",NULL);
	if(s==NULL){
		printf("You must add a line corresponding to pointing_field_extension in the parser file : commons:pointing_field_extension\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("pointing_field_extension : [%s]\n",s);
		pextension=s;
	}else{
		printf("you must specify pointing_field_extension\n");
		return -1 ;
	}//pointing_field_extension = _def ;
	 */

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

	/*char * pPath;
	pPath = getenv ("TMPBATCH");
	if (pPath!=NULL){
		noiseSppreffile=pPath;
		printf ("The current path is: %s\n",pPath);
	}else{
		s = iniparser_getstring(ini, "commons:tmp_dir",NULL);
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
		}//noise_suffixe = ./RCW_120_M/ ;
	}
	 */
	s = iniparser_getstring(ini, "commons:output_dir",NULL);
	if(s==NULL){
		printf("Warning : The line corresponding to output directory in the ini file has been erased : commons:output_dir\n");
		cout << "Using default output directory : " << dirfile << endl;
		outdir=dirfile;
	}else{
		str=(string)s;
		if(str.size()!=0){
			printf("output_dir : [%s]\n",s);
			outdir=s;
		}else{
			cout << "Using default output directory : " << dirfile << endl;
			outdir=dirfile;
		}//output_dir = ./RCW_120_M/ ;
	}

	//#ifdef USE_MPI
	// for poutdir default value
	char * pPath;
	pPath = getenv ("TMPBATCH");
	pPath=NULL; // TODO : Remove apres test !!!!
	if (pPath!=NULL){
		tmp_dir=pPath;
		printf ("The current path is: %s\n",pPath);
	}else{
		//#else

		s = iniparser_getstring(ini, "commons:temp_dir",NULL);
		if(s==NULL){
			printf("Warning : The line corresponding to temporary directory in the ini file has been erased : commons:output_dir\n");
			cout << "Using default output directory : " << dirfile << endl;
			tmp_dir=dirfile;
		}else{
			str=(string)s;
			if(str.size()!=0){
				printf("temp_dir : [%s]\n",s);
				tmp_dir=s;
			}else{
				cout << "Using default output directory : " << dirfile << endl;
				tmp_dir=dirfile;
			}//output_dir = ./RCW_120_M/ ;
		}
		//#endif
	}

	s = iniparser_getstring(ini, "commons:out_file_str",NULL);
	if(s==NULL){
		printf("You must add a line corresponding to a prefixe for generated files in the ini file : commons:out_file_str\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("out_file_str : [%s]\n",s);
		termin=s;
	}else{
		printf("You must specify a prefixe for generated files : out_file_str\n");
		return -1;
	}//out_file_str = sanepic ;


	i = iniparser_getint(ini, "commons:time_offset", 0);
	//if(isnan((double)i)){
	if(i!=0){
		printf("time_offset :      [%d]\n", i);
		shift_data_to_point=i;
	}//time_offset =  ;*/



	// crossing constraint removal box coordinates
	s = iniparser_getstring(ini, (char*)"commons:box_coord_x1", NULL); // faire un data file
	if(s==NULL){
		printf("Warning : The line corresponding to box_coord_x1 in the ini file has been erased : commons:box_coord_x1\n");
	}else{
		str=(string)s;
		if(str.size()!=0){
			printf("box_coord_x1 :      [%ld]\n", atol(s));

			if (xxi.size() != yyf.size()){
				printf("box_coord_x1 requires at least box_coord_x2, box_coord_y1, box_coord_y2. Exiting.\n");
				return -1;
			}
			xxi.push_back(atol(s));//box_coord_x1 =  ;
		}
	}

	s = iniparser_getstring(ini, (char*)"commons:box_coord_x2", NULL);
	if(s==NULL){
		printf("Warning : The line corresponding to box_coord_x2 in the ini file has been erased : commons:box_coord_x2\n");
	}else{
		str=(string)s;
		if(str.size()!=0){
			printf("box_coord_x2 :      [%ld]\n", atol(s));
			if (xxf.size() != xxi.size()-1){
				printf("box_coord_x2 requires at least box_coord_x1, box_coord_y1, box_coord_y2. Exiting. \n");
				return -1;
			}
			xxf.push_back(atol(s));//box_coord_x2 =  ;
		}
	}

	s = iniparser_getstring(ini, (char*)"commons:box_coord_y1", NULL);
	if(s==NULL){
		printf("Warning : The line corresponding to box_coord_y1 in the ini file has been erased : commons:box_coord_y1\n");
	}else{
		str=(string)s;
		if(str.size()!=0){
			printf("box_coord_y1 :      [%ld]\n", atol(s));
			if (yyi.size() != xxi.size()-1){
				printf("box_coord_y1 requires at least box_coord_x1, box_coord_x2, box_coord_y2. Exiting.\n");
				return -1;
			}
			yyi.push_back(atol(s));//box_coord_y1 =  ;
		}
	}

	s = iniparser_getstring(ini, (char*)"commons:box_coord_y2", NULL);
	if(s==NULL){
		printf("Warning : The line corresponding to box_coord_y2 in the ini file has been erased : commons:box_coord_y2\n");
	}else{
		str=(string)s;
		if(str.size()!=0){
			printf("box_coord_y2 :      [%ld]\n", atol(s));
			if (yyf.size() != xxi.size()-1){
				printf("box_coord_y2 requires at least box_coord_x1, box_coord_x2, box_coord_y1. Exiting. \n");
				return -1;
			}
			yyf.push_back(atol(s));//box_coord_y2 =  ;
		}
	}

	b = iniparser_getboolean(ini, "commons:map_flagged_data", -1);
	if(b!=-1){
		printf("map_flagged_data:    [%d]\n", b);
		flgdupl=b;
	}
	//flgdupl = False ;

	b = iniparser_getboolean(ini, "sanepic_preprocess:no_baseline", -1);
	if(b!=-1){
		printf("no_baseline:    [%d]\n", b);
		NORMLIN=b;
	}
	//NORMLIN = False ;

	b = iniparser_getboolean(ini, "sanepic_preprocess:correlation", -1);
	if(b!=-1){
		printf("correlation:    [%d]\n", b);
		CORRon=b;
	}//CORRon = True

	b = iniparser_getboolean(ini, "sanepic_conjugate_gradient:project_gaps", -1);
	if(b!=-1){
		printf("projgaps:    [%d]\n", b);
		projgaps=b;
	}//projgaps= False


	b = iniparser_getboolean(ini, "sanepic_preprocess:remove_poly", -1);
	if(b!=-1){
		printf("remove_poly:    [%d]\n", b);
		remove_polynomia=b;
	}//remove_poly = True

	i = iniparser_getint(ini, "sanepic_conjugate_gradient:iterW", 0);
	//if(isnan((double)i)){
	if(i>0){
		printf("iterw :      [%d]\n", i);
		iterw=i;
	}//iterw =  ;

	/*
	b = iniparser_getboolean(ini, "sanepic_conjugate_gradient:Noise_PS_estimation", -1);
	if(b!=-1){
		printf("doInitPS:    [%d]\n", b);
		doInitPS=b;
	}//doInitPS = False
	 */

	s = iniparser_getstring(ini, (char*)"sanepic_estim_PS:noise_estim", NULL);
	if(s==NULL){
		printf("You must add a line corresponding to the mixing matrix of noise components in the ini file : sanepic_conjugate_gradient:noise_estim\n");
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
	/*if (fframes_vec.size() != nsamples_vec.size()) {
		cerr << "Must give at least one first frame number. Exiting.\n";
		return -1 ;
	}*/
	if (xxi.size() != xxf.size() || xxi.size() != yyi.size() || xxi.size() != yyf.size()) {
		cerr << "box_coord_x1 box_coord_x2 box_coord_y1 box_coord_y2 must have the same size. Exiting.\n";
		return -1 ;
	}

	//ntotscan = fframes_vec.size();
	ndet = bolonames.size();

	int nnf;
	nnf = (int)extentnoiseSP.size();
	//nnf=1; // Debug
	if (nnf != 1 && nnf != ntotscan){
		cerr << "ERROR: There should be one noise power spectrum file per scan, or a single one for all the scans. Check -K options" << endl;
		return -1;
	}
	if (nnf == 1 && ntotscan > 1)
		extentnoiseSP.resize(ntotscan, extentnoiseSP[0]);

	//  printf("%d\n",nnf);

	if (((int)fcut.size()!=nnf)&&((long)fcut.size()!=ntotscan)){
		cerr << "Please give a correct number of noise cut frequency : 1 or 1 per scan\n";
		exit(0);
	}

	if(fcut.size()==1)
		fcut.resize(ntotscan, fcut[0]);



	iniparser_freedict(ini);
	return 0 ;
}
