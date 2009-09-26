/*
 * parsePre.cpp
 *
 *  Created on: 12 juin 2009
 *      Author: matthieu
 */

#include "inputFileIO.h"
#include "parsePre.h"

//int parse_sanePos_ini_file(char * ini_name)
int parse_sanePre_ini_file(char * ini_name, int  &shift_data_to_point, long  &napod,double &fsamp, bool &NOFILLGAP,bool &NORMLIN, bool &flgdupl,
		bool &CORRon, long &ntotscan, long &ndet, int &nnf,	double &f_lp, double &f_lp_Nk, string &dirfile, string &outdir, string &poutdir, string &bextension,
		string &fextension, string &pextension, string &termin, string &noiseSppreffile,
		int &coordsyst, std::vector<string> &bolonames,std::vector<long> &fframes_vec, std::vector<long> &nsamples_vec, string &fname, double &pixdeg)
/*int parse_sanePos_ini_file(char * ini_name,int &bfixc, int  &shift_data_to_point, long  &napod, bool &NOFILLGAP, bool &flgdupl,
		double * srccoord, double * coordscorner, double &radius, long &ntotscan, long &ndet, long &nnf,
		double &pixdeg, double * tancoord, double * tanpix, string &dirfile, string &outdir, string &poutdir, string &bextension,
		string &fextension, string &pextension, string &file_offsets, string &file_frame_offsets, string &termin,
		int &coordsyst, std::vector<string> &bolonames, std::vector<long> &fframes_vec, std::vector<long> &nsamples_vec,
		std::vector<long> &xxi, std::vector<long> &xxf, std::vector<long> &yyi, std::vector<long> &yyf, std::vector<double> &fcut, std::vector<string> &extentnoiseSP)
 */{
	dictionary	*	ini ;

	/* Some temporary variables to hold query results */
	int				b ;
	int				i ;
	double			d ;
	char		*	s ;
	long l;
	string str;

	// for poutdir default value
	char * pPath;
	pPath = getenv ("TMPBATCH");
	if (pPath!=NULL)
		printf ("The current path is: %s",pPath);


	std::vector<long> xxi, xxf, yyi, yyf; // box for crossing constraints removal coordinates lists (left x, right x, top y, bottom y)




	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	// printf dictionnary to stderr for debugging
	//iniparser_dump(ini, stderr);


	s = iniparser_getstring(ini, "sanepic_parallel_scheme:fname",NULL);
	str=(string)s;
	if((s!=NULL)&&(str.size()!=0)){
		printf("fname : [%s]\n",s);
		fname=s;
	}else{
		printf("You need to run Find_best_frame_order first !\n");
	}

	/* Get sanepic_preprocess attributes */
	printf("sanepic_preprocess:\n");

	s = iniparser_getstring(ini, "commons:data_directory", NULL);
	if(s!=NULL){
		printf("data_directory : [%s]\n",s);
		dirfile = s;
	}else{
		printf("You must specify a data directory : commons:data_directory\n");
		exit(0);
	}//./RCW_120_M/

	s = iniparser_getstring(ini, "commons:channel",NULL);
	if(s!=NULL){
		printf("channel file : [%s]\n",s);
		read_strings((string)s, bolonames);
	}else{
		printf("You must specify a bolometer file : commons:channel\n");
		exit(0);
	}//	channel =./RCW_120_M/bolos_commons.txt ;

	const char *temp;

	s = iniparser_getstring(ini, "commons:frame_file",NULL);
	if(s!=NULL){
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
		exit(0);
	}//frame_file =./RCW_120_M/frame_file.txt ;

	i = iniparser_getint(ini, "commons:coord_syst", -1);
	if((i==1)||(i==2)||(i==3)){
		printf("Coordinate system :      [%d]\n", i);
		coordsyst=i;
	}else{
		printf("Choose a coordinate system between 1 and 3 : commons:coord_syst\n");
		exit(0);
	}//coord_syst = 1 ;

	s = iniparser_getstring(ini, "commons:bolofield_extension",NULL);
	if(s!=NULL){
		printf("bolofield extension : [%s]\n",s);
		bextension=s;
	}else{
		printf("You must specify a bolo extension : commons:bolofield_extension\n");
		exit(0);
	}//_data

	i = iniparser_getint(ini, "commons:apodize_Nsamples", -1);
	if(i>0){
		printf("apodize_Nsamples :      [%d]\n", i);
		napod=i;
	}else{
		printf("You must choose a number of samples to apodize commons:apodize_Nsamples\n");
		exit(0);
	}//apodize_Nsamples = 100 ;

	d = iniparser_getdouble(ini,(char*)"sanepic_preprocess:sampling_frequency", -1.0);
	if (d<=0.0){
		printf("sampling_frequency cannot be negative or 0 ! Or maybe you forgot to mention sampling frequency \n");
		exit(0);
	}else{
		printf("sampling_frequency  :   [%g]\n", d);
		fsamp=d;
	}//sampling_frequency = 25.0 ;

	d = iniparser_getdouble(ini,(char*)"sanepic_preprocess:filter_frequency", -1.0);
	if (d<=0.0){
		printf("filter_frequency cannot be negative ! or maybe you have to mention filter frequency \n");
		exit(0);
	}else{
		printf("filter_frequency  :   [%g]\n", d);
		f_lp=d;
	}//filter_frequency = 0.005 ;

	d = iniparser_getdouble(ini,(char*)"sanepic_preprocess:noise_cut_frequency", 0.0);
	if (d<0.0){
		printf("noise_cut_frequency cannot be negative !\n");
		exit(0);
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
	}

	s = iniparser_getstring(ini, "commons:flag_field_extension",NULL);
	if(s!=NULL){
		printf("flag_field_extension : [%s]\n",s);
		fextension=s;
	}else{
		printf("you must specify flag_field_extension\n");
		exit(0);
	}//flag_field_extension = _flag ;

	s = iniparser_getstring(ini, "commons:pointing_field_extension",NULL);
	if(s!=NULL){
		printf("pointing_field_extension : [%s]\n",s);
		pextension=s;
	}else{
		printf("you must specify pointing_field_extension\n");
		exit(0);
	}//pointing_field_extension = _def ;

	s = iniparser_getstring(ini, "sanepic_preprocess:noise_prefixe",NULL);
	if(s!=NULL){
		printf("noise_prefixe : [%s]\n",s);
		//extentnoiseSp.push_back(s);
	}else{
		printf("you must specify noise_prefixe\n");
		exit(0);
	}//noise_prefixe = NoisePS ;

	s = iniparser_getstring(ini, "commons:noise_suffixe",NULL);
	if(s!=NULL){
		printf("noise_suffixe : [%s]\n",s);
		noiseSppreffile=s;
	}else{
		printf("you must specify a noise_suffixe\n");
		exit(0);
	}//noise_suffixe = ./RCW_120_M/ ;

	s = iniparser_getstring(ini, "commons:output_dir",NULL);
	if(s!=NULL){
		printf("output_dir : [%s]\n",s);
		outdir=s;
	}else{
		printf("Using default output directory : [%s]\n",pPath);
	}//output_dir = ./RCW_120_M/ ;

	s = iniparser_getstring(ini, "commons:out_file_str",NULL);
	if(s!=NULL){
		printf("out_file_str : [%s]\n",s);
		termin=s;
	}else{
		printf("You must specify a prefixe for generated files : out_file_str\n");
	}//out_file_str = sanepic ;

	/*s = iniparser_getstring(ini, "sanepic_preprocess:offset_file",NULL);
	if(s!=NULL){
		printf("offset_file : [%s]\n",s);
		file_offsets=s;
	}else{
		printf("You must specify bolometers position : offset_file\n");
		exit(0);
	}//offset_file = ./RCW_120_M/bolo_positions.txt ;*/

	/*s = iniparser_getstring(ini, "sanepic_preprocess:file_frame_offsets",NULL);
	str=(string)s;
	if((s!=NULL)&&(str.size()!=0)){
		printf("file_frame_offsets : [%s]\n",s);
		file_frame_offsets=s;
	}//file_frame_offsets =  ;*/

	i = iniparser_getint(ini, "commons:time_offset", 0);
	//if(isnan((double)i)){
	if(i!=0){
		printf("time_offset :      [%d]\n", i);
		shift_data_to_point=i;
	}
	/*}else{
		//printf("test was good\n");
	}//time_offset =  ;*/



	/*	i = iniparser_getint(ini, (char*)"sanepic_preprocess:box_coord_x1", -1); // faire un data file
	if((i!=-1)&&(i!=0)){
			printf("box_coord_x1_1 :      [%d]\n", i);

			if (xxi.size() != yyf.size()){
				printf("box_coord_x1 requires at least box_coord_x2, box_coord_y1, box_coord_y2. Exiting.\n");
				exit(1);
			}
			xxi.push_back(i);//box_coord_x1 =  ;
		}*/

	s = iniparser_getstring(ini, "sanepic_compute_positions:pixsize",NULL);
	if(s==NULL){
		printf("Pixel size is not referred in ini file : add a line : sanepic_compute_positions:pixsize\n");
		return -1;
	}
	str=(string)s;
	if((str.size()!=0)){
		temp=str.c_str();
		pixdeg=atof(temp);
		cout << "Pixsize : " << setprecision(str.size()) << pixdeg << endl;
	}else{
		printf("Pixsize cannot be negative ! or you forgot to mention pixel size\n");
		return -1 ;
	}


	// crossing constraint removal box coordinates
	s = iniparser_getstring(ini, (char*)"commons:box_coord_x1", NULL); // faire un data file
	str=(string)s;
	if((s!=NULL)&&(str.size()!=0)){
		printf("box_coord_x1 :      [%ld]\n", atol(s));

		if (xxi.size() != yyf.size()){
			printf("box_coord_x1 requires at least box_coord_x2, box_coord_y1, box_coord_y2. Exiting.\n");
			exit(1);
		}
		xxi.push_back(atol(s));//box_coord_x1 =  ;
	}


	s = iniparser_getstring(ini, (char*)"commons:box_coord_x2", NULL);
	str=(string)s;
	if((s!=NULL)&&(str.size()!=0)){
		printf("box_coord_x2 :      [%ld]\n", atol(s));
		if (xxf.size() != xxi.size()-1){
			printf("box_coord_x2 requires at least box_coord_x1, box_coord_y1, box_coord_y2. Exiting. \n");
			exit(1);
		}
		xxf.push_back(atol(s));//box_coord_x2 =  ;
	}

	s = iniparser_getstring(ini, (char*)"commons:box_coord_y1", NULL);
	str=(string)s;
	if((s!=NULL)&&(str.size()!=0)){
		printf("box_coord_y1 :      [%ld]\n", atol(s));
		if (yyi.size() != xxi.size()-1){
			printf("box_coord_y1 requires at least box_coord_x1, box_coord_x2, box_coord_y2. Exiting.\n");
			exit(1);
		}
		yyi.push_back(atol(s));//box_coord_y1 =  ;
	}

	s = iniparser_getstring(ini, (char*)"commons:box_coord_y2", NULL);
	str=(string)s;
	if((s!=NULL)&&(str.size()!=0)){
		printf("box_coord_y2 :      [%ld]\n", atol(s));
		if (yyf.size() != xxi.size()-1){
			printf("box_coord_y2 requires at least box_coord_x1, box_coord_x2, box_coord_y1. Exiting. \n");
			exit(1);
		}
		yyf.push_back(atol(s));//box_coord_y2 =  ;
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


	// path in which data are written
	if (pPath != NULL){
		poutdir = pPath;
	} else {
		poutdir = outdir;
	}



	// Set default parameter values
	if (fframes_vec.size() == 0) {
		fframes_vec.push_back(0);
		nsamples_vec.push_back(-1);
	}
	if (fframes_vec.size() == 1 && nsamples_vec.size() == 0)
		nsamples_vec.push_back(-1);

	// Check improper usage
	if (bolonames.size() == 0) {
		cerr << "Must provide at least one channel.\n\n";
		exit(1);
		//usage(argv[0]);
	}
	if (fframes_vec.size() != nsamples_vec.size()) {
		cerr << "Must give at least one first frame number. Exiting.\n";
		exit(1);
	}
	if (xxi.size() != xxf.size() || xxi.size() != yyi.size() || xxi.size() != yyf.size()) {
		cerr << "box_coord_x1 box_coord_x2 box_coord_y1 box_coord_y2 must have the same size. Exiting.\n";
		exit(1);
	}

	ntotscan = fframes_vec.size();
	ndet = bolonames.size();

	//nnf = extentnoiseSp.size();
	nnf=1; // Temporarily
	if (nnf != 1 && nnf != ntotscan){
		cerr << "ERROR: There should be one noise power spectrum file per scan, or a single one for all the scans. Check -K options" << endl;
		exit(1);
	}
	if (nnf == 1 && ntotscan > 1)
		//extentnoiseSp.resize(ntotscan, extentnoiseSp[0]);

		//  printf("%d\n",nnf);


		iniparser_freedict(ini);
	return 0 ;
}
