/*
 * parsePos.cpp
 *
 *  Created on: 11 juin 2009
 *      Author: matthieu
 */


#include "parsePos.h"

int parse_sanePos_ini_file(char * ini_name,bool &bfixc, int  &shift_data_to_point, long  &napod, bool &NOFILLGAP, bool &flgdupl,
		double * srccoord, double * coordscorner, double &radius, long &ntotscan, long &ndet, int &nnf,
		double &pixdeg, string &dirfile, string &outdir, string &poutdir, string &bextension,
		string &fextension, string &pextension, string &file_offsets, string &file_frame_offsets, /*string &termin,*/
		int &coordsyst, std::vector<string> &bolonames,std::vector<long> &fframes_vec, std::vector<long> &nsamples_vec,
		std::vector<long> &xxi,std::vector<long> &xxf, std::vector<long> &yyi, std::vector<long> &yyf)
{
	dictionary	*	ini ;

	/* Some temporary variables to hold query results */
	int				b ;
	int				i ;
	double			d ;
	char		*	s ;
	long l;
	string str;
	const char *temp;


	//std::vector<long> xxi, xxf, yyi, yyf; // box for crossing constraints removal coordinates lists (left x, right x, top y, bottom y)



	/* COMMAND-LINE PARAMETER PROCESSING */
	int tmpcount = 0; // parser variable (check if parsing was performed well)
	int tmpcount2 = 0; // parser variable (check if parsing was performed well)

	// load dictionnary
	ini = iniparser_load(ini_name);

	if (ini==NULL) {
		fprintf(stderr, "cannot parse file: %s\n", ini_name);
		return -1 ;
	}

	// printf dictionnary to stderr for debugging
	//iniparser_dump(ini, stderr);



	/* Get sanepic_compute_positions attributes */
	printf("sanepic_compute_positions:\n");

	//#ifdef USE_MPI
	// for poutdir default value
	char * pPath;
	pPath = getenv ("TMPBATCH");
	if (pPath!=NULL){
		outdir=pPath;
		printf ("The current path is: %s\n",pPath);
	}else{
		//#else

		s = iniparser_getstring(ini, "commons:temp_dir",NULL);
		if(s==NULL){
			printf("Warning : The line corresponding to temporary directory in the ini file has been erased : commons:output_dir\n");
			cout << "Using default output directory : " << dirfile << endl;
			outdir=dirfile;
		}else{
			str=(string)s;
			if(str.size()!=0){
				printf("temp_dir : [%s]\n",s);
				outdir=s;
			}else{
				cout << "Using default output directory : " << dirfile << endl;
				outdir=dirfile;
			}//output_dir = ./RCW_120_M/ ;
		}
	}
	//#endif


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
		read_bolofile((string)s, bolonames);
	}else{
		printf("You must specify a bolometer file : commons:channel\n");
		return -1 ;
	}//	channel =./RCW_120_M/bolos_commons.txt ;



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
		read_bolofile((string)s,dummy);
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

	/*d = iniparser_getdouble(ini,(char*)"sanepic_compute_positions:pixsize", -1.0);
	if (d<=0){
		printf("Pixsize cannot be negative ! or you forgot to mention pixel size\n");
		return -1 ;
	}else{
		printf("Pixsize :   [%g]\n", d);
		pixdeg=d;
	}//pixsize =  0.00168725828819 ;*/

	i = iniparser_getint(ini, "commons:coord_syst", -1);
	if((i==1)||(i==2)||(i==3)){
		printf("Coordinate system :      [%d]\n", i);
		coordsyst=i;
	}else{
		printf("Choose a coordinate system between 1 and 3 : commons:coord_syst\n");
		return -1 ;
	}//coord_syst = 1 ;

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

	i = iniparser_getint(ini, "commons:apodize_Nsamples", -1);
	if(i>0){
		printf("apodize_Nsamples :      [%d]\n", i);
		napod=i;
	}else{
		printf("You must choose a number of samples to apodize commons:apodize_Nsamples\n");
		return -1 ;
	}//apodize_Nsamples = 100 ;


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


	/*
	s = iniparser_getstring(ini, "sanepic_compute_positions:noise_prefixe",NULL);
	if(s!=NULL){
		printf("noise_prefixe : [%s]\n",s);
		//extentnoiseSp.push_back(s);
	}else{
		printf("you must specify noise_prefixe\n");
		exit(0);
	}//noise_prefixe = NoisePS ;*/

	/*s = iniparser_getstring(ini, "sanepic_compute_positions:noise_suffixe",NULL);
	if(s!=NULL){
		printf("noise_suffixe : [%s]\n",s);
		noiseSppreffile=s;
	}else{
		printf("you must specread_bolofileify a noise_suffixe\n");
		exit(0);
	}//noise_suffixe = ./RCW_120_M/ ;*/



	/*s = iniparser_getstring(ini, "commons:out_file_str",NULL);
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
		return -1 ;
	}//out_file_str = sanepic ;*/

	s = iniparser_getstring(ini, "sanepic_compute_positions:offset_file",NULL);
	if(s==NULL){
		printf("You must add a line corresponding to bolo positions in the ini file : sanepic_compute_positions:offset_file\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("offset_file : [%s]\n",s);
		file_offsets=s;
	}else{
		printf("You must specify bolometers position : offset_file\n");
		return -1 ;
	}//offset_file = ./RCW_120_M/bolo_positions.txt ;

	s = iniparser_getstring(ini, "sanepic_compute_positions:file_frame_offsets",NULL);
	if(s==NULL){
		printf("Warning : The line corresponding to the file containing offsets for different frame range has been erased : sanepic_compute_positions:file_frame_offsets\n");
	}else{
		str=(string)s;
		if(str.size()!=0){
			printf("file_frame_offsets : [%s]\n",s);
			file_frame_offsets=s;
		}
	}//file_frame_offsets =  ;


	b = iniparser_getboolean(ini, "commons:nofill_gap", -1);
	if(b!=-1){
		printf("nofill_gap:    [%d]\n", b);
		NOFILLGAP=b;
	}
	//NOFILLGAP = 0 ;

	i = iniparser_getint(ini, "sanepic_compute_positions:time_offset", 0);
	//if(isnan((double)i)){
	if(i!=0){
		printf("time_offset :      [%d]\n", i);
		shift_data_to_point=i;
	}
	/*}else{
		//printf("test was good\n");
	}//time_offset =  ;*/

	//TODO: faire un fichier + centre + radius au lieu de carré
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


	d = iniparser_getdouble(ini, (char*)"sanepic_compute_positions:RA_source", 0.0);
	if(d!=0.0){
		printf("RA_source :      [%g]\n", d);
		srccoord[0]=d;
		tmpcount2 += 1;
	}//RA_source =  ;

	d = iniparser_getdouble(ini, (char*)"sanepic_compute_positions:DEC_source", 0.0);
	if(d!=0.0){
		printf("DEC_source :      [%g]\n", d);
		srccoord[1]=d;
		tmpcount2 += 1;
	}//DEC_source =  ;

	d = iniparser_getdouble(ini, (char*)"sanepic_compute_positions:map_radius", 0.0);
	if(d!=0.0){
		printf("map_radius :      [%g]\n", d);
		radius=d;
		tmpcount2 += 1;
	}//map_radius =  ;

	d = iniparser_getdouble(ini, (char*)"sanepic_compute_positions:RA_min", 0.0);
	if(d!=0.0){
		printf("RA_min :      [%g]\n", d);
		coordscorner[0] = d;
		tmpcount += 1;
	}//RA_min =  ;

	d = iniparser_getdouble(ini, (char*)"sanepic_compute_positions:RA_max", 0.0);
	if(d!=0.0){
		printf("RA_max :      [%g]\n", d);
		coordscorner[1] = d;
		tmpcount += 1;
	}//RA_max =  ;

	d = iniparser_getdouble(ini, (char*)"sanepic_compute_positions:DEC_min", 0.0);
	if(d!=0.0){
		printf("DEC_min :      [%g]\n", d);
		coordscorner[2] = d;
		tmpcount += 1;
	}
	//DEC_min =  ;

	d = iniparser_getdouble(ini, (char*)"sanepic_compute_positions:DEC_max", 0.0);
	if(d!=0.0){
		printf("DEC_max :      [%g]\n", d);
		coordscorner[3] = d;
		tmpcount += 1;
	}//DEC_max =  ;


	b = iniparser_getboolean(ini, "commons:map_flagged_data", -1);
	if(b!=-1){
		printf("map_flagged_data:    [%d]\n", b);
		flgdupl=b;
	}//flgdupl = False ;


	if (tmpcount == 1 || tmpcount == 2 || tmpcount == 3){
		cerr << "ERROR: None or all the following keywords must be set: RA_min RA_max DEC_min DEC_max. Exiting. \n";
		return -1 ;
	}
	if (tmpcount2 == 1 || tmpcount2 == 2){
		cerr << "ERROR: None or all the following keywords must be set: RA_source DEC_source map_radius. Exiting. \n";
		return -1 ;
	}



	if (tmpcount == 4)
		bfixc = 1;

	if (tmpcount2 == 3){
		if (tmpcount == 4){
			cerr << "ERROR: Conflicting input parameter: RA_min RA_max DEC_min DEC_max keywords are not compatible with RA_source DEc_source map_radius keywords . Exiting. \n";
			return -1 ;
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
		cerr << "ERROR: You must provide coordinates of the source in RA/DEC for telescope coordinates, use RA_source DEC_source map_radius\n";
		return -1 ;
	}



	// Set default parameter values
	if (fframes_vec.size() == 0) {
		fframes_vec.push_back(0);
		nsamples_vec.push_back(-1);
	}
	if (fframes_vec.size() == 1 && nsamples_vec.size() == 0)
		nsamples_vec.push_back(-1);

	// Check improper usage
	//if (dirfile == "") usage(argv[0]);
	if (bolonames.size() == 0) {
		cerr << "Must provide at least one channel.\n\n";
		return -1 ;
		//usage(argv[0]);
	}
	if (fframes_vec.size() != nsamples_vec.size()) {
		cerr << "Must give at least one first frame number. Exiting.\n";
		return -1 ;
	}
	if (xxi.size() != xxf.size() || xxi.size() != yyi.size() || xxi.size() != yyf.size()) {
		cerr << "box_coord_x1 box_coord_x2 box_coord_y1 box_coord_y2 must have the same size. Exiting.\n";
		return -1 ;
	}

	ntotscan = fframes_vec.size();
	ndet = bolonames.size();

	printf("Number of scans : %ld\n",ntotscan);
	printf("%ld bolometers will be used\n",ndet);


	iniparser_freedict(ini);
	return 0 ;
}

