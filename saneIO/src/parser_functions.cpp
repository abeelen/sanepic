/*
 * parser_functions.cpp
 *
 *  Created on: 20 oct. 2009
 *      Author: matthieu
 */


#include <iostream>
#include <iomanip>
//#include <fstream>
//#include <cstdlib>
//#include <cstdio>
//#include <string>
//#include <unistd.h>
//#include <string>
//#include <vector>
//#include <algorithm>
//
//#include "dataIO.h"
//#include "mpi_architecture_builder.h"

extern "C"{
#include "iniparser.h"
#include "dictionary.h"
}

#include "inputFileIO.h"
#include "parser_functions.h"

using namespace std;

int read_dirfile(dictionary	*ini, struct directories &dir, int rank){

	char *s;
	string str;

	s = iniparser_getstring(ini, "commons:data_directory", NULL);
	if(s==NULL){
		if(rank==0)
			printf("You must add a line in ini file specifying a data directory : commons:data_directory\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		//printf("data_directory : [%s]\n",s);
		dir.dirfile = s;
	}else{
		if(rank==0)
			printf("You must specify a data directory : commons:data_directory\n");
		return -1 ;
	}//./RCW_120_M/

	return 0;
}


int read_tmpdir(dictionary	*ini, struct directories &dir, int rank){



	char *s, *pPath;
	string str;


	pPath = getenv ("TMPBATCH");
	if (pPath!=NULL){
		dir.tmp_dir=pPath;
		if(rank==0)
			printf ("The current path is: %s\n",pPath);
	}else{
		//#else

		s = iniparser_getstring(ini, "commons:temp_dir",NULL);
		if(s==NULL){
			if(rank==0){
				printf("Warning : The line corresponding to temporary directory in the ini file has been erased : commons:output_dir\n");
				cout << "Using default output directory : " << dir.dirfile << endl;
			}
			dir.tmp_dir=dir.dirfile;
		}else{
			str=(string)s;
			if(str.size()!=0){
				//printf("temp_dir : [%s]\n",s);
				dir.tmp_dir=s;
			}else{
				if(rank==0)
					cout << "Warning : The line corresponding to temporary directory in the ini file is clear : Using default output directory : " << dir.dirfile << endl;
				dir.tmp_dir=dir.dirfile;
			}//output_dir = ./RCW_120_M/ ;
		}
	}


	return 0;
}


int read_outdir(dictionary	*ini, struct directories &dir, int rank){

	char *s;
	string str;


	s = iniparser_getstring(ini, "commons:output_dir",NULL);
	if(s==NULL){
		if(rank==0){
			printf("Warning : The line corresponding to output directory in the ini file has been erased : commons:output_dir\n");
			cout << "Using default output directory : " << dir.dirfile << endl;
		}
		dir.outdir=dir.dirfile;
	}else{
		str=(string)s;
		if(str.size()!=0){
			//printf("output_dir : [%s]\n",s);
			dir.outdir=s;
		}else{
			if(rank==0)
				cout << "Using default output directory : " << dir.dirfile << endl;
			dir.outdir=dir.dirfile;
		}//output_dir = ./RCW_120_M/ ;
	}

	return 0;

}


int read_channel_list(dictionary	*ini, std::vector<string> &bolonames, int rank){

	string str;
	char *s;

	s = iniparser_getstring(ini, "commons:channel",NULL);
	if(s==NULL){
		if(rank==0)
			printf("You must add a line in ini file corresponding to a bolometer file : commons:channel\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		if(rank==0)
			printf("channel file : [%s]\n",s);
		read_strings((string)s, bolonames);
	}else{
		if(rank==0)
			printf("You must specify a bolometer file : commons:channel\n");
		return -1 ;
	}//	channel =./RCW_120_M/bolos_commons.txt ;

	return 0;

}


int read_fits_file_list(dictionary	*ini, struct directories &dir, struct samples &samples_str, int rank){

	char *s;
	string str;

	s = iniparser_getstring(ini, "commons:fits_filelist",NULL);
	if(s==NULL){
		if(rank==0)
			printf("You must add a line in the ini file corresponding to a frame file : commons:fits_filelist\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		samples_str.filename=str;
		read_fits_list(str, samples_str.fitsvect, samples_str.noisevect, samples_str.scans_index, samples_str.framegiven);
		cout << "fitsvect " << samples_str.fitsvect[0] << endl; //" " << samples_str.fitsvect[1] << " " << samples_str.fitsvect[2] << " " << samples_str.fitsvect[3] << endl;
		cout << "noisevect " << samples_str.noisevect[0] << endl; //" " << samples_str.noisevect[1] << " " << samples_str.noisevect[2] << " " << samples_str.noisevect[3] << endl;
		cout << "scans_index " << samples_str.scans_index[0] << endl; //" " << samples_str.scans_index[1] << " " << samples_str.scans_index[2] << " " << samples_str.scans_index[3] << endl;

		for(int ii=0;ii<(int)((samples_str.fitsvect).size());ii++){
			cout << dir.dirfile + samples_str.fitsvect[ii] << endl;
			samples_str.fitsvect[ii] = dir.dirfile + samples_str.fitsvect[ii];
		}

		readFrames(samples_str.fitsvect, samples_str.nsamples);

		cout << "after read frames\n";

		if((int)(samples_str.noisevect).size()==0){ // if no noise file is given in fits_filelist
			// read the noise file name in the ini file
			s = iniparser_getstring(ini, "commons:noise_prefixe",NULL);
			if(s==NULL){
				if(rank==0)
					printf("You must add a line in the ini file corresponding to a frame file : commons:noise_prefixe\n");
				return -1;
			}
			str=(string)s;
			if(str.size()!=0){
				samples_str.noisevect.push_back(str);
				// meme fichier de bruit pour tous les scans
				(samples_str.noisevect).resize(samples_str.fitsvect.size(),samples_str.noisevect[0]);

			}else{
				if(rank==0)
					printf("You must specify a fits noise filename : commons:noise_prefixe\n");
				return -1 ;
			}//frame_file =./RCW_120_M/fits_files.txt ;
		}

	}else{
		if(rank==0)
			printf("You must specify a fits filelist : commons:fits_filelist\n");
		return -1 ;
	}//frame_file =./RCW_120_M/fits_files.txt ;


	return 0;

}


int read_pixel_size(dictionary	*ini, struct input_commons &com, int rank){

	char *s;
	string str;
	const char* temp;

	s = iniparser_getstring(ini, "sanepic_compute_positions:pixsize",NULL);
	if(s==NULL){
		if(rank==0)
			printf("You must add a line in the ini file corresponding to pixel size : sanepic_compute_positions:pixsize\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		temp=str.c_str();
		com.pixdeg=atof(temp);
		//cout << "Pixsize : " << setprecision(str.size()) << u_opt.pixdeg << endl;
	}else{
		if(rank==0)
			printf("Pixsize cannot be negative ! or you forgot to mention pixel size\n");
		return -1 ;
	}

	return 0;

}

int read_apodize_samples(dictionary	*ini, struct input_commons &com, int rank){

	int i;

	i = iniparser_getint(ini, "commons:apodize_Nsamples", -1);
	if(i>0){
		//printf("apodize_Nsamples :      [%d]\n", i);
		com.napod=i;
	}else{
		if(rank==0)
			printf("You must choose a number of samples to apodize commons:apodize_Nsamples\n");
		return -1 ;
	}//apodize_Nsamples = 100 ;

	return 0;

}

int read_nofillgap(dictionary	*ini, struct input_commons &com, int rank){

	bool b;

	b = iniparser_getboolean(ini, "commons:nofill_gap", 0);
	//if(b!=0){
	//printf("nofill_gap:    [%d]\n", b);
	com.NOFILLGAP=b;
	//}
	//NOFILLGAP = 0 ;

	return 0;

}


int read_box_coord(dictionary	*ini, std::vector<struct box> &boxFile, int rank){

	char *s;
	string str;

	s = iniparser_getstring(ini, (char*)"commons:box_coord_file", NULL);
	if(s==NULL){
		if(rank==0)
			printf("Warning : The line corresponding to box_coord_file in the ini file has been erased : commons:box_coord_file\n");
	}else{
		str=(string)s;
		if(str.size()!=0){
			//printf("box_coord_file :      [%ld]\n", atol(s));

			//std::vector<box> boxList;
			readBoxFile(str, boxFile);
		}
	}

	return 0;
}

/*
int read_RA_DEC_min_max(dictionary	*ini, struct user_options_sanepos &u_opt, int &tmpcount){

	double d;

	d = iniparser_getdouble(ini, (char*)"sanepic_compute_positions:RA_min", 0.0);
	if(d!=0.0){
		//printf("RA_min :      [%g]\n", d);
		u_opt.coordscorner[0] = d;
		tmpcount += 1;
	}//RA_min =  ;

	d = iniparser_getdouble(ini, (char*)"sanepic_compute_positions:RA_max", 0.0);
	if(d!=0.0){
		//printf("RA_max :      [%g]\n", d);
		u_opt.coordscorner[1] = d;
		tmpcount += 1;
	}//RA_max =  ;

	d = iniparser_getdouble(ini, (char*)"sanepic_compute_positions:DEC_min", 0.0);
	if(d!=0.0){
		//printf("DEC_min :      [%g]\n", d);
		u_opt.coordscorner[2] = d;
		tmpcount += 1;
	}
	//DEC_min =  ;

	d = iniparser_getdouble(ini, (char*)"sanepic_compute_positions:DEC_max", 0.0);
	if(d!=0.0){
		//printf("DEC_max :      [%g]\n", d);
		u_opt.coordscorner[3] = d;
		tmpcount += 1;
	}//DEC_max =  ;

	return 0;
}


int read_RA_DEC_radius_source(dictionary	*ini, struct user_options_sanepos &u_opt, int tmpcount){

	int tmpcount2=0;
	double d;

	d = iniparser_getdouble(ini, (char*)"sanepic_compute_positions:RA_source", 0.0);
	if(d!=0.0){
		//printf("RA_source :      [%g]\n", d);
		u_opt.srccoord[0]=d;
		tmpcount2 += 1;
	}//RA_source =  ;

	d = iniparser_getdouble(ini, (char*)"sanepic_compute_positions:DEC_source", 0.0);
	if(d!=0.0){
		//printf("DEC_source :      [%g]\n", d);
		u_opt.srccoord[1]=d;
		tmpcount2 += 1;
	}//DEC_source =  ;

	d = iniparser_getdouble(ini, (char*)"sanepic_compute_positions:map_radius", 0.0);
	if(d!=0.0){
		//printf("map_radius :      [%g]\n", d);
		u_opt.radius=d;
		tmpcount2 += 1;
	}//map_radius =  ;


	if (tmpcount2 == 1 || tmpcount2 == 2){
		cerr << "ERROR: None or all the following keywords must be set: RA_source DEC_source map_radius. Exiting. \n";
		return -1 ;
	}


	if (tmpcount2 == 3){
		if (tmpcount == 4){
			cerr << "ERROR: Conflicting input parameter: RA_min RA_max DEC_min DEC_max keywords are not compatible with RA_source DEc_source map_radius keywords . Exiting. \n";
			return -1 ;
		}
		//u_opt.bfixc = 1;
		u_opt.coordscorner[0] = u_opt.srccoord[0];
		u_opt.coordscorner[1] = u_opt.srccoord[0];
		u_opt.coordscorner[2] = u_opt.srccoord[1];
		u_opt.coordscorner[3] = u_opt.srccoord[1];
	}

	return 0;

}
 */

int read_map_flagged_data(dictionary	*ini, struct input_commons &com, int rank){

	bool b;

	b = iniparser_getboolean(ini, "commons:map_flagged_data", 0);
	//if(b!=0){
	//printf("map_flagged_data:    [%d]\n", b);
	com.flgdupl=b;
	//}//flgdupl = False ;

	return 0;
}


int read_sampling_frequency(dictionary	*ini, struct user_options &u_opt, int rank){

	double d;

	d = iniparser_getdouble(ini,(char*)"sanepic_preprocess:sampling_frequency", -1.0);
	if (d<=0.0){
		if(rank==0)
			printf("sampling_frequency cannot be negative or 0 ! Or maybe you forgot to mention sampling frequency \n");
		return -1;
	}else{
		//printf("sampling_frequency  :   [%g]\n", d);
		u_opt.fsamp=d;
	}//sampling_frequency = 25.0 ;

	return 0;
}


int read_filter_frequency(dictionary	*ini, struct user_options &u_opt, int rank){

	double d;

	d = iniparser_getdouble(ini,(char*)"sanepic_preprocess:filter_frequency", -1.0);
	if (d<=0.0){
		if(rank==0)
			printf("filter_frequency cannot be negative ! or maybe you have to mention filter frequency \n");
		return -1;
	}else{
		//printf("filter_frequency  :   [%g]\n", d);
		u_opt.f_lp=d;
	}//filter_frequency = 0.005 ;

	return 0;
}


int read_noise_cut_freq(dictionary	*ini, std::vector<double> &fcut, int rank){

	string str;
	char *s;

	s = iniparser_getstring(ini, "sanepic_preprocess:fcut_file",NULL);
	if(s==NULL){
		if(rank==0)
			printf("You must add a line corresponding to noise cut frequency file in the parser file : sanepic_preprocess:fcut_file\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		//printf("fcut file : [%s]\n",s);
		//read frame file function
		std::vector<string> dummy2;
		read_strings(str,dummy2);
		if(((int)dummy2.size())==0){
			if(rank==0)
				printf("You must provide at least one number of noise cut frequency (or one per scan) in fcut_file !\n");
			return -1;}
		for(int ii=0; ii<(int)dummy2.size(); ii++)
			fcut.push_back(atof(dummy2[ii].c_str()));

	}else{
		if(rank==0)
			printf("You must specify a noise cut frequency file : sanepic_preprocess:fcut_file\n");
		return -1;
	}//frame_file =./RCW_120_M/fcut_file.txt ;

	return 0;

}

/*
int read_noise_file_list(dictionary	*ini, std::vector<string> &extentnoiseSP){


	string str;
	char *s;

	s = iniparser_getstring(ini, "commons:noise_prefixe_file",NULL);
	if(s==NULL){
		printf("You must add a line corresponding to noise_prefixe file in the ini file : commons:noise_prefixe_file\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		//printf("noise_prefixe_file : [%s]\n",s);
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

	return 0;

}*/

int read_baseline(dictionary	*ini, struct user_options &u_opt, int rank){

	bool b;

	b = iniparser_getboolean(ini, "sanepic_preprocess:no_baseline", 0);
	//if(b!=0){
	//printf("no_baseline:    [%d]\n", b);
	u_opt.NORMLIN=b;
	//}
	//NORMLIN = False ;

	return 0;

}

int read_correlation(dictionary	*ini, struct user_options &u_opt, int rank){

	bool b;

	b = iniparser_getboolean(ini, "sanepic_preprocess:correlation", 1);
	//if(b!=1){
	//printf("correlation:    [%d]\n", b);
	u_opt.CORRon=b;
	//}//CORRon = True

	return 0;
}

int read_remove_poly(dictionary	*ini, struct user_options &u_opt, int rank){

	bool b;

	b = iniparser_getboolean(ini, "sanepic_preprocess:remove_poly", 1);
	//if(b!=1){
	//printf("remove_poly:    [%d]\n", b);
	u_opt.remove_polynomia=b;
	//}//remove_poly = True

	return 0;

}

int read_projgaps(dictionary	*ini, struct user_options &u_opt, int rank){

	bool b;

	b = iniparser_getboolean(ini, "sanepic_conjugate_gradient:project_gaps", 0);
	//if(b!=0){
	//printf("projgaps:    [%d]\n", b);
	u_opt.projgaps=b;
	//}//projgaps= False

	return 0 ;
}

int read_iter(dictionary	*ini, int &iterw, int rank){

	int i;

	i = iniparser_getint(ini, "sanepic_conjugate_gradient:iterW", 0);
	//if(isnan((double)i)){
	if(i>0){
		//printf("iterw :      [%d]\n", i);
		iterw=i;
	}//iterw =  ;

	return 0;
}

int read_ell_file(dictionary	*ini, string &ellFile, int rank){


	string str;
	char *s;

	s = iniparser_getstring(ini, "sanepic_estim_PS:ell_file",NULL);
	if(s==NULL){
		if(rank==0)
			printf("You must add a line in ini file corresponding to a ell file :sanepic_estim_PS:ell_file\n");
		return -1;
	}
	ellFile=(string)s;
	if(ellFile.size()!=0){
		//printf("ell File: [%s]\n",s);
	}else{
		if(rank==0)
			printf("You must specify a ell file : sanepic_estim_PS:ell_file\n");
		return -1 ;
	}//	channel =./RCW_120_M/bolos_commons.txt ;

	return 0;
}


int read_map_file(dictionary	*ini, string &signame, int rank){


	string str;
	char *s;

	s = iniparser_getstring(ini, (char*)"sanepic_estim_PS:map_file", NULL);
	if(s==NULL){
		if(rank==0)
			printf("You must add a line corresponding to the fits map file (created by sanePic) in the ini file : sanepic_estim_PS:map_file\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		if(rank==0)
			printf("noise_estim :      [%s]\n", s);
		signame=s;

	}else{
		signame="NOSIGFILE";
	}//

	return 0;
}

int read_cov_matrix_file(dictionary	*ini, string &fname, int rank){

	string str;
	char*s;

	s = iniparser_getstring(ini, "sanepic_inv_matrix:cov_matrix_file",NULL);
	if(s!=NULL){
		if(rank==0)
			printf("cov_matrix_file: [%s]\n",s);
		fname=s;
		//read_strings((string)s, bolonames);
	}else{
		if(rank==0)
			printf("You must specify a noise covariance matrix file to invert : sanepic_inv_matrix:cov_matrix_file\n");
		return(-1);
	}//	fname = ./RCW_120_M/BoloPS0sanepic_binary.psd

	return 0;

}

int read_mixmatfile(dictionary	*ini, string &MixMatfile, int rank){

	string str;
	char *s;

	s = iniparser_getstring(ini, (char*)"sanepic_estim_PS:noise_estim", NULL);
	if(s==NULL){
		if(rank==0)
			printf("You must add a line corresponding to the mixing matrix of noise components in the ini file : sanepic_estim_PS:noise_estim\n");
		return -1;
	}
	str=(string)s;
	if(str.size()!=0){
		printf("noise_estim :      [%s]\n", s);
		MixMatfile=s;

	}else{
		if(rank==0)
			printf("You must give filename containing the mixing matrix of noise components : noise_estim\n");
		return(-1);
	}// MixMatfile = Mixlaboca

	return 0;
}

int read_directories(dictionary	*ini, struct directories &dir, int rank){

	if(read_dirfile(ini, dir, rank)==-1)
		return -1;

	if(read_tmpdir(ini, dir, rank)==-1)
		return -1;

	if(read_outdir(ini, dir, rank)==-1)
		return -1;

	return 0;
}

int read_commons(dictionary	*ini, struct input_commons &commons, int rank){

	if(read_map_flagged_data(ini,  commons, rank))
		return -1;

	if(read_pixel_size(ini,  commons, rank))
		return -1;

	if(read_apodize_samples(ini, commons, rank))
		return -1;

	if(read_nofillgap(ini, commons, rank))
		return -1;

	return 0;
}

int read_user_options(dictionary *ini,struct user_options &u_opt, int rank){


	if(read_sampling_frequency(ini, u_opt, rank)==-1)
		return -1;

	if(read_filter_frequency(ini, u_opt, rank)==-1)
		return -1;

	if(read_baseline(ini, u_opt, rank)==-1)
		return -1;

	if(read_correlation(ini, u_opt, rank)==-1)
		return -1;

	if(read_remove_poly(ini, u_opt, rank)==-1)
		return -1;

	if(read_projgaps(ini, u_opt, rank)==-1)
		return -1;



	return 0;
}

void print_commons(struct input_commons commons){

	/*if((commons.shift_data_to_point)>0)
		cout << "A time offset will be subsrtacted to the data to match the pointing : [" << commons.shift_data_to_point << "]\n";
	 */
	if(commons.napod>0)
		cout << "You have specified a number of samples to apodize : [" << commons.napod << "]\n";

	if(commons.NOFILLGAP)
		cout << "You have set nofill_gap to True : the gaps in data timeline will NOT be filled\n";
	else
		cout << "You have set nofill_gap to False (default) : the gaps in data timeline will be filled\n";

	if(commons.flgdupl)
		cout << "Flagged data are put in a separate map : map_flagged_data = True\n";

	cout << "You have specified a pixel size : [" << setprecision(14) << commons.pixdeg << " deg]\n";

}

void print_directories(struct directories dir){

	cout << "Data directory : [" << dir.dirfile << "]\n";
	cout << "Temporary directory : [" << dir.tmp_dir << "]\n";
	cout << "Output directory : [" << dir.outdir << "]\n";

}

/*
void print_parser_sanepos(struct user_options_sanepos u_opt){


	printf("sanePos parser operations completed :\n\n");
	cout << "You have specified the following options : \n";

	///////
	 print_directories(u_opt.dir);
	//cout << "Data directory : [" << u_opt.dirfile << "]\n";
	//cout << "Temporary directory : [" << u_opt.tmp_dir << "]\n";
	//cout << "Output directory : [" << u_opt.outdir << "]\n";

	///////
	 print_commons(u_opt.commons);


	//////////

	//if(u_opt.srccoord[0]>-1000){
	//	cout << "RA of the tangent point and of the source for telescope coordinate maps : [" << u_opt.srccoord[0] <<"]\n";
	//	cout << "DEC of the tangent point and of the source for telescope coordinate maps : [" << u_opt.srccoord[1] <<"]\n";
	//	cout << "fixed radius (half a side) of the map in degrees : [" << u_opt.radius << "]\n";
	//}

}*/


void print_parser(struct user_options u_opt){


	cout << "You have chosen a sampling frequency equal to : [" << u_opt.fsamp << " Hz]\n";

	if(u_opt.NORMLIN)
		cout << "No baseline will be removed from the data\n";
	else
		cout << "A baseline will be removed from the data (default)\n";

	if(u_opt.remove_polynomia)
		cout << "Remove a polynomia fitted to the data to reduce fluctuations on timescales larger than the length of the considered segment\n";
	else
		cout << "No polynomia will be used\n";

	if(u_opt.CORRon)
		cout << "Correlations between detectors are included in the analysis\n";
	else
		cout << "Correlations between detectors are NOT included in the analysis\n";

	cout << "Frequency of the high pass filter applied to the data : [" << u_opt.f_lp << " Hz]\n";

	//cout << "Noise power spectrum file prefixe : [" << u_opt.noiseSppreffile << "]\n";

	if(u_opt.projgaps)
		cout << "Gaps are projected to a pixel in the map, gap filling of noise only is performed iteratively\n";
	else
		cout <<  "Gaps are NOT projected to a pixel in the map (default)\n";


}
