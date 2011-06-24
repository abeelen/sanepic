#include <iostream>
#include <stdio.h>
#include <string.h>



#include "imageIO.h"
#include "struct_definition.h"
#include "inputFileIO.h"
#include "utilities.h"

#include <sstream>
#include <string>
#include <vector>

extern "C" {
#include <fitsio.h>
#include <nrutil.h>
#include "wcslib/wcslib.h"
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
}

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "Unknown"
#endif

using namespace std;


int get_fits_META(string fname, std::vector<string> &key, std::vector<int> &datatype, std::vector<string> &val, std::vector<string> &com){

	fitsfile *fp;
	int fits_status = 0; // MUST BE initialized... otherwise it fails on the call to the function...
	char *keylist[11]={(char *)"EQUINOX", (char *)"TIMESYS",(char *)"TYPE",(char *)"CREATOR", (char *)"INSTRUME",
			(char *)"DATE-OBS",(char *)"DATE-END", (char *)"OBJECT",
			(char *)"RADESYS", (char *)"TELESCOP", (char *)"OBSERVER"};
	int keynum=11;

	// default commentaries values
	com.push_back("Equinox of celestial coordinate system");
	com.push_back("All dates are in UTC time");
	com.push_back("Product Type Identification");
	com.push_back("Generator of this product");
	com.push_back("Instrument attached to this product");
	com.push_back("Start date of this product");
	com.push_back("End date of this product");
	com.push_back("Target name");
	com.push_back("Coordinate reference frame for the RA and DEC");
	com.push_back("Name of telescope");
	com.push_back("Observer name");


	for(long kk=0; kk< keynum-1; kk++){
		key.push_back(keylist[kk]);
		datatype.push_back(TSTRING);

		// generating default values
		val.push_back("Unknown");
	}

	datatype[0] = TINT; // TINT for equinox
	val[0]="2000"; // equinox set to 2000 by default

	if (fits_open_file(&fp, fname.c_str(), READWRITE, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	char comment[80];
	char value[80];

	int erased=0;
	for(long ii=0;ii<keynum;ii++){
		fits_read_keyword(fp, keylist[ii], value, comment, &fits_status);

		if(fits_status==0){
			val[ii]=value;
			com[ii]=comment;
		}else{
			// reset status
			fits_status=0; // keyword not found, keeping default value

			//			key.erase(key.begin()+ii-erased);
			//			datatype.erase(datatype.begin()+ii-erased);
			erased++;
		}
	}

	if(erased==12)
		return 0;

	// Change string in due form
	for(long kk=0; kk<(long)key.size(); kk++){
		string tmp =  val[kk];
		string tmp2;
		if((int)tmp[0]==39){
			tmp.erase(tmp.begin());
			int marker=0;
			for(int ii=0; ii<(int)tmp.size(); ii++)
				if((int)tmp[ii]==32 || (int)tmp[ii]==39){ // remove ' character and white spaces
					marker = ii;
					break;
				}

			if(marker>0)
				tmp2.insert(tmp2.begin(), tmp.begin(), tmp.begin()+marker);
			else
				tmp2=tmp;

			val[kk] = tmp2;
		}
	}


	// close file
	if(fits_close_file(fp, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	return 0;

}





int write_fits_wcs(string fname, struct wcsprm * wcs, long NAXIS1, long NAXIS2,  char dtype, void *data, string table_name ,bool fits_already_exist, std::vector<string> key, std::vector<int> datatype, std::vector<string> val, std::vector<string> com)
{
	// all angles in degrees

	fitsfile *fp;
	int fits_status = 0; // MUST BE initialized... otherwise it fails on the call to the function...
	int start=1;
	int naxis = 2;                  // number of dimensions
	long naxes[] = {NAXIS1, NAXIS2}; // size of dimensions
	long fpixel[] = {1, 1};          // index for write_pix
	long long ndata = NAXIS1 * NAXIS2;            // number of data points

	char *header, *hptr;
	int nkeyrec;

	if(fits_already_exist){
		if (fits_open_file(&fp, fname.c_str(), READWRITE, &fits_status)){
			fits_report_error(stderr, fits_status);
			return 1;
		}
	}else{
		// create fits file
		if ( fits_create_file(&fp, fname.c_str(), &fits_status) ){
			fits_report_error(stderr, fits_status);
			return 1;
		}
		start=0;
	}

	for(int ii=start;ii<2;ii++){
		naxis=2*ii;
		naxes[0]=NAXIS1*ii;
		naxes[1]=NAXIS2*ii;

		// create fits image (switch on data type)
		switch (dtype) {
		case 'd':    // double
			if ( fits_create_img(fp, DOUBLE_IMG, naxis, naxes, &fits_status) ){
				fits_report_error(stderr, fits_status);
				return 1;
			}
			break;
		case 'l':    // long
			if ( fits_create_img(fp, LONG_IMG, naxis, naxes, &fits_status) ){
				fits_report_error(stderr, fits_status);
				return 1;
			}
			break;
		default:
			printf("write_fits: data type %c not supported. Exiting.\n",dtype);
			return 1;
		}

		if((int)val.size()>0){
			// add META DATA
			// add equinox as an integer if equinox was found in META DATA
			char* value_char = (char*)(val[0].c_str());
			char * pEnd;
			int value_int = strtol(value_char,&pEnd,10);
			string comment_equi = com[0];

			if(fits_write_key(fp, datatype[0], (char*)key[0].c_str(), &value_int, (char*)comment_equi.c_str(), &fits_status)){
				fits_report_error(stderr, fits_status);
				fits_status=0;
			}


			// add the other META DATA as STRING
			for(long kk=1; kk<(long)key.size(); kk++){
				string value = val[kk];
				string comment = com[kk];

				if(fits_write_key(fp, datatype[kk], (char*)key[kk].c_str(), (char*)value.c_str(), (char*)comment.c_str(), &fits_status)){
					fits_report_error(stderr, fits_status);
					fits_status=0;
				}
			}
		}


		// Transform wcsprm struture to header
		if ( (fits_status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header)) ){
			printf("wcshdo ERROR %d: %s.\n", fits_status, wcs_errmsg[fits_status]);
			return 1;
		}

		hptr = header;
		// write it to the fits file
		for (int keyrec = 0; keyrec < nkeyrec; keyrec++, hptr += 80){
			if ( fits_write_record(fp, (const char*) hptr, &fits_status)){
				fits_report_error(stderr, fits_status);
				return 1;
			}
		}
		free(header);

		// write date to file
		if ( fits_write_date(fp, &fits_status) ){
			fits_report_error(stderr, fits_status);
			return 1;
		}

	}


	if(fits_update_key(fp, TSTRING, (char *)"EXTNAME", (void*)(table_name.c_str()),
			(char *) "table name", &fits_status))
		return 1;

	// write map data
	switch (dtype) {
	case 'd':    // double
		if ( fits_write_pix(fp, TDOUBLE, fpixel, ndata, (double*) data, &fits_status) ){
			fits_report_error(stderr, fits_status);
			return 1;
		}
		break;
	case 'l':    // long
		if ( fits_write_pix(fp, TLONG, fpixel, ndata, (long*) data, &fits_status) ){
			fits_report_error(stderr, fits_status);
			return 1;
		}
		break;
	}


	// Add Comments to the header
	std::vector<std::string> comments;
	comments.push_back(" ");
	comments.push_back("This fits file was generated by SANEPIC version " + (string)PACKAGE_VERSION);
	comments.push_back("For more informations about SANEPIC and for SANEPIC updates, please");
	comments.push_back("check our dedicated website : http://www.ias.u-psud.fr/sanepic ");
	comments.push_back(" ");

	for (int ii=0; ii< comments.size(); ii++){
		if(fits_write_comment(fp, (char *)comments[ii].c_str(),  &fits_status))
			return 1;
	}

	if (fits_write_chksum(fp, &fits_status)){
		cout << "error checksum !\n";
		return 1;
	}

	// close file
	if(fits_close_file(fp, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	return 0;
}

int write_fits_history2(std::string fname,long NAXIS1, long NAXIS2, struct param_common dir, struct param_sanePre proc_param, struct param_sanePos pos_param, std::vector<double> fcut, struct samples samples_struct, struct param_sanePS PS_param, struct param_sanePic Pic_param, struct param_saneInv Inv_param)
{

	fitsfile *fptr;
	int fits_status = 0;
	//	long naxes[2] = { 1, 1 };
	std::vector<string> key_vect;
	std::vector<string> comm_vect;
	std::vector<string> value_vect;

	//fill the key_vect
	// [commons]
	key_vect.push_back("DATADIR");
	key_vect.push_back("INPUTDIR");
	for(long iframe=0;iframe < samples_struct.ntotscan; iframe ++)
		key_vect.push_back("SOURCES");
	key_vect.push_back("BOLO_GLOBAL_FILE");
	key_vect.push_back("BOLO_SUFFIX");
	key_vect.push_back("TMPDIR");
	key_vect.push_back("OUTPUTDIR");

	// [sanePos]
	key_vect.push_back("PIXSIZE");
	key_vect.push_back("MAP_DUPLICATION");
	key_vect.push_back("FILE_FORMAT");
	key_vect.push_back("GAPS_PROJECTED");
	key_vect.push_back("MASK_FILE");
	key_vect.push_back("RA_NOM");
	key_vect.push_back("DEC_NOM");
	key_vect.push_back("PROJ_TYPE");

	// [sanePre]
	key_vect.push_back("NAPOD");
	key_vect.push_back("SAMPLING_FREQ");
	key_vect.push_back("FILTER_FREQ");
	key_vect.push_back("FCUT_FILE");
	key_vect.push_back("POLY_ORDER");
	key_vect.push_back("POLY_SUBTRACTION");
	key_vect.push_back("NORMLIN");
	key_vect.push_back("CORRELATION");
	key_vect.push_back("FILLGAPS");

	// [saneInv]
	key_vect.push_back("COV_MATRIX_FILE");
	key_vect.push_back("NOISEDIR");
	key_vect.push_back("COV_MATRIX_SUFFIX");

	// [sanePS]
	key_vect.push_back("ELL_GLOBAL_FILE");
	key_vect.push_back("ELL_SUFFIX");
	key_vect.push_back("MIXING_MATRIX_GLOBAL_FILE");
	key_vect.push_back("MIXING_MATRIX_SUFFIX");
	key_vect.push_back("NOISE COMPONENT");
	key_vect.push_back("MAP_FILE");
	key_vect.push_back("SAVE_DATA");

	// [sanePic]
	key_vect.push_back("ITERW");
	key_vect.push_back("ITERMAX");


	/*--------------------------------- fill comment vetor -------------------------------------*/

	// [commons]
	comm_vect.push_back("Sources data path");
	comm_vect.push_back("input path");
	for(long iframe=0;iframe < samples_struct.ntotscan; iframe ++)
		comm_vect.push_back("Data Source Fits File");
	comm_vect.push_back("every scans have the same detector list which name is filled in this field");
	comm_vect.push_back("bolometers filelist suffix : can be void");
	comm_vect.push_back("temporary data path");
	comm_vect.push_back("output path");

	// [sanePos]
	comm_vect.push_back("SIZE OF THE PIXEL (deg)");
	comm_vect.push_back("flagged data are put in a separate map");
	comm_vect.push_back("gaps are projected to a pixel in the map, if so gap filling of noise only is performed iteratively");
	comm_vect.push_back("Sources files format");
	comm_vect.push_back("name of the fits file that was used to mask radiant sources");
	comm_vect.push_back("RA nominal");
	comm_vect.push_back("DEC nominal");
	comm_vect.push_back("Projection type");

	// [sanePre]
	comm_vect.push_back("Number of samples to apodize");
	comm_vect.push_back("sampling frequency (Hz)");
	comm_vect.push_back("Butterworth filter frequency (Hz)");
	comm_vect.push_back("file containing frequency under which noise power spectra are thresholded");
	comm_vect.push_back("polynomia order");
	comm_vect.push_back("Remove a polynomia fitted to the data to reduce fluctuations on timescales larger than the length of the considered segment ?");
	comm_vect.push_back("Do we remove a baseline from the data ?");
	comm_vect.push_back("Correlations between detectors are not included in the analysis ?");
	comm_vect.push_back("Do we fill the gaps in the timeline with White noise + baseline ?");

	// [saneInv]
	comm_vect.push_back("this file contains the matrix you want to invert (same for all scans)");
	comm_vect.push_back("covariance matrices path");
	comm_vect.push_back("Each scan has his own covmat named basename(scan_filename) + cov_suffix");

	// [sanePS]
	comm_vect.push_back("ell file (the same for all the scans)");
	comm_vect.push_back("each scan has an ell file named : basename(scan_filename) + ell_suffix");
	comm_vect.push_back("MixingMatrix file is the same for all the scans");
	comm_vect.push_back("each scan has a mixmat file named : basename(scan_filename) + MixingMatrix_suffix");
	comm_vect.push_back("number of noise component estimated in sanePS");
	comm_vect.push_back("A map whose signal will be subtracted from the data for the noise estimation step");
	comm_vect.push_back("Do we save the session after each step ?");

	// [sanePic]
	comm_vect.push_back("iterW = -1 -> no save on disk will be done. iterW > 0 -> save on disk every iterW iterations");
	comm_vect.push_back("Maximum number of conjugate gradient iterations");


	//  (commons)
	value_vect.push_back(dir.dirfile);
	value_vect.push_back(dir.input_dir);
	for(int num=0;num<(int)samples_struct.ntotscan;num++)
		value_vect.push_back(samples_struct.fitsvect[num]);
	value_vect.push_back(dir.bolo_global_filename);
	value_vect.push_back(dir.bolo_suffix);
	value_vect.push_back(dir.tmp_dir);
	value_vect.push_back(dir.output_dir);

	// [sanePos]
	value_vect.push_back(StringOf(pos_param.pixdeg));
	value_vect.push_back(StringOf(pos_param.flgdupl ? "True" : "False"));
	value_vect.push_back(StringOf(pos_param.fileFormat ? "HIPE" : "SANEPIC"));
	value_vect.push_back(StringOf(pos_param.projgaps ? "True" : "False"));
	value_vect.push_back(pos_param.maskfile);
	value_vect.push_back(StringOf(pos_param.ra_nom));
	value_vect.push_back(StringOf(pos_param.dec_nom));
	value_vect.push_back(pos_param.projtype);


	// [sanePre]
	value_vect.push_back(StringOf(proc_param.napod));
	value_vect.push_back(StringOf(proc_param.fsamp));
	value_vect.push_back(StringOf(proc_param.f_lp));
	value_vect.push_back(StringOf(proc_param.fcut_file));
	value_vect.push_back(StringOf(proc_param.poly_order));
	value_vect.push_back(StringOf(proc_param.remove_polynomia ? "True" : "False"));
	value_vect.push_back(StringOf(proc_param.NORMLIN ? "False" : "True"));
	value_vect.push_back(StringOf(proc_param.CORRon ? "True" : "False"));
	value_vect.push_back(StringOf(proc_param.NOFILLGAP ? "False" : "True"));

	// [saneInv]
	value_vect.push_back(Inv_param.cov_matrix_file);
	value_vect.push_back(Inv_param.noise_dir);
	value_vect.push_back(Inv_param.cov_matrix_suffix);

	// [sanePS]
	value_vect.push_back(PS_param.ell_global_file);
	value_vect.push_back(PS_param.ell_suffix);
	value_vect.push_back(PS_param.mix_global_file);
	value_vect.push_back(PS_param.mix_suffix);
	value_vect.push_back(StringOf(PS_param.ncomp));
	value_vect.push_back(StringOf(PS_param.signame));
	value_vect.push_back(StringOf(PS_param.save_data));

	// [sanePic]
	value_vect.push_back(StringOf(Pic_param.iterw));
	value_vect.push_back(StringOf(Pic_param.itermax));


	if (fits_open_file(&fptr, fname.c_str(), READWRITE, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}
	// ---------------------------------------------
	// write the Channel List

	char *ttype[] = { (char*) "KEY", (char*) "values", (char*) "comments" };
	char *tform[] = { tableFormat(key_vect), tableFormat(value_vect), tableFormat(comm_vect) };
	char *tunit[] = { (char*) "None", (char*) "None", (char*) "None" };
	char **data, **data2, **data_value;

	data = vString2carray(key_vect);
	data2 = vString2carray(comm_vect);
	data_value = vString2carray(value_vect);

	if (fits_create_tbl(fptr, BINARY_TBL, key_vect.size(), 3, ttype, tform, tunit,
			(char*)"History", &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 1, 1, 1, key_vect.size(), data, &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 2, 1, 1, value_vect.size(), data_value, &fits_status))
		return 1;
	if (fits_write_col(fptr, TSTRING, 3, 1, 1, comm_vect.size(), data2, &fits_status))
		return 1;
	if (fits_write_key(fptr, TSTRING, (char *) "TUNIT1", (char *) "NONE",
			(char *) "physical unit of the field", &fits_status))
		return 1;

	if (fits_write_chksum(fptr, &fits_status)){
		cout << "error checksum !\n";
		return 1;
	}

	// close file
	if(fits_close_file(fptr, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}


	for(long ii=0;ii<(long)key_vect.size();ii++){
		delete [] data[ii];
		delete [] data2[ii];
		delete [] data_value[ii];
	}

	delete [] data2;
	delete [] data_value;
	delete [] data;
	delete [] tform[0];
	delete [] tform[1];
	delete [] tform[2];


	return 0;
}


int write_fits_META(string fname, long NAXIS1, long NAXIS2, string path, struct param_sanePre proc_param, struct param_sanePos pos_param, std::vector<double> fcut, struct samples samples_struct, long ncomp){

	fitsfile *fp;
	int fits_status = 0; // MUST BE initialized... otherwise it fails on the call to the function...
	std::ostringstream oss;

	if (fits_open_file(&fp, fname.c_str(), READWRITE, &fits_status))
		fits_report_error(stderr, fits_status);

	fits_create_img(fp, 8, 0, 0, &fits_status);

	fits_update_key(fp, TSTRING, (char *)"EXTNAME", (void*)"History",
			(char *) "table name", &fits_status);

	string keyname = "PATHNAME";
	string value = path;
	string comm = "Source data path";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	for(int num=0;num<(int)samples_struct.ntotscan;num++){
		oss << "Source" << num;
		string keyname = oss.str();
		oss.str("");
		string value = samples_struct.fitsvect[num];
		string comm = "Data Source Fits File";
		if (fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
			fits_report_error(stderr, fits_status);
			return 1;
		}

	}


	oss << proc_param.napod;
	keyname = "NAPOD";
	value = oss.str();
	comm = "Number of samples to apodize";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}


	oss << proc_param.poly_order;
	keyname = "POLYORDR";
	value = oss.str();
	comm = "Fitted"; // polynomia order";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}


	// use only 8 characters for keyname to avoid HIERARCH keyword addition by cfitsio...
	oss << proc_param.fsamp;
	keyname = "SAMPFREQ";
	value = oss.str();
	comm = "sampling frequency (Hz)";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	oss << proc_param.f_lp;
	keyname = "FILTFREQ";
	value = oss.str();
	comm = "Butterworth filter frequency (Hz)";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}


	keyname = "FILLGAPS";
	if(proc_param.NOFILLGAP)
		value = "no";
	else
		value = "yes";
	comm = "Do we fill the gaps in the timeline with White noise + baseline ?";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	keyname = "NORMLIN";
	if(proc_param.NORMLIN)
		value = "no";
	else
		value = "yes";
	comm = "Do we remove a baseline from the data ?";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	keyname = "CORREL";
	if(proc_param.CORRon)
		value = "yes";
	else
		value = "no";
	comm = "Correlations between detectors are not included in the analysis ?";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	keyname = "POLYNSUB";
	if(proc_param.remove_polynomia)
		value = "yes";
	else
		value = "no";
	comm = "Remove a polynomia fitted to the data to reduce fluctuations on timescales larger than the length of the considered segment ?";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	oss << pos_param.pixdeg;
	keyname = "PIXSIZE";
	value = oss.str();
	comm = "SIZE OF THE PIXEL (deg)";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	keyname = "DUPLMAP";
	if(pos_param.flgdupl)
		value = "yes";
	else
		value = "no";
	comm = "flagged data are put in a separate map";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	keyname = "GAPSPROJ";
	if(pos_param.projgaps)
		value = "yes";
	else
		value = "no";
	comm = "gaps are projected to a pixel in the map, if so gap filling of noise only is performed iteratively";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	keyname = "FORMAT";
	if(pos_param.fileFormat)
		value = "HIPE";
	else
		value = "SANEPIC";
	comm = "Sources file format";

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	keyname = "MASKFILE";
	oss << pos_param.maskfile;
	value = oss.str();
	comm = "name of the fits file that was used to mask radiant sources";
	oss.str("");

	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}


	//	oss << det.ndet;
	//	keyname = "NUMDET";
	//	value = oss.str();
	//	comm = "Number of detectors that were used for the analysis";
	//	oss.str("");
	//
	//	if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
	//		fits_report_error(stderr, fits_status);
	//		return 1;
	//	}

	if(ncomp>0){
		oss << ncomp;
		keyname = "COMPONEN";
		value = oss.str();
		comm = "number of noise component to estimate in sanePS";
		oss.str("");

		if ( fits_write_key(fp, TSTRING, (char*) keyname.c_str(), (char*)value.c_str(), (char*)comm.c_str(), &fits_status)){
			fits_report_error(stderr, fits_status);
			return 1;
		}
	}


	// close file
	if(fits_close_file(fp, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	return 0;
}




int write_fits_mask(std::string fname, std::string maskfile)
{

	fitsfile *fptr, *outfptr;
	int fits_status = 0;
	long naxes[2] = { 1, 1 };

	if (fits_open_file(&fptr, maskfile.c_str(), READWRITE, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	if (fits_open_file(&outfptr, fname.c_str(), READWRITE, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	if(fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", 0, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}


	if(fits_copy_header(fptr, outfptr, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	if (fits_write_chksum(fptr, &fits_status)){
		cout << "error checksum !\n";
		return 1;
	}

	// Retrieve the image size
	if (fits_get_img_size(fptr, 2, naxes, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}


	//	long NAXIS1 = naxes[0];

	if(fits_movnam_hdu(outfptr, IMAGE_HDU, (char*) "mask", 0, &fits_status))
		if (fits_copy_data(fptr, outfptr, &fits_status)){
			fits_report_error(stderr, fits_status);
			return 1;
		}


	// close files
	if(fits_close_file(fptr, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	if(fits_close_file(outfptr, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	return 0;

}

int read_mask_wcs(string fname, string extname, struct wcsprm *& wcs, long &NAXIS1, long &NAXIS2,  short *& data)
/*
 * Read the extension 'extname' from the 'fname' fits file, extension must be a 2D image
 * Return a wcs structure, the size of the image, the image data type and image itself cast to int, float or double
 */
{
	fitsfile *fptr;
	int status = 0, anynul, wcsstatus[NWCSFIX];
	char *header;
	int nkeyrec, nwcs, nreject;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };

	// Open the fits file...
	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ... and move to the 'extname'
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) extname.c_str(), 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// Retrieve the image size
	if (fits_get_img_size(fptr, 2, naxes, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	NAXIS1 = naxes[0];
	NAXIS2 = naxes[1];
	// Allocate the image container and read its depending on the type
	//	switch (dtype) {
	//	case 's':
	data = new short[NAXIS1*NAXIS2];
	if (fits_read_pix(fptr, TSHORT, fpixel, (long long) NAXIS1*NAXIS2, 0, data, &anynul, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// retrieve the wcs header
	if (fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	if (( status = wcspih(header, nkeyrec, WCSHDR_all, 2, &nreject, &nwcs, &wcs))) {
		fprintf(stderr, "wcspih ERROR %d: %s.\n", status,wcshdr_errmsg[status]);
	}

	if ((status = wcsfix(7, 0, wcs, wcsstatus))) {
		for (long ii = 0; ii < NWCSFIX; ii++) {
			if (wcsstatus[ii] > 0) {
				fprintf(stderr, "wcsfix ERROR %d: %s.\n", status, wcsfix_errmsg[wcsstatus[ii]]);
			}
		}
	}

	if ((status= wcsset(wcs))) {
		printf("wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
	}

	if (fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	free(header);
	return(0);
}


int read_fits_signal(string fname, double *S, long long* indpix, long NAXIS1, long NAXIS2, struct wcsprm * wcs)
/*
 * This function read the sanePic generated map and converts it into S (only seen pixels)
 */
{
	fitsfile *fptr;
	int status = 0, anynul, wcsstatus[NWCSFIX];
	char *header;
	int nkeyrec, nwcs, nreject;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	long imNAXIS1, imNAXIS2;
	long mi;
	double *map;
	struct wcsprm *wcs_fits;

	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*)"Image", 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// retrieve the wcs header
	if (fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	if (( status = wcspih(header, nkeyrec, WCSHDR_all, 2, &nreject, &nwcs, &wcs_fits))) {
		fprintf(stderr, "wcspih ERROR %d: %s.\n", status,wcshdr_errmsg[status]);
	}

	if ((status = wcsfix(7, 0, wcs_fits, wcsstatus))) {
		for (long ii = 0; ii < NWCSFIX; ii++) {
			if (wcsstatus[ii] > 0) {
				fprintf(stderr, "wcsfix ERROR %d: %s.\n", status, wcsfix_errmsg[wcsstatus[ii]]);
			}
		}
	}

	if ((status= wcsset(wcs_fits))) {
		printf("wcsset ERROR %d: %s.\n", status, wcs_errmsg[status]);
	}


	if(fits_get_img_size(fptr, 2, naxes, &status))
		return 1;

	imNAXIS1=(long)naxes[0];
	imNAXIS2=(long)naxes[1];

	if(compare_wcs(fname, wcs, wcs_fits, NAXIS1, NAXIS2, imNAXIS1, imNAXIS2))
		return 1;

	// Initialize the data container
	map = new double [imNAXIS1*imNAXIS2];

	if (fits_read_pix(fptr, TDOUBLE, fpixel, (long long) imNAXIS1*imNAXIS2, 0, map, &anynul, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// work correctly
	for (long ii=0; ii<imNAXIS1; ii++) {
		for (long jj=0; jj<imNAXIS2; jj++) {
			mi = jj*imNAXIS1 + ii;
			if (indpix[mi] >= 0){
				S[indpix[mi]]= map[mi];
			}
		}
	}

	if (fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	//	nwcs=1;
	wcsvfree(&nwcs, &wcs_fits);
	free(header);
	return 0;
}


int save_keyrec(string tmpdir, struct wcsprm * wcs, long NAXIS1, long NAXIS2){

	FILE *fout;
	int nkeyrec, status;
	char *header, *hptr;

	if ((status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header))) {
		printf("%4d: %s.\n", status, wcs_errmsg[status]);
		return 1;
	}


	tmpdir=tmpdir + "mapHeader.keyrec";
	fout = fopen(tmpdir.c_str(),"w");
	if (fout==NULL) {fputs ("Creation error : File error on mapHeader.keyrec\n",stderr); return (1);}

	fprintf(fout,"NAXIS1  = %20ld / %-47s\n",NAXIS1,"length of data axis 1");
	fprintf(fout,"NAXIS2  = %20ld / %-47s\n",NAXIS2,"length of data axis 2");

	hptr = header;
	for (int i = 0; i < nkeyrec; i++, hptr += 80) {
		fprintf(fout,"%.80s\n", hptr);
	}
	fclose(fout);
	free(header);

	return 0;
}

int print_MapHeader(struct wcsprm *wcs){

	int nkeyrec;
	char * header, *hptr ;
	if (int status = wcshdo(WCSHDO_all, wcs, &nkeyrec, &header)) {
		printf("%4d: %s.\n", status, wcs_errmsg[status]);
		return 1;
	}
	hptr = header;
	printf("\n\n Map Header :\n");
	for (int ii = 0; ii < nkeyrec; ii++, hptr += 80) {
		printf("%.80s\n", hptr);
	}
	free(header);

	return 0;
}

int read_keyrec(string tmpdir, struct wcsprm * & wcs, long * NAXIS1, long * NAXIS2, int rank){

	char *memblock=NULL;
	int nkeyrec=0, nreject, nwcs, status;

	if(rank==0){
		tmpdir = tmpdir + "mapHeader.keyrec";

		FILE *fin;
		size_t result;
		int size;

		fin = fopen(tmpdir.c_str(),"r");
		if (fin==NULL) {fputs ("Read error : File error on mapHeader.keyrec",stderr); return 1;}

		fseek(fin, 0L, SEEK_END);     /* Position to end of file */
		size = ftell(fin);            /* Get file length */
		rewind(fin);                  /* Back to start of file */


		nkeyrec = size/81;

		char comment[47];

		// Read the two first lines, NAXIS1/NAXIS2
		result = fscanf(fin,"NAXIS1  = %20ld / %47c\n",NAXIS1,(char *) &comment);
		result = fscanf(fin,"NAXIS2  = %20ld / %47c\n",NAXIS2,(char *) &comment);


		memblock = new char [(nkeyrec-2)*80];

		for (int ii = 0; ii < nkeyrec; ii++) {
			result = fread(&memblock[ii*80], 80, sizeof(char), fin);
			fseek(fin, 1, SEEK_CUR); // skip newline char
		}
		fclose (fin);
	}

#ifdef USE_MPI
	//	int position=0;
	MPI_Bcast(NAXIS1,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(NAXIS2,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(&nkeyrec,1,MPI_INT,0,MPI_COMM_WORLD);

	if(rank!=0)
		memblock = new char [(nkeyrec-2)*80];

	MPI_Bcast(memblock, (nkeyrec-2)*80, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

	/* Parse the primary header of the FITS file. */
	/* -2 to handle the firts two NAXIS? keyword */
	if ((status = wcspih(memblock, nkeyrec-2, WCSHDR_all, 2, &nreject, &nwcs, &wcs))) {
		fprintf(stderr, "wcspih ERROR %d: %s.\n", status,wcshdr_errmsg[status]);
		return 1;
	}
	delete[] memblock;
	//	free(comment);

	return 0;
}

int compare_wcs(std::string fname, struct wcsprm *wcs, struct wcsprm *wcs_fits, long NAXIS1, long NAXIS2, long imNAXIS1, long imNAXIS2){

	// compatibility verifications !

	if(wcs_fits->crpix[0]!=wcs->crpix[0] || wcs_fits->crpix[1]!=wcs->crpix[1]){
		cout << "CRPIX are different between mapheader.keyrec and the map_file : " <<  fname << endl;
		cout << "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}
	if(wcs_fits->crval[0]!=wcs->crval[0] || wcs_fits->crval[1]!=wcs->crval[1]){
		cout << "CRVAL are different between mapheader.keyrec and the map_file : " <<  fname << endl;
		cout << "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}
	if(wcs_fits->cdelt[0]!=wcs->cdelt[0] || wcs_fits->cdelt[1]!=wcs->cdelt[1]){
		cout << "CDELT are different between mapheader.keyrec and the map_file : " <<  fname << endl;
		cout << "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}
	if(wcs_fits->lonpole!=wcs->lonpole || wcs_fits->latpole!=wcs->latpole){
		cout << "LONPOLE and/or LATPOLE are different between mapheader.keyrec and the map_file : " <<  fname << endl;
		cout << "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}
	if(NAXIS1!=imNAXIS1 || NAXIS2!=imNAXIS2){
		cout << "NAXIS are different between mapheader.keyrec and the map_file : " <<  fname << endl;
		cout << "Please re-run sanePos with the correct data/ini_file. Exiting...\n";
		return 1;
	}


	return 0;

}

