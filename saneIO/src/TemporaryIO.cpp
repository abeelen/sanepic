#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "StructDefinition.h"
#include "InputFileIO.h"
#include "TemporaryIO.h"
#include "DataIO.h"
#include "Utilities.h"
#include "ErrorCode.h"

extern "C" {
#include "getdata.h"
#include <wcslib/cel.h>
#include <wcslib/wcs.h>
#include <wcslib/sph.h>
#include <wcslib/wcsmath.h>
#include <wcslib/wcstrig.h>
#include <wcslib/prj.h>
#include <fitsio.h>
}

#ifdef USE_MPI
#include "mpi.h"
#endif

using namespace std;

int writeDataFlagToDirfile(struct param_common dir,
		struct samples samples_struct) {

	double *signal;
	int *flag;
	long ns;
	DIRFILE* D, *H;
	std::vector<string> det_vect;

	fitsfile *fptr;
	int status = 0;
	int rowIndex = 0, naxis = 0, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };

	for (long iframe = samples_struct.iframe_min;
			iframe < samples_struct.iframe_max; iframe++) {

		det_vect = samples_struct.bolo_list[iframe];

		string base_name = samples_struct.basevect[iframe];
		string fits_name = samples_struct.fitsvect[iframe];

		string datadir = dir.tmp_dir + "dirfile/" + base_name + "/Data";
		string flagdir = dir.tmp_dir + "dirfile/" + base_name + "/Flag";

		D = gd_open((char *) datadir.c_str(),
				GD_RDWR | GD_CREAT | GD_VERBOSE | GD_UNENCODED);
		H = gd_open((char *) flagdir.c_str(),
				GD_RDWR | GD_CREAT | GD_VERBOSE | GD_UNENCODED);

		if (fits_open_file(&fptr, fits_name.c_str(), READONLY, &status)) {
			fits_report_error(stderr, status);
			return 1;
		}

		//TODO: if the dirfile structure is set by sub_rank 0 in a first loop, then this loop can be parallelized by sub_rank
		for (long idet = 0; idet < (long) det_vect.size(); idet++) {

			string field = det_vect[idet];
			string data_outfile = "Data_" + base_name + "_" + field;
			string flag_outfile = "Flag_" + base_name + "_" + field;

			// Retrieve the row of the specified channel
			rowIndex = find_channel_index(fptr, field.c_str());

			// Move ptr to signal hdu
			if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", 0,
					&status)) {
				fits_report_error(stderr, status);
				return 1;
			}

			// Retrieve the size of the signal
			if (fits_get_img_dim(fptr, &naxis, &status)) {
				fits_report_error(stderr, status);
				return 1;
			}
			if (naxis != 2) {
				fits_report_error(stderr, status);
				return 1;
			}
			if (fits_get_img_size(fptr, 2, naxes, &status)) {
				fits_report_error(stderr, status);
				return 1;
			}

			// Allocate Memory
			ns = naxes[0];
			signal = new double[ns];
			flag = new int[ns];

			// Retrieve the corresponding row
			fpixel[0] = 1;
			fpixel[1] = rowIndex;
			if (fits_read_pix(fptr, TDOUBLE, fpixel, ns, 0, signal, &anynul,
					&status)) {
				fits_report_error(stderr, status);
				return 1;
			}

			// Move ptr to mask hdu ...
			if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", 0, &status)) {
				fits_report_error(stderr, status);
				return 1;
			}

			// Retrieve the size of the mask
			if (fits_get_img_dim(fptr, &naxis, &status)) {
				fits_report_error(stderr, status);
				return 1;
			}
			if (naxis != 2) {
				fits_report_error(stderr, status);
				return 1;
			}
			if (fits_get_img_size(fptr, 2, naxes, &status)) {
				fits_report_error(stderr, status);
				return 1;
			}
			if (ns != naxes[0])
				return 1;

			if (fits_read_pix(fptr, TINT, fpixel, ns, 0, flag, &anynul,
					&status)) {
				fits_report_error(stderr, status);
				return 1;
			}

			//configure dirfile field
			gd_entry_t E;
			E.field = (char*) data_outfile.c_str();
			E.field_type = GD_RAW_ENTRY;
			E.fragment_index = 0;
			E.spf = 1;
			E.data_type = GD_DOUBLE;
			E.scalar[0] = NULL;

			// add to the dirfile
			gd_add(D, &E);
			//			gd_flush(D,NULL);

			// write binary file on disk
			int n_write = gd_putdata(D, (char*) data_outfile.c_str(), 0, 0, 0,
					ns, GD_DOUBLE, signal);
			if (gd_error(D) != 0) {
				cout << "error putdata in write_data_flag : wrote " << n_write
						<< " and expected " << ns << endl;
				return 1;
			}

			//configure dirfile field
			gd_entry_t F;
			F.field = (char*) flag_outfile.c_str();
			F.field_type = GD_RAW_ENTRY;
			F.fragment_index = 0;
			F.spf = 1;
			F.data_type = GD_INT32;
			F.scalar[0] = NULL;

			// add to the dirfile
			gd_add(H, &F);

			// write binary file on disk
			n_write = gd_putdata(H, (char*) flag_outfile.c_str(), 0, 0, 0, ns,
					GD_INT32, flag);
			if (gd_error(H) != 0) {
				cerr << "EE - error putdata in write_data_flag : wrote "
						<< n_write << " and expected " << ns << endl;
				return 1;
			}

			delete[] flag;
			delete[] signal;

			gd_flush(D, data_outfile.c_str());
			gd_flush(H, flag_outfile.c_str());
		}

		// close fits file
		if (fits_close_file(fptr, &status)) {
			fits_report_error(stderr, status);
			return 1;
		}

		// close dirfile
		if (gd_close(D)) {
			cout << "Dirfile gd_close error in write_data_flag for : "
					<< datadir << endl;
			return 1;
		}

		// close dirfile
		if (gd_close(H)) {
			cout << "Dirfile gd_close error in write_data_flag for : "
					<< flagdir << endl;
			return 1;
		}

	}

	return 0;
}

int writeLonLatToDirfile(struct param_common dir,
		struct samples samples_struct) {

	double *lon;
	double *lat;
	long ns;
	DIRFILE* D, *H;

	fitsfile *fptr;
	int status = 0;
	int rowIndex = 0, naxis = 0, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };

	std::vector<string> det_vect;

	for (long iframe = samples_struct.iframe_min;
			iframe < samples_struct.iframe_max; iframe++) {

		det_vect = samples_struct.bolo_list[iframe];

		string base_name = samples_struct.basevect[iframe];
		string fits_name = samples_struct.fitsvect[iframe];

		string LONdir = dir.tmp_dir + "dirfile/" + base_name + "/Lon";
		string LATdir = dir.tmp_dir + "dirfile/" + base_name + "/Lat";

		D = gd_open((char *) LONdir.c_str(),
				GD_RDWR | GD_CREAT | GD_VERBOSE | GD_UNENCODED); // | GD_TRUNC |
		H = gd_open((char *) LATdir.c_str(),
				GD_RDWR | GD_CREAT | GD_VERBOSE | GD_UNENCODED); // | GD_TRUNC |

		if (fits_open_file(&fptr, fits_name.c_str(), READONLY, &status)) {
			fits_report_error(stderr, status);
			return 1;
		}

		//TODO: if the dirfile structure is set by sub_rank 0 in a first loop, then this loop can be parallelized by sub_rank
		for (long idet = 0; idet < (long) det_vect.size(); idet++) {

			string field = det_vect[idet];
			string lon_outfile = "Lon_" + base_name + "_" + field;
			string lat_outfile = "Lat_" + base_name + "_" + field;

			// Retrieve the row of the specified channel
			rowIndex = find_channel_index(fptr, field.c_str());

			// Move ptr to lon hdu
			if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "lon", 0, &status)) {
				fits_report_error(stderr, status);
				return 1;
			}

			// Retrieve the size of the signal
			if (fits_get_img_dim(fptr, &naxis, &status)) {
				fits_report_error(stderr, status);
				return 1;
			}
			if (naxis != 2) {
				fits_report_error(stderr, status);
				return 1;
			}
			if (fits_get_img_size(fptr, 2, naxes, &status)) {
				fits_report_error(stderr, status);
				return 1;
			}

			// Allocate Memory
			ns = naxes[0];
			lon = new double[ns];
			lat = new double[ns];

			// ---------------------------------------------
			// Retrieve the corresponding row
			fpixel[0] = 1;
			fpixel[1] = rowIndex;
			if (fits_read_pix(fptr, TDOUBLE, fpixel, ns, 0, lon, &anynul,
					&status)) {
				fits_report_error(stderr, status);
				return 1;
			}

			// ---------------------------------------------
			// Move ptr to lat hdu
			if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "lat", 0, &status)) {
				fits_report_error(stderr, status);
				return 1;
			}

			// Retrieve the size of the signal
			if (fits_get_img_dim(fptr, &naxis, &status)) {
				fits_report_error(stderr, status);
				return 1;
			}
			if (naxis != 2) {
				fits_report_error(stderr, status);
				return 1;
			}
			if (fits_get_img_size(fptr, 2, naxes, &status)) {
				fits_report_error(stderr, status);
				return 1;
			}
			if (ns != naxes[0])
				return 1;

			// ---------------------------------------------
			// Retrieve the corresponding row
			if (fits_read_pix(fptr, TDOUBLE, fpixel, ns, 0, lat, &anynul,
					&status)) {
				fits_report_error(stderr, status);
				return 1;
			}

			//configure dirfile field
			gd_entry_t E;
			E.field = (char*) lon_outfile.c_str();
			E.field_type = GD_RAW_ENTRY;
			E.fragment_index = 0;
			E.spf = 1;
			E.data_type = GD_DOUBLE;
			E.scalar[0] = NULL;

			// add to the dirfile
			gd_add(D, &E);

			// write binary file on disk
			int n_write = gd_putdata(D, (char*) lon_outfile.c_str(), 0, 0, 0,
					ns, GD_DOUBLE, lon);
			if (gd_error(D) != 0) {
				cout << "error putdata in write_lon : wrote " << n_write
						<< " and expected " << ns << endl;
				return 1;
			}

			//configure dirfile field
			gd_entry_t F;
			F.field = (char*) lat_outfile.c_str();
			F.field_type = GD_RAW_ENTRY;
			F.fragment_index = 0;
			F.spf = 1;
			F.data_type = GD_DOUBLE;
			F.scalar[0] = NULL;

			// add to the dirfile
			gd_add(H, &F);

			// write binary file on disk
			n_write = gd_putdata(H, (char*) lat_outfile.c_str(), 0, 0, 0, ns,
					GD_DOUBLE, lat);
			if (gd_error(H) != 0) {
				cout << "error putdata in write_lat : wrote " << n_write
						<< " and expected " << ns << endl;
				return 1;
			}

			delete[] lon;
			delete[] lat;

			gd_flush(D, E.field);
			gd_flush(H, F.field);
		}

		// close fits file
		if (fits_close_file(fptr, &status)) {
			fits_report_error(stderr, status);
			return 1;
		}

		// close dirfile
		if (gd_close(D)) {
			cout << "Dirfile gd_close error  " << LONdir << endl;
			return 1;
		}

		// close dirfile
		if (gd_close(H)) {
			cout << "Dirfile gd_close error for : " << LATdir << endl;
			return 1;
		}

	}

	return 0;
}

int exportLonLatToDirfile(struct param_common dir,
		struct samples samples_struct) {

	DIRFILE* D, *H;

	std::vector<string> det_vect;

	for (long iframe = samples_struct.iframe_min;
			iframe < samples_struct.iframe_max; iframe++) {

		det_vect = samples_struct.bolo_list[iframe];

		string base_name = samples_struct.basevect[iframe];
		string fits_name = samples_struct.fitsvect[iframe];

		string LONdir = dir.tmp_dir + "dirfile/" + base_name + "/Lon";
		string LATdir = dir.tmp_dir + "dirfile/" + base_name + "/Lat";

		D = gd_open((char *) LONdir.c_str(),
				GD_RDWR | GD_CREAT | GD_VERBOSE | GD_UNENCODED); // | GD_TRUNC |
		H = gd_open((char *) LATdir.c_str(),
				GD_RDWR | GD_CREAT | GD_VERBOSE | GD_UNENCODED); // | GD_TRUNC |

		double *reflon, *reflat, *phi, **offsets;

		long ns = samples_struct.nsamples[iframe];

		// read bolo offsets
		// TODO : This function should also return the PRJCODE to be used below...
		if (read_all_bolo_offsets_from_fits(fits_name, det_vect, offsets))
			return 1;

		// read reference position
		long test_ns;
		if (read_ReferencePosition_from_fits(fits_name, reflon, reflat, phi,
				test_ns))
			return 1;

		if (test_ns != ns) {
			cerr << "Read position does not correspond to frame position"
					<< endl;
			cerr << "Check !!" << endl;
			return 1;
		}

		// Save some computation time...
		double *cosphi, *sinphi;
		cosphi = new double[ns];
		sinphi = new double[ns];

		for (long ii = 0; ii < ns; ii++) {
			cosphi[ii] = cos(phi[ii] / 180.0 * M_PI);
			sinphi[ii] = sin(phi[ii] / 180.0 * M_PI);
		}

		// find the pointing solution at each time stamp for each detector
		struct celprm arrayProj;
		celini(&arrayProj);

		// TODO: use the PRJCODE read from the file...
		tanset(&arrayProj.prj);

		//TODO Looping on detector is absolutely inefficient here...
		//     as one need to initialize projection center by time stamp
		//     and therefore we can not parallelize call to wcslib

		double *lon, *lat;
		lon = new double[ns];
		lat = new double[ns];

		for (long idet = 0; idet < (long) det_vect.size(); idet++) {

			string field = det_vect[idet];
			string lon_outfile = "Lon_" + base_name + "_" + field;
			string lat_outfile = "Lat_" + base_name + "_" + field;

			double offxx, offyy, dummy1, dummy2, lon_deg, lat_deg;
			int status;

			// Deproject the detector offset.
			for (long ii = 0; ii < ns; ii++) {

				arrayProj.ref[0] = reflon[ii];
				arrayProj.ref[1] = reflat[ii];

				if (celset(&arrayProj))
					cerr << "problem celset\n";

				offxx = (cosphi[ii] * offsets[idet][0]
				                                    - sinphi[ii] * offsets[idet][1]) * -1;
				offyy = sinphi[ii] * offsets[idet][0]
				                                   + cosphi[ii] * offsets[idet][1];

				// Projection away from the detector plane, into the spherical sky...
				if (celx2s(&arrayProj, 1, 0, 1, 1, &offxx, &offyy, &dummy1,
						&dummy2, &lon_deg, &lat_deg, &status) == 1) {
					printf("   TAN(X2S) ERROR 1: %s\n", prj_errmsg[1]);
					continue;
				}
				lon[ii] = lon_deg;
				lat[ii] = lat_deg;

			}

			//configure dirfile field
			gd_entry_t E;
			E.field = (char*) lon_outfile.c_str();
			E.field_type = GD_RAW_ENTRY;
			E.fragment_index = 0;
			E.spf = 1;
			E.data_type = GD_DOUBLE;
			E.scalar[0] = NULL;

			// add to the dirfile
			gd_add(D, &E);

			// write binary file on disk
			int n_write = gd_putdata(D, (char*) lon_outfile.c_str(), 0, 0, 0,
					ns, GD_DOUBLE, lon);
			if (gd_error(D) != 0) {
				cout << "error putdata in write_lon : wrote " << n_write
						<< " and expected " << ns << endl;
				return 1;
			}

			//configure dirfile field
			gd_entry_t F;
			F.field = (char*) lat_outfile.c_str();
			F.field_type = GD_RAW_ENTRY;
			F.fragment_index = 0;
			F.spf = 1;
			F.data_type = GD_DOUBLE;
			F.scalar[0] = NULL;

			// add to the dirfile
			gd_add(H, &F);

			// write binary file on disk
			n_write = gd_putdata(H, (char*) lat_outfile.c_str(), 0, 0, 0, ns,
					GD_DOUBLE, lat);
			if (gd_error(H) != 0) {
				cout << "error putdata in write_lat : wrote " << n_write
						<< " and expected " << ns << endl;
				return 1;
			}

			gd_flush(D, lon_outfile.c_str());
			gd_flush(H, lat_outfile.c_str());
		}

		delete[] lon;
		delete[] lat;

		// close dirfile
		if (gd_close(D)) {
			cout << "Dirfile gd_close error  " << LONdir << endl;
			return 1;
		}

		// close dirfile
		if (gd_close(H)) {
			cout << "Dirfile gd_close error for : " << LATdir << endl;
			return 1;
		}

	}

	return 0;
}

int readDataFromDirfile(DIRFILE* D, string scan_name, string field,
		double *data, long ns) {

	// set dirfile name and binary name
	string outfile = "Data_" + scan_name + "_" + field;
	const char * field_code;
	field_code = outfile.c_str();

	// fill data array
	int nget = gd_getdata(D, field_code, 0, 0, 0, ns, GD_DOUBLE, data);
	if (gd_error(D) != 0) {
		cout << "error getdata in read_data_from_dirfile : read " << nget
				<< endl;
		return 1;
	}

	gd_flush(D, field_code);

	return 0;
}

int readFlagFromDirfile(DIRFILE* H, string scan_name, string field, int *mask,
		long ns) {

	// set dirfile name and binary name
	string outfile = "Flag_" + scan_name + "_" + field;
	const char * field_code;
	field_code = outfile.c_str();

	// fill mask array
	int nget = gd_getdata(H, field_code, 0, 0, 0, ns, GD_INT32, mask);
	if (gd_error(H) != 0) {
		cout << "error getdata in read_flag_from_dirfile : read " << nget
				<< endl;
		return 1;
	}

	gd_flush(H, field_code);

	return 0;
}

int readLonFromDirfile(DIRFILE* D, string scan_name, string field, double *lon,
		long ns) {

	// set dirfile name and binary name
	string outfile = "Lon_" + scan_name + "_" + field;
	const char * field_code;
	field_code = outfile.c_str();

	// fill ra array
	int nget = gd_getdata(D, field_code, 0, 0, 0, ns, GD_DOUBLE, lon);
	if (gd_error(D) != 0) {
		cout << "error getdata in read_LON_from_dirfile : read " << nget
				<< endl;
		return 1;
	}

	gd_flush(D, field_code);

	return 0;
}

int readLatFromDirfile(DIRFILE* D, string scan_name, string field, double *lat,
		long ns) {

	// set dirfile name and binary name
	string outfile = "Lat_" + scan_name + "_" + field;
	const char * field_code;
	field_code = outfile.c_str();

	// fill lat array
	int nget = gd_getdata(D, field_code, 0, 0, 0, ns, GD_DOUBLE, lat);
	if (gd_error(D) != 0) {
		cout << "error getdata in read_LAT_from_dirfile : read " << nget
				<< endl;
		return 1;
	}

	gd_flush(D, field_code);

	return 0;
}

int writeSampToPix(DIRFILE *D, string scan_name, std::string boloname, long ns,
		long long *&samptopix)
/*!  write a sample to pixel vector to disk  */
{
	string outfile = "Indexes_" + scan_name + "_" + boloname;

	// write binary file on disk
	int n_write = gd_putdata(D, (char*) outfile.c_str(), 0, 0, 0, ns, GD_INT64,
			samptopix);
	if (gd_error(D) != 0) {
		cout << "error putdata in write_samptopix : wrote " << n_write
				<< " and expected " << ns << endl;
		return 1;
	}

	gd_flush(D, NULL);

	return 0;
}

int readSampToPix(DIRFILE* D, string scan_name, std::string boloname,
		long long *samptopix, long ns)
/*!  read a sample to pixel vector from disk  */
{
	// set binary file name
	string outfile = "Indexes_" + scan_name + "_" + boloname;
	const char * field_code;
	field_code = outfile.c_str();

	// fill samptopix array
	int nget = gd_getdata(D, field_code, 0, 0, 0, ns, GD_INT64, samptopix);
	if (gd_error(D) != 0) {
		cout << "error getdata in read_samptopix : read " << nget << endl;
		return 1;
	}

	gd_flush(D, field_code);

	return 0;
}

int writeIndexCCR(long long map_size, long long npixsrc, long long * indpsrc,
		std::string outdir)
/*!  write the bright sources index to disk  */
{
	FILE *fp;
	string testfile = outdir + "indpsrc.bin";

	if ((fp = fopen(testfile.c_str(), "w")) != NULL) {
		fwrite(&map_size, sizeof(long long), 1, fp);
		fwrite(&npixsrc, sizeof(long long), 1, fp);
		fwrite(indpsrc, sizeof(long long), map_size, fp);
		fclose(fp);
	} else {
		cerr << "EE - Could not open " << testfile << endl;
		return 1;
	}

	return 0;
}

int readIndexCCR(long long &map_size, long long &npixsrc, long long *&indpsrc, 	std::string outdir)
/*!  read the bright sources index from disk  */
{
	FILE *fp;
	string testfile;

	testfile = outdir + "indpsrc.bin";
	if ((fp = fopen(testfile.c_str(), "r")) != NULL) {
		if (1 != fread(&map_size, sizeof(long long), 1, fp)) {
			cerr << endl << "EE - failed reading map_size in indpsrc.bin file" << endl;
		}
		if (1 != fread(&npixsrc, sizeof(long long), 1, fp)) {
			cerr << endl << "EE - failed reading npixsrc in indpsrc.bin file" << endl;
		}
		indpsrc = new long long[map_size];
		if ((unsigned int) map_size != fread(indpsrc, sizeof(long long), map_size, fp)) {
			cerr << endl << "EE - failed reading indpsrc in indpsrc.bin file" << endl;
		}
		fclose(fp);
	} else {
		cerr << endl << "EE - cannot find indpsrc.bin file at " << testfile << endl;
		return 1;
	}

	return 0;
}

int write_indpix(long long ind_size, long long npix, long long *indpix,
		string outdir, int flagon)
/*!  write the map index to disk  */
{
	FILE *fp;
	string testfile2;

	testfile2 = outdir + "Indpix_for_conj_grad.bin";

	if ((fp = fopen(testfile2.c_str(), "w")) != NULL) {
		fwrite(&flagon, sizeof(int), 1, fp); // mat 04/06
		fwrite(&npix, sizeof(long long), 1, fp);
		fwrite(&ind_size, sizeof(long long), 1, fp);
		fwrite(indpix, sizeof(long long), ind_size, fp);
		fclose(fp);
	} else {
		cerr << "ERROR : Could not open " << testfile2 << endl;
		return 1;
	}

#ifdef DEBUG
	testfile2 = outdir + "Indpix_for_conj_grad.txt";
	if((fp = fopen(testfile2.c_str(),"w"))) {
		fprintf(fp,"%d\n",flagon);
		fprintf(fp,"%lld\n",npix);
		fprintf(fp,"%lld\n",ind_size);
		for(long ii =0;ii<ind_size;ii++)
			fprintf(fp,"%lld ",indpix[ii]);
		fclose(fp);
	} else {
		cerr << "ERROR : Could not find " << testfile2 << endl;
		return 1;
	}

#endif

	return 0;
}

int read_indpix(long long &ind_size, long long &npix, long long *&indpix,
		string outdir, int &flagon)
/*!  read the map index from disk  */
{
	FILE *fp;
	string testfile2;

	testfile2 = outdir + "Indpix_for_conj_grad.bin";
	if ((fp = fopen(testfile2.c_str(), "r")) != NULL) {
		if (1 != fread(&flagon, sizeof(int), 1, fp))
			cerr << "EE - failed reading flagon in Indpix file" << endl;

		if (1 != fread(&npix, sizeof(long long), 1, fp))
			cerr << "EE - failed reading npix in Indpix file" << endl;

		if (1 != fread(&ind_size, sizeof(long long), 1, fp))
			cerr << "EE - failed reading ind_size in Indpix file" << endl;

		indpix = new long long[ind_size];
		if ((unsigned int) ind_size
				!= fread(indpix, sizeof(long long), ind_size, fp))
			cerr << "EE - failed reading indpix in Indpix file" << endl;

		fclose(fp);
	} else {
		cerr << "Error : cannot find Indpix file " << testfile2 << endl;
		return 1;
	}

	return 0;
}

int write_PNd(double *PNd, long long npix, string outdir, std::string filename)
/*!  write the map preconditioner to disk  */
{
	FILE *fp;
	string testfile2;

	testfile2 = outdir + filename;

	if ((fp = fopen(testfile2.c_str(), "w")) != NULL) {
		fwrite(&npix, sizeof(long long), 1, fp);
		fwrite(PNd, sizeof(double), npix, fp);
		fclose(fp);
	} else {
		cerr << "ERROR : Could not find " << testfile2 << endl;
		return 1;
	}

#ifdef DEBUG
	testfile2 = outdir + "PNdCorr.txt";

	if((fp = fopen(testfile2.c_str(),"w"))) {
		fprintf(fp,"%lld\n",npix);
		for(int ii= 0; ii<npix;ii++)
			fprintf(fp,"%lf ",PNd[ii]);
		fclose(fp);
	} else {
		cerr << "ERROR : Could not find " << testfile2 << endl;
		return 1;
	}

#endif
	return 0;
}

// correlation between npix here and npix in Indpix file is done in sanePic (main.cpp)
int read_PNd(double *&PNdtot, long long &npix, string outdir,
		std::string filename)
/*!  read the map preconditioner from disk  */
{
	FILE *fp;
	string testfile2;

	testfile2 = outdir + filename;

	if ((fp = fopen(testfile2.c_str(), "r")) != NULL) {
		if (1 != fread(&npix, sizeof(long long), 1, fp))
			cerr << "EE - failed reading npix in file " << testfile2 << endl;
		PNdtot = new double[npix];
		if ((unsigned) npix != fread(PNdtot, sizeof(double), npix, fp))
			cerr << "EE - failed reading PNdtot in file " << testfile2 << endl;
		fclose(fp);
	} else {
		cerr << "EE - Unable to read file : " << testfile2 << endl;
		return 1;
	}

	return 0;
}

int writeFdata(DIRFILE *D, long ns, fftw_complex *fdata, string prefixe,
		long idet, string filename, std::vector<std::string> bolonames)
/*! write Fourier data file to disk */
{

	// set dirfile and binary names
	string outfile = prefixe + filename + "_" + bolonames[idet];
	const char * field_code;
	field_code = outfile.c_str();

	// write binary file on disk
	int n_write = gd_putdata(D, field_code, 0, 0, 0, ns / 2 + 1, GD_COMPLEX128,
			fdata);
	if (gd_error(D) != 0) {
		cout << "error putdata in write_fdata : wrote " << n_write
				<< " and expected " << ns << endl;
		return 1;
	}

	gd_flush(D, field_code);

	return 0;
}

int readFdata(DIRFILE* D, string filename, string boloname, string prefixe,
		fftw_complex *fdata, long ns)
/*!  read the map preconditioner from disk  */
{

	// set dirfile and bin names
	string outfile = prefixe + filename + "_" + boloname;
	const char * field_code;
	field_code = outfile.c_str();

	// fill fdata with binary
	int nget = gd_getdata(D, field_code, 0, 0, 0, ns / 2 + 1, GD_COMPLEX128,
			fdata);
	if (gd_error(D) != 0) {
		cout << "error getdata in read_fdata : read " << nget << endl;
		return 1;
	}

	gd_flush(D, field_code);

	return 0;
}

uint16_t readFramesFromDirfile(std::string tmp_dir, struct samples &samples_struct) {

	uint16_t returnCode = 0;

	DIRFILE* H;
	string scan_name, filedir;
	unsigned int nframe;

	// Need to go and open the particular dirfile to avoid cascading to the first non null RAW item...

	for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++) {

		// Open temporary dirfile to read noise size...
		scan_name = samples_struct.basevect[iframe];

		// dirfile path
		filedir = tmp_dir + "dirfile/" + scan_name + "/Data/";

		// open dirfile
		H = gd_open((char *) filedir.c_str(), GD_RDWR | GD_VERBOSE | GD_UNENCODED);

		nframe = gd_nframes(H);

		// close dirfile
		if (gd_close(H))
			cerr << " EE - Dirfile gd_close error in readFramesFrpmDirfile for : " << filedir << endl;

		if (nframe == 0)
			returnCode |= NBINS_WRONG_VALUE;

		// get nbins value
		samples_struct.nsamples[iframe] = (long) nframe;
	}

	return returnCode;
}

uint16_t readNoiseBinSizeFromDirfile(std::string tmp_dir, struct samples &samples_struct) {

	uint16_t returnCode = 0;

	DIRFILE* H;
	string scan_name, filedir;
	unsigned int nframe;

	samples_struct.nbins.assign(samples_struct.ntotscan, 0);
	samples_struct.ndet.assign(samples_struct.ntotscan, 0);

	for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++) {

		// Open temporary dirfile to read noise size...
		scan_name = samples_struct.basevect[iframe];

		// dirfile path
		filedir = tmp_dir + "dirfile/" + scan_name + "/Noise_data/Ell/";

		// open dirfile
		H = gd_open((char *) filedir.c_str(), GD_RDWR | GD_VERBOSE | GD_UNENCODED);

		nframe = gd_nframes(H);


		// close dirfile
		if (gd_close(H))
			cerr << " EE - Dirfile gd_close error in get_noise_bin_sizes for : " << filedir << endl;

		if (nframe == 0)
			returnCode |= NBINS_WRONG_VALUE;

		// get nbins value
		samples_struct.nbins[iframe] = (long) (nframe - 1);

		// get ndet value
		filedir = tmp_dir + "dirfile/" + scan_name + "/Noise_data/";
		H = gd_open((char *) filedir.c_str(), GD_RDWR | GD_VERBOSE | GD_UNENCODED);
		nframe = gd_nframes(H);

		// close dirfile
		if (gd_close(H))
			cerr << "EE - Dirfile gd_close error in get_noise_bin_sizes for : " << filedir << endl;

		if (nframe == 0)
			returnCode |= NDET_WRONG_VALUE;

		samples_struct.ndet[iframe] = ((long) nframe / samples_struct.nbins[iframe]);

	}

	return returnCode;
}


uint16_t init_dirfile(std::string tmp_dir, struct samples & samples_struct,
		int sub_rank) {

	uint16_t returnCode = 0;

	// Close previously openened dirfile
	for (long iframe = 0; iframe < samples_struct.ntotscan; iframe++) {
		if (samples_struct.dirfile_pointers[iframe]) {
			if (gd_close(samples_struct.dirfile_pointers[iframe])) {
				cerr << "EE - error closing dirfile...";
				returnCode |= TMP_PATH_PROBLEM;
			} else {
				samples_struct.dirfile_pointers[iframe] = NULL;
			}
		}
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	//TODO This should be done on iframe_min iframe_max in case of subranks and by subrank 0
	if (sub_rank == 0) {

		DIRFILE *temp;

		for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++) {

			string dirfile_basename = tmp_dir + "dirfile/" + samples_struct.basevect[iframe];

			// create folders
			samples_struct.dirfile_pointers[iframe] = gd_open(
					(char *) dirfile_basename.c_str(),
					GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE | GD_UNENCODED);

			const char * subdirs[] = { "Data", "Flag", "Lon", "Lat", "fData",
					"Indexes", "Noise_data", "Noise_data/Ell" };
			for (unsigned long ii = 0; ii < Elements_in(subdirs); ii++) {
				string dirfile_name = dirfile_basename + "/" + subdirs[ii];
				temp = gd_open((char *) dirfile_name.c_str(),
						GD_RDWR | GD_CREAT | GD_TRUNC | GD_VERBOSE
						| GD_UNENCODED);
				returnCode |= gd_close(temp);
				string dirfile_format = string(subdirs[ii]) + "/format";
				gd_include(samples_struct.dirfile_pointers[iframe],
						(char *) dirfile_format.c_str(), 0,
						GD_RDWR | GD_CREAT | GD_UNENCODED);
			}

			returnCode |= gd_flush(samples_struct.dirfile_pointers[iframe],
					NULL);
			returnCode |= gd_close(samples_struct.dirfile_pointers[iframe]);
		}

	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// and reopen it for furter use....
	for (long iframe = samples_struct.iframe_min;
			iframe < samples_struct.iframe_max; iframe++) {
		string dirfile_basename = tmp_dir + "dirfile/"
				+ samples_struct.basevect[iframe];
		samples_struct.dirfile_pointers[iframe] = gd_open(
				(char *) dirfile_basename.c_str(),
				GD_RDWR | GD_VERBOSE | GD_UNENCODED);
	}
	return returnCode;
}

uint16_t cleanup_dirfile_sanePos(std::string tmp_dir,
		struct samples & samples_struct, int sub_rank) {

	uint16_t returnCode = 0;
	std::vector<string> det_vect;

	// Close previously opened dirfile
	for (long iframe = samples_struct.iframe_min;
			iframe < samples_struct.iframe_max; iframe++) {
		if (samples_struct.dirfile_pointers[iframe]) {
			if (gd_close(samples_struct.dirfile_pointers[iframe])) {
				cerr << "EE - error closing dirfile...";
				returnCode |= TMP_PATH_PROBLEM;
			} else {
				samples_struct.dirfile_pointers[iframe] = NULL;
			}
		}
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if (sub_rank == 0) {

		for (long iframe = samples_struct.iframe_min;
				iframe < samples_struct.iframe_max; iframe++) {

			det_vect = samples_struct.bolo_list[iframe];

			string scan_name = samples_struct.basevect[iframe];
			string index_path = tmp_dir + "dirfile/" + scan_name + "/Indexes";

			DIRFILE *S = gd_open((char *) index_path.c_str(),
					GD_RDWR | GD_TRUNC | GD_VERBOSE | GD_UNENCODED);

			// then generate binaries and fill format file
			const char * prefixes[] = { "Indexes_" };
			for (unsigned long ip = 0; ip < Elements_in(prefixes); ip++) {
				for (long idet = 0; idet < (long) det_vect.size(); idet++) {
					string outfile = prefixes[ip] + scan_name + "_" + det_vect[idet];

					//configure dirfile field
					gd_entry_t E;
					E.field = (char*) outfile.c_str();
					E.field_type = GD_RAW_ENTRY;
					E.fragment_index = 0;
					E.spf = 1;
					E.data_type = GD_INT64;
					E.scalar[0] = NULL;

					// add to the dirfile
					returnCode |= gd_add(S, &E);
				}
				// flush all detector at once..
				returnCode |= gd_flush(S, NULL);

			}

			if (gd_close(S))
				cerr << "EE - Error closing " << index_path << "-> memory leaks ..."
				<< endl;
		}
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// and reopen it for furter use....
	for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++) {
		string dirfile_basename = tmp_dir + "dirfile/" + samples_struct.basevect[iframe];
		samples_struct.dirfile_pointers[iframe] = gd_open( (char *) dirfile_basename.c_str(), GD_RDWR | GD_VERBOSE | GD_UNENCODED);
	}

	return returnCode;
}

int cleanup_dirfile_saneInv(std::string tmp_dir, struct samples & samples_struct, int sub_rank) {

	uint16_t returnCode = 0;

	// Close previously opened dirfile
	for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++) {
		if (samples_struct.dirfile_pointers[iframe]) {
			if (gd_close(samples_struct.dirfile_pointers[iframe])) {
				cerr << "EE - error closing dirfile...";
				returnCode |= TMP_PATH_PROBLEM;
			} else {
				samples_struct.dirfile_pointers[iframe] = NULL;
			}
		}
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if (sub_rank == 0) {
		for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++) {

			string scan_name = samples_struct.basevect[iframe];

			string noise_path = tmp_dir + "dirfile/" + scan_name + "/Noise_data";
			string ell_path = noise_path + "/Ell";

			DIRFILE *S = gd_open((char *) noise_path.c_str(),  GD_RDWR | GD_TRUNC | GD_UNENCODED);
			DIRFILE *D = gd_open((char *) ell_path.c_str(),    GD_RDWR | GD_TRUNC | GD_UNENCODED);

			gd_add_string(S, (const char *) "WARN_Noise",
					(const char *) "The size of this dirfile is different from the main dirfile", 0);
			gd_add_string(D, (const char *) "WARN_Ell",
					(const char *) "The size of this dirfile is different from the main dirfile", 0);

			std::vector<string> det_vect = samples_struct.bolo_list[iframe];

			string outfile;
			for (long idet = 0; idet < (long) det_vect.size(); idet++) {

				// ell binary filename
				outfile = "Ell_InvNoisePS_" + scan_name + "_" + det_vect[idet];

				// configure dirfile field for ell
				gd_entry_t E;
				E.field = (char*) outfile.c_str();
				E.field_type = GD_RAW_ENTRY;
				E.fragment_index = 0;
				E.spf = 1;
				E.data_type = GD_DOUBLE;
				E.scalar[0] = NULL;

				// add to the dirfile
				returnCode |= gd_add(D, &E);

				// spectra filename
				outfile = "InvNoisePS_" + scan_name + "_" + det_vect[idet];

				// set field information for spectra
				E.field = (char*) outfile.c_str();

				// add to the dirfile
				returnCode |= gd_add(S, &E);
			}

			// flush all detector at once..
			returnCode |= gd_flush(D, NULL);
			returnCode |= gd_flush(S, NULL);

			if (gd_close(S))
				cerr << "EE - Error closing " << noise_path << "-> memory leaks ..." << endl;
			if (gd_close(D))
				cerr << "EE - Error closing " << ell_path << "-> memory leaks ..." << endl;

		}

	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// and reopen it for furter use....
	for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++) {
		string dirfile_basename = tmp_dir + "dirfile/" + samples_struct.basevect[iframe];
		samples_struct.dirfile_pointers[iframe] = gd_open( 	(char *) dirfile_basename.c_str(), GD_RDWR | GD_VERBOSE | GD_UNENCODED);
	}

	return returnCode;
}

uint16_t cleanup_dirfile_sanePic(std::string tmp_dir, struct samples & samples_struct, int sub_rank) {
	uint16_t returnCode = 0;

	// Close previously opened dirfile
	for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++) {
		if (samples_struct.dirfile_pointers[iframe]) {
			if (gd_close(samples_struct.dirfile_pointers[iframe])) {
				cerr << "EE - error closing dirfile...";
				returnCode |= TMP_PATH_PROBLEM;
			} else {
				samples_struct.dirfile_pointers[iframe] = NULL;
			}
		}
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if (sub_rank == 0) {

		for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++) {

			std::vector<string> det_vect = samples_struct.bolo_list[iframe];

			//get fourier transform dirfile names !
			string scan_name = samples_struct.basevect[iframe];
			string fdata_path = tmp_dir + "dirfile/" + scan_name + "/fData";

			// clean up the dirfiles with TRUNC option
			DIRFILE *S = gd_open((char *) fdata_path.c_str(), GD_RDWR | GD_TRUNC | GD_UNENCODED);

			if ( gd_error(S) == GD_E_OK ) {

				// then generate binaries and fill format file
				const char * prefixes[] = { "fData_", "fPs_" };
				for (unsigned long ip = 0; ip < Elements_in(prefixes); ip++)
					for (long idet = 0; idet < (long) det_vect.size(); idet++) {
						string outfile = prefixes[ip] + scan_name + "_"
								+ det_vect[idet];
						//configure dirfile field
						gd_entry_t E;
						E.field = (char*) outfile.c_str();
						E.field_type = GD_RAW_ENTRY;
						E.fragment_index = 0;
						E.spf = 1;
						E.data_type = GD_COMPLEX128;
						E.scalar[0] = NULL;

						// add to the dirfile
						returnCode |= gd_add(S, &E);
					}

				returnCode |= gd_flush(S, NULL);


				if (gd_close(S))
					cerr << "EE - Error closing " << fdata_path << "-> memory leaks ..." << endl;
			} else {
				returnCode |= TMP_PATH_PROBLEM;
			}

			//TODO: Probably the last chance to test the consistency of sizes....

			//			// check sizes in Indexes, data and flag format
			//			string indexes_path = tmp_dir + "dirfile/" + scan_name + "/Indexes";
			//			string data_path = tmp_dir + "dirfile/" + scan_name + "/data";
			//			DIRFILE *I = gd_open((char *) indexes_path.c_str(), GD_RDWR | GD_VERBOSE | GD_UNENCODED);
			//			DIRFILE *D = gd_open((char *) data_path.c_str(), GD_RDWR | GD_VERBOSE | GD_UNENCODED);
			//
			//			long nframeI = gd_nframes(I);
			//			long nframeD = gd_nframes(I);
			//
			//			//			long ns = samples_struct.nsamples[iframe];
			//			gd_close(I);
			//			gd_close(D);
			//
			//			if (nframeI != nframeD ) {
			//				cout << "Error... Dirfile data or Indexes has incorrect size !!" << endl;
			//				cout << indexes_path << " : " << nframeI << endl;
			//				cout << data_path    << " : " << nframeD << endl;
			//				return 1;
			//			}

		}
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// and reopen it for furter use....
	for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++) {
		string dirfile_basename = tmp_dir + "dirfile/" + samples_struct.basevect[iframe];
		samples_struct.dirfile_pointers[iframe] = gd_open( (char *) dirfile_basename.c_str(), GD_RDWR | GD_VERBOSE | GD_UNENCODED);
	}

	return returnCode;
}

