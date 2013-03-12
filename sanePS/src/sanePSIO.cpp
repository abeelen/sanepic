#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>

#include "sanePSIO.h"
#include "ErrorCode.h"
#include "CovMatrixIO.h"
#include "StructDefinition.h"
#include "InputFileIO.h"
#include "ErrorCode.h"


extern "C" {
#include <fitsio.h>
}

using namespace std;

int restore_session(string tmp_dir, string filename, int &goto_step, double **commonm2, double **N,
		double **P, double **Rellth, double **Rellexp, double *SPref, long ndet, int ncomp, long nbins, long ns)
{

	cout << "Restoring : " << filename << endl;
	FILE* fp;
	string file = tmp_dir + "data_saved_sanePS_" + filename + ".bi";
	size_t len=0;

	if (NULL == (fp = fopen(file.c_str(), "rb")))
	{
		cout << "Unable to open " << file << " for reading\n";
		return 0;
	}

	len+=fread(&goto_step,sizeof(int),1,fp);

	switch(goto_step)
	{
	case 2:
		// load commonm2
		for (int ii = 0; ii < ncomp; ii++)
			for(long in=0; in<ns; in++)
				len+=fread(&commonm2[ii][in], sizeof(double), 1, fp);
		break;

	case 3:
		for (int idet = 0; idet < ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				//				fprintf(fp,"%10.15g\n",N[idet][ibin]);
				len+=fread(&N[idet][ibin], sizeof(double), 1, fp);

		for (int ii = 0; ii < ncomp; ii++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&P[ii][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&Rellth[idet][ibin], sizeof(double), 1, fp);
		break;

	case 4:
	case 5:
		// load N, P, Rellexp, Rellth, SPref
		for (int idet = 0; idet < ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&N[idet][ibin], sizeof(double), 1, fp);

		for (int ii = 0; ii < ncomp; ii++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&P[ii][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&Rellth[idet][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fread(&Rellexp[idet][ibin], sizeof(double), 1, fp);

		len+=fread(SPref, sizeof(double), nbins, fp);
		break;

	case 6:
		return 0;

	default :
		cerr << "EE - completed_step has an incorrect value in : " << file << ". Exiting ...\n" << endl;
		return 1;

	}

	fclose(fp);

	return 0;
}


int save_session(string tmp_dir, string filename, int goto_step, double **commonm2, double **N,
		double **P, double **Rellth, double **Rellexp, double *SPref, long ndet, int ncomp, long nbins, long ns){

	FILE* fp;
	string file = tmp_dir + "data_saved_sanePS_" + filename + ".bi";
	size_t len=0;

	if (NULL == (fp = fopen(file.c_str(), "w+")))
	{
		cout << "Unable to open " << file << " for writing\n";
		return 1;
	}

	len+=fwrite(&goto_step,sizeof(int),1,fp);

	switch(goto_step)
	{
	case 2:
		// load commonm2
		for (int ii = 0; ii < ncomp; ii++)
			for(long in=0; in<ns; in++)
				len+=fwrite(&commonm2[ii][in], sizeof(double), 1, fp);
		break;

	case 3:
		// write N, P, Rellth
		for (int idet = 0; idet < ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&N[idet][ibin], sizeof(double), 1, fp);

		for (int ii = 0; ii < ncomp; ii++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&P[ii][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&Rellth[idet][ibin], sizeof(double), 1, fp);
		break;

	case 4:
	case 5:
		// write N, P, Rellexp, SPref
		for (int idet = 0; idet < ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&N[idet][ibin], sizeof(double), 1, fp);

		for (int ii = 0; ii < ncomp; ii++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&P[ii][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&Rellth[idet][ibin], sizeof(double), 1, fp);

		for (int idet = 0; idet < ndet*ndet; idet++)
			for (int ibin = 0; ibin < nbins; ibin++)
				len+=fwrite(&Rellexp[idet][ibin], sizeof(double), 1, fp);

		len+=fwrite(SPref, sizeof(double), nbins, fp);
		break;

	case 6:
		return 0;

	default :
		cerr << " EE - completed_step has an incorrect value in : " << file << ". Exiting ...\n" << endl;
		return 1;

	}

	fclose(fp);

	return 0;
}


int write_to_disk(string outdirSpN, string fits_filename, struct param_sanePS structPS, std::vector<std::string> det,  long nbins, double *ell, double **mixmat,
		double **Rellth, double **Rellexp, double **N, double *SPref, double **P)
{


	ostringstream temp_stream;
	string testfile;
	string tempstr1, tempstr2;
	string nameSpfile;
	string basename;

	long ndet = (long)det.size();

	FILE *fp;
	double *data1d;

	basename = FitsBasename(fits_filename);

	temp_stream << outdirSpN << basename << structPS.ell_suffix;
	nameSpfile = temp_stream.str();
	temp_stream.str("");

	fp = fopen(nameSpfile.c_str(),"w");
	for (long ii=0;ii<nbins;ii++)
		fprintf(fp,"%g\n",ell[ii]);
	fclose(fp);

	temp_stream << "!" << outdirSpN << basename << structPS.cov_matrix_suffix;
	nameSpfile= temp_stream.str();
	temp_stream.str("");

	if(write_CovMatrix(nameSpfile, det, nbins, ell, Rellth))
		return FILE_PROBLEM;


	temp_stream << "!" << outdirSpN << basename << "_exp" << structPS.cov_matrix_suffix;
	nameSpfile= temp_stream.str();
	temp_stream.str("");

	if(write_CovMatrix(nameSpfile, det, nbins, ell, Rellexp))
		return FILE_PROBLEM;


#ifdef DEBUG

	temp_stream << outdirSpN << basename << "_exp.psd";
	nameSpfile= temp_stream.str();
	temp_stream.str("");
	fp = fopen(nameSpfile.c_str(),"w");

	for (long idet1=0;idet1<ndet;idet1++){
		for (long idet2=0;idet2<ndet;idet2++){

			///// write power spectrum to disk
			tempstr1 = det[idet1];
			tempstr2 = det[idet2];
			fprintf(fp,"%s%s%s\n",tempstr1.c_str(),"-",tempstr2.c_str());
			fprintf(fp,"%d\n",(int)nbins);
			for (long ii=0;ii<nbins;ii++){
				fprintf(fp,"%g\t",ell[ii]);
				fprintf(fp,"%10.15g\n",(Rellexp[idet1*(ndet)+idet2][ii]+Rellexp[idet2*(ndet)+idet1][ii])/2.0*SPref[ii]);
			}
			fprintf(fp,"%g\n",ell[nbins]);
		}
	}
	fclose(fp);

#endif

	temp_stream << outdirSpN  << basename << structPS.mix_suffix;
	testfile= temp_stream.str();
	temp_stream.str("");

	fp = fopen(testfile.c_str(),"w");
	fprintf(fp,"%i\n",structPS.ncomp);
	for (long idet=0;idet<ndet;idet++){
		for (long iComp=0;iComp<structPS.ncomp;iComp++){
			fprintf(fp,"%10.15g ",mixmat[idet][iComp]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);


#ifdef DEBUG

	//**************** Write component power spectra to disk
	for (long idet1=0;idet1<ndet;idet1++){

		tempstr1 = det[idet1];
		temp_stream << outdirSpN + tempstr1 + "_uncnoise_" << basename << ".psd";
		nameSpfile= temp_stream.str();
		temp_stream.str("");

		fp = fopen(nameSpfile.c_str(),"w");
		for (long ii=0;ii<nbins;ii++){
			fprintf(fp,"%10.15g\n",N[idet1][ii]*SPref[ii]);
		}
		fprintf(fp,"\n");
		fclose(fp);
	}

	data1d = new double[ndet*nbins];
	for (long i=0; i< ndet; i++)
		for (long j=0; j<nbins; j++)
			data1d[i*nbins+j] = N[i][j]*SPref[j];

	temp_stream << "!" + outdirSpN + "Nfinal_" << basename << "_uncnoise.fits";
	testfile= temp_stream.str();
	temp_stream.str("");
	write_psd_tofits(testfile,nbins,ndet,'d',data1d);
	delete [] data1d;

	for (long jj=0;jj<structPS.ncomp;jj++){

		temp_stream << outdirSpN + "Comp_" << jj << "_uncnoise_" << basename << ".psd";
		nameSpfile= temp_stream.str();
		temp_stream.str("");

		fp = fopen(nameSpfile.c_str(),"w");
		for (long ii=0;ii<nbins;ii++){
			fprintf(fp,"%10.15g\n",P[jj][ii]*SPref[ii]);
		}
		fprintf(fp,"\n");
		fclose(fp);
	}

#endif


	data1d = new double[structPS.ncomp*nbins];
	for (long i=0; i< structPS.ncomp; i++)
		for (long j=0; j<nbins; j++)
			data1d[i*nbins+j] = P[i][j]*SPref[j];

	temp_stream << "!" + outdirSpN + "Nfinal_" << basename << "_cnoise.fits";

	// get filename
	testfile= temp_stream.str();
	temp_stream.str("");
	write_psd_tofits(testfile,nbins,structPS.ncomp,'d',data1d);



	delete [] data1d;

	return EXIT_SUCCESS;

}


int write_psd_tofits(string fname, long nx, long ny,
		char dtype, void * psd1d) {

	fitsfile *fp;
	int fits_status = 0;

	long naxis = 2;           // number of dimensions
	long naxes[] = {nx, ny};  // size of dimensions
	long fpixel[] = {1, 1};   // index for write_pix
	long ndata = nx * ny;     // number of data points

	// create fits file
	if ( fits_create_file(&fp, fname.c_str(), &fits_status) ){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	// create fits image (switch on data type)
	switch (dtype) {
	case 'd':    // double
		if ( fits_create_img(fp, DOUBLE_IMG, naxis, naxes, &fits_status) )
		{
			fits_report_error(stderr, fits_status);
			return 1;
		}
		if ( fits_write_pix(fp, TDOUBLE, fpixel, ndata, (double*) psd1d, &fits_status) ){
			fits_report_error(stderr, fits_status);
			return 1;
		}
		break;
	case 'l':    // long
		if ( fits_create_img(fp, LONG_IMG, naxis, naxes, &fits_status) ){
			fits_report_error(stderr, fits_status);
			return 1;
		}
		if ( fits_write_pix(fp, TLONG, fpixel, ndata, (long*) psd1d, &fits_status) ){
			fits_report_error(stderr, fits_status);
			return 1;
		}
		break;
	default:
		cerr << "write_fits: data type '" << dtype << "' not supported. Exiting.\n";
		return 1;
	}

	// close file
	if (fits_close_file(fp, &fits_status)){
		fits_report_error(stderr, fits_status);
		return 1;
	}

	return 0;

}


