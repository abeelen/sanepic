/*
 * Corr_preprocess.cpp
 *
 *  Created on: 20 juil. 2009
 *      Author: matthieu
 */

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

#include "Corr_preprocess.h"
#include "covMatrixIO.h"
#include "dataIO.h"
#include "inline_IO2.h"
#include "inputFileIO.h"
#include "todprocess.h"
#include "map_making.h"
#include "struct_definition.h"
#include <time.h>


//temp
#include <fstream>
#include <iostream>
#include <sstream>

#include <gsl/gsl_math.h>
#include <fftw3.h>
//#include <time.h>

using namespace std;

#if defined(PARA_BOLO) && ! defined(USE_MPI)
#define USE_MPI
#endif

//void write_tfAS(double *S, long long *indpix, long NAXIS1, long NAXIS2, long long npix,
//		bool flgdupl, string dir, long ns, long ndet, long iframe, std::vector<string> bolonames)

void write_tfAS(double *S, struct detectors det,long long *indpix, long NAXIS1, long NAXIS2, long long npix,
		bool flgdupl, string dir, long ns, long iframe, int rank, int size)
{



	double *Ps;
	long long *samptopix;

	fftw_plan fftplan;
	fftw_complex *fdata;

	samptopix = new long long[ns];
	Ps = new double[ns];
	fdata = new fftw_complex[ns/2+1];

	fill(samptopix,samptopix+ns,0);
	fill(Ps,Ps+ns,0.0);
	for(long ii =0; ii<ns/2+1;ii++){
		fdata[ii][0]=0.0;
		fdata[ii][1]=0.0;
	}

	int factdupl = 1;
	if(flgdupl==1)  factdupl = 2;


#ifdef PARA_BOLO
	for (long idet1=rank*det.ndet/size;idet1<(rank+1)*det.ndet/size;idet1++){
#else
		for (long idet1=0;idet1<det.ndet;idet1++){
#endif



			//Read pointing data
			read_samptopix(ns, samptopix, dir, idet1, iframe, det.boloname);

			//		cout << "samptopix : " << endl;
			//		cout << samptopix[0] << " " << samptopix[1] << " " << samptopix[2] << endl;

			// temporary down
			deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,flgdupl,factdupl);

			//Fourier transform of the data
			fftplan = fftw_plan_dft_r2c_1d(ns, Ps, fdata, FFTW_ESTIMATE);
			fftw_execute(fftplan);
			fftw_destroy_plan(fftplan);


			write_fPs(ns, fdata, dir, idet1, iframe, det.boloname);

		}

		delete[] samptopix;
		delete[] Ps;
		delete[] fdata;

	}

#ifdef LARGE_MEMORY
	void write_ftrProcesdata(double *S, struct param_process proc_param, struct samples samples_struct, struct param_positions pos_param,
			string tmp_dir,	struct detectors det, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2,
			long long npix,	long long npixsrc, long long addnpix, double f_lppix, long ns, long iframe, int rank, int size, std::ofstream &file, fftw_complex *&fdatas)
#else
	void write_ftrProcesdata(double *S, struct param_process proc_param, struct samples samples_struct, struct param_positions pos_param,
			string tmp_dir,	struct detectors det, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2,
			long long npix,	long long npixsrc, long long addnpix, double f_lppix, long ns, long iframe, int rank, int size, std::ofstream &file)
#endif
	{



		double *data, *bfilter, *data_lp, *Ps;
		short *flag;
		long long *samptopix;

		fftw_plan fftplan;
		fftw_complex *fdata;


		string field1, fits_filename;


		data_lp = new double[ns];

		samptopix = new long long[ns];
		Ps = new double[ns];
		bfilter = new double[ns/2+1];
		fdata = new fftw_complex[ns/2+1];

		fill(data_lp,data_lp+ns,0.0);
		fill(Ps,Ps+ns,0.0);
		fill(bfilter,bfilter+(ns/2+1),0.0);
		fill(samptopix,samptopix+ns,0);

		for (long ii=0;ii<ns/2+1;ii++){
			fdata[ii][0] = 0.0;
			fdata[ii][1] = 0.0;
		}


		int factdupl = 1;
		if(pos_param.flgdupl==1)		factdupl = 2;

		fits_filename = samples_struct.fits_table[iframe];
		cout << "fits file : " << fits_filename << endl;


#ifdef PARA_BOLO
		for (long idet1=rank*det.ndet/size;idet1<(rank+1)*det.ndet/size;idet1++){
			cout << "[ " << rank << " ] progression write_ftr : " << 100.0*(1.0-((double)(rank+1)-(double)idet1*(double)size/(double)det.ndet)) << " %" << endl;
			std::ostringstream oss;
			oss << tmp_dir + "fdata_" << iframe << "_" << det.boloname[idet1] << ".bi";
			time_t rawtime;
			struct tm * timeinfo;
			time ( &rawtime );
			timeinfo = localtime ( &rawtime );
			file << "Writing file : " << oss.str() << " at " << asctime (timeinfo) << endl;
#else
			for (long idet1=0;idet1<det.ndet;idet1++){
				cout << "[ " << rank << " ] progression write_ftr : " << 100.0*(1.0-((double)(rank+1)-(double)idet1*(double)size/(double)det.ndet)) << " %" << endl;
#endif

				field1 = det.boloname[idet1];
				//				cout << field1 << endl;

				//				if (rank==1)
				//					cout << idet1 << endl;



				for (long ii=0;ii<ns/2+1;ii++){
					fdata[ii][0] = 0.0;
					fdata[ii][1] = 0.0;
				}
				//		cout << field1 << "  apres ALLOC "  << endl;



				long test_ns;
				read_signal_from_fits(fits_filename, field1, data, test_ns);
				if (test_ns != ns) {
					cerr << "Read signal does not correspond to frame size : Check !!" << endl;
					exit(-1);
				}

				read_flag_from_fits(fits_filename , field1, flag, test_ns);
				if (test_ns != ns) {
					cerr << "Read flag does not correspond to frame size : Check !!" << endl;
					exit(-1);
				}






				//TODO : Ps should not be here...  remove the signal before or make the deproject inside MapMakePreProcess
				//TODO : write fdata inside MapMakePreProcess.. or create a function same is true in sanePS


				if (S != NULL){
					//// Read pointing
					read_samptopix(ns, samptopix, tmp_dir, idet1, iframe, det.boloname);

					//TODO : Fix that... same number of argument... not the same calling as in sanePS
					if (addnpix){
						deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,2,factdupl,samples_struct.ntotscan,indpsrc,npixsrc);
					} else {
						deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,2,factdupl);
					}
				}


				//TODO : Ps should not be here...  remove the signal before or make the deproject inside MapMakePreProcess
				//TODO : write fdata inside MapMakePreProcess.. or create a function same is true in sanePS

				if (S != NULL){
					//********************  pre-processing of data ********************//
					MapMakPreProcessData(data,flag,ns,proc_param.napod,proc_param.poly_order,f_lppix,data_lp,bfilter, // default poly order = 4
							proc_param.NORMLIN,proc_param.NOFILLGAP,proc_param.remove_polynomia,Ps);
				}
				else {
					MapMakPreProcessData(data,flag,ns,proc_param.napod,proc_param.poly_order,f_lppix,data_lp,bfilter, // default poly order = 4
							proc_param.NORMLIN,proc_param.NOFILLGAP,proc_param.remove_polynomia);
				}

				//Fourier transform of the data
				fftplan = fftw_plan_dft_r2c_1d(ns, data_lp, fdata, FFTW_ESTIMATE);
				fftw_execute(fftplan);
				fftw_destroy_plan(fftplan);

#ifdef LARGE_MEMORY
				for(long ii=0;ii<(ns/2+1);ii++){
					fdatas[((ns/2+1)*idet1)+ii][0]=fdata[ii][0];
					fdatas[((ns/2+1)*idet1)+ii][1]=fdata[ii][1];
				}
#else
				//write fourier transform to disk
				write_fdata(ns, fdata, tmp_dir, idet1, iframe, det.boloname);
#endif


				delete [] flag;
				delete [] data;
			} // idet1



			delete[] data_lp;
			delete[] samptopix;
			delete[] Ps;
			delete[] bfilter;
			delete[] fdata;

		}

#ifdef LARGE_MEMORY
		void do_PtNd(double *PNd, string *noise_table, string dir, string prefixe,
				struct detectors det, double f_lppix, double fsamp, long ns, int rank, int size,
				long long *indpix, long NAXIS1, long NAXIS2, long long npix, long iframe,
				double *Mp, long *hits,std::ofstream &file,fftw_complex *fdatas)
#else
		void do_PtNd(double *PNd, string *noise_table, string dir, string prefixe,
				struct detectors det, double f_lppix, double fsamp, long ns, int rank, int size,
				long long *indpix, long NAXIS1, long NAXIS2, long long npix, long iframe,
				double *Mp, long *hits,std::ofstream &file)
#endif
		{

			long  nbins;
			string field1, field2;
			string extentNoiseSp;

			string nameSpfile;

			long long *samptopix;
			double *ell, *SpN, *bfilter, *bfilter_, *Nk, *Nd;

			double powered;

			fftw_plan fftplan;
			fftw_complex *fdata, *Ndf;


			samptopix = new long long[ns];
			Nd = new double[ns];
			bfilter = new double[ns/2+1];
			bfilter_ = new double[ns/2+1];
			Nk = new double[ns/2+1];
			fdata = new fftw_complex[ns/2+1];
			Ndf = new fftw_complex[ns/2+1];

			double **SpN_all;

			//TODO : This is a butterworth filter.... why not use butterworth()
			for (long ii=0;ii<ns/2+1;ii++){
				powered=gsl_pow_int(double(ii)/f_lppix,16);
				bfilter[ii] = powered /(1.0+powered);
			}
			for (long ii=0;ii<ns/2+1;ii++)
				bfilter_[ii] = 1.0/(bfilter[ii]+0.000001);

			fill(Nd,Nd+ns,0.0);
			fill(Nk,Nk+(ns/2+1),0.0);
			fill(samptopix,samptopix+ns,0);


			//			cout << "do_ptnd " << rank << endl;



#ifdef PARA_BOLO
			for (long idet1=rank*det.ndet/size;idet1<(rank+1)*det.ndet/size;idet1++){
				cout << "[ " << rank << " ] progression do_ptNd : " << 100.0*(1.0-((double)(rank+1)-(double)idet1*(double)size/(double)det.ndet)) << " %" << endl;
				std::ostringstream oss;
				oss << "frame : " << iframe << " bolo : " << det.boloname[idet1];
				time_t rawtime;
				struct tm * timeinfo;
				time ( &rawtime );
				timeinfo = localtime ( &rawtime );
				file << "do_ptnd : " << oss.str() << " at " << asctime (timeinfo) << endl;
#else
				for (long idet1=0;idet1<det.ndet;idet1++){
					cout << "[ " << rank << " ] progression do_ptNd : " << 100.0*(1.0-((double)(rank+1)-(double)idet1*(double)size/(double)det.ndet)) << " %" << endl;
#endif
					field1 = det.boloname[idet1];
					//		cout << field1 << endl;



					//Read pointing data
					read_samptopix(ns, samptopix, dir, idet1, iframe, det.boloname);


					//**************************************** Noise power spectrum
					string extname = "_InvNoisePS";
					//string suffix = Basename(noise_table[iframe]) + extname;
					string suffix = FitsBasename(noise_table[iframe]) + extname;

					//read noise PS file
					long ndet2;
					read_InvNoisePowerSpectra(dir, field1,  suffix, &nbins, &ndet2, &ell, &SpN_all);
					if(det.ndet!=ndet2) cout << "Error. The number of detector in noisePower Spectra file must be egal to input bolofile number\n";

					SpN = new double[nbins];
					fill(SpN,SpN+nbins,0.0);

					//Init N-1d
					for (long ii=0;ii<ns/2+1;ii++){
						Ndf[ii][0] = 0.0;
						Ndf[ii][1] = 0.0;
					}


					for (long idet2=0;idet2<det.ndet;idet2++){
						field2 = det.boloname[idet2];

						fill(Nd,Nd+ns,0.0);
						fill(Nk,Nk+(ns/2+1),0.0);

#ifdef LARGE_MEMORY
						for(long ii=0;ii<(ns/2+1);ii++){
							fdata[ii][0]=fdatas[((ns/2+1)*idet2)+ii][0];
							fdata[ii][1]=fdatas[((ns/2+1)*idet2)+ii][1];
						}
#else
						//read Fourier transform of the data
						read_fdata(ns, fdata, prefixe, dir, idet2, iframe, det.boloname);
#endif



						//****************** Cross power spectrum of the noise  ***************//
						for (int ii=0;ii<nbins;ii++){
							SpN[ii] = SpN_all[idet2][ii];
						}



						// TODO : Why do we need to reinterpolate the noise power spectrum here ?
						// interpolate logarithmically the noise power spectrum
						InvbinnedSpectrum2log_interpol(ell,SpN,bfilter_,nbins,ns,fsamp,Nk);
						//InvbinnedSpectrum2bis(ell,SpN,bfilter_,nbins,ns,fsamp,Nk);

						//


						for (long jj=0;jj<ns/2+1;jj++){
							if (isnan(Nk[jj])) {
								printf("isnan has been found : iframe %ld, det1 %ld, det2 %ld\n",iframe, idet1, idet2);
								exit(1);
							}
						}

						//********************************* compute N^-1 d  ***********************//
						for (long ii=0;ii<ns/2+1;ii++){
							Ndf[ii][0] += (fdata[ii][0]*Nk[ii]);
							Ndf[ii][1] += (fdata[ii][1]*Nk[ii]);
						}




						//Compute weight map for preconditioner
						if ((Mp != NULL) && (idet2 == idet1))
							compute_diagPtNPCorr(Nk,samptopix,ns,NAXIS1, NAXIS2,indpix,npix,f_lppix,Mp);
						//


					}// end of idet2 loop

					// dEBUG


					fftplan = fftw_plan_dft_c2r_1d(ns, Ndf, Nd, FFTW_ESTIMATE);
					fftw_execute(fftplan);
					fftw_destroy_plan(fftplan);

					for (long ii=0;ii<ns;ii++){
						PNd[indpix[samptopix[ii]]] += Nd[ii]; // Nd real
					}

					//compute hit counts
					if (hits != NULL){
						for (long ii=0;ii<ns;ii++){
							hits[indpix[samptopix[ii]]] += 1;
						}
					}


					delete[] ell;
					delete[] SpN;
					free_dmatrix(SpN_all,0,det.ndet-1,0,nbins-1);


				}// end of idet1 loop


				delete[] samptopix;
				delete[] Nd;
				delete[] bfilter;
				delete[] bfilter_;
				delete[] Nk;
				delete[] fdata;
				delete[] Ndf;


			}

