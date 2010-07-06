/*
 * conjugate_gradient.cpp
 *
 *  Created on: 3 juil. 2009
 *      Author: matthieu
 */


#include "conjugate_gradient.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

#include "imageIO.h"
#include "Corr_preprocess.h"
#include "NoCorr_preprocess.h"


#if defined(USE_MPI) || defined(PARA_BOLO)
#include "mpi.h"
#endif

using namespace std;


void sanepic_conjugate_gradient(struct samples samples_struct,  struct param_positions pos_param,
		struct detectors det, struct common dir, struct param_process proc_param, long long npix, double* &S,
		long iframe_min, long iframe_max, std::vector<double> fcut, long long *indpix, struct wcsprm * wcs,
		long NAXIS1, long NAXIS2, int iterw, long long *indpsrc, long long npixsrc, int flagon, int rank, int size,
		double *&PNdtot, long long addnpix)
{


	FILE *fp;
	//bool fru;
	string fname, testfile;
	ostringstream temp_stream;
	double *PtNPmatS,  *PtNPmatStot=NULL, *r, *q, *qtot=NULL, *d, *Mp, *Mptot=NULL, *s; // =NULL to avoid warnings
	// Mp = M in the paper = preconditioner
	long *hits, *hitstot=NULL;
	double *PNd;

	double var0 = 0.0, var_n = 0.0, delta0 = 0.0, delta_n = 0.0, alpha = 0.0;
	double delta_o, rtq, beta;

	long mi;
	double *map1d;

	int iter;
	double  f_lppix_Nk,f_lppix;
	long ns;
	long long npixeff;

	int factdupl = 1;
	if(pos_param.flgdupl==1) factdupl = 2;

	// memory allocs
	r           = new double[npix];
	q           = new double[npix];
	d           = new double[npix];
	Mp          = new double[npix];
	s           = new double[npix];
	PtNPmatS    = new double[npix];
	hits        = new long[npix];
	PNd         = new double[npix];

	map1d       = new double[NAXIS1*NAXIS2];

	for (int idupl = 0;idupl<=pos_param.flgdupl;idupl++){



		//Conjugate gradien Inversion
		if (pos_param.projgaps || !flagon){
			npixeff = npix;
		} else {
			npixeff = npix-1;
		}



		printf("[%2.2i] npix = %lld, npixeff = %lld\n", rank, npix, npixeff);


		//t1 = time(0);

		fill(PtNPmatS,PtNPmatS+npix,0.0);
		fill(Mp,Mp+npix,0.0);
		fill(hits,hits+npix,0);
		fill(r,r+npix,0.0); // useless
		fill(d,d+npix,0.0); // useless
		fill(s,s+npix,0.0);
//		fill(PNd,PNd+npix,0.0);



		for (long iframe=iframe_min;iframe<iframe_max;iframe++){
			ns = samples_struct.nsamples[iframe];
			f_lppix_Nk = fcut[iframe]*double(ns)/proc_param.fsamp;


			//    cout << "[" << rank << "] " << iframe << "/" << iframe_max << endl;


			// preconditioner computation : Mp
			if (proc_param.CORRon){
				write_tfAS(S,det,indpix,NAXIS1, NAXIS2,npix,pos_param.flgdupl, dir.tmp_dir,ns,iframe, rank, size);
				// read pointing + deproject + fourier transform

				//				do_PtNd(PtNPmatS,extentnoiseSp_all,noiseSppreffile,tmp_dir,prefixe,bolonames, // todo : replace all do_PtNd par la new version avec les struct
				//						f_lppix_Nk,fsamp,ns,ndet,/*size_det,rank_det,*/indpix,NAXIS1, NAXIS2,npix,iframe,Mp,hits);

#ifdef PARA_BOLO
				//				cout << "rank " << rank << " a fini et attend ! \n";
				MPI_Barrier(MPI_COMM_WORLD);
#endif

				do_PtNd(PtNPmatS, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
						proc_param.fsamp,ns, rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,Mp,hits);
				//				do_PtNd(PNd, samples_struct.noise_table,dir.tmp_dir,prefixe,det,f_lppix_Nk,
				//								u_opt.fsamp,ns/*,size_det,rank_det*/,indpix,NAXIS1, NAXIS2,npix,iframe,Mp,hits/*,fdata_buffer*/);
				// return Pnd = At N-1 d
			} else {

				//				do_PtNPS_nocorr(S,extentnoiseSp_all,noiseSppreffile,tmp_dir,dirfile,bolonames,
				//						f_lppix_Nk,fsamp,flgdupl,factdupl,ns,ndet,indpix,
				//						NAXIS1, NAXIS2,npix,iframe,PtNPmatS,Mp,hits);

				do_PtNPS_nocorr(S, samples_struct.noise_table, dir, det,f_lppix_Nk,
						proc_param.fsamp, pos_param.flgdupl, ns, indpix, NAXIS1, NAXIS2, npix,
						iframe, PtNPmatS, Mp, hits);
			}

		} // end of iframe loop

		//cout << "do_ptnd" << endl;


#if defined(USE_MPI) || defined(PARA_BOLO)

		if(rank==0){
			PtNPmatStot = new double[npix];
			hitstot=new long[npix];
			Mptot = new double[npix];
			qtot = new double[npix];

			// todo : Is it really necessary to initialize those variables ?
			fill(PtNPmatStot,PtNPmatStot+npix,0.0);
			fill(hitstot,hitstot+npix,0);
			fill(Mptot,Mptot+npix,0.0);
			fill(qtot,qtot+npix,0.0);
		}

		MPI_Reduce(PtNPmatS,PtNPmatStot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(hits,hitstot,npix,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(Mp,Mptot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else

		PtNPmatStot=PtNPmatS; // ajout Mat 02/06
		hitstot=hits;
		Mptot=Mp;

#endif

		// intitialisation of the Conjugate gradient with preconditioner
		// see : http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
		if (rank == 0) {

			for (long ii=0;ii<npixeff;ii++)
				if (Mptot[ii] == 0)
					printf("ERROR: Mp[%ld] has elements = 0\n",ii);


			for (long ii=0;ii<npixeff;ii++)
				Mptot[ii] = 1.0/Mptot[ii]; // M : preconditioner


			for (long ii=0;ii<npixeff;ii++)
				r[ii] = PNdtot[ii] - PtNPmatStot[ii]; // r = b - Ax

			//			cout << r[0] << " " << r[1] << " " << r[2] << " " << endl;

			for (long ii=0;ii<npixeff;ii++)
				d[ii] =  Mptot[ii] * r[ii]; // d = M-1 * r


			delta_n = 0.0;
			for (long ii=0;ii<npixeff;ii++)
				delta_n += r[ii]*d[ii]; // delta_new = rT * d

			//			cout << "delta_n : " << delta_n << endl;

			var_n = 0.0;
			for (long ii=0;ii<npixeff;ii++)
				var_n += r[ii]*r[ii];


			delta0 = delta_n; // delta_0 <= delta_new
			var0 = var_n;
			printf("[%2.2i] var0 = %lf\n",rank, var0);

		}

#if defined(USE_MPI) || defined(PARA_BOLO)
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&var0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&var_n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(d,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

		printf("[%2.2i] Main Conjugate gradient loop started\n",rank);


		//start loop
		iter = 0; // max iter = 2000, but ~100 iterations are required to achieve convergence

		// while i<imax and var_new > epsilon² * var_0 : epsilon² = 1e-10 => epsilon = 1e-5
		while(((iter < 2000) && (var_n/var0 > 1e-10) && (idupl || !pos_param.flgdupl)) || (!idupl && var_n/var0 > 1e-4)){ // 2000
			// added brackets in order to avoid warning, mat-27/05
			//			if(iter==2)
			//				break;
			fill(q,q+npixeff,0.0); // q <= A*d
			//fill(qtot,qtot+npixeff,0.0);



			for (long iframe=iframe_min;iframe<iframe_max;iframe++){
				ns = samples_struct.nsamples[iframe];
				//				ff = fframes[iframe];
				f_lppix_Nk = fcut[iframe]*double(ns)/proc_param.fsamp;

				if (proc_param.CORRon){
					write_tfAS(d,det,indpix,NAXIS1, NAXIS2,npix,pos_param.flgdupl, dir.tmp_dir,ns,iframe, rank, size);
					// read pointing + deproject + fourier transform

					//					do_PtNd(q,extentnoiseSp_all,noiseSppreffile,tmp_dir,prefixe,bolonames,f_lppix_Nk,
					//							fsamp,ns,ndet,/*size_det,rank_det,*/indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL);

#ifdef PARA_BOLO
					//					cout << "rank " << rank << " a fini et attend ! \n";
					MPI_Barrier(MPI_COMM_WORLD);
#endif
					do_PtNd(q, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
							proc_param.fsamp,ns, rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL);
					// return Pnd = At N-1 d
				} else {

					//					do_PtNPS_nocorr(d,extentnoiseSp_all,noiseSppreffile,tmp_dir,dirfile,bolonames,
					//							f_lppix_Nk,fsamp,flgdupl,factdupl,ns,ndet,/*size_det,rank_det,*/indpix,
					//							NAXIS1, NAXIS2,npix,iframe,q,NULL,NULL);

					do_PtNPS_nocorr(d, samples_struct.noise_table, dir, det,f_lppix_Nk,
							proc_param.fsamp, pos_param.flgdupl, ns, indpix, NAXIS1, NAXIS2, npix,
							iframe, q, NULL, NULL);
				}
			} // end of iframe loop




#if defined(USE_MPI) || defined(PARA_BOLO)
			//cout << " q reduction\n";
			MPI_Reduce(q,qtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
			qtot=q;
#endif



			if (rank == 0){
				rtq= 0.0;
				for (long ii=0;ii<npixeff;ii++)
					rtq += qtot[ii] * d[ii]; // rtq = (dT * q)

				alpha = delta_n/rtq; // alpha <= delta_new / (dT * q)
				//				cout << "rtq : " << rtq << endl;
				//				cout << "alpha : " << alpha << endl;


				for (long ii=0;ii<npixeff;ii++)
					S[ii] += alpha*d[ii]; // x = x + alpha * d, x = S = signal
			}

#if defined(USE_MPI) || defined(PARA_BOLO)
			MPI_Barrier(MPI_COMM_WORLD);
			//cout << rank << " S bcast\n";
			MPI_Bcast(S ,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif


			// every 10 iterations do ....
			if ((iter % 10) == 0){ // if iter is divisible by 10, recompute PtNPmatStot

				fill(PtNPmatS,PtNPmatS+npixeff,0.0);
				//fill(PtNPmatStot,PtNPmatStot+npixeff,0.0);


				for (long iframe=iframe_min;iframe<iframe_max;iframe++){

					ns = samples_struct.nsamples[iframe];
					//					ff = fframes[iframe];
					f_lppix_Nk = fcut[iframe]*double(ns)/proc_param.fsamp;

					if (proc_param.CORRon){
						write_tfAS(S,det,indpix,NAXIS1, NAXIS2,npix,pos_param.flgdupl, dir.tmp_dir,ns,iframe, rank, size);
						// read pointing + deproject + fourier transform

						//						do_PtNd(PtNPmatS,extentnoiseSp_all,noiseSppreffile,tmp_dir,prefixe,bolonames,
						//								f_lppix_Nk,fsamp,ns,ndet,/*size_det,rank_det,*/indpix,NAXIS1, NAXIS2,npix,iframe,
						//								NULL,NULL);

#ifdef PARA_BOLO
						//						cout << "rank " << rank << " a fini et attend ! \n";
						MPI_Barrier(MPI_COMM_WORLD);
#endif
						do_PtNd(PtNPmatS, samples_struct.noise_table,dir.tmp_dir,"fPs_",det,f_lppix_Nk,
								proc_param.fsamp,ns,rank,size,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL);
						// return Pnd = At N-1 d
					} else {

						//						do_PtNPS_nocorr(S,extentnoiseSp_all,noiseSppreffile,tmp_dir,dirfile,bolonames,
						//								f_lppix_Nk,fsamp,flgdupl,factdupl,ns,ndet,/*size_det,rank_det,*/
						//								indpix,NAXIS1, NAXIS2,npix,iframe,PtNPmatS,NULL,NULL);

						do_PtNPS_nocorr(S, samples_struct.noise_table, dir, det,f_lppix_Nk,
								proc_param.fsamp, pos_param.flgdupl, ns, indpix, NAXIS1, NAXIS2, npix,
								iframe, PtNPmatS, NULL, NULL);
					}
				} // end of iframe loop


#if defined(USE_MPI) || defined(PARA_BOLO)
				//				cout << rank << " PtNPmatS reduction\n";
				MPI_Reduce(PtNPmatS,PtNPmatStot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
				PtNPmatStot=PtNPmatS;
#endif


				if (rank == 0){
					for (long ii=0;ii<npixeff;ii++)
						r[ii] = PNdtot[ii] - PtNPmatStot[ii]; //r = b - Ax
				}


			} else {

				if (rank == 0){
					for (long ii=0;ii<npixeff;ii++)
						r[ii] -= alpha*qtot[ii]; // else r = r - alpha * q
				}
			}

			if (rank == 0){

				for (long ii=0;ii<npixeff;ii++)
					s[ii] = Mptot[ii]*r[ii]; // s = M-1 * r


				delta_o = delta_n; // delta_0 <= delta_new

				delta_n = 0.0;
				for (long ii=0;ii<npixeff;ii++)
					delta_n += r[ii]*s[ii]; // delta_new = rT * s

				var_n = 0.0;
				for (long ii=0;ii<npixeff;ii++)
					var_n += r[ii]*r[ii];



				beta = delta_n/delta_o; // beta = delta_new / delta_0
				for (long ii=0;ii<npixeff;ii++)
					d[ii] = s[ii] + beta*d[ii]; // d = s + beta * d


				cout << "iter = " << iter;
				cout << ", crit  = " << setiosflags(ios::scientific) << setiosflags(ios::floatfield) << var_n/var0;
				cout << ", crit2 = " << setiosflags(ios::scientific) << setiosflags(ios::floatfield) << delta_n/delta0;
				cout << "\r " << flush;

				/*sprintf(testfile,"%s%s%s%s",outdir.c_str(),"Iteration_file_",termin.c_str(),".txt");
					      fp = fopen(testfile,"w");
					      if(fp!=NULL){
						for(ii=0;ii<npix;ii++)
						  fprintf(fp,"[%2.2i] iter = %d, crit = %10.15g, crit2 = %10.15g  \n",rank, iter,var_n/var0,delta_n/delta0);
						fclose(fp);}*/

				//      printf("[%2.2i] iter = %d, crit = %10.15g, crit2 = %10.15g     \n",rank, iter,var_n/var0,delta_n/delta0);


				//todo : Do we really need that ? Should be covered by sanePre
				if (iter == 0){
					for (long ii=0; ii<NAXIS1; ii++) {
						for (long jj=0; jj<NAXIS2; jj++) {
							mi = jj*NAXIS1 + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = Mptot[indpix[mi]];
							} else {
								map1d[mi] = NAN;
							}
						}
					}


					fname = '!' + dir.output_dir + "binMap_noisevar.fits"; // write preconditioner
					//					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
					write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);


					for (long ii=0; ii<NAXIS1; ii++) {
						for (long jj=0; jj<NAXIS2; jj++) {
							mi = jj*NAXIS1 + ii;
							map1d[mi] = 0.0;
						}
					}


					//todo : Treat the 2 pixels properly
					if (addnpix){
						for (long iframe = 0;iframe<samples_struct.ntotscan;iframe++){
							for (long ii=0; ii<NAXIS1; ii++) {
								for (long jj=0; jj<NAXIS2; jj++) {
									mi = jj*NAXIS1 + ii;
									long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
									if ((indpsrc[mi] != -1) && (indpix[ll] != -1))
										map1d[mi] += 1.0/Mptot[indpix[ll]];
								}
							}
						}

						fname = '!' + dir.output_dir + "optimMap_" + "_invnoisevaruncpix.fits";
						//						write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
						write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);
					}

				} // end of if (iter == 0)

				if (iterw && (iter % iterw) == 0){

					// make the map
					for (long ii=0; ii<NAXIS1; ii++) {
						for (long jj=0; jj<NAXIS2; jj++) {
							mi = jj*NAXIS1 + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = S[indpix[mi]];
							} else {
								map1d[mi] = NAN;
							}
						}
					}

					temp_stream << "!" + dir.output_dir + "optimMap_" + "flux" << iter << "b.fits";


					// Transform into string
					fname= temp_stream.str();
					// Clear ostringstream buffer
					temp_stream.str("");
					//					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
					write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);

					if (pos_param.flgdupl){
						for (long ii=0; ii<NAXIS1; ii++) {
							for (long jj=0; jj<NAXIS2; jj++) {
								mi = jj*NAXIS1 + ii;
								if (indpix[mi] >= 0){
									map1d[mi] = S[indpix[mi+NAXIS1*NAXIS2]]; //-finalmap[ii][jj];
								} else {
									map1d[mi] = 0.0;
								}
							}
						}


						temp_stream << "!" + dir.output_dir + "optimMap_fluxflags_" << iter << "b.fits";

						// Transform into string
						fname= temp_stream.str();
						// Clear ostringstream buffer
						temp_stream.str("");
						//						write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
						write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);
					}



					if (addnpix){
						// initialize the container
						for (long jj=0; jj<NAXIS2 ; jj++){
							for (long ii=0; ii<NAXIS1 ; ii++){
								mi = jj*NAXIS1 + ii;
								map1d[mi] = 0.0;
							}
						}
						// loop thru frame to coadd all pixels
						for (long iframe = 0;iframe<samples_struct.ntotscan;iframe++){
							for (long ii=0; ii<NAXIS1; ii++) {
								for (long jj=0; jj<NAXIS2; jj++) {
									mi = jj*NAXIS1 + ii;
									long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
									if ((indpsrc[mi] != -1) && (indpix[ll] != -1))
										map1d[mi] += S[indpix[ll]]/Mptot[mi];;
								}
							}
						}


						temp_stream << "!" + dir.output_dir + "optimMap_fluxuncpix_" << iter << "b.fits";

						// Transform into string
						fname= temp_stream.str();
						// Clear ostringstream buffer
						temp_stream.str("");
						//						write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
						write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);

					}
				} // end of if (iterw && (iter % iterw) == 0)


				temp_stream << dir.output_dir + "ConvFile.txt";

				// Transform into string
				testfile= temp_stream.str();
				// Clear ostringstream buffer
				temp_stream.str("");
				fp = fopen(testfile.c_str(),"a");
				fprintf(fp,"iter = %d, crit = %10.15g, crit2 = %10.15g\n",iter,var_n/var0, delta_n/delta0);
				fclose(fp);

			} // end of if (rank == 0)


#if defined(USE_MPI) || defined(PARA_BOLO)
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&var_n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(d ,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

			iter++; // i = i +1

		} // end of while loop
		printf("\n");








		if  ((pos_param.projgaps || (pos_param.flgdupl)) && !idupl){


			fill(PNd,PNd+npix,0.0);

			for (long iframe=iframe_min;iframe<iframe_max;iframe++){

				ns = samples_struct.nsamples[iframe];
				//				ff = fframes[iframe];
				f_lppix = proc_param.f_lp*double(ns)/proc_param.fsamp;
				f_lppix_Nk = fcut[iframe]*double(ns)/proc_param.fsamp;

				if (proc_param.CORRon){

					//					write_ftrProcesdata(S,indpix,indpsrc,NAXIS1, NAXIS2,npix,npixsrc,ntotscan,addnpix,flgdupl,factdupl,2,
					//							tmp_dir,dirfile,bolonames,fits_table,f_lppix,ns,
					//							napod,ndet,NORMLIN,NOFILLGAP,remove_polynomia,iframe);

					write_ftrProcesdata(S,proc_param,samples_struct,pos_param,dir.tmp_dir,det,indpix,indpsrc,NAXIS1, NAXIS2,npix,
							npixsrc,addnpix,f_lppix,ns,	iframe, rank, size);
					// fillgaps + butterworth filter + fourier transform
					// "fdata_" files generation (fourier transform of the data)

					//					do_PtNd(PNd,extentnoiseSp_all,noiseSppreffile,tmp_dir,prefixe,bolonames,f_lppix_Nk,
					//							fsamp,ns,ndet,/*size_det,rank_det,*/indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL);
#ifdef PARA_BOLO
					//cout << "rank " << rank << " a fini et attend ! \n";
					MPI_Barrier(MPI_COMM_WORLD);
#endif

					do_PtNd(PNd, samples_struct.noise_table,dir.tmp_dir,"fdata_",det,f_lppix_Nk,
							proc_param.fsamp,ns, rank, size,indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL);

					// return Pnd = At N-1 d
				} else {

					do_PtNd_nocorr(PNd,dir.tmp_dir,proc_param, pos_param, samples_struct,det, f_lppix, f_lppix_Nk,
							addnpix, ns,indpix, indpsrc, NAXIS1, NAXIS2, npix, npixsrc, iframe, S);
				}
			} // end of iframe loop




#if defined(USE_MPI) || defined(PARA_BOLO)
			MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
			PNdtot=PNd; // ajout Mat 02/07
#endif
		}



	}// end of idupl loop



#if defined(USE_MPI) || defined(PARA_BOLO)
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	temp_stream << dir.output_dir + "testfile.txt";

	// Transform into string
	testfile= temp_stream.str();
	// Clear ostringstream buffer
	temp_stream.str("");
	fp = fopen(testfile.c_str(),"a");
	fprintf(fp,"test avant ecriture \n");
	fclose(fp);


	//todo : Should be outside of this function
	// write map function : NAXIS1, NAXIS2, S, indpix, outdir, termin, addnpix, ntotscan,
	// pixdeg, tancoord, tanpix, coordsyst, Mptot, indpsrc, npixsrc, factdupl,
	//
	//******************************  write final map in file ********************************



	if (rank == 0){

		printf(" after CC INVERSION %lld\n",npix*(npix+1)/2);


		//todo : writing the maps should NOT be here...
		//todo : In general separate the reading/writing from the computationfor (long ii=0; ii<NAXIS1; ii++) {
		for (long jj=0; jj<NAXIS2; jj++) {
			mi = jj*NAXIS1 + ii;
			if (indpix[mi] >= 0){
				map1d[mi] = S[indpix[mi]];
			} else {
				map1d[mi] = NAN;
			}
		}
	}

	fname = '!' + dir.output_dir + "optimMap_flux.fits";
	write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);


	for (long ii=0; ii<NAXIS1; ii++) {
		for (long jj=0; jj<NAXIS2; jj++) {
			mi = jj*NAXIS1 + ii;
			if (indpix[mi] >= 0){
				map1d[mi] = Mptot[indpix[mi]];
			} else {
				map1d[mi] = NAN;
			}
		}
	}


	fname = '!' + dir.output_dir + "optimMap_noisevar.fits";
	//		write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
	write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);

	if (addnpix){
		for (long iframe = 0;iframe<samples_struct.ntotscan;iframe++){
			for (long ii=0; ii<NAXIS1; ii++) {
				for (long jj=0; jj<NAXIS2; jj++) {
					mi = jj*NAXIS1 + ii;
					long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
					if ((indpsrc[mi] != -1)  && (indpix[ll] != -1)){
						map1d[mi] = S[ll];
					} else {
						map1d[mi] = 0.0;
					}
				}
			}



			temp_stream << "!" + dir.output_dir + "optimMap_flux_fr" << iframe << ".fits";
			// Transform into string
			fname= temp_stream.str();
			// Clear ostringstream buffer
			temp_stream.str("");
			write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);

			for (long ii=0; ii<NAXIS1; ii++) {
				for (long jj=0; jj<NAXIS2; jj++) {
					mi = jj*NAXIS1 + ii;
					long long ll = factdupl*NAXIS1*NAXIS2 + iframe*npixsrc + indpsrc[mi];
					if ((indpsrc[mi] != -1)  && (indpix[ll] != -1)){
						map1d[mi] = Mptot[indpix[ll]];
					} else {
						map1d[mi] = 0.0;
					}
				}
			}

			//fname = '!' + outdir + "optimMap_" + termin + "_noisevar_fr" + iframestr + ".fits";
			temp_stream << "!" + dir.output_dir + "optimMap_noisevar_fr" << iframe << ".fits";

			// Transform into string
			fname= temp_stream.str();
			// Clear ostringstream buffer
			temp_stream.str("");
			//					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
			write_fits_wcs(fname, wcs, NAXIS1, NAXIS2, 'd', (void *)map1d);
		}
	}
}

// end of write map function

#if defined(USE_MPI) || defined(PARA_BOLO)
delete [] qtot;
delete [] Mptot;
delete [] PtNPmatStot;
delete [] hitstot;
#endif

delete [] r;
delete [] q;

delete [] d;
delete [] Mp;

delete [] s;
delete [] PtNPmatS;

delete [] hits;

delete [] map1d;
delete [] PNd;

}
