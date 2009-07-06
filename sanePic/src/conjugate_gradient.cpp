/*
 * conjugate_gradient.cpp
 *
 *  Created on: 3 juil. 2009
 *      Author: matthieu
 */


#include "conjugate_gradient.h"


void sanepic_conjugate_gradient(bool flgdupl, int npix, double* &S,long iframe_min, long iframe_max,
		long *nsamples, long *fframes,std::vector<double> fcut,double fsamp,
		long *indpix, int nn, int factdupl, string poutdir, string termin, long ndet,
		string *extentnoiseSp_all,string noiseSppreffile, std::vector<string> bolonames, int size_det,
		int rank_det, int iterw, double pixdeg, double *tancoord, double *tanpix,int coordsyst,
		long *indpsrc, long npixsrc, int flagon, bool projgaps, int rank, bool CORRon,
		string dirfile, double *&PNdtot,double *&PNd, long ntotscan,long addnpix,bool NORMLIN,bool NOFILLGAP,
		long napod,int shift_data_to_point,bool remove_polynomia,string fextension,string bextension,
		string flpoint_field,string scerr_field, string outdir){



	FILE *fp;
	char testfile[100];


	char iterchar[30];
	char iframechar[30];
	string iterstr, iframestr;
	bool fru;
	string fname;
	double *PtNPmatS,  *PtNPmatStot, *r, *q, *qtot, *d, *Mp, *Mptot, *s;
	long *hits, *hitstot;

	double var0, var_n, delta0, delta_n, delta_o, rtq, alpha, beta;

	long mi;
	double *map1d;
	string prefixe;
	int iter;
	double  f_lppix_Nk,f_lp,f_lppix;
	long ns,ff;
	int npixeff;

	double errarcsec = 15.0; // rejection criteria : scerr[ii] > errarcsec, sample is rejected
	// source error

	//time t1,t2;

	// memory allocs
	r = new double[npix];
	q = new double[npix];
	qtot = new double[npix];
	d = new double[npix];
	Mp = new double[npix];
	Mptot = new double[npix];
	s = new double[npix];
	PtNPmatS = new double[npix];
	PtNPmatStot = new double[npix];
	hits = new long[npix];
	hitstot = new long[npix];

	map1d = new double[nn*nn];

	for (int idupl = 0;idupl<=flgdupl;idupl++){

		for (long ii=0;ii<npix;ii++) S[ii] = 0.0;//PNd[ii];

		//Conjugate gradien Inversion
		if (projgaps || !flagon){
			npixeff = npix;
		} else {
			npixeff = npix-1;
		}



		printf("[%2.2i] npix = %d, npixeff = %d\n", rank, npix, npixeff);


		//t1 = time(0);

		fill(PtNPmatS,PtNPmatS+npix,0.0);
		fill(PtNPmatStot,PtNPmatStot+npix,0.0);
		fill(Mp,Mp+npix,0.0);
		fill(Mptot,Mptot+npix,0.0);
		fill(hits,hits+npix,0);
		fill(hitstot,hitstot+npix,0);


		prefixe = "fPs";

		for (long iframe=iframe_min;iframe<iframe_max;iframe++){
			ns = nsamples[iframe];
			ff = fframes[iframe];
			f_lppix_Nk = fcut[iframe]*double(ns)/fsamp;


			//    cout << "[" << rank << "] " << iframe << "/" << iframe_max << endl;

			if (CORRon){
				write_tfAS(S,indpix,nn,npix,flgdupl,factdupl, poutdir,termin,ff,ns,ndet,iframe);
				// read pointing + deproject + fourier transform

				do_PtNd(PtNPmatS,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,
						f_lppix_Nk,fsamp,ff,ns,ndet,size_det,rank_det,indpix,nn,npix,iframe,Mp,hits);
				// return Pnd = At N-1 d
			} else {

				do_PtNPS_nocorr(S,extentnoiseSp_all,noiseSppreffile,poutdir,termin,dirfile,bolonames,
						f_lppix_Nk,fsamp,flgdupl,factdupl,ff,ns,ndet,size_det,rank_det,indpix,
						nn,npix,iframe,PtNPmatS,Mp,hits);
			}

		} // end of iframe loop




#ifdef USE_MPI
		MPI_Reduce(PtNPmatS,PtNPmatStot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(hits,hitstot,npix,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(Mp,Mptot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
		for(long ii=0;ii<npix;ii++){
			PtNPmatStot[ii]=PtNPmatS[ii]; // ajout Mat 02/06
			hitstot[ii]=hits[ii];
			Mptot[ii]=Mp[ii];
		}
#endif


		if (rank == 0) {

			for (long ii=0;ii<npixeff;ii++)
				if (Mptot[ii] == 0)
					printf("ERROR: Mp[%ld] has elements = 0\n",ii);


			for (long ii=0;ii<npixeff;ii++)
				Mptot[ii] = 1.0/Mptot[ii];


			for (long ii=0;ii<npixeff;ii++)
				r[ii] = PNdtot[ii] - PtNPmatStot[ii];

			for (long ii=0;ii<npixeff;ii++)
				d[ii] =  Mptot[ii] * r[ii];


			delta_n = 0.0;
			for (long ii=0;ii<npixeff;ii++)
				delta_n += r[ii]*d[ii];

			var_n = 0.0;
			for (long ii=0;ii<npixeff;ii++)
				var_n += r[ii]*r[ii];


			delta0 = delta_n;
			var0 = var_n;
			printf("[%2.2i] var0 = %lf\n",rank, var0);

		}

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&var0,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&var_n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(d,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

		printf("[%2.2i] Main Conjugate gradient loop started\n",rank);


		//start loop
		iter = 0;
		while(((iter < 2000) && (var_n/var0 > 1e-10) && (idupl || !flgdupl)) || (!idupl && var_n/var0 > 1e-4)){
			// added brackets in order to avoid warning, mat-27/05

			fill(q,q+npixeff,0.0);
			fill(qtot,qtot+npixeff,0.0);
			cout << "dans while\n";
			exit(0);


			prefixe = "fPs";

			for (long iframe=iframe_min;iframe<iframe_max;iframe++){
				ns = nsamples[iframe];
				ff = fframes[iframe];
				f_lppix_Nk = fcut[iframe]*double(ns)/fsamp;

				if (CORRon){
					write_tfAS(d,indpix,nn,npix,flgdupl,factdupl, poutdir,termin,ff,ns,ndet,iframe);
					// read pointing + deproject + fourier transform

					do_PtNd(q,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,f_lppix_Nk,
							fsamp,ff,ns,ndet,size_det,rank_det,indpix,nn,npix,iframe,NULL,NULL);
					// return Pnd = At N-1 d
				} else {

					do_PtNPS_nocorr(d,extentnoiseSp_all,noiseSppreffile,poutdir,termin,dirfile,bolonames,
							f_lppix_Nk,fsamp,flgdupl,factdupl,ff,ns,ndet,size_det,rank_det,indpix,
							nn,npix,iframe,q,NULL,NULL);
				}
			} // end of iframe loop




#ifdef USE_MPI
			MPI_Reduce(q,qtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
			for(long ii=0;ii<npix;ii++)
				qtot[ii]=q[ii]; // ajout mat 02/06
#endif



			if (rank == 0){
				rtq= 0.0;
				for (long ii=0;ii<npixeff;ii++)
					rtq += qtot[ii] * d[ii];

				alpha = delta_n/rtq;


				for (long ii=0;ii<npixeff;ii++)
					S[ii] += alpha*d[ii];
			}

#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(S ,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif



			if ((iter % 10) == 0){

				fill(PtNPmatS,PtNPmatS+npixeff,0.0);
				fill(PtNPmatStot,PtNPmatStot+npixeff,0.0);
				prefixe = "fPs";

				for (long iframe=iframe_min;iframe<iframe_max;iframe++){
					ns = nsamples[iframe];
					ff = fframes[iframe];
					f_lppix_Nk = fcut[iframe]*double(ns)/fsamp;

					if (CORRon){
						write_tfAS(S,indpix,nn,npix,flgdupl,factdupl, poutdir,termin,ff,ns,ndet,iframe);
						// read pointing + deproject + fourier transform

						do_PtNd(PtNPmatS,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,
								f_lppix_Nk,fsamp,ff,ns,ndet,size_det,rank_det,indpix,nn,npix,iframe,
								NULL,NULL);
						// return Pnd = At N-1 d
					} else {

						do_PtNPS_nocorr(S,extentnoiseSp_all,noiseSppreffile,poutdir,termin,dirfile,bolonames,
								f_lppix_Nk,fsamp,flgdupl,factdupl,ff,ns,ndet,size_det,rank_det,
								indpix,nn,npix,iframe,PtNPmatS,NULL,NULL);
					}
				} // end of iframe loop



#ifdef USE_MPI
				MPI_Reduce(PtNPmatS,PtNPmatStot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
				for(long ii=0;ii<npixeff;ii++)
					PtNPmatStot[ii]=PtNPmatS[ii];
#endif


				if (rank == 0){
					for (long ii=0;ii<npixeff;ii++)
						r[ii] = PNdtot[ii] - PtNPmatStot[ii];
				}


			} else {

				if (rank == 0){
					for (long ii=0;ii<npixeff;ii++)
						r[ii] -= alpha*qtot[ii];
				}
			}





			if (rank == 0){

				for (long ii=0;ii<npixeff;ii++)
					s[ii] = Mptot[ii]*r[ii];


				delta_o = delta_n;

				delta_n = 0.0;
				for (long ii=0;ii<npixeff;ii++)
					delta_n += r[ii]*s[ii];

				var_n = 0.0;
				for (long ii=0;ii<npixeff;ii++)
					var_n += r[ii]*r[ii];



				beta = delta_n/delta_o;
				for (long ii=0;ii<npixeff;ii++)
					d[ii] = s[ii] + beta*d[ii];


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



				if (iter == 0){
					for (long ii=0; ii<nn; ii++) {
						for (long jj=0; jj<nn; jj++) {
							mi = jj*nn + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = Mptot[indpix[mi]];
							} else {
								map1d[mi] = 0.0;
							}
						}
					}

					fname = '!' + outdir + "optimMap_" + termin + "_noisevar.fits";
					write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


					for (long ii=0; ii<nn ; ii++){
						for (long jj=0; jj<nn; jj++){
							mi = jj*nn + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = Mptot[indpix[mi]] * PNdtot[indpix[mi]];
							} else {
								map1d[mi] = 0.0;
							}
						}
					}
					fname = '!' + outdir + "binMap_" + termin + "_flux.fits";
					write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);



					for (long ii=0; ii<nn; ii++) {
						for (long jj=0; jj<nn; jj++) {
							mi = jj*nn + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = hitstot[indpix[mi]];
							} else {
								map1d[mi] = 0.0;
							}
						}
					}

					if (addnpix){
						for (long iframe = 0;iframe<ntotscan;iframe++){
							for (long ii=0; ii<nn; ii++) {
								for (long jj=0; jj<nn; jj++) {
									mi = jj*nn + ii;
									if ((indpsrc[mi] != -1) && (indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc] != -1))
										map1d[mi] += hitstot[indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc]];
								}
							}
						}
					}

					fname = '!' + outdir + "optimMap_" + termin + "_hits.fits";
					write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


					for (long ii=0; ii<nn; ii++) {
						for (long jj=0; jj<nn; jj++) {
							mi = jj*nn + ii;
							map1d[mi] = 0.0;
						}
					}

					if (addnpix){
						for (long iframe = 0;iframe<ntotscan;iframe++){
							for (long ii=0; ii<nn; ii++) {
								for (long jj=0; jj<nn; jj++) {
									mi = jj*nn + ii;
									if ((indpsrc[mi] != -1) && (indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc] != -1))
										map1d[mi] += 1.0/Mptot[indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc]];
								}
							}
						}

						fname = '!' + outdir + "optimMap_" + termin + "_invnoisevaruncpix.fits";
						write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
					}

				}


				if (iterw && (iter % iterw) == 0){

					// make the map
					for (long ii=0; ii<nn; ii++) {
						for (long jj=0; jj<nn; jj++) {
							mi = jj*nn + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = -S[indpix[mi]];
							} else {
								map1d[mi] = 0.0;
							}
						}
					}

					sprintf(iterchar,"%d",iter);
					iterstr = iterchar;
					fname = '!' + outdir + "optimMap_" + termin + "_flux" + iterstr + "b.fits";
					write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


					if (flgdupl){
						for (long ii=0; ii<nn; ii++) {
							for (long jj=0; jj<nn; jj++) {
								mi = jj*nn + ii;
								if (indpix[mi] >= 0){
									map1d[mi] = -S[indpix[mi+nn*nn]]; //-finalmap[ii][jj];
								} else {
									map1d[mi] = 0.0;
								}
							}
						}


						sprintf(iterchar,"%d",iter);
						iterstr = iterchar;
						fname = '!' + outdir + "optimMap_" + termin + "_fluxflags" + iterstr + "b.fits";
						write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
					}
					if (addnpix){
						for (long ii=0; ii<nn ; ii++){
							for (long jj=0; jj<nn ; jj++){
								mi = jj*nn + ii;
								map1d[mi] = 0.0;
							}
						}
						for (long iframe = 0;iframe<ntotscan;iframe++){
							for (long ii=0; ii<nn; ii++) {
								for (long jj=0; jj<nn; jj++) {
									mi = jj*nn + ii;
									if ((indpsrc[mi] != -1)  && (indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc] != -1)){
										map1d[mi] += -S[indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc]]/Mptot[indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc]];
									}
								}
							}
						}


						sprintf(iterchar,"%d",iter);
						iterstr = iterchar;
						fname = '!' + outdir + "optimMap_" + termin + "_fluxuncpix_" + iterstr + "b.fits";
						write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);

					}
				}



				sprintf(testfile,"%s%s%s%s",outdir.c_str(),"ConvFile_",termin.c_str(),".txt");
				fp = fopen(testfile,"a");
				fprintf(fp,"iter = %d, crit = %10.15g, crit2 = %10.15g\n",iter,var_n/var0, delta_n/delta0);
				fclose(fp);

			}

#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&var_n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(d ,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

			iter++;

		} // end of while loop
		printf("\n");




		/*	sprintf(testfile,"%s%s%s%s",outdir.c_str(),"testfile_",termin.c_str(),".txt");
		fp = fopen(testfile,"a");
		fprintf(fp,"test en sortie de la boucle while \n");
		fclose(fp);*/





		if  ((projgaps || (flgdupl)) && !idupl){


			fill(PNd,PNd+npix,0.0);
			fill(PNdtot,PNdtot+npix,0.0);
			prefixe = "fdata";

			for (long iframe=iframe_min;iframe<iframe_max;iframe++){

				ns = nsamples[iframe];
				ff = fframes[iframe];
				f_lppix = f_lp*double(ns)/fsamp;
				f_lppix_Nk = fcut[iframe]*double(ns)/fsamp;

				if (CORRon){

					write_ftrProcesdata(S,indpix,indpsrc,nn,npix,npixsrc,ntotscan,addnpix,flgdupl,factdupl,2,
							poutdir,termin,errarcsec,dirfile,scerr_field,flpoint_field,bolonames,
							bextension,fextension,shift_data_to_point,f_lppix,ff,ns,
							napod,ndet,NORMLIN,NOFILLGAP,remove_polynomia,iframe);
					// fillgaps + butterworth filter + fourier transform
					// "fdata_" files generation (fourier transform of the data)

					do_PtNd(PNd,extentnoiseSp_all,noiseSppreffile,poutdir,prefixe,termin,bolonames,f_lppix_Nk,
							fsamp,ff,ns,ndet,size_det,rank_det,indpix,nn,npix,iframe,NULL,NULL);
					// return Pnd = At N-1 d
				} else {

					do_PtNd_nocorr(PNd,extentnoiseSp_all,noiseSppreffile,poutdir,termin,errarcsec,dirfile,
							scerr_field,flpoint_field,bolonames,bextension,fextension,
							shift_data_to_point,f_lppix,f_lppix_Nk,fsamp,ntotscan,addnpix,
							flgdupl,factdupl,2,ff,ns,napod,ndet,size_det,rank_det,indpix,indpsrc,
							nn,npix,npixsrc,NORMLIN,NOFILLGAP,remove_polynomia,iframe,S);
				}
			} // end of iframe loop




#ifdef USE_MPI
			MPI_Reduce(PNd,PNdtot,npix,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
			for(long ii=0;ii<npix;ii++)
				PNdtot[ii]=PNd[ii]; // ajout Mat 02/07
#endif
		}



	}// end of idupl loop



#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	sprintf(testfile,"%s%s%s%s",outdir.c_str(),"testfile_",termin.c_str(),".txt");
	fp = fopen(testfile,"a");
	fprintf(fp,"test avant ecriture \n");
	fclose(fp);



	// write map function : nn, S, indpix, outdir, termin, addnpix, ntotscan,
	// pixdeg, tancoord, tanpix, coordsyst, Mptot, indpsrc, npixsrc, factdupl,
	//

	//******************************  write final map in file ********************************



	if (rank == 0){

		printf(" after CC INVERSION %d\n",npix*(npix+1)/2);





		for (long ii=0; ii<nn; ii++) {
			for (long jj=0; jj<nn; jj++) {
				mi = jj*nn + ii;
				if (indpix[mi] >= 0){
					map1d[mi] = -S[indpix[mi]];
				} else {
					map1d[mi] = 0.0;
				}
			}
		}


		fname = '!' + outdir + "optimMap_" + termin + "_flux.fits";
		write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


		for (long ii=0; ii<nn; ii++) {
			for (long jj=0; jj<nn; jj++) {
				mi = jj*nn + ii;
				if (indpix[mi] >= 0){
					map1d[mi] = Mptot[indpix[mi]];
				} else {
					map1d[mi] = 0.0;
				}
			}
		}


		fname = '!' + outdir + "optimMap_" + termin + "_noisevar.fits";
		write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);

		if (addnpix){
			for (long iframe = 0;iframe<ntotscan;iframe++){
				fru = 0;
				for (long ii=0; ii<nn; ii++) {
					for (long jj=0; jj<nn; jj++) {
						mi = jj*nn + ii;
						if ((indpsrc[mi] != -1)  && (indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc] != -1)){
							fru = 1;
							map1d[mi] = -S[indpix[indpsrc[mi]+factdupl*nn*nn+npixsrc*iframe]];
						} else {
							map1d[mi] = 0.0;
						}
					}
				}


				if (fru){
					sprintf(iframechar,"%ld",iframe);
					iframestr = iframechar;
					fname = '!' + outdir + "optimMap_" + termin + "_flux_fr" + iframestr + ".fits";
					write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


					for (long ii=0; ii<nn; ii++) {
						for (long jj=0; jj<nn; jj++) {
							mi = jj*nn + ii;
							if ((indpsrc[mi] != -1)  && (indpix[indpsrc[mi]+factdupl*nn*nn+iframe*npixsrc] != -1)){
								map1d[mi] = Mptot[indpix[indpsrc[mi]+factdupl*nn*nn+npixsrc*iframe]];
							} else {
								map1d[mi] = 0.0;
							}
						}
					}

					fname = '!' + outdir + "optimMap_" + termin + "_noisevar_fr" + iframestr + ".fits";
					write_fits(fname, pixdeg, nn, nn, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
				}
			}
		}
	}

	// end of write map function



	delete [] r;
	delete [] q;
	delete [] qtot;
	delete [] d;
	delete [] Mp;
	delete [] Mptot;
	delete [] s;
	delete [] PtNPmatS;
	delete [] PtNPmatStot;
	delete [] hits;
	delete [] hitstot;
	delete [] map1d;

}
