/*
 * conjugate_gradient.cpp
 *
 *  Created on: 3 juil. 2009
 *      Author: matthieu
 */


#include "conjugate_gradient.h"


void sanepic_conjugate_gradient(bool flgdupl, int npix, double* &S,long iframe_min, long iframe_max,
		long *nsamples, long *fframes,std::vector<double> fcut,double f_lp,double fsamp,
		long *indpix, int NAXIS1, int NAXIS2, int factdupl, string tmp_dir, string termin, string termin_internal, long ndet,
		string *extentnoiseSp_all,string noiseSppreffile, std::vector<string> bolonames,/* int size_det,
		int rank_det,*/ int iterw, double pixdeg, double *tancoord, double *tanpix,int coordsyst,
		long *indpsrc, long npixsrc, int flagon, bool projgaps, int rank, bool CORRon,
		string dirfile, double *&PNdtot, long ntotscan,long addnpix,bool NORMLIN,bool NOFILLGAP,
		long napod,int shift_data_to_point,bool remove_polynomia,/*string fextension,string bextension,
		string flpoint_field,string scerr_field,*/ string outdir, string *fits_table){



	FILE *fp;
	//char testfile[100];


	//char iterchar[30];
	//char iframechar[30];
	//string iterstr, iframestr;
	bool fru;
	string fname, testfile;
	ostringstream temp_stream;
	double *PtNPmatS,  *PtNPmatStot, *r, *q, *qtot, *d, *Mp, *Mptot, *s;
	// Mp = M in the paper = preconditioner
	long *hits, *hitstot;
	double *PNd;

	double var0 = 0.0, var_n = 0.0, delta0 = 0.0, delta_n = 0.0, alpha = 0.0;
	double delta_o, rtq, beta;

	long mi;
	double *map1d;
	string prefixe;
	int iter;
	double  f_lppix_Nk,f_lppix;
	long ns,ff;
	int npixeff;

	//double errarcsec = 15.0; // rejection criteria : scerr[ii] > errarcsec, sample is rejected
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
	PNd=new double[npix];

	map1d = new double[NAXIS1*NAXIS2];

	for (int idupl = 0;idupl<=flgdupl;idupl++){

		//for (long ii=0;ii<npix;ii++) S[ii] = 0.0;//PNd[ii]; // useless car fill dans main

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


		prefixe = "fPs_";

		for (long iframe=iframe_min;iframe<iframe_max;iframe++){
			ns = nsamples[iframe];
			ff = fframes[iframe];
			f_lppix_Nk = fcut[iframe]*double(ns)/fsamp;


			//    cout << "[" << rank << "] " << iframe << "/" << iframe_max << endl;


			// preconditioner computation : Mp
			if (CORRon){
				write_tfAS(S,indpix,NAXIS1, NAXIS2,npix,flgdupl,factdupl, tmp_dir,ff,ns,ndet,iframe);
				// read pointing + deproject + fourier transform

				do_PtNd(PtNPmatS,extentnoiseSp_all,noiseSppreffile,tmp_dir,prefixe,bolonames,
						f_lppix_Nk,fsamp,ff,ns,ndet,/*size_det,rank_det,*/indpix,NAXIS1, NAXIS2,npix,iframe,Mp,hits);
				// return Pnd = At N-1 d
			} else {

				do_PtNPS_nocorr(S,extentnoiseSp_all,noiseSppreffile,tmp_dir,dirfile,bolonames,
						f_lppix_Nk,fsamp,flgdupl,factdupl,ff,ns,ndet,/*size_det,rank_det,*/indpix,
						NAXIS1, NAXIS2,npix,iframe,PtNPmatS,Mp,hits);
			}

		} // end of iframe loop

		//cout << "do_ptnd" << endl;


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

			for (long ii=0;ii<npixeff;ii++)
				d[ii] =  Mptot[ii] * r[ii]; // d = M-1 * r


			delta_n = 0.0;
			for (long ii=0;ii<npixeff;ii++)
				delta_n += r[ii]*d[ii]; // delta_new = rT * d

			var_n = 0.0;
			for (long ii=0;ii<npixeff;ii++)
				var_n += r[ii]*r[ii];


			delta0 = delta_n; // delta_0 <= delta_new
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
		iter = 0; // max iter = 2000, but ~100 iterations are required to achieve convergence

		// while i<imax and var_new > epsilon² * var_0 : epsilon² = 1e-10 => epsilon = 1e-5
		while(((iter < 2000) && (var_n/var0 > 1e-10) && (idupl || !flgdupl)) || (!idupl && var_n/var0 > 1e-4)){
			// added brackets in order to avoid warning, mat-27/05

			fill(q,q+npixeff,0.0); // q <= A*d
			fill(qtot,qtot+npixeff,0.0);

			prefixe = "fPs_";

			for (long iframe=iframe_min;iframe<iframe_max;iframe++){
				ns = nsamples[iframe];
				ff = fframes[iframe];
				f_lppix_Nk = fcut[iframe]*double(ns)/fsamp;

				if (CORRon){
					write_tfAS(d,indpix,NAXIS1, NAXIS2,npix,flgdupl,factdupl, tmp_dir,ff,ns,ndet,iframe);
					// read pointing + deproject + fourier transform

					do_PtNd(q,extentnoiseSp_all,noiseSppreffile,tmp_dir,prefixe,bolonames,f_lppix_Nk,
							fsamp,ff,ns,ndet,/*size_det,rank_det,*/indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL);
					// return Pnd = At N-1 d
				} else {

					do_PtNPS_nocorr(d,extentnoiseSp_all,noiseSppreffile,tmp_dir,dirfile,bolonames,
							f_lppix_Nk,fsamp,flgdupl,factdupl,ff,ns,ndet,/*size_det,rank_det,*/indpix,
							NAXIS1, NAXIS2,npix,iframe,q,NULL,NULL);
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
					rtq += qtot[ii] * d[ii]; // rtq = (dT * q)

				alpha = delta_n/rtq; // alpha <= delta_new / (dT * q)


				for (long ii=0;ii<npixeff;ii++)
					S[ii] += alpha*d[ii]; // x = x + alpha * d, x = S = signal
			}

#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(S ,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif



			if ((iter % 10) == 0){ // if iter is divisible by 10, recompute PtNPmatStot

				fill(PtNPmatS,PtNPmatS+npixeff,0.0);
				fill(PtNPmatStot,PtNPmatStot+npixeff,0.0);
				prefixe = "fPs_";

				for (long iframe=iframe_min;iframe<iframe_max;iframe++){
					ns = nsamples[iframe];
					ff = fframes[iframe];
					f_lppix_Nk = fcut[iframe]*double(ns)/fsamp;

					if (CORRon){
						write_tfAS(S,indpix,NAXIS1, NAXIS2,npix,flgdupl,factdupl, tmp_dir,ff,ns,ndet,iframe);
						// read pointing + deproject + fourier transform

						do_PtNd(PtNPmatS,extentnoiseSp_all,noiseSppreffile,tmp_dir,prefixe,bolonames,
								f_lppix_Nk,fsamp,ff,ns,ndet,/*size_det,rank_det,*/indpix,NAXIS1, NAXIS2,npix,iframe,
								NULL,NULL);
						// return Pnd = At N-1 d
					} else {

						do_PtNPS_nocorr(S,extentnoiseSp_all,noiseSppreffile,tmp_dir,dirfile,bolonames,
								f_lppix_Nk,fsamp,flgdupl,factdupl,ff,ns,ndet,/*size_det,rank_det,*/
								indpix,NAXIS1, NAXIS2,npix,iframe,PtNPmatS,NULL,NULL);
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



				if (iter == 0){
					for (long ii=0; ii<NAXIS1; ii++) {
						for (long jj=0; jj<NAXIS2; jj++) {
							mi = jj*NAXIS2 + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = Mptot[indpix[mi]];
							} else {
								map1d[mi] = 0.0;
							}
						}
					}

					fname = '!' + outdir + "optimMap_" + termin + "_noisevar.fits"; // write preconditioner
					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


					for (long ii=0; ii<NAXIS1 ; ii++){
						for (long jj=0; jj<NAXIS2; jj++){
							mi = jj*NAXIS1 + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = Mptot[indpix[mi]] * PNdtot[indpix[mi]];
							} else {
								map1d[mi] = 0.0;
							}
						}
					}
					fname = '!' + outdir + "binMap_" + termin + "_flux.fits";
					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);



					for (long ii=0; ii<NAXIS1; ii++) {
						for (long jj=0; jj<NAXIS2; jj++) {
							mi = jj*NAXIS1 + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = hitstot[indpix[mi]];
							} else {
								map1d[mi] = 0.0;
							}
						}
					}

					if (addnpix){
						for (long iframe = 0;iframe<ntotscan;iframe++){
							for (long ii=0; ii<NAXIS1; ii++) {
								for (long jj=0; jj<NAXIS2; jj++) {
									mi = jj*NAXIS1 + ii;
									if ((indpsrc[mi] != -1) && (indpix[indpsrc[mi]+factdupl*NAXIS1*NAXIS2+iframe*npixsrc] != -1))
										map1d[mi] += hitstot[indpix[indpsrc[mi]+factdupl*NAXIS2*NAXIS2+iframe*npixsrc]];
								}
							}
						}
					}

					fname = '!' + outdir + "optimMap_" + termin + "_hits.fits";
					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


					for (long ii=0; ii<NAXIS1; ii++) {
						for (long jj=0; jj<NAXIS2; jj++) {
							mi = jj*NAXIS1 + ii;
							map1d[mi] = 0.0;
						}
					}

					if (addnpix){
						for (long iframe = 0;iframe<ntotscan;iframe++){
							for (long ii=0; ii<NAXIS1; ii++) {
								for (long jj=0; jj<NAXIS2; jj++) {
									mi = jj*NAXIS1 + ii;
									if ((indpsrc[mi] != -1) && (indpix[indpsrc[mi]+factdupl*NAXIS1*NAXIS2+iframe*npixsrc] != -1))
										map1d[mi] += 1.0/Mptot[indpix[indpsrc[mi]+factdupl*NAXIS1*NAXIS2+iframe*npixsrc]];
								}
							}
						}

						fname = '!' + outdir + "optimMap_" + termin + "_invnoisevaruncpix.fits";
						write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
					}

				}


				if (iterw && (iter % iterw) == 0){

					// make the map
					for (long ii=0; ii<NAXIS1; ii++) {
						for (long jj=0; jj<NAXIS2; jj++) {
							mi = jj*NAXIS1 + ii;
							if (indpix[mi] >= 0){
								map1d[mi] = S[indpix[mi]];
							} else {
								map1d[mi] = 0.0;
							}
						}
					}

					/*sprintf(iterchar,"%d",iter);
					iterstr = iterchar;
					fname = '!' + outdir + "optimMap_" + termin + "_flux" + iterstr + "b.fits";*/
					temp_stream << "!" + outdir + "optimMap_" + termin + "_flux" << iter << "b.fits";

					// récupérer une chaîne de caractères
					fname= temp_stream.str();
					// Clear ostringstream buffer
					temp_stream.str("");
					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


					// TODO : check indices... seems wrong... (OK for me : mat 30/09)
					if (flgdupl){
						for (long ii=0; ii<NAXIS1; ii++) {
							for (long jj=0; jj<NAXIS2; jj++) {
								mi = jj*NAXIS1 + ii;
								if (indpix[mi] >= 0){
									map1d[mi] = -S[indpix[mi+NAXIS1*NAXIS2]]; //-finalmap[ii][jj];
								} else {
									map1d[mi] = 0.0;
								}
							}
						}


						/*sprintf(iterchar,"%d",iter);
						iterstr = iterchar;
						fname = '!' + outdir + "optimMap_" + termin + "_fluxflags" + iterstr + "b.fits";*/
						temp_stream << "!" + outdir + "optimMap_" + termin + "_fluxflags" << iter << "b.fits";

						// récupérer une chaîne de caractères
						fname= temp_stream.str();
						// Clear ostringstream buffer
						temp_stream.str("");
						write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
					}
					if (addnpix){
						for (long ii=0; ii<NAXIS1 ; ii++){
							for (long jj=0; jj<NAXIS2 ; jj++){
								mi = jj*NAXIS1 + ii;
								map1d[mi] = 0.0;
							}
						}
						for (long iframe = 0;iframe<ntotscan;iframe++){
							for (long ii=0; ii<NAXIS1; ii++) {
								for (long jj=0; jj<NAXIS2; jj++) {
									mi = jj*NAXIS1 + ii;
									if ((indpsrc[mi] != -1)  && (indpix[indpsrc[mi]+factdupl*NAXIS1*NAXIS2+iframe*npixsrc] != -1)){
										map1d[mi] += -S[indpix[indpsrc[mi]+factdupl*NAXIS1*NAXIS2+iframe*npixsrc]]/Mptot[indpix[indpsrc[mi]+factdupl*NAXIS1*NAXIS2+iframe*npixsrc]];
									}
								}
							}
						}


						/*sprintf(iterchar,"%d",iter);
						iterstr = iterchar;
						fname = '!' + outdir + "optimMap_" + termin + "_fluxuncpix_" + iterstr + "b.fits";*/
						temp_stream << "!" + outdir + "optimMap_" + termin + "_fluxuncpix_" << iter << "b.fits";

						// récupérer une chaîne de caractères
						fname= temp_stream.str();
						// Clear ostringstream buffer
						temp_stream.str("");
						write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);

					}
				}



				//	sprintf(testfile,"%s%s%s%s",outdir.c_str(),"ConvFile_",termin.c_str(),".txt");
				temp_stream << outdir + "ConvFile_" + termin + "txt";

				// récupérer une chaîne de caractères
				testfile= temp_stream.str();
				// Clear ostringstream buffer
				temp_stream.str("");
				fp = fopen(testfile.c_str(),"a");
				fprintf(fp,"iter = %d, crit = %10.15g, crit2 = %10.15g\n",iter,var_n/var0, delta_n/delta0);
				fclose(fp);

			}

#ifdef USE_MPI
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Bcast(&var_n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(d ,npix,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

			iter++; // i = i +1

		} // end of while loop
		printf("\n");




		/*	sprintf(testfile,"%s%s%s%s",outdir.c_str(),"testfile_",termin.c_str(),".txt");
		fp = fopen(testfile,"a");
		fprintf(fp,"test en sortie de la boucle while \n");
		fclose(fp);*/





		if  ((projgaps || (flgdupl)) && !idupl){


			fill(PNd,PNd+npix,0.0);
			fill(PNdtot,PNdtot+npix,0.0);
			prefixe = "fdata_";

			for (long iframe=iframe_min;iframe<iframe_max;iframe++){

				ns = nsamples[iframe];
				ff = fframes[iframe];
				f_lppix = f_lp*double(ns)/fsamp;
				f_lppix_Nk = fcut[iframe]*double(ns)/fsamp;

				if (CORRon){

					write_ftrProcesdata(S,indpix,indpsrc,NAXIS1, NAXIS2,npix,npixsrc,ntotscan,addnpix,flgdupl,factdupl,2,
							tmp_dir,dirfile,bolonames,fits_table,shift_data_to_point,f_lppix,ff,ns,
							napod,ndet,NORMLIN,NOFILLGAP,remove_polynomia,iframe);
					// fillgaps + butterworth filter + fourier transform
					// "fdata_" files generation (fourier transform of the data)

					do_PtNd(PNd,extentnoiseSp_all,noiseSppreffile,tmp_dir,prefixe,bolonames,f_lppix_Nk,
							fsamp,ff,ns,ndet,/*size_det,rank_det,*/indpix,NAXIS1, NAXIS2,npix,iframe,NULL,NULL);
					// return Pnd = At N-1 d
				} else {

					do_PtNd_nocorr(PNd,extentnoiseSp_all,noiseSppreffile,tmp_dir,dirfile,
							bolonames, fits_table, shift_data_to_point,f_lppix,f_lppix_Nk,fsamp,ntotscan,addnpix,
							flgdupl,factdupl,2,ff,ns,napod,ndet,/*size_det,rank_det,*/indpix,indpsrc,
							NAXIS1, NAXIS2,npix,npixsrc,NORMLIN,NOFILLGAP,remove_polynomia,iframe,S);
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

	//sprintf(testfile,"%s%s%s%s",outdir.c_str(),"testfile_",termin.c_str(),".txt");
	temp_stream << outdir + "testfile_" + termin_internal + "txt";

	// récupérer une chaîne de caractères
	testfile= temp_stream.str();
	// Clear ostringstream buffer
	temp_stream.str("");
	fp = fopen(testfile.c_str(),"a");
	fprintf(fp,"test avant ecriture \n");
	fclose(fp);



	// write map function : NAXIS1, NAXIS2, S, indpix, outdir, termin, addnpix, ntotscan,
	// pixdeg, tancoord, tanpix, coordsyst, Mptot, indpsrc, npixsrc, factdupl,
	//

	//******************************  write final map in file ********************************



	if (rank == 0){

		printf(" after CC INVERSION %d\n",npix*(npix+1)/2);





		for (long ii=0; ii<NAXIS1; ii++) {
			for (long jj=0; jj<NAXIS2; jj++) {
				mi = jj*NAXIS1 + ii; // mat 30/09
				if (indpix[mi] >= 0){
					map1d[mi] = S[indpix[mi]];
				} else {
					map1d[mi] = 0.0;
				}
			}
		}


		//////// test pour sanePS
		fp = fopen("test_signal_pic.txt","w");
		for (int i =0;i<npix;i++)
			fprintf(fp,"%lf\n",S[i]);
		fclose(fp);
		////////////////////////////////////

		fname = '!' + outdir + "optimMap_" + termin + "_flux.fits";
		write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


		for (long ii=0; ii<NAXIS1; ii++) {
			for (long jj=0; jj<NAXIS2; jj++) {
				mi = jj*NAXIS1 + ii;
				if (indpix[mi] >= 0){
					map1d[mi] = Mptot[indpix[mi]];
				} else {
					map1d[mi] = 0.0;
				}
			}
		}


		fname = '!' + outdir + "optimMap_" + termin + "_noisevar.fits";
		write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);

		if (addnpix){
			for (long iframe = 0;iframe<ntotscan;iframe++){
				fru = 0;
				for (long ii=0; ii<NAXIS1; ii++) {
					for (long jj=0; jj<NAXIS2; jj++) {
						mi = jj*NAXIS1 + ii;
						if ((indpsrc[mi] != -1)  && (indpix[indpsrc[mi]+factdupl*NAXIS1*NAXIS2+iframe*npixsrc] != -1)){
							fru = 1;
							map1d[mi] = -S[indpix[indpsrc[mi]+factdupl*NAXIS1*NAXIS2+npixsrc*iframe]];
						} else {
							map1d[mi] = 0.0;
						}
					}
				}


				if (fru){
					/*sprintf(iframechar,"%ld",iframe);
					iframestr = iframechar;
					fname = '!' + outdir + "optimMap_" + termin + "_flux_fr" + iframestr + ".fits";*/
					temp_stream << "!" + outdir + "optimMap_" + termin + "_flux_fr" << iframe << ".fits";

					// récupérer une chaîne de caractères
					fname= temp_stream.str();
					// Clear ostringstream buffer
					temp_stream.str("");
					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);


					for (long ii=0; ii<NAXIS1; ii++) {
						for (long jj=0; jj<NAXIS2; jj++) {
							mi = jj*NAXIS1 + ii;
							if ((indpsrc[mi] != -1)  && (indpix[indpsrc[mi]+factdupl*NAXIS1*NAXIS2+iframe*npixsrc] != -1)){
								map1d[mi] = Mptot[indpix[indpsrc[mi]+factdupl*NAXIS1*NAXIS2+npixsrc*iframe]];
							} else {
								map1d[mi] = 0.0;
							}
						}
					}

					//fname = '!' + outdir + "optimMap_" + termin + "_noisevar_fr" + iframestr + ".fits";
					temp_stream << "!" + outdir + "optimMap_" + termin + "_noisevar_fr" << iframe << ".fits";

					// récupérer une chaîne de caractères
					fname= temp_stream.str();
					// Clear ostringstream buffer
					temp_stream.str("");
					write_fits(fname, pixdeg, NAXIS1, NAXIS2, tancoord, tanpix, coordsyst, 'd', (void *)map1d);
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
