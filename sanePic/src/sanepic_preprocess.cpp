/*
 * sanepic_preprocess.cpp
 *
 *  Created on: 2 juil. 2009
 *      Author: matthieu
 */


#include "sanepic_preprocess.h"

void sanepic_preprocess(int nn, std::vector<long> xxi, std::vector<long> xxf,
		std::vector<long> yyi, std::vector<long> yyf, long *&indpsrc, long &npixsrc,
		long ntotscan, long &addnpix,int &npix, int &factdupl, bool flgdupl, string termin,
		string outdir, double *&PNdtot, long *&indpix, int &flagon){

	//char testfile[100];
	unsigned char *mask; // samples flags, pointing flags, rejected samples list
	//FILE *fp;
	int npix2;

	//************************************* Deal with masking the point sources
	// define the mask
	mask = new unsigned char[nn*nn];
	for (long ii=0;ii<nn*nn;ii++)
		mask[ii] = 1;


	if (xxi.size() != 0){
		for (long ib = 0;ib < (long)xxi.size(); ib++){ // to avoid warning, mat-27/05
			// for each box crossing constraint removal
			for (long ii=xxi[ib];ii<xxf[ib];ii++)
				for (long ll=yyi[ib];ll<yyf[ib];ll++)
					mask[ll*nn + ii] = 0;  // mask is initialised to 0
		}
	}



	//long npixsrc = 0;

	for (long ii=0;ii<nn*nn;ii++){
		if (mask[ii] == 0){
			indpsrc[ii] = npixsrc;
			npixsrc += 1;
		} else {
			indpsrc[ii] = -1;
		}
	}
	addnpix = ntotscan*npixsrc; // addnpix = number of pix to add in pixon = number of scans * number of pix in box crossing constraint removal

	cout  << "addnpix : " << addnpix << endl;
	//*************************************************************************************************//

	if (flgdupl) factdupl = 2; // -M =1, default 0 : if flagged data are put in a duplicated map

	// read npix, PNdtot from file
	read_PNd(PNdtot, npix, termin, outdir);
	/*for (ii=0;ii<20;ii++)
		cout << PNdtot[ii] << " ";
	cout << endl << "avant read indpix\n";
	exit(0);*/


	indpix=new long[factdupl*nn*nn+2 + addnpix];
	// read indpix

	read_indpix(factdupl*nn*nn+2 + addnpix, npix2, indpix, termin, outdir, flagon);

	if (npix!=npix2)
		cout << "Warning ! Indpix_for_conj_grad.bi and PNdCorr_*.bi are not compatible, npix!=npix2" << endl;



	/*cout << "flagon : " << flagon << endl;
	for (ii=0;ii<nn*nn;ii=ii+nn+1)
		cout << indpix[ii] << " ";
	cout << endl;
	exit(1);*/
	/*
	if (rank == 0){
		sprintf(testfile,"%s%s%s%s",outdir.c_str(),"testfile_",termin.c_str(),".txt");
		fp = fopen(testfile,"a");
		for (ii=0;ii<=ntotscan;ii++) fprintf(fp,"frnum[%ld] = %ld \n",ii,frnum[ii]);
		fclose(fp);
	}*/

	//end of PREPROC

	delete [] mask;

}
