/*
 * sanepic_preprocess.cpp
 *
 *  Created on: 2 juil. 2009
 *      Author: matthieu
 */


#include "sanepic_preprocess.h"

void sanepic_preprocess(int nn, std::vector<long> xxi, std::vector<long> xxf,
		std::vector<long> yyi, std::vector<long> yyf, long *&indpsrc, long &npixsrc,
		long ntotscan, long &addnpix){

	//char testfile[100];
	unsigned char *mask; // samples flags, pointing flags, rejected samples list
	//FILE *fp;


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

	delete [] mask;

}

//end of PREPROC

