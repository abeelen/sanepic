/*
 * sanePre_preprocess.cpp
 *
 *  Created on: 25 juin 2009
 *      Author: matthieu
 */

#include "sanePre_preprocess.h"

#include "todprocess.h"
#include "map_making.h"

#include <iostream>
#include <vector>
#include <string>


using namespace std;

//
//
////TODO : Not used anymore...
//long Compute_indpsrc_addnpix(int NAXIS1, int NAXIS2, long ntotscan, std::vector<long> xxi, std::vector<long> xxf, std::vector<long> yyi,
//		std::vector<long> yyf, long* &indpsrc,long &npixsrc)
//{
//
//	long addnpix=0;
//
//	//************************************* Deal with masking the point sources
//	// define the mask
//	unsigned char *mask;
//
//	mask = new unsigned char[NAXIS1*NAXIS2];
//	for (long ii=0;ii<NAXIS1*NAXIS2;ii++)
//		mask[ii] = 1;
//
//
//	if (xxi.size() != 0){
//		for (long ib = 0;ib < (long)xxi.size(); ib++){ // to avoid warning, mat-27/05
//			// for each box crossing constraint removal
//			for (long ii=xxi[ib];ii<xxf[ib];ii++)
//				for (long ll=yyi[ib];ll<yyf[ib];ll++)
//					mask[ll*NAXIS1 + ii] = 0;  // mask is initialised to 0
//		}
//	}
//
//
//
//
//	indpsrc = new long[NAXIS1*NAXIS2];
//	for (long ii=0;ii<NAXIS1*NAXIS2;ii++){
//		if (mask[ii] == 0){
//			indpsrc[ii] = npixsrc;
//			npixsrc += 1;
//		} else {
//			indpsrc[ii] = -1;
//		}
//	}
//	addnpix = ntotscan*npixsrc; // addnpix = number of pix to add in pixon = number of scans * number of pix in box crossing constraint removal
//
//	//debug
//	cout  << "addnpix : " << addnpix << endl;
//
//	//clean up
//	delete [] mask;
//
//	return addnpix;
//}
