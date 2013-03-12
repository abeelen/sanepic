/*
 * frame_order.cpp
 *
 *  Created on: Dec 6, 2012
 *      Author: abeelen
 */


#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
#include <cstdio>

#include <assert.h>
#include <stdint.h>

using namespace std;

#include "FrameOrder.h"
#include "ErrorCode.h"
#include "InputFileIO.h"
#include "Utilities.h"
#include "MPIConfiguration.h"

uint16_t writeParallelScheme(string outdir, vector<long> order, struct samples samples_struct)
// Write the Parallelization Scheme for further use.
{

	ofstream file;
	string fname = outdir + parallelScheme_filename;

	cout << endl << "II - Scheme written in file : " << fname << endl;

	file.open(fname.c_str(), ios::out);
	if(!file.is_open()){
		cerr << "Error writing file [" << fname << "]" << endl;
		return FILE_PROBLEM;
	}

	for (long ii=0;ii<samples_struct.ntotscan;ii++){
		file << Basename(samples_struct.fitsvect[ii]) << " " << order[ii] << endl;
	}

	file.close();

	return EXIT_SUCCESS;
}

void distributeFrames(vector<long> nsamples, vector<long> int_weight, vector<float> float_weight, vector<long> & order){
  /**
   *  Distribute frames of [nsamples_i] between nodes with integer
   * (nCPU per node for example) weight and float weight (relative
   *  speed of nodes/per CPU)
   */

  assert(int_weight.size() == float_weight.size());

  long nFrames = nsamples.size();
  long nNode   = int_weight.size();

  order.clear();
  order.resize(nFrames);

  for (long ii=0; ii< nFrames; ii++)
    order[ii] = ii;

  // pseudoSample_perProc contains the weighted number of samples per CPU
  vector<float> pseudoSample_perProc(nNode);
  for (long ii=0; ii< nNode; ii++)
    pseudoSample_perProc[ii] = 0;

  vector<size_t> index(nFrames);
  // Find indexes of sorted nsamples
  paired_sort(nsamples, index );
  // inverse it (bigger first)
  reverse(index.begin(), index.end());

  vector<size_t> index_float(nNode);
  paired_sort(float_weight, index_float);
  reverse(index_float.begin(), index_float.end());

  // assign the biggest frames first, onto the most weighted node...
  for (long iFrame = 0; iFrame < min(nNode, nFrames); iFrame++){
    order[index[iFrame]] = index_float[iFrame];
    pseudoSample_perProc[iFrame] = ceil(nsamples[index[iFrame]]*1.0/int_weight[index_float[iFrame]])*1.0/float_weight[index_float[iFrame]];
  }

  long idCPU;
  vector<float>::iterator it;
  // ... then adjust with smaller frames by...
  for (long iFrame = nNode; iFrame < nFrames; iFrame++){
    // ... findind the index of the minimum value in pseudoSample_perProc
    it =   find( pseudoSample_perProc.begin(), pseudoSample_perProc.end(), *( min_element(pseudoSample_perProc.begin(), pseudoSample_perProc.end()) ) );
    idCPU = distance(pseudoSample_perProc.begin(), it);
    // ... and assign it to it...
    order[index[iFrame]] = idCPU;
    pseudoSample_perProc[idCPU] += ceil(nsamples[index[iFrame]]*1.0/int_weight[idCPU])*1.0/float_weight[idCPU];
  }

}

void  AssignNodeFloatWeight(vector<string> nodeName, map<string, float> map_nodeWeight, vector<float> & nodeWeight){
  /**
   * Take nodeName corresponding weight in nodeWeight associative array, store into float_weight
   *
   */

  // default value of 1:
  nodeWeight.clear();
  nodeWeight.assign(nodeName.size(), 1);

  for( map<string, float>::iterator ii=map_nodeWeight.begin(); ii !=map_nodeWeight.end(); ++ii){
    string ii_nodeName = (*ii).first;
    float ii_weight = (*ii).second;

    for (vector<string>::iterator jj=nodeName.begin(); jj != nodeName.end(); ++jj){
      if ( *jj == ii_nodeName ){
	long id = distance(nodeName.begin(), jj);
	nodeWeight[id] = ii_weight;
      }
    }
  }

}

uint16_t removeProcName(std::string outdir){
	string filename;
	filename = outdir+processorName_filename;
	return remove(filename.c_str());

}

uint16_t readNodeWeight(string & output, string pathIn, map<string, float> & nodeWeight){
  /**
   * read node weight into an associative hash
   */

  string fname("node.weight");

  vector<string> output_string;
  vector<float>  output_float;

  // if not found, no big deal...
  if ( read_file_2col(output, pathIn+fname,  output_string, output_float) ){
	  output += "WW - no node weight found, using default 1\n";
  }

  nodeWeight.clear();
  for (unsigned int ii=0; ii< output_string.size(); ii++)
    nodeWeight[output_string[ii]] = output_float[ii];

  return 0;

}
