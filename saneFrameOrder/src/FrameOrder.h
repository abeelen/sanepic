/*
 * frameOrder.h
 *
 *  Created on: Dec 6, 2012
 *      Author: abeelen
 */

#ifndef FRAMEORDER_H_
#define FRAMEORDER_H_

#include <string>
#include <vector>
#include <map>

#include <stdint.h>

using namespace std;

void distributeFrames(vector<long> nsamples, vector<long> int_weight, vector<float> float_weight, vector<long> & order);
void  AssignNodeFloatWeight(vector<string> nodeName, map<string, float> map_nodeWeight, vector<float> & nodeWeight);

//! Writes the parallel scheme file to disk
/*!
  \param outdir Output directory : outdir/parallel_scheme
  \param position The processor indice array. This array will be filled by find_best_order_frames
  \param frnum Frame number (proc 0 has to do scans which frame's index goes from frnum(0) to frnum(1), proc 1 from frnum(1) to frnum(2), etc...)
  \param size The number of processors to be used
  \param samples_struct A samples structure that contains everything about frames, noise files and frame processing order
  \return An integer >0 if there were a problem, or 0 if everything went OK
 */
uint16_t writeParallelScheme(std::string outdir, vector<long> order, struct samples samples_struct);

uint16_t removeProcName(std::string outdir);

uint16_t readNodeWeight(string & output, string pathIn, map<string, float> & nodeWeight);

#endif /* FRAMEORDER_H_ */
