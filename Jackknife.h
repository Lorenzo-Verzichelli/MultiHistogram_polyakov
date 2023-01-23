#pragma once
#include <iostream>
#include <vector>

#include "weihted_mean.h"

void jackknife_mean_susc(std::vector<double>* obses, int nrun, std::vector<double>* weights, unsigned int block,
	double& mean, double& variance, double& err_mean, double& err_variance);
/*Computes mean and variance with jeckknife errors of the observable obses at emperature 1/beta, blocked with size block

	results will be stored in the arguments passed by reference
	passed pointers must point to nrun long array, wich will NOT be modified

	the basic idea comes from Bonati's python script for Multi-histogram
	(hopfully I have understud it correctly)
*/

void compute_weights(double beta, int nrun, double* betas, std::vector<double>* energies, double* logZs,
	std::vector<double>* weights);
/*Computes the LOG of weights for the jeckknife procedure
	w_{i, s} ^{-1} (\beta) = \sum_j n_j Z_j ^{-1} \exp( (\beta - \beta_j) E_{i, s} )

	The idea is that if we need to repeat jeckknife procedure chanfìging the block, we compute once the weights
	??? Should we compute just once also Z(\beta) = \sum_{i, s} w_{i, s} (or even better the stat of the whol sample) ???

*/