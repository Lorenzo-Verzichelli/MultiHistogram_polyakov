#pragma once

#include <iostream>

#include <vector>
#include <cmath>

void compute_zetas(int nrun, double* betas, std::vector<double>* energies, double* logZs);
/*Compute log(Z(beta)) for each SIMULATED beta
	Named "MultiHistRw" in Maio's source

	Results stored in logZs (can be zeros as passed)
	betas, energies and logZs must poin to preallocated memory of len nrun
	betas and energies will NOT be modified

	Recursively calls logZ_update until results are stable (Delta2 < 1.e-14)
	Newman-Barkema eq. 8.35
	\Delta^2 = \sum_k ( (Z' _k - Z_k) / (Z' _k) )^2
*/

double logZ_update(double beta, int nrun, double* betas, std::vector<double>* energies, double* logZs);
/*Updates the zetas for the iterative search, ancillary to compute_zetas
	Named "LnZ" in Maio's source

	same requests as for compute_zetas
	NO array passed will be changed

	Implements equation (8.36) from Newman-Barkema chap 8:
	Z_k = \sum_{i, s} ( \sum_j n_j Z_j ^{-1} \exp((\beta_k - \beta_j) E_{i, s}) )^{-1}
	i, j, k enumerate the simulated temperatures (0, ... , nrun-1)
	s		enumerates the samples simulated at a given temperature
*/

double log_obs_no_div_Z(std::vector<double>* obses, double power, double beta, int nrun, double* betas, std::vector<double>* energies, double* logZs);
/*Computes log(Z(\beta) obs)
	Named "LnO_n" in Maio's source

*/

double extrap_obs_pwr(std::vector<double>* obses, double power, double beta, int nrun, double* betas, std::vector<double>* energies, double* logZs);
/*Compute the expectation value of obs ^ power at beta
	NOTE that this version returns obs ^ power, already divided by Z(beta)

	pointers must point to nrun long array, wich will NOT be modified
	for each i = 0, ... nrun - 1 should hold:
		energies[i].size() == obses[i].size()
		NAN will be returned otherwise
*/