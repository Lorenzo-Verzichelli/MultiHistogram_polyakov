#include "MultiHistogram.h"

void compute_zetas(int nrun, double* betas, std::vector<double>* energies, double* logZs)
{
	double* new_logZs = new double[nrun];
	double Delta2 = 1.;	//distance from prior iteration

	double A;		//shift in Z
	double term;	//to update Delta2
	do {
		A = logZs[nrun / 2]; //not quite like newman-barkema (?)
		for (int k = 0; k < nrun; k++) logZs[k] -= A;

		for (int k = 0; k < nrun; k++) new_logZs[k] = logZ_update(betas[k], nrun, betas, energies, logZs);
		// \Delta^2 = \sum_k((Z' _k - Z_k) / (Z' _k)) ^ 2
		Delta2 = 0.;
		for (int k = 0; k < nrun; k++) {
			term = expm1(new_logZs[k] - logZs[k]);
			Delta2 += term * term;
		}
		
		for (int k = 0; k < nrun; k++)
			logZs[k] = new_logZs[k];
	} while (Delta2 > 1.e-14);

	delete[] new_logZs;
}

double logZ_update(double beta, int nrun, double* betas, std::vector<double>* energies, double* logZs)
/*
Formula: \sum_{i, s} ( \sum_j n_j Z_j ^{-1} \exp((\beta_k - \beta_j) E_{i, s}) )^{-1}
*/
{
	double res = 0.;
	double* logNs = new double[nrun];
	for (int i = 0; i < nrun; i++) {
		logNs[i] = log(static_cast<double>(energies[i].size()));
	}
	bool init = true;
	double E;			//will loop over measured energies
	double LogDen;		//will cumulate the log of the denominator
	double LogDen_term; //to update LogDen
	double diff;		//will contain difference of logs
	for (int i = 0; i < nrun; i++) {
		std::vector<double>::iterator cur = energies[i].begin();
		std::vector<double>::iterator last = energies[i].end();
		for (; cur != last; cur++) {
			E = *cur;
			LogDen = logNs[0] - logZs[0] + (beta - betas[0]) * E;
			for (int j = 1; j < nrun; j++) {
				LogDen_term = logNs[j] - logZs[j] + (beta - betas[j]) * E;
				diff = LogDen - LogDen_term; //avoid big exp
				if (diff > 0.) LogDen += log1p(exp(-diff));
				else LogDen = LogDen_term + log1p(exp(diff));
			}
			if (init) {
				res -= LogDen;
				init = false;
			}
			else {
				diff = res + LogDen;
				if (diff > 0) res += log1p(exp(-diff));
				else res = -LogDen + log1p(exp(diff));
			}
		}
	}
	delete[] logNs;
	return res;
}

double log_obs_no_div_Z(std::vector<double>* obses, double power, double beta, int nrun, double* betas, std::vector<double>* energies, double* logZs)
{
	//check input
	for (int i = 0; i < nrun; i++) {
		if (energies[i].size() != obses[i].size()) {
			std::cerr << "ERROR in extrap_obs_pwr: energies and obses for run " << i << " have incosistent size ("
				<< energies[i].size() << " and " << obses[i].size() << " respectively)\n";
			exit(EXIT_FAILURE);
		}
	}

	double res = 0.;
	double* logNs = new double[nrun];
	for (int k = 0; k < nrun; k ++) 
		logNs[k] = log(static_cast<double>(energies[k].size()));
	
	bool init = true;
	double E;			//will loop over energies
	double LogDen;		//will cumulate the log of the denominator
	double LogDen_term;	//to update LogDen
	double diff;		//difference of logs
	double LogO;		//will loop over log(obses ^ power)
	for (int i = 0; i < nrun; i++) {
		std::vector<double>::iterator cur = energies[i].begin();
		std::vector<double>::iterator last = energies[i].end();
		std::vector<double>::iterator curO = obses[i].begin();
		for (; cur != last; cur++, curO++) {
			E = *cur;
			LogDen = logNs[0] - logZs[0] + (beta - betas[0]) * E;
			for (int j = 1; j < nrun; j++) {
				LogDen_term = logNs[j] - logZs[j] + (beta - betas[j]) * E;
				diff = LogDen - LogDen_term;
				if (diff > 0) LogDen += log1p(exp(-diff));
				else LogDen = LogDen_term + log1p(exp(diff));
			}
			LogO = power * log((*curO));
			LogDen -= LogO;
			if (init) {
				res = -LogDen;
				init = false;
			}
			else {
				diff = res + LogDen;
				if (diff > 0) res += log1p(exp(-diff));
				else res = -LogDen + log1p(exp(diff));
			}
		}
	}
	delete[] logNs;
	return res;
}

double extrap_obs_pwr(std::vector<double>* obses, double power, double beta, int nrun, double* betas, std::vector<double>* energies, double* logZs)
{
	//check input
	for (int i = 0; i < nrun; i++) {
		if (energies[i].size() != obses[i].size()) {
			std::cerr << "ERROR in extrap_obs_pwr: energies and obses for run " << i << " have incosistent size ("
				<< energies[i].size() << " and " << obses[i].size() << " respectively)\n";
			exit(EXIT_FAILURE);
		}
	}

	double res = 0.;
	double* logNs = new double[nrun];
	for (int k = 0; k < nrun; k++)
		logNs[k] = log(static_cast<double>(energies[k].size()));

	bool init = true;
	double E;			//will loop over energies
	double LogDen;		//will cumulate the log of the denominator
	double LogDen_term;	//to update LogDen
	double diff;		//difference of logs
	double LogO;		//will loop over log(obses ^ power)
	double LogZ = 0;	//will cumulate the log of Z(beta)
	for (int i = 0; i < nrun; i++) {
		std::vector<double>::iterator cur = energies[i].begin();
		std::vector<double>::iterator last = energies[i].end();
		std::vector<double>::iterator curO = obses[i].begin();
		for (; cur != last; cur++, curO++) {
			E = *cur;
			LogDen = logNs[0] - logZs[0] + (beta - betas[0]) * E;
			for (int j = 1; j < nrun; j++) {
				LogDen_term = logNs[j] - logZs[j] + (beta - betas[j]) * E;
				diff = LogDen - LogDen_term;
				if (diff > 0) LogDen += log1p(exp(-diff));
				else LogDen = LogDen_term + log1p(exp(diff));
			}
			LogO = power * log((*curO));
			LogO -= LogDen;
			if (init) {
				LogZ = -LogDen;
				res = LogO;
				init = false;
			}
			else {
				diff = res - LogO;
				if (diff > 0) res += log1p(exp(-diff));
				else res = LogO + log1p(exp(diff));

				diff = LogZ + LogDen;
				if (diff > 0) LogZ += log1p(exp(diff));
				else LogZ = -LogDen + log1p(exp(diff));
			}
		}
	}
	delete[] logNs;
	return exp(res - LogZ);
}
