#include "Jackknife.h"

void jackknife_mean_susc(std::vector<double>* obses, int nrun, std::vector<double>* weights, unsigned int block, double& mean, double& variance, double& err_mean, double& err_variance)
{
	//check input
	for (int i = 0; i < nrun; i++) {
		if (weights[i].size() != obses[i].size()) {
			std::cerr << "ERROR in jackknife_mean_susc: weights and obses for run " << i << " have incosistent size ("
				<< weights[i].size() << " and " << obses[i].size() << " respectively)\n";
			exit(EXIT_FAILURE);
		}
	}

	int* block_nums = new int[nrun];
	int* start_blocks = new int[nrun];
	for (int j = 0; j < nrun; j++) {
		block_nums[j] = (int) obses[j].size() / block;
		start_blocks[j] = (int) obses[j].size() - block * block_nums[j];
	}

	weighted_mean stat; //statistics of the whole sample
	for (int i = 0; i < nrun; i++) {
		for (unsigned int s = start_blocks[i]; s < obses[i].size(); s++) {
			stat.add_log_weig(obses[i][s], weights[i][s]);
		}
	}

	int block_num_tot = 0;
	for (int j = 0; j < nrun; j++) block_num_tot += block_nums[j];
	double* jack_means = new double[block_num_tot];
	double* jack_vars = new double[block_num_tot];

	weighted_mean jack_stat;
	int block_pos, block_end;
	int block_count = 0;
	//loop over removed block
	for (int i = 0; i < nrun; i++) {
		for (int k = 0; k < block_nums[i]; k++) { //for each (i, k): one block
			jack_stat = stat; //set to statistics of the whole sample
			block_pos = k * block + start_blocks[i]; //where the block starts
			block_end = block_pos + block;
			for (int s = block_pos; s < block_end; s++) { //loop inside the block
				jack_stat.remove_log_weig(obses[i][s], weights[i][s]);
			}
			jack_means[block_count] = jack_stat.mean();
			jack_vars[block_count] = jack_stat.variance();

			block_count++;
		}
	}

	mean = 0;
	for (block_count = 0; block_count < block_num_tot; block_count++) mean += jack_means[block_count];
	mean /= block_num_tot;

	err_mean = 0;
	for (block_count = 0; block_count < block_num_tot; block_count++)
		err_mean += (jack_means[block_count] - mean) * (jack_means[block_count] - mean);
	err_mean *= block_num_tot - 1;
	err_mean /= block_num_tot;
	err_mean = sqrt(err_mean);

	variance = 0;
	for (block_count = 0; block_count < block_num_tot; block_count++) variance += jack_vars[block_count];
	variance /= block_num_tot;

	err_variance = 0;
	for (block_count = 0; block_count < block_num_tot; block_count++)
		err_variance += (jack_vars[block_count] - variance) * (jack_vars[block_count] - variance);
	err_variance *= block_num_tot - 1;
	err_variance /= block_num_tot;
	err_variance = sqrt(err_variance);

	delete[] jack_means;
	delete[] jack_vars;
	delete[] start_blocks;
	delete[] block_nums;
}

void compute_weights(double beta, int nrun, double* betas, std::vector<double>* energies, double* logZs, std::vector<double>* weights)
{
	//w_{ i, s }^ {-1} (\beta) = \sum_j n_j Z_j^ { -1 } \exp( (\beta - \beta_j) E_ { i, s } )
	for (int i = 0; i < nrun; i++) {
		weights[i].reserve(energies[i].size());
	}

	double* logNs = new double[nrun];
	for (int i = 0; i < nrun; i++) {
		logNs[i] = log(static_cast<double>(energies[i].size()));
	}

	double LogDen;
	double LogDen_term;
	double diff;
	for (int i = 0; i < nrun; i++) {
		for (double& E : energies[i]) {
			LogDen = logNs[0] - logZs[0] + (beta - betas[0]) * E;
			for (int j = 1; j < nrun; j++) {
				LogDen_term = logNs[j] - logZs[j] + (beta - betas[j]) * E;
				diff = LogDen - LogDen_term; //avoid big exp
				if (diff > 0.) LogDen += log1p(exp(-diff));
				else LogDen = LogDen_term + log1p(exp(diff));
			}
			weights[i].push_back(-LogDen);
		}
	}

	delete[] logNs;
}
