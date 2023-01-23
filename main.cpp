#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>

#include "MultiHistogram.h"
#include "Jackknife.h"

#define DATA_COL 4
/*
The input file must contain the name of each data file followed by the corresponding beta

data file should look like:
{Re Pi_s} {Re Pi_t} {Re L} {Im L}
(columns separated by ' ')
*/

#define ARG_INPUT 1
#define	ARG_OUTPUT 2
#define ARG_BETA_MIN 3
#define ARG_BETA_MAX 4
#define ARG_BETA_STEPS 5
#define ARG_BLOCK 6
#define ARG_SPACE 7
#define ARG_TIME 8
#define ARG_DIM 9

#define ARGS_NUM 10

int main(int argc, char* argv[]) {

	if (argc != ARGS_NUM) {
		std::cout << "This program computes chi_L throu resampling in beta (multi histogram)\n"
			<< "Compiled from file " << __FILE__ << "\n"
			<< "Authothor Lorenzo Verzichelli (lorenzo.verzichelli@gmail.com) \n";
#ifdef __INTEL_COMPILER
		std::cout << "Compiled with icc\n";
#elif defined( __GNUC__ )
		std::cout << "Compiled with gcc " << __GNUC__ << " " << __GNUC_MINOR__ << " " << " " << __GNUC_PATCHLEVEL__ << "\n";
#endif // __INTEL_COMPILER
		std::cout << "Usage: " << argv[0] << " input_file output_file beta_min beta_max beta_step_num block L T space_time_dim\n";
		exit(EXIT_FAILURE);
	}

	std::ifstream input(argv[ARG_INPUT]); //open input file for reading
	//count data file
	std::string data_fn;
	int nrun = 0;					//number of input data file (one for each beta)
	double beta;
	while (input >> data_fn) {
		nrun++;
		if (!(input >> beta)) {
			std::cerr << "Expected value for beta after " << data_fn << " in " << argv[ARG_INPUT] << "\n";
			exit(EXIT_FAILURE);
		}
	}

	if (nrun == 0) {
		std::cerr << "Problems in reading input (" << argv[ARG_INPUT] << ")\n";
		exit(EXIT_FAILURE);
	}

	double* betas = new double[nrun];
	std::vector<double>* energies = new std::vector<double>[nrun];
	std::vector<double>* polyakovs = new std::vector<double>[nrun];

	input.clear();
	input.seekg(0);

	double beta_min = atof(argv[ARG_BETA_MIN]);
	double beta_max = atof(argv[ARG_BETA_MAX]);
	int beta_step_num = atoi(argv[ARG_BETA_STEPS]);
	int block = atoi(argv[ARG_BLOCK]);
	int space_ext = atoi(argv[ARG_SPACE]);
	int time_ext = atoi(argv[ARG_TIME]);
	int st_dim = atoi(argv[ARG_DIM]);
	int quad_vol = time_ext;
	for (int dim = 1; dim < st_dim; dim++) {
		quad_vol *= space_ext;
	}
/*	
	std::cout << "Input file: " << argv[1] << "\n"
		<< "beta min: " << beta_min << "\n"
		<< "beta max: " << beta_max << "\n"
		<< "beta num: " << beta_step_num << "\n"
		<< "space ext: " << space_ext << "\n"
		<< "tim ext: " << time_ext << "\n"
		<< "space time dim: " << st_dim << "\n"; 
*/
	int t_plaq_num = (st_dim - 1) * quad_vol;
	int s_plaq_num = t_plaq_num * (st_dim - 2) / 2;

	std::ifstream data_in;
	double read[DATA_COL];	//data just read
	int read_row;	//row reading
	int read_col;	//col reading
	int crun = 0;
	while (input >> data_fn) {
		if (crun == nrun) {
			std::cerr << "Somthing strange happend in " << argv[ARG_INPUT] << ": expected " << nrun << " lines, but more were found \n";
			break;
		}
		if (!(input >> beta)) {
			std::cerr << "Expected value for beta after " << data_fn << " in " << argv[ARG_INPUT] << "\n";
			exit(EXIT_FAILURE);
		}
		betas[crun] = beta;

		data_in.open(data_fn);
		read_row = 0;
		while (data_in >> read[0]) {
			read_row++;

			for (read_col = 1; read_col < DATA_COL; read_col++) {
				if (!(data_in >> read[read_col])) {
					std::cerr << "Expected data at " << data_fn << " row " << read_row << " col " << read_col << "\n";
					exit(EXIT_FAILURE);
				}
			}

			energies[crun].push_back(- s_plaq_num * read[0] - t_plaq_num * read[1]);
			polyakovs[crun].push_back(sqrt(read[2] * read[2] + read[3] * read[3]));
		}

		data_in.close();
		crun++;

		if (read_row == 0) std::cerr << "No data found in " << data_fn << "\n";
	}
	if (crun < nrun)
		std::cerr << "Somthing strange happend in " << argv[ARG_INPUT] << ": expected " << nrun << " lines, but only " << crun << " were found\n";

	unsigned int total_stat = 0;
	for (int crun = 0; crun < nrun; crun++) {
		total_stat += (unsigned int) energies[crun].size();
	}
	if (total_stat == 0) {
		std::cerr << "No data was read! Unable to proceed\n";
		exit(EXIT_FAILURE);
	}

	input.close();

	std::cout << "beta \t #ener \t #poly \n";
	for (crun = 0; crun < nrun; crun++)
		std::cout << betas[crun] << "\t" << energies[crun].size() << "\t" << polyakovs[crun].size() << "\n"; 

	double* target_betas = new double[beta_step_num];
	beta_step_num--;
	target_betas[0] = beta_min;
	target_betas[beta_step_num] = beta_max;
	double delta_beta = (beta_max - beta_min) / beta_step_num;
	for (int j = 1; j < beta_step_num; j++)
		target_betas[j] = beta_min + delta_beta * j;

	std::time_t zeta_time = std::time(nullptr);
	std::cout << "Computing zetas..."; std::cout.flush(); 
	double* logZs = new double[nrun];
	compute_zetas(nrun, betas, energies, logZs);

	std::time_t obs_time = std::time(nullptr);
	std::cout << " done! in " << std::difftime(obs_time, zeta_time) <<" seconds\n"
		<< "Computing obsevables..."; std::cout.flush();

	double poly_mean, poly_susc;
/*
	for (int j = 0; j <= beta_step_num; j++) {
		poly_mean = extrap_obs_pwr(polyakovs, 1., target_betas[j], nrun, betas, energies, logZs);
		poly_susc = extrap_obs_pwr(polyakovs, 2., target_betas[j], nrun, betas, energies, logZs);
		poly_susc -= poly_mean * poly_mean;

		std::cout << std::fixed << std::setprecision(15) << target_betas[j] << "\t" << poly_mean << "\t" << poly_susc <<"\n";
	}
	
	double logZ_extr;
	std::cout << "a' la Maio" << std::endl;
	for (int j = 0; j <= beta_step_num; j++) {
		logZ_extr = logZ_update(target_betas[j], nrun, betas, energies, logZs);
    poly_mean = exp(log_obs_no_div_Z(polyakovs, 1., target_betas[j], nrun, betas, energies, logZs) - logZ_extr);
    poly_susc = exp(log_obs_no_div_Z(polyakovs, 2., target_betas[j], nrun, betas, energies, logZs) - logZ_extr);
    poly_susc -= poly_mean * poly_mean;

    std::cout << std::fixed << std::setprecision(15) << target_betas[j] << "\t" << poly_mean << "\t" << poly_susc <<"\n";
	}
*/
	std::ofstream output(argv[ARG_OUTPUT]);
	if (output.fail()) {
		std::cerr << "Unable to open output file" << std::endl;
		exit(EXIT_FAILURE);
	}

	double err_mean, err_susc;
	std::vector<double>* weights = new std::vector<double>[nrun];
	for (int j = 0; j <= beta_step_num; j++) {
		compute_weights(target_betas[j], nrun, betas, energies, logZs, weights);
	
		jackknife_mean_susc(polyakovs, nrun, weights, block, poly_mean, poly_susc, err_mean, err_susc);
		
		output << target_betas[j] << " "
			<< poly_mean << " " << err_mean << " "
			<< poly_susc << " " << err_susc << std::endl;
		
		for (int k = 0; k < nrun; k++) weights[k].clear();
	}
	std::cout << " done! in " << std::difftime(std::time(nullptr), obs_time) << " seconds" << std::endl;

	delete[] weights;
	delete[] logZs;
	delete[] target_betas;
	delete[] betas;
	delete[] polyakovs;
	delete[] energies;

	return 0;
}
