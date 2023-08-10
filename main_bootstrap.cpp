#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>
#include <cstdlib>

#include "MultiHistogram.h"
#include "Jackknife.h"

#define DATA_COL 4
/*
The input file must contain the name of each data file followed by the corresponding beta

data file should look like:
{Re Pi_s} {Re Pi_t} {Re L} {Im L}
(columns separated by ' ')
*/

#define ARG_INPUT		1
#define	ARG_BLOCK		2
#define ARG_BETA_MIN	3
#define ARG_BETA_MAX	4
#define ARG_BETA_STEPS	5
#define ARG_SPACE		6
#define ARG_TIME		7
#define ARG_DIM			8
#define ARG_BOOT_NUM    9

#define ARGS_NUM 		10

#define POLY_F_NAME "out_MH_poly_mean.dat"
#define SUSC_F_NAME "out_MH_poly_susc.dat"

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
		std::cout << "Usage: " << argv[0] << " input_file block_size beta_min beta_max beta_step_num block L T space_time_dim bootstrap_iterations\n"
            << "output in " POLY_F_NAME " and " SUSC_F_NAME << std::endl;
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
	std::vector<double>* energies_read = new std::vector<double>[nrun];
    std::vector<double>* energies = new std::vector<double>[nrun];
	std::vector<double>* polyakovs_read = new std::vector<double>[nrun];
	std::vector<double>* polyakovs = new std::vector<double>[nrun];

	input.clear();
	input.seekg(0);

	double beta_min = atof(argv[ARG_BETA_MIN]);
	double beta_max = atof(argv[ARG_BETA_MAX]);
    int block = atoi(argv[ARG_BLOCK]);
	int beta_step_num = atoi(argv[ARG_BETA_STEPS]);
	int space_ext = atoi(argv[ARG_SPACE]);
	int time_ext = atoi(argv[ARG_TIME]);
	int st_dim = atoi(argv[ARG_DIM]);
    int boot_num = atoi(argv[ARG_BOOT_NUM]);
	int quad_vol = time_ext;
	for (int dim = 1; dim < st_dim; dim++) {
		quad_vol *= space_ext;
	}

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

			energies_read[crun].push_back(- s_plaq_num * read[0] - t_plaq_num * read[1]);
			polyakovs_read[crun].push_back(sqrt(read[2] * read[2] + read[3] * read[3]));
		}

		data_in.close();
		crun++;

		if (read_row == 0) std::cerr << "No data found in " << data_fn << "\n";
	}
	if (crun < nrun)
		std::cerr << "Somthing strange happend in " << argv[ARG_INPUT] << ": expected " << nrun << " lines, but only " << crun << " were found\n";

	input.close();

    std::cout << "beta \t #ener \t #poly \n";
	for (crun = 0; crun < nrun; crun++)
		std::cout << betas[crun] << "\t" << energies[crun].size() << "\t" << polyakovs[crun].size() << "\n"; 

    srand((int) time(nullptr));

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
	compute_zetas(nrun, betas, energies_read, logZs);

	std::time_t obs_time = std::time(nullptr);
	std::cout << " done! in " << std::difftime(obs_time, zeta_time) <<" seconds\n"
		<< "Computing obsevables..."; std::cout.flush();

	double poly_mean, poly_susc;
    int num_blocks, start, tot_blocks;
    unsigned int samples, total_stat;

  	double err_mean, err_susc;

    std::vector<double>* weights = new std::vector<double>[nrun];
	std::ofstream output_poly(POLY_F_NAME);
    std::ofstream output_susc(SUSC_F_NAME);

    total_stat = 0;
    for (int crun = 0; crun < nrun; crun++) {
        samples = (unsigned int) energies_read[crun].size();
        num_blocks = samples / block;
        tot_blocks = num_blocks * block;
        energies[crun].reserve(tot_blocks);
        polyakovs[crun].reserve(tot_blocks);
        total_stat += tot_blocks;
    }

    if (total_stat == 0) {
        std::cerr << "No data was read! Unable to proceed\n";
		exit(EXIT_FAILURE);
    }

    for (int boot = 0; boot < boot_num; boot++) {
	    for (int crun = 0; crun < nrun; crun++) {
            samples = (unsigned int) energies_read[crun].size();
            num_blocks = samples / block;
            for (int j = 0; j < num_blocks; j++) {
                start = (rand() % num_blocks) * block;
                for (int k = 0; k < block; k++) {
                    energies[crun].push_back(energies_read[crun][start + k]);
                    polyakovs[crun].push_back(polyakovs_read[crun][start + k]);
                }
            }
        }

		for (int j = 0; j <= beta_step_num; j++) {
			compute_weights(target_betas[j], nrun, betas, energies, logZs, weights);
			jackknife_mean_susc(polyakovs, nrun, weights, block, poly_mean, poly_susc, err_mean, err_susc);
			output_poly << boot << target_betas[j] << " "
				<< poly_mean << " " << err_mean << "\n";
            output_susc << boot << target_betas[j] << " "
				<< poly_susc << " " << err_susc << "\n";
			
			for (int k = 0; k < nrun; k++) weights[k].clear();
		}

        for (int crun = 0; crun < nrun; crun ++) {
            energies[crun].clear();
            polyakovs[crun].clear();
        }

        std::cout << " iteration " << boot << " done! in " << std::difftime(std::time(nullptr), obs_time)
            << " seconds" << std::endl << "                      ";
    }
	
    output_poly.close();
    output_susc.close();

	delete[] weights;
	delete[] logZs;
	delete[] target_betas;
	delete[] betas;
	delete[] polyakovs;
	delete[] energies;
    delete[] polyakovs_read;
	delete[] energies_read;
    
	return 0;
}
