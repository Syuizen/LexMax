#include <ilcplex/ilocplex.h>
#include "Algorithms.hpp"
#include "test_alg.hpp"
#include <functional>
#include <vector>
#include <numeric>
#include <ctime>

int main(int argc, const char** argv) {
	/*
	Usage: LexAlg file_name, output_file_name, algorithm_choices (one of {cplex, bisection, greedy, rounding, packing})
	Example: LexAlg instance.mps output.txt cplex
	*/

	if (argc < 3) {
		std::cout << "Please speficify the input parameters" << std::endl;
		std::cout << "Parameters: file_name, output_file_name, algorithm_type" << std::endl;
		return 0;
	}

	_Algorithm_Choice choice;
	
	/*
	if (string(argv[1]) == "_STAB_" && argc == 4) {
		int num_of_nodes_arr[1] = {150};
		double prob_arr[6] = { 5, 30, 35, 40, 45, 50 };
		//int seed_arr[100]
		const char* output_file = argv[2];
		string _type = string(argv[3]);
		if (_type == "CPLEX") {
			choice = cplex;
			return 0;
		}
		else if (_type == "Bisection") {
			choice = bisection;
		}
		else if (_type == "Greedy") {
			choice = greedy;
		}
		else if (_type == "Rounding") {
			choice = rounding;
		}
		else if (_type == "Packing") {
			choice = packing;
		}
		else if (_type == "ProbRounding") {
			choice = prob_rounding;
		}
		else if (_type == "ProbPacking") {
			choice = prob_packing;
		}
		for (IloInt i = 0; i < 1; ++i) {
			IloInt num_of_nodes = num_of_nodes_arr[i];
			for (IloInt j = 0; j < 6; ++j) {
				IloInt edge_prob = prob_arr[j];	
				for (IloInt k = 1; k <= 100; ++k) {
					//IloInt rand_seed = seed_arr[k];
					IloInt rand_seed = k;
					test_stab(output_file, choice, _type, num_of_nodes, edge_prob, rand_seed);
				}
			}
		}
		return 0; 
	}
	*/

	if (string(argv[1]) == "_QCP_" && argc == 10) {
		std::cout << "start QCP test" << std::endl;
		const char* output_file = argv[2];
		string _type = string(argv[3]);
		IloInt seed = atoi(argv[4]);
		IloBool binary_QCP = atoi(argv[5]);
		IloInt sparsity = atoi(argv[6]);
		IloInt num_qc = atoi(argv[7]);
		IloInt num_linear = atoi(argv[8]);
		IloInt num_var = atoi(argv[9]);

		if (_type == "CPLEX") {
			choice = cplex;
		}
		else if (_type == "Bisection") {
			choice = bisection;
		}
		else if (_type == "Greedy") {
			choice = greedy;
		}
		else if (_type == "Rounding") {
			choice = rounding;
		}
		else if (_type == "Packing") {
			choice = packing;
		}
		else if (_type == "ProbRounding") {
			choice = prob_rounding;
		}
		else if (_type == "ProbPacking") {
			choice = prob_packing;
		}
		test_QCP(output_file, choice, _type, seed, binary_QCP, sparsity, num_qc, num_linear, num_var);
		return 0;
	}
	

	if (argc == 3) {
		const char* file_address = argv[1];
		const char* output_file = argv[2];
		cout << "TRY ALL METHODS" << endl;
		choice = cplex;
		test(file_address, output_file, choice, string("CPLEX"), 3);
		//choice = bisection;
		//test(file_address, output_file, choice, string("Bisection"), 3);
		//choice = greedy;
		//test(file_address, output_file, choice, string("Greedy"), 3);
		choice = rounding;
		test(file_address, output_file, choice, string("Rounding"), 5);
		choice = packing;
		test(file_address, output_file, choice, string("Packing"), 5);
		choice = roundingend;
		test(file_address, output_file, choice, string("RoundingEnd"), 10);
		choice = roundingendint;
		test(file_address, output_file, choice, string("RoundingEndInt"), 10);
		//choice = prob_rounding;
		//test(file_address, output_file, choice, string("ProbRounding"), 3);
		//choice = prob_packing;
		//test(file_address, output_file, choice, string("ProbPacking"), 3);
		choice = roundinglast;
		test(file_address, output_file, choice, string("RoundingLast"), 10);
		choice = packinglast;
		test(file_address, output_file, choice, string("PackingLast"), 10);
	}
	if (argc == 4) {
		const char* file_address = argv[1];
		const char* output_file = argv[2];
		//IloInt freq = atoi(argv[4]);
		string _type = string(argv[3]);
		if (_type == "CPLEX") {
			choice = cplex;
		}
		else if (_type == "Bisection") {
			choice = bisection;
		}
		else if (_type == "Greedy") {
			choice = greedy;
		}
		else if (_type == "Rounding") {
			choice = rounding;
		}
		else if (_type == "Packing") {
			choice = packing;
		}
		else if (_type == "ProbRounding") {
			choice = prob_rounding;
		}
		else if (_type == "ProbPacking") {
			choice = prob_packing;
		}
		else if (_type == "RoundingLast") {
			choice = roundinglast;
		}
		else if (_type == "PackingLast") {
			choice = packinglast;
		}
		else if (_type == "RoundingEnd") {
			choice = roundingend;
		}
		else if (_type == "RoundingEndInt") {
			choice = roundingendint;
		}
		if (choice != cplex) {
			test(file_address, output_file, choice, _type, 5);
		}
		else {
			test(file_address, output_file, choice, _type, 1);
		}
		
	}

	if (argc == 5) {
		const char* file_address = argv[1];
		const char* output_file = argv[2];
		string _type = string(argv[3]);
		IloInt freq = atoi(argv[4]);
		if (_type == "CPLEX") {
			choice = cplex;
		}
		else if (_type == "Bisection") {
			choice = bisection;
		}
		else if (_type == "Greedy") {
			choice = greedy;
		}
		else if (_type == "Rounding") {
			choice = rounding;
		}
		else if (_type == "Packing") {
			choice = packing;
		}
		else if (_type == "ProbRounding") {
			choice = prob_rounding;
		}
		else if (_type == "ProbPacking") {
			choice = prob_packing;
		}
		else if (_type == "RoundingLast") {
			choice = roundinglast;
		}
		else if (_type == "PackingLast") {
			choice = packinglast;
		}
		else if (_type == "RoundingEnd") {
			choice = roundingend;
		}
		else if (_type == "RoundingEndInt") {
			choice = roundingendint;
		}
		test(file_address, output_file, choice, _type, freq);
	}

    return 0;
}
