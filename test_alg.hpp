#include <ilcplex/ilocplex.h>
#include "Algorithms.hpp"
#include <functional>
#include <vector>
#include <numeric>
#include <ctime>

class test_param {
public:
	IloInt callback_freq = 200;
	IloNum best_lex_obj = -IloInfinity;
	IloInt num_lex_alg_call = 0;
	IloBool prob_rounding = IloFalse;
	vector<IloNum> result_list;
	IloBool new_heur_sol = IloFalse;
};


void test_stab(const char* output_file, _Algorithm_Choice choice, string _type, IloInt num_of_nodes, IloInt edge_prob, IloInt rand_seed);

void test(const char* file_address, const char* output_file, _Algorithm_Choice choice, string _type, IloInt freq, double runtimelimit = 3600);

void test_QCP(const char* output_file, _Algorithm_Choice choice, string _type, IloInt seed, IloBool binary_QCP, IloInt sparsity, IloInt num_quad_cons = 30, IloInt num_linear_cons = 1000, IloInt num_var = 200);