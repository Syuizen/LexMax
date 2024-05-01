#ifndef Algorithms_hpp
#define Algorithms_hpp

#include <ilcplex/ilocplex.h>
#include <vector>
#include <stdio.h>

ILOSTLBEGIN

enum _Problem_Type {Read_Only, Random_QC};

enum _Algorithm_Choice {cplex, bisection, greedy, rounding, packing, prob_rounding, prob_packing, roundinglast, packinglast, cplexet, roundingend, roundingendint};


class Quadratic_Constraint {
public:
	IloArray<IloNumArray> quad_coeff_mat;
	IloNumArray linear_coeff;
	IloNum ub;
	Quadratic_Constraint(IloEnv env, IloNum var_num) {
		quad_coeff_mat = IloArray<IloNumArray>(env, var_num);
		for (IloInt i = 0; i < var_num; i++) {
			quad_coeff_mat[i] = IloNumArray(env, var_num);
		}
		linear_coeff = IloNumArray(env, var_num);
	}
};

class LexAlg {
private:
	// IloInt cpu_limit = 20;
	// The following 3 variables are used to generate the stable set model based on a random graph
protected:
    const char* name; // the file containing model information
	IloEnv env; // invoking environment
	IloNumArray all_sol_array; // solution vector of all variables in the invoking model which may contain some auxiliary variables defined in some methods
	IloNumVarArray packing_var_array; // variables of the model after converting to packing type
	IloRangeArray packing_range_array; // constraints of the model after converting to packing type
	IloNumArray h_array; // h vector used to convert model to packing type
	IloNumArray var_change_ub; // the upper bound of all variables x_i needed to be replaced by u_i-y_i where u_i is the upper bound of x_i
	IloNumArray var_change_lb; // the upper bound of all variables x_i needed to be replaced by y_i-l_i where l_i is the lower bound of x_i
	IloNumArray var_change_indi; // indicator vector. If x_i is replaced by u_i-y_i, then var_change_indi[i] = -1; if x_i is replaced by y_i + l_i, then var_change_indi[i] = 1; Otherwise, it's 0
	IloBool packing_convertor = IloFalse; // indicate whether the convert_to_packing() has been called or not
	IloInt num_of_nodes; // total number of nodes in the graph
	IloInt edge_prob; // the probability of existance of each edge
	IloInt rand_seed; // the seed of random number generator
	IloBool generate_rand_graph = IloFalse; // indicate whether we should generate a random graph or not
	IloIntArray sign_of_cols; // indicate the sign of elements in each column. 1 if all elements are nonnegative, -1 if all elements are nonpositive, 0 if otherwise
	IloArray<IloIntArray> nonzero_coeff_index_mat; // store all indices of nonzero entry in constr_matrix
	IloInt num_quadratic_cons; // the number of quadratic type constraints
	IloInt num_linear_cons; // the number of linear type constraints
	IloInt num_var; // the number of variables
	IloInt sparsity_QCP; // the expected sparsity of QCP

public:
	const IloNum block_cut_off = 10000; // used to determine the size of the block
	IloInt h_0; // value used to convert the model to packing type
	IloBool linear_prob = IloTrue; // whether the invoking model is linear
	IloBool not_packing = IloFalse; //If packing algorithm outputs a new model then not_packing = TRUE 
	IloNum obj_extra_const; // equal to the const in the objective function
	IloNumVarArray var_array; // variables of the input model
	IloModel model; // invoking model
	IloBool MSE_relaxed_flag = IloTrue; // indicate if relax variables when call getMSEsolution()
	IloNumArray callback_sol_arr; // the solution returned from callback methods
	IloBool run_MSE = IloTrue; // indicate if call getMSEsolution() or proj()
	IloArray<IloNumArray> constr_matrix; // the matrix of coefficients of constraints in the invoking model
	IloArray<Quadratic_Constraint> quad_constr_arr; // the array of coefficients of quadratic constraints in the invoking model
	IloNumArray constr_matrix_bound; // the bound of linear constraints in the invoking model
	IloNumArray new_obj_coeff; // the coefficients of the objective function of the model after modifying 
	IloObjective obj_func; // the invoking objective function
	IloRangeArray range_array; // constraints of the model after modifying
	IloNumArray sol_array; // solution vector of the variables after replacement
	vector<IloInt> index_array; // permutation vector
	IloNumArray ori_coeff_array; // coefficients of the objective function of input model
	IloNumArray ori_sol_array; // solution vector of the variables in the input model
	IloCplex cplex; // invoking cplex object
	IloBool generate_QCP = IloFalse; // if generate a QCP instance randomly
	IloBool binary_QCP = IloFalse; // if generate a Binary QCP instance
	IloBool callback_sol_is_int = IloFalse; // if we use the integral solution found by heuristic method instead of fractional solution of relaxed model
	IloNum mse_running_time = 0; // store the total running time of inner cplex solving mse in packing (not the whole method)
	vector<IloNum> infeasible_packing_sol_array; // store how infeasible the packing solution is 
	void convert_to_packing();
	void check_packing();
	void define_var();
	LexAlg(IloEnv &_env, const char* _name) { name = _name; env = _env; define_var(); readproblem(); check_packing(); env.out() << "Successfully read problem." << endl;};
	LexAlg(IloEnv &_env, IloInt _num_of_nodes, IloInt _edge_prob, IloInt _rand_seed) {
		generate_rand_graph = IloTrue;
		num_of_nodes = _num_of_nodes; edge_prob = _edge_prob; rand_seed = _rand_seed;
		env = _env; define_var(); readproblem(); check_packing(); env.out() << "Successfully read problem." << endl; 
	};
	LexAlg(IloEnv &_env, IloInt _num_quad_cons, IloInt _num_linear_cons, IloInt _num_var, IloInt _rand_seed, IloBool _binary_QCP, IloInt _QCP_sparsity) {
		linear_prob = IloFalse;
		generate_QCP = IloTrue;
		binary_QCP = _binary_QCP;
		num_quadratic_cons = _num_quad_cons;
		num_linear_cons = _num_linear_cons;
		num_var = _num_var;
		rand_seed = _rand_seed;
		sparsity_QCP = _QCP_sparsity;
		env = _env; define_var(); readproblem(); env.out() << "Successfully generate QCP problem" << endl;
	}
	void random_QCP_generator(IloInt num_quadratic_cons, IloInt num_linear_cons, IloInt num_var, IloInt rand_seed, IloInt expected_sparsity);
    void readproblem();
    virtual void solve() = 0;
	void packing_matrix();	
	void nonzero_index_mat_generate();
	void get_ori_sol();
	void random_graph_generator(IloInt num_of_nodes, IloInt prob, IloInt seed);
};

class CPLEX_solver : public LexAlg {
public:
	void solve();
	CPLEX_solver(IloEnv &_env, const char* _name) : LexAlg(_env, _name) {};
	CPLEX_solver(IloEnv &_env, IloInt _num_of_nodes, IloInt _edge_prob, IloInt _rand_seed) : LexAlg(_env, _num_of_nodes,_edge_prob, _rand_seed) {};
	CPLEX_solver(IloEnv &_env, IloInt _num_quad_cons, IloInt _num_linear_cons, IloInt _num_var, IloInt _rand_seed, IloBool _binary_QCP, IloInt _QCP_sparsity) : LexAlg(_env, _num_quad_cons, _num_linear_cons, _num_var, _rand_seed, _binary_QCP, _QCP_sparsity) {};
};

class Bisection: public LexAlg {
public:
    void solve();
	Bisection(IloEnv &_env, const char* _name): LexAlg(_env, _name) {};
	Bisection(IloEnv &_env, IloInt _num_of_nodes, IloInt _edge_prob, IloInt _rand_seed) : LexAlg(_env, _num_of_nodes, _edge_prob, _rand_seed) {};
	Bisection(IloEnv &_env, IloInt _num_quad_cons, IloInt _num_linear_cons, IloInt _num_var, IloInt _rand_seed, IloBool _binary_QCP, IloInt _QCP_sparsity) : LexAlg(_env, _num_quad_cons, _num_linear_cons, _num_var, _rand_seed, _binary_QCP, _QCP_sparsity) {};
};

class Greedy: public LexAlg {
public:
    void solve();
	Greedy(IloEnv &_env, const char* _name) : LexAlg(_env, _name) {};
	Greedy(IloEnv &_env, IloInt _num_of_nodes, IloInt _edge_prob, IloInt _rand_seed): LexAlg(_env, _num_of_nodes, _edge_prob, _rand_seed) {};
	Greedy(IloEnv &_env, IloInt _num_quad_cons, IloInt _num_linear_cons, IloInt _num_var, IloInt _rand_seed, IloBool _binary_QCP, IloInt _QCP_sparsity) : LexAlg(_env, _num_quad_cons, _num_linear_cons, _num_var, _rand_seed, _binary_QCP, _QCP_sparsity) {};
};

class RoundingLP: public LexAlg {
public:
	IloInt size_of_block;
	IloInt max_size_of_block;
    void solve();
	RoundingLP(IloEnv &_env, const char* _name, IloInt block_size) : LexAlg(_env, _name) { max_size_of_block = block_size; 	RLP_preprocessor();};
	RoundingLP(IloEnv &_env, IloInt _num_of_nodes, IloInt _edge_prob, IloInt _rand_seed, IloInt block_size): LexAlg(_env, _num_of_nodes, _edge_prob, _rand_seed) { max_size_of_block = block_size; 	RLP_preprocessor(); };
	RoundingLP(IloEnv &_env, IloInt _num_quad_cons, IloInt _num_linear_cons, IloInt _num_var, IloInt _rand_seed, IloBool _binary_QCP, IloInt block_size, IloInt _QCP_sparsity) : LexAlg(_env, _num_quad_cons, _num_linear_cons, _num_var, _rand_seed, _binary_QCP, _QCP_sparsity) { max_size_of_block = block_size; };
	void clear_all_block_constrains();
	void RLP_preprocessor();
	void getMSEsolution(IloBool relaxed);
private:
	void proj();
	void proj_along(IloNumArray proj_direction);
	void packing_greedy_solver();
	IloRangeArray mse_range_array_lb; // constraints of the model to obtain MSE solution
	IloRangeArray mse_range_array_ub; // constraints of the model to obtain MSE solution
	IloNumVarArray mse_var_array; // extra varaibles in the MSE model
	IloInt num_of_var; // total number of variables in the invoking model
	IloInt block_ind; // the index of current working block. 
	IloNumArray pi; // the power series used in the block solver
	IloNumArray coeff_array; // the coefficients of objective function of the invoking model in the block solver
	IloArray<IloRangeArray> block_constr_array; // the new costraints added by the block solver
	IloRangeArray conv_lexi_array; // the inequalities defining the lexi order
	IloBool if_min; // indicate whether the invoking model is minimization or maximization
	IloInt block_start_ind; // the index of the first variable in the block solver
	IloInt block_end_ind; // the index of the variable after current block. In other words, current block only optimizes variables with indices from block_start_ind to block_end_ind (but not include block_end_ind)
	void block_model_solver(IloBool &if_optimal, IloRangeArray &block_constr);
	void conv_lexi_order_generator(IloRangeArray &_range_array, const IloNumArray &_h, IloNum ub_val);
	IloNum getQuadraticBestSolPerVar(IloNumArray &current_sol, IloInt var_ind, IloNum ub);
	IloNum getQuadraticBestSolPerVarPerCons(Quadratic_Constraint &quad_Cons, IloNumArray &current_sol, IloInt var_ind, IloNum ub);
};

#endif /* Algorithms_hpp */
