#include <ilcplex/ilocplex.h>
#include "Algorithms.hpp"
#include "test_alg.hpp"
#include <functional>
#include <vector>
#include <numeric>
#include <ctime>
#include <algorithm>

/*
Probabilistic Rounding.
x^: x_i = 1 with prob x^_i. x_i = 0 with prob 1- x^_i
if x is feasible, then take it as a solution
if not, find x's permutation and run LexAlg

Check packing as a preprosessor

Write solution into separate columns
*/

ILOSTLBEGIN

ILOUSERCUTCALLBACK2(add_lex_cut, LexAlg*, LexAlg_obj, IloBool*, new_sol_avail) {
	if (*new_sol_avail == IloTrue) {
		*new_sol_avail = IloFalse;
		IloRange new_constr(getEnv(),0,0);
		IloInt prev_ind;
		IloInt N = LexAlg_obj->var_array.getSize();		
		IloInt sub_ind;
		IloInt main_ind;
		IloNumArray phi_array(getEnv(), N);
		IloNum temp_num;
		for (IloInt _ind_var = N - 1; _ind_var >= 0; --_ind_var) {
			// Construct the phi array for current variable
			main_ind = LexAlg_obj->index_array[_ind_var];
			prev_ind = -1;
			for (IloInt _ind = _ind_var + 1; _ind < N; ++_ind) {
				sub_ind = LexAlg_obj->index_array[_ind];
				if (LexAlg_obj->sol_array[sub_ind] > 0) {
					if (prev_ind > 0) {
						phi_array[sub_ind] = phi_array[prev_ind] * (LexAlg_obj->var_array[prev_ind].getUB() - LexAlg_obj->sol_array[prev_ind] + 1);
					}
					else {
						prev_ind = sub_ind;
						phi_array[sub_ind] = LexAlg_obj->var_array[main_ind].getUB() - LexAlg_obj->sol_array[main_ind];
					}
				}
			}
			//check if the fractional solution violate the constriant
			temp_num = LexAlg_obj->callback_sol_arr[main_ind] + IloScalProd(phi_array, LexAlg_obj->callback_sol_arr) - IloScalProd(phi_array, LexAlg_obj->sol_array);
			if (temp_num > LexAlg_obj->sol_array[main_ind]) {
				new_constr.setExpr(LexAlg_obj->var_array[main_ind] + IloScalProd(phi_array, LexAlg_obj->var_array) - IloScalProd(phi_array, LexAlg_obj->sol_array));
				new_constr.setUB(LexAlg_obj->sol_array[main_ind]);
				add(new_constr, IloCplex::CutManagement::UseCutFilter);
			}
			//_range_array.add(LexAlg_obj->var_array[main_ind] + IloScalProd(phi_array, LexAlg_obj->var_array) - IloScalProd(phi_array, LexAlg_obj->sol_array) <= LexAlg_obj->sol_array[main_ind]);
		}
	}
}


ILOHEURISTICCALLBACK2(store_frac_sol, LexAlg*, LexAlg_obj, IloBool*, sol_received) {
	if (getNnodes() <= 0) {
		getValues(LexAlg_obj->callback_sol_arr, LexAlg_obj->var_array);
		*sol_received = IloTrue;
	}
}


ILOMIPINFOCALLBACK2(store_int_sol, LexAlg*, LexAlg_obj, IloBool*, sol_received) {
	if (getNnodes() <= 0 && hasIncumbent()) {
		getIncumbentValues(LexAlg_obj->callback_sol_arr, LexAlg_obj->var_array);
		*sol_received = IloTrue;
	}
}


ILOHEURISTICCALLBACK4(callback_get_info, LexAlg*, LexAlg_obj, IloNum*, rel_gap, IloNum*, best_result, IloBool*, new_sol_avail) {
	if (getNnodes() > 0) {
		*best_result = getIncumbentObjValue();
		*rel_gap = getMIPRelativeGap();
		//abort(); // stop CPLEX at root node
	} else {
		*best_result = getIncumbentObjValue();
		*rel_gap = getMIPRelativeGap();
		if (*new_sol_avail) {
			setSolution(LexAlg_obj->var_array, LexAlg_obj->sol_array);
			*new_sol_avail = IloFalse;
		}
	}
}

ILOINCUMBENTCALLBACK6(LexSolver_int, LexAlg*, LexAlg_obj, test_param*, parm_setting, IloInt*, external_num, IloNum*, rel_gap, IloNum*, best_result, IloBool*, new_sol_avail) {
	if (hasIncumbent() && (getNnodes() <= 0)) {
		getEnv().out() << "******FIND INTEGRAL SOLUTION********" << endl;
		// Apply rounding ONLY at the root node
		// get the heuristic solutions
		// solve it to better results
		*external_num += 1;
		if (*external_num % parm_setting->callback_freq == 0) {
			getEnv().out() << "freq: " << parm_setting->callback_freq << endl;
			getEnv().out() << "num: " << *external_num << endl;

			IloInt N = LexAlg_obj->var_array.getSize();
			try {
				getIncumbentValues(LexAlg_obj->callback_sol_arr, LexAlg_obj->var_array);  // Get fractional solution for LP relaxation
				// Sort the values in lpvals (not necessary cause we cannot sort variables with 0 values)
				// sort(LexAlg_obj->index_array.begin(), LexAlg_obj->index_array.end(), [&](IloInt a, IloInt b) {return parm_setting->callback_sol_arr[a] > parm_setting->callback_sol_arr[b]; });
				LexAlg_obj->solve();
				parm_setting->num_lex_alg_call += 1;
				IloNum lex_result = IloScalProd(LexAlg_obj->new_obj_coeff, LexAlg_obj->sol_array) + LexAlg_obj->obj_extra_const;
				parm_setting->result_list.push_back(lex_result);
				if (lex_result > parm_setting->best_lex_obj) {
					parm_setting->best_lex_obj = lex_result;
					*new_sol_avail = IloTrue;
				}
			}
			catch (...) {
				throw;
			}
		}
	}
}



ILOHEURISTICCALLBACK6(LexSolver, LexAlg*, LexAlg_obj, test_param*, parm_setting, IloInt*, external_num, IloNum*, rel_gap, IloNum*, best_result, IloBool*, new_sol_avail) {
	if (getNnodes() > 0) { 
		*best_result = getIncumbentObjValue();
		*rel_gap = getMIPRelativeGap();
		abort(); // stop CPLEX at root node
	}
	if (getNnodes() <= 0) {
		*external_num += 1;
		if (*external_num % parm_setting->callback_freq == 0) {
			getEnv().out() << "freq: " << parm_setting->callback_freq << endl;
			getEnv().out() << "num: " << *external_num << endl;

			// Apply rounding ONLY at the root node
			// Round the fractional node solution by first sorting it in increasing values
			// and then finding the lexmax solution in that permutation
			IloInt N = LexAlg_obj->var_array.getSize();
			//getEnv().out() << N << endl;
			//IloNumArray sol_array(getEnv(), N);
			try {
				getValues(LexAlg_obj->callback_sol_arr, LexAlg_obj->var_array);            // Get fractional solution for LP relaxation

				if (parm_setting->prob_rounding) {
					IloRandom rand_generator(getEnv(), 0);
					for (IloInt t = 0; t < N; t++) {
						if (rand_generator.getInt(100) > LexAlg_obj->callback_sol_arr[t] * 100) {
							LexAlg_obj->callback_sol_arr[t] = 0;
						}
						else {
							LexAlg_obj->callback_sol_arr[t] = 1;
						}
					}
				}

				IloBool run_lex = IloFalse;
				if (parm_setting->num_lex_alg_call > 1) {
					for (IloInt i = 0; i < N - 1; ++i) {
						if (LexAlg_obj->sol_array[LexAlg_obj->index_array[i]] > LexAlg_obj->sol_array[LexAlg_obj->index_array[i + 1]]) {
							run_lex = IloTrue;
						}
					}
				}
				else {
					run_lex = IloTrue;
				}
				if (run_lex) {
					// Sort the values in lpvals
					sort(LexAlg_obj->index_array.begin(), LexAlg_obj->index_array.end(), [&](IloInt a, IloInt b) {return LexAlg_obj->callback_sol_arr[a] > LexAlg_obj->callback_sol_arr[b]; });
					LexAlg_obj->solve();
					parm_setting->num_lex_alg_call += 1;
					IloNum lex_result = IloScalProd(LexAlg_obj->new_obj_coeff, LexAlg_obj->sol_array) + LexAlg_obj->obj_extra_const;
					parm_setting->result_list.push_back(lex_result);
					if (lex_result > parm_setting->best_lex_obj) {
						parm_setting->best_lex_obj = lex_result;
						setSolution(LexAlg_obj->var_array, LexAlg_obj->sol_array, lex_result);
					}
					*new_sol_avail = IloTrue;
				}
			}
			catch (...) {
				//sol_array.end();
				throw;
			}
		}
	}
}


ILOMIPINFOCALLBACK3(get_root_node_obj_val, IloNum*, root_node_objective_value, IloNum*, rel_gap, IloInt*, Nnode) {
	try {
		if (getNnodes() <= 0) {
			*root_node_objective_value = getIncumbentObjValue();
			*rel_gap = getMIPRelativeGap();
		}
		else {
			*Nnode = getNnodes();
		}
	}
	catch (...) {
		throw;
	}
}

/*
void test_stab(const char* output_file, _Algorithm_Choice choice, string _type, IloInt num_of_nodes, IloInt edge_prob, IloInt rand_seed) {
	try {
		
		IloEnv env;
		ofstream time_record;
		time_record.open(output_file, std::ofstream::out | std::ofstream::app);
		std::clock_t s_start;

		test_param parm_setting;
		parm_setting.callback_freq = 200;
		parm_setting.prob_rounding = IloFalse;

		LexAlg* algorithm = new CPLEX_solver(env, num_of_nodes, edge_prob, rand_seed);
		switch (choice) {
		case cplex:
			algorithm = new CPLEX_solver(env, num_of_nodes, edge_prob, rand_seed);
			break;
		case bisection:
			algorithm = new Bisection(env, num_of_nodes, edge_prob, rand_seed);
			break;
		case greedy:
			algorithm = new RoundingLP(env, num_of_nodes, edge_prob, rand_seed, 14);
			break;
		case rounding:
			algorithm = new RoundingLP(env, num_of_nodes, edge_prob, rand_seed, 14);
			//algorithm->model.add(IloConversion(env, algorithm->var_array, ILOFLOAT)); % The relaxed version is very slow to solve. % TODO: Alternatively, in branch-and-cut callback, add Lexi-constraint in each cut.
			break;
		case packing:
			algorithm = new RoundingLP(env, num_of_nodes, edge_prob, rand_seed, 14);
			algorithm->packing_matrix();
			//algorithm->run_MSE = IloFalse; //use proj method.
			algorithm->MSE_relaxed_flag = IloFalse; //Run IP MSE problem
			break;
		case prob_rounding:
			algorithm = new RoundingLP(env, num_of_nodes, edge_prob, rand_seed, 14);
			algorithm->model.add(IloConversion(env, algorithm->var_array, ILOFLOAT));
			break;
		case prob_packing:
			algorithm = new RoundingLP(env, num_of_nodes, edge_prob, rand_seed, 14);
			algorithm->packing_matrix();
			algorithm->model.add(IloConversion(env, algorithm->var_array, ILOFLOAT));
			break;
		}
		if (choice != cplex) {
			IloCplex cplex(env);
			IloNum best_result = -IloInfinity;
			IloObjective obj_func(env);
			IloModel model(env);
			obj_func = IloMaximize(env, IloScalProd(algorithm->new_obj_coeff, algorithm->var_array) + algorithm->obj_extra_const);
			model.add(algorithm->var_array);
			model.add(obj_func);
			model.add(algorithm->range_array);
			cplex.extract(model);
			IloInt callback_freq = 3;
			IloInt external_num=-1;
			IloNum rel_gap;
			IloBool new_sol_avail = IloFalse;
			env.out() << "Ready" << endl;
			if (choice == prob_rounding || choice == prob_packing) {
				parm_setting.prob_rounding = IloTrue;
				cplex.use(LexSolver_int(env, algorithm, &parm_setting, &external_num, &rel_gap, &best_result, &new_sol_avail));
				cplex.use(callback_get_info(env, algorithm, &rel_gap, &best_result, &new_sol_avail));
				cplex.use(add_lex_cut(env, algorithm, &new_sol_avail));
			}
			else {
				cplex.use(LexSolver_int(env, algorithm, &parm_setting, &external_num, &rel_gap, &best_result, &new_sol_avail));
				cplex.use(callback_get_info(env, algorithm, &rel_gap, &best_result, &new_sol_avail));
				cplex.use(add_lex_cut(env, algorithm, &new_sol_avail));
			}
			//cplex.setParam(IloCplex::Param::MIP::Limits::Nodes, 1);
			//cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);
			//cplex.setParam(IloCplex::Param::MIP::Strategy::RINSHeur, -1);
			cplex.setParam(IloCplex::Param::TimeLimit, 10800);
			cplex.setParam(IloCplex::Param::Threads, 1);
			algorithm->cplex.setParam(IloCplex::Param::Threads, 1);
			s_start = std::clock();
			cplex.solve();

			time_record << "Instance: " << "Stable set, " << "num of nodes: " << num_of_nodes << ", " << "edge prob: " << edge_prob << ", " << "rand seed: " << rand_seed << ", ";
			time_record << "ALG: " << _type << ", " << "Time: " << (std::clock() - s_start) / (double)CLOCKS_PER_SEC << " sec." << ", " << "Best OBJ: " << best_result;
			time_record << ", Permutations: " << parm_setting.num_lex_alg_call;
			time_record << ", " << "LP bound: " << cplex.getBestObjValue() << ", " << "Gap: " << rel_gap << endl;

			cplex.end();
		}
		else {
			//set up a new environment and cplex object.
			//Get the root node objective value by cplex default setting
			IloNum root_node_val = 0;
			IloNum rel_gap;
			algorithm->cplex.use(get_root_node_obj_val(env, &root_node_val, &rel_gap));
			//algorithm->cplex.setParam(IloCplex::Param::MIP::Limits::Nodes, 1);
			//algorithm->cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);
			algorithm->cplex.setParam(IloCplex::Param::TimeLimit, 10800);
			algorithm->cplex.setParam(IloCplex::Param::Threads,1);
			//algorithm->cplex.setParam(IloCplex::Param::MIP::Strategy::RINSHeur, -1);
			s_start = std::clock();
			algorithm->solve();
			time_record << "Instance: " << "Stable set, " << "num of nodes: "<< num_of_nodes <<", " << "edge prob: " << edge_prob << ", " << "rand seed: " << rand_seed << ", ";
			time_record << "ALG: " << _type << ", " << "Time: " << (std::clock() - s_start) / (double)CLOCKS_PER_SEC << " sec." << ", " << "Best OBJ: " << root_node_val;
			time_record << ", " << "LP bound: " << algorithm->cplex.getBestObjValue()  << ", " << "Gap: " << rel_gap << endl;
		}
		env.end();
		time_record.close();
	}
	catch (IloException& ex) {
		cerr << "Concert exception caught: " << ex << endl;
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}

	return;
}

**/

void test(const char* file_address, const char* output_file, _Algorithm_Choice choice, string _type, IloInt freq, double runtimelimit) {
	try {
		LexAlg* algorithm;
		IloEnv env;
		
		
		/*
		TEST PARAMETERS SETTING
		*/
		test_param parm_setting;
		parm_setting.callback_freq = freq;
		parm_setting.prob_rounding = IloFalse;
		//algorithm = new CPLEX_solver(env, file_address);

		switch (choice) {
		case cplex:
			algorithm = new CPLEX_solver(env, file_address);
			break;
		case cplexet:
			algorithm = new CPLEX_solver(env, file_address);
			break;
		case bisection:
			algorithm = new Bisection(env, file_address);
			break;
		case greedy:
			algorithm = new RoundingLP(env, file_address, 200);
			break;
		case rounding:
			algorithm = new RoundingLP(env, file_address, 200);
			//algorithm->model.add(IloConversion(env, algorithm->var_array, ILOFLOAT)); % The relaxed version is very slow to solve. % TODO: Alternatively, in branch-and-cut callback, add Lexi-constraint in each cut.
			break;
		case roundinglast:
			algorithm = new RoundingLP(env, file_address, 200);
			break;
		case packinglast:
			algorithm = new RoundingLP(env, file_address, 200);
			algorithm->packing_matrix();
			algorithm->MSE_relaxed_flag = IloFalse; //Run IP MSE problem
			break;
		case packing:
			algorithm = new RoundingLP(env, file_address, 200);
			//algorithm->cplex.exportModel("init.lp");
			algorithm->packing_matrix();
			//algorithm->cplex.exportModel("after.lp");
			//algorithm->run_MSE = IloFalse; //use proj method.
			algorithm->MSE_relaxed_flag = IloFalse; //Run IP MSE problem
			break;
		case prob_rounding:
			algorithm = new RoundingLP(env, file_address, 200);
			algorithm->model.add(IloConversion(env, algorithm->var_array, ILOFLOAT));
			break;
		case prob_packing:
			algorithm = new RoundingLP(env, file_address, 200);
			algorithm->packing_matrix();
			algorithm->model.add(IloConversion(env, algorithm->var_array, ILOFLOAT));
			break;
		case roundingend:
			algorithm = new RoundingLP(env, file_address, 200);
			break;
		case roundingendint:
			algorithm = new RoundingLP(env, file_address, 200);
			break;
		
		default:
			algorithm = new CPLEX_solver(env, file_address);
		}

		std::clock_t s_start;
		IloBool sol_received = IloFalse;

		if (choice != cplex && choice != cplexet) {
			IloCplex cplex(env);
			IloNum best_result = -IloInfinity;
			IloObjective obj_func(env);
			IloModel model(env);
			obj_func = IloMaximize(env, IloScalProd(algorithm->new_obj_coeff, algorithm->var_array) + algorithm->obj_extra_const);
			model.add(algorithm->var_array);
			model.add(obj_func);
			model.add(algorithm->range_array);
			cplex.extract(model);
			IloInt external_num = -1;
			IloNum rel_gap = IloInfinity;
			IloInt Nnode = 0;
			env.out() << "Ready" << endl;
			IloBool new_sol_avail = IloFalse;
			if (choice == prob_rounding || choice == prob_packing) {
				parm_setting.prob_rounding = IloTrue;
				algorithm->callback_sol_is_int = IloTrue;
				//cplex.use(LexSolver(env, algorithm, &parm_setting, &external_num, &rel_gap, &best_result, &new_sol_avail));
				//cplex.use(add_lex_cut(env, algorithm, &new_sol_avail));
				//cplex.use(get_root_node_obj_val(env, &best_result, &rel_gap));
				cplex.use(LexSolver_int(env, algorithm,  &parm_setting, &external_num, &rel_gap, &best_result, &new_sol_avail));
				cplex.use(callback_get_info(env, algorithm, &rel_gap, &best_result, &new_sol_avail));
			}
			else {

				if (choice == roundinglast || choice == packinglast){
					/*
						Callback: fractional solutions 
					*/
					env.out() << "Fractional Callback" << endl;
					cplex.use(LexSolver(env, algorithm, &parm_setting, &external_num, &rel_gap, &best_result, &new_sol_avail));
					cplex.use(add_lex_cut(env, algorithm, &new_sol_avail));
					cplex.use(get_root_node_obj_val(env, &best_result, &rel_gap, &Nnode));
				}
				if (choice == rounding || choice == packing) {
					/*
						Callback: incumbent solutions
					*/
					env.out() << "Integral Callback" << endl;
					algorithm->callback_sol_is_int = IloTrue;
					cplex.use(LexSolver_int(env, algorithm, &parm_setting, &external_num, &rel_gap, &best_result, &new_sol_avail));
					cplex.use(callback_get_info(env, algorithm, &rel_gap, &best_result, &new_sol_avail));
				}
				if (choice == roundingend) {
					/*
						Callback: store fractional solutions withouth performing any computations
					*/
					env.out() << "Fractional Callback" << endl;
					cplex.use(store_frac_sol(env, algorithm, &sol_received));
				}
				if (choice == roundingendint) {
					/*
						Callback: store incumbent solutions withouth performing any computations
					*/
					env.out() << "Incumbenet Callback" << endl;
					algorithm->callback_sol_is_int = IloTrue;
					cplex.use(store_int_sol(env, algorithm, &sol_received));
				}
			}
			cplex.setParam(IloCplex::Param::MIP::Strategy::Search, 1); // turn off dynamic search
			cplex.setParam(IloCplex::Param::MIP::Limits::Nodes, 1);
			if (choice == roundinglast || choice == packinglast) {
				cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);
				cplex.setParam(IloCplex::Param::MIP::Strategy::RINSHeur, -1);
			}
			cplex.setParam(IloCplex::Param::TimeLimit, runtimelimit);
			cplex.setParam(IloCplex::Param::Threads, 1);
			//algorithm->cplex.setParam(IloCplex::Param::Threads, 1);
			algorithm->cplex.setParam(IloCplex::Param::TimeLimit, runtimelimit);
			s_start = std::clock();
			cplex.solve();

			if (cplex.getStatus() == IloAlgorithm::Status::Feasible || cplex.getStatus() == IloAlgorithm::Status::Optimal) {
				if (best_result < cplex.getObjValue()) {
					best_result = cplex.getObjValue();
				}
				rel_gap = cplex.getMIPRelativeGap();
			}
			if (choice == roundingend || choice == roundingendint) {
				env.out() << "Rounding END Method" << endl;
				if (sol_received) {
					// perform tuning
					env.out() << "Rounding END Start" << endl;
					sort(algorithm->index_array.begin(), algorithm->index_array.end(), [&](IloInt a, IloInt b) {return algorithm->callback_sol_arr[a] > algorithm->callback_sol_arr[b]; });
					parm_setting.num_lex_alg_call += 1;
					algorithm->solve();
					parm_setting.best_lex_obj = IloScalProd(algorithm->new_obj_coeff, algorithm->sol_array) + algorithm->obj_extra_const;
				}
			}

			double time_gap = (double)(std::clock() - s_start) / CLOCKS_PER_SEC;
			if (choice == rounding) {
				_Algorithm_Choice cplexet_ag;
				cplexet_ag = cplexet;
				std::string alg_type = "CPLEXET";
				test(file_address, output_file, cplexet_ag, alg_type, 1, time_gap);
			}

			

			ofstream time_record;
			time_record.open(output_file, std::ofstream::out | std::ofstream::app);
			time_record << endl << "----------------------------------------------------------" << endl;
			time_record << "Instance: " << file_address << ", " << "if packing: ";
			if (algorithm->not_packing) { time_record << " NOT packing, "; }
			else { time_record << "IS packing, "; }
			time_record << "ALG: " << _type << ", " << "Time: " << time_gap << " sec." << ", " << "Best INT OBJ (CPLEX): " << best_result;
			//time_record << ", BEST INT (CPLEX RETURN):" << cplex.getObjValue();
			time_record << ", " << "BEST LEX OBJ: " << parm_setting.best_lex_obj << ", " << "TOTAL CALL LEX TIMES: " << parm_setting.num_lex_alg_call; 
			time_record << ", " << "LP bound: " << cplex.getBestObjValue() << ", " << "Gap: " << rel_gap;
			if (choice == packing || choice == packinglast) {
				time_record << ", " << "Proj_time: " << algorithm->mse_running_time;
				time_record << ", " << "h_0: " << algorithm->h_0;
			}
			time_record << endl;
			time_record << "------------------- List of lex max results ------------- " << endl;
			for (auto it = parm_setting.result_list.begin(); it < parm_setting.result_list.end(); ++it) {
				time_record << *it << ", ";
			}
			time_record << endl;
			
			if (choice == packing || choice == packinglast) {
				time_record << "------------------- List of packing Error results ------------- " << endl;
				for (auto it = algorithm->infeasible_packing_sol_array.begin(); it < algorithm->infeasible_packing_sol_array.end(); ++it) {
					time_record << *it << ", ";
				}
			}
			time_record << endl;
			time_record.close();
			cplex.end();
		}
		else {
			//set up a new environment and cplex object.
			//Get the root node objective value by cplex default setting
			IloNum root_node_val = 0;
			IloNum rel_gap = 1e75;
			IloInt Nnode = 0;
			algorithm->cplex.use(get_root_node_obj_val(env, &root_node_val, &rel_gap, &Nnode));
			if (choice != cplexet) {
				algorithm->cplex.setParam(IloCplex::Param::MIP::Limits::Nodes, 1);
			}
			algorithm->cplex.setParam(IloCplex::Param::MIP::Strategy::Search, 1);
			algorithm->cplex.setParam(IloCplex::Param::TimeLimit, runtimelimit);
			algorithm->cplex.setParam(IloCplex::Param::Threads,1);
			s_start = std::clock();
			algorithm->solve();
			double obj_value = -1e75;
			if (algorithm->cplex.getStatus() == IloAlgorithm::Status::Feasible || algorithm->cplex.getStatus() == IloAlgorithm::Status::Optimal) {
				obj_value = algorithm->cplex.getObjValue();
			}
			ofstream time_record;
			time_record.open(output_file, std::ofstream::out | std::ofstream::app);
			time_record << endl << "----------------------------------------------------------" << endl;
			time_record << "Instance: " << file_address << ", ";
			time_record << "ALG: " << _type << ", " << "Time: " << (std::clock() - s_start) / (double)CLOCKS_PER_SEC << " sec.";
			time_record << ", " << "LP bound: " << algorithm->cplex.getBestObjValue() << ", " << "Gap: " << rel_gap;
			time_record << ", " << "Best OBJ: " << obj_value << ", " << "Nodes: " << Nnode << endl;
			time_record.close();
		}
		env.end();
		
	}
	catch (IloException& ex) {
		cerr << "Concert exception caught: " << ex << endl;
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}

	return;
}


void test_QCP(const char* output_file, _Algorithm_Choice choice, string _type, IloInt rand_seed, IloBool binary_QCP, IloInt sparsity,
				IloInt num_quad_cons, IloInt num_linear_cons, IloInt num_var) {
	try {
		LexAlg* algorithm;
		IloEnv env;
		ofstream time_record;
		time_record.open(output_file, std::ofstream::out | std::ofstream::app);
		std::clock_t s_start;

		test_param parm_setting;
		parm_setting.callback_freq = 200;
		parm_setting.prob_rounding = IloFalse;
		algorithm = new CPLEX_solver(env, num_quad_cons, num_linear_cons, num_var, rand_seed, binary_QCP, sparsity);
		switch (choice) {
		case cplex:
			algorithm = new CPLEX_solver(env, num_quad_cons, num_linear_cons, num_var, rand_seed, binary_QCP, sparsity);
			break;
		case bisection:
			algorithm = new Bisection(env, num_quad_cons, num_linear_cons, num_var, rand_seed, binary_QCP, sparsity);
			break;
		case greedy:
			algorithm = new Greedy(env, num_quad_cons, num_linear_cons, num_var, rand_seed, binary_QCP, sparsity);
			break;
		case rounding:
			algorithm = new RoundingLP(env, num_quad_cons, num_linear_cons, num_var, rand_seed, binary_QCP, 200, sparsity);
			//algorithm->model.add(IloConversion(env, algorithm->var_array, ILOFLOAT)); % The relaxed version is very slow to solve. % TODO: Alternatively, in branch-and-cut callback, add Lexi-constraint in each cut.
			break;
		case packing:
			algorithm = new RoundingLP(env, num_quad_cons, num_linear_cons, num_var, rand_seed, binary_QCP, 200, sparsity);
			algorithm->packing_matrix();
			//algorithm->run_MSE = IloFalse; //use proj method.
			algorithm->MSE_relaxed_flag = IloFalse; //Run IP MSE problem
			break;
		case prob_rounding:
			algorithm = new RoundingLP(env, num_quad_cons, num_linear_cons, num_var, rand_seed, binary_QCP, 200, sparsity);
			algorithm->model.add(IloConversion(env, algorithm->var_array, ILOFLOAT));
			break;
		case prob_packing:
			algorithm = new RoundingLP(env, num_quad_cons, num_linear_cons, num_var, rand_seed, binary_QCP, 200, sparsity);
			algorithm->packing_matrix();
			algorithm->model.add(IloConversion(env, algorithm->var_array, ILOFLOAT));
			break;
		}
		if (choice != cplex) {
			IloCplex cplex(env);
			IloNum best_result = -IloInfinity;
			IloObjective obj_func(env);
			IloModel model(env);
			obj_func = IloMaximize(env, IloScalProd(algorithm->new_obj_coeff, algorithm->var_array) + algorithm->obj_extra_const);
			model.add(algorithm->var_array);
			model.add(obj_func);
			model.add(algorithm->range_array);
			cplex.extract(model);
			IloInt callback_freq = 3;
			IloInt external_num = -1;
			IloNum rel_gap;
			
			IloBool new_sol_avail = IloFalse;
			env.out() << "Ready" << endl;
			if (choice == prob_rounding || choice == prob_packing) {
				parm_setting.prob_rounding = IloTrue;
				cplex.use(LexSolver_int(env, algorithm, &parm_setting, &external_num, &rel_gap, &best_result, &new_sol_avail));
				cplex.use(callback_get_info(env, algorithm, &rel_gap, &best_result, &new_sol_avail));
				cplex.use(add_lex_cut(env, algorithm, &new_sol_avail));
			}
			else {
				algorithm->callback_sol_is_int = IloTrue;
				cplex.use(LexSolver_int(env, algorithm, &parm_setting, &external_num, &rel_gap, &best_result, &new_sol_avail));
				cplex.use(callback_get_info(env, algorithm, &rel_gap, &best_result, &new_sol_avail));
				//cplex.use(LexSolver(env, algorithm, &parm_setting, &external_num, &rel_gap, &best_result, &new_sol_avail));
				//cplex.use(get_root_node_obj_val(env, &best_result, &rel_gap));
			}
			cplex.setParam(IloCplex::Param::MIP::Limits::Nodes, 1);
			//cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);
			//cplex.setParam(IloCplex::Param::MIP::Strategy::RINSHeur, -1);
			cplex.setParam(IloCplex::Param::TimeLimit, 600);
			cplex.setParam(IloCplex::Param::Threads, 1);
			algorithm->cplex.setParam(IloCplex::Param::Threads, 1);
			s_start = std::clock();
			cplex.solve();

			time_record << "Instance: " << "QCP, " << "rand seed: " << rand_seed << ", " << "sparsity (expected): " << sparsity << ", ";
			time_record << "ALG: " << _type << ", " << "Time: " << (std::clock() - s_start) / (double)CLOCKS_PER_SEC << " sec." << ", " << "Best OBJ: " << best_result;
			time_record << ", Permutations: " << parm_setting.num_lex_alg_call;
			time_record << ", " << "LP bound: " << cplex.getBestObjValue() << ", " << "Gap: " << rel_gap << endl;

			cplex.end();
		}
		else {
			//set up a new environment and cplex object.
			//Get the root node objective value by cplex default setting
			IloNum root_node_val = 0;
			IloNum rel_gap;
			IloInt Nnode;
			algorithm->cplex.use(get_root_node_obj_val(env, &root_node_val, &rel_gap, &Nnode));
			algorithm->cplex.setParam(IloCplex::Param::MIP::Limits::Nodes, 1);
			//algorithm->cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);
			algorithm->cplex.setParam(IloCplex::Param::TimeLimit, 600);
			algorithm->cplex.setParam(IloCplex::Param::Threads, 1);
			//algorithm->cplex.setParam(IloCplex::Param::MIP::Strategy::RINSHeur, -1);
			s_start = std::clock();
			algorithm->solve();
			time_record << "Instance: " << "QCP, "<< "rand seed: " << rand_seed << ", " << "sparsity (expected): " << sparsity << ", ";
			time_record << "ALG: " << _type << ", " << "Time: " << (std::clock() - s_start) / (double)CLOCKS_PER_SEC << " sec." << ", " << "Best OBJ: " << root_node_val;
			time_record << ", " << "LP bound: " << algorithm->cplex.getBestObjValue() << ", " << "Gap: " << rel_gap << endl;
		}
		env.end();
		time_record.close();
	}
	catch (IloException& ex) {
		cerr << "Concert exception caught: " << ex << endl;
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}

	return;
}
