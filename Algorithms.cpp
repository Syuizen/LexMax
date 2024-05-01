#include "Algorithms.hpp"
#include <string>
#include <ilcplex/ilocplex.h>
#include "BigInteger/BigIntegerLibrary.hh"
#include <vector>


ILOMIPINFOCALLBACK2(get_root_node_sol, IloNumVarArray*, var_array, IloNumArray*, sol_array) {
	getEnv().out() << "--------------------------RUNNING MSE------------------" << endl;
	try {
		if (hasIncumbent()){ // if has feasible integer solution
			getIncumbentValues(*sol_array, *var_array);
			getEnv().out() << "--------------------------MSE SOLUTION------------------" << endl;
			getEnv().out() << *sol_array << endl;
			getEnv().out() << "--------------------------MSE SOLUTION------------------" << endl;
			abort(); //terminate current cplex solver process
		}
	}
	catch (...) {
		throw;
	}
}

void LexAlg::define_var() {
	cplex = IloCplex(env);
	model = IloModel(env);
	var_array = IloNumVarArray(env);
	range_array = IloRangeArray(env, 0);
	obj_func = IloObjective(env);
	cplex.extract(model);
}

void LexAlg::nonzero_index_mat_generate() {
	// create the nonzero_coeff_index_mat
	nonzero_coeff_index_mat = IloArray<IloIntArray>(env, var_array.getSize());
	// the row index corresponds to variables
	// Note that, in constr_matrix, the COLUMN index corresponds to variables
	for (IloInt i = 0; i < nonzero_coeff_index_mat.getSize(); ++i) {
		nonzero_coeff_index_mat[i] = IloIntArray(env, 0);
	}
	for (IloInt i = 0; i < constr_matrix.getSize(); ++i) {
		for (IloInt var_ind = 0; var_ind < var_array.getSize(); ++var_ind) {
			if (constr_matrix[i][var_ind] != 0) {
				nonzero_coeff_index_mat[var_ind].add(i);
			}
		}
	}
}

void LexAlg::random_graph_generator(IloInt num_of_nodes, IloInt prob, IloInt seed){
	env.out() << "start generating random graph" << endl;
	IloRandom rand_generator(env, seed);
	IloInt rand_num;
	var_array = IloNumVarArray(env, num_of_nodes, 0, 1, IloNumVar::Int);
	//model.add(var_array);
	for (IloInt i = 0; i < num_of_nodes; i++) {
		for (IloInt j = i + 1; j < num_of_nodes; j++) {
			rand_num = rand_generator.getInt(100);
			if (rand_num < prob) {
				range_array.add(var_array[i] + var_array[j] <= 1);
				model.add(range_array[range_array.getSize() -1]);
			}
		}
	}
	obj_func = IloMaximize(env, IloSum(var_array));
	model.add(obj_func);
}

void LexAlg::random_QCP_generator(IloInt _num_quad_cons, IloInt _num_linear_cons, IloInt _num_var, IloInt _rand_seed, IloInt _expected_sparsity) {
	quad_constr_arr = IloArray<Quadratic_Constraint>(env, num_quadratic_cons);
	
	IloRandom random_generator(env, _rand_seed);
	if (binary_QCP) {
		var_array = IloNumVarArray(env, _num_var, 0, 1, IloNumVar::Int);
	}
	else {
		var_array = IloNumVarArray(env, _num_var, 0, 100, IloNumVar::Int);
	}
	
	quad_constr_arr = IloArray<Quadratic_Constraint>(env, _num_quad_cons);
	constr_matrix = IloArray<IloNumArray>(env, _num_linear_cons);


	// generate quadratic constraints
	for (IloInt i = 0; i < _num_quad_cons; i++) {
		Quadratic_Constraint single_quad_constr(env, _num_var);
		for (IloInt j = 0; j < _num_var; j++) {
			for (IloInt z = 0; z < _num_var; z++) {
				if (random_generator.getInt(100) < _expected_sparsity) {
					single_quad_constr.quad_coeff_mat[j][z] = 1; // random_generator.getInt(2); // bug the densiy further reduced by 1/2
				}
				else {
					single_quad_constr.quad_coeff_mat[j][z] = 0;
				}
			}
			single_quad_constr.linear_coeff[j] = random_generator.getInt(2);
		}
		if (!binary_QCP) {
			// make quadratic matrix be diagonally-dominant
			for (IloInt j = 0; j < _num_var; j++) {
				for (IloInt z = 0; z < _num_var; z++) {
					if (j != z)
						single_quad_constr.quad_coeff_mat[j][j] += single_quad_constr.quad_coeff_mat[j][z];
				}
			}
		}

		single_quad_constr.ub = random_generator.getInt(1000);
		quad_constr_arr[i] = single_quad_constr;

		// add into the model
		// build a single quadratic constraint
		IloNumExprArray quad_constr(env, _num_var);
		for (IloInt j = 0; j < _num_var; j++) {
			quad_constr[j] = IloScalProd(single_quad_constr.quad_coeff_mat[j], var_array);
		}
		range_array.add(IloScalProd(quad_constr, var_array) + IloScalProd(single_quad_constr.linear_coeff, var_array) <= single_quad_constr.ub);
		model.add(range_array[range_array.getSize() - 1]);
	}

	
	// generate linear constraints
	constr_matrix_bound = IloNumArray(env, _num_linear_cons);
	constr_matrix = IloArray<IloNumArray>(env, _num_linear_cons);
	for (IloInt i = 0; i < _num_linear_cons; i++) {
		IloNumArray single_linear_constr(env, _num_var);
		IloNum rand_num;
		for (IloInt j = 0; j < _num_var; j++) {
			single_linear_constr[j] = random_generator.getInt(1); // bug here, this always give 0
		}
		constr_matrix[i] = single_linear_constr;
		rand_num = random_generator.getInt(100);
		constr_matrix_bound[i] = rand_num;
		// add into the model
		range_array.add(IloScalProd(single_linear_constr, var_array) <= rand_num);
		model.add(range_array[range_array.getSize() - 1]);
	}

	// generate obj coeff array
	ori_coeff_array = IloNumArray(env, _num_var);
	for (IloInt i = 0; i < _num_var; i++) {
		ori_coeff_array[i] = random_generator.getInt(100);
	}
	obj_func = IloMaximize(env, IloScalProd(ori_coeff_array, var_array));
	model.add(obj_func);
}

void LexAlg::readproblem() {
	if (generate_QCP) {
		random_QCP_generator(num_quadratic_cons, num_linear_cons, num_var, rand_seed, sparsity_QCP);
		obj_extra_const = 0;
		cplex.extract(model);
		new_obj_coeff = ori_coeff_array.copy();
		sol_array = IloNumArray(env, var_array.getSize());
		callback_sol_arr = IloNumArray(env, var_array.getSize());
		sign_of_cols = IloIntArray(env, var_array.getSize());
		for (IloInt i = 0; i < var_array.getSize(); i++) {
			index_array.push_back(i);
			sign_of_cols[i] = 1;
		}
		nonzero_index_mat_generate();
		all_sol_array = IloNumArray(env, var_array.getSize());
		return;
	}
	else {
		if (generate_rand_graph) {
			random_graph_generator(num_of_nodes, edge_prob, rand_seed);
			obj_extra_const = 0;
			env.out() << "random graph is generated" << endl;
			cplex.extract(model);
		}
		else {
			cplex.importModel(model, name, obj_func, var_array, range_array);
			cplex.extract(model);
		}
	}

	/*
	***************TEST SECTION*************** 
	*/

	sol_array = IloNumArray(env, var_array.getSize());
	var_change_ub = IloNumArray(env, var_array.getSize());
	var_change_lb = IloNumArray(env, var_array.getSize());
	var_change_indi = IloNumArray(env, var_array.getSize());
	callback_sol_arr = IloNumArray(env, var_array.getSize());
	for (IloInt i = 0; i < var_array.getSize(); i++) {
		index_array.push_back(i);
	}

	//check if any variable has nonzero lower bound
	for (IloInt i = 0; i < var_array.getSize(); ++i) {
		if (var_array[i].getLB() != 0) {
			var_change_indi[i] = 1;
			var_change_lb[i] = var_array[i].getLB();
		}
	}
	
	
	
	ori_coeff_array = IloNumArray(env, var_array.getSize());
	new_obj_coeff = IloNumArray(env, var_array.getSize());
	new_obj_coeff = ori_coeff_array.copy();
	for (IloExpr::LinearIterator it = IloExpr(obj_func.getExpr()).getLinearIterator(); it.ok(); ++it) {
		for (IloInt i = 0; i < var_array.getSize(); i++) {
			if (var_array[i].getId() == it.getVar().getId()) {
				ori_coeff_array[i] = it.getCoef();
				break;
			}
		}
	}
	


	if (obj_func.getSense() == IloObjective::Maximize) {
		//the model is maximization problem
		for (IloInt i = 0; i < ori_coeff_array.getSize(); i++) {
			if (ori_coeff_array[i] < 0) { 
				var_change_ub[i] = var_array[i].getUB();
				var_change_indi[i] = -1;
				var_change_lb[i] = 0; // note that, if we replace x_i by u_i - y_i, then y_i must have zero as its lower bound and hence we no longer need to store the lower bound
				new_obj_coeff[i] = -ori_coeff_array[i];
			}
			else {
				var_change_ub[i] = 0;
				new_obj_coeff[i] = ori_coeff_array[i];
			}
		}
	}
	else {
		//the model is minimization problem
		for (IloInt i = 0; i < ori_coeff_array.getSize(); i++) {
			if (ori_coeff_array[i] > 0) { 
				var_change_ub[i] = var_array[i].getUB();
				var_change_indi[i] = -1;
				var_change_lb[i] = 0; // note that, if we replace x_i by u_i - y_i, then y_i must have zero as its lower bound and hence we no longer need to store the lower bound
				new_obj_coeff[i] = ori_coeff_array[i];
			}
			else {
				var_change_ub[i] = 0;
				new_obj_coeff[i] = -ori_coeff_array[i];
			}
		}
	}
	obj_extra_const = IloScalProd(new_obj_coeff, var_change_lb) -IloScalProd(new_obj_coeff, var_change_ub)+obj_func.getConst();
	obj_func.setExpr( IloScalProd(new_obj_coeff, var_array) + obj_extra_const);
	obj_func.setSense(IloObjective::Maximize);
}

void LexAlg::check_packing() {
	/*
	This function converts the invoking model to the following form

	(P1)

		max		c1 x1 + ... + cn xn
		s.t.	A11 x1 + ... + A1n xn <= b1
				A12 x1 + ... + A2n xn <= b2
				...
				Am1 x1 + ... + Amn xn <= bm
				0 <= xi <= ui
	
	This function checks if constraint matrix of given model_file is packing or not.
	
	This function also initializes constr_matrix which store coefficients of each constraint.
	*/
	/***********************************
	ALERT: Assume all variables are nonnegative
	************************************/
	try {
		env.out() << "--------------start packing check----------------" << endl;

		constr_matrix = IloArray<IloNumArray>(env, range_array.getSize());
		constr_matrix_bound = IloNumArray(env, range_array.getSize());
		for (IloInt i = 0; i < constr_matrix.getSize(); i++) {
			constr_matrix[i] = IloNumArray(env, var_array.getSize());
		}
		IloBool finite_upper_bound_flag = IloFalse;
		IloBool finite_lower_bound_flag = IloFalse;
		IloBool packing_flag = IloTrue;
		IloRange temp_range;
		IloInt j = 0;
		IloNum temp_bound;

		// change bounds of all variables which have been replaced in readproblem()
		for (IloInt k = 0; k < var_array.getSize(); ++k) {
			if (var_change_indi[k] != 0) {
				// k-th var x_k is replaced by y_k = u_k - x_k
				// OR
				// k-th var x_k is replaced by y_k = x_k - l_k
				IloInt l_bound = var_array[k].getLB();
				IloInt u_bound = var_array[k].getUB();
				var_array[k].setUB(IloInfinity);
				var_array[k].setLB(0);
				var_array[k].setUB(u_bound - l_bound);
			}
		}
		IloInt num_of_range = range_array.getSize();
		for (IloInt i = 0; i < num_of_range; ++i) {
			finite_upper_bound_flag = IloFalse;
			finite_lower_bound_flag = IloFalse;
			temp_range = range_array[i];
			j = 0;
			for (IloExpr::LinearIterator it = temp_range.getLinearIterator(); it.ok(); ++it) {
				while (it.getVar().getId() != var_array[j].getId()) {
					constr_matrix[i][j] = 0;
					j = j + 1;
				}
				constr_matrix[i][j] = it.getCoef(); //Store the coefficients of the current constraint into constr_matrix
				j = j + 1;
			}
			// Update the i-th constraint if we have changed variables in the readproblem()
			for (IloInt k = 0; k < var_array.getSize(); ++k) {
				if (var_change_indi[k] == -1) {
					// k-th var x_k is replaced by y_k = u_k - x_k
					// Change coefficients of current constraint
					constr_matrix[i][k] = constr_matrix[i][k] * (-1); // a_k is replaced by -a_k
				}
			}

			// The next part are used to update the bound of constraints and break equality into two inequalities
			if (temp_range.getUB() < IloInfinity) {
				finite_upper_bound_flag = IloTrue;
			}
			if (temp_range.getLB() > -IloInfinity) {
				finite_lower_bound_flag = IloTrue;
			}

			/*
			Exactly one of the following 3 if-blocks will be run in each iteration

			1. b1 <= Ax <= b2 
			2. Ax <= b
			3. b <= Ax

			*/

			if (finite_upper_bound_flag == IloTrue && finite_lower_bound_flag == IloTrue) {
				// This case b1 <= Ax <= b2, then we need to split this constraint into two 
				temp_bound = temp_range.getLB() + IloScalProd(var_change_ub, constr_matrix[i]) - IloScalProd(var_change_lb, constr_matrix[i]); // new bound: b - SUM_{k: x_k is changed by u_k - y_k} a_k * u_k + SUM_{k: x_k is changed by l_k + y_k} a_k * l_k
				
				IloNumArray temp_constr_row(env, var_array.getSize());
				for (IloInt _temp_ind = 0; _temp_ind < var_array.getSize(); ++_temp_ind) {
					temp_constr_row[_temp_ind] = -constr_matrix[i][_temp_ind];
				}
				range_array.add(IloScalProd(temp_constr_row, var_array) <= -temp_bound); // b <= Ax ----> - Ax <= -b
				model.add(range_array[range_array.getSize() - 1]);
				constr_matrix.add(temp_constr_row);
				range_array[i].setLB(-IloInfinity);
				constr_matrix_bound.add(-temp_bound);
				// now change the Ax <= b2 part
				temp_bound = temp_range.getUB() + IloScalProd(var_change_ub, constr_matrix[i]) - IloScalProd(var_change_lb, constr_matrix[i]); // new bound: b - SUM_{k: x_k is changed by u_k - y_k} a_k * u_k + SUM_{k: x_k is changed by l_k + y_k} a_k * l_k
				range_array[i].setUB(temp_bound);
				range_array[i].setExpr(IloScalProd(constr_matrix[i], var_array));
				constr_matrix_bound[i] = temp_range.getUB();
			}
			if (finite_upper_bound_flag == IloTrue && finite_lower_bound_flag == IloFalse) {
				// Constriant type: Ax <= b ---> Ax <= b - SUM_{k: x_k is changed by u_k - x_k} a_k * u_k
				// This is the upper bound for current constraint after variable replacements
				// Note that, constr_matrix[i][k] = - a_k ( it has been changed ) if x_k is replaced
				temp_bound = temp_range.getUB() + IloScalProd(var_change_ub, constr_matrix[i]) - IloScalProd(var_change_lb, constr_matrix[i]); // new bound: b - SUM_{k: x_k is changed by u_k - y_k} a_k * u_k + SUM_{k: x_k is changed by l_k + y_k} a_k * l_k
				range_array[i].setExpr(IloScalProd(constr_matrix[i], var_array));
				range_array[i].setUB(temp_bound);
				constr_matrix_bound[i] = temp_range.getUB();
			}
			if (finite_upper_bound_flag == IloFalse && finite_lower_bound_flag == IloTrue) {
				// Constraint type: b <= Ax ----> b - SUM_{k: x_k is changed by u_k - x_k} a_k * u_k <= Ax
				// This is the lower bound for current constraint after variable replacements
				// Note that, constr_matrix[i][k] = - a_k ( it has been changed ) if x_k is replaced
				// Then convert b' <= Ax ----> -Ax <= -b'
				temp_bound = temp_range.getLB() + IloScalProd(var_change_ub, constr_matrix[i]) - IloScalProd(var_change_lb, constr_matrix[i]); // new bound: b - SUM_{k: x_k is changed by u_k - y_k} a_k * u_k + SUM_{k: x_k is changed by l_k + y_k} a_k * l_k
				for (IloInt _temp_ind = 0; _temp_ind < var_array.getSize(); ++_temp_ind) {
					constr_matrix[i][_temp_ind] = -constr_matrix[i][_temp_ind];
				}
				range_array[i].setExpr(IloScalProd(constr_matrix[i], var_array));
				range_array[i].setLB(-IloInfinity);
				range_array[i].setUB(-temp_bound);
				constr_matrix_bound[i] = temp_range.getUB();
			}

			// if we have a constraint like b <= Ax <= c, then current model cannot be packing
			if (finite_lower_bound_flag && finite_upper_bound_flag) { packing_flag = IloFalse; }
		}
		// Check if constraint matrix is packing (column by column)
		sign_of_cols = IloIntArray(env, var_array.getSize());
		IloBool positive_sign_flag = IloFalse;
		IloBool negative_sign_flag = IloFalse;
		if (packing_flag) {
			for (IloInt k = 0; k < var_array.getSize(); ++k) {
				positive_sign_flag = IloFalse;
				negative_sign_flag = IloFalse;
				for (IloInt i = 0; i < constr_matrix.getSize(); ++i) {
					if (constr_matrix[i][k] > 0) { positive_sign_flag = IloTrue; 
					}
					if (constr_matrix[i][k] < 0) { negative_sign_flag = IloTrue; 
					}
					if (positive_sign_flag && negative_sign_flag) {
						packing_flag = IloFalse;
					}
				}
				if (positive_sign_flag) { sign_of_cols[k] = 1; }
				if (negative_sign_flag) { sign_of_cols[k] = -1; }
				if (positive_sign_flag == IloTrue && negative_sign_flag == IloTrue) { sign_of_cols[k] = 9; }
			}
		}
		
		if (packing_flag) {
			not_packing = IloFalse;
			//string outputfilename =  name;
			//outputfilename += "_yes.lp";
			//cplex.exportModel(outputfilename.c_str());
			nonzero_index_mat_generate();
			env.out() << "---------Current model satifies packing property------------" << endl;
			return;
		}
		else {
			not_packing = IloTrue;
			//string outputfilename = name;
			//outputfilename += "_not.lp";
			//cplex.exportModel(outputfilename.c_str());

			env.out() << "---------Current model does NOT satisfy packing property------------" << endl;
		}
	}
	catch (...) {
		cplex.exportModel("error_model.lp");
		throw;
	}
}


void LexAlg::convert_to_packing() {
	/*
	This function converts the invoking model to the packing type
	Suppose the invoking model is the following
	
	(P1)
	
		max		c1 x1 + ... + cn xn
		s.t.	A11 x1 + ... + A1n xn <= b1
				A12 x1 + ... + A2n xn <= b2
				...
				Am1 x1 + ... + Amn xn <= bm
				0 <= xi <= ui

	Then the packing model is: (we introduce a new variable called s)

	(P2)

		max		c1 x1 + ... + cn xn
		s.t.	(A11 + h1) x1 + ... + (A1n + hn) xn + s <= b1 + h0
				(A12 + h1) x1 + ... + (A2n + hn) xn + s <= b2 + h0
				...
				(Am1 + h1) x1 + ... + (Amn + hn) xn + s <= bm + h0
				h1 x1 + h2 x2 + ... + hn xn + s <= h0  (****)
				0 <= xi <= ui
				0 <= s

	where
		hj = max { 0, max_i ( -Aij ) } = max { 0, - (min_i Aij) }, the negative of minimum entry of j-th column if the minimum entry is negative
		h0 = max { 0, max_i ( -bi ), OPTh } and OPTh is the optimal value of { max hx: Ax <= b, 0 <= x <= u }
	
	Note that (P1) and (P2) are NOT equivalent. We should project the solution of (P2) on the hyperplane h1 x1 + h2 x2 + ... + hn xn + s = h0 to get a feasible solution to (P1)
	*/
	
	if (not_packing == IloFalse) { return; } // if the invoking model is the packing model, then we are done
	// If the invoking model is not packing, then start converting

	packing_convertor = IloTrue;
	h_array = IloNumArray(env, var_array.getSize()); // To construct vector h, we need to find minimum entry of each column
	// Note that the initial values of h_array are 0
	h_0 = 0; 	// To construct h0, we need to find the minimum value of bound bi's

	// Find the minimum value of each column and bound
	for (IloInt ind_col = 0; ind_col < var_array.getSize(); ++ind_col) {
		for (IloInt ind_row = 0; ind_row < constr_matrix.getSize(); ++ind_row) {
			// If the constraint type is Ax <= b
			if (range_array[ind_row].getUB() < IloInfinity && constr_matrix[ind_row][ind_col] < h_array[ind_col]) {
				h_array[ind_col] = constr_matrix[ind_row][ind_col];
				if (h_0 > range_array[ind_row].getUB()) { h_0 = range_array[ind_row].getUB(); }
			}
			// If the constraint type is Ax => b, we should consider -Ax <= -b
			if (range_array[ind_row].getLB() > -IloInfinity && (-constr_matrix[ind_row][ind_col]) < h_array[ind_col]) {
				h_array[ind_col] = -constr_matrix[ind_row][ind_col];
				if (h_0 > (-range_array[ind_row].getLB())) { h_0 = -range_array[ind_row].getLB(); }
			}
		}
	}

	// h_array[j] = min{ 0, min_i A_{ij} }
	// We now construct hj = max { 0, - (min_i Aij) } = - h_array[j]
	for (IloInt _ind = 0; _ind < h_array.getSize(); ++_ind) {
		h_array[_ind] = -h_array[_ind];
	}
	// h0 = min_i b_i and we need to let h0 = max { 0, max_i ( -bi ), OPTh } and OPTh is the optimal value of { max hx: Ax <= b, 0 <= x <= u }
	h_0 = -h_0;
	
	// Find OPT h
		// Change the objective function from cx to hx
	obj_func.setExpr(IloScalProd(var_array, h_array));
	IloConversion temp_conver = IloConversion(env, var_array, ILOFLOAT);
	model.add(temp_conver);
	cplex.solve();
	IloNum z_h = cplex.getObjValue(); // z_h = OPTh
	h_0 = IloMax(h_0, IloCeil(z_h)); // h_0 = max { 0, max_i ( -bi ), OPTh }
	model.remove(temp_conver); // Remove relaxation
	temp_conver.end(); // Release memory

	/******************
		Start convert the model to (P2) form
	*******************/
	packing_var_array = IloNumVarArray(env, 0); // the array of all variables appearing in the packing model
	packing_var_array.add(var_array); // packing_var_array = var_array + s
	IloNumVar s(env, 0, h_0, ILOINT); // Create a new variable
	packing_var_array.add(s); // The last variable of packing_var_array is s
	packing_range_array = IloRangeArray(env, 0); // Store the new constraints appering when constructing packing model
	IloNum temp_ub;
	for (IloInt index_row = 0; index_row < constr_matrix.getSize(); index_row++) {
		
		/*
		Change constraints. There are mainly 3 types of constraints.
			1. Ax <= b -----> (A.+ h) x + s <= b + h0
			2. Ax => b -----> (-A. + h) x + s <= -b + h0
			3. b1 <= Ax <= b2 ------> two constraints: (A.+ h) x + s <= b2 + h0    and   (-A. + h) x + s <= -b1 + h0
		*/

		if (range_array[index_row].getUB() < IloInfinity) {
			temp_ub = range_array[index_row].getUB();
			for (IloInt index_col = 0; index_col < constr_matrix[index_row].getSize(); index_col++) {
				constr_matrix[index_row][index_col] += h_array[index_col];
			}
			constr_matrix_bound[index_row] = temp_ub + h_0;
			packing_range_array.add(IloScalProd(constr_matrix[index_row], var_array) + s <= constr_matrix_bound[index_row]);
			model.add(packing_range_array[packing_range_array.getSize() - 1]);
			//range_array[index_row].setExpr(IloScalProd(constr_matrix[index_row], var_array) + s <= constr_matrix_bound[index_row]);
		}
		if (range_array[index_row].getLB() > -IloInfinity) {
			IloNumArray new_constr_row(env, var_array.getSize());
			for (IloInt index_col = 0; index_col < constr_matrix[index_row].getSize(); index_col++) {
				if (range_array[index_row].getUB() < IloInfinity) {
					new_constr_row[index_col] = -constr_matrix[index_row][index_col] + 2 * h_array[index_col];
				}
				else {
					constr_matrix[index_row][index_col] = -constr_matrix[index_row][index_col] + h_array[index_col];
				}
			}
			if (range_array[index_row].getUB() < IloInfinity) {
				constr_matrix.add(new_constr_row);
				constr_matrix_bound.add(-range_array[index_row].getLB() + h_0);
				packing_range_array.add(IloScalProd(new_constr_row, var_array) + s <= -range_array[index_row].getLB() + h_0);
			}
			else {
				constr_matrix_bound[index_row] = -range_array[index_row].getLB() + h_0;
				packing_range_array.add(IloScalProd(constr_matrix[index_row], var_array) + s <= -range_array[index_row].getLB() + h_0);
			}
			model.add(packing_range_array[packing_range_array.getSize() - 1]);
		}
	}
	packing_range_array.add(IloScalProd(h_array, var_array) + s <= h_0); // Add the constraint hx + s <= h0
	model.add(packing_range_array[packing_range_array.getSize() - 1]);
	constr_matrix.add(h_array);
	constr_matrix_bound.add(h_0);
	// Set the objective function back to original
	obj_func.setExpr(IloScalProd(new_obj_coeff, var_array) + obj_extra_const);
	obj_func.setSense(IloObjective::Maximize);
	not_packing = IloFalse; // Label current model as packing type
	model.remove(range_array); // Remove the old constriants
	for (IloInt k = 0; k < var_array.getSize(); ++k) {
		sign_of_cols[k] = 1; // all entries in every column are non-negative
	}
}


void LexAlg::packing_matrix() {
	/*
	If matrix is not packing, then convert it to a matrix of packing type
	*/
	/***********************************
	ALERT: Assume all variables are nonnegative
	************************************/
	//check_packing();
	if (not_packing) {
		convert_to_packing();
		nonzero_index_mat_generate();
	}
}

void LexAlg::get_ori_sol() {
	// Get the solution to the original model (the model specified in the input file)
	// Before running this function, sol_array must have been assgined values
	ori_sol_array = IloNumArray(env, sol_array.getSize());
	for (IloInt _ind = 0; _ind < var_array.getSize(); ++_ind) {
		if (var_change_indi[_ind] != 0) {
			// current variable xi has been replaced by ui - yi
			ori_sol_array[_ind] = var_change_ub[_ind] - sol_array[_ind];
		}
		else {
			ori_sol_array[_ind] = sol_array[_ind];
		}
	}
}

void Bisection::solve() {
	IloInt num_of_var = var_array.getSize();
	IloNumArray ub_array(env, num_of_var);
	IloNumArray lb_array(env, num_of_var);
	IloRangeArray equality_array(env, 0);
	for (IloInt i = 0; i < var_array.getSize(); i++) {
		lb_array[i] = var_array[i].getLB();
	}
	for (IloInt i = 0; i < ub_array.getSize(); i++) {
		ub_array[i] = var_array[i].getUB();
	}

    IloInt start_index = var_array.getSize() - 1;
    IloBool loop_flag = IloTrue;
    IloBool func_flag = IloTrue;
    IloInt max_index;
    IloInt loop_upper;
    IloBool if_max_problem = IloTrue;
    
    if (obj_func.getSense() == IloObjective::Maximize) {
        if_max_problem = IloTrue;
    }
    else {
        if_max_problem = IloFalse;
    }
    
    IloNumArray coef_bisection(env, var_array.getSize());
    //coef_bisection[0] = 0;
    obj_func.setLinearCoefs(var_array, coef_bisection);
    try {
        cplex.solve();
    }
    catch (IloException& ex) {
        cerr << "Problem is infeasible" << endl;
        return;
    }
    
    while (func_flag) {
        env.out() << " LOOP START " << endl;
        loop_flag = IloTrue;
        while (loop_flag) {
            //cplex.exportModel("temp1.lp");
            for (IloInt i = start_index; i >= 0; i--) {
                if (lb_array[index_array[i]] != ub_array[index_array[i]]) {
                    loop_flag = IloFalse;
                    max_index = i;
                    break;
                }
            }
            start_index = max_index;
            
            env.out() << "update index: " << max_index << endl;
            
            if (loop_flag && start_index == 0) {
                //two solutions are identical
                loop_flag = IloFalse;
                func_flag = IloFalse;
            }
            
        }
        
        if (func_flag == IloFalse) {
            break;
        }
        
        loop_upper = IloCeil((lb_array[index_array[max_index]] + ub_array[index_array[max_index]]) / 2);
        if (if_max_problem) {
            cplex.addCut(var_array[index_array[max_index]] >= loop_upper);
            try {
                cplex.solve();
				cplex.clearCuts();
            }
            catch (IloException& ex) {
                //do nothing
            }
            if (cplex.getStatus() == IloAlgorithm::Infeasible) {
                //var_array[index_array[max_index]].setUB(loop_upper - 1);
				equality_array.add(var_array[index_array[max_index]] <= loop_upper - 1);
				model.add(equality_array[equality_array.getSize() - 1]);
				ub_array[index_array[max_index]] = loop_upper - 1;
			}
            else {
                //var_array[index_array[max_index]].setLB(loop_upper);
				equality_array.add(var_array[index_array[max_index]] >= loop_upper);
				model.add(equality_array[equality_array.getSize() - 1]);
				lb_array[index_array[max_index]] = loop_upper;
            }
        }
        else {
            cplex.addCut(var_array[index_array[max_index]] <= loop_upper - 1);
            try {
                cplex.solve();
				cplex.clearCuts();
            }
            catch (IloException& ex) {
                //do nothing
            }
            if (cplex.getStatus() == IloAlgorithm::Infeasible) {
                //var_array[index_array[max_index]].setLB(loop_upper);
				equality_array.add(var_array[index_array[max_index]] >= loop_upper);
				model.add(equality_array[equality_array.getSize()-1]);
				lb_array[index_array[max_index]] = loop_upper;
            }
            else {
                //var_array[index_array[max_index]].setUB(loop_upper - 1);
				equality_array.add(var_array[index_array[max_index]] <= loop_upper - 1);
				model.add(equality_array[equality_array.getSize() - 1]);
				ub_array[index_array[max_index]] = loop_upper - 1;
            }
        }
    }
    env.out() << "OKAY" << endl;
    cplex.solve();
    cplex.getValues(var_array, sol_array);
	model.remove(equality_array);
	return;
}

void Greedy::solve() {
	/*
	For a given permutation PER, 
	Try the greedy method by maximize each variable in the order PER(n), PER(n-1), ..., PER(1)
	Each iteration, we solve an IP problem
	*/

	try {
		env.out() << "Start Greedy Solver" << endl;
		IloNumArray cofarray(env, var_array.getSize());
		IloInt numofvar = var_array.getSize();;
		cofarray = IloNumArray(env, numofvar);
		for (IloInt i = 0; i < numofvar; i++) {
			cofarray[i] = 0;
		}
		cofarray[index_array[numofvar-1]] = 1;
		obj_func.setLinearCoefs(var_array, cofarray);
		cplex.solve();
		IloNum objval;
		objval = cplex.getObjValue();
		for (IloInt i = numofvar - 2; i > -1; i--) {
			cofarray[index_array[i]] = 1;
			cofarray[index_array[i+1]] = 0;
			obj_func.setLinearCoefs(var_array, cofarray);
			cplex.addCut(var_array[index_array[i+1]] == objval);
			cplex.solve();
			objval = cplex.getObjValue();
			env.out() << "current index: " << i << endl;
		}
		cplex.solve();
		cplex.getValues(var_array, sol_array);
		cplex.clearCuts();
	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	return;
}

void RoundingLP::block_model_solver(IloBool &if_optimal, IloRangeArray &block_constr) {
	try {
		//IloIntArray index_array(env, var_array.getSize());
		while (true) {
			// solve LP relaxation
			cplex.solve();

			if (cplex.getStatus() == IloAlgorithm::InfeasibleOrUnbounded || cplex.getStatus() == IloAlgorithm::Infeasible) {
				env.out() << "infeasible" << endl;
				if_optimal = IloFalse;
				return;
			}
			else {
				cplex.getValues(sol_array, var_array);
				IloNumArray block_sol_array(env, size_of_block);

				for (IloInt i = 0; i < size_of_block; ++i) {
					block_sol_array[i] = sol_array[index_array[i + block_start_ind]];
				}

				IloBool all_are_int = IloTrue;
				IloNum tol = IloPower(10, -5);
				for (IloInt i = block_start_ind; i < block_end_ind; i++) {
					if (sol_array[index_array[i]] - IloFloor(sol_array[index_array[i]]) > tol && IloCeil(sol_array[index_array[i]]) - sol_array[index_array[i]] > tol) {
						all_are_int = IloFalse;
						break;
					}
				}
				if (all_are_int == IloTrue) {
					env.out() << "block is optimal" << endl;
					for (IloInt i = block_end_ind - 1; i >= block_start_ind; i--) {
						block_constr.add(var_array[index_array[i]] == IloRound(sol_array[index_array[i]]));
						model.add(block_constr[block_constr.getSize() - 1]);
					}
					if_optimal = IloTrue;
					return;
				}
				else {
					model.remove(conv_lexi_array);
					conv_lexi_array.clear();
					conv_lexi_order_generator(conv_lexi_array, pi, IloScalProd(coeff_array, sol_array));
					model.add(conv_lexi_array);
				}
			}
		}
	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	return;
}

void RoundingLP::RLP_preprocessor() {
	mse_range_array_lb = IloRangeArray(env, 0);
	mse_range_array_ub = IloRangeArray(env, 0);
	mse_var_array = IloNumVarArray(env);
	num_of_var = var_array.getSize();
	pi = IloNumArray(env, max_size_of_block);
	coeff_array = IloNumArray(env, num_of_var);
	//block_constr_array = IloArray<IloRangeArray>(env, IloCeil((num_of_var) / 3) + 1);
	block_constr_array = IloArray<IloRangeArray>(env, 0);
	conv_lexi_array = IloRangeArray(env, 0);
	if_min = IloFalse;
	if (obj_func.getSense() == IloObjective::Minimize) {
		if_min = IloTrue;
	}

	//initialize block_constr_array
	/*
	for (IloInt i = 0; i < block_constr_array.getSize(); i++) {
		block_constr_array[i] = IloRangeArray(env, 0);
	}
	*/

	//initial coefficent is 0
	for (IloInt i = 0; i < num_of_var; i++) {
		coeff_array[i] = 0;
	}

	/*
	Let all variables with all-nonpositive-coefficents in constrains to attain their upperbound
	*/

	if (packing_convertor) {
		// If the packing convertor function has been called, then there is an auxiliary variable s in the invoking model
		all_sol_array = IloNumArray(env, packing_var_array.getSize());
	}
	else {
		// If the packing convertor function has not been called, then the invoking model contains exactly same variables as the input model
		all_sol_array = IloNumArray(env, var_array.getSize());
	}

	IloInt nonzero_row;
	// if packing_convertor is true, then coefficent matrix cannot have any negative entries
	if (not_packing == IloFalse) {
		for (IloInt var_ind = 0; var_ind < var_array.getSize(); ++var_ind) {
			if (sign_of_cols[var_ind] <= 0) {
				// For each non-positive column, let the associated variables attain their upper bounds
				all_sol_array[var_ind] = var_array[var_ind].getUB();
				for (IloInt _ind_2 = 0; _ind_2 < nonzero_coeff_index_mat[var_ind].getSize(); ++_ind_2) {
					// Then update the upper bound of each constraint
					nonzero_row = nonzero_coeff_index_mat[var_ind][_ind_2];
					constr_matrix_bound[nonzero_row] -= all_sol_array[var_ind] * constr_matrix[nonzero_row][var_ind];
				}
				cplex.addCut(var_array[var_ind] == var_array[var_ind].getUB());
				//add extra cut to enforce the variable to be its maximum
			}
		}
	}
	env.out() << "Preprossor finished" << endl;
}

void RoundingLP::solve() {
	try {
		//env.out() << "Current permutation: " << index_array << endl;
		// TODO: add the selection about whether use Greedy_solver() or not
		if (not_packing == IloFalse) {
			packing_greedy_solver();
			if (packing_convertor){
				if (run_MSE) {
					env.out() << "RUN MSE" << endl;
					proj_along(new_obj_coeff);
					//env.out() << sol_array << endl;
					getMSEsolution(MSE_relaxed_flag);
					//env.out() << "after MSE" << endl;
					//env.out() << sol_array << endl;
				}
				else {
					proj();
				}
			}
			return;
		}
		IloRangeArray lb_var_range_arr(env, 0);
		if (callback_sol_is_int) {
			//if the callback solution is an integral solution, then we check if we can increase any variables
			for (IloInt i = 0; i < var_array.getSize(); ++i) {
				//var_array[i].setLB(callback_sol_arr[i]);
				lb_var_range_arr.add(var_array[i] >= callback_sol_arr[i]);
			}
		}
		model.add(lb_var_range_arr);
		IloBool block_opt_flag;
		block_ind = 0;
		block_end_ind = var_array.getSize();
		size_of_block = 14;
		IloNum _temp = 1;
		while (_temp < block_cut_off && block_end_ind > size_of_block && size_of_block < max_size_of_block) {
			size_of_block += 1;
			_temp  *= (var_array[index_array[block_end_ind - size_of_block]].getUB() + 1);
		}
		block_start_ind = IloMax(block_end_ind - size_of_block, 0);
		block_constr_array.add(IloRangeArray(env, 0));
		while (block_start_ind >= 0) {
			env.out() << "current block start index: " << block_start_ind << endl;
			//_temp = 1;
			pi[0] = 1;
			for (IloInt i = block_start_ind; i < block_end_ind - 1; i++) {
				pi[i - block_start_ind + 1] = (var_array[index_array[i]].getUB() + 1)*pi[i - block_start_ind];
			}
			for (IloInt i = block_start_ind; i < block_end_ind; i++) {
				coeff_array[index_array[i]] = pi[i - block_start_ind];
			}
			obj_func.setLinearCoefs(var_array, coeff_array);
			block_opt_flag = IloTrue;

			env.out() << "start block solver" << endl;
			block_model_solver(block_opt_flag, block_constr_array[block_ind]);
			env.out() << "block solver finished" << endl;
			if (block_opt_flag == IloFalse) {
				/*
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				WARNING : THIS BLOCK METHOD IS NOT USED NOW
				TO MAKE IT WORK, ONE NEED TO ADD THE FEATURE TO DETERMINE size_of_block WHEN GO TO PREVIOUS BLOCK
				*/
				env.out() << "cannot optimize block" << endl;
				/*
				Backwards tracking method:
				Delete all cuts which have been added for this block.
				Then go to previous block.
				Cut previous optimal solution and resolve previous block.
				*/
				model.remove(block_constr_array[block_ind]);
				block_constr_array[block_ind].clear();
				if (block_ind == 0) {
					env.out() << "problem is infeasible" << endl;
					return;
				}
				//Go back to previous block
				env.out() << "go to previous block" << endl;
				block_ind = block_ind - 1;
				for (IloInt i = block_start_ind; i < block_end_ind; i++) {
					coeff_array[index_array[i]] = 0;
				}
				block_start_ind = block_end_ind;
				block_end_ind = block_start_ind + size_of_block;
				for (IloInt i = block_start_ind; i < block_end_ind; i++) {
					coeff_array[index_array[i]] = pi[i - block_start_ind];
				}
				model.remove(block_constr_array[block_ind]);
				block_constr_array[block_ind].clear();
				/*
				Backwards track new constraints
				*/
				model.remove(conv_lexi_array);
				conv_lexi_array.clear();
				conv_lexi_order_generator(conv_lexi_array, pi, IloScalProd(coeff_array, sol_array) - 1);
				model.add(conv_lexi_array);
				env.out() << "new model created" << endl;
			}
			else {
				for (IloInt i = block_start_ind; i < block_end_ind; i++) {
					coeff_array[index_array[i]] = 0;
				}
				
				block_end_ind = block_start_ind;
				if (block_end_ind == 0) {
					block_start_ind = -1;
				}
				else {
					//Determine the size of next block
					size_of_block = 0;
					_temp = 1;
					while (_temp < block_cut_off && block_end_ind > size_of_block && size_of_block < max_size_of_block) {
						size_of_block += 1;
						_temp *= (var_array[index_array[block_end_ind - size_of_block]].getUB() + 1);
					}
					block_start_ind = IloMax(block_end_ind - size_of_block, 0);
					block_constr_array.add(IloRangeArray(env, 0));
					block_ind += 1; //Go to next block
				}
			}
		}
		cplex.solve();
		cplex.getValues(var_array, sol_array);
		//-------------------------------------------
		clear_all_block_constrains(); 
		model.remove(lb_var_range_arr);
		//-------------------------------------------
		env.out() << "Lex solver finished" << endl;
	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	return;
}

void RoundingLP::conv_lexi_order_generator(IloRangeArray &_range_array, const IloNumArray &_h, IloNum ub_val) {
	/*
	This algorithm is used to generate linear constraints of lexi_order set Conv(X) where
	X = {x: hx < ub_val} where h is a power series. Index_array defines a lexicographical order.
	if if_min = False, then let X = {x: hx > ub_val}
	*/
	IloNumArray theta(env, _h.getSize());
	IloNumArray phi_array(env, _h.getSize());
	IloNumVarArray block_var_array(env, 0);

	for (IloInt i = 0; i < size_of_block; ++i) {
		block_var_array.add(var_array[index_array[i + block_start_ind]]);
	}

	IloNum temp_num;
	if (if_min) {
		//if model is minimization
		temp_num = IloCeil(ub_val);
	}
	else {
		//if model is maximization
		temp_num = IloFloor(ub_val);
	}
	for (IloInt i = size_of_block -1 ; i >= 0; --i) {
		theta[i] = IloMin(IloFloor((temp_num) / _h[i]), block_var_array[i].getUB());
		//env.out() << var_array[i] << " : theta: " << theta[i] << endl;
		temp_num -= theta[i] * _h[i];
	}

	// generate linear constraints for conv(X)
	// ASSUME the model is MAX
	IloInt prev_ind;
	for (IloInt _ind_var = size_of_block - 1; _ind_var >= 0; --_ind_var) {
		// Construct the phi array for current variable
		prev_ind = -1;
		for (IloInt _ind = _ind_var + 1; _ind < size_of_block; ++_ind) {
			if (theta[_ind] > 0) {
				if (prev_ind > 0) {
					phi_array[_ind] = phi_array[prev_ind] * (block_var_array[prev_ind].getUB() - theta[prev_ind] + 1);
				}
				else {
					prev_ind = _ind;
					phi_array[_ind] = block_var_array[_ind_var].getUB() - theta[_ind_var];
				}
			}
		}
	_range_array.add(block_var_array[_ind_var] + IloScalProd(phi_array, block_var_array) - IloScalProd(phi_array, theta) <= theta[_ind_var]);
	}
	theta.end();
	phi_array.end();
	return;
}


void RoundingLP::proj_along(IloNumArray proj_direction) {
	IloNum theta;
	IloNum top;
	IloNum bottom;

	env.out() << "--- projection ---" << endl;

	//env.out() << "--- call back solution ---" << endl;
	//env.out() << callback_sol_arr << endl;

	bottom = 0;
	top = h_0 - all_sol_array[all_sol_array.getSize() - 1];
	for (IloInt i = 0; i < h_array.getSize(); i++) {
		top -= h_array[i] * all_sol_array[i];
		bottom += h_array[i] * proj_direction[i];
	}
	
	if (bottom != 0) {
		// projected line does not parallel to the h_array
		theta = top / bottom;
		for (IloInt i = 0; i < h_array.getSize(); i++) {
			all_sol_array[i] += theta * proj_direction[i];
		}
	}
	//env.out() << "--- all_solution ---" << endl;
	//env.out() <<  all_sol_array << endl;

	for (IloInt _ind = 0; _ind < var_array.getSize(); _ind++) {
		sol_array[_ind] = all_sol_array[_ind];
	}
	//env.out() << "--- solution ---" << endl;
	//env.out() << sol_array << endl;
	
	return;
}



void RoundingLP::proj() {
	//TODO: handle arbitrary direction vector
	BigInteger biginteger_theta;
	BigInteger top((int)(h_0 - all_sol_array[all_sol_array.getSize() - 1])); //top =  h_0 - s
	BigInteger bottom(1);
	BigInteger _precision((long)1e15);
	BigInteger power_2(1);
	BigInteger num_2(2);
	vector<BigInteger> biginteger_h;
	vector<BigInteger> biginteger_x;
	BigInteger temp;
	IloNum theta;
	for (IloInt i = 0; i < h_array.getSize(); i++) {
		temp = (int)IloCeil(h_array[i]);
		biginteger_h.push_back(temp);
		temp = (int)IloCeil(all_sol_array[i]);
		biginteger_x.push_back(temp);
	}
	for (IloInt i = 0; i < h_array.getSize(); i++) {
		bottom += biginteger_h[i] * power_2;
		top -= biginteger_h[i] * biginteger_x[i];
		power_2 *= num_2;
	}
	top *= _precision;
	biginteger_theta = top / bottom;
	theta = biginteger_theta.toLong() / (1e15);
	all_sol_array[all_sol_array.getSize() - 1] -= theta;
	power_2 = power_2 / num_2;
	for (IloInt i = all_sol_array.getSize() - 2; i >= 0; i--) {
		all_sol_array[i] += (biginteger_theta*power_2).toLong() / (1e15);
		power_2 = power_2 / num_2;
	}
	for (IloInt _ind = 0; _ind < var_array.getSize(); _ind++) {
		sol_array[_ind] = all_sol_array[_ind];
	}
	return;

}

IloNum RoundingLP::getQuadraticBestSolPerVarPerCons(Quadratic_Constraint &quad_Cons, IloNumArray &current_sol, IloInt var_ind, IloNum ub){
	IloArray<IloNumArray> quad_coeff = quad_Cons.quad_coeff_mat;
	IloNumArray linear_coeff = quad_Cons.linear_coeff;
	IloNumArray quad_coeff_row(env, current_sol.getSize()); 
	IloNumArray quad_coeff_col(env, current_sol.getSize());
	for (IloInt i = 0; i < current_sol.getSize(); i++) {
		quad_coeff_row[i] = quad_coeff[var_ind][i];
		quad_coeff_col[i] = quad_coeff[i][var_ind];
	}
	IloNum b = IloScalProd(quad_coeff_col, current_sol) + IloScalProd(quad_coeff_row, current_sol) + linear_coeff[var_ind];
	IloNum d = IloScalProd(linear_coeff, current_sol);
	IloNumArray temp_sum_arr(env, current_sol.getSize());
	for (IloInt i = 0; i < current_sol.getSize(); i++) {
		temp_sum_arr[i] = IloScalProd(quad_coeff[i], current_sol);
	}
	d += IloScalProd(temp_sum_arr, current_sol);
	IloNum a = quad_coeff[var_ind][var_ind];
	IloNum theta;
	if (a > 0) {
		theta = IloMin(ub, IloFloor((b / 2 + sqrt(quad_Cons.ub + pow(b, 2) / 4 - d)) / sqrt(a)));
	}
	if (a == 0 && b > 0) {
		theta = IloMin(ub, IloFloor((quad_Cons.ub - d) / b));
	}
	else {
		theta = ub;
	}
	return theta;
}

IloNum RoundingLP::getQuadraticBestSolPerVar(IloNumArray &current_sol, IloInt var_ind, IloNum ub) {
	IloNum theta = ub;
	for (IloInt j = 0; j < quad_constr_arr.getSize(); j++) {
		theta = IloMin(theta, getQuadraticBestSolPerVarPerCons(quad_constr_arr[j], current_sol, var_ind, ub));
	}
	return theta;
}



void RoundingLP::packing_greedy_solver() {
	//Assume our problem is maximization
	IloNumArray temp_constr_matrix_bound = constr_matrix_bound.copy();

	env.out() << "start packing greedy solver" << endl;

	// Initial values of all_sol_array are 0

	/*
	We optimize the current model by the following strategy:

		1. Find all non-positive columns and let the corresponding variable to be its upper bound
		(This part has been moved to the initialized part, since objective coefficients are all 
		non-negative, no matter under which permutation lex max point must force these varibles to attain
		their upper bound)
		2. Greedily solve the model for variables corresponding to the positive columns

	*/
	IloNumArray quad_best_sol(env, var_array.getSize());
	if (linear_prob != IloTrue) {
		// the invoking model contains quadratic constraints
		for (IloInt i = var_array.getSize() - 1; i>0; i--) {
			IloInt var_ind = index_array[i];
			quad_best_sol[i] = getQuadraticBestSolPerVar(quad_best_sol, var_ind, var_array[var_ind].getUB());
		}		
	}

	if (callback_sol_is_int == IloTrue) {
		// If we get an integral solution from the callback function
		// this time we compute the increament vector 
		for (IloInt i = 0; i < var_array.getSize(); ++i) {
			if (sign_of_cols[i] <= 0) {
				callback_sol_arr[i] = 0;
			}
		}

		for (IloInt i = 0; i < constr_matrix_bound.getSize(); ++i) {
			temp_constr_matrix_bound[i] -= IloScalProd(constr_matrix[i], callback_sol_arr);
		}
	}

	IloInt var_ind;
	IloInt index_row;
	IloNum temp_val;
	for (IloInt _ind = var_array.getSize()-1; _ind >= 0; -- _ind) {
		var_ind = index_array[_ind];
		if (sign_of_cols[var_ind] == 1) {
			// constraint coefficients are positive for this variable
			// Let current variable to be as large as possible
			if (callback_sol_is_int == IloTrue) {
				all_sol_array[var_ind] = var_array[var_ind].getUB() - callback_sol_arr[var_ind];
				//if (all_sol_array[var_ind] == 0) {
				//	continue;
				//}
			}
			else {
				all_sol_array[var_ind] = var_array[var_ind].getUB();
			}
			for (IloInt i = 0; i < nonzero_coeff_index_mat[var_ind].getSize(); ++i) {
				index_row = nonzero_coeff_index_mat[var_ind][i];
				temp_val = temp_constr_matrix_bound[index_row] / constr_matrix[index_row][var_ind];
				if (all_sol_array[var_ind] > temp_val){
					all_sol_array[var_ind] = IloFloor(temp_val);
				}
			}
			// update the upper bound of each constraint
			for (IloInt i = 0; i < nonzero_coeff_index_mat[var_ind].getSize(); ++i) {
				index_row = nonzero_coeff_index_mat[var_ind][i];
				temp_constr_matrix_bound[index_row] -= all_sol_array[var_ind] * constr_matrix[index_row][var_ind];
			}
		}
	}

	if (callback_sol_is_int) {
		for (IloInt i = 0; i < var_array.getSize(); ++i) {
			all_sol_array[i] += callback_sol_arr[i];
		}
	}

	if (packing_convertor) {
		all_sol_array[all_sol_array.getSize() - 1] = IloMin(temp_constr_matrix_bound);
		for (IloInt _ind=0; _ind < var_array.getSize(); ++_ind){
			sol_array[_ind] = all_sol_array[_ind];
		}
	} else {
		for (IloInt _ind = 0; _ind < var_array.getSize(); _ind++) {
			sol_array[_ind] = all_sol_array[_ind];
		}
	}
	if (linear_prob == IloFalse) {
		for (IloInt _ind = 0; _ind < var_array.getSize(); _ind++) {
			sol_array[_ind] = IloMin(quad_best_sol[_ind], sol_array[_ind]);
		}
	}

	env.out() << "-- after packing greedy solver ---- " << endl;
	env.out() << all_sol_array << endl;

	return;
}

void RoundingLP::clear_all_block_constrains() {
	//Should call this function only after running RoundingLP.solve()
	for (IloInt i = 0; i < block_constr_array.getSize(); i++) {
		model.remove(block_constr_array[i]);
		block_constr_array[i].clear();
	}
	model.remove(conv_lexi_array);
}

void RoundingLP::getMSEsolution(IloBool relaxed) {
	/*
	Projection method: project the optimal solution to (P2) to a feasible solution to (P1) 
	In this function, the projected point is derived by finding the minimizer of square distance.
	*/
	try {
		if (mse_var_array.getSize() == 0) {
			for (IloInt i = 0; i < var_array.getSize(); i++) {
				IloNumVar abs_var(env);
				mse_var_array.add(abs_var);
			}
			model.add(mse_var_array);
		}

		if (mse_range_array_lb.getSize() == 0) {
			for (IloInt i = 0; i < var_array.getSize(); i++) {
				IloRange upper_range(env, var_array[i] - mse_var_array[i], sol_array[i]);
				mse_range_array_ub.add(upper_range);
				IloRange lower_range(env, sol_array[i], mse_var_array[i] + var_array[i]);
				mse_range_array_lb.add(lower_range);
			}
			model.add(mse_range_array_lb);
			model.add(mse_range_array_ub);
		}
		else {
			for (IloInt i = 0; i < var_array.getSize(); i++) {
				mse_range_array_lb[i].setLb(sol_array[i]);
				mse_range_array_ub[i].setUb(sol_array[i]);
			}
		}

		obj_func.setExpr(IloSum(mse_var_array));
		// Minimize L1 norm
		cplex.getEnv().out() << "Start MSE" << endl;


		// add range

		//obj_func.setExpr(var_array.getSize()*IloAbs(IloScalProd(var_array,sol_array)- IloSum(sol_array))-IloSum(var_array));
		obj_func.setSense(IloObjective::Minimize);

		packing_range_array[packing_range_array.getSize() - 1].setLB(h_0);
		infeasible_packing_sol_array.push_back(IloScalProd(h_array, sol_array));
		cplex.extract(model);

		cplex.setParam(IloCplex::Param::Emphasis::MIP, 1); // emphsize on the feasibility
		
		for (IloInt i = 0; i < var_array.getSize(); ++i) {
			sol_array[i] = -IloInfinity; // Initialize the solution array as negative infinity.
		}

		if (relaxed) {
			IloConversion temp_cov = IloConversion(env, packing_var_array, ILOFLOAT);
			model.add(temp_cov);
			std::clock_t s_start = std::clock();
			cplex.solve();
			mse_running_time += (double)(std::clock() - s_start) / CLOCKS_PER_SEC;
			cplex.getValues(var_array, sol_array);
			model.remove(temp_cov);
			temp_cov.end();
		}
		else {
			cplex.use(get_root_node_sol(env, &var_array, &sol_array));
			std::clock_t s_start = std::clock();
			cplex.solve();
			mse_running_time += (double)(std::clock() - s_start) / CLOCKS_PER_SEC;
			if (cplex.getStatus() == IloAlgorithm::Status::Feasible || cplex.getStatus() == IloAlgorithm::Status::Optimal) {
				cplex.getValues(var_array, sol_array);
			}
			//env.out() << sol_array << endl;
		}

		cplex.setParam(IloCplex::Param::Emphasis::MIP, 0); // set back to the default parameter
	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	return;
}

void CPLEX_solver::solve() {
	cplex.solve();
}
