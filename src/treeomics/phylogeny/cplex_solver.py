#!/usr/bin/python
import logging
import math
from collections import defaultdict, Counter
from random import sample
import numpy as np
from scipy.special import logsumexp
import cplex as cp
import heapq
from itertools import islice
from phylogeny.solution import Solution

"""Find maximal subset of compatible mutation patterns weighted by reliability scores via CPLEX ILP solver"""
__author__ = 'Johannes REITER'
__date__ = 'April, 2014'


# get logger for application
logger = logging.getLogger('treeomics')

NO_PLOTTED_SOLUTIONS = 5


def solve_conflicting_phylogeny(cf_graph, no_muts, pool_size, no_plotted_solutions, time_limit=None, n_max_threads=0):
    """
    Translates given conflict graph into a integer linear program and
    solves the ILP for the minimum number of mutation patterns (set of identical mutation patterns)
    which need to be ignored
    :param cf_graph: Conflict graph: nodes correspond to mutation patterns and edges to their conflicts
    :param no_muts: Number of processed mutations
    :param pool_size: number of best solutions explored by ILP solver to estimate confidence
    :param time_limit: time limit for MILP solver in seconds
    :param n_max_threads: Sets the default maximal number of parallel threads that will be invoked by CPLEX
      (0: default, let CPLEX decide; 1: single threaded; N: uses up to N threads)
      https://www.ibm.com/support/knowledgecenter/en/SS9UKU_12.5.0/com.ibm.cplex.zos.help/Parameters/topics/Threads.html
    :param no_plotted_solutions: number of best solutions from the solution pool that will be plotted
    :return list of top solutions, dictionary with calculated node likelihoods based on likelihood of
            each solution in the solution pool
    """

    logger.info('Build linear program (cplex) for finding the minimal number of conflicting mutations.')

    # the number of columns in the ILP is given by the number of nodes in the conflict graph
    # weighting of the mutation patterns corresponds to the number of mutation which are conflicting
    objective_function = [data['weight'] for _, data in cf_graph.nodes(data=True)]

    logger.debug('Objective function: ' + ', '.join(
        '{}: {:.1e}'.format(var_idx, weight) for var_idx, weight in enumerate(objective_function, 1)))

    # look at the fourth highest reliability score value to get an impression of the magnitude of these values
    ex_rs = heapq.nlargest(min(5, cf_graph.order()), objective_function)[-1]
    # scale all values to avoid numerical issues with CPLEX
    Solution.SCALING_FACTOR = 1e20 / ex_rs
    scaled_obj_func = [Solution.SCALING_FACTOR * x for x in objective_function]

    # column types
    ctypes = ['B' for _ in range(len(objective_function))]

    # add column names to the ILP
    cnames = []
    ilp_col_mps = []
    for col_idx, node in enumerate(cf_graph.nodes(), 0):
        cnames.append(str(node))
        ilp_col_mps.append(node)

    lp = cp.Cplex()

    # set objective function which is the minimum number of positions with a sequencing error
    lp.objective.set_sense(lp.objective.sense.minimize)
    lp.parameters.threads.set(n_max_threads)
    # see more information at:
    # http://www.ibm.com/support/knowledgecenter/en/SS9UKU_12.4.0/com.ibm.cplex.zos.help/Parameters/topics/SolnPoolGap.html?view=embed
    # lp.parameters.mip.tolerances.absmipgap.set(1e-15)
    # lp.parameters.mip.tolerances.mipgap.set(1e-15)
    # lp.parameters.simplex.tolerances.optimality.set(1e-09)

    # set time limit for MILP solver
    if time_limit is not None:
        lp.parameters.timelimit.set(time_limit)

    lp.variables.add(obj=scaled_obj_func, types=ctypes, names=cnames)

    # add evolutionary constraints
    constraints = []    # LHS (left hand side) of the rows in the ILP
    row_names = []      # names of the rows (constraints)
    for constraint_idx, (source, sink) in enumerate(cf_graph.edges(), 1):
        constraint = [[str(source), str(sink)], [1, 1]]

        constraints.append(constraint)
        row_names.append(str(source)+'-'+str(sink))

        # logger.debug('Add constraint {}: {}'.format(constraint_idx, str(source)+'-'+str(sink)))

    row_rhss = [1 for _ in range(len(constraints))]     # 1 is the RHS in all constraints
    row_senses = ['G' for _ in range(len(constraints))]     # greater equal is used in all constraints
    lp.linear_constraints.add(lin_expr=constraints, senses=row_senses, rhs=row_rhss, names=row_names)

    logger.debug('Added {} constraints.'.format(len(constraints)))

    # # 3... Node file on disk and compressed
    # lp.parameters.mip.strategy.file.set(3)
    # lp.parameters.workmem.set(2048)

    # ################ explore the solution space by keeping a pool of the best solutions ###############
    # more information at:
    # https://www.ibm.com/support/knowledgecenter/SS9UKU_12.5.0/com.ibm.cplex.zos.help/Parameters/topics/PopulateLim.html
    # and populate.py an example within the CPLEX installation
    if pool_size > 1:
        lp.parameters.mip.limits.populate.set(pool_size)
        # strategy for replacing a solution in the solution pool when the solution pool has reached its capacity
        lp.parameters.mip.pool.replace.set(1)       # 1...replace the solution which has the worst objective

        # Controls the trade-off between the number of solutions generated for the solution pool and
        # the amount of time or memory consumed.
        lp.parameters.mip.pool.intensity.set(4)     # 4...very aggressive: enumerate all practical solutions
        # lp.parameters.mip.pool.capacity.set(pool_size)

        # set the solution pool relative gap parameter to obtain solutions
        # of objective value within 10% of the optimal
        # lp.parameters.mip.pool.relgap.set(5)
        try:
            lp.populate_solution_pool()         # solve the Integer Linear Program (ILP)
        except cp.exceptions.CplexSolverError as e:
            logger.error("Exception raised during populate")
            raise e
    else:
        lp.solve()

    # assess obtained solutions
    solutions, weighted_node_lh = assess_solutions(
        lp.solution, objective_function, cf_graph, ilp_col_mps, no_muts, pool_size, no_plotted_solutions)

    return solutions, weighted_node_lh


def assess_solutions(sol, objective_function, cf_graph, ilp_col_mps, no_muts, pool_size, no_plotted_solutions):
    """
    Assess solution pool and infer the compatible and incompatible mutation patterns as well as their support
    :param sol: solutions obtained by the ILP solver
    :param objective_function: objective function provided to the ILP solver
    :param cf_graph: conflict graph
    :param ilp_col_mps: list with mapping from column ids to nodes (mutation patterns)
    :param no_muts: Number of processed mutations
    :param pool_size: number of best solutions explored by ILP solver to estimate confidence
    :param no_plotted_solutions: number of solutions where the inferred trees will be plotted
    :return list of top solutions, dictionary with calculated node likelihoods based on likelihood of
            each solution in the solution pool
    """

    assert pool_size >= no_plotted_solutions, \
        'Solution pool size needs to be larger than the number of plotted solutions!'
    assert 1 <= no_plotted_solutions <= 20, 'At least 1 and not more than 20 solutions can be plotted: {}'.format(
        no_plotted_solutions)

    no_sols = sol.pool.get_num()
    logger.debug("The solution pool contains {} solutions.".format(no_sols))
    solve_stat = sol.get_status()

    # column solution values
    solution_values = sol.get_values()

    # proportional to incompatible mutations (depending on the weight definition)
    obj_value = sol.get_objective_value()

    opt_sol = Solution(1, objective_function, obj_value, solution_values, cf_graph, no_muts)

    # get likelihood of optimal solution
    opt_lh = math.exp(opt_sol.llh)
    logger.debug('Likelihood of optimal solution: {:.2e}'.format(opt_lh))

    logger.debug('Minimum vertex cover is of weight (objective value) {:.3e} (original weight: {:3e}).'.format(
        opt_sol.get_unscaled_obj_value(), sum(val for val in objective_function)))
    logger.info('Solution status: {}'.format(sol.status[solve_stat]))
    logger.debug('Column solution values: ' +
                 ', '.join('{}: {}'.format(var_idx, status) for var_idx, status in enumerate(solution_values, 1)))

    # meanobjval = lp.solution.pool.get_mean_objective_value() / scaling_factor
    # logger.debug("The average objective value of the solutions is {:.3e}".format(meanobjval))

    total_obj_value = sum(sol.pool.get_objective_value(i) for i in range(no_sols)) / Solution.SCALING_FACTOR
    min_obj_value = min(sol.pool.get_objective_value(i) for i in range(no_sols)) / Solution.SCALING_FACTOR
    max_obj_value = max(sol.pool.get_objective_value(i) for i in range(no_sols)) / Solution.SCALING_FACTOR
    logger.debug('Best solution: {:.2e}; worst solution in pool: {:.2e}; sum of all solutions: {:.2e}'.format(
        min_obj_value, max_obj_value, total_obj_value))

    # print the objective value of each solution and its difference to the incumbent
    # names = lp.solution.pool.get_names()
    # print("Solution     Objective      Number of variables    log weight   weight   Likelihood")
    # print("             value          that differ compared")
    # print("                            to the incumbent")

    # list of the top ranked solutions
    solutions = list()
    # log likelihoods of all solutions in the pool
    llhs = list()

    if pool_size > 1:
        # record the chosen patterns in the solution space (pool)
        node_frequencies = Counter()
        # weighted and summed likelihood of each node across the explored solution pool
        weighted_node_lh = defaultdict(float)
        weighted_node_llh = defaultdict(list)
        # weighted_node_counts = defaultdict(float)
        # total_log_sol_weight = 0.0

        rank = 1
        for i in sorted(range(no_sols), key=lambda k: sol.pool.get_objective_value(k)):

            if rank == 1:
                solution = opt_sol

            else:
                if all(prev_vals == cur_vals for prev_vals, cur_vals in zip(solution_values, sol.pool.get_values(i))):
                    continue    # same solution exists twice in solution pool

                # list of almost binary variables of a mutation pattern is present or absent in the solution
                solution_values = sol.pool.get_values(i)

                # get objective function value of the i'th solution
                obj_value = sol.pool.get_objective_value(i)

                solution = Solution(rank, objective_function, obj_value, solution_values, cf_graph, no_muts)

            llhs.append(solution.llh)

            if rank <= no_plotted_solutions:
                solutions.append(solution)

            # calculate a solution weight depending on the best and the worst values
            # log_sol_weight = 1.0 - ((objval_i - min_obj_value) / (max_obj_value - min_obj_value))
            # sol_weight = 1.0 - ((math.exp(-objval_i) - math.exp(-min_obj_value)) /
            #                     (math.exp(-max_obj_value) - math.exp(-min_obj_value)))
            # total_log_sol_weight += log_sol_weight

            for ilp_col_idx, val in enumerate(solution_values):
                if round(val, 5) == 0:  # mutation pattern was selected
                    node_frequencies[ilp_col_mps[ilp_col_idx]] += 1
                    weighted_node_lh[ilp_col_mps[ilp_col_idx]] += math.exp(solution.llh)
                    weighted_node_llh[ilp_col_mps[ilp_col_idx]].append(solution.llh)
                    # weighted_node_counts[ilp_col_mps[ilp_col_idx]] += log_sol_weight

                    # compute the number of variables that differ in solution i and in the incumbent
                    # no_differences = 0
                    # for j in range(len(objective_function)):
                    #     if round(solution_values[j] - solution_values[j], 6) > 0:
                    #         no_differences += 1

                    # print every solution
                    # print("{:>5s}   {:13g}       {:>3} / {:>3} {:>22.2%} {:>9.2%} {:>10.2e}".format(
                    #     names[i], objval_i, no_differences, len(objective_function), log_sol_weight, sol_weight, lh))

            rank += 1

        for node in node_frequencies.keys():
            # weighted_node_counts[node] /= total_log_sol_weight
            node_frequencies[node] /= len(llhs)     # normalize by the number of produced solutions
            # weighted_node_lh[node] *= sum(math.exp(llh) for llh in llhs)

            # sum likelihoods of solutions where node was present
            llh_sum_pres_sols = logsumexp(weighted_node_llh[node])
            # sum likelihood of all solutions
            llh_sum_all_sols = logsumexp(llhs)
            # likelihood of node across all evaluated solutions
            weighted_node_lh[node] = math.exp(llh_sum_pres_sols - llh_sum_all_sols)

        logger.debug('Mutation patterns sorted by their weighted likelihood:')
        for node, weighted_lh in islice(sorted(weighted_node_lh.items(), key=lambda k: -k[1]), 0, 50):
            logger.debug('Weighted lh {:.4%} (freq: {:.2%}): {}'.format(
                weighted_lh, node_frequencies[node], ', '.join(str(n) for n in node)))
            # logger.debug('opt lh {:.3e}; total lh {:.3e}'.format(opt_lh, sum(lhs)))

    else:   # only one solution was inferred
        solutions.append(opt_sol)
        weighted_node_lh = None

    return solutions, weighted_node_lh


def bootstrapping_solving(cf_graph, mp_weights, idx_to_mp, no_samples):
    """
    Generate and solve MILP of the down-sampled data-set and track the identified
    mutation pattern occurrences
    :param cf_graph: Conflict graph: nodes correspond to mutation patterns and edges to their conflicts
    :param mp_weights: 2-dimensional array with log probability that this variant has this mutation pattern
    :param idx_to_mp: dictionary from mutation patterns to the column ids used in the mp_weights array
    :param no_samples: Number of samples with replacement for the bootstrapping
    :return: observed occurrences of mutation patterns
    """

    # record the chosen patterns in the down-sampled data set
    node_frequencies = Counter()

    logger.debug('Build linear programs (cplex) for the robustness analysis through bootstrapping.')

    # add column names to the ILP
    ilp_col_names = []
    ilp_col_mps = []
    ilp_cols = dict()
    for col_idx, node in enumerate(cf_graph.nodes(), 0):
        ilp_col_names.append(str(node))
        ilp_cols[node] = col_idx
        ilp_col_mps.append(node)

    # column types
    var_types = ['B' for _ in range(len(ilp_col_names))]

    # add evolutionary constraints
    constraints = []    # LHS (left hand side) of the rows in the ILP
    row_names = []      # names of the rows (constraints)
    for constraint_idx, (source, sink) in enumerate(cf_graph.edges(), 1):
        constraint = [[str(source), str(sink)], [1, 1]]
        constraints.append(constraint)
        row_names.append(str(source)+'-'+str(sink))

        # logger.debug('Add constraint {}: {}'.format(constraint_idx, str(source)+'-'+str(sink)))

    row_rhss = [1 for _ in range(len(constraints))]     # 1 is the RHS in all constraints
    row_senses = ['G' for _ in range(len(constraints))]     # greater equal is used in all constraints
    logger.debug('Generated {} constraints.'.format(len(constraints)))
    logger.info('Do bootstrapping with {} samples.'.format(no_samples))

    m = len(mp_weights)   # number of variants
    for rep in range(no_samples):

        # obtain sample of used variants
        used_muts = np.random.choice(m, m, replace=True)
        # the number of columns in the ILP is given by the number of nodes in the conflict graph
        # weighting of the mutation patterns corresponds to the number of mutation which are conflicting
        objective_function = np.zeros(cf_graph.order())

        # update objective function (mutation pattern scores)
        # decrease objective function values according to the removed patterns
        for used_mut in used_muts:
            for col_idx, log_ml in mp_weights[used_mut].items():
                # add the (negative log probability) part of the reliability score of this mutation in this pattern
                # note we are in log space
                if idx_to_mp[col_idx] in ilp_cols.keys():
                    objective_function[ilp_cols[idx_to_mp[col_idx]]] -= math.log(-math.expm1(log_ml))
                # else: for parsimony-uninformative mutation patterns nothing needs to be done

        # logger.debug('Update objective function: ' + ', '.join(
        #     '{}: {:.3f}'.format(var_idx, weight) for var_idx, weight in enumerate(objective_function, 1)))

        # generate new MILP
        lp = cp.Cplex()
        # lp.set_error_stream(None)
        # lp.set_warning_stream(None)
        lp.set_results_stream(None)
        lp.set_log_stream(None)
        lp.variables.add(obj=objective_function, types=var_types, names=ilp_col_names)
        lp.objective.set_sense(lp.objective.sense.minimize)
        lp.linear_constraints.add(lin_expr=constraints, senses=row_senses, rhs=row_rhss, names=row_names)

        # solve the Integer Linear Program (ILP)
        lp.solve()
        sol = lp.solution       # obtain solution
        # solve_stat = sol.get_status()
        # logger.debug('Solution status: {}'.format(sol.status[solve_stat]))
        #
        # # proportional to incompatible mutations (depending on the weight definition)
        # objective_value = sol.get_objective_value()
        # logger.debug('Minimum vertex cover is of weight (objective value) {:.4f} (original weight: {:4f}).'
        #               .format(objective_value, sum(val for val in objective_function)))

        # column solution values
        # solution_values = sol.get_values()
        # logger.debug('Column solution values: ' + ', '.join(
        #     '{}: {}'.format(var_idx, status) for var_idx, status in enumerate(solution_values, 1)))

        solution_values = sol.get_values()
        for ilp_col_idx, val in enumerate(solution_values):
            if round(val, 5) == 0:
                node_frequencies[ilp_col_mps[ilp_col_idx]] += 1

        if no_samples >= 100 and rep > 0 and rep % (no_samples/100) == 0:
            logger.debug('{:.0%} of bootstrapping completed.'.format(1.0*rep/no_samples))

    return node_frequencies


def solve_downsampled_nodes(cf_graph, mp_weights, col_ids_mp, no_replications):
    """
    Generate and solve MILP of the down-sampled data-set and track the identified
    mutation pattern occurrences
    :param cf_graph: Conflict graph: nodes correspond to mutation patterns and edges to their conflicts
    :param mp_weights: 2-dimensional array with log probability that this variant has this mutation pattern
    :param col_ids_mp: dictionary from mutation patterns to the column ids used in the mp_weights array
    :param no_replications: Number of replications per used fraction of variants
    :return: observed occurrences of mutation patterns per variant fraction
    """

    # record the chosen patterns in the down-sampled data set
    node_frequencies = defaultdict(Counter)

    logger.debug('Build linear programs (cplex) for the robustness analysis through down-sampling.')

    # the number of columns in the ILP is given by the number of nodes in the conflict graph
    # weighting of the mutation patterns corresponds to the number of mutation which are conflicting
    objective_function = []

    # add column names to the ILP
    var_names = []
    for col_idx, mp in sorted(col_ids_mp.items(), key=lambda k: k[0]):
        var_names.append(str(mp))
        objective_function.append(cf_graph.node[mp]['weight'])

    logger.debug('Objective function: ' + ', '.join(
        '{}: {:.3f}'.format(var_idx, weight) for var_idx, weight in enumerate(objective_function, 1)))

    # column types
    var_types = ['B' for _ in range(len(objective_function))]

    # add evolutionary constraints
    constraints = []    # LHS (left hand side) of the rows in the ILP
    row_names = []      # names of the rows (constraints)
    for constraint_idx, (source, sink) in enumerate(cf_graph.edges(), 1):
        constraint = [[str(source), str(sink)], [1, 1]]
        constraints.append(constraint)
        row_names.append(str(source)+'-'+str(sink))

        # logger.debug('Add constraint {}: {}'.format(constraint_idx, str(source)+'-'+str(sink)))

    row_rhss = [1 for _ in range(len(constraints))]     # 1 is the RHS in all constraints
    row_senses = ['G' for _ in range(len(constraints))]     # greater equal is used in all constraints
    logger.debug('Generated {} constraints.'.format(len(constraints)))

    mut_ids = [i for i in range(len(mp_weights))]
    for removed_fraction in range(90, 0, -10):
        for rep in range(no_replications):

            # obtain sample of used shared mutations
            removed_muts = sample(mut_ids, int(round(0.01*removed_fraction*len(mut_ids))))

            # update objective function (mutation pattern scores)
            # decrease objective function values according to the removed patterns
            for col_idx, mp in sorted(col_ids_mp.items(), key=lambda k: k[0]):
                for removed_mut in removed_muts:
                    # detract the part of the reliability score of this mutation in this pattern
                    # note we are in log space
                    objective_function[col_idx] += math.log(-math.expm1(mp_weights[removed_mut][col_idx]))

            # logger.debug('Update objective function: ' + ', '.join(
            #     '{}: {:.3f}'.format(var_idx, weight) for var_idx, weight in enumerate(objective_function, 1)))

            # generate new MILP
            lp = cp.Cplex()
            # lp.set_error_stream(None)
            # lp.set_warning_stream(None)
            lp.set_results_stream(None)
            lp.set_log_stream(None)
            lp.variables.add(obj=objective_function, types=var_types, names=var_names)
            lp.objective.set_sense(lp.objective.sense.minimize)
            lp.linear_constraints.add(lin_expr=constraints, senses=row_senses, rhs=row_rhss, names=row_names)

            # solve the Integer Linear Program (ILP)
            lp.solve()
            sol = lp.solution       # obtain solution
            # solve_stat = sol.get_status()
            # logger.debug('Solution status: {}'.format(sol.status[solve_stat]))

            # proportional to incompatible mutations (depending on the weight definition)
            # objective_value = sol.get_objective_value()
            # logger.debug('Minimum vertex cover is of weight (objective value) {:.4f} (original weight: {:4f}).'
            #               .format(objective_value, sum(val for val in objective_function)))

            # column solution values
            # solution_values = sol.get_values()
            # logger.debug('Column solution values: ' + ', '.join(
            #     '{}: {}'.format(var_idx, status) for var_idx, status in enumerate(solution_values, 1)))

            solution_values = sol.get_values()
            for col_idx, val in enumerate(solution_values):
                if round(val, 5) == 0:
                    node_frequencies[100-removed_fraction][col_ids_mp[col_idx]] += 1

            # increase objective function values according to the removed values to its original value
            for col_idx, mp in sorted(col_ids_mp.items(), key=lambda k: k[0]):
                for removed_mut in removed_muts:
                    # detract the part of the reliability score of this mutation in this pattern
                    # note we are in log space
                    objective_function[col_idx] -= math.log(-math.expm1(mp_weights[removed_mut][col_idx]))

        logger.debug('Finished sampling and solving {:.0%} used fraction of variants.'.format(0.01*removed_fraction))

    return node_frequencies


def solve_downsampled_binary_nodes(cf_graph, mut_pattern_scores, shared_mutations, no_replications, no_samples):
    """
    Generate and solve MILP of the down-sampled data-set and track the identified
    mutation pattern occurrences when each variant has exactly one mutation pattern
    :param cf_graph: Conflict graph: nodes correspond to mutation patterns and edges to their conflicts
    :param mut_pattern_scores: Mutation pattern score of each mutation (key: mut_idx)
    :param shared_mutations: List of the shared (parsimony-informative) mutations (mut_idx)
    :param no_replications: Number of replications per used fraction of variants
    :param no_samples: number of samples
    :return: observed occurrences of mutation patterns per variant fraction
    """

    # record the chosen patterns in the down-sampled data set
    node_frequencies = defaultdict(Counter)

    logger.debug('Build linear programs (cplex) for the robustness analysis through down-sampling.')

    # the number of columns in the ILP is given by the number of nodes in the conflict graph
    # weighting of the mutation patterns corresponds to the number of mutation which are conflicting
    objective_function = []

    # build index from mutations to patterns
    mutations = dict()

    # add column names to the ILP
    var_names = []
    node_indices = dict()   # map nodes (mutation patterns) to column id in the ILP
    for col_idx, (node, data) in enumerate(cf_graph.nodes(data=True), 0):

        var_names.append(str(node))
        objective_function.append(data['weight'])
        # nodes are given by a frozenset of samples (mutation patterns)
        node_indices[node] = col_idx
        for mut_idx in data['muts']:
            mutations[mut_idx] = node

    logger.debug('Objective function: ' + ', '.join(
        '{}: {:.3f}'.format(var_idx, weight) for var_idx, weight in enumerate(objective_function, 1)))

    # column types
    var_types = ['B' for _ in range(len(objective_function))]

    # add evolutionary constraints
    constraints = []    # LHS (left hand side) of the rows in the ILP
    row_names = []      # names of the rows (constraints)
    for constraint_idx, (source, sink) in enumerate(cf_graph.edges(), 1):
        constraint = [[str(source), str(sink)], [1, 1]]
        constraints.append(constraint)
        row_names.append(str(source)+'-'+str(sink))

        # logger.debug('Add constraint {}: {}'.format(constraint_idx, str(source)+'-'+str(sink)))

    row_rhss = [1 for _ in range(len(constraints))]     # 1 is the RHS in all constraints
    row_senses = ['G' for _ in range(len(constraints))]     # greater equal is used in all constraints
    logger.debug('Generated {} constraints.'.format(len(constraints)))

    for removed_fraction in range(95, 0, -5):
        for rep in range(no_replications):

            # obtain sample of used shared mutations
            removed_muts = sample(shared_mutations, int(round(0.01*removed_fraction*len(shared_mutations))))

            # update objective function (mutation pattern scores)
            # decrease objective function values according to the removed patterns
            updated_nodes = set()
            for removed_mut in removed_muts:
                updated_nodes.add(mutations[removed_mut])
                objective_function[node_indices[mutations[removed_mut]]] -= mut_pattern_scores[removed_mut]

            # logger.debug('Update objective function: ' + ', '.join(
            #     '{}: {:.3f}'.format(var_idx, weight) for var_idx, weight in enumerate(objective_function, 1)))

            # generate new MILP
            lp = cp.Cplex()
            # lp.set_error_stream(None)
            # lp.set_warning_stream(None)
            lp.set_results_stream(None)
            lp.set_log_stream(None)
            lp.variables.add(obj=objective_function,
                             types=var_types, names=var_names)
            lp.objective.set_sense(lp.objective.sense.minimize)
            lp.linear_constraints.add(lin_expr=constraints, senses=row_senses, rhs=row_rhss, names=row_names)

            # solve the Integer Linear Program (ILP)
            lp.solve()
            sol = lp.solution       # obtain solution
            # solve_stat = sol.get_status()
            # logger.debug('Solution status: {}'.format(sol.status[solve_stat]))

            # proportional to incompatible mutations (depending on the weight definition)
            # objective_value = sol.get_objective_value()
            # logger.debug('Minimum vertex cover is of weight (objective value) {:.4f} (original weight: {:4f}).'
            #               .format(objective_value, sum(val for val in objective_function)))

            # column solution values
            # solution_values = sol.get_values()
            # logger.debug('Column solution values: ' + ', '.join(
            #     '{}: {}'.format(var_idx, status) for var_idx, status in enumerate(solution_values, 1)))

            for node in cf_graph.nodes():
                if round(sol.get_values(str(node)), 5) == 0 and 1 < len(node) < no_samples:
                    node_frequencies[100-removed_fraction][node] += 1

            # increase objective function values again to the initial values
            for removed_mut in removed_muts:
                objective_function[node_indices[mutations[removed_mut]]] += mut_pattern_scores[removed_mut]

        logger.debug('Finished sampling and solving {:.0%} used fraction of variants.'.format(0.01*removed_fraction))

    return node_frequencies
