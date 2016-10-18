#!/usr/bin/python
import logging
import math
from collections import defaultdict, Counter
from random import sample
import numpy as np
import cplex as cp
import heapq

"""Find maximal subset of compatible mutation patterns weighted by reliability scores via CPLEX ILP solver"""
__author__ = 'Johannes REITER'
__date__ = 'April, 2014'


# get logger for application
logger = logging.getLogger('treeomics')


def solve_conflicting_phylogeny(cf_graph, time_limit=None):
    """
    Translates given conflict graph into a integer linear program and
    solves the ILP for the minimum number of mutation patterns (set of identical mutation patterns)
    which need to be ignored
    :param cf_graph: Conflict graph: nodes correspond to mutation patterns and edges to their conflicts
    :param time_limit: time limit for MILP solver in seconds
    :return set of conflicting mutation patterns, set of compatible mutation patterns
    """

    logger.info('Build linear program (cplex) for finding the minimal number of conflicting mutations.')

    # the number of columns in the ILP is given by the number of nodes in the conflict graph
    # weighting of the mutation patterns corresponds to the number of mutation which are conflicting
    objective_function = [data['weight'] for _, data in cf_graph.nodes_iter(data=True)]

    logger.debug('Objective function: ' + ', '.join(
        '{}: {:.1e}'.format(var_idx, weight) for var_idx, weight in enumerate(objective_function, 1)))

    # look at the fourth highest reliability score value to get an impression of the magnitude of these values
    ex_rs = heapq.nlargest(4, objective_function)[-1]
    # scale all values to avoid numerical issues with CPLEX
    scaling_factor = 1e20 / ex_rs
    scaled_obj_func = [scaling_factor * x for x in objective_function]

    # column types
    ctypes = ['B' for _ in range(len(objective_function))]

    # add column names to the ILP
    cnames = []
    for col_idx, node in enumerate(cf_graph.nodes_iter(), 0):
        cnames.append(str(node))

    lp = cp.Cplex()

    # set objective function which is the minimum number of positions with a sequencing error
    lp.objective.set_sense(lp.objective.sense.minimize)
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
    for constraint_idx, (source, sink) in enumerate(cf_graph.edges_iter(), 1):
        constraint = [[str(source), str(sink)], [1, 1]]

        constraints.append(constraint)
        row_names.append(str(source)+'-'+str(sink))

        # logger.debug('Add constraint {}: {}'.format(constraint_idx, str(source)+'-'+str(sink)))

    row_rhss = [1 for _ in range(len(constraints))]     # 1 is the RHS in all constraints
    row_senses = ['G' for _ in range(len(constraints))]     # greater equal is used in all constraints
    lp.linear_constraints.add(lin_expr=constraints, senses=row_senses, rhs=row_rhss, names=row_names)

    logger.debug('Added {} constraints.'.format(len(constraints)))

    # solve the Integer Linear Program (ILP)
    lp.solve()
    sol = lp.solution       # obtain solution

    solve_stat = sol.get_status()
    # column solution values
    solution_values = sol.get_values()
    # proportional to incompatible mutations (depending on the weight definition)
    objective_value = sol.get_objective_value()

    logger.info('Minimum vertex cover is of weight (objective value) {:.3e} (original weight: {:3e}).'.format(
        objective_value/scaling_factor, sum(val for val in objective_function)))

    logger.info('Solution status: {}'.format(sol.status[solve_stat]))

    logger.debug('Column solution values: ' +
                 ', '.join('{}: {}'.format(var_idx, status) for var_idx, status in enumerate(solution_values, 1)))

    # translate solution of the ILP back to the phylogeny problem
    # removing the minimum vertex cover (conflicting mutation patterns) gives the maximum compatible set of mps
    compatible_nodes = set()
    incompatible_nodes = set()
    conflicting_mutations_weight = 0

    for col_idx, (node, data) in enumerate(cf_graph.nodes_iter(data=True)):

        if round(sol.get_values(str(node)), 5) == 0:
            compatible_nodes.add(node)
        else:
            incompatible_nodes.add(node)
            conflicting_mutations_weight += data['weight']

    # assert round(conflicting_mutations_weight, 4) == round(objective_value, 4), \
    #     "As long as weights are given by the number of mutations: {} == {}".format(conflicting_mutations_weight,
    #                                                                                objective_value)
    logger.debug('Compatible mutation patterns after the conflicting mps have been removed: {}'.format(
        compatible_nodes))

    return incompatible_nodes, compatible_nodes


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
    for col_idx, node in enumerate(cf_graph.nodes_iter(), 0):
        ilp_col_names.append(str(node))
        ilp_cols[node] = col_idx
        ilp_col_mps.append(node)

    # column types
    var_types = ['B' for _ in range(len(ilp_col_names))]

    # add evolutionary constraints
    constraints = []    # LHS (left hand side) of the rows in the ILP
    row_names = []      # names of the rows (constraints)
    for constraint_idx, (source, sink) in enumerate(cf_graph.edges_iter(), 1):
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
                objective_function[ilp_cols[idx_to_mp[col_idx]]] -= math.log(-math.expm1(log_ml))

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

        # proportional to incompatible mutations (depending on the weight definition)
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

    logger.debug('Finished bootstrapping.')

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
    for constraint_idx, (source, sink) in enumerate(cf_graph.edges_iter(), 1):
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
    for col_idx, (node, data) in enumerate(cf_graph.nodes_iter(data=True), 0):

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
    for constraint_idx, (source, sink) in enumerate(cf_graph.edges_iter(), 1):
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

            for node in cf_graph.nodes_iter():
                if round(sol.get_values(str(node)), 5) == 0 and 1 < len(node) < no_samples:
                    node_frequencies[100-removed_fraction][node] += 1

            # increase objective function values again to the initial values
            for removed_mut in removed_muts:
                objective_function[node_indices[mutations[removed_mut]]] += mut_pattern_scores[removed_mut]

        logger.debug('Finished sampling and solving {:.0%} used fraction of variants.'.format(0.01*removed_fraction))

    return node_frequencies
