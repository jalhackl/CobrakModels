# Import toy model and pretty-print function and pathway-using ecTFBA
from cobrak.printing import print_dict, print_optimization_result
from cobrak.lps import perform_lp_optimization, perform_lp_variability_analysis

#used for local
from cobrak.nlps import perform_nlp_irreversible_optimization_with_active_reacs_only

from cobrak.evolution import perform_nlp_evolutionary_optimization



def run_cobrak_basic(cobrak_model, pop_size=None, evolution_num_gens=10, objective_target="ATP_Consumption", objective_sense=+1, variability_dict_fixed_params=None, mode="evo", print_result=True,
    after_evo_perform_postprocessing=False, evo_return_all_results=False):


    variability_dict = perform_lp_variability_analysis(
        cobrak_model,
        with_enzyme_constraints=True,
        with_thermodynamic_constraints=True,
    )



    if variability_dict_fixed_params is not None:
        variability_dict = {k: variability_dict_fixed_params[k] if k in variability_dict_fixed_params else v for k, v in variability_dict.items()}



    if mode == "evo":

        all_results = perform_nlp_evolutionary_optimization(
            cobrak_model=cobrak_model,
            objective_target=objective_target,
            objective_sense=objective_sense,
            variability_dict=variability_dict,
            with_kappa=True,
            with_gamma=True,
            with_alpha=False,
            with_iota=False,
            sampling_wished_num_feasible_starts=2,
            #objvalue_json_path="blabla.json",
            evolution_num_gens=evolution_num_gens,
            pop_size=pop_size
        )

        result = all_results[list(all_results.keys())[0]][0] # 0->The first element is the best



        if after_evo_perform_postprocessing:
            from cobrak.evolution import postprocess

            postprocess_results, best_postprocess_result = postprocess(
                cobrak_model=cobrak_model,
                opt_dict=result,
                objective_target=objective_target,
                objective_sense=objective_sense,
                variability_data={}, # No variability dict given -> An ecTFVA is automatically run for us
            )

            print_optimization_result(cobrak_model, best_postprocess_result)

            #if postprocessing and True, return all results (including postprocessing results) for further analysis; otherwise, only return the main result of the evolutionary optimization
            if evo_return_all_results:
                return all_results, result, postprocess_results, best_postprocess_result
        
        #if not performing postprocessing, but still want to return all results of the evolutionary optimization:
        if evo_return_all_results:
            return all_results, result
        

    elif mode == "local":

        ectfba_result = perform_lp_optimization(
        cobrak_model=cobrak_model,
        objective_target=objective_target,
        objective_sense=objective_sense,
        # Set enzyme constraints as they are also used in the NLP
        with_enzyme_constraints=True,
        # This following setting is important to find thermodynamically active reactions
        with_thermodynamic_constraints=True,
        )

        # Run (local and fast) NLP (by default, with the IPOPT solver)
        nlp_result = perform_nlp_irreversible_optimization_with_active_reacs_only(
            cobrak_model,
            objective_target=objective_target,
            objective_sense=objective_sense,
            # Set the suitable set of thermodynamically active reactions
            optimization_dict=ectfba_result,
            # Set the variable bounds from our preparatory variability analysis
            variability_dict=variability_dict,
            # We use the saturation term constraint (otherwise, κ is set to 1 for all reactions);
            # default is True anyway
            with_kappa=True,
            # We use the thermodynamic term constraint (otherwise, γ is set to 1 for all reactions);
            # default is True anyway
            with_gamma=True,
        )

        result = nlp_result



    if print_result:

        print_optimization_result(cobrak_model, result)
    
    return result

    


def run_ectfba_basic(cobrak_model, objective_target="ATP_Consumption", objective_sense=+1, 
                     variability_dict_fixed_params=None, 
                     mode="ectfba", use_variability_dict=True, print_result=True):

    
    if use_variability_dict:
        variability_dict = perform_lp_variability_analysis(
            cobrak_model,
            with_enzyme_constraints=True,
            with_thermodynamic_constraints=True,
        )
    else:
        variability_dict=None
    

    if use_variability_dict:
        if variability_dict_fixed_params is not None:
            variability_dict = {k: variability_dict_fixed_params[k] if k in variability_dict_fixed_params else v for k, v in variability_dict.items()}

    if mode == "ectfba":

        # Perform ecTFBA
        result = perform_lp_optimization(
            cobrak_model=cobrak_model,
            objective_target=objective_target,
            objective_sense=objective_sense,
            with_thermodynamic_constraints=True,
            with_enzyme_constraints=True,

            variability_dict=variability_dict
        )



    if print_result:
        # Pretty print result as dictionary
        print_dict(result)
        print_optimization_result(cobrak_model, result)
    
    return result

    


def remove_reaction(model, reaction_id):
    """
    Remove a reaction and clean up any unused enzymes and metabolites.
    Returns a report of what was removed.
    """

    report = {
        "reaction_removed": None,
        "enzymes_removed": [],
        "metabolites_removed": [],
        "constraints_cleaned": 0,
    }

    # --- Remove reaction ---
    removed_reaction = model.reactions.pop(reaction_id, None)
    if removed_reaction is not None:
        report["reaction_removed"] = reaction_id

    # --- Track used enzymes ---
    used_enzymes = {
        enz_id
        for r in model.reactions.values()
        if r.enzyme_reaction_data
        for enz_id in r.enzyme_reaction_data.identifiers
    }

    old_enzymes = set(model.enzymes.keys())
    removed_enzymes = old_enzymes - used_enzymes

    model.enzymes = {
        k: v for k, v in model.enzymes.items() if k in used_enzymes
    }

    report["enzymes_removed"] = sorted(removed_enzymes)

    # --- Track used metabolites ---
    used_metabolites = set()
    for r in model.reactions.values():
        used_metabolites.update(r.stoichiometries.keys())

    old_metabolites = set(model.metabolites.keys())
    removed_metabolites = old_metabolites - used_metabolites

    model.metabolites = {
        k: v for k, v in model.metabolites.items() if k in used_metabolites
    }

    report["metabolites_removed"] = sorted(removed_metabolites)

    # --- constraints ---
    cleaned_constraints = 0
    for constraint in model.extra_linear_constraints:
        for met in list(constraint.stoichiometries.keys()):
            if met not in model.metabolites:
                constraint.stoichiometries.pop(met)
                cleaned_constraints += 1

    report["constraints_cleaned"] = cleaned_constraints

    return report



import numpy as np

def make_even_grid_with_mean(n_points, mean, p_low=0.5, p_high=1.5, dtype=float):
    """
    Create an evenly spaced grid with 'mean' included, roughly within [p_low*mean, p_high*mean].
    """
    # Determine middle index
    mid_index = n_points // 2

    # Compute approximate low/high
    low_approx = mean * p_low
    high_approx = mean * p_high

    # Compute step so that the total span fits in the approximate range
    total_span = max(mean - low_approx, high_approx - mean) * 2
    step = total_span / (n_points - 1)

    # Start value so mean is at mid_index
    start = mean - step * mid_index

    # Generate evenly spaced grid
    grid = np.array([start + i * step for i in range(n_points)], dtype=dtype)
    
    return grid


def grid_around_value(fixp, step, n_lower, n_upper, dtype=float):
    """
    Create an evenly spaced grid around a fixed value.

    Parameters:
        fixp: the value to include exactly in the grid
        step: spacing between points
        n_lower: number of points below fixp
        n_upper: number of points above fixp
        dtype: int or float

    Returns:
        np.ndarray: evenly spaced grid including fixp
    """
    start_point = fixp - n_lower * step
    end_point = fixp + (n_upper + 1) * step  # +1 to include the uppermost point

    lower_grid = np.arange(start_point, fixp, step, dtype=dtype)
    upper_grid = np.arange(fixp, end_point, step, dtype=dtype)

    grid = np.concatenate([lower_grid, upper_grid])
    return grid