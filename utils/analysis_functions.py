import pandas as pd


def compute_total_enzyme_weight(solution_dict, enzyme_dict):
    total_weight = 0.0
    
    for key, value in solution_dict.items():
        # Only look at enzyme variables
        if key.startswith("enzyme_") and "_of_" in key:
            
            # Extract enzyme name
            enzyme_name = key[len("enzyme_"):key.index("_of_")]
            
            # Check that enzyme exists
            if enzyme_name in enzyme_dict:
                mw = enzyme_dict[enzyme_name].molecular_weight
                total_weight += value * mw
            else:
                raise KeyError(f"Enzyme '{enzyme_name}' not found in enzyme_dict.")
    
    return total_weight


def compute_grouped_enzyme_mass(solution_dict, enzyme_dict, groups, other_group_name="Other"):
    enzyme_totals = {}
    total_mass = 0.0

    # accumulate per enzyme
    for key, concentration in solution_dict.items():
        if key.startswith("enzyme_") and "_of_" in key:
            
            enzyme_name = key[len("enzyme_"):key.index("_of_")]

            if enzyme_name not in enzyme_dict:
                raise KeyError(f"Enzyme '{enzyme_name}' not found in enzyme_dict.")

            mw = enzyme_dict[enzyme_name].molecular_weight
            contribution = concentration * mw

            enzyme_totals.setdefault(enzyme_name, 0.0)
            enzyme_totals[enzyme_name] += contribution
            total_mass += contribution

    # group construction
    group_totals = {group: 0.0 for group in groups}
    group_totals[other_group_name] = 0.0

    for enzyme_name, mass in enzyme_totals.items():
        assigned = False
        for group_name, enzyme_list in groups.items():
            if enzyme_name in enzyme_list:
                group_totals[group_name] += mass
                assigned = True
                break

        if not assigned:
            group_totals["Other"] += mass

    return total_mass, group_totals, enzyme_totals



def grouped_wrapper(row, groups, enzyme_dict):
    total_mass, group_totals, enzyme_totals = compute_grouped_enzyme_mass(
        row.to_dict(),
        enzyme_dict,
        groups
    )
    
    return pd.Series({
        "total_mass": total_mass,
        "group_totals": group_totals,
        "enzyme_totals": enzyme_totals
    })



def append_analysis_rows(df, enzyme_dict, groups=None):
    if groups is None:
        groups = dict()

    df["total_enzyme_weight"] = df.apply(
        lambda row: compute_total_enzyme_weight(row.to_dict(), enzyme_dict),
        axis=1
    )

    '''
    def grouped_wrapper(row, groups, enzyme_dict):
        total_mass, group_totals, enzyme_totals = compute_grouped_enzyme_mass(
            row.to_dict(),
            enzyme_dict,
            groups
        )
        
        return pd.Series({
            "total_mass": total_mass,
            "group_totals": group_totals,
            "enzyme_totals": enzyme_totals
        })
    '''

    df[["total_mass", "group_totals", "enzyme_totals"]] = df.apply(
        grouped_wrapper,
        axis=1,
        args=(groups, enzyme_dict)
    )

    return df



import numpy as np

def add_atp_adp_columns(df, mode="full"):
    """
    Adds ATP/ADP concentration and ratio columns.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe.
    mode : str
        "full"  -> expects *_c and *_m columns (compartmentalized)
        "simple" -> uses non-compartmental columns only

    Returns
    -------
    pandas.DataFrame
        Updated dataframe.
    """

    df = df.copy()

    if mode == "full":
        if "x_ATP" in df.columns:
            df["ATP conc., c"] = np.exp(df["x_ATP"])
        if "x_ATP_m" in df.columns:
            df["ATP conc., m"] = np.exp(df["x_ATP_m"])
        if "x_ADP" in df.columns:
            df["ADP conc., c"] = np.exp(df["x_ADP"])
        if "x_ADP_m" in df.columns:
            df["ADP conc., m"] = np.exp(df["x_ADP_m"])

        if "ATP conc., c" in df.columns and "ADP conc., c" in df.columns:
            df["ATP/ADP, c"] = df["ATP conc., c"] / df["ADP conc., c"]

        if "ATP conc., m" in df.columns and "ADP conc., m" in df.columns:
            df["ATP/ADP, m"] = df["ATP conc., m"] / df["ADP conc., m"]

    elif mode == "simple":
        if "x_ATP" in df.columns:
            df["ATP conc."] = np.exp(df["x_ATP"])
        if "x_ADP" in df.columns:
            df["ADP conc."] = np.exp(df["x_ADP"])

        if "ATP conc." in df.columns and "ADP conc." in df.columns:
            df["ATP/ADP"] = df["ATP conc."] / df["ADP conc."]

    return df