# IMPORT SECTION
from math import log

from cobrak.constants import STANDARD_R, STANDARD_T
from cobrak.dataclasses import (
    Enzyme,
    EnzymeReactionData,
    ExtraLinearConstraint,
    Metabolite,
    Model,
    Reaction,
)



def initialize_model_wo_comp(k_cat_glycolysis=140_000, 
                                    k_cat_respiration=140_000,
                                    k_cat_overflow=140_000,
                                    k_cat_biomass=142,
                                    biomass_atp=4,
                                    biomass_enzyme_W=1500,
                                    Km_Biomass_dict = None,
                                    atp_concentration_constraint=False,
                                    
                                    glycolysis_enzyme_W=1000,
                                    respiration_enzyme_W=2500,
                                    overflow_enzyme_W=500,
                                    
                                    biomass_dG0=0):
    
    if not Km_Biomass_dict:
        Km_Biomass_dict={  # Michaelis-Menten constants in M=mol⋅l⁻¹; Default is {}
                        "B": 0.0001,  
                        "ADP": 0.0001,
                        "M": 0.0001,
                        "ATP": 0.0001,
                    }
        

    # EXAMPLE MODEL DEFINITION SECTION
    toy_model = Model(
        reactions={
            # Metabolic reactions
            "Glycolysis": Reaction(
                # Stoichiometrically relevant member variables
                stoichiometries={
                    "S": -1,  # Negative stoichiometry → S is consumed by Glycolysis
                    "ADP": -2,
                    "M": +1,  # Positive stoichiometry → M is produced by Glycolysis
                    "ATP": +2,  # Two ATP molecules are produced by Glycolysis
                },
                min_flux=0.0,  # Minimal flux in mmol⋅gDW⁻¹⋅h⁻¹; should be ≥0 for most analyses
                max_flux=1000.0,  # Maximal flux in mmol⋅gDW⁻¹⋅h⁻¹
                # Thermodynamically relevant member variables
                # (only neccessary if thermodynamic constraints are used)
                dG0=-10.0,  # Standard Gibb's free energy ΔG'° in kJ⋅mol⁻¹; Default is None (no ΔG'°)
                dG0_uncertainty=None,  # ΔG'° uncertainty in kJ⋅mol⁻¹; Default is None (no uncertainty)
                # Let's set the variable for enzyme-kinetic parameters
                # of the dataclass EnzymeReactionData
                # (Default is None, i.e. no enzyme parameters given)
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=[
                        "E_glyc"
                    ],  # Subunit(s) which constitute the reaction's catalyst
                    k_cat=k_cat_glycolysis,  # Turnover number in h⁻¹
                    k_ms={  # Michaelis-Menten constants in M=mol⋅l⁻¹; Default is {}
                        "S": 0.0001,  # e.g., K_m of reaction Glycolysis regarding metabolite A
                        "ADP": 0.0001,
                        "M": 0.0001,
                        "ATP": 0.0001,
                    },
                    special_stoichiometries={},  # No special stoichiometry, all subunits occur once
                ),
                # Extra information member variables
                annotation={"description": "This is reaction Glycolysis"},  # Default is {}
                name="Reaction Glycolysis",  # Default is ""
            ),
            "Respiration": Reaction(
                stoichiometries={
                    "M": -1,
                    "ADP": -4,
                    "C": +1,
                    "ATP": +4,
                },
                min_flux=0.0,
                max_flux=1_000.0,
                dG0=-10.0,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=["E_resp"],
                    k_cat=k_cat_respiration,
                    k_ms={
                        "ADP": 0.00027,
                        "M": 0.00027,
                        "C": 0.0001,
                        "ATP": 0.0001,
                    },
                ),
            ),
            "Overflow": Reaction(
                stoichiometries={
                    "M": -1,
                    "P": +1,
                },
                min_flux=0.0,
                max_flux=1_000.0,
                dG0=-10.0,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=["E_over"],
                    k_cat=k_cat_overflow,
                    k_ms={
                        "M": 0.001,
                        "P": 0.0001,
                    },
                ),
            ),
            # Exchange reactions
            "EX_S": Reaction(
                stoichiometries={
                    "S": +1,
                },
                min_flux=0.0,
                max_flux=1_000.0,
            ),
            "EX_C": Reaction(
                stoichiometries={
                    "C": -1.0,
                },
                min_flux=0.0,
                max_flux=1_000.0,
            ),
            "EX_P": Reaction(
                stoichiometries={
                    "P": -1,
                },
                min_flux=0.0,
                max_flux=1_000.0,
            ),
            
            "ATP_Consumption": Reaction(
                stoichiometries={
                    "ATP": -1,
                    "ADP": +1,
                },
                min_flux=0.0,
                max_flux=1_000.0,
            ),
            


            # biomass part
            "EX_B": Reaction(
                stoichiometries={
                    "B": -1,
                },
                min_flux=0.0,
                max_flux=1_000.0,
            ),

            "Biomass_Rea": Reaction(
                stoichiometries={
                    "ATP": -biomass_atp,
                    "ADP": +biomass_atp,
                    "M": -1,
                    "B": +1
                },
                min_flux=0.0,
                max_flux=1_000.0,
                dG0=biomass_dG0,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=["Biomass_enz"],
                    k_cat=k_cat_biomass,
                ),
                k_ms=Km_Biomass_dict
            ),


        },
        metabolites={
            "S": Metabolite(
                log_min_conc=log(
                    1e-6
                ),  # optional, minimal ln(concentration); Default is ln(1e-6 M)
                log_max_conc=log(
                    0.02
                ),  # optional, maximal ln(concentration); Default is ln(0.02 M)
                annotation={
                    "description": "This is metabolite S"
                },  # optional, default is ""
                name="Metabolite S",  # optional, default is ""
                formula="X",  # optional, default is ""
                charge=0,  # optional, default is 0
            ),
            "M": Metabolite(),
            "C": Metabolite(),
            "P": Metabolite(),
            "ATP": Metabolite(),
            "ADP": Metabolite(),

            "B": Metabolite(),
        },
        enzymes={
            "E_glyc": Enzyme(
                molecular_weight=glycolysis_enzyme_W,  # Molecular weight in kDa
                min_conc=None,  # Optional concentration in mmol⋅gDW⁻¹; Default is None (minimum is 0)
                max_conc=None,  # Optional maximal concentration in mmol⋅gDW⁻¹; Default is None (only protein pool restricts)
                annotation={"description": "Enzyme of Glycolysis"},  # Default is {}
                name="Glycolysis enzyme",  # Default is ""
            ),
            "E_resp": Enzyme(molecular_weight=respiration_enzyme_W),
            "E_over": Enzyme(molecular_weight=overflow_enzyme_W),

            "Biomass_enz": Enzyme(molecular_weight=biomass_enzyme_W)

        },
        extra_linear_constraints = [],

        #comment, because we want to minimize and do not set a fixed value
        #max_prot_pool=max_prot_pool,  # In g⋅gDW⁻¹; This value is used for our analyses with enzyme constraints
        
        
        # We set the following two constraints:
        # 1.0 * EX_A - 1.0 * Glycolysis ≤ 0.0
        # and
        # 1.0 * EX_A + 1.0 * Glycolysis ≥ 0.0
        # in other words, effectively,
        # 1.0 * EX_A = 1.0 * Glycolysis
        #extra_linear_constraints=[
        #    ExtraLinearConstraint(
        #        stoichiometries={
        #            "EX_S": -1.0,
        #            "Glycolysis": 1.0,
        #        },
        #        lower_value=0.0,
        #        upper_value=0.0,
        #    )
        #],  # Keep in mind that this is a list as multiple extra flux constraints are possible
        kinetic_ignored_metabolites=[],
        R=STANDARD_R,
        T=STANDARD_T,
        max_conc_sum=float("inf"),
        annotation={"description": "COBRA-k toy model"},
    )

    if atp_concentration_constraint:
        toy_model.extra_linear_constraints = [
            ExtraLinearConstraint(
                stoichiometries={
                    "x_ATP": 1.0,
                    "x_ADP": -1.0,
                },
                lower_value=log(3.0),
            )
        ]

    return toy_model