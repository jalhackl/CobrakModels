# IMPORT SECTION
from math import log
import re

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
                                    Km_Biomass_dict=None,
                                    atp_concentration_constraint=False,
                                    
                                    glycolysis_enzyme_W=1000,
                                    respiration_enzyme_W=2500,
                                    overflow_enzyme_W=500,
                                    
                                    biomass_dG0=0,
                                    atp_consumption_dG=0,

                                    Km_Glycolysis_dict=None,
                                    Km_Overflow_dict=None,
                                    Km_Respiration_dict=None,
                                    remove_compartment_info=True,

                                    # dG0 values
                                    glycolysis_dG0=-10.0,
                                    respiration_dG0=-10.0,
                                    overflow_dG0=-10.0,
                                    transport_dg0=None):

    if not Km_Biomass_dict:
        Km_Biomass_dict = {}

    if not Km_Glycolysis_dict:
        Km_Glycolysis_dict = {}

    if not Km_Overflow_dict:
        Km_Overflow_dict = {}

    if not Km_Respiration_dict:
        Km_Respiration_dict = {}

    if remove_compartment_info:

        def remove_suffix(d):
            return {
                re.sub(r'_[A-Za-z]$', '', key): value
                for key, value in d.items()
            }

        Km_Biomass_dict = remove_suffix(Km_Biomass_dict)
        Km_Glycolysis_dict = remove_suffix(Km_Glycolysis_dict)
        Km_Overflow_dict = remove_suffix(Km_Overflow_dict)
        Km_Respiration_dict = remove_suffix(Km_Respiration_dict)

    if transport_dg0 is not None:
        glycolysis_dG0 = glycolysis_dG0
        respiration_dG0 = respiration_dG0
        overflow_dG0 = overflow_dG0
        biomass_dG0 = biomass_dG0

    # EXAMPLE MODEL DEFINITION SECTION
    toy_model = Model(
        reactions={
            "Glycolysis": Reaction(
                stoichiometries={
                    "S": -1,
                    "ADP": -2,
                    "M": +1,
                    "ATP": +2,
                },
                min_flux=0.0,
                max_flux=1000.0,
                dG0=glycolysis_dG0,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=["E_glyc"],
                    k_cat=k_cat_glycolysis,
                    k_ms=Km_Glycolysis_dict,
                    special_stoichiometries={},
                ),
                annotation={"description": "This is reaction Glycolysis"},
                name="Reaction Glycolysis",
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
                dG0=respiration_dG0,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=["E_resp"],
                    k_cat=k_cat_respiration,
                    k_ms=Km_Respiration_dict,
                ),
            ),

            "Overflow": Reaction(
                stoichiometries={
                    "M": -1,
                    "P": +1,
                },
                min_flux=0.0,
                max_flux=1_000.0,
                dG0=overflow_dG0,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=["E_over"],
                    k_cat=k_cat_overflow,
                    k_ms=Km_Overflow_dict,
                ),
            ),

            "EX_S": Reaction(
                stoichiometries={"S": +1},
                min_flux=0.0,
                max_flux=1_000.0,
            ),

            "EX_C": Reaction(
                stoichiometries={"C": -1.0},
                min_flux=0.0,
                max_flux=1_000.0,
            ),

            "EX_P": Reaction(
                stoichiometries={"P": -1},
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
                dG0=atp_consumption_dG
            ),
        },

        metabolites={
            "S": Metabolite(
                log_min_conc=log(1e-6),
                log_max_conc=log(0.02),
                annotation={"description": "This is metabolite S"},
                name="Metabolite S",
                formula="X",
                charge=0,
            ),
            "M": Metabolite(),
            "C": Metabolite(),
            "P": Metabolite(),
            "ATP": Metabolite(),
            "ADP": Metabolite(),
        },

        enzymes={
            "E_glyc": Enzyme(
                molecular_weight=glycolysis_enzyme_W,
                min_conc=None,
                max_conc=None,
                annotation={"description": "Enzyme of Glycolysis"},
                name="Glycolysis enzyme",
            ),
            "E_resp": Enzyme(molecular_weight=respiration_enzyme_W),
            "E_over": Enzyme(molecular_weight=overflow_enzyme_W),
        },

        extra_linear_constraints=[],
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

