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



def initialize_model_w_comp(k_cat_glycolysis=140_000, 
                                    k_cat_respiration=140_000,
                                    k_cat_overflow=140_000,
                                    k_cat_transport=140_000,

                                    # dG0 parameters
                                    glycolysis_dG0=-10.0,
                                    respiration_dG0=-10.0,
                                    overflow_dG0=-10.0,
                                    transport_dg0=None,
                                    biomass_dG0=0,

                                    k_cat_biomass=142,
                                    biomass_atp=4,
                                    biomass_enzyme_W=1500,


                                    M_transporter_enzyme_W=50,
                                    C_transporter_enzyme_W=50,
                                    ATP_transporter_enzyme_W=50,
                                    ADP_transporter_enzyme_W=50,

                                    atp_concentration_constraint=False,
                                    atp_m_concentration_constraint=False,
                                    
                                    # transporter Kcat
                                    kcat_transport_dict=None,
                                    Km_M_dict=None, Km_C_dict=None, Km_ATP_dict=None, Km_ADP_dict=None, 
                                    Km_Biomass_dict=None,
                                    glycolysis_enzyme_W=1000,
                                    respiration_enzyme_W=2500,
                                    overflow_enzyme_W=500,

                                    Km_Glycolysis_dict=None,
                                    Km_Overflow_dict=None,
                                    Km_Respiration_dict=None
                                    ):
    
    kcat_transport_default_value = 1e20
    if k_cat_transport is not None:
        kcat_transport_default_value = k_cat_transport

    if not kcat_transport_dict:
        kcat_transport_dict = {"kcat_M": None, "kcat_C": None, "kcat_ATP": None, "kcat_ADP": None}

        for key, value in kcat_transport_dict.items():
            if value is None:
                kcat_transport_dict[key] = kcat_transport_default_value

    if not Km_M_dict:
        Km_M_dict = {}

    if not Km_C_dict:
        Km_C_dict = {}

    if not Km_ATP_dict:
        Km_ATP_dict = {}

    if not Km_ADP_dict:
        Km_ADP_dict = {}

    if not Km_Biomass_dict:
        Km_Biomass_dict = {}

    if not Km_Glycolysis_dict:
        Km_Glycolysis_dict = {}

    if not Km_Overflow_dict:
        Km_Overflow_dict = {}

    if not Km_Respiration_dict:
        Km_Respiration_dict = {}

    # apply global transport dG0 if provided
    if transport_dg0 is not None:
        glycolysis_dG0 = glycolysis_dG0
        respiration_dG0 = respiration_dG0
        overflow_dG0 = overflow_dG0
        biomass_dG0 = biomass_dG0

    toy_model_compartments = Model(
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
                    "M_m": -1,
                    "ADP_m": -4,
                    "C_m": +1,
                    "ATP_m": +4,
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

            "Transport_M": Reaction(
                stoichiometries={
                    "M": -1,
                    "M_m": +1
                },
                min_flux=0.0,
                max_flux=1_000.0,
                dG0=0.0,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=["M_transporter"],
                    k_cat=kcat_transport_dict["kcat_M"],
                    k_ms=Km_C_dict
                ),
            ),

            "Transport_C": Reaction(
                stoichiometries={
                    "C_m": -1,
                    "C": +1
                },
                min_flux=0.0,
                max_flux=1_000.0,
                dG0=0.0,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=["C_transporter"],
                    k_cat=kcat_transport_dict["kcat_C"],
                    k_ms=Km_C_dict
                ),
            ),

            "Transport_ATP": Reaction(
                stoichiometries={
                    "ATP_m": -1,
                    "ATP": +1
                },
                min_flux=0.0,
                max_flux=1_000.0,
                dG0=0.0,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=["ATP_transporter"],
                    k_cat=kcat_transport_dict["kcat_ATP"],
                    k_ms=Km_ATP_dict
                ),
            ),

            "Transport_ATP_back": Reaction(
                stoichiometries={
                    "ATP": -1,
                    "ATP_m": +1
                },
                min_flux=0.0,
                max_flux=1_000.0,
                dG0=0.0,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=["ATP_transporter"],
                    k_cat=kcat_transport_dict["kcat_ATP"],
                    k_ms=Km_ATP_dict
                ),
            ),

            "Transport_ADP": Reaction(
                stoichiometries={
                    "ADP": -1,
                    "ADP_m": +1,
                },
                min_flux=0.0,
                max_flux=1_000.0,
                dG0=0.0,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=["ADP_transporter"],
                    k_cat=kcat_transport_dict["kcat_ADP"],
                    k_ms=Km_ADP_dict
                ),
            ),

            "Transport_ADP_back": Reaction(
                stoichiometries={
                    "ADP_m": -1,
                    "ADP": +1,
                },
                min_flux=0.0,
                max_flux=1_000.0,
                dG0=0.0,
                enzyme_reaction_data=EnzymeReactionData(
                    identifiers=["ADP_transporter"],
                    k_cat=kcat_transport_dict["kcat_ADP"],
                    k_ms=Km_ADP_dict
                ),
            ),

            "EX_S": Reaction({"S": +1}, 0.0, 1_000.0),
            "EX_C": Reaction({"C": -1.0}, 0.0, 1_000.0),
            "EX_P": Reaction({"P": -1}, 0.0, 1_000.0),

            "ATP_Consumption": Reaction(
                stoichiometries={"ATP": -1, "ADP": +1},
                min_flux=0.0,
                max_flux=1_000.0,
            ),

            "EX_B": Reaction(
                stoichiometries={"B": -1},
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
                    k_ms=Km_Biomass_dict
                ),
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
            "ATP_m": Metabolite(compartment="m"),
            "ADP_m": Metabolite(compartment="m"),
            "M_m": Metabolite(compartment="m"),
            "C_m": Metabolite(compartment="m"),
            "B": Metabolite(),
        },

        enzymes={
            "E_glyc": Enzyme(molecular_weight=glycolysis_enzyme_W),
            "E_resp": Enzyme(molecular_weight=respiration_enzyme_W),
            "E_over": Enzyme(molecular_weight=overflow_enzyme_W),
            "ATP_transporter": Enzyme(molecular_weight=ATP_transporter_enzyme_W),
            "ADP_transporter": Enzyme(molecular_weight=ADP_transporter_enzyme_W),
            "M_transporter": Enzyme(molecular_weight=M_transporter_enzyme_W),
            "C_transporter": Enzyme(molecular_weight=C_transporter_enzyme_W),
            "Biomass_enz": Enzyme(molecular_weight=biomass_enzyme_W),
        },

        kinetic_ignored_metabolites=[],
        R=STANDARD_R,
        T=STANDARD_T,
        max_conc_sum=float("inf"),
        annotation={"description": "COBRA-k toy model"},
        extra_linear_constraints=[],
    )

    if atp_concentration_constraint:
        toy_model_compartments.extra_linear_constraints = [
            ExtraLinearConstraint(
                stoichiometries={
                    "x_ATP": 1.0,
                    "x_ADP": -1.0,
                },
                lower_value=log(3.0),
            )
        ]

    if atp_m_concentration_constraint:
        toy_model_compartments.extra_linear_constraints.append(
            ExtraLinearConstraint(
                stoichiometries={
                    "x_ATP_m": 1.0,
                    "x_ADP_m": -1.0,
                },
                lower_value=log(3.0),
            )
        )

    return toy_model_compartments
