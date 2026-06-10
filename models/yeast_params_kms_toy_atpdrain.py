k_cat_glycolysis=1.19 * 10**6
k_cat_respiration=4.35 * 10**6
k_cat_overflow=1.68 * 10**6

k_cat_bio=126
k_cat_transport=7.2 * 10**5
#p.16 metabolicshiftsyeast
k_cat_glycolysis_cal = k_cat_glycolysis * 0.43
k_cat_respiration_cal = k_cat_respiration * 0.0077
k_cat_overflow_cal = k_cat_overflow * 1.1
k_cat_bio_cal = 142

k_cat_glycolysis = k_cat_glycolysis_cal
k_cat_respiration = k_cat_respiration_cal
k_cat_overflow = k_cat_overflow_cal
k_cat_bio = k_cat_bio_cal



#transporter enzyme weight values
biomass_enzyme_W=500

ATP_transporter_enzyme_W=50
ADP_transporter_enzyme_W=50

M_transporter_enzyme_W=50
C_transporter_enzyme_W=50

glycolysis_enzyme_W=500
respiration_enzyme_W=1000-4*M_transporter_enzyme_W

fermentation_enzyme_W=100
#or
#M_transporter_enzyme_W=150
#C_transporter_enzyme_W=50
#ATP_transporter_enzyme_W=300
#ADP_transporter_enzyme_W=300



Km_Glycolysis_dict={  # Michaelis-Menten constants in M=mol⋅l⁻¹; Default is {}
                        "S": 0.0001,  # e.g., K_m of reaction Glycolysis regarding metabolite A
                        "ADP": 0.0001,
                        "M": 0.0001,
                        "ATP": 0.0001,
                    }



Km_Respiration_dict={
                        "ADP_m": 0.00027,
                        "M_m": 0.00027,
                        "C_m": 0.0001,
                        "ATP_m": 0.0001,
                    }



Km_Overflow_dict={
                        "M": 0.001,
                        "P": 0.0001,
                    }


Km_Biomass_dict={  # Michaelis-Menten constants in M=mol⋅l⁻¹; Default is {}
                        "B": 0.0001,  
                        "ADP": 0.0001,
                        "M": 0.0001,
                        "ATP": 0.0001,
                    }





