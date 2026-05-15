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

k_cat_bio = k_cat_bio_cal
k_cat_glycolysis = k_cat_glycolysis_cal
k_cat_respiration = k_cat_respiration_cal
k_cat_overflow = k_cat_overflow_cal 



#km transporter values
# Transport_M (pyruvate)
k_cat_transport_M = 3e5
Km_M = 5e-4
Km_M_dict = {"M": 5.0e-4, "M_m": 5.0e-4}
#Km_M_m_dict = {"M_m": 5.0e-4}

# Transport_C (CO2)
k_cat_transport_C = 1e7
Km_C = 1.0
Km_C_dict = {"C": 1, "C_m": 1}
Km_C_m_dict = {"C_m": 1}


# Transport_ATP_ADP
##k_cat_transport_ATP_ADP = 5e5
##Km_ADP = 5e-5
##Km_ATP = 5e-5
##Km_ADP_dict = {"ADP": 5.0e-5, "ADP_m": 5.0e-5}
#Km_ADP_m_dict = {"ADP_m": 5.0e-5}
##Km_ATP_dict = {"ATP": 5.0e-5, "ATP_m": 5.0e-5}
#Km_ATP_m_dict = {"ATP_m": 5.0e-5}

#transporter enzyme weight values
biomass_enzyme_W=1500

ATP_transporter_enzyme_W=500
ADP_transporter_enzyme_W=500

M_transporter_enzyme_W=500
C_transporter_enzyme_W=500

#or
#M_transporter_enzyme_W=150
#C_transporter_enzyme_W=50
#ATP_transporter_enzyme_W=300
#ADP_transporter_enzyme_W=300

