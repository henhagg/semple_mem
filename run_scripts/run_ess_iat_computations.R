library("LaplacesDemon")
library("mcmcse")

source(file.path("algorithm", "compute_ess_and_iat.R"), local = TRUE)

##### mRNA model simulated data

# SeMPLE
compute_ess_and_iat_semple(
    input_dir = file.path("results", "mrna_fix_tzero", "num_observation_40", "50k_R4"),
    semple_round_index = 4
)

# PEPSDI
compute_ess_and_iat_pepsdi(input_dir = file.path("results", "PEPSDI", "40ind", "attempt9_50ksamples", "Npart1000_nsamp50000_corr0.99_exp_id1_run1"))
