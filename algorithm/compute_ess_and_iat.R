library("LaplacesDemon")
library("mcmcse")

compute_ess_and_iat_semple <- function(input_dir, semple_round_index, output_significant_digits = 3) {
    # read samples from file
    samples_eta <- read.csv(file.path(input_dir, paste0("eta_round", semple_round_index, ".csv")))
    samples_kappa_xi <- read.csv(file.path(input_dir, paste0("kappa_xi_round", semple_round_index, ".csv")))
    samples <- cbind(samples_eta, samples_kappa_xi)

    # compute univariate ESS
    ess <- LaplacesDemon::ESS(samples)
    ess["min_ESS_univariate"] <- min(ess)

    # compute multivariate ESS
    ess["ESS_multivariate"] <- mcmcse::multiESS(samples)

    # compute IAT
    iat <- t(sapply(1:ncol(samples), function(i) LaplacesDemon::IAT(samples[, i])))
    colnames(iat) <- colnames(samples)

    # read runtime from file
    elapsed_time <- read.csv(file.path(input_dir, "elapsed_time.csv"))

    if (!is.null(elapsed_time$round)) { # if there are multiple rounds in elapsed time
        runtime <- elapsed_time$seconds[elapsed_time$round == semple_round_index]
    } else {
        runtime <- elapsed_time$seconds
    }

    # compute ess per second
    ess_per_second <- ess / runtime

    # round to X significant digits
    ess <- signif(ess, digits = output_significant_digits)
    ess_per_second <- signif(ess_per_second, digits = output_significant_digits)
    iat <- signif(iat, digits = output_significant_digits)

    # write results to files
    write.csv(t(ess), file = file.path(input_dir, "ess.csv"), quote = FALSE, row.names = FALSE)
    write.csv(t(ess_per_second), file = file.path(input_dir, "ess_per_second.csv"), quote = FALSE, row.names = FALSE)
    write.csv(iat, file = file.path(input_dir, "iat.csv"), quote = FALSE, row.names = FALSE)
}


compute_ess_and_iat_pepsdi <- function(input_dir, output_significant_digits = 3) {
    # read samples from csv
    samples_mu <- read.csv(file.path(input_dir, "Mean.csv"))
    samples_scale <- read.csv(file.path(input_dir, "Scale.csv"))
    samples_kappa_sigma <- read.csv(file.path(input_dir, "Kappa_sigma.csv"))

    # Match PEPSDI parameter names with SeMPLE naming
    colnames(samples_scale) <- paste0("tau", 1:ncol(samples_scale))
    colnames(samples_kappa_sigma)[which(names(samples_kappa_sigma) == "sigma1")] <- "xi1"

    samples <- cbind(samples_mu, samples_scale, samples_kappa_sigma)

    # compute univariate ESS
    ess <- LaplacesDemon::ESS(samples)
    ess["min_ess"] <- min(ess)

    # compute multivariate ESS
    ess["ESS_multivariate"] <- mcmcse::multiESS(samples)

    # compute IAT
    iat <- t(sapply(1:ncol(samples), function(i) LaplacesDemon::IAT(samples[, i])))
    colnames(iat) <- colnames(samples)

    # read runtime from file
    runtime_milliseconds <- read.csv(file.path(input_dir, "Run_time.csv"))$Run_time

    # compute ess per second
    runtime_seconds <- runtime_milliseconds / 1000
    ess_per_second <- ess / runtime_seconds

    # round to X significant digits
    ess <- signif(ess, digits = output_significant_digits)
    ess_per_second <- signif(ess_per_second, digits = output_significant_digits)
    iat <- signif(iat, digits = output_significant_digits)

    # write results to files
    write.csv(t(ess), file = file.path(input_dir, "ess.csv"), quote = FALSE, row.names = FALSE)
    write.csv(t(ess_per_second), file = file.path(input_dir, "ess_per_second.csv"), quote = FALSE, row.names = FALSE)
    write.csv(iat, file = file.path(input_dir, "iat.csv"), quote = FALSE, row.names = FALSE)
}
