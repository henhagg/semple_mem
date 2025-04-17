library(ggplot2)
library(gridExtra)
library(latex2exp)
library(rjson)
library(ggpubr)
library(tidyr)
library(dplyr)

source(file.path("plot_scripts", "plot_ppc.R"), local = TRUE)
source(file.path("plot_scripts", "plot_kde_kalman_and_semple_with_prior.R"), local = TRUE)
source(file.path("plot_scripts", "plot_kde_pepsdi_and_semple.R"), local = TRUE)

plot_mc_separate_ind_param_from_csv = function(input_dir,
                                               round_index,
                                               individual_index,
                                               fontsize = 15,
                                               pdf_width = 5,
                                               pdf_height = 3) {
  
  param = read.csv(file.path(input_dir, paste0("ind_param_round", round_index, ".csv")))
  settings = rjson::fromJSON(file = file.path(input_dir, "settings.json"))
  
  true_ind_param_file_path = file.path(dirname(settings$true_param_file_name), "ind_param.csv")
  if(file.exists(true_ind_param_file_path)){
    true_ind_param = read.csv(true_ind_param_file_path)
  } else {
    true_ind_param = NULL
  }

  param_names = colnames(param)
  
  output_dir = file.path(input_dir, "plots")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  param_names_file_path = file.path("models", settings$model_name, "param_names.csv")
  if(file.exists(param_names_file_path)) {
    axis_labels = read.csv(param_names_file_path)
  } else {
    axis_labels = param_names
    names(axis_labels) = param_names
  }
  
  for (i in 1:(length(param_names) - 2)) {
    param_name = sym(param_names[i])
    
    
    plot = ggplot(data = param %>% filter(id == individual_index)) +
      geom_line(mapping = aes(x = gibbs_cycle, y = !!param_name)) +
      labs(x = "Iteration", y = TeX(as.character(axis_labels[param_names[i]]))) +
      theme(text = element_text(size = fontsize))
    if(!is.null(true_ind_param)){
      plot = plot + geom_hline(aes(yintercept = true_ind_param[i, individual_index]),
                               color = "red",
                               linetype = "dashed",
                               linewidth = 1)
    }
    
    output_file_name = paste0("mc_ind",
                              individual_index,
                              "_c",
                              i,
                              "_round",
                              round_index,
                              ".pdf")
    ggsave(
      filename = file.path(output_dir, output_file_name),
      width = pdf_width,
      height = pdf_height
    )
  }

}

param_name_to_tex = function(param_name) {
  texified_string = sub("(.*)(\\d)", "\\1_\\2", param_name)
  TeX(paste0("$\\", texified_string,"$"))
}

plot_mc_shared_param = function(input_dir,
                                param_type,
                                round_index,
                                fontsize = 15,
                                pdf_width = 5,
                                pdf_height = 3) {
  settings = rjson::fromJSON(file = file.path(input_dir, "settings.json"))
  true_param = settings$true_param
  if(true_param[[1]] == "NA") {
    true_param_exists = FALSE
  } else {
    true_param_exists = TRUE
  }
  
  output_dir = file.path(input_dir, "plots")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  param = read.csv(file.path(input_dir, paste0(param_type, "_round", round_index, ".csv")))
  param_names = names(param)
  param$gibbs_cycle = 1:nrow(param)
  
  param_names_file_path = file.path("models", settings$model_name, "param_names.csv")
  if(file.exists(param_names_file_path)) {
    axis_labels = read.csv(param_names_file_path)
  } else {
    axis_labels = param_names
    names(axis_labels) = param_names
  }
  
  for (i in 1:length(param_names)) {
    param_name_sym = sym(param_names[i])
    
    plot = ggplot(param) +
      geom_line(mapping = aes(x = gibbs_cycle, y = !!param_name_sym)) +
      labs(x = "Iteration", y = TeX(as.character(axis_labels[param_names[i]]))) +
      theme(text = element_text(size = fontsize))
    if (true_param_exists) {
      plot = plot + geom_hline(aes(yintercept = as.numeric(true_param[param_names[i]])),
                               color = "red",
                               linetype = "dashed",
                               linewidth = 1)
    }
    
    
    output_file_name = paste0("mc_", param_names[i], "_round", round_index, ".pdf")
    ggsave(
      filename = file.path(output_dir, output_file_name),
      width = pdf_width,
      height = pdf_height
    )
  }
}

plot_kde = function(param, true_param) {
  param_names = names(param)
  kde = lapply(1:length(param_names), function(i) {
    x = sym(param_names[i])
    ggplot(param, aes(x = !!x)) +
      geom_density() +
      labs(x = param_name_to_tex(param_names[i]), y = "") +
      geom_vline(aes(xintercept = true_param[[param_names[i]]]),
                 color = "red",
                 linewidth = 0.5,
                 linetype = "dashed")
  })
}

plot_kde_multiple_rounds = function(input_dir,
                                    param_type,
                                    rounds,
                                    axis_labels,
                                    true_param = NULL,
                                    xlims = NULL,
                                    trim = TRUE) {
  param = read.csv(file.path(input_dir, paste0(param_type, "_round", rounds[1], ".csv")))
  param_names = names(param)
  param$round = rep(rounds[1], nrow(param))
  for (i in 2:length(rounds)) {
    param_next = read.csv(file.path(input_dir, paste0(param_type, "_round", rounds[i], ".csv")))
    param_next$round = rep(rounds[i], nrow(param_next))
    param = rbind(param, param_next)
  }
  param$round = cut(param$round, breaks = length(rounds), labels = rounds)

  kde_plot_list = lapply(1:length(param_names), function(i) {
    x = sym(param_names[i])
    round = sym("round")
    p = ggplot(param) +
      geom_density(aes(
        x = !!x,
        group = !!round,
        color = !!round
      ), trim = trim) +
      labs(x = TeX(as.character(axis_labels[param_names[i]])), y = "", color = "") +
      scale_color_manual(labels = c("Prior samples", "Posterior samples"), values = c("#1B9E77", "#D95F02"))
    if(!is.null(xlims)){
      p = p + xlim(xlims[[i]][1], xlims[[i]][2])
    }
    if(!is.null(true_param)){
      p = p + geom_vline(
        aes(xintercept = true_param[[param_names[i]]]),
        color = "black",
        linewidth = 0.5,
        linetype = "dashed"
      )
    }
    p = p + theme(plot.margin = margin(0, 0.2, 0, 0, "cm"),
                  axis.title = element_text(size = 15, face = "bold"),
                  legend.text = element_text(size = 10),
                  axis.text = element_text(size = 10))
    return(p)
  })
  return(kde_plot_list)
}

plot_kde_multiple_param_multiple_rounds = function(input_dir,
                                                   param_types,
                                                   rounds,
                                                   xlims = NULL,
                                                   ncol = NULL,
                                                   nrow = NULL,
                                                   pdf_width = 7,
                                                   pdf_height = 7) {
  
  settings = rjson::fromJSON(file = file.path(input_dir, "settings.json"))
  true_param = settings$true_param
  if(true_param[[1]] == "NA") {true_param = NULL}
  
  output_dir = file.path(input_dir, "plots")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  axis_labels = read.csv(file.path("models", settings$model_name, "param_names.csv"))
  
  plots = lapply(1:length(param_types), function(i) {
    plot_kde_multiple_rounds(
      input_dir,
      param_types[i],
      rounds,
      true_param,
      axis_labels = axis_labels,
      xlims = xlims[[i]]
    )
  })
  
  pdf(file.path(output_dir, paste0("kde_multiround.pdf")),
      width = pdf_width,
      height = pdf_height)
  print(ggarrange(
    plotlist = unlist(plots, recursive = F),
    common.legend = T,
    ncol = ncol,
    nrow = nrow
  ))
  dev.off()
}

plot_kde_ind_param_multiple_rounds = function(input_dir,
                                              individual_indices,
                                              xlims = NULL,
                                              round_index = 2,
                                              pdf_width = 7,
                                              pdf_height = 7) {
  settings = rjson::fromJSON(file = file.path(input_dir, "settings.json"))
  output_dir = file.path(input_dir, "plots")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  true_ind_param_file_path = file.path(dirname(settings$true_param_file_name), "ind_param.csv")
  if (file.exists(true_ind_param_file_path)) {
    true_ind_param = read.csv(true_ind_param_file_path)
  } else {
    true_ind_param = NULL
  }
  
  param_prior = read.csv(file.path(input_dir, "ind_param_round0.csv"))
  param_names = names(param_prior)
  param_prior$round = rep(0, nrow(param_prior))
  num_ind_param = length(param_names)
  
  param_names_file_path = file.path("models", settings$model_name, "param_names.csv")
  if(file.exists(param_names_file_path)) {
    axis_labels = read.csv(param_names_file_path)
  } else {
    axis_labels = param_names
    names(axis_labels) = param_names
  }
  
  for (individual_index in individual_indices) {
    ind_param = read.csv(file.path(input_dir, paste0(
      "ind_param_round", round_index, ".csv"
    )))
    param_posterior = ind_param %>% dplyr::filter(id == individual_index) %>% dplyr::select(!c(gibbs_cycle, id))
    param_posterior$round = rep(2, nrow(param_posterior))
    
    param = rbind(param_prior, param_posterior)
    param$round = cut(param$round, breaks = 2, labels = c(0, 2))
    
    kde_plot_list = lapply(1:length(param_names), function(i) {
      x = sym(param_names[i])
      round = sym("round")
      p = ggplot(param) +
        geom_density(aes(
          x = !!x,
          group = !!round,
          color = !!round
        )) +
        labs(x = TeX(as.character(axis_labels[param_names[i]])),
             y = "",
             color = "") +
        scale_color_manual(
          labels = c("Prior samples", "Posterior samples"),
          values = c("#1B9E77", "#D95F02")
        )
      if (!is.null(true_ind_param)) {
        true_ind_param_list = as.list(true_ind_param[, individual_index])
        names(true_ind_param_list) = param_names
        p = p + geom_vline(
          aes(xintercept = true_ind_param_list[[param_names[i]]]),
          color = "black",
          linewidth = 0.5,
          linetype = "dashed"
        )
      }
      if(!is.null(xlims)){
        p = p + xlim(xlims[[i]][1], xlims[[i]][2])
      }
      p = p + theme(plot.margin = margin(0, 0.2, 0, 0, "cm"),
                    axis.title = element_text(size = 12, face = "bold"),
                    legend.text = element_text(size = 10),
                    axis.text = element_text(size = 8))
      return(p)
    })
    
    pdf(file.path(
      output_dir,
      paste0("kde_multiround_ind", individual_index, ".pdf")
    ),
    width = pdf_width,
    height = pdf_height)
    print(ggarrange(
      plotlist = kde_plot_list,
      common.legend = T,
      ncol = num_ind_param
    ))
    dev.off()
  }
}

plot_full_analysis = function(input_dir) {
  settings = rjson::fromJSON(file = file.path(input_dir, "settings.json"))
  
  plot_kde_multiple_param_multiple_rounds(
    input_dir,
    param_types = c("eta", "kappa_xi"),
    rounds = c(0, 2:settings$num_rounds)
  )
  
  for (r in 2:settings$num_rounds) {
    plot_mc_shared_param(input_dir, param_type = "kappa_xi", round_index =  r)
    plot_mc_shared_param(input_dir, param_type = "eta", round_index =  r)
    plot_mc_separate_ind_param_from_csv(input_dir,
                                        round_index = r,
                                        individual_index = 1)
  }
  
  plot_ppc(input_dir, round_index = 2)
  for (i in 1:min(settings$num_observation, 5)) {
    plot_ppc_individual(
      input_dir = input_dir,
      round_index = 2,
      individual_index = i
    )
  }
}

plot_analysis_only_individual_param = function(input_dir) {
  settings = rjson::fromJSON(file = file.path(input_dir, "settings.json"))
  
  plot_kde_multiple_param_multiple_rounds(
    input_dir,
    param_types = c("eta"),
    rounds = c(0, 2:settings$num_rounds)
  )
  
  for (r in 2:settings$num_rounds) {
    plot_mc_shared_param(input_dir, param_type = "eta", round_index =  r)
    plot_mc_separate_ind_param_from_csv(input_dir,
                                        round_index = r,
                                        individual_index = 1)
  }
  
  plot_ppc(input_dir, round_index = 2)
  for (i in 1:min(settings$num_observation, 5)) {
    plot_ppc_individual(
      input_dir = input_dir,
      round_index = 2,
      individual_index = i
    )
  }
}

plot_observed_data = function(input_dir) {
  observed_data = read.csv(file.path(input_dir, "observations.csv"))
  observed_data_gathered = observed_data %>% pivot_longer(!c(time), names_to = "observation", values_to = "value")
  ggplot(
    observed_data_gathered,
    aes(x = time, y = value, colour = observation)) +
    labs(x = "Time (hours)", y = "Measurement") + geom_line() + theme_minimal() + theme(legend.position = "none") 
}

plot_bic_from_file = function(input_file,
                              save_to_file = FALSE,
                              pdf_width = NA,
                              pdf_height = NA) {
  bic_data = read.csv(file = input_file)
  ggplot(bic_data, aes(x = K, y = BIC)) + geom_line() + geom_point()
  
  if (save_to_file) {
    ggsave(
      filename = file.path(dirname(input_file), "bic.pdf"),
      width = pdf_width,
      height = pdf_height
    )
  }
}

######################## BIC PLOTS ##########################
# plot_bic_from_file(
#   input_file = "results/ornstein_uhlenbeck_unperturbed_noise/bic/10k/bic.csv",
#   save_to_file = TRUE,
#   pdf_width = 3.5,
#   pdf_height = 2
# )
# 
# plot_bic_from_file(
#   input_file = "results/mrna_fix_tzero/bic/50k/bic.csv",
#   save_to_file = TRUE,
#   pdf_width = 3.5,
#   pdf_height = 2
# )
# 
# plot_bic_from_file(
#   input_file = "results/mrna_indep_prior/bic/10k/bic.csv",
#   save_to_file = TRUE,
#   pdf_width = 3.5,
#   pdf_height = 2
# )

########### ORNSTEIN-UHLENBECK #############
# plot_kde_kalman_and_semple_with_prior(
#   input_dir_semple = "results/ornstein_uhlenbeck_unperturbed_noise/num_observation_40/10k_samples_hmc",
#   kalman_samples_file = "models/ornstein_uhlenbeck_kalman/results_priormeanxi-1_2/kalman_post_samples_full.csv",
#   burnin_kalman = 100000,
#   thinning_kalman = 10,
#   xlims = list(c(-2,2), c(0, 3), c(-2, 2), NULL, NULL, NULL, c(-3,0)),
#   ylims = list(NULL, NULL, NULL, NULL, NULL, NULL, c(0,5)),
#   pdf_width = 4,
#   pdf_height = 7,
#   font_size_axis = 15,
#   font_size_legend = 10,
#   ncol = 2,
#   nrow = 4
# )
# plot_ppc(input_dir = "results/ornstein_uhlenbeck_unperturbed_noise/num_observation_40/10k_samples_hmc",
#          round_index = 2,
#          pdf_width = 3,
#          pdf_height = 2,
#          font_size_axis = 9,
#          font_size_ticks = 9)
# plot_mc_shared_param(
#   input_dir = "results/ornstein_uhlenbeck_unperturbed_noise/num_observation_40/10k_samples_hmc",
#   param_type = "kappa_xi",
#   round_index =  2
# )
# plot_mc_shared_param(
#   input_dir = "results/ornstein_uhlenbeck_unperturbed_noise/num_observation_40/10k_samples_hmc",
#   param_type = "eta",
#   round_index =  2
# )
# plot_mc_separate_ind_param_from_csv(
#   input_dir = "results/ornstein_uhlenbeck_unperturbed_noise/num_observation_40/10k_samples_hmc",
#   round_index = 2,
#   individual_index = 1
# )
# plot_kde_ind_param_multiple_rounds(
#   input_dir = "results/ornstein_uhlenbeck_unperturbed_noise/num_observation_40/10k_samples_hmc",
#   individual_indices = 1:5,
#   xlims = NULL,
#   round_index = 2,
#   pdf_width = 4.5,
#   pdf_height = 2
# )

############## MRNA REAL DATA KDE ONLY INDIVIDUAL PARAM ###############
# plot_kde_multiple_param_multiple_rounds(
#   input_dir = "results/mrna_indep_prior_only_individual_param/egfp_40ind/burnin10",
#   param_types = c("eta"),
#   xlims = list(eta = list(
#     c(-4, 2),
#     c(-8, 0),
#     c(-4, 4),
#     c(-2, 2),
#     c(2, 8),
#     c(-4, 5),
#     c(1, 5),
#     c(-4, 2),
#     c(0, 25),
#     c(0, 25),
#     c(0, 25),
#     c(0, 25),
#     c(0, 25),
#     c(0, 25),
#     c(0, 40),
#     c(0, 25)
#   )),
#   rounds = c(0, 2),
#   ncol = 2,
#   nrow = 8,
#   pdf_width = 4,
#   pdf_height = 10
# )

############## MRNA FIX TZERO WITH FIXED EFFECTS #############
# plot_kde_pepsdi_and_semple(
#   input_dir_semple = "results/mrna_fix_tzero/num_observation_40/10k_hmc_randomseed3_K10",
#   semple_round_index = 2,
#   input_dir_pepsdi = "../PEPSDI_results/40ind/attempt9_50ksamples/Npart1000_nsamp50000_corr0.99_exp_id1_run1",
#   sample_indices_pepsdi = 1:50000,
#   thinning_pepsdi = 5,
#   include_semple_prior = TRUE,
#   xlims = list(c(-2,1), c(-4,-2), c(-4,4), c(0,25), c(0,25), c(0,25), c(2,8), c(-4,5), c(1.5,2.5), c(-2,-1)),
#   font_size_axis_label = 15,
#   font_size_ticks = 8,
#   pdf_width = 3.5,
#   pdf_height = 7,
#   ncol = 2,
#   nrow = 5
# )
# plot_kde_ind_param_multiple_rounds(
#   input_dir = "results/mrna_fix_tzero/num_observation_40/10k_hmc_randomseed3_K10",
#   individual_indices = 1:5,
#   xlims = list(c(-4, 2.5), c(-5, -1), c(-2,2)),
#   round_index = 2,
#   pdf_width = 4.5,
#   pdf_height = 2
# )
# plot_mc_separate_ind_param_from_csv(
#   input_dir = "results/mrna_fix_tzero/num_observation_40/10k_hmc_randomseed3_K10",
#   round_index = 2,
#   individual_index = 1
# )
# plot_mc_shared_param(
#   input_dir = "results/mrna_fix_tzero/num_observation_40/10k_hmc_randomseed3_K10",
#   param_type = "kappa_xi",
#   round_index =  2
# )
# plot_mc_shared_param(
#   input_dir = "results/mrna_fix_tzero/num_observation_40/10k_hmc_randomseed3_K10",
#   param_type = "eta",
#   round_index =  2
# )

############# MRNA REAL DATA FIXED EFFECTS ####################
# plot_kde_multiple_param_multiple_rounds(
#   input_dir = "results/mrna_indep_prior/egfp_40ind/K5",
#   param_types = c("eta", "kappa_xi"),
#   xlims = list(
#   eta = list(c(-4,2), c(-8,0), c(-4,4), c(-2,2), c(0, 25), c(0, 25), c(0, 25), c(0, 25)),
#   kappa_xi = list(c(2, 8), c(-4, 5), c(1, 5), c(-4, 2))
#   ),
#   rounds = c(0, 2),
#   ncol = 2,
#   nrow = 6,
#   pdf_width = 4,
#   pdf_height = 8
# )
# plot_ppc(
#   input_dir = "results/mrna_indep_prior/egfp_40ind/K5",
#   round_index = 2,
#   pdf_width = 3,
#   pdf_height = 2,
#   font_size_axis = 9,
#   font_size_ticks = 9
# )
# for(i in 1:5) {
#   plot_ppc_individual(
#     input_dir = "results/mrna_indep_prior/egfp_40ind/K5",
#     round_index = 2,
#     individual_index = i,
#     pdf_width = 3,
#     pdf_height = 2,
#     xlab = "Time (hours)",
#     ylab = "Measurement"
#   )
# }
# plot_mc_separate_ind_param_from_csv(
#   input_dir = "results/mrna_indep_prior/egfp_40ind/K5",
#   round_index = 2,
#   individual_index = 1
# )
# plot_mc_shared_param(
#   input_dir = "results/mrna_indep_prior/egfp_40ind/K5",
#   param_type = "kappa_xi",
#   round_index =  2
# )
# plot_mc_shared_param(
#   input_dir = "results/mrna_indep_prior/egfp_40ind/K5",
#   param_type = "eta",
#   round_index =  2
# )
