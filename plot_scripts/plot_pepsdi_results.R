plot_full_pepsdi_analysis <- function(input_dir,
                                      true_pop_param_list,
                                      observed_data_file,
                                      num_plot_columns = 3) {
  # read samples from csv
  samples_mu <- read.csv(file.path(input_dir, "Mean.csv"))
  samples_scale <- read.csv(file.path(input_dir, "Scale.csv"))
  samples_kappa_sigma <- read.csv(file.path(input_dir, "Kappa_sigma.csv"))
  #   samples_ind_param = read.csv(file.path(input_dir, "Ind_param.csv"))

  # match PEPSDI parameter names with SeMPLE naming
  colnames(samples_scale) <- paste0("tau", 1:ncol(samples_scale))
  colnames(samples_kappa_sigma)[which(names(samples_kappa_sigma) == "sigma1")] <- "xi1"

  # read parameter names to use as axis labels from file
  param_names_file_path <- file.path("models", "mrna_fix_tzero", "param_names.csv")
  if (file.exists(param_names_file_path)) {
    axis_labels <- read.csv(param_names_file_path)
  } else {
    axis_labels <- param_names
    names(axis_labels) <- param_names
  }

  # create output dir for plots
  output_dir <- file.path(input_dir, "plots")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # plot traceplots
  plot_traceplots_pop_param(samples_mu, axis_labels, true_pop_param_list, output_dir)
  plot_traceplots_pop_param(samples_scale, axis_labels, true_pop_param_list, output_dir)
  plot_traceplots_pop_param(samples_kappa_sigma, axis_labels, true_pop_param_list, output_dir)
  # plot_traceplots_ind_param(samples_ind_param, ind_index = 1.0, output_dir)
}

plot_pvc_pepsdi <- function(input_dir,
                            observed_data_file,
                            pdf_width = 5,
                            pdf_height = 3) {
  pvc_ind_param <- read.csv(file.path(input_dir, "Pvc_ind_param.csv"))
  pvc_quant <- read.csv(file.path(input_dir, "Pvc_quant.csv"))
  observed_data <- read.csv(observed_data_file)

  ggplot() +
    geom_ribbon(
      data = pvc_quant,
      mapping = aes(x = time, ymin = y1_qu05_low, ymax = y1_qu95_upp),
      alpha = 0.2
    ) +
    geom_line(
      data = observed_data,
      mapping = aes(
        x = time,
        y = obs,
        group = id,
        colour = id
      )
    ) +
    ylab("") +
    theme(legend.position = "none")

  ggsave(
    filename = file.path(input_dir, "plots", "pvc.pdf"),
    width = pdf_width,
    height = pdf_height
  )
}

plot_traceplots_pop_param <- function(samples,
                                      axis_labels,
                                      true_pop_param_list,
                                      output_dir,
                                      fontsize = 15,
                                      pdf_width = 5,
                                      pdf_height = 3) {
  param_names <- names(samples)
  samples$gibbs_cycle <- 1:nrow(samples)

  for (i in 1:length(param_names)) {
    param_name_sym <- sym(param_names[i])
    ggplot(samples) +
      geom_line(mapping = aes(x = gibbs_cycle, y = !!param_name_sym)) +
      geom_hline(aes(yintercept = as.numeric(true_pop_param_list[param_names[i]])),
        color = "red",
        linetype = "dashed"
      ) +
      labs(x = "Iteration", y = TeX(as.character(axis_labels[param_names[i]]))) +
      theme_bw() +
      theme(text = element_text(size = fontsize))

    output_file_name <- paste0("traceplot_", param_names[i], ".pdf")
    ggsave(
      filename = file.path(output_dir, output_file_name),
      width = pdf_width,
      height = pdf_height
    )
  }
}

plot_traceplots_ind_param <- function(samples,
                                      ind_index,
                                      output_dir,
                                      fontsize = 15,
                                      pdf_width = 5,
                                      pdf_height = 3) {
  samples <- samples %>% filter(id == ind_index)
  samples <- subset(samples, select = -c(id))
  param_names <- names(samples)
  samples$gibbs_cycle <- 1:nrow(samples)

  for (i in 1:length(param_names)) {
    param_name_sym <- sym(param_names[i])
    ggplot(samples) +
      geom_line(mapping = aes(x = gibbs_cycle, y = !!param_name_sym)) +
      labs(x = "Iteration", y = param_name_to_tex(paste0("log ", param_names[i]))) +
      theme_bw() +
      theme(text = element_text(size = fontsize))

    output_file_name <- paste0("traceplot_", param_names[i], "_ind", ind_index, ".pdf")
    ggsave(
      filename = file.path(output_dir, output_file_name),
      width = pdf_width,
      height = pdf_height
    )
  }
}

plot_kde_pepsdi_and_semple_multiround <- function(input_dir_semple,
                                                  semple_rounds = c(2),
                                                  input_dir_pepsdi,
                                                  sample_indices_pepsdi,
                                                  thinning_pepsdi,
                                                  xlims = NULL,
                                                  ncol = NULL,
                                                  nrow = NULL,
                                                  pdf_width = 7,
                                                  pdf_height = 7,
                                                  font_size_axis_label = 10,
                                                  font_size_ticks = 10,
                                                  font_size_legend = 10) {
  settings <- rjson::fromJSON(file = file.path(input_dir_semple, "settings.json"))
  true_param <- settings$true_param
  output_dir <- file.path(input_dir_semple, "plots")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  axis_labels <- read.csv(file.path("models", settings$model_name, "param_names.csv"))

  # Read and combine samples from each SeMPLE round
  samples_list <- lapply(semple_rounds, function(r) {
    eta <- read.csv(file.path(input_dir_semple, paste0("eta_round", r, ".csv")))
    kx <- read.csv(file.path(input_dir_semple, paste0("kappa_xi_round", r, ".csv")))
    samples <- cbind(eta, kx)
    samples$method <- paste0("SeMPLE round r=", r)
    return(samples)
  })
  semple_samples <- do.call(rbind, samples_list)

  # Read PEPSDI samples
  samples_mu <- read.csv(file.path(input_dir_pepsdi, "Mean.csv"))
  samples_scale <- read.csv(file.path(input_dir_pepsdi, "Scale.csv"))
  samples_kappa_sigma <- read.csv(file.path(input_dir_pepsdi, "Kappa_sigma.csv"))

  # Match PEPSDI parameter names with SeMPLE naming
  colnames(samples_scale) <- paste0("tau", 1:ncol(samples_scale))
  colnames(samples_kappa_sigma)[which(names(samples_kappa_sigma) == "sigma1")] <- "xi1"

  # Select and thin PEPSDI samples
  samples_pepsdi <- cbind(samples_mu, samples_scale, samples_kappa_sigma)
  samples_pepsdi <- samples_pepsdi[sample_indices_pepsdi, ]
  samples_pepsdi <- samples_pepsdi[seq(1, nrow(samples_pepsdi), thinning_pepsdi), ]
  samples_pepsdi$method <- "PEPSDI"

  # Combine all samples
  all_samples <- rbind(semple_samples, samples_pepsdi)
  # Set method order: SeMPLE rounds in order, then PEPSDI
  method_levels <- c(paste0("SeMPLE round r=", semple_rounds), "PEPSDI")
  all_samples$method <- factor(all_samples$method, levels = method_levels)

  param_names <- setdiff(names(all_samples), "method")

  kde_plot_list <- lapply(1:length(param_names), function(i) {
    x <- sym(param_names[i])
    p <- ggplot(all_samples) +
      geom_density(aes(
        x = !!x,
        group = method,
        color = method
      )) +
      labs(x = TeX(as.character(axis_labels[param_names[i]])), y = "", color = "") +
      scale_color_brewer(palette = "Dark2")

    # Add true parameter value line if available
    if (!is.null(true_param) && (param_names[i] %in% names(true_param))) {
      tp_val <- as.numeric(true_param[[param_names[i]]])
      if (!is.na(tp_val)) {
        p <- p + geom_vline(
          xintercept = tp_val,
          color = "black",
          linewidth = 0.5,
          linetype = "dashed"
        )
      }
    }

    if (!is.null(xlims) && !is.null(xlims[[i]])) {
      p <- p + xlim(xlims[[i]][1], xlims[[i]][2])
    }

    p <- p + theme_bw() + theme(
      axis.title = element_text(size = font_size_axis_label),
      axis.text = element_text(size = font_size_ticks),
      legend.text = element_text(size = font_size_legend),
      plot.margin = margin(0, 0.2, 0, 0, "cm")
    )
    return(p)
  })

  pdf(file.path(output_dir, "kde_pepsdi_semple_multiround.pdf"),
    width = pdf_width,
    height = pdf_height
  )
  print(ggarrange(
    plotlist = kde_plot_list,
    common.legend = TRUE,
    ncol = ncol,
    nrow = nrow
  ))
  dev.off()
}

plot_kde_pepsdi_and_semple <- function(input_dir_semple,
                                       semple_round_index,
                                       input_dir_pepsdi,
                                       sample_indices_pepsdi,
                                       thinning_pepsdi,
                                       include_semple_prior = FALSE,
                                       xlims = NULL,
                                       ncol = NULL,
                                       nrow = NULL,
                                       pdf_width = 7,
                                       pdf_height = 7,
                                       font_size_axis_label = 10,
                                       font_size_ticks = 10,
                                       font_size_legend = 10) {
  settings <- rjson::fromJSON(file = file.path(input_dir_semple, "settings.json"))
  true_param <- settings$true_param
  output_dir <- file.path(input_dir_semple, "plots")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  axis_labels <- read.csv(file.path("models", settings$model_name, "param_names.csv"))

  # read semple posterior samples
  semple_eta_posterior <- read.csv(file.path(input_dir_semple, paste0("eta_round", semple_round_index, ".csv")))
  semple_kappa_xi_posterior <- read.csv(file.path(input_dir_semple, paste0("kappa_xi_round", semple_round_index, ".csv")))
  semple_posterior_samples <- cbind(semple_eta_posterior, semple_kappa_xi_posterior)

  # read pepsdi samples from csv
  samples_mu <- read.csv(file.path(input_dir_pepsdi, "Mean.csv"))
  samples_scale <- read.csv(file.path(input_dir_pepsdi, "Scale.csv"))
  samples_kappa_sigma <- read.csv(file.path(input_dir_pepsdi, "Kappa_sigma.csv"))

  # match parameter names with the naming in semple
  colnames(samples_scale) <- paste0("tau", 1:ncol(samples_scale))
  colnames(samples_kappa_sigma)[which(names(samples_kappa_sigma) == "sigma1")] <- "xi1"

  # select and thin out a subset of the pepsdi samples
  samples_pepsdi <- cbind(samples_mu, samples_scale, samples_kappa_sigma)
  samples_pepsdi <- samples_pepsdi[sample_indices_pepsdi, ]
  samples_pepsdi <- samples_pepsdi[seq(1, nrow(samples_pepsdi), thinning_pepsdi), ]

  if (!include_semple_prior) {
    all_samples <- rbind(semple_posterior_samples, samples_pepsdi)
    all_samples$method <- c(rep("SeMPLE", nrow(semple_posterior_samples)), rep("PEPSDI", nrow(samples_pepsdi)))
  } else {
    # read semple prior samples
    semple_eta_prior <- read.csv(file.path(input_dir_semple, "eta_round0.csv"))
    semple_kappa_xi_prior <- read.csv(file.path(input_dir_semple, "kappa_xi_round0.csv"))
    semple_prior_samples <- cbind(semple_eta_prior, semple_kappa_xi_prior)

    all_samples <- rbind(
      semple_prior_samples,
      semple_posterior_samples,
      samples_pepsdi
    )
    all_samples$method <- c(
      rep("Prior", nrow(semple_prior_samples)),
      rep("SeMPLE", nrow(semple_posterior_samples)),
      rep("PEPSDI", nrow(samples_pepsdi))
    )
    all_samples$method <- factor(all_samples$method, levels = c("Prior", "SeMPLE", "PEPSDI"))
  }
  param_names <- names(semple_posterior_samples)

  kde_plot_list <- lapply(1:length(param_names), function(i) {
    x <- sym(param_names[i])
    p <- ggplot(all_samples) +
      geom_density(aes(
        x = !!x,
        group = method,
        color = method
      )) +
      labs(x = TeX(as.character(axis_labels[param_names[i]])), y = "", color = "") +
      geom_vline(
        aes(xintercept = true_param[[param_names[i]]]),
        color = "black",
        linewidth = 0.5,
        linetype = "dashed"
      ) +
      scale_color_brewer(palette = "Dark2")
    if (!is.null(xlims[[i]])) {
      p <- p + xlim(xlims[[i]][1], xlims[[i]][2])
    }
    p <- p + theme_bw() + theme(
      axis.title = element_text(size = font_size_axis_label, face = "bold"),
      axis.text = element_text(size = font_size_ticks),
      legend.text = element_text(size = font_size_legend),
      plot.margin = margin(0, 0.3, 0, -0.3, "cm")
    )
    return(p)
  })

  pdf(
    file = file.path(output_dir, paste0("kde_comp_pepsdi.pdf")),
    width = pdf_width,
    height = pdf_height
  )
  print(ggarrange(
    plotlist = kde_plot_list,
    common.legend = T,
    ncol = ncol,
    nrow = nrow
    # align = "none"
    # hjust = 20
    # vjust = 3
  ))
  dev.off()
}
