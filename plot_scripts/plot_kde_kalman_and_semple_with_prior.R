plot_kde_kalman_and_semple_with_prior = function(input_dir_semple,
                                                 semple_round = 2,
                                                 kalman_samples_file,
                                                 burnin_kalman,
                                                 thinning_kalman,
                                                 xlims = NULL,
                                                 ylims = NULL,
                                                 pdf_width = 7,
                                                 pdf_height = 7,
                                                 ncol = NULL,
                                                 nrow = NULL,
                                                 font_size_axis = 10,
                                                 font_size_legend = 10) {
  settings = rjson::fromJSON(file = file.path(input_dir_semple, "settings.json"))
  true_param = settings$true_param
  output_dir = file.path(input_dir_semple, "plots")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  axis_labels = read.csv(file.path("models", settings$model_name, "param_names.csv"))
  
  # read semple prior samples
  semple_eta_prior = read.csv(file.path(input_dir_semple, "eta_round0.csv"))
  semple_xi_prior = read.csv(file.path(input_dir_semple, "kappa_xi_round0.csv"))
  semple_prior_samples = cbind(semple_eta_prior, semple_xi_prior)
  
  # read semple posterior samples
  semple_eta_posterior = read.csv(file.path(input_dir_semple, paste0("eta_round", semple_round, ".csv")))
  semple_xi_posterior = read.csv(file.path(input_dir_semple, paste0("kappa_xi_round", semple_round, ".csv")))
  semple_posterior_samples = cbind(semple_eta_posterior, semple_xi_posterior)
  
  # read kalman filter samples
  kalman_samples = read.csv(kalman_samples_file)
  kalman_samples_without_burnin = kalman_samples[burnin_kalman:nrow(kalman_samples), ]
  kalman_samples_thinned = kalman_samples_without_burnin[seq(1, nrow(kalman_samples_without_burnin), thinning_kalman), ]
  
  all_samples = rbind(semple_prior_samples,
                      semple_posterior_samples,
                      kalman_samples_thinned)
  all_samples$method = c(rep("Prior", nrow(semple_prior_samples)),
                         rep("SeMPLE posterior", nrow(semple_posterior_samples)),
                         rep("MCMC posterior", nrow(kalman_samples_thinned)))
  
  all_samples$method = factor(all_samples$method, levels = c("Prior", "SeMPLE posterior", "MCMC posterior"))
  
  param_names = names(kalman_samples_thinned)
  
  kde_plot_list = lapply(1:length(param_names), function(i) {
    x = sym(param_names[i])
    p = ggplot(all_samples) +
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
    if(!is.null(xlims[[i]])){
      p = p + xlim(xlims[[i]][1], xlims[[i]][2])
    }
    if(!is.null(ylims[[i]])){
      p = p + coord_cartesian(ylim = c(ylims[[i]][1], ylims[[i]][2]))
    }
    p = p + theme_bw()
    p = p + theme(
      axis.title = element_text(size = font_size_axis, face = "bold"),
      legend.text = element_text(size = font_size_legend),
      plot.margin = margin(0, 0.5, 0, 0, "cm")
    )
    return(p)
  })
  
  pdf(file.path(output_dir, paste0("kde_comp_kalman_with_prior.pdf")),
      width = pdf_width,
      height = pdf_height)
  print(ggarrange(
    plotlist = kde_plot_list,
    common.legend = T,
    ncol = ncol,
    nrow = nrow
  ))
  dev.off()
}


plot_kde_only_noise = function(input_dir_semple,
                               kalman_samples_file,
                               burnin_kalman,
                               thinning_kalman,
                               xlims = NULL,
                               ylims = NULL,
                               pdf_width = 7,
                               pdf_height = 7,
                               ncol = NULL,
                               nrow = NULL,
                               font_size_axis = 10,
                               font_size_legend = 10) {
  
  settings = rjson::fromJSON(file = file.path(input_dir_semple, "settings.json"))
  true_param = settings$true_param
  output_dir = file.path(input_dir_semple, "plots")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  axis_labels = read.csv(file.path("models", settings$model_name, "param_names.csv"))
  
  # read semple prior samples
  semple_eta_prior = read.csv(file.path(input_dir_semple, "eta_round0.csv"))
  semple_eta_prior = semple_eta_prior["mu3"]
  semple_xi_prior = read.csv(file.path(input_dir_semple, "kappa_xi_round0.csv"))
  semple_prior_samples = cbind(semple_eta_prior, semple_xi_prior)
  semple_prior_samples["mu3xi1"] = semple_prior_samples["mu3"] * semple_prior_samples["xi1"]
  
  # read semple posterior samples
  semple_eta_posterior = read.csv(file.path(input_dir_semple, "eta_round2.csv"))
  semple_eta_posterior = semple_eta_posterior["mu3"]
  semple_xi_posterior = read.csv(file.path(input_dir_semple, "kappa_xi_round2.csv"))
  semple_posterior_samples = cbind(semple_eta_posterior, semple_xi_posterior)
  semple_posterior_samples["mu3xi1"] = semple_posterior_samples["mu3"] * semple_posterior_samples["xi1"]
  
  # read kalman filter samples
  kalman_samples = read.csv(kalman_samples_file)
  kalman_samples_without_burnin = kalman_samples[burnin_kalman:nrow(kalman_samples), ]
  kalman_samples_thinned = kalman_samples_without_burnin[seq(1, nrow(kalman_samples_without_burnin), thinning_kalman), ]
  kalman_samples_thinned = kalman_samples_thinned[c("mu3","xi1")]
  
  kalman_samples_thinned["mu3xi1"] = kalman_samples_thinned["mu3"] * kalman_samples_thinned["xi1"]
  
  axis_labels["mu3xi1"] = "$\\mu_{c_3} $\\log $\\xi"
  
  all_samples = rbind(semple_prior_samples,
                      semple_posterior_samples,
                      kalman_samples_thinned)
  all_samples$method = c(rep("Prior", nrow(semple_prior_samples)),
                         rep("SeMPLE posterior", nrow(semple_posterior_samples)),
                         rep("MCMC posterior", nrow(kalman_samples_thinned)))
  
  all_samples$method = factor(all_samples$method, levels = c("Prior", "SeMPLE posterior", "MCMC posterior"))
  
  param_names = names(all_samples)
  param_names = param_names[-length(param_names)] # remove method from parameter names
  true_param$mu3xi1 = true_param$mu3 * true_param$xi1
  
  kde_plot_list = lapply(1:length(param_names), function(i) {
    x = sym(param_names[i])
    p = ggplot(all_samples) +
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
    if(!is.null(xlims[[i]])){
      p = p + xlim(xlims[[i]][1], xlims[[i]][2])
    }
    if(!is.null(ylims[[i]])){
      p = p + coord_cartesian(ylim = c(ylims[[i]][1], ylims[[i]][2]))
    }
    p = p + theme(
      axis.title = element_text(size = font_size_axis, face = "bold"),
      legend.text = element_text(size = font_size_legend),
      plot.margin = margin(0, 0.5, 0, 0, "cm")
    )
    return(p)
  })
  
  pdf(file.path(output_dir, paste0("kde_comp_kalman_with_prior_only_sigma.pdf")),
      width = pdf_width,
      height = pdf_height)
  print(ggarrange(
    plotlist = kde_plot_list,
    common.legend = T,
    ncol = ncol,
    nrow = nrow
  ))
  dev.off()
}


plot_kde_kalman_and_semple_multiround = function(input_dir_semple,
                                                  kalman_samples_file,
                                                  burnin_kalman,
                                                  thinning_kalman,
                                                  semple_rounds = c(2),
                                                  xlims = NULL,
                                                  ylims = NULL,
                                                  pdf_width = 7,
                                                  pdf_height = 7,
                                                  ncol = NULL,
                                                  nrow = NULL,
                                                  font_size_axis = 10,
                                                  font_size_legend = 10) {
  # Read settings and axis labels
  settings = rjson::fromJSON(file = file.path(input_dir_semple, "settings.json"))
  true_param = settings$true_param
  output_dir = file.path(input_dir_semple, "plots")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  axis_labels = read.csv(file.path("models", settings$model_name, "param_names.csv"))

  # Read semple posterior samples for each requested round and label them
  semple_list = lapply(semple_rounds, function(r) {
    eta = read.csv(file.path(input_dir_semple, paste0("eta_round", r, ".csv")))
    kx = read.csv(file.path(input_dir_semple, paste0("kappa_xi_round", r, ".csv")))
    samples = cbind(eta, kx)
    samples$method = paste0("SeMPLE round r=", r)
    return(samples)
  })
  semple_all = do.call(rbind, semple_list)

  # Read kalman samples, apply burnin and thinning
  kalman_samples = read.csv(kalman_samples_file)
  kalman_samples_without_burnin = kalman_samples[burnin_kalman:nrow(kalman_samples), ]
  kalman_samples_thinned = kalman_samples_without_burnin[seq(1, nrow(kalman_samples_without_burnin), thinning_kalman), ]
  kalman_samples_thinned$method = "MCMC posterior"

  # Combine all samples (only posterior methods)
  all_samples = rbind(semple_all, kalman_samples_thinned)

  # Make method a factor with an order: semple rounds in increasing order, then MCMC
  method_levels = c(paste0("SeMPLE round r=", semple_rounds), "MCMC posterior")
  # Ensure unique and preserve order
  method_levels = unique(method_levels)
  all_samples$method = factor(all_samples$method, levels = method_levels)

  # Parameter names: take from kalman samples if present, otherwise from semple
  if (nrow(kalman_samples_thinned) > 0) {
    param_names = names(kalman_samples_thinned)
  } else {
    param_names = names(semple_all)
  }
  # Remove 'method' if present
  param_names = setdiff(param_names, "method")

  kde_plot_list = lapply(1:length(param_names), function(i) {
    pn = param_names[i]
    x = sym(pn)
    p = ggplot(all_samples) +
      geom_density(aes(
        x = !!x,
        group = method,
        color = method
      )) +
      labs(x = TeX(as.character(axis_labels[pn])), y = "", color = "") +
      scale_color_brewer(palette = "Dark2") +
      theme_bw() +
      theme(
        axis.title = element_text(size = font_size_axis, face = "bold"),
        legend.text = element_text(size = font_size_legend),
        plot.margin = margin(0, 0.5, 0, 0, "cm")
      )
    # Add true value vertical line if available
    if (!is.null(true_param) && (pn %in% names(true_param))) {
      tp_val = as.numeric(true_param[[pn]])
      if (!is.na(tp_val)) {
        p = p + geom_vline(xintercept = tp_val, color = "black", linetype = "dashed", linewidth = 0.5)
      }
    }
    if (!is.null(xlims) && !is.null(xlims[[i]])) {
      p = p + xlim(xlims[[i]][1], xlims[[i]][2])
    }
    if (!is.null(ylims) && !is.null(ylims[[i]])) {
      p = p + coord_cartesian(ylim = c(ylims[[i]][1], ylims[[i]][2]))
    }
    return(p)
  })

  pdf(file.path(output_dir, paste0("kde_comp_kalman_semple_multiround.pdf")), width = pdf_width, height = pdf_height)
  print(ggarrange(plotlist = kde_plot_list, common.legend = T, ncol = ncol, nrow = nrow))
  dev.off()
}

