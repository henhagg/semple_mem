plot_kde_pepsdi_and_semple = function(input_dir_semple,
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
                                      font_size_ticks = 10) {
  
  settings = rjson::fromJSON(file = file.path(input_dir_semple, "settings.json"))
  true_param = settings$true_param
  output_dir = file.path(input_dir_semple, "plots")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  axis_labels = read.csv(file.path("models", settings$model_name, "param_names.csv"))
  
  # read semple posterior samples 
  semple_eta_posterior = read.csv(file.path(input_dir_semple, paste0("eta_round", semple_round_index, ".csv")))
  semple_kappa_xi_posterior = read.csv(file.path(input_dir_semple, paste0("kappa_xi_round", semple_round_index, ".csv")))
  semple_posterior_samples = cbind(semple_eta_posterior, semple_kappa_xi_posterior)
  
  # read pepsdi samples from csv
  samples_mu = read.csv(file.path(input_dir_pepsdi, "Mean.csv"))
  samples_scale = read.csv(file.path(input_dir_pepsdi, "Scale.csv"))
  samples_kappa_sigma = read.csv(file.path(input_dir_pepsdi, "Kappa_sigma.csv"))
  
  # match parameter names with the naming in semple
  colnames(samples_scale) = paste0("tau", 1:ncol(samples_scale))
  colnames(samples_kappa_sigma)[which(names(samples_kappa_sigma) == "sigma1")] = "xi1"
  
  # select and thin out a subset of the pepsdi samples
  samples_pepsdi = cbind(samples_mu, samples_scale, samples_kappa_sigma)
  samples_pepsdi = samples_pepsdi[sample_indices_pepsdi, ]
  samples_pepsdi = samples_pepsdi[seq(1, nrow(samples_pepsdi), thinning_pepsdi), ]
  
  if (!include_semple_prior) {
    all_samples = rbind(semple_posterior_samples, samples_pepsdi)
    all_samples$method = c(rep("SeMPLE", nrow(semple_posterior_samples)), rep("PEPSDI", nrow(samples_pepsdi)))
  } else {
    # read semple prior samples
    semple_eta_prior = read.csv(file.path(input_dir_semple, "eta_round0.csv"))
    semple_kappa_xi_prior = read.csv(file.path(input_dir_semple, "kappa_xi_round0.csv"))
    semple_prior_samples = cbind(semple_eta_prior, semple_kappa_xi_prior)
    
    all_samples = rbind(semple_prior_samples,
                        semple_posterior_samples,
                        samples_pepsdi)
    all_samples$method = c(rep("Prior", nrow(semple_prior_samples)),
                           rep("SeMPLE", nrow(semple_posterior_samples)),
                           rep("PEPSDI", nrow(samples_pepsdi)))
    all_samples$method = factor(all_samples$method, levels = c("Prior", "SeMPLE", "PEPSDI"))
  }
  param_names = names(semple_posterior_samples)
  
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
    p = p + theme(
      axis.title = element_text(size = font_size_axis_label, face = "bold"),
      axis.text = element_text(size = font_size_ticks),
      plot.margin = margin(0, 0.2, 0, 0, "cm")
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
  ))
  dev.off()
}
