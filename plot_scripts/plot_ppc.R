plot_ppc = function(input_dir,
                    round_index,
                    pdf_width = 5,
                    pdf_height = 3,
                    font_size_axis = 10,
                    font_size_ticks = 10) {
  ppc_quantiles = read.csv(file.path(input_dir, paste0("ppc_round", round_index, ".csv")))
  settings = rjson::fromJSON(file = file.path(input_dir, "settings.json"))
  observed_data = read.csv(
    file.path(
      "models",
      settings$model_name,
      settings$observation_name,
      "observations.csv"
    )
  )
  
  observed_data_gathered = observed_data %>% pivot_longer(!c(time), names_to = "Observation", values_to = "value")
  ppc_quantiles_gathered = ppc_quantiles %>% pivot_longer(!c(time), names_to = "quantile", values_to = "value")
  
  ggplot() +
    geom_ribbon(
      data = ppc_quantiles,
      mapping = aes(x = time, ymin = qu0025, ymax = qu0975),
      alpha = 0.2
    ) +
    geom_line(data = observed_data_gathered,
              mapping = aes(x = time, y = value, colour = Observation)) +
    labs(x = "Time", y = "") +
    theme(
      legend.position = "none",
      axis.title = element_text(size = font_size_axis),
      axis.text = element_text(size = font_size_ticks)
    )
  
  ggsave(
    filename = file.path(input_dir, "plots", paste0("ppc_round", round_index, ".pdf")),
    width = pdf_width,
    height = pdf_height
  )
}

plot_ppc_individual = function(input_dir,
                               round_index,
                               individual_index,
                               pdf_width = 5,
                               pdf_height = 3,
                               xlab = NULL,
                               ylab = NULL) {
  ppc_quantiles = read.csv(file.path(input_dir, paste0("ppc_individual_round", round_index, ".csv")))
  settings = rjson::fromJSON(file = file.path(input_dir, "settings.json"))
  observed_data = read.csv(
    file.path(
      "models",
      settings$model_name,
      settings$observation_name,
      "observations.csv"
    )
  )
  
  obs_column_name = sym(paste0("obs", individual_index))
  
  p = ggplot() +
    geom_ribbon(
      data = ppc_quantiles %>% filter(id == individual_index),
      mapping = aes(x = time, ymin = qu0025, ymax = qu0975),
      alpha = 0.2
    ) +
    geom_line(data = observed_data,
              mapping = aes(x = time, y = !!obs_column_name))
  if(!is.null(ylab)){
    p = p + ylab(ylab)
  } else {
    p = p + ylab("")
  }
  if(!is.null(xlab)){
    p = p + xlab(xlab)
  }
  # ggtitle(paste0("Posterior predictive check observation ", individual_index, ", algorithm round ", round_index))
  
  ggsave(
    filename = file.path(
      input_dir,
      "plots",
      paste0("ppc_ind", individual_index, "_round", round_index, ".pdf")
    ),
    width = pdf_width,
    height = pdf_height
  )
}
