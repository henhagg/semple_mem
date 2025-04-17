library(ggplot2)
library(dplyr)

plot_runtime_comparison = function(input_dir_fe,
                                   num_individuals_fe,
                                   input_dir_re,
                                   num_individuals_re,
                                   algorithm_round_index) {
  df = data.frame(matrix(ncol = 3, nrow = 0))
  
  N_fe = length(input_dir_fe)
  for (i in 1:N_fe) {
    input_dir = input_dir_fe[i]
    elapsed_time = read.csv(file.path(input_dir, "elapsed_time.csv"))
    elapsed_time_seconds = elapsed_time[algorithm_round_index + 1, 1]
    elapsed_time_minutes = elapsed_time_seconds / 60
    df = rbind(df,
               c(
                 num_individuals_fe[i],
                 elapsed_time_minutes,
                 "With fixed-effects"
               ))
  }
  
  N_re = length(input_dir_re)
  for (i in 1:N_re) {
    input_dir = input_dir_re[i]
    elapsed_time = read.csv(file.path(input_dir, "elapsed_time.csv"))
    elapsed_time_seconds = elapsed_time[algorithm_round_index + 1, 1]
    elapsed_time_minutes = elapsed_time_seconds / 60
    df = rbind(df,
               c(
                 num_individuals_re[i],
                 elapsed_time_minutes,
                 "Only random-effects"
               ))
  }
  
  colnames(df) = c("M", "runtime", "alg_version")
  print(df)
  df$M = as.numeric(df$M)
  df$runtime = as.numeric(df$runtime)
  
  ggplot(data = df) +
    geom_line(mapping = aes(
      x = M,
      y = runtime,
      group = alg_version,
      colour = alg_version
    )) +
    geom_point(mapping = aes(
      x = M,
      y = runtime,
      group = alg_version,
      colour = alg_version
    )) +
    scale_x_continuous(breaks = num_individuals_re) +
    scale_y_continuous(breaks = round(seq(0, max(df$runtime), by = 40), 1)) +
    labs(x = "Number of individuals (M)", y = "Run-time [minutes]", colour = "Algorithm version") +
    theme_classic() +
    theme(legend.position = c(0.28, 0.8),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10))
  
  ggsave("results/runtime_comparison.pdf", width = 3.8, height = 2.2)
  
  # plot ratio
  ratio = df$runtime[1:N_fe]/df$runtime[(N_fe+1):(N_fe+N_re)]
  M = df$M[1:N_fe]
  df_ratio = data.frame(M, ratio)
  ggplot(data = df_ratio) +
    geom_line(mapping = aes(
      x = M,
      y = ratio,
    )) +
    geom_point(mapping = aes(
      x = M,
      y = ratio,
    )) +
    scale_x_continuous(breaks = num_individuals_re) +
    # scale_y_continuous(breaks = round(seq(0, max(df$runtime), by = 40), 1)) +
    labs(x = "Number of individuals (M)", y = "", title = "Run-time ratio fixed-effects/random-effects") +
    theme_classic() +
    theme(legend.position = c(0.28, 0.8),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.title = element_text(size=12))
  
  ggsave("results/runtime_ratio.pdf", width = 3.8, height = 2.2)
}

plot_runtime_comparison(
  input_dir_fe = c(
    "results/mrna_indep_prior/10ind/K10",
    "results/mrna_indep_prior/40ind/K10",
    "results/mrna_indep_prior/100ind/K10",
    "results/mrna_indep_prior/200ind/K10",
    "results/mrna_indep_prior/400ind/K10"
  ),
  num_individuals_fe = c(10, 40, 100, 200, 400),
  input_dir_re = c(
    "results/mrna_indep_prior_only_individual_param/10ind/K10",
    "results/mrna_indep_prior_only_individual_param/40ind/K10",
    "results/mrna_indep_prior_only_individual_param/100ind/K10",
    "results/mrna_indep_prior_only_individual_param/200ind/K10",
    "results/mrna_indep_prior_only_individual_param/400ind/K10"
  ),
  num_individuals_re = c(10, 40, 100, 200, 400),
  algorithm_round_index = 2
)