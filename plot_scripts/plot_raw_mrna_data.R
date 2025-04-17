library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)

plot_data = function(data, ind_indices = NULL) {
  data_gathered = data %>% pivot_longer(!c(time), names_to = "ind", values_to = "value")
  if(!is.null(ind_indices)){
    plot = ggplot(data_gathered %>% dplyr::filter(ind %in% ind_indices), aes(x = time, y = value, colour = ind)) + geom_line(show.legend = F)
  } else {
    plot = ggplot(data_gathered, aes(x = time, y = value, colour = ind)) + geom_line(show.legend = F)
  }
  print(plot)
}