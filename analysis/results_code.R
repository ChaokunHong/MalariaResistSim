library(AMRSpreadModel)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(cowplot)
library(tidyr)
library(broom)
library(scales)
library(mblm)
library(MASS)
library(mgcv)

# Figure2. The change in resistance over time for one ft and EIR for different values of rTR_true.
EIR_range <- c(1, 10, 50)
ft_range <- c(0.2, 0.4, 0.4)

range_results <- range_model(
  param_ranges = list(EIR = EIR_range, ft = ft_range, rTR_true = c(0.02, 0.05, 0.1, 0.15, 0.2)),
  use_eir_ft = TRUE,
  time_scale = "year",
  run_time = 7,
  time_points = c(0:12),
  ton = 2, toff = 999999, res_time = 1,
  init_res = 0.1, day0_res = 0.01
)
range_results <- as.data.frame(range_results)

plot_list <- list()
for (eir in EIR_range) {
  for (ft in ft_range) {
    subset_results <- range_results[which(range_results$EIR == eir & range_results$ft == ft), ]
    p <- plot_range(
      subset_results,
      x_var = "t",
      y_var = "prevalence_res",
      color_var = "rTR_true",
      title = paste("EIR =", eir, ", ft =", ft)
    ) + theme(legend.position = "none") +
      xlab("Time (Years)") +
      ylab("Prevalence (resistant)") +
      scale_y_continuous(breaks = c (0, 0.1, 0.2, 0.3),
                         limits = c(0, 0.3))
    plot_list[[paste("EIR", eir, "ft", ft)]] <- p
  }
}

g_legend <- function(a.gplot) {
  tmp <- ggplotGrob(a.gplot)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend_plot <- ggplot(range_results, aes(x = t, y = prevalence_res, color = as.factor(rTR_true))) +
  geom_line() +
  labs(color = "rTR_true") +
  theme_minimal() +
  theme(legend.position = "bottom")

legend <- g_legend(legend_plot)
plot1 <- do.call(grid.arrange, c(plot_list, ncol = 3, list(bottom = legend)))
































# Figure 3. Comparisons of resistance over time for different ft and EIR and TR_true combinations
EIR_range2 <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 30, 50, 100, 200)
ft_range2 <- seq(0.1, 0.9, 0.1)

range_results <- range_model(
  param_ranges = list(EIR = EIR_range2, ft = ft_range2, rTR_true = c(0.02, 0.05, 0.1, 0.15, 0.2)),
  use_eir_ft = TRUE,
  run_time = 3650,
  time_points = c(3650),
  ton = 200, toff = 999999, res_time =100,
  init_res = 0.1, day0_res = 0.01
)

range_results$rTR_true <- as.factor(range_results$rTR_true)

combined_plot <- ggplot(range_results, aes(x = ft, y = prevalence_res, color = rTR_true)) +
  geom_line() +
  facet_wrap(~ EIR, scales = "fixed") +
  scale_y_continuous(breaks = seq(0.1, 1, length.out = 5),
                     limits = c(0, 1)) +
  labs(color = "rTR_true", x = "ft", y = "Prevalence (resistant)", title = "Comparisons of resistance over time for different ft and EIR and TR_true combinations") +
  theme_minimal()

print(combined_plot)














# Supplementary Figures: Model Checking figures
EIR_range3 <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 30, 50, 100, 200)
ft_range3 <- seq(0.1, 0.9, 0.1)
start_results <- starting_params( EIR = EIR_range3, ft = ft_range3)
start_results$prevalence <- start_results$A + start_results$D + start_results$T

ggplot(start_results, aes(x = EIR, y = prevalence, color = as.factor(ft), group = as.factor(ft))) +
  geom_line() +
  labs(title = "Malaria Prevalence vs EIR",
       x = "EIR",
       y = "Prevalence",
       color = "ft") +
  ylim(0, 1) +
  theme_minimal()



















# Uganda study
library(ggplot2)
library(dplyr)

first_line <- readLines("/Users/hongchaokun/Documents/BMR/Marial Drug Resistance Modenling/data/Uganda_allele_frequency.txt", n = 1)
print(first_line)

column_names <- c("District", "year", "Locus", "n", "min_year", "adj_year", "nobs", "lrsmed", "freq", "wt_med")

data <- read.table("/Users/hongchaokun/Documents/BMR/Marial Drug Resistance Modenling/data/Uganda_allele_frequency.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)

print(colnames(data))

filtered_data <- data %>% filter(Locus == "K13")

filtered_data <- filtered_data %>% mutate(freq = freq * 100)

filtered_data$freq_initial <- filtered_data$freq / 100

ggplot(filtered_data, aes(x = year, y = freq, size = n)) +
  geom_point(aes(color = "K13"), alpha = 0.6) +
  geom_smooth(method = "lm", aes(color = "K13"), se = TRUE) +
  facet_wrap(~ District) +
  scale_size_continuous(range = c(1, 4), breaks = c(20, 40, 60, 80)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
  labs(x = "Year", y = "Prevalence of Resistance (%)", color = "Mutations", size = "Sample Size") +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "lightgrey"),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "right")





# One District
library(dplyr)

filtered_data <- data %>% filter(Locus == "K13" & District == "Agago")

plot2 <- ggplot(filtered_data, aes(x = year, y = freq, size = n)) +
  geom_point(aes(color = "K13"), alpha = 0.6) +
  geom_smooth(method = "lm", aes(color = "K13"), se = TRUE) +
  facet_wrap(~ District) +
  scale_size_continuous(range = c(1, 4), breaks = c(20, 40, 60, 80)) +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.2), limits = c(0, 0.6)) +
  labs(x = "Year", y = "Prevalence of Resistance", color = "Mutations", size = "Sample Size") +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "lightgrey"),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "right")
print(plot2)


plot2_smaller <- ggdraw() +
  draw_plot(plot2, width = 0.4, height = 0.8)
combined_plot <- plot_grid(plot1, plot2_smaller, labels = c("A", "B"), ncol = 1, rel_heights = c(2, 0.4))
print(combined_plot)





































########################################################################################################################################################################################
### Compare the model results with the Uganda data
# Load required libraries
library(dplyr)
library(ggplot2)
library(broom)
library(scales)

# Load data
uganda_raw <- read.table("/Users/hongchaokun/Documents/BMR/Marial Drug Resistance Modenling/data/Uganda_allele_frequency.txt", header = TRUE, stringsAsFactors = FALSE)

# Calculate slope function
calculate_slope <- function(x, y) {
  if(length(y) < 2) return(NA)
  lm_model <- lm(y ~ x)
  return(coef(lm_model)[2])
}

# Process Uganda data to calculate the slope for each District
uganda_data <- uganda_raw %>%
  dplyr::filter(Locus == "K13") %>%
  dplyr::group_by(District) %>%
  dplyr::summarise(actual_slope = calculate_slope(year, freq),
                   start_year = min(year),
                   end_year = max(year))

# Extract initial resistance values for each District
initial_resistances <- uganda_raw %>%
  dplyr::filter(Locus == "K13") %>%
  dplyr::group_by(District) %>%
  dplyr::arrange(year) %>%
  dplyr::summarise(day0_res = ifelse(first(freq) == 0, 0.005, first(freq)),
                   start_year = first(year))

# Assume range_results is already generated

# Process model results and adjust time axis
model_slopes <- range_results %>%
  dplyr::mutate(year = t + 2015) %>%
  dplyr::group_by(EIR, ft, rTR_true) %>%
  dplyr::summarise(model_slope = calculate_slope(year, prevalence_res), .groups = "drop")

# Compare slopes
slope_comparison <- model_slopes %>%
  dplyr::cross_join(uganda_data) %>%
  dplyr::mutate(slope_diff = abs(model_slope - actual_slope)) %>%
  dplyr::arrange(District, slope_diff)

# Select the closest N model results for each District
N <- 10
top_matches <- slope_comparison %>%
  dplyr::group_by(District) %>%
  dplyr::slice_head(n = N)

for (district in unique(top_matches$District)) {
  district_data <- initial_resistances %>%
    dplyr::filter(District == district)
  start_year <- district_data$start_year
  end_year <- uganda_data %>%
    dplyr::filter(District == district) %>%
    dplyr::pull(end_year)

  model_params <- top_matches %>%
    dplyr::filter(District == district) %>%
    dplyr::select(EIR, ft, rTR_true) %>%
    distinct()

  plot_data <- range_results %>%
    dplyr::semi_join(model_params, by = c("EIR", "ft", "rTR_true")) %>%
    dplyr::mutate(
      year = t + start_year - 1,
      model_label = sprintf("EIR=%.0f, ft=%.1f, rTR=%.2f", EIR, ft, rTR_true)
    ) %>%
    dplyr::filter(year >= start_year, year <= end_year)

  actual_data <- uganda_raw %>%
    dplyr::filter(District == district, Locus == "K13") %>%
    dplyr::select(year, freq) %>%
    dplyr::mutate(
      prevalence_res = freq,
      model_label = "Actual Data"
    )

  combined_data <- bind_rows(
    plot_data %>% dplyr::select(year, prevalence_res, model_label),
    actual_data %>% dplyr::select(year, prevalence_res, model_label)
  )

  fitted_models <- plot_data %>%
    dplyr::group_by(model_label) %>%
    dplyr::do(broom::tidy(lm(prevalence_res ~ year, data = .))) %>%
    dplyr::filter(term == "year") %>%
    dplyr::select(model_label, slope = estimate) %>%
    dplyr::ungroup()

  actual_fit <- broom::tidy(lm(freq ~ year, data = actual_data)) %>%
    dplyr::filter(term == "year") %>%
    dplyr::pull(estimate)

  fitted_models <- fitted_models %>%
    dplyr::mutate(label_with_slope = sprintf("%s (%.4f)", model_label, slope)) %>%
    dplyr::add_row(model_label = "Actual Data", slope = actual_fit, label_with_slope = sprintf("Actual Data (%.4f)", actual_fit))

  n_colors <- nrow(fitted_models) - 1
  color_palette <- scales::hue_pal()(n_colors)
  color_palette <- c(color_palette, "black")

  # Set colors and linetypes
  color_mapping <- setNames(color_palette, fitted_models$model_label)
  linetype_mapping <- setNames(c(rep("solid", n_colors), "dashed"), fitted_models$model_label)

  p <- ggplot(combined_data, aes(x = year, y = prevalence_res, color = model_label, linetype = model_label)) +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x, linewidth = 1, show.legend = TRUE) +
    scale_color_manual(name = "Model Parameters and Actual Data (Slope)",
                       values = color_mapping) +
    scale_linetype_manual(name = "Model Parameters and Actual Data (Slope)",
                          values = linetype_mapping) +
    scale_x_continuous(breaks = seq(start_year, end_year, by = 1)) +
    coord_cartesian(xlim = c(start_year, end_year)) +
    labs(title = paste("Trend Comparison for", district),
         subtitle = paste("Data range:", start_year, "-", end_year),
         x = "Year",
         y = "Prevalence / Frequency") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(t = 10, b = 10),
      legend.key.width = unit(1.5, "cm"),
      legend.text = element_text(size = 8),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(ncol = 3), linetype = guide_legend(ncol = 3))

  print(p)
}


########################################################################################################################################################################################




















































































########################################################################################################################################################################################
### Compare the model results with the Uganda data (scatter plot)
calculate_slope <- function(x, y) {
  if (length(y) < 2) return(NA)
  lm_model <- lm(y ~ x)
  return(coef(lm_model)[2])
}

# Load data
uganda_raw <- read.table("/Users/hongchaokun/Documents/BMR/Marial Drug Resistance Modenling/data/Uganda_allele_frequency.txt", header = TRUE, stringsAsFactors = FALSE)

# Process Uganda data to calculate the slope for each District
uganda_data <- uganda_raw %>%
  dplyr::filter(Locus == "K13") %>%
  dplyr::group_by(District) %>%
  dplyr::summarise(actual_slope = calculate_slope(year, freq),
                   start_year = min(year),
                   end_year = max(year))

# Extract initial resistance values for each District
initial_resistances <- uganda_raw %>%
  dplyr::filter(Locus == "K13") %>%
  dplyr::group_by(District) %>%
  dplyr::arrange(year) %>%
  dplyr::summarise(day0_res = ifelse(first(freq) == 0, 0.005, first(freq)),
                   start_year = first(year))

# Create an empty list to store plots
plots <- list()

for (district in unique(uganda_data$District)) {
  # Get initial resistance value and start year
  district_data <- initial_resistances %>%
    dplyr::filter(District == district)
  day0_res <- district_data$day0_res
  start_year <- district_data$start_year

  # Get end year
  end_year <- uganda_data %>%
    dplyr::filter(District == district) %>%
    dplyr::pull(end_year)

  # Calculate the number of years to model
  years_to_model <- end_year - start_year + 1

  # Generate model results
  EIR_range <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300)
  ft_range <- seq(0.1, 0.9, 0.1)
  range_results <- range_model(
    param_ranges = list(EIR = EIR_range, ft = ft_range, rTR_true = c(0.02, 0.05, 0.1, 0.15, 0.2)),
    use_eir_ft = TRUE,
    time_scale = "year",
    run_time = years_to_model,
    time_points = c(0:years_to_model),
    ton = 2, toff = 8, res_time = 0,
    init_res = 0, day0_res = day0_res
  )

  # Check if range_results is empty or NULL
  if (is.null(range_results) || nrow(range_results) == 0) {
    warning(paste("No successful model runs for district:", district))
    next
  }

  range_results <- as.data.frame(range_results)

  # Process model results and adjust time axis
  model_slopes <- range_results %>%
    dplyr::mutate(year = t + start_year - 1) %>%
    dplyr::filter(year >= start_year) %>%
    dplyr::group_by(EIR, ft, rTR_true) %>%
    dplyr::summarise(model_slope = calculate_slope(year, prevalence_res), .groups = "drop")

  # Compare slopes
  slope_comparison <- model_slopes %>%
    dplyr::cross_join(uganda_data %>% dplyr::filter(District == district)) %>%
    dplyr::mutate(slope_diff = abs(model_slope - actual_slope)) %>%
    dplyr::arrange(slope_diff)

  # Select the closest N model results for each District
  N <- 20
  top_matches <- slope_comparison %>%
    dplyr::slice_head(n = N)

  model_params <- top_matches %>%
    dplyr::select(EIR, ft, rTR_true) %>%
    dplyr::distinct()

  plot_data <- range_results %>%
    dplyr::semi_join(model_params, by = c("EIR", "ft", "rTR_true")) %>%
    dplyr::mutate(
      year = t + start_year - 1,
      model_label = sprintf("EIR=%.0f, ft=%.1f, rTR=%.2f", EIR, ft, rTR_true)
    ) %>%
    dplyr::filter(year >= start_year)

  actual_data <- uganda_raw %>%
    dplyr::filter(District == district, Locus == "K13") %>%
    dplyr::select(year, freq) %>%
    dplyr::mutate(
      prevalence_res = freq,
      model_label = "Actual Data"
    )

  combined_data <- dplyr::bind_rows(
    plot_data %>% dplyr::select(year, prevalence_res, model_label),
    actual_data %>% dplyr::select(year, prevalence_res, model_label)
  )

  fitted_models <- plot_data %>%
    dplyr::group_by(model_label) %>%
    dplyr::do(broom::tidy(lm(prevalence_res ~ year, data = .))) %>%
    dplyr::filter(term == "year") %>%
    dplyr::select(model_label, slope = estimate) %>%
    dplyr::ungroup()

  actual_fit <- broom::tidy(lm(freq ~ year, data = actual_data)) %>%
    dplyr::filter(term == "year") %>%
    dplyr::pull(estimate)

  fitted_models <- fitted_models %>%
    dplyr::mutate(label_with_slope = sprintf("%s (%.4f)", model_label, slope)) %>%
    dplyr::add_row(model_label = "Actual Data", slope = actual_fit, label_with_slope = sprintf("Actual Data (%.4f)", actual_fit))

  n_colors <- nrow(fitted_models) - 1
  color_palette <- scales::hue_pal()(n_colors)
  color_palette <- c(color_palette, "black")

  # Update combined_data with label_with_slope
  combined_data <- combined_data %>%
    dplyr::left_join(fitted_models %>% dplyr::select(model_label, label_with_slope), by = "model_label") %>%
    dplyr::mutate(label_with_slope = factor(label_with_slope, levels = fitted_models$label_with_slope))

  # Set colors and linetypes
  color_mapping <- setNames(color_palette, levels(combined_data$label_with_slope))
  linetype_mapping <- setNames(c(rep("solid", n_colors), "dashed"), levels(combined_data$label_with_slope))

  # Plot
  p <- ggplot() +
    geom_line(data = combined_data %>% dplyr::filter(model_label != "Actual Data"), aes(x = year, y = prevalence_res, color = label_with_slope), linewidth = 1) +
    geom_point(data = actual_data, aes(x = year, y = prevalence_res), size = 2, color = "black") +
    geom_smooth(data = actual_data, aes(x = year, y = prevalence_res), color = "black",
                method = "lm", se = TRUE, formula = y ~ x, linewidth = 1, linetype = "dashed") +
    scale_color_manual(name = "Model Parameters and Actual Data (Slope)",
                       values = color_mapping) +
    scale_x_continuous(breaks = seq(start_year, end_year, by = 1)) +
    coord_cartesian(xlim = c(start_year, end_year)) +
    labs(title = paste("Trend Comparison for", district),
         subtitle = paste("Data range:", start_year, "-", end_year),
         x = "Year",
         y = "Prevalence of Resistance") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(t = 10, b = 10),
      legend.key.width = unit(1.5, "cm"),
      legend.text = element_text(size = 8),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(ncol = 5))

  # Add plot to the list
  plots[[district]] <- p

  print(p)
}

combined_plot <- marrangeGrob(grobs = plots, ncol = 2, nrow = 2)
ggsave("/Users/hongchaokun/Desktop/Uganda/trend_scatter_combined.pdf", combined_plot, width = 40, height = 20)


























































########################################################################################################################################################################################
### Compare the model results with the Uganda data (slope plot)
calculate_slope <- function(x, y) {
  if (length(y) < 2) return(NA)
  lm_model <- lm(y ~ x)
  return(coef(lm_model)[2])
}

# Load data
uganda_raw <- read.table("/Users/hongchaokun/Documents/BMR/Marial Drug Resistance Modenling/data/Uganda_allele_frequency.txt", header = TRUE, stringsAsFactors = FALSE)

# Process Uganda data to calculate the slope for each District
uganda_data <- uganda_raw %>%
  dplyr::filter(Locus == "K13") %>%
  dplyr::group_by(District) %>%
  dplyr::summarise(actual_slope = calculate_slope(year, freq),
                   start_year = min(year),
                   end_year = max(year))

# Extract initial resistance values for each District
initial_resistances <- uganda_raw %>%
  dplyr::filter(Locus == "K13") %>%
  dplyr::group_by(District) %>%
  dplyr::arrange(year) %>%
  dplyr::summarise(day0_res = ifelse(first(freq) == 0, 0.005, first(freq)),
                   start_year = first(year))

# Create an empty list to store plots
plots_slope <- list()

for (district in unique(uganda_data$District)) {
  # Get initial resistance value and start year
  district_data <- initial_resistances %>%
    dplyr::filter(District == district)
  day0_res <- district_data$day0_res
  start_year <- district_data$start_year

  # Get end year
  end_year <- uganda_data %>%
    dplyr::filter(District == district) %>%
    dplyr::pull(end_year)

  # Calculate the number of years to model
  years_to_model <- end_year - start_year + 1

  # Generate model results
  EIR_range <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300)
  ft_range <- seq(0.1, 0.9, 0.1)
  range_results <- range_model(
    param_ranges = list(EIR = EIR_range, ft = ft_range, rTR_true = c(0.02, 0.05, 0.1, 0.15, 0.2)),
    use_eir_ft = TRUE,
    time_scale = "year",
    run_time = years_to_model,
    time_points = c(0:years_to_model),
    ton = 1, toff = 999999, res_time = 0,
    init_res = 0, day0_res = day0_res
  )

  # Check if range_results is empty or NULL
  if (is.null(range_results) || nrow(range_results) == 0) {
    warning(paste("No successful model runs for district:", district))
    next
  }

  range_results <- as.data.frame(range_results)

  # Process model results and adjust time axis
  model_slopes <- range_results %>%
    dplyr::mutate(year = t + start_year - 1) %>%
    dplyr::filter(year >= start_year) %>%
    dplyr::group_by(EIR, ft, rTR_true) %>%
    dplyr::summarise(model_slope = calculate_slope(year, prevalence_res), .groups = "drop")

  # Compare slopes
  slope_comparison <- model_slopes %>%
    dplyr::cross_join(uganda_data %>% dplyr::filter(District == district)) %>%
    dplyr::mutate(slope_diff = abs(model_slope - actual_slope)) %>%
    dplyr::arrange(slope_diff)

  # Select the closest N model results for each District
  N <- 20
  top_matches <- slope_comparison %>%
    dplyr::slice_head(n = N)

  model_params <- top_matches %>%
    dplyr::select(EIR, ft, rTR_true) %>%
    dplyr::distinct()

  plot_data <- range_results %>%
    dplyr::semi_join(model_params, by = c("EIR", "ft", "rTR_true")) %>%
    dplyr::mutate(
      year = t + start_year - 1,
      model_label = sprintf("EIR=%.0f, ft=%.1f, rTR=%.2f", EIR, ft, rTR_true)
    ) %>%
    dplyr::filter(year >= start_year)

  actual_data <- uganda_raw %>%
    dplyr::filter(District == district, Locus == "K13") %>%
    dplyr::select(year, freq) %>%
    dplyr::mutate(
      prevalence_res = freq,
      model_label = "Actual Data"
    )

  combined_data <- dplyr::bind_rows(
    plot_data %>% dplyr::select(year, prevalence_res, model_label),
    actual_data %>% dplyr::select(year, prevalence_res, model_label)
  )

  fitted_models <- plot_data %>%
    dplyr::group_by(model_label) %>%
    dplyr::do(broom::tidy(lm(prevalence_res ~ year, data = .))) %>%
    dplyr::filter(term == "year") %>%
    dplyr::select(model_label, slope = estimate) %>%
    dplyr::ungroup()

  actual_fit <- broom::tidy(lm(freq ~ year, data = actual_data)) %>%
    dplyr::filter(term == "year") %>%
    dplyr::pull(estimate)

  fitted_models <- fitted_models %>%
    dplyr::mutate(label_with_slope = sprintf("%s (%.4f)", model_label, slope)) %>%
    dplyr::add_row(model_label = "Actual Data", slope = actual_fit, label_with_slope = sprintf("Actual Data (%.4f)", actual_fit))

  n_colors <- nrow(fitted_models) - 1
  color_palette <- scales::hue_pal()(n_colors)
  color_palette <- c(color_palette, "black")

  # Update combined_data with label_with_slope
  combined_data <- combined_data %>%
    dplyr::left_join(fitted_models %>% dplyr::select(model_label, label_with_slope), by = "model_label") %>%
    dplyr::mutate(label_with_slope = factor(label_with_slope, levels = fitted_models$label_with_slope))

  # Set colors and linetypes
  color_mapping <- setNames(color_palette, levels(combined_data$label_with_slope))
  linetype_mapping <- setNames(c(rep("solid", n_colors), "dashed"), levels(combined_data$label_with_slope))

  # Plot
  p_slope <- ggplot(combined_data, aes(x = year, y = prevalence_res, color = label_with_slope, linetype = label_with_slope)) +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x, linewidth = 1) +
    geom_point(data = combined_data %>% filter(model_label == "Actual Data"),
               aes(x = year, y = prevalence_res),
               size = 2, color = "black") +
    scale_color_manual(name = "Model Parameters and Actual Data (Slope)",
                       values = color_mapping) +
    scale_linetype_manual(name = "Model Parameters and Actual Data (Slope)",
                          values = linetype_mapping) +
    scale_x_continuous(breaks = seq(start_year, end_year, by = 1)) +
    coord_cartesian(xlim = c(start_year, end_year)) +
    labs(title = paste("Trend Comparison for", district),
         subtitle = paste("Data range:", start_year, "-", end_year),
         x = "Year",
         y = "Prevalence / Frequency") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(t = 10, b = 10),
      legend.key.width = unit(1.5, "cm"),
      legend.text = element_text(size = 8),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(ncol = 3), linetype = guide_legend(ncol = 3))

  # Add plot to the list
  plots_slope[[district]] <- p_slope

  print(p_slope)
}

combined_plot <- marrangeGrob(grobs = plots_slope, ncol = 2, nrow = 2)
ggsave("/Users/hongchaokun/Desktop/Uganda/trend_linear_combined.pdf", combined_plot, width = 40, height = 20)
























































































































########################################################################################################################################################################################
#### Log
# Modified slope calculation function using log method
calculate_slope <- function(x, y) {
  if (length(y) < 2) return(NA)
  # Add a small constant to y to avoid log(0)
  y_adjusted <- y + 1e-10
  log_model <- lm(log(y_adjusted) ~ x)
  # Return the exponential of the coefficient minus 1 for percentage change
  return(exp(coef(log_model)[2]) - 1)
}

# Load data
uganda_raw <- read.table("/Users/hongchaokun/Documents/BMR/Marial Drug Resistance Modenling/data/Uganda_allele_frequency.txt", header = TRUE, stringsAsFactors = FALSE)

# Process Uganda data to calculate the slope for each District
uganda_data <- uganda_raw %>%
  dplyr::filter(Locus == "K13") %>%
  dplyr::group_by(District) %>%
  dplyr::summarise(actual_slope = calculate_slope(year, freq),
                   start_year = min(year),
                   end_year = max(year))

# Extract initial resistance values for each District
initial_resistances <- uganda_raw %>%
  dplyr::filter(Locus == "K13") %>%
  dplyr::group_by(District) %>%
  dplyr::arrange(year) %>%
  dplyr::summarise(day0_res = ifelse(first(freq) == 0, 0.005, first(freq)),
                   start_year = first(year))

# Create an empty list to store plots
plots <- list()

for (district in unique(uganda_data$District)) {
  # Get initial resistance value and start year
  district_data <- initial_resistances %>%
    dplyr::filter(District == district)
  day0_res <- district_data$day0_res
  start_year <- district_data$start_year

  # Get end year
  end_year <- uganda_data %>%
    dplyr::filter(District == district) %>%
    dplyr::pull(end_year)

  # Calculate the number of years to model
  years_to_model <- end_year - start_year + 1

  # Generate model results
  EIR_range <- c(280, 290, 300)
  ft_range <- seq(0.7, 0.9, 0.1)
  range_results <- range_model(
    param_ranges = list(EIR = EIR_range, ft = ft_range, rTR_true = c(0.02, 0.05, 0.1, 0.15, 0.2)),
    use_eir_ft = TRUE,
    time_scale = "year",
    run_time = years_to_model,
    time_points = c(0:years_to_model),
    ton = 1, toff = 999999, res_time = 0,
    init_res = 0, day0_res = day0_res
  )

  # Check if range_results is empty or NULL
  if (is.null(range_results) || nrow(range_results) == 0) {
    warning(paste("No successful model runs for district:", district))
    next
  }

  range_results <- as.data.frame(range_results)

  # Process model results and adjust time axis
  model_slopes <- range_results %>%
    dplyr::mutate(year = t + start_year - 1) %>%
    dplyr::filter(year >= start_year) %>%
    dplyr::group_by(EIR, ft, rTR_true) %>%
    dplyr::summarise(model_slope = calculate_slope(year, prevalence_res), .groups = "drop")

  # Compare slopes
  slope_comparison <- model_slopes %>%
    dplyr::cross_join(uganda_data %>% dplyr::filter(District == district)) %>%
    dplyr::mutate(slope_diff = abs(model_slope - actual_slope)) %>%
    dplyr::arrange(slope_diff)

  # Select the closest N model results for each District
  N <- 20
  top_matches <- slope_comparison %>%
    dplyr::slice_head(n = N)

  model_params <- top_matches %>%
    dplyr::select(EIR, ft, rTR_true) %>%
    dplyr::distinct()

  plot_data <- range_results %>%
    dplyr::semi_join(model_params, by = c("EIR", "ft", "rTR_true")) %>%
    dplyr::mutate(
      year = t + start_year - 1,
      model_label = sprintf("EIR=%.0f, ft=%.1f, rTR=%.2f", EIR, ft, rTR_true)
    ) %>%
    dplyr::filter(year >= start_year)

  actual_data <- uganda_raw %>%
    dplyr::filter(District == district, Locus == "K13") %>%
    dplyr::select(year, freq) %>%
    dplyr::mutate(
      prevalence_res = freq,
      model_label = "Actual Data"
    )

  combined_data <- dplyr::bind_rows(
    plot_data %>% dplyr::select(year, prevalence_res, model_label),
    actual_data %>% dplyr::select(year, prevalence_res, model_label)
  )

  fitted_models <- plot_data %>%
    dplyr::group_by(model_label) %>%
    dplyr::do(broom::tidy(lm(prevalence_res ~ year, data = .))) %>%
    dplyr::filter(term == "year") %>%
    dplyr::select(model_label, slope = estimate) %>%
    dplyr::ungroup()

  actual_fit <- broom::tidy(lm(freq ~ year, data = actual_data)) %>%
    dplyr::filter(term == "year") %>%
    dplyr::pull(estimate)

  fitted_models <- fitted_models %>%
    dplyr::mutate(label_with_slope = sprintf("%s (%.4f)", model_label, slope)) %>%
    dplyr::add_row(model_label = "Actual Data", slope = actual_fit, label_with_slope = sprintf("Actual Data (%.4f)", actual_fit))

  n_colors <- nrow(fitted_models) - 1
  color_palette <- scales::hue_pal()(n_colors)
  color_palette <- c(color_palette, "black")

  # Update combined_data with label_with_slope
  combined_data <- combined_data %>%
    dplyr::left_join(fitted_models %>% dplyr::select(model_label, label_with_slope), by = "model_label") %>%
    dplyr::mutate(label_with_slope = factor(label_with_slope, levels = fitted_models$label_with_slope))

  # Set colors and linetypes
  color_mapping <- setNames(color_palette, levels(combined_data$label_with_slope))
  linetype_mapping <- setNames(c(rep("solid", n_colors), "dashed"), levels(combined_data$label_with_slope))

  # Plot
  p <- ggplot() +
    geom_line(data = combined_data %>% dplyr::filter(model_label != "Actual Data"), aes(x = year, y = prevalence_res, color = label_with_slope), linewidth = 1) +
    geom_point(data = actual_data, aes(x = year, y = prevalence_res), size = 2, color = "black") +
    geom_smooth(data = actual_data, aes(x = year, y = prevalence_res), color = "black",
                method = "lm", se = TRUE, formula = y ~ x, linewidth = 1, linetype = "dashed") +
    scale_color_manual(name = "Model Parameters and Actual Data (Slope)",
                       values = color_mapping) +
    scale_x_continuous(breaks = seq(start_year, end_year, by = 1)) +
    coord_cartesian(xlim = c(start_year, end_year)) +
    labs(title = paste("Trend Comparison for", district),
         subtitle = paste("Data range:", start_year, "-", end_year),
         x = "Year",
         y = "Prevalence of Resistance") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(t = 10, b = 10),
      legend.key.width = unit(1.5, "cm"),
      legend.text = element_text(size = 8),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(ncol = 5))

  # Add plot to the list
  plots[[district]] <- p

  print(p)
}
