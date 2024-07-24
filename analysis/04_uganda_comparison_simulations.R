library(tidyverse)

# -------------------------------------------------------- o
# 1. First lets work out suitable EIR and ft ranges --------
# -------------------------------------------------------- o

# Grab our covariate parameter ranges
chaokun_dat <-  read.csv("https://raw.githubusercontent.com/OJWatson/hrpup/main/analysis/data_derived/chaokun_uganda.csv")

# simple plots to show what this is:
chaokun_dat %>% filter(year>2015) %>% ggplot(aes(year, pfpr210, color = name_1)) + geom_line() +
  theme_bw() + xlab("Year") + ylab("Malaria Prevalence") +
  scale_color_discrete(name = "Admin")

chaokun_dat %>% filter(year>2015) %>% ggplot(aes(year, ft, color = name_1)) + geom_line() +
  theme_bw() + xlab("Year") + ylab("Treatment Coverage (ft)") +
  scale_color_discrete(name = "Admin")

# Next we need to get an ft range and an EIR range for these malaraia prevalence
pfpr_to_eir_heuristic <- function(ft, pfpr, sdd = 0.03){

  if(sdd < 0.03) stop("sdd must not be less than 0.03 or likelihood may fail")

  ll_function <- function(params, data) {

    # get eq
    start <- phi_eir_rel(exp(params[["eir"]]), data$ft)

    # what is the prev
    prev <- start$D + start$T + start$A

    # calculate the log likelihood using a truncated normal
    ll <- log(truncnorm::dtruncnorm(prev, a = 0, b = 1, mean = data$pfpr, sd = data$sdd))
    return(ll)
  }

  # starting conditions
  # we use a log to make it easier for the solver
  start <- c("eir" = log(20))

  # create the data list
  data <- list("ft" = ft, "pfpr" = pfpr, "sdd" = sdd)

  # create bounds
  lower <- log(c(0.01))
  upper <- log(c(500))

  # fit our best model
  fit <- optim(
    par = start,
    fn = ll_function,
    data = data,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(
      trace = TRUE,
      fnscale = -1,
      maxit = 10000
    ),
    hessian = TRUE
  )

  return(exp(as.numeric(fit$par)))

}

# what are the ranges in pfpr and ft for each region
uga_ranges <- chaokun_dat %>%
  filter(year > 2015) %>%
  group_by(name_1) %>%
  summarise(pfpr_low = min(pfpr210),
            pfpr_high = max(pfpr210),
            ft_low = min(ft),
            ft_high = max(ft))

# now work out what the eir should be
uga_eir_ranges <- uga_ranges %>%
  split(.$name_1) %>%
  map(.f = function(x){
    x$eir_high <- pfpr_to_eir_heuristic(x$ft_low, x$pfpr_high)
    x$eir_low <- pfpr_to_eir_heuristic(x$ft_high, x$pfpr_low)
    return(x)
  }) %>%
  do.call(rbind, .) %>%
  rename(district = name_1)

# -------------------------------------------------------- o
# 2. Create data for range function --------
# -------------------------------------------------------- o

# first we need to get the actual freq data
uga_data <- read.csv("https://raw.githubusercontent.com/bailey-lab/selmar/main/analysis/data/data-derived/Uganda_allele_frequency.txt", sep = " ")
uga_data <- uga_data %>% filter(Locus == "K13") %>%
  select(district = District, year, n, freq)

# work out suitable starting res
start_freq <- function(x){
  if(x[1]>0) {
    return(x[1])
  } else {
    return(mean(x[seq_len(which(x>0)[1])]))
  }
}

# find out starting freq and join with eir, ft data
uga_start_res <- uga_data %>% group_by(district) %>%
  summarise(f1 = start_freq(freq),
            years = max(year)-min(year),
            year0 = year[1])

# not to just fix a few spelling and admin name differences
uga_eir_ranges <- uga_eir_ranges %>%
  mutate(district = replace(district, district == "Kabale", "Rukiga"))
uga_eir_ranges <- uga_eir_ranges %>%
  mutate(district = replace(district, district == "Amolatar", "Amoleta"))

# and merge this data together
uga_param_data <- left_join(uga_eir_ranges %>% mutate(district = replace(district, district == "Amolatar", "Amoleta")), uga_start_res)

# -------------------------------------------------------- o
# 3. Run range model --------
# -------------------------------------------------------- o

# let's do this in parallel

# set up cluster
future::plan(future::multisession, workers = 8)

# a function to run our param ranges function and grabbing just the prevalence and summarising
generate_sims <- function(uga_param_data) {

  sim_res <- uga_param_data %>% split(.$district) %>%
    furrr::future_map(function(x){

      param_ranges <- list(
        EIR = seq(x$eir_low, x$eir_high, length.out = 10),
        ft = seq(x$ft_low, x$ft_high, length.out = 10),
        rTR_true = c(0.01, 0.02, 0.05, 0.1, 0.2)
      )

      out <- range_model(param_ranges, use_eir_ft = TRUE,
                         run_time = (x$years+2)*365,
                         time_scale = "day",
                         ton = 2*365,
                         toff = (x$years+2)*365,
                         day0_res = x$f1,
                         init_res = 0,
                         res_time = 0,
                         rTR_true = 0.1,
                         verbose = FALSE)

      clean <- out %>% filter(t >= 2*365) %>%
        mutate(t = t-(2*365)) %>%
        mutate(year = x$year0 + t/365) %>%
        select(year, EIR, ft, rTR_true, prevalence_res) %>%
        group_by(year, rTR_true) %>%
        summarise(res_low = quantile(prevalence_res, 0.025),
                  res_med = quantile(prevalence_res, 0.5),
                  res_high = quantile(prevalence_res, 0.975))


      return(clean %>% mutate(district = x$district))

    }, .options = furrr_options(seed = TRUE))

  return(sim_res)
}

# generate our results now
sim_res <- generate_sims(uga_param_data)

# quick plot to look at these
comparison_plot <- function(sim_res, uga_data) {
  do.call(rbind, sim_res) %>%
    ggplot(aes(year, res_med, ymin = res_low, ymax = res_high, color = as.factor(1/rTR_true), fill = as.factor(1/rTR_true))) +
    geom_line() +
    geom_ribbon(alpha=0.2) +
    geom_point(data = uga_data, aes(year, freq, size = n), inherit.aes = FALSE) +
    # This geom smooth is fitting a binomial regression so it goes between 0 and 1
    geom_smooth(
      method="glm",
      data = uga_data,
      aes(pos = round(n*freq), total = n, x = year, y = freq, group = district),
      method.args=list(family="binomial"),
      formula = cbind(pos, total-pos) ~ x,
      fullrange = FALSE,
      span = 10,
      se = FALSE,
      show.legend = FALSE,
      color = "black",
      inherit.aes = FALSE
    ) +
    facet_wrap(~district, scales = "free_y") +
    theme_bw() +
    xlab("Year") +
    ylab("Resistance Prevalence") +
    scale_color_viridis_d(name = "Parasite Clearance \nDuration (days)") +
    scale_fill_viridis_d(name = "Parasite Clearance \nDuration (days)") +
    scale_size_binned(name = "Sample Size\n", breaks = seq(0,200,40), limits = c(0, 200), range = c(0,4))
}

# now make the plot
gg1 <- comparison_plot(sim_res, uga_data)
gg1

# these plots look good but interestingly the starting value of the binomial model smooth
# starts quite differntly. Perhaps we should just start the model with the assumed starting
# value from the binomial model

# the warning is just for Jinja which is very hard to fit too
new_uga_start_res <- uga_data %>%
  split(.$district) %>%
  map(function(x){
    model <- glm(cbind(pos, n-pos) ~ year, data = x %>% mutate(pos = round(freq*n)), family = "binomial")
    f1 <- exp(predict(model))[1]
    data.frame("district" = x$district[1], f1 = f1)
  }) %>%
  do.call(rbind, .)

# and merge this data together
uga_new_param_data <- left_join(uga_eir_ranges %>% mutate(district = replace(district, district == "Amolatar", "Amoleta")),
                            uga_start_res %>% select(-f1) %>% left_join(new_uga_start_res))

# generate our results now
sim_res2 <- generate_sims(uga_new_param_data)

# now make the plot
gg2 <- comparison_plot(sim_res2, uga_data)
gg2

# great that looks better so let's save the plot
# N.B. You may want to add this to your R package if you use it a lot
save_figs <- function(name,
                      fig,
                      width = 6,
                      height = 6,
                      plot_dir = file.path(here::here(), "analysis/plots"),
                      pdf_plot = TRUE,
                      font_family = "Helvetica") {

  if(!is.null(font_family)) {
    fig <- fig + ggplot2::theme(text = ggplot2::element_text(family = font_family))
  }

  dir.create(plot_dir, showWarnings = FALSE)
  fig_path <- function(name) {paste0(plot_dir, "/", name)}

  cowplot::save_plot(filename = fig_path(paste0(name,".png")),
                     plot = fig,
                     base_height = height,
                     base_width = width)

  if(pdf_plot) {
    pdf(file = fig_path(paste0(name,".pdf")), width = width, height = height)
    print(fig)
    dev.off()
  }

}

# and save this figure
save_figs("uganda_comparison", gg2, width = 10, height = 6)
