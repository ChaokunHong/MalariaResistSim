# AMRSpreadModel

**AMRSpreadModel** is an R package designed to create and run deterministic models of the spread of antimalarial drug resistance. This package provides a range of functionalities to simulate malaria transmission dynamics and visualize the results.

## Installation

To install the package, you can use the `devtools` package to install directly from GitHub:

```{r}
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("ChaokunHong/MalariaResistSim")
```

## Dependencies

The package imports several other R packages, including: **odin**, **ICDMM**, **tidyr**, **stats**, **ggplot2**, **data.table**, **utils**, **dplyr**, **Rcpp**, **magrittr**. These packages will be installed automatically when you install the **AMRSpreadModel** package.

## Usage

### Creating a Malaria Model

The **malaria_model** function is used to create a model of malaria transmission dynamics. You can either provide a list of parameters directly or generate parameters based on the Entomological Inoculation Rate (EIR) and treatment rate (ft).

#### Example 1: Create a model using parameters

```{r}
library(AMRSpreadModel)

# Define parameters
params <- list(
  S = 1000, D = 10, A = 5, T = 3, Sv = 200, Ev = 10, Iv = 5,
  m = 0.1, a = 0.3, b = 0.55, phi = 0.3, fT = 0.5, rD = 0.1, rA = 0.1,
  rT_S = 0.1, rT_R = 0.1, e = 0.2, mu = 0.1, n = 10, cA = 0.1, cD = 0.2, cT = 0.3,
  ton = 5000, toff = 50000, res_time = 3000, init_res = 0.01, rTR_true = 0.1
)

# Create the model
model <- malaria_model(params = params)

# Print model summary
print(model)
```

#### Example 2: Create a model using EIR and treatment rate

```{r}
# Create the model using EIR and ft
model <- malaria_model(EIR = 10, ft = 0.5)

# Print model summary
print(model)
```

### Running the Model

Once the model is created, you can run it over a specified period and visualize the results. \#### Example: Run the model for 1000 days

```{r}
# Run the model for 1000 days
results <- model$run(0:1000)

# Plot the results
plot <- plot_model(results, res_time = 3000, ton = 5000, toff = 50000)
print(plot)
```

### Running the Model Over Ranges of Parameters

The **range_model** function allows you to run the model over a range of parameter values and summarize the results. \#### Example: Run the model over parameter ranges

```{r}
# Define parameter ranges
param_ranges <- list(EIR = seq(1, 20, by = 1), ft = seq(0.1, 0.9, by = 0.1))

# Run the model over the parameter ranges
range_results <- range_model(param_ranges, run_time = 1000, time_points = c(100, 200, 300, 400, 500))

# Plot the results for the parameter ranges
plot <- plot_range(range_results, x_var = "EIR", y_var = "prevalence", color_var = "ft")
print(plot)
```

### Using **modelrun** for Custom Simulations

The **modelrun** function allows you to run the model with custom parameters and visualize the results. \#### Example: Run the model with custom parameters

```{r}
# Define parameters
params <- list(
  S = 1000, D = 10, A = 5, T = 3, Sv = 200, Ev = 10, Iv = 5,
  m = 0.1, a = 0.3, b = 0.55, phi = 0.3, fT = 0.5, rD = 0.1, rA = 0.1,
  rT_S = 0.1, rT_R = 0.1, e = 0.2, mu = 0.1, n = 10, cA = 0.1, cD = 0.2, cT = 0.3,
  ton = 5000, toff = 50000, res_time = 3000, init_res = 0.01, rTR_true = 0.1
)

# Create the model
model <- malaria_model(params = params)

# Define custom time points
time_points <- c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)

# Run the model with custom time points
results <- model$run(time_points)

# Print the results
print(results)
```

### Creating Equilibrium Initial Conditions

The **equilibrium_init_create** function generates initial conditions for the model based on equilibrium assumptions. \#### Example: Create equilibrium initial conditions

```{r}
# Define parameters
params <- list(EIR = 10, ft = 0.5, phi = 0.3, b = 0.55, m = 0.1, a = 0.3, cA = 0.1, cD = 0.2, cT = 0.3, n = 10, mu = 0.1, rD = 0.1, rA = 0.1, rT_S = 0.1, rT_R = 0.1)

# Create equilibrium initial conditions
equilibrium_init <- equilibrium_init_create(params)

# Print equilibrium initial conditions
print(equilibrium_init)
```

### Generating Parameters Based on EIR and ft

The **phi_eir_rel** function generates a list of model parameters based on specified values of EIR and ft. \#### Example: Generate parameters based on EIR and ft

```{r}
# Generate parameters
params <- phi_eir_rel(EIR = 10, ft = 0.5)

# Print generated parameters
print(params)
```

## Contributing

Contributions to the package are welcome. You can contribute by reporting issues, suggesting new features, or submitting pull requests.

## License

This package is licensed under the GPL-3 license.

## Contact

For any questions or issues, please contact Chaokun Hong at [chaokunhong\@gmail.com](mailto:chaokunhong@gmail.com)
