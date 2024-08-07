% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria_model.R
\name{malaria_model}
\alias{malaria_model}
\title{Create Malaria Model}
\usage{
malaria_model(
  params = NULL,
  EIR = NULL,
  ft = NULL,
  ton = 5000,
  toff = 50000,
  day0_res = 0,
  init_res = 0.01,
  res_time = 3000,
  rTR_true = 0.1,
  verbose = FALSE
)
}
\arguments{
\item{params}{A list of parameters to pass to the model. If NULL, parameters will be generated based on EIR and ft.}

\item{EIR}{Entomological Inoculation Rate (used if params is NULL)}

\item{ft}{Treatment rate (used if params is NULL)}

\item{ton}{Time at which treatment is turned on}

\item{toff}{Time at which treatment is turned off}

\item{day0_res}{Resistant at Day 0. Default = 0.01}

\item{init_res}{Initial resistance level at res_time}

\item{res_time}{Time at which resistance is introduced}

\item{rTR_true}{True treatment rate for resistant parasites}

\item{verbose}{Logical. If TRUE, prints detailed logs. Default is FALSE.}
}
\value{
An object of class \code{odin_model}.
}
\description{
Create Malaria Model
}
