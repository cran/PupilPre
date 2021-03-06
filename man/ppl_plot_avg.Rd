% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{ppl_plot_avg}
\alias{ppl_plot_avg}
\title{Plots average Pupil.}
\usage{
ppl_plot_avg(data, xlim = NA, Column = NULL, Averaging = "Event",
  Condition1 = NULL, Condition2 = NULL, Cond1Labels = NA,
  Cond2Labels = NA, ErrorBar = TRUE, PupilPreTheme = TRUE,
  ConfLev = 95, CItype = "simultaneous", ErrorBand = FALSE,
  ErrorType = "SE")
}
\arguments{
\item{data}{A data table object output after having executed
\code{\link{create_time_series}}.}

\item{xlim}{A vector of two integers specifying the limits of the x-axis.}

\item{Column}{A character string specifying the desired column.}

\item{Averaging}{A character string indicating how the averaging should
be done. "Event" (default) will produce the overall mean in the data, while
"Subject" or "Item" (or, in principle, any other column name) will
calculate the grand mean by that factor.}

\item{Condition1}{A string containing the column name corresponding to the
first condition, if available.}

\item{Condition2}{A string containing the column name corresponding to the
second condition, if available.}

\item{Cond1Labels}{A named character vector specifying the desired custom
labels of the levels of the first condition.}

\item{Cond2Labels}{A named character vector specifying the desired custom
labels of the levels of the second condition.}

\item{ErrorBar}{A logical indicating whether error bars should be
included in the plot.}

\item{PupilPreTheme}{A logical indicating whether the theme included with the
function should be applied, or ggplot2's base theme (to which any other
custom theme could be added)..}

\item{ConfLev}{A number indicating the confidence level of the CI.}

\item{CItype}{A string indicating "simultaneous" or "pointwise". Simultaneous
performs a Bonferroni correction for the interval.}

\item{ErrorBand}{A logical indicating whether error bands should be
included in the plot.}

\item{ErrorType}{A string indicating "SE" (Standard Error) or "CI"
(Confidence Interval).}
}
\description{
\code{ppl_plot_avg} calculates the grand or conditional averages with
standard error. It then plots the results.
N.B.: This function will work for data with a maximum of 2 conditions.
}
\examples{
# Load example data
data("Pupilex7")

ppl_plot_avg(data = Pupilex7, xlim = c(0, 1900), Column = "Pupil",
             Condition1 = NULL, Condition2 = NULL, Cond1Labels = NA,
             Cond2Labels = NA, ErrorBar = TRUE, PupilPreTheme = TRUE)

# Please see the vignettes for detailed example usage.
# vignette("PupilPre_Plotting", package="PupilPre")

}
