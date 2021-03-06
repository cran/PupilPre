% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process.R
\name{apply_cleanup_change}
\alias{apply_cleanup_change}
\title{Applies user-selected changes to auto cleanup}
\usage{
apply_cleanup_change(data = data, LogFile = NULL)
}
\arguments{
\item{data}{A data frame object created from \code{downsample}.}

\item{LogFile}{A character string indicating the name (and location) of the
log file.}
}
\value{
An object of type data table as described in \link[tibble]{tibble}.
}
\description{
\code{apply_cleanup_change} applies to each event the user-selected
changes to the automatic cleanup based on information stored in the RDS file
using \code{verify_cleanup_app} which was created using either
\code{blink_cleanup} or \code{artifact_cleanup}.
}
\examples{
if (interactive()) {

# Load example data
data("Pupilex3")

# Ensure the log file exists by running cleanup
# Writing log file to temporary folder for the example
dat <- clean_blink(Pupilex3, BlinkPadding = c(100, 100), Delta = 5,
                   MaxValueRun = 5, NAsAroundRun = c(2,2),
                   LogFile = paste0(tempdir(),"/BlinkCleanupLog.rds"))

# Read log file from temporary folder
verify_cleanup_app(dat, LogFile = paste0(tempdir(),"/BlinkCleanupLog.rds"))

# Make verification via the app interface

# Read log file from the temporary folder
dat <- apply_user_cleanup(dat,
                          LogFile = paste0(tempdir(),"/BlinkCleanupLog.rds"))
}

# Please see the vignettes for detailed example usage.
# vignette("PupilPre_Cleanup", package="PupilPre")


}
