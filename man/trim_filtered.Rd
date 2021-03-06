% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process.R
\name{trim_filtered}
\alias{trim_filtered}
\title{Trim the beginning and end of filtered events.}
\usage{
trim_filtered(data = data, RmSkipped = NULL, RmEdges = NULL)
}
\arguments{
\item{data}{A data frame object created from \code{auto_cleanup}.}

\item{RmSkipped}{A logical value indicating whether or not to remove events
that were skipped during the filtering process (due to NAs).}

\item{RmEdges}{A numeric vector of length 2 indicating the number of
milliseconds to remove from the beginning and end of each event.}
}
\value{
An object of type data table as described in \link[tibble]{tibble}.
}
\description{
\code{trim_filtered} removes events skipped by the filter as well as a
specified number of milliseconds from the beginning and the end of each
filtered event (as to remove artifacts created by the filter).
}
\examples{
# Load example data
data("Pupilex6")

dat <- trim_filtered(data = Pupilex6, RmSkipped = TRUE,
                     RmEdges = c(75, 75))

# Please see the vignettes for detailed example usage.
vignette("PupilPre_Interpolation_and_Filtering", package="PupilPre")

}
