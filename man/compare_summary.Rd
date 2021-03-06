% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{compare_summary}
\alias{compare_summary}
\title{A utility function to compare pupil size data before and after applying
the cleanup}
\usage{
compare_summary(data, ReturnData = FALSE)
}
\arguments{
\item{data}{A data table object output by \code{\link{create_time_series}}.}

\item{ReturnData}{A logical indicating whether to return a data table
containing the summary information.}
}
\value{
Summary information.
}
\description{
\code{compare_summary} is a utility function to compare pupil size data
before and after applying the cleanup and summarizes a comparison between
Pupil and Pupil_Previous by Event
}
\examples{
# Load example data
data("Pupilex4")

compare_summary(Pupilex4)

# Please see the vignettes for detailed example usage.
# vignette("PupilPre_Basic_Preprocessing", package="PupilPre")

}
