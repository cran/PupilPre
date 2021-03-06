% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{NA_summary}
\alias{NA_summary}
\title{Check missing data}
\usage{
NA_summary(data, Summary = "Event", PupilColumn = NULL,
  ReturnData = FALSE)
}
\arguments{
\item{data}{A data table object output by \code{\link{create_time_series}}.}

\item{Summary}{A character string indicating the type of summary.}

\item{PupilColumn}{A character string indicating which column to use for the summary.}

\item{ReturnData}{A logical indicating whether to return a data table
containing the summary information.}
}
\value{
Summary information.
}
\description{
\code{NA_summary} summarizes missing data by Event, Subject, or Item.
}
\examples{
# Load example data
data("Pupilex3")

NA_summary(Pupilex3, Summary = "Event", PupilColumn = "Pupil")

# Please see the vignettes for detailed example usage.
# vignette("PupilPre_Basic_Preprocessing", package="PupilPre")

}
