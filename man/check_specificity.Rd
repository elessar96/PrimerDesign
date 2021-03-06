% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primer_design_functions.R
\name{check_specificity}
\alias{check_specificity}
\title{Primer design}
\usage{
check_specificity(primers, reference, cutoff = 0.005)
}
\arguments{
\item{primers}{A data.frame as produced by design_primers()}

\item{reference}{Input sequences (DNAStringSet object) that primers should be compared to.}

\item{cutoff}{Threshold of the Hamming distance below which sequences are treated as part of the same cluster.}
}
\value{
data.frame as from design_primers() with additional information on matching sequences.
}
\description{
This function searches for matches with other sequences than the target sequence.
}
\details{
This function searches for sequences in the data set that match with the primer-probe set and stores the sequence name and minimum alignment score of any mismatching sequence in the primers data.frame. To this end it first splits the input data into multiple clusters and compares the primer sequences to each of the cluster consensus sequences.
}
\keyword{primer}
