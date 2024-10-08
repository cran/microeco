% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{tidy_taxonomy}
\alias{tidy_taxonomy}
\title{Clean up the taxonomic table to make taxonomic assignments consistent.}
\format{
\code{\link{data.frame}} object.
}
\usage{
tidy_taxonomy(
  taxonomy_table,
  column = "all",
  pattern = c(".*unassigned.*", ".*uncultur.*", ".*unknown.*", ".*unidentif.*",
    ".*unclassified.*", ".*No blast hit.*", ".*Incertae.sedis.*"),
  replacement = "",
  ignore.case = TRUE,
  na_fill = ""
)
}
\arguments{
\item{taxonomy_table}{a data.frame with taxonomic information (rows are features; columns are taxonomic levels); 
or a microtable object with \code{tax_table} in it.}

\item{column}{default "all"; "all" or a number; 'all' represents cleaning up all the columns; a number represents cleaning up this specific column.}

\item{pattern}{default c(".*unassigned.*", ".*uncultur.*", ".*unknown.*", ".*unidentif.*", ".*unclassified.*", ".*No blast hit.*", ".*Incertae.sedis.*"); 
the characters (regular expressions) to be removed or replaced; removed when parameter \code{replacement = ""}, 
replaced when parameter replacement has something; Note that the capital and small letters are not distinguished when \code{ignore.case = TRUE}.}

\item{replacement}{default ""; the characters used to replace the character in \code{pattern} parameter.}

\item{ignore.case}{default TRUE; if FALSE, the pattern matching is case sensitive and if TRUE, case is ignored during matching.}

\item{na_fill}{default ""; used to replace \code{NA}.}
}
\value{
data.frame
}
\description{
Clean up the taxonomic table to make taxonomic assignments consistent.
}
\examples{
data("taxonomy_table_16S")
tidy_taxonomy(taxonomy_table_16S)
}
