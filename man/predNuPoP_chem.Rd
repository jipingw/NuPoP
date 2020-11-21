\name{predNuPoP_chem}
\alias{predNuPoP_chem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{R function for nucleosome positioning prediction, occupancy score and nucleosome binding affinity score calculation using chemical map profile}
\description{This function invokes Fortran codes to compute the Viterbi prediction of nucleosome positioning, nucleosome occupancy score and nucleosome binding affinity score . A pre-trained linker DNA length distribution for the current species is used in a duration Hidden Markov model. Nucleosome profile is trained based on chemical maps.}

\usage{predNuPoP_chem(file,species=7,model=4)}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{file}{a string for the path and name of a DNA sequence file in FASTA format.  This sequence file can be located in any directory. It must contain only one sequence of any length. By FASTA format, we require each line to be of the same length (the last line can be shorter; the first line should be '>sequenceName'). The length of each line should be not longer than 10 million bp.}
\item{species}{an integer from 0 to 11 as the label for a species indexed as follows: 1 = Human; 2 = Mouse; 3 = Rat; 4 = Zebrafish; 5 = D. melanogaster; 6 = C. elegans; 7 = S. cerevisiae; 8 = C. albicans; 9 = S. pombe; 10 = A. thaliana;  11 = Maize; 0 = Other. The default is 7 = S. cerevisiae . If \code{species=0} is specified, NuPoP will identify a species from 1-11 that has most similar base composition to the input sequence, and then use the models from the selected species for prediction.}
\item{model}{an integer = 4 or 1. NuPoP has two models integrated. One is the first order Markov chain for both nucleosome and linker DNA states. The other is 4th order (default). The latter distinguishes nucleosome/linker in up to 5-mer usage, and thus is slightly more effective in prediction, but runs slower. The time used by 4th order model is about 2.5 times of the 1st order model.}
}
\value{\code{predNuPoP} outputs the prediction results into the current working directory. The output file is named after the input file with an added extension \code{_Prediction1.txt} or \code{_Prediction4.txt}, where 1 or 4 stands for the order of Markov chain models specified. The output file has five columns, \code{Position}, \code{P-start}, \code{Occup}, \code{N/L}, \code{Affinity}:
  \item{Position}{position in the input DNA sequence}
  \item{P-start}{probability that the current position is the start of a nucleosome}
\item{Occup}{nucleosome occupancy score}
\item{N/L}{nucleosome (1) or linker (0) for each position based on Viterbi prediction}
\item{Affinity}{nucleosome binding affinity score}
}

\examples{
library(NuPoP)
predNuPoP_chem(system.file("extdata", "test.seq", package="NuPoP"),species=7,model=4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}

