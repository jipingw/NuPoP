\name{plotNuPoP}
\alias{plotNuPoP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{R function for plotting the predicted nucleosome positioning map and nucleosome occupancy map}
\description{This function produces two plots from a specified region based on the prediction results from function \code{predNuPoP}. The first plot is the nucleosome occupancy (grey color). In the second plot, in addition to the occupancy,  Viterbi prediction (red rectangle) and the posterior probability for a position to be the start of a nucleosome (blue color) are superimposed.
}
\usage{plotNuPoP(predNuPoPResults)}

\arguments{
  \item{predNuPoPResults}{NuPoP prediction results from predNuPoP function. It must be a data frame read in by readNuPoP function.}
}
\value{\code{plotNuPoP} outputs two plots: the nucleosome occupancy score map and Viterbi optimal nucleosome positioning map together with posterior probability for a position to be the start of a nucleosome.
}


\examples{
library(NuPoP)

## the prediction results are stored in the current working directory
## the user should replace "system.file("extdata","test.seq_Prediction4.txt",package="NuPoP")"
## by the actual path and file name generated from prediction.
temp=readNuPoP(system.file("extdata","test.seq_Prediction4.txt",package="NuPoP"),startPos=1,endPos=5000)
plotNuPoP(temp)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}

