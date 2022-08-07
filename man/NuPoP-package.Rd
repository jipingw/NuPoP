\name{NuPoP-package}
\alias{NuPoP-package}
\alias{NuPoP}
\docType{package}
\title{
 An R package for nucleosome positioning prediction
}
\description{
\code{NuPoP} is an R package for \code{Nu}cleosome \code{Po}sitioning \code{P}rediction.
This package is built upon a duration hidden Markov model proposed in Xi et al 2010 and
Wang et al 2008. The core of the package was written in Fortran. \code{NuPoP} has prediction
function trained based on both MNase maps and chemical maps of nucleosome positioning data,
all generated from Ji-Ping Wang's lab. Some detailed information about this package can be
found at Ji-Ping Wang's github: https://github.com/jipingw/NuPoP, and its performatics relative to another derived package nuCpos can be found at https://github.com/jipingw/NuPoP_doc.
\code{NuPoP} also provides additional functions for visualization of nucleosome predictions.
Four functions  including \code{predNuPoP},\code{predNuPoP_chem}, \code{readNuPoP},
and \code{plotNuPoP} are provided for nucleosome positioning prediction based on MNase map and chemical map respecively, prediction results readin, and prediction results  visualization respectively.
The input DNA sequence can be of any length.
}
\details{
\tabular{ll}{
Package: \tab NuPoP\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-06-24\cr
License: \tab  GPL-2\cr
}
\code{predNuPoP}: R function invoking Fortran codes to predict nucleosome positioning, nucleosome occupancy and binding affinity.

\code{readNuPoP}: R function to read in the prediction results by \code{predNuPoP}.

\code{plotNuPoP}: R function to visualize predictions.

\code{predNuPoP_chem}: R function invoking Fortran codes to predict nucleosome positioning, nucleosome occupancy and binding affinity using profiles trained based on chemically mapped nucleosomes in yeast, S.pombe, mouse and human.
}
\author{
Ji-Ping Wang, Liqun Xi

Maintainer: Ji-Ping Wang<jzwang@northwestern.edu>
}
\references{
Xi, L., Fondufe-Mittendorf, Y., Xia, L., Flatow, J., Widom, J. and Wang, J.-P. (2010), Predicting nucleosome positioning using a duration Hidden Markov Model, BMC Bioinformatics , doi:10.1186/1471-2105-11-346

Wang, J.-P., Fondufe-Mittendorf, Y., Xi, L., Tsai, G., Segal, E. and Widom, J.(2008), Preferentially quantized linker DNA lengths in Saccharomyces cerevisiae, PLoS Computational Biology, 4(9) e1000175
}

\examples{
\dontrun{
predNuPoP(system.file("extdata", "test.seq", package="NuPoP"),species=7,model=4)

## the prediction results are stored in the current working directory
## the user should replace "system.file("extdata","test.seq_Prediction4.txt",package="NuPoP")"
## by the actual path and file name generated from prediction.

temp=readNuPoP(system.file("extdata","test.seq_Prediction4.txt",package="NuPoP"),startPos=1,endPos=5000)
plotNuPoP(temp)
}
}
