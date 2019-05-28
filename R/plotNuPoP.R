plotNuPoP=function(predNuPoPResults)
{
  ## ====================================================================================================
  ## Purpose: This function plots the optimal nucleosome positioning map and nucleosome occupancy map
  ## input:   predNuPoPResults--- NuPoP prediction results from predNuPoP function. It must be
  ##          a data frame read in by readNuPoP function.
  ## output:  on-screen graphs
  ## ==================================================================================================

  par(mfrow=c(2,1))

  nop = nrow(predNuPoPResults)
  plot(predNuPoPResults[,1],predNuPoPResults[,3],type='n',ylim=c(-0.05,1.0),xlab='position',ylab='occupancy')
  title("Nucleosome occupancy plot")
  polygon(c(predNuPoPResults[1,1],predNuPoPResults[,1],predNuPoPResults[nop,1]),c(0,predNuPoPResults[,3],0),col=8)

  plot(predNuPoPResults[,1],predNuPoPResults[,3],type='n',ylim=c(-0.05,1.0),xlab='position',ylab='probability/occupancy')
  title("Occupancy(grey)/probability(blue)/Viterbi(red)")
  polygon(c(predNuPoPResults[1,1],predNuPoPResults[,1],predNuPoPResults[nop,1]),c(0,predNuPoPResults[,3],0),col=8)
  points(predNuPoPResults[,1],predNuPoPResults[,4],type='l',col=2)
  points(predNuPoPResults[,1],predNuPoPResults[,2],type='h',col=4)
}
