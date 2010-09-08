readNuPoP <-function(file, startPos, endPos)
  ## ====================================================================================================
  ## Purpose: This function reads in the optimal nucleosome positioning map and nucleosome occupancy map
  ## input:   file    --- a string specifying the output file name of the predictions results from
  ##                      function predNpred
  ##          startPos--- an integer specifying the start position of the DNA sequence for visualizing
  ##                      of the nucleosome positioning map
  ##          endPos  --- an integer specifying the end position of the DNA sequence for visualizing
  ##                      of the nucleosome positioning map
  ## ==================================================================================================

{
    skp = startPos
    nop = endPos - startPos + 1
    file = as.character(file)
    results=read.table(file, skip = skp, nrows = nop)
    colnames(results)=c("Position","P.start", "Occup", "N/L", "Affinity")
    return(results)
}
