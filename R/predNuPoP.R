
predNuPoP=function(file,species=7,model=4)
{

  ## ====================================================================================================
  ## Purpose: This function is a wrapper of Fortran codes for prediction of nucleosome positioning
  ## input:   file    --- a string specifying the file name of DNA sequence in FASTA format
  ##          species --- an integer between 0 and 10 to specify which organism the DNA sequence is from.
  ##                      The default is set as 7 for Yeast. 
  ##          model   --- 1 or 4, specifying the order of Markov models to be used in the
  ##                      duration Hidden Markov model.
  ## output:  a result file in text format is output to the current working directory.
  ## ==================================================================================================

  file=as.character(file); n=nchar(file); maxlen=500; maxlen=as.integer(maxlen)
  species=as.integer(species); model=as.integer(model)
  rep=1; rep=as.integer(rep); ind=0
  file_name_num=as.integer(charToRaw(file))

  
  if(model==1){
    results=.Fortran("vtbfb",n,as.integer(file_name_num),freqL,tranL,freqN,tranN,maxlen,rep,species,Pd,ind=as.integer(ind),PACKAGE = "NuPoP")
    ind=results$ind
    if(ind==0){
      FilePath=getwd()
      parts=strsplit(file,"/")[[1]]
      cat(paste("Prediction output: '"), FilePath, "/", parts[length(parts)],"_Prediction1.txt'",sep="")
    }
  } else if(model==4){
    results=.Fortran("vtbfbNL4",n,as.integer(file_name_num),freqL,tranL,tranL2,tranL3,tranL4,freqN4,tranN4,maxlen,rep,species,Pd,ind=as.integer(ind),PACKAGE = "NuPoP")
    ind=results$ind
    if(ind==0){
      FilePath=getwd()
      parts=strsplit(file,"/")[[1]]
      cat(paste("Prediction output: '"), FilePath, "/", parts[length(parts)],"_Prediction4.txt'",sep="")
    }
  } else if(model!=1 && model!=4){
    stop("'model' should be 1 or 4 only; exit!")
  }

  if(ind==1){
    print("In current directory, the input file does not exist, stop!")
  } else if(ind==2){
    print("The input file is not in FASTA format, or contains characters other than A/a, C/c, G/g, T/t, N/n, stop!")
  } else if(ind==3){
    print("The input species label is incorrectly specificed. It should be an integer from 0 to 10!")
  }
}
