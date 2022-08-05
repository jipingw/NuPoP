predNuPoP_chem=function(file,species=7,model=4)
{
  ## ====================================================================================================
  ## Purpose: This function is a wrapper of Fortran codes for prediction of nucleosome positioning
  ## input:   file    --- a string specifying the file name of DNA sequence in FASTA format
  ##          species --- an integer between 0 and 10 to specify which organism the DNA sequence is from.
  ##                      The default is set as 7 for Yeast. 
  ##          model   --- 1 or 4, specifying the order of Markov models to be used in the
  ##                      duration Hidden Markov model.
  ## output:  a result file in text format is output to the current working directory.
  ##
  ## ==================================================================================================

  file=as.character(file); n=nchar(file); maxlen=500; maxlen=as.integer(maxlen)
  species=as.integer(species); model=as.integer(model)
  rep=1; rep=as.integer(rep); ind=0
  file_name_num=as.integer(charToRaw(file))

  if(model==1){
    if(species==1){ #huamn, species=7 was passed to fortran codes such that re-scaling is not needed
      results=.Fortran("cvtbfb",n,as.integer(file_name_num),cfreqL_h,ctranL_h,cfreqN_h,ctranN_h,maxlen,rep,species=species,Pd,ind=as.integer(ind),PACKAGE = 
"NuPoP")
    }else if(species==2){ #mouse
      results=.Fortran("cvtbfb",n,as.integer(file_name_num),cfreqL_m,ctranL_m,cfreqN_m,ctranN_m,maxlen,rep,species=species,Pd,ind=as.integer(ind),PACKAGE = 
"NuPoP")
    }else if (species==7){ #yeast
      results=.Fortran("cvtbfb",n,as.integer(file_name_num),cfreqL_y,ctranL_y,cfreqN_y,ctranN_y,maxlen,rep,species=species,Pd,ind=as.integer(ind),PACKAGE = 
"NuPoP")
    }else if(species==9){ #pombe
      results=.Fortran("cvtbfb",n,as.integer(file_name_num),cfreqL_p,ctranL_p,cfreqN_p,ctranN_p,maxlen,rep,species=species,Pd,ind=as.integer(ind),PACKAGE = 
"NuPoP")
    }else{
      results=.Fortran("cvtbfb",n,as.integer(file_name_num),cfreqL_y,ctranL_y,cfreqN_y,ctranN_y,maxlen,rep,species=species,Pd,ind=as.integer(ind),PACKAGE = 
"NuPoP")
    }
    ind=results$ind
    if(ind==0){
      FilePath=getwd()
      parts=strsplit(file,"/")[[1]]
      cat(paste("Prediction output: '"), FilePath, "/", parts[length(parts)],"_Prediction1.txt'",sep="")
    }
  }else if(model==4){
    if(species==1){#species=7 was passed to fortran codes such that re-scaling is not needed
      
results=.Fortran("cvtbfbNL4",n,as.integer(file_name_num),cfreqL_h,ctranL_h,ctranL2_h,ctranL3_h,ctranL4_h,cfreqN4_h,ctranN4_h,maxlen,rep,species=species,Pd,ind=as.integer(ind),PACKAGE = "NuPoP")
    }else if(species==2){
      
results=.Fortran("cvtbfbNL4",n,as.integer(file_name_num),cfreqL_m,ctranL_m,ctranL2_m,ctranL3_m,ctranL4_m,cfreqN4_m,ctranN4_m,maxlen,rep,species=species,Pd,ind=as.integer(ind),PACKAGE = "NuPoP")
     }else if (species==7){
      
results=.Fortran("cvtbfbNL4",n,as.integer(file_name_num),cfreqL_y,ctranL_y,ctranL2_y,ctranL3_y,ctranL4_y,cfreqN4_y,ctranN4_y,maxlen,rep,species=species,Pd,ind=as.integer(ind),PACKAGE = "NuPoP")
     }else if(species==9){
    
results=.Fortran("cvtbfbNL4",n,as.integer(file_name_num),cfreqL_p,ctranL_p,ctranL2_p,ctranL3_p,ctranL4_p,cfreqN4_p,ctranN4_p,maxlen,rep,species=species,Pd,ind=as.integer(ind),PACKAGE = "NuPoP")
     }else{
      
results=.Fortran("cvtbfbNL4",n,as.integer(file_name_num),cfreqL_y,ctranL_y,ctranL2_y,ctranL3_y,ctranL4_y,cfreqN4_y,ctranN4_y,maxlen,rep,species,Pd,ind=as.integer(ind),PACKAGE = "NuPoP")
    }
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
