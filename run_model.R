### SET VARS ################################
setwd("/Users/armourc/Documents/Sharpton_Lab/MetaAnalysis/new_database/data/linear_modeling/")
#############################################
library(cplm)
#library(cplm,,lib.loc="/usr/lib64/R/library/")


args <- commandArgs(TRUE)

file      <- args[1]
dbtype    <- args[2]
db        <- args[3]
level     <- args[4]
name      <- args[5]
phenotype <- args[6]

df <- read.table(file,header=T,stringsAsFactors = F,sep="\t")

outpath <- paste("./model_output/",db,sep="")

if( !dir.exists(outpath) ){
  dir.create(outpath,recursive = T)
}

#################
### FUNCTIONS ###
#################

run_model <- function(df){
#   abund   <- df$abund
#   status  <- factor(df$status)
#   study   <- factor(df$study)
#   subject <- df$subject_id
#   sample  <- df$sample_id
#   obesity <- factor(df$obesity_status)
  
  if( phenotype == "obesity" ){
    model   <- cpglm(abund ~ status + study,data=df)
    model   <- cpglm(abund ~ status,data=df)
    info    <- summary(model)
    pvalue  <- info$coefficients["status","Pr(>|t|)"]
  }else if( phenotype == "arthritis" ){
    #model   <- cpglm(abund ~ factor(status,levels = c("none","moderate","high")) + obesity_status,data=df)
    model   <- cpglm(abund ~ factor(status,levels = c("none","moderate","high")),data=df)
    info    <- summary(model)
    p.mod   <- info$coefficients['factor(status, levels = c("none", "moderate", "high"))moderate',"Pr(>|t|)"]
    p.high  <- info$coefficients['factor(status, levels = c("none", "moderate", "high"))high',"Pr(>|t|)"]
    pvalue  <- c(p.mod,p.high)
  }else if( phenotype == "carcinoma" ){
    #model   <- cpglm(abund ~ factor(status,levels = c("controls","advanced adenoma","carcinoma")) + obesity_status,data=df)
    model   <- cpglm(abund ~ factor(status,levels = c("controls","advanced adenoma","carcinoma")),data=df)
    info    <- summary(model)
    p.ad    <- info$coefficients['factor(status, levels = c("controls", "advanced adenoma", "carcinoma"))advanced adenoma',"Pr(>|t|)"]
    p.car   <- info$coefficients['factor(status, levels = c("controls", "advanced adenoma", "carcinoma"))carcinoma',"Pr(>|t|)"]
    pvalue  <- c(p.ad,p.car)
  }else if( length(unique(df$study)) == 1){   
    #model   <- cpglm(abund ~ status + obesity_status,data=df)
    model   <- cpglm(abund ~ status,data=df)
    info    <- summary(model)
    pvalue  <- info$coefficients["status","Pr(>|t|)"]
  }else{
    #model   <- cpglm(abund ~ status + study + obesity_status,data=df)
    model   <- cpglm(abund ~ status + study,data=df)
    info    <- summary(model)
    pvalue  <- info$coefficients["status","Pr(>|t|)"]
  }
  
  return(pvalue)
}

adjust_name <- function(name){
  if( length(strsplit(name,"\\s")[[1]]) > 1 ){
    name <- gsub("\\s","--",name)
  }
  if( length(strsplit(name,"\\/")[[1]]) > 1 ){
    name <- gsub("\\/","..",name)
  }
  return(paste("[[",name,"]]",sep=""))
}

############
### MAIN ###
############
output <- rep(NA,12)
pvalue <- NA
plevel <- NA

if(data.class(try(run_model(df),silent=T)) == "try-error"){
  pvalue <- NA
}else{
  pvalue <- run_model(df)
}
studies   <- paste(unique(df$study),sep="",collapse = ";")

if(phenotype == "arthritis"){
  plevels   <- c("moderate","high")
  case_n    <- c(length(which(df$status == "moderate")),length(which(df$status == "high"))) 
  ctrl_n    <- length(which(df$status == "none"))
  case_mean <- c(mean(df$abund[df$status == "moderate"], na.rm=T),mean(df$abund[df$status == "high"], na.rm=T))
  ctrl_mean <- mean(df$abund[df$status == "none"], na.rm=T)

  output1   <- paste(dbtype,db,level,name,phenotype,plevels[1],studies,case_n[1],ctrl_n,case_mean[1],ctrl_mean,pvalue[1],sep = "\t")
  output2   <- paste(dbtype,db,level,name,phenotype,plevels[2],studies,case_n[2],ctrl_n,case_mean[2],ctrl_mean,pvalue[2],sep = "\t")
  
  output <- paste(output1,output2,sep="\n")
}else if(phenotype == "carcinoma"){
  plevels   <- c("advanced_adenoma","carcinoma")
  case_n    <- c(length(which(df$status == "advanced adenoma")),length(which(df$status == "carcinoma")))
  ctrl_n    <- length(which(df$status == "controls"))
  case_mean <- c(mean(df$abund[df$status == "advanced adenoma"], na.rm=T),mean(df$abund[df$status == "carcinoma"], na.rm=T))
  ctrl_mean <- mean(df$abund[df$status == "controls"], na.rm=T)
  
  output1   <- paste(dbtype,db,level,name,phenotype,plevels[1],studies,case_n[1],ctrl_n,case_mean[1],ctrl_mean,pvalue[1],sep = "\t")
  output2   <- paste(dbtype,db,level,name,phenotype,plevels[2],studies,case_n[2],ctrl_n,case_mean[2],ctrl_mean,pvalue[2],sep = "\t")
  
  output <- paste(output1,output2,sep="\n")
}else{
  case_n    <- length(which(df$status == 1))
  ctrl_n    <- length(which(df$status == 0))
  case_mean <- mean(df$abund[df$status == 1], na.rm=T)
  ctrl_mean <- mean(df$abund[df$status == 0], na.rm=T)
  output <- c(dbtype,db,level,name,phenotype,plevel,studies,case_n,ctrl_n,case_mean,ctrl_mean,pvalue)
}

writeLines(output,sep="\t")

mod_name <- adjust_name(name)

outfile_name <- paste(dbtype,db,level,mod_name,phenotype,sep="_")
outfile <- paste(outpath,"/",outfile_name,"_result.txt",sep="")

write(output,outfile,sep="\t",ncolumns = length(output))

