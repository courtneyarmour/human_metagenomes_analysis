
##### SET THESE VARS ##########################################################
setwd("/Users/armourc/Documents/Sharpton_Lab/MetaAnalysis/new_database/")     #
database <- "metaphlan2"                                                            #
outdir <- paste("./data/linear_modeling/",database,"_model_tables/",sep="")   #
infodir <- paste("./data/linear_modeling/",database,"_summary_tables/",sep="")#
total_dir_num <- 100                                                          #
###############################################################################

if( !dir.exists(outdir) ){
  dir.create(outdir,recursive = T)
}
if( !dir.exists(infodir) ){
  dir.create(infodir,recursive = T)
}

##############
### INPUTS ###
##############
source("./scripts/metagenome_analysis_functions.R")

### metadata
all_attributes <- read.table("./data/all_attributes.txt",header=T,stringsAsFactors = F,sep="\t")
# convert to table with ids as rows and attibutes as columns
att.map <- convert_attibutes(all_attributes)

### samples to exclude 
exclude_samples <- read.table("./data/exclude_samples.txt",header=T,stringsAsFactors = F,sep="\t")

#sample_to_subject mapping table
sample_to_subject <- read.table("./data/sample_to_subject.txt",header=T,stringsAsFactors = F,sep="\t")

### sample mapping
#sample_mapping <- read.table("./data/sample_mapping.txt",header=T,stringsAsFactors = F,sep="\t")

#representative samples table
#rep_samps <- read.table("./data/representative_samples.txt",header=T,stringsAsFactors = F,sep="\t")

### abundances
# database type var
if(database == "mOTU" | database == "metaphlan2"){
  dbtype <- "taxa"
}else if(database == "kegg" | database == "eggnog"){
  dbtype <- "function"
}
#automate input of table
filename <- paste("./data/abundance_",dbtype,"_",database,".txt",sep="")
abunds <- read.table(file = filename,
                     header = T,
                     stringsAsFactors = F,
                     sep = "\t")

# remove info for samples to exclude
abunds <- exclude_samp_abunds(abunds,exclude_samples)

### objects in the databases
all_names     <- read.table("./data/all_names.txt",header=T,stringsAsFactors = F,sep="\t")
# limit to database of interest
all_db_names  <- all_names[which(all_names[,"database"] == database),]
# only test objects prevalent in > 0.5% of samples
db_names_0.5  <- all_db_names[which(all_db_names[,"prevalence"] > 0.5),]

############
### MAIN ###
############

# calculate small value to adjust zeros for log transformation
adjust_value <- get_adjust_value(abunds,dbtype)

phenotypes <- c('carcinoma', 'cirrhosis', 't2d', 'arthritis',
                'crohns','ulcerative_colitis', 'obesity')

max_item_num <- find_dir_size(total_dir_num,length(db_names_0.5$name))

item_counter   <- 1
dir_counter    <- 1
summaries      <- NULL
data_summaries <- NULL

for( i in 1:nrow(db_names_0.5) ){
  if( item_counter > max_item_num ){
    dir_counter  <- dir_counter + 1
    item_counter <- 1
  }
  
  table <- as.data.frame(matrix(nrow = length(phenotypes),ncol=11)); row <- 0
  colnames(table) <- c('type','database','level','name','phenotype','studies','case_n',
                       'ctrl_n','case_mean','ctrl_mean','p_value')
  
  type      <- db_names_0.5[i,'type']
  db        <- db_names_0.5[i,'database']
  level     <- db_names_0.5[i,'level']
  
  if(level == "gene_family"){
    mod_level <- "gene"
  }else{
    mod_level <- level
  }

  name      <- db_names_0.5[i,'name']
  mod_name  <- adjust_name(name)
  
  abun <- get_abund(abunds,type,db,level,name)
    
  for( phenotype in phenotypes ){
    #include_comorb <- "obesity"
    include_comorb <- NULL
    subj_ids       <- get_ids(att.map,phenotype,include_comorb)
    
    #samp_ids       <- sample_to_subject[which(sample_to_subject[,"subject_id"] %in% subj_ids),"sample_id"]
    #samp_ids       <- samp_ids[which(samp_ids %in% row.names(abun))]
    #subj_ids       <- sample_to_subject[which(sample_to_subject[,"sample_id"] %in% samp_ids),"subject_id"] 
    
    sub_att.map            <- att.map[which(att.map[,"subject_id"] %in% subj_ids),] 
    row.names(sub_att.map) <- sub_att.map[,"subject_id"]
    
    if( phenotype == "carcinoma" ){
      samp_ids       <- sample_to_subject[which(sample_to_subject[,"subject_id"] %in% subj_ids),"sample_id"]
      samp_ids       <- samp_ids[which(samp_ids %in% row.names(abun))]
      subj_ids       <- sample_to_subject[which(sample_to_subject[,"sample_id"] %in% samp_ids),"subject_id"] 
      
      status <- factor(sub_att.map[subj_ids,"carcinoma_type"],levels = c("controls","advanced adenoma","carcinoma"))
    }else if( phenotype == "arthritis"){
      sub_att.map  <- sub_att.map[which(sub_att.map[,"arthritis_type"] %in% c("none","moderate","high")),]
      subj_ids     <- subj_ids[which(subj_ids %in% sub_att.map[,"subject_id"])]
      samp_ids     <- sample_to_subject[which(sample_to_subject[,"subject_id"] %in% subj_ids),"sample_id"]
      samp_ids     <- samp_ids[which(samp_ids %in% row.names(abun))]
      subj_ids     <- sample_to_subject[which(sample_to_subject[,"sample_id"] %in% samp_ids),"subject_id"]
      
      status <- factor(sub_att.map[subj_ids,"arthritis_type"],levels = c("none","moderate","high"))
    }else if( phenotype == "t2d" ){
      sub_att.map  <- sub_att.map[which(sub_att.map[,"glucose_tolerance"] %in% c("NGT","T2D",NA)),]
      subj_ids     <- subj_ids[which(subj_ids %in% sub_att.map[,"subject_id"])]
      samp_ids     <- sample_to_subject[which(sample_to_subject[,"subject_id"] %in% subj_ids),"sample_id"]
      samp_ids     <- samp_ids[which(samp_ids %in% row.names(abun))]
      subj_ids     <- sample_to_subject[which(sample_to_subject[,"sample_id"] %in% samp_ids),"subject_id"]
      
      status <- sub_att.map[subj_ids,phenotype]
    }else{
      samp_ids       <- sample_to_subject[which(sample_to_subject[,"subject_id"] %in% subj_ids),"sample_id"]
      samp_ids       <- samp_ids[which(samp_ids %in% row.names(abun))]
      if( phenotype %in% c("crohns","ulcerative_colitis") ){
        samp_ids <- remove_bioreps(samp_ids)
      }
      subj_ids       <- sample_to_subject[which(sample_to_subject[,"sample_id"] %in% samp_ids),"subject_id"]
      
      status <- sub_att.map[subj_ids,phenotype]
    }
    
    abund          <- abun[samp_ids,3]
    abund_log      <- log(abund+adjust_value)
    study          <- abun[samp_ids,"study_id"]
    
    df <- data.frame(subj_ids,samp_ids,abund,abund_log,status,study,stringsAsFactors = F) 
    colnames(df) <- c("subject_id","sample_id","abund","abund_log","status","study")
    
#     if( phenotype != "obesity"){
#       obesity_status <- sub_att.map[subj_ids,"obesity"]
#       df <- cbind(df,obesity_status)
#     }
    
    if( !(phenotype %in% summaries) ){
      pheno_summary <- get_summary_stats(df,sub_att.map,phenotype)
      summaries     <- c(summaries, phenotype)
      
      out <- paste(infodir,phenotype,"_attribute_summary.txt",sep="")
      write.table(pheno_summary,file=out,quote=F,sep="\t")
    }
    if( !(phenotype %in% data_summaries ) ){
      if(phenotype == "obesity"){
        cols <- colnames(df[which(!(colnames(df) %in% c("abund","abund_log")))])
        data_summary <- df[,cols]
      }else{
        cols <- colnames(df[which(!(colnames(df) %in% c("abund","abund_log")))])
        data_summary <- df[,cols]
        
      }
      data_summaries <- c(data_summaries,phenotype)
      
      out <- paste(infodir,phenotype,"_data_summary.txt",sep="")
      write.table(data_summary,file=out,quote = F,sep = "\t")
    }
   
     table_name <- paste(type,db,mod_level,mod_name,phenotype,sep="_")
     
     sub_outdir <- paste(outdir,"dir",dir_counter,sep="")
     if( !dir.exists(sub_outdir) ){
       dir.create(sub_outdir)
     }
     
     outfile <- paste(sub_outdir,"/",table_name,".tab",sep="")
     write.table(df,file=outfile,quote=F,sep="\t")
  }
  item_counter <- item_counter+1
}



####
# ra_summary  <- read.table("./data/linear_modeling/mOTU_summary_tables/arthritis_attribute_summary.txt",header=T,stringsAsFactors = F,sep="\t")
# cro_summary <- read.table("./data/linear_modeling/mOTU_summary_tables/crohns_attribute_summary.txt",header=T,stringsAsFactors = F,sep="\t")
# uc_summary  <- read.table("./data/linear_modeling/mOTU_summary_tables/ulcerative_colitis_attribute_summary.txt",header=T,stringsAsFactors = F,sep="\t")
# cir_summary <- read.table("./data/linear_modeling/mOTU_summary_tables/cirrhosis_attribute_summary.txt",header=T,stringsAsFactors = F,sep="\t")
# ob_summary  <- read.table("./data/linear_modeling/mOTU_summary_tables/obesity_attribute_summary.txt",header=T,stringsAsFactors = F,sep="\t")
# t2d_summary <- read.table("./data/linear_modeling/mOTU_summary_tables/t2d_attribute_summary.txt",header=T,stringsAsFactors = F,sep="\t")
# car_summary <- read.table("./data/linear_modeling/mOTU_summary_tables/carcinoma_attribute_summary.txt",header=T,stringsAsFactors = F,sep="\t")
