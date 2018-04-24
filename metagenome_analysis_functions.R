################
### PACKAGES ###
################
library(reshape2)
library(vegan)
library(ggplot2)
library(cowplot)
library(gplots)
library(scales)

#################
### FUNCTIONS ###
#################

convert_attibutes <- function(attributes){
  # convert table
  att.map <- dcast(attributes,subject_id~attribute,value.var = "value",fill=NA,drop=F)
  
  #alter numeric fields
  numeric_fields <- c("age","bmi")
  for(num_field in numeric_fields){
    att.map[,num_field] <- as.numeric(att.map[,num_field])
  }
  
  #fix incorrect zero entries
  att.map[which(att.map[,"bmi"] == 0),"bmi"] <- NA
  att.map[which(att.map[,"age"] == 0),"age"] <- NA
  
  #add missing classes
  att.map[!is.na(att.map[,'bmi']) & att.map[,'bmi'] >= 30, 'obesity'] <- 1
  att.map[!is.na(att.map[,'bmi']) & att.map[,'bmi'] < 30, 'obesity'] <- 0
  att.map[!is.na(att.map[,'ibd_type']) & att.map[,'ibd_type'] == 'crohns disease', 'crohns'] <- 1
  att.map[!is.na(att.map[,'ibd_type']) & att.map[,'ibd_type'] != 'crohns disease', 'crohns'] <- 0
  att.map[!is.na(att.map[,'ibd_type']) & att.map[,'ibd_type'] == 'ulcerative colitis', 'ulcerative_colitis'] <- 1
  att.map[!is.na(att.map[,'ibd_type']) & att.map[,'ibd_type'] != 'ulcerative colitis', 'ulcerative_colitis'] <- 0
  att.map[!is.na(att.map[,'glucose_tolerance']) & att.map[,'glucose_tolerance'] == 'T2D', 't2d'] <- 1
  att.map[!is.na(att.map[,'glucose_tolerance']) & att.map[,'glucose_tolerance'] == 'NGT', 't2d'] <- 0
  #alter factor fields
  factor_fields <- c("arthritis","arthritis_type","bmi_type","carcinoma","carcinoma_type","cirrhosis",
                     "continent","country","glucose_tolerance","hypertension","ibd","ibd_type",
                     "lifestyle","metformin","t1d","t2d","sex","obesity","crohns","ulcerative_colitis")
  for(field in factor_fields){
    att.map[,field] <- as.factor(att.map[,field])
  }
  return(att.map)
}

modify_mapping <- function(mapping){
  status <- mapping$status
  status_mod <- ifelse(status %in%  c(0,"controls","none"),0,1)
  mapping$status <- status_mod
  
  mapping <- cbind(mapping,pheno.status=paste(mapping$phenotype,mapping$status,sep="."))
  
  return(mapping)
}

reduce_abunds <- function(abunds,mapping){
  abund_samps <- paste0(abunds[,"sample_id"],";",abunds[,"study_id"])
  map_samps    <- paste0(mapping[,"sample_id"],";",mapping[,"study"])
  
  abunds_reduced <- abunds[which(abund_samps %in% map_samps),]
  return(abunds_reduced)
}

reduce_mapping <- function(abunds,mapping){
  abund_samps <- paste0(abunds[,"sample_id"],";",abunds[,"study_id"])
  map_samps    <- paste0(mapping[,"sample_id"],";",mapping[,"study"])
  
  mapping_reduced <- mapping[which(map_samps %in% abund_samps),]
  return(mapping_reduced)
} 

split_def <- function(def){
  return(strsplit(def," \\[")[[1]][1])
}

subset_cplm_results <- function(result_table,level,phenotype){
  subset <- result_table[which(result_table$level == level &
                               result_table$phenotype == phenotype),,drop=F]
  return(subset)
}

check_if_empty <- function(matrix){
  if(length(rownames(matrix)) == 0){
    return(1)
  }else{
    return(0)
  }
}

add_correction <- function(results_table){
  phenotypes  <- unique(results_table[,"phenotype"])
  levels      <- unique(results_table[,"level"])
  
  new_table <- NULL
  for(phenotype in phenotypes){
    for(level in levels){
      subset <- subset_cplm_results(results_table,level,phenotype)
      if( check_if_empty(subset) ){
        next
      }
      pvalues <- subset[,"pvalue"]
      
      fdr <- p.adjust(pvalues,method="fdr")
      
      corrected_matrix <- NULL
      if( is.null(new_table) ){
        new_table <- cbind(subset,fdr)
        colnames(new_table) <- c(colnames(subset),"fdr")
        row.names(new_table) <- row.names(subset)
      }else{
        corrected_matrix <- cbind(subset,fdr)
        colnames(corrected_matrix) <- c(colnames(subset),"fdr")
        #check that the tables can be merged (ie same number of columns and in same order)
        if( length(colnames(corrected_matrix)) == length(colnames(new_table)) ){
          for( i in 1:length(colnames(new_table)) ){
            if( colnames(corrected_matrix)[i] != colnames(new_table)[i] ){
              stop("cannot merge tables, columns in different order")
            }
          }
        }else{
          stop("cannot merge tables, unequal number of columns")
        }
        #now combine all_phenotest_wqval table with new values
        new_table <- rbind(new_table,corrected_matrix)
      }
    }
  }
  return(new_table)
}

count_biomarkers <- function(results_wcorrection,cutoffs){
  rnames  <- NULL
  cnames  <- c("Database","Level","Min","Max",cutoffs)
  df      <- NULL
  
  phenotypes <- unique(results_wcorrection[,"phenotype"])
  levels <- unique(results_wcorrection[,"level"])
  for(phenotype in phenotypes){
    for(level in levels){
      subset <- subset_cplm_results(results_wcorrection,level,phenotype)
      if( check_if_empty(subset) ){
        next
      }
      
      #set up table
      rnames <- c(rnames,phenotype)
      values <- NULL    # will hold (database,level,min,max,and counts for each cutoff)
      
      if( length(unique(subset[,"fdr"])) == 1 ){
        if(is.na(unique(subset[,"fdr"]))){
          range <- c(NA,NA)
        }
      }else{
        range <- range(na.omit(subset[,"fdr"]))
      }
      min <- range[1]
      max <- range[2]
      
      values[1] <- "kegg"
      values[2] <- level
      values[3] <- min
      values[4] <- max
      
      for( i in 1:length(cutoffs) ){
        j <- i + 4     # position in value vector
        if( is.na(min) ){
          values[j] <- NA  
        }else if( min > cutoffs[i] ){
          values[j] <- 0
        }else{
          if( phenotype %in% c("carcinoma","arthritis") ){
            below_cutoff <- subset[which(subset[,"fdr"] < cutoffs[i]),]
            values[j] <- length(unique(below_cutoff[,"name"]))
          }else{
            below_cutoff <- subset[which(subset[,"fdr"] < cutoffs[i]),]
            values[j] <- length(row.names(below_cutoff))
          }
        }
      }
      df <- rbind(df,values)
      colnames(df) <- cnames
      row.names(df) <- rnames
    }
  }
  return(df)
}

find_common_biomarkers <- function(result_table,cutoff,detail_results){
  biomarker_table <- as.data.frame(matrix(nrow=0,ncol=8)); row <- 0
  colnames(biomarker_table) <- c("type","database","level","name","phenotypes","fdr","change_in_case","phenotype_count")
  all_names <- unique(result_table[,"name"])
  for( name in all_names ){
    levels <- unique(result_table[which(result_table[,"name"] == name),"level"])
    for(level in levels){
      type <- unique(result_table[,"type"])
      database <- unique(result_table[,"database"])
      subset <- result_table[which(result_table[,"type"] == type &
                                     result_table[,"database"] == database &
                                     result_table[,"level"] == level &
                                     result_table[,"name"] == name),,drop=F]
      #keep only those that are below the cutoff 
      subset_cutoff <- subset[which(subset[,"fdr"] < cutoff),]
      if( nrow(subset_cutoff) > 0){
        #get phenotypes
        phenos <- subset_cutoff[,"phenotype"]
        #get q or fdr values
        vals <- subset_cutoff[,"fdr"]
        # modify to only include arthritis or carcinoma level with lower q or fdr value
        plevel <- NULL
        if("arthritis" %in% phenos){
          arth_data <- subset_cutoff[which(subset_cutoff$phenotype == "arthritis"),]
          arth_vals <- arth_data[,"fdr"]
          if(length(arth_vals) > 1){
            if(arth_vals[1] <= arth_vals[2]){
              plevel <- c(plevel,arthritis=arth_data[,"phenotype_level"][1])
              exclude <- arth_data[,"phenotype_level"][2]
            }else{
              plevel <- c(plevel,arthritis=arth_data[,"phenotype_level"][2])
              exclude <- arth_data[,"phenotype_level"][1]
            }
            subset_cutoff <- subset_cutoff[-which(subset_cutoff$phenotype_level == exclude),]
            phenos           <- subset_cutoff[,"phenotype"]
            vals             <- subset_cutoff[,"fdr"]
          }else{
            plevel <- c(plevel,arthritis=arth_data[,"phenotype_level"])
          }
        }
        if("carcinoma" %in% phenos){
          car_data <- subset_cutoff[which(subset_cutoff$phenotype == "carcinoma"),]
          car_vals <- car_data[,"fdr"]
          if(length(car_vals) > 1){
            if(car_vals[1] <= car_vals[2]){
              plevel <- c(plevel,carcinoma=car_data[,"phenotype_level"][1])
              exclude <- car_data[,"phenotype_level"][2]
            }else{
              plevel <- c(plevel,carcinoma=car_data[,"phenotype_level"][2])
              exclude <- car_data[,"phenotype_level"][1]
            }
            subset_cutoff <- subset_cutoff[-which(subset_cutoff$phenotype_level == exclude),]
            phenos           <- subset_cutoff[,"phenotype"]
            vals             <- subset_cutoff[,"fdr"]
          }else{
            plevel <- c(plevel,carcinoma=car_data[,"phenotype_level"])
          }
          if(plevel["carcinoma"] == "advanced_adenoma"){
            plevel["carcinoma"] <- "adenoma"
          }
        }
        
        
        
        #get direction of change in cases
        change <- subset_cutoff$case_mean - subset_cutoff$ctrl_mean
        #mod_change <- ifelse(change > 0,"U","D")
        
        # get direction of change based on slope estimate
        # change <- NULL
        # for( i in 1:nrow(subset_cutoff) ){
        #   pheno <- subset_cutoff$phenotype[i] 
        #   comp <- "status"
        #   if(pheno == "arthritis"){
        #     comp <- paste0("status_",plevel["arthritis"])
        #   }
        #   if(pheno == "carcinoma"){
        #     comp <- paste0("status_",plevel["carcinoma"])
        #   }
        #   slope <- detail_results[which(detail_results$database == database &
        #                                   detail_results$level == level &
        #                                   detail_results$name == name &
        #                                   detail_results$phenotype == pheno &
        #                                   detail_results$component == comp),"estimate"]
        #   change <- c(change,slope)
        # }
        
        mod_change <- ifelse(change > 0,"U","D")
        if(level == "gene"){
          mod_level <- "gene_family"
        }else{
          mod_level <- level
        }
        
        #fill in table with information
        row <- row + 1
        biomarker_table[row,'type']            <- type
        biomarker_table[row,'database']        <- database
        biomarker_table[row,'level']           <- level
        biomarker_table[row,'name']            <- name
        biomarker_table[row,'phenotypes']      <- paste(phenos,collapse = ";")
        biomarker_table[row,"fdr"]             <- paste(vals,collapse = ";")
        biomarker_table[row,"change_in_case"]  <- paste(mod_change,collapse = ";")
        biomarker_table[row,'phenotype_count'] <- length(phenos)
      }
    }
  }
  return(biomarker_table)
}

count_shared_biomarkers <- function(shared_markers_table){
  databases <- unique(shared_markers_table[,"database"])
  levels <- unique(shared_markers_table[,"level"])
  phenotypes <- unique(unlist(strsplit(shared_markers_table[,"phenotypes"],";")))
  
  count_table <- as.data.frame(matrix(nrow=0,ncol = (length(phenotypes)+4) ));row <- 1
  colnames(count_table) <- c("type","database","level","phenotype",phenotypes)
  
  for(database in databases){
    for(level in levels){
      subset <- shared_markers_table[which(shared_markers_table[,"database"] == database &
                                             shared_markers_table[,"level"] == level),]
      for(phenotypeRow in phenotypes){
        count_table[row,"type"] <- unique(subset[,"type"]) 
        count_table[row,"database"] <- database
        count_table[row,"level"] <- level
        count_table[row,"phenotype"] <- phenotypeRow
        for(phenotypeCol in phenotypes){
          count <- 0
          for( i in 1:length(row.names(subset)) ){
            marker_phenos <- strsplit(subset[i,"phenotypes"],";")[[1]]
            if(phenotypeRow %in% marker_phenos &
               phenotypeCol %in% marker_phenos){
              count <- count + 1
            }
          }
          count_table[row,phenotypeCol] <- count
        }
        row <- row + 1
      }
    }
  }
  return(count_table)
}

create_abund_map <- function(abund_long,database,level,metric){
  abund_subset <- abund_long[which(abund_long[,"db"] == database &
                                     abund_long[,"level"] == level),,drop=F]
  if( check_if_empty(abund_subset) ){
    stop("no data matching those parameters")
  }
  abund.map <- acast(abund_subset,sample_id ~ name,value.var = metric,fill=0) 
}

get_p_name <- function(phenotype){
  if(phenotype == 'carcinoma'){
    name <- "Colorectal cancer"
  }else if(phenotype == 'cirrhosis'){
    name <- "Liver cirrhosis"
  }else if(phenotype == 't2d'){
    name <- "Type II diabetes"
  }else if(phenotype == 'arthritis'){
    name <- "Rheumatoid arthritis"
  }else if(phenotype == 'crohns'){
    name <- "Crohn's disease" 
  }else if(phenotype == 'ulcerative_colitis'){
    name <- "Ulcerative colitis"
  }else if(phenotype == 'obesity'){
    name <- "Obesity"
  }
  return(name)
}

get_p_abbrev <- function(phenotype){
  if(phenotype == 'carcinoma'){
    abbrev <- "CC"
  }else if(phenotype == 'cirrhosis'){
    abbrev <- "LC"
  }else if(phenotype == 't2d'){
    abbrev <- "T2D"
  }else if(phenotype == 'arthritis'){
    abbrev <- "RA"
  }else if(phenotype == 'crohns'){
    abbrev <- "CD" 
  }else if(phenotype == 'ulcerative_colitis'){
    abbrev <- "UC"
  }else if(phenotype == 'obesity'){
    abbrev <- "OB"
  }
  return(abbrev)
}

unadjust_name <- function(name){
  if( length(strsplit(name,"--")[[1]]) > 1 ){
    name <- gsub("--"," ",name)
  }
  if( length(strsplit(name,"\\.\\.")[[1]]) > 1 ){
    name <- gsub("\\.\\.","/",name)
  }
  if( length(strsplit(name,"+")[[1]]) > 1 ){
    name <- gsub("\\+(.*)\\+","(\\1)",name,perl=T)
  }
  return(name)
}

make_binary <- function(x){
  if(is.na(x)){
    x <- NA
  }else if(x > 0){
    x <- 1
  }else if(x < 0){
    x <- -1
  }else{
    x <- x
  }
  return(x)
}

exclude_zero <- function(abund.map){
  mod.abund.map <- abund.map[,colSums(abund.map) != 0]
  return(mod.abund.map)
}

### calculates a small value used to adjust before log transformation (since there are zeros)
get_adjust_value <- function(abunds,db_type){
  if(db_type == "taxa"){
    min_val <- min(abunds$relative_abundance[which(abunds$relative_abundance > 0 )])
  }else if(db_type == "function"){
    min_val <- min(abunds$copy_number[which(abunds$copy_number > 0 )])
  }
  adj_value <- min_val/10
  return(adj_value)
}

### used in building tables to run the tests
### calculates the number of files to include in each directory
find_dir_size <- function(dir_num,test_num){
  if(test_num > dir_num){
    items_per_dir <- as.integer(test_num/dir_num + 0.5)
  }else{
    items_per_dir <- 1
  }
  return(items_per_dir)
}

### used in building tables to run the tests
### adjusts the names of items to test so they can be used in the filename
adjust_name <- function(name){
  if( length(strsplit(name,"\\s")[[1]]) > 1 ){
    name <- gsub("\\s","--",name)
  }
  if( length(strsplit(name,"\\/")[[1]]) > 1 ){
    name <- gsub("\\/","..",name)
  }
  if( length(strsplit(name,"[()]")[[1]]) > 1 ){
    name <- gsub("[()]","+",name)
  }
  return(paste("[[",name,"]]",sep=""))
}

### used in building tables to run the tests
### subsets the abundance table 
get_abund <- function(abund,type,db,level,name){
  sub <- abund[which(abund[,"type"] == type & 
                       abund[,"db"] == db &
                       abund[,"level"] == level &
                       abund[,"name"] == name),,drop=F]
  if( type == 'taxa' ){
    df <- sub[,c("sample_id","study_id","relative_abundance")]
    colnames(df) <- c("sample_id","study_id","relabun")
  }else if( type == 'function'){
    df <- sub[,c("sample_id","study_id","copy_number")]
    colnames(df) <- c("sample_id","study_id","copynum")
  }
  row.names(df) <- df$sample_id
  return(df)
}

# used to grab the subjects for the desired phenotype without any other phenotypes
# if you want other phenotypes to be included, provide a vector with the names of those 
# phenotypes to exclude_pheno ( ie. Obesity )
get_ids <- function(att.map,phenotype,exclude_pheno){
  morb   <- c("arthritis","carcinoma","cirrhosis","metformin","t2d",
              "obesity","crohns","ulcerative_colitis","t1d")
  comorb <- morb[!(morb %in% phenotype) & !(morb %in% exclude_pheno)]
  
  with_pheno <- !is.na(att.map[,phenotype])
  no_comorb  <- apply(is.na(att.map[,comorb]) | att.map[,comorb] == 0,1,all)
  return(att.map[with_pheno & no_comorb,"subject_id"])
}

get_obesity_ids <- function(att.map,phenotype,exclude_pheno){
  morb   <- c("arthritis","carcinoma","cirrhosis","metformin","t2d",
              "obesity","crohns","ulcerative_colitis","t1d")
  comorb <- morb[!(morb %in% phenotype) & !(morb %in% exclude_pheno)]
  
  with_pheno <- !is.na(att.map[,phenotype])
  no_comorb  <- apply(is.na(att.map[,comorb]) | att.map[,comorb] == 0,1,all)
  sub_att.map <- att.map[with_pheno & no_comorb,]
  lean_obese <- sub_att.map[which(sub_att.map[,"bmi_type"] %in% c("normal_weight","obese")),]
  return(lean_obese[,"subject_id"])
}

exclude_samp_abunds <- function(abunds,exclude_samples){
  #   exclude_rows <- NULL
  #   for( i in 1:nrow(exclude_samples) ){
  #     samp_id  <- exclude_samples[i,"sample_id"]
  #     study_id <- exclude_samples[i,"study_id"]
  # 
  #     rows <- row.names(abunds[which(abunds[,"sample_id"] == samp_id & abunds[,"study_id"] == study_id),])
  #     exclude_rows <- c(exclude_rows,rows)
  #   }
  #   keep_rows <- row.names(abunds)[!(row.names(abunds) %in% exclude_rows)]
  #   return(abunds[keep_rows,])
  
  ex_names    <- paste(exclude_samples[,"sample_id"],exclude_samples[,"study_id"],sep=".")
  abund_names <- paste(abunds[,"sample_id"],abunds[,"study_id"],sep=".")
  
  return(abunds[!(abund_names %in% ex_names),])
}
