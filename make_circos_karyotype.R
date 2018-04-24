#setwd("/Users/armourc/Documents/Sharpton_Lab/MetaAnalysis/new_database/")
setwd("/Users/armourc/Box Sync/Sharpton_Lab/MetaAnalysis/new_database/")

### FUNCTIONS
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

### INPUTS
results_wcorrection <- read.table("./stats_out/manuscript_figures/kegg_0.2/kegg_output_wcorrection.txt",
                                  header=T,stringsAsFactors = F,sep="\t")
detail_results <- read.table("./model_output/kegg_output_detail.txt",header=F,stringsAsFactors = F,sep="\t")
colnames(detail_results) <- c("type","database","level","name","phenotype","component","estimate","std.err","t_value","p_value")

shared_markers <- read.table("./stats_out/manuscript_figures/kegg_0.2/kegg_shared_markers.txt",
                             header=T,stringsAsFactors = F,sep="\t",quote = "")

phenotypes <- c("arthritis","carcinoma","cirrhosis","crohns",
                "obesity","t2d","ulcerative_colitis")

cutoff <- 0.2
level  <- "module"

colors <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7")
names(colors) <- phenotypes

### MAIN

karyotype <- NULL
for(phenotype in phenotypes){
  sub_results <- results_wcorrection[which(results_wcorrection$phenotype == phenotype &
                                           results_wcorrection$level == level),]
  
  sig_results <- sub_results[which(sub_results$fdr < cutoff),c("level","name","phenotype")]
  sig_results <- unique(sig_results)
  row.names(sig_results) <- sig_results$name
  
  # if(phenotype %in% c("arthritis","carcinoma")){
  #   #n_sig <- length(unique(sig_results$name))
  #   sig_results <- sig_results[which]
  # }
  # 
  n_sig <- nrow(sig_results)
  
  
  
  sub_shared <- shared_markers[which(shared_markers$name %in% sig_results$name),]
  name_order <- sub_shared$name[order(sub_shared$phenotype_count,decreasing = T)]
  
  sig_results <- sig_results[name_order,]
  
  
  data <- data.frame(band=rep("band",n_sig),
                     chr_id=get_p_abbrev(phenotype),
                     id=sig_results$name,
                     label=sig_results$name,
                     start=seq(0,n_sig-1,1),
                     end=seq(1,n_sig,1),
                     color=colors[phenotype])
  
  karyotype <- rbind(karyotype,data)
}

write.table(karyotype,"/Users/armourc/Desktop/karyotype.disease.moduleV2.txt",
            sep=" ",row.names = F,col.names = F,quote=F)

#### table for heatmap showing increased or decreased in cases
heatmap_info <- NULL
for(phenotype in phenotypes){
  sig_names <- karyotype[which(karyotype$chr_id == get_p_abbrev(phenotype)),"id"]
  pheno_detail_results <- detail_results[which(detail_results$phenotype == phenotype),]
  sub_pheno_detail_results <- pheno_detail_results[which(pheno_detail_results$name %in% sig_names),]
  if(phenotype %in% c("arthritis","carcinoma")){
    sub_pheno_detail_results <- sub_pheno_detail_results[which(sub_pheno_detail_results$component != "intercept"),]
    pheno_slope_est <- NULL
    for(name in unique(sub_pheno_detail_results[,"name"])){
      name_fdrs <- results_wcorrection[which(results_wcorrection$phenotype == phenotype &
                                             results_wcorrection$name == name &
                                             results_wcorrection$fdr < cutoff),]
      if(nrow(name_fdrs) == 1){
        plevel <- name_fdrs$phenotype_level
      }else{
        if(nrow(name_fdrs) > 2){
          stop("why so many rows?!?")
        }
        if(name_fdrs$fdr[1] < name_fdrs$fdr[2]){
          plevel <- name_fdrs$phenotype_level[1]
        }else{
          plevel <- name_fdrs$phenotype_level[2]
        }
      }
      if(plevel == "advanced_adenoma"){
        plevel <- "adenoma"
      }
      info <- sub_pheno_detail_results[which(sub_pheno_detail_results$name == name &
                                               sub_pheno_detail_results$component == paste0("status_",plevel)),c("name","estimate")]
      pheno_slope_est <- rbind(pheno_slope_est,info)
    }
  }else{
    pheno_slope_est <- sub_pheno_detail_results[which(sub_pheno_detail_results$component == "status"),c("name","estimate")]
  }
  colnames(pheno_slope_est) <- c("id","value")
  pheno_slope_est$value <- ifelse(pheno_slope_est$value > 0,1,-1)
  
  sub_karyotype <- karyotype[which(karyotype$chr_id == get_p_abbrev(phenotype)),]
  
  data <- merge(sub_karyotype[,c("chr_id","id","start","end")],pheno_slope_est,by="id")
  #cols <- ifelse(data$value == 1,"color=blue","color=red")
  cols <- ifelse(data$value == 1,"color=0,0,255,0.8","color=255,0,0,0.8")
  data <- data.frame(data[,c("chr_id","start","end","value")],options=cols)
  
  heatmap_info <- rbind(heatmap_info,data)
}
write.table(heatmap_info,"/Users/armourc/Desktop/slope_infoV3.txt",
            sep=" ",row.names = F,col.names = F,quote=F)


### info for links

link_markers <- shared_markers[which(shared_markers$phenotype_count > 1 &
                                     shared_markers$level == "module"),]
table <- NULL
for(i in 1:nrow(link_markers)){
  marker<- link_markers[i,"name"]
  phenos <- unlist(strsplit(link_markers[i,"phenotypes"],";"))
  
  if(length(phenos) == 5){
    option <- "color=black,thickness=5p" 
  }else if(length(phenos) == 4){
    option <- "color=dgrey,thickness=4p"
  }else if(length(phenos) == 3){
    option <- "color=lgrey,thickness=3p"
  }else if(length(phenos) == 2){
    option <- "color=vlgrey,thickness=3p"
  }else{
    stop("ERROR!!")
  }
  
  if(length(phenos) == 2){
    info <- c(get_p_abbrev(phenos[1]),
              karyotype[which(karyotype$id == marker & karyotype$chr_id == get_p_abbrev(phenos[1])),"start"],
              karyotype[which(karyotype$id == marker & karyotype$chr_id == get_p_abbrev(phenos[1])),"end"],
              get_p_abbrev(phenos[2]),
              karyotype[which(karyotype$id == marker & karyotype$chr_id == get_p_abbrev(phenos[2])),"start"],
              karyotype[which(karyotype$id == marker & karyotype$chr_id == get_p_abbrev(phenos[2])),"end"],
              option)
    table <- rbind(table,info)
  }else{
    if(length(phenos) == 2){
      info <- c(get_p_abbrev(phenos[1]),
                karyotype[which(karyotype$id == marker & karyotype$chr_id == get_p_abbrev(phenos[1])),"start"],
                karyotype[which(karyotype$id == marker & karyotype$chr_id == get_p_abbrev(phenos[1])),"end"],
                get_p_abbrev(phenos[2]),
                karyotype[which(karyotype$id == marker & karyotype$chr_id == get_p_abbrev(phenos[2])),"start"],
                karyotype[which(karyotype$id == marker & karyotype$chr_id == get_p_abbrev(phenos[2])),"end"],
                option)
      table <- rbind(table,info)
    }else{
      for(j in 1:(length(phenos)-1)){
        phenoA <- phenos[j]
        for(k in 2:length(phenos)){
          if(j == k){next()}
          phenoB <- phenos[k]
          
          info <- c(get_p_abbrev(phenoA),
                    karyotype[which(karyotype$id == marker & karyotype$chr_id == get_p_abbrev(phenoA)),"start"],
                    karyotype[which(karyotype$id == marker & karyotype$chr_id == get_p_abbrev(phenoA)),"end"],
                    get_p_abbrev(phenoB),
                    karyotype[which(karyotype$id == marker & karyotype$chr_id == get_p_abbrev(phenoB)),"start"],
                    karyotype[which(karyotype$id == marker & karyotype$chr_id == get_p_abbrev(phenoB)),"end"],
                    option)
          table <- rbind(table,info)
        }
      }
    }
  }
    
}
#write.table(table,"/Users/armourc/Desktop/marker_link_info.txt",
#            row.names = F,col.names = F,quote=F,sep=" ")

table <- as.data.frame(table)

table5 <- table[which(table$V7 == "color=black,thickness=5p"),]
table4 <- table[which(table$V7 == "color=dgrey,thickness=4p"),]
table3 <- table[which(table$V7 == "color=lgrey,thickness=3p"),]
table2 <- table[which(table$V7 == "color=vlgrey,thickness=3p"),]

write.table(table5,"/Users/armourc/Desktop/marker_link_info5.txt",
            row.names = F,col.names = F,quote=F,sep=" ")
write.table(table4,"/Users/armourc/Desktop/marker_link_info4.txt",
            row.names = F,col.names = F,quote=F,sep=" ")
write.table(table3,"/Users/armourc/Desktop/marker_link_info3.txt",
            row.names = F,col.names = F,quote=F,sep=" ")
write.table(table2,"/Users/armourc/Desktop/marker_link_info2.txt",
            row.names = F,col.names = F,quote=F,sep=" ")


