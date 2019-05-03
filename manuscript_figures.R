#####################
### SET VARIABLES ###
#####################
#setwd("/Users/armourc/Documents/Sharpton_Lab/MetaAnalysis/new_database/")
setwd("/Users/armourc/Box Sync/Sharpton_Lab/MetaAnalysis/new_database/")

cutoff  <- 0.2
outdir  <- paste0("./stats_out/manuscript_figures/kegg_",cutoff,"/")
if( !dir.exists(outdir) ){ dir.create(outdir,recursive = TRUE) }

source("./scripts/metagenome_analysis_functions.R")

set.seed(98)
##############
### INPUTS ###
##############

### METADATA
all_attributes <- read.table("./data/all_attributes.txt",header=T,stringsAsFactors = F,sep="\t")
# convert to table with ids as rows and attibutes as columns
att.map <- convert_attibutes(all_attributes)

samp_to_subj <- read.table("./data/sample_to_subject.txt",header = T,stringsAsFactors = F,sep="\t")

mod_att.map <- merge(samp_to_subj,att.map,by="subject_id")

### MAPPING
disease_to_sample <- read.table("./data/disease_to_sample.txt",header = T,stringsAsFactors = F,sep = "\t")

### ABUNDANCE
abunds <- read.table(file = "./data/abundance_function_kegg.txt",
                     header = T,
                     stringsAsFactors = F,
                     sep = "\t")
abunds <- reduce_abunds(abunds,disease_to_sample)

#remove repeat samples
if( length(which(!(disease_to_sample[,"sample_id"] %in% abunds[,"sample_id"]))) >= 1 ){
  disease_to_sample <- reduce_mapping(abunds,disease_to_sample)
}
mod_disease_to_sample <- modify_mapping(disease_to_sample)

### MODEL OUTPUT
results <- read.table("./model_output/kegg_output.tab",header=F,stringsAsFactors = F,sep="\t")
colnames(results) <- c("type","database","level","name","phenotype","phenotype_level","studies","case_n","ctrl_n","case_mean","ctrl_mean","pvalue")

#detailed model output
detail_results <- read.table("./model_output/kegg_output_detail.txt",header=F,stringsAsFactors = F,sep="\t")
colnames(detail_results) <- c("type","database","level","name","phenotype","component","estimate","std.err","t_value","p_value")

### MISC
db_levels  <- c("gene_family","module","pathway")
phenotypes <- c("arthritis","carcinoma","cirrhosis","crohns",
                "obesity","t2d","ulcerative_colitis")

kegg_definitions <- read.table("./data/kegg_definitions.txt",header=T,stringsAsFactors = F,sep="\t",quote="",fill=F)
mod_kegg_definitions <- kegg_definitions
mod_kegg_definitions$definition <- unlist(lapply(kegg_definitions$definition,function(x) split_def(x)))
rownames(mod_kegg_definitions) <- mod_kegg_definitions$name

################
### ANALYSIS ###
################
# add fdr correction 
results_wcorrection <- add_correction(results)
write.table(results_wcorrection,file=paste0(outdir,"/kegg_output_wcorrection.txt"),quote=F,sep="\t")

# count biomarkers
cutoffs            <- c(0.05,0.10,0.15,0.20)
biomarker_counts   <- count_biomarkers(results_wcorrection,cutoffs)
write.table(biomarker_counts,file=paste0(outdir,"kegg_biomarker_counts.txt"),quote=F,sep="\t")

# find shared biomarkers at fdr cutoff
shared_markers <- find_common_biomarkers(results_wcorrection,cutoff,detail_results)
shared_markers <- merge(kegg_definitions[,c("name","definition")],shared_markers,by="name")
shared_markers <- shared_markers[,c("type","database","level","name","definition","phenotypes","fdr","change_in_case","phenotype_count")]
write.table(shared_markers,file=paste0(outdir,"/kegg_shared_markers.txt"),quote=F,sep="\t")

# count marker overlap between phenotypes
shared_counts <- count_shared_biomarkers(shared_markers)
write.table(shared_counts,file=paste0(outdir,"/kegg_shared_counts.txt"),quote=F,sep="\t")

### ALPHA DIVERSITY
out <- paste0(outdir,"alpha_diversity/")

#richness
outpath <- paste0(out,"richness/")
if(!dir.exists(outpath)){dir.create(outpath,recursive = T)}


gene_abund.map <- create_abund_map(abunds,"kegg","gene_family","relative_abundance")

alpha.div <- specnumber(gene_abund.map,MARGIN = 1)

max_a <- max(alpha.div)
min_a <- min(alpha.div)

df <- NULL
ttest_ps <- NULL
kstest_ps <- NULL
kwtest_ps <- NULL
table     <- as.data.frame(matrix(nrow=length(phenotypes),ncol=8,
                                  dimnames=list(phenotypes,
                                                c("control_mean","case_mean","ttest_pvalue","t_sig",
                                                  "kstest_pvalue","ks_sig","kwtest_pvalue","kw_sig"))))
plotlist  <- list()
plotlist1 <- list()
for(phenotype in phenotypes){
  sub_mapping <- mod_disease_to_sample[which(mod_disease_to_sample[,"phenotype"] == phenotype),]
  case_ids <- sub_mapping[which(sub_mapping[,"status"] == 1),"sample_id"]
  ctrl_ids <- sub_mapping[which(sub_mapping[,"status"] == 0),"sample_id"]
  
  case_div <- alpha.div[case_ids]
  ctrl_div <- alpha.div[ctrl_ids]
  
  case_mean <- mean(case_div)
  ctrl_mean <- mean(ctrl_div)
  
  case_data <- as.data.frame(cbind(rep(phenotype,length(case_ids)),rep(1,length(case_ids)),ids=case_ids,alpha=case_div))
  case_data$alpha<- as.numeric(as.character(case_data$alpha))
  
  ctrl_data <- as.data.frame(cbind(rep(phenotype,length(ctrl_ids)),rep(0,length(ctrl_ids)),ids=ctrl_ids,alpha=ctrl_div))
  ctrl_data$alpha <- as.numeric(as.character(ctrl_data$alpha))
  
  pheno_df <- as.data.frame(rbind(case_data,ctrl_data))
  colnames(pheno_df) <- c("phenotype","status","sample_id","alpha")
  pheno_df$alpha <- as.numeric(as.character(pheno_df$alpha))
  
  df <- as.data.frame(rbind(df,pheno_df))
  colnames(df) <- c("phenotype","status","sample_id","alpha")
  
  #t-test
  t_test <- t.test(alpha ~ status,data=pheno_df)
  ttest_ps <- c(ttest_ps,t_test$p.value)
  
  #ks test
  ks_test <- ks.test(case_data$alpha,ctrl_data$alpha)
  kstest_ps <- c(kstest_ps,ks_test$p.value)
  
  #kruskal wallis
  kw_test <- kruskal.test(alpha ~ status, data=pheno_df)
  kwtest_ps <- c(kwtest_ps,kw_test$p.value)
  
  table[phenotype,"case_mean"]     <- case_mean 
  table[phenotype,"control_mean"]  <- ctrl_mean
  table[phenotype,"ttest_pvalue"]  <- t_test$p.value
  table[phenotype,"kstest_pvalue"] <- ks_test$p.value
  table[phenotype,"kwtest_pvalue"] <- kw_test$p.value
  
  if(t_test$p.value < 0.05){
    t_significance <- "*"
    if(t_test$p.value < 0.01){
      t_significance <- "**"
      if(t_test$p.value < 0.001){
        t_significance <- "***" 
      }
    }
  }else{
    t_significance <- "NS"
  }
  table[phenotype,"t_sig"] <- t_significance
    
  if(ks_test$p.value < 0.05){
    k_significance <- "*"
    if(ks_test$p.value < 0.01){
      k_significance <- "**"
      if(ks_test$p.value < 0.001){
        k_significance <- "***" 
      }
    }
  }else{
    k_significance <- "NS"
  }
  table[phenotype,"ks_sig"] <- k_significance
  
  if(kw_test$p.value < 0.05){
    kw_significance <- "*"
    if(kw_test$p.value < 0.01){
      kw_significance <- "**"
      if(kw_test$p.value < 0.001){
        kw_significance <- "***" 
      }
    }
  }else{
    kw_significance <- "NS"
  }
  table[phenotype,"kw_sig"] <- kw_significance
  
  if(t_significance == "NS"){
    mod_sig <- ""
  }else{
    mod_sig <- t_significance    
  }
  
  #histogram
  plot <- ggplot(pheno_df,aes(x=alpha,colour=status,fill=status)) + geom_histogram(alpha=0.3) +
    scale_fill_manual(values = c("red","blue"),breaks=c(0,1),labels=c("control","case")) +
    scale_colour_manual(values = c("red","blue"),breaks=c(0,1),labels=c("control","case")) +
    ggtitle(paste0(get_p_name(phenotype)," ",mod_sig)) + 
    xlab(paste0("Number of KOs")) + ylab("Number of subjects") +
    scale_x_continuous(limits=c(3500,6200))
  plotlist[[phenotype]] <- plot
  
  plot <- ggplot(pheno_df,aes(x=alpha,clour=status,fill=status)) + geom_density(alpha=0.4) +
    scale_fill_manual(values = c("red","blue"),breaks=c(0,1),labels=c("control","case")) +
    scale_colour_manual(values = c("red","blue"),breaks=c(0,1),labels=c("control","case")) +
    ggtitle(paste0(get_p_name(phenotype)," ",mod_sig)) + 
    xlab(paste0("Number of KOs")) +
    scale_x_continuous(limits=c(3500,6200))
  plotlist1[[phenotype]] <- plot
}
write.table(table,file=paste0(outpath,"kegg_gene_richness_table.txt"),quote=F,sep="\t",row.names = T)

names(ttest_ps) <- phenotypes
t_sig <- ifelse(ttest_ps<0.05,"*","")

sig_pheno <- paste0(unlist(lapply(df$phenotype,function(x) get_p_abbrev(x))),t_sig[df$phenotype])
df$phenotype <- sig_pheno

colnames(df) <- c("phenotype","status","sample_id","alpha")
df <- as.data.frame(df)
df$alpha <- as.numeric(as.character(df$alpha))
df$status <- factor(df$status,levels=c(0,1),labels=c("control","case"))

#boxplot
outfile <- paste0(outpath,"kegg_gene_richness_boxplot.pdf")
pdf(file=outfile)
plot <- ggplot(data=df,aes(x=phenotype,y=alpha,fill=status)) + geom_boxplot(alpha=0.75) +
  ylab("number of KOs") + ggtitle("KEGG gene family richness") + xlab("") +
  scale_fill_manual(values=c("blue","red"))
print(plot)
dev.off()
#violin plot
outfile <- paste0(outpath,"kegg_gene_richness_violinplot.pdf")
pdf(file=outfile)
plot <- ggplot(data=df,aes(x=phenotype,y=alpha,fill=status)) + geom_violin(alpha=0.75) +
  stat_summary(fun.y=median, geom="point", size=1, color="black",position=position_dodge(1)) +
  scale_fill_manual(values=c("blue","red")) + xlab("") + ylab("number of KOs") +
  ggtitle("KEGG gene family richness")
print(plot)
dev.off()
#histogram
outfile <- paste0(outpath,"kegg_gene_richness_histogram.pdf") 
pdf(outfile,width=17)
plot <- plot_grid(plotlist = plotlist,nrow = 2)  
print(plot)
dev.off()
#density1
outfile <- paste0(outpath,"kegg_gene_richness_density.pdf")
pdf(outfile,width=16,height=7)
plot <- plot_grid(plotlist = plotlist1,nrow=2)
print(plot)
dev.off()
#density2
outfile <- paste0(outpath,"kegg_gene_richness_density2.pdf")
pdf(outfile,width=13,height=14)
plot <- plot_grid(plotlist = plotlist1,nrow=3)
print(plot)
dev.off()

### BETA DIVERSTIY
level <- "gene_family"
outpath <- paste0(outdir,"beta_diversity/")
if(!(dir.exists(outpath))){dir.create(outpath,recursive=T)}

dims <- c(1,2,3,4)

abund.map <- create_abund_map(abunds,"kegg","gene_family","copy_number")
dist <- metaMDS(abund.map,distance="bray",pc=TRUE,k=length(dims))

#adonis
# a_table <- as.data.frame(matrix(data=NA,nrow=7,ncol=3,dimnames = list(phenotypes,c("level","p_value","rsq"))))
# 
# for(phenotype in phenotypes){
#   mapping_subset <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == phenotype),]
#   
#   ids <- mapping_subset[,"sample_id"]
#   sub_abunds <- exclude_zero(abund.map[ids,])
#   a_test <- adonis(sub_abunds~status,data=mapping_subset,permutations = 1000)
#   a_pval <- a_test$aov.tab["status","Pr(>F)"]
#   rsq <- a_test$aov.tab["status","R2"]
#   
#   a_table[phenotype,] <- c(level,a_pval,rsq)
# }
#write.table(a_table,paste0(outpath,level,"_adonis.txt"),quote=F,sep="\t")

#adonis with metadata
a_table <- NULL
for(phenotype in phenotypes){
  mapping_subset <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == phenotype),]
  
  ids <- mapping_subset[,"sample_id"]
  sub_abunds <- exclude_zero(abund.map[ids,])
  
  metadata_subset <- mod_att.map[which(mod_att.map$sample_id %in% ids),
                                 c("subject_id","sample_id","age","bmi","country","sex")]
  mapping_subset2 <- merge(mapping_subset,metadata_subset,by="sample_id")
  mapping_subset2 <- na.omit(mapping_subset2)
  sub_abunds2 <- sub_abunds[mapping_subset2$sample_id,]
  
  a_test <- NA
  #print(phenotype)
  if(phenotype %in% c("obesity","t2d")){
    a_test <- adonis(sub_abunds2~status+age+bmi+sex+study,data=mapping_subset2,permutations=1000)
  }else{
    a_test <- adonis(sub_abunds2~status+age+bmi+sex,data=mapping_subset2,permutations = 1000)
  }
  variable <- row.names(a_test$aov.tab)[which(!row.names(a_test$aov.tab) %in% c("Residuals","Total"))]
  p_value <- a_test$aov.tab[variable,"Pr(>F)"]
  r_squared <- a_test$aov.tab[variable,"R2"]
  db_level <- rep(level,length(variable))
  disease <- rep(phenotype,length(variable))

  pheno_info <- cbind(db_level,disease,variable,p_value,r_squared)

  a_table <- rbind(a_table,pheno_info)
}
colnames(a_table) <- c("level","disease","variable","p_value","r_squared")
#write.table(a_table,paste0(outpath,level,"_adonis_with_metadata.txt"),quote=F,sep="\t")
a_table <- read.table(paste0(outpath,level,"_adonis_with_metadata.txt"),header=T,sep="\t")

#ob and t2d adonis with study
# sub_a_table <- NULL
# 
# for(phenotype in c("obesity","t2d")){
#   mapping_subset <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == phenotype),]
#   ids <- mapping_subset[,"sample_id"]
#   sub_abunds <- exclude_zero(abund.map[ids,])
#   a_test <- adonis(sub_abunds~status+study,data=mapping_subset,permutations = 1000)
#   
#   metadata_subset <- mod_att.map[which(mod_att.map$sample_id %in% ids),
#                                  c("subject_id","sample_id","age","bmi","country","sex")]
#   mapping_subset2 <- merge(mapping_subset,metadata_subset,by="sample_id")
#   mapping_subset2 <- na.omit(mapping_subset2)
#   sub_abunds2 <- sub_abunds[mapping_subset2$sample_id,]
#   a_test2 <- adonis(sub_abunds2~status+age+bmi+sex+country,data=mapping_subset2,permutations=1000)
#   a_test3 <- adonis(sub_abunds2~status+age+bmi+sex+country+study,data=mapping_subset2)
#   
#   stat_pval <- a_test$aov.tab["status","Pr(>F)"]
#   stat_rsq <- a_test$aov.tab["status","R2"]
#   stud_pval <- a_test$aov.tab["study","Pr(>F)"]
#   stud_rsq <- a_test$aov.tab["study","R2"]
#     
#   #sub_a_table[phenotype,] <- c(level,a_pval,rsq)
#   sub_a_table <- rbind(sub_a_table,c(phenotype,level,"status",stat_pval,stat_rsq))
#   sub_a_table <- rbind(sub_a_table,c(phenotype,level,"study",stud_pval,stud_rsq))
#   
# }
# colnames(sub_a_table) <- c("phenotype","level","component","p_value","rsq")
# write.table(sub_a_table,paste0(outpath,level,"_ob-t2d_adonis.txt"),quote=F,sep="\t")

#NMDS
for(a in 1:length(dims)){
  for(b in 1:length(dims)){
    dim.1 <- dims[a]
    dim.2 <- dims[b]
    if(dim.1 >= dim.2){
      next()
    }
    
    plist <- list()
    plist2 <- list()
    plistG <- list()
    #a_table <- as.data.frame(matrix(data=NA,nrow=7,ncol=3,dimnames = list(phenotypes,c("level","p_value","rsq"))))
    #a_table <- as.data.frame(matrix(data=NA,nrow=0,ncol=4))
    #alist <- list()
    for(phenotype in phenotypes){
      local({
        i <- phenotype
        
        filename <- paste0( outpath,"/nmds-",level,"_","_dims_", dim.1,"_", dim.2,".pdf")
        
        mapping_subset <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == phenotype),]
        sub_data <- dist$points
        sub_data <- sub_data[which(row.names(sub_data) %in% mapping_subset[,"sample_id"]),]
        
        #build dataframe
        df <- sub_data[,c(dim.1,dim.2)]
        df <- cbind(sample_id=row.names(df),df)
        df <- merge(df,mapping_subset,by="sample_id")
        
        #fix data classes
        df[,2] <- as.numeric(as.character(df[,2]))
        df[,3] <- as.numeric(as.character(df[,3]))
        df[,"status"] <- factor(df$status,labels = c("control","case"))
        
        #adonis
        # ids <- mapping_subset[,"sample_id"]
        # sub_abunds <- exclude_zero(abund.map[ids,])
        # a_test <- adonis(sub_abunds~status,data=mapping_subset,permutations = 5000)
        # a_pval <- a_test$aov.tab["status","Pr(>F)"]
        # rsq <- a_test$aov.tab["status","R2"]
        # 
        # alist[[i]] <<- c(level,a_pval,rsq)
        # a_table <<- rbind(a_table,c(phenotype,level,a_pval,rsq))
        # a_table[phenotype,] <- c(phenotype,level,a_pval,rsq)
        
        a_pval <- a_table[which(a_table$disease == phenotype & 
                                a_table$variable == "status"),"p_value"]
        
        if(a_pval < 0.05){
          significance <- "*"
          if(a_pval < 0.01){
            significance <- "**"
            if(a_pval < 0.001){
              significance <- "***" 
            }
          }
        }else{
          significance <- ""
        }
        
        #title
        title <- paste0(get_p_name(phenotype)," ",significance)
        
        plot <- ggplot(data=df,aes(x=df[,2],y=df[,3],colour=status,shape=status,fill=status)) + 
          geom_point(size=1,alpha=0.6) +
          stat_ellipse(aes(colour=status,lty=status)) +
          scale_shape_manual(values = c(21,21)) + 
          scale_color_manual(values = c("blue","red")) +
          scale_fill_manual(values = c("blue","red")) +
          scale_linetype_manual(values = c(1,1)) +
          xlab(colnames(df)[2]) + ylab(colnames(df)[3]) +
          ggtitle(title) +
          guides(colour=guide_legend(paste0("Status")),
                 shape=guide_legend(paste0("Status")),
                 linetype=guide_legend(paste0("Status")),
                 fill=guide_legend("Status"))
        
        pdf(file=paste0(outpath,level,"_",phenotype,"_",dim.1,"-",dim.2,".pdf"))
        print(plot)
        dev.off()
        
        plot <- plot + theme(legend.position="none")
        plist[[i]] <<- plot
        
        plot <- plot + theme(axis.ticks=element_blank(),
                             axis.text=element_blank())
       
        plist2[[i]] <<- plot
        
        plot <- plot + scale_colour_manual(values=c("black","Grey")) +
          scale_fill_manual(values=c("black","Grey"))
        plistG[[i]] <<- plot
      })
      
    }
    
    pdf(file=paste0(outpath,level,"_all_dim",dim.1,"_dim",dim.2,".pdf"),width = 10,height = 13)
    all_plot <- plot_grid(plotlist = plist,nrow = 4)
    print(all_plot)
    dev.off()
    
    pdf(file=paste0(outpath,level,"_all-wide_dim",dim.1,"_dim",dim.2,".pdf"),width = 12,height = 7)
    #png(filename=paste0(outpath,level,"_all-wide_dim",dim.1,"_dim",dim.2,".png"),width=650,height = 350,res=200)
    wide_plot <- plot_grid(plotlist = plist2,nrow=2)
    print(wide_plot)
    dev.off()
    
    pdf(file=paste0(outpath,level,"_all-wideG1_dim",dim.1,"_dim",dim.2,".pdf"),width = 12,height = 7)
    plot1 <- plot_grid(plotlist = list(plistG[["arthritis"]],plist2[["carcinoma"]],plist2[["cirrhosis"]],
                                       plist2[["crohns"]],plist2[["obesity"]],plist2[["t2d"]],plist2[["ulcerative_colitis"]])
                       ,nrow=2)
    print(plot1)
    dev.off()
    
    pdf(file=paste0(outpath,level,"_all-wideG2_dim",dim.1,"_dim",dim.2,".pdf"),width = 12,height = 7)
    plot2 <- plot_grid(plotlist = list(plistG[["arthritis"]],plistG[["carcinoma"]],plistG[["cirrhosis"]],
                                       plist2[["crohns"]],plistG[["obesity"]],plistG[["t2d"]],plistG[["ulcerative_colitis"]])
                       ,nrow=2)
    print(plot2)
    dev.off()
    
    pdf(file=paste0(outpath,level,"_all-wideG3_dim",dim.1,"_dim",dim.2,".pdf"),width = 12,height = 7)
    plot3 <- plot_grid(plotlist = list(plistG[["arthritis"]],plistG[["carcinoma"]],plistG[["cirrhosis"]],
                                       plist2[["crohns"]],plist2[["obesity"]],plistG[["t2d"]],plistG[["ulcerative_colitis"]])
                       ,nrow=2)
    print(plot3)
    dev.off()
    
  }
}



### BETA DISPERSION

level <- "gene_family"
sub_abunds <- abunds[which(abunds$level == level),]
abund.map <- create_abund_map(sub_abunds,"kegg",level,"copy_number")

results_table <- matrix(nrow=length(phenotypes),ncol=2,dimnames = list(phenotypes,c("anova_p","ptest_p")))
p_list <- list()
p_list2 <- list()
type <- "centroid"
for(pheno in phenotypes){
  local({
    i <- pheno
    
    ids <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == pheno),"sample_id"]
    case_ids <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == pheno &
                                              mod_disease_to_sample$status == 1),"sample_id"]
    sub_abund.map <- abund.map[ids,]
    
    dist2 <- vegdist(sub_abund.map,method="bray")
    group <- as.factor(ifelse(row.names(sub_abund.map) %in% case_ids,"case","control"))
    
    ## Calculate multivariate dispersions
    mod <- betadisper(dist2,group,type=type,bias.adjust=TRUE)
    
    ## Perform test
    anova_result <- anova(mod)
    a_pval <- anova_result$`Pr(>F)`[1]
    
    ## Permutation test for F
    ptest_result <- permutest(mod, pairwise = TRUE)
    p_pval <- ptest_result$tab["Groups","Pr(>F)"]
    
    ## add values to table
    results_table[pheno,] <<- c(a_pval,p_pval)
    
    outpath <- paste0(outdir,"bias_adj-",type,"/",pheno,"/")
    if(!dir.exists(outpath)){dir.create(outpath,recursive = T)}
    
    ## Tukey's Honest Significant Differences
    mod.HSD <- TukeyHSD(mod)
    
    pdf(paste0(outpath,level,"_",pheno,"_plot1.pdf"))
    print(plot(mod.HSD))
    dev.off()
    
    ## Plot the groups and distances to centroids on the
    ## first two PCoA axes
    pdf(paste0(outpath,level,"_",pheno,"_plot2.pdf"))
    print(plot(mod))
    dev.off()
    
    ## Draw a boxplot of the distances to centroid for each group
    pdf(paste0(outpath,level,"_",pheno,"_plot3.pdf"))
    print(boxplot(mod))
    dev.off()      
    
    dist_to_centroid <- mod$distances 
    classification <- factor(mod$group,levels=c("control","case"),labels=c("control","case"))
    
    df <- as.data.frame(matrix(nrow=length(dist_to_centroid),ncol=3))
    colnames(df) <- c("sampleID","distance","status")
    df$sampleID <- names(dist_to_centroid)
    df$distance <- dist_to_centroid
    df$status <- classification
    
    if(p_pval < 0.05){
      significance <- "*"
      if(p_pval < 0.01){
        significance <- "**"
        if(p_pval < 0.001){
          significance <- "***" 
        }
      }
    }else{
      significance <- ""
    } 
    
    title <- paste0(get_p_name(pheno)," ",significance)
    
    pdf(paste0(outpath,level,"_",pheno,"_boxplot.pdf"),height = 4,width = 4)
    plot <- ggplot(df,aes(x=status,y=distance,fill=status,colour=status)) + geom_boxplot(alpha=0.5) +
      ylab("Dispersion") + xlab("") +
      scale_fill_manual(values = c("blue","red")) +
      scale_colour_manual(values = c("blue","red")) +
      ggtitle(title)
    print(plot)
    dev.off()
    
    plot <- plot + theme(legend.position="none") +
      ggtitle("")
    
    p_list[[i]] <<- plot
    
    plot <- plot + theme(axis.ticks=element_blank(),
                         axis.text=element_blank())
    p_list2[[i]] <<- plot
  })
}
pdf(file=paste0(outdir,"bias_adj-",type,"/",level,"_all.pdf"),width=4,height=12)
all_plot <- plot_grid(plotlist = p_list,nrow = 4)
print(all_plot)
dev.off()

pdf(file=paste0(outdir,"bias_adj-",type,"/",level,"_all-wide.pdf"),width=5,height=7)
wide_plot <- plot_grid(plotlist=p_list,nrow=3)
print(wide_plot)
dev.off()

write.table(results_table,paste0(outdir,"bias_adj-",type,"/",level,"_pvals.txt"),quote=F,sep="\t")

full_list <- list()
for(pheno in phenotypes){
  full_list[[paste0(pheno,"_bdiv")]] <- plist2[[pheno]]
  full_list[[paste0(pheno,"_bdisp")]] <- p_list2[[pheno]]
}
pdf(file="./stats_out/manuscript_figures/kegg_0.2/beta-div-disp-wide.pdf",height = 6,width=13)
wide_plot <- plot_grid(plotlist = full_list,nrow=3,ncol=6,rel_widths = c(1.5,1,1.5,1,1.5,1))
print(wide_plot)
dev.off()

### HEATMAPS
outpath <- paste0(outdir,"heatmap/")
if(!dir.exists(outpath)){dir.create(outpath,recursive = T)}

#make matrix of phenotype by module model coef
level <- "module" 
sub_detail_results <- detail_results[which(detail_results[,"level"] == level),]

table <- sub_detail_results[which(grepl(sub_detail_results[,"component"],pattern = "status.*",perl=T)),c("name","phenotype","component","estimate")]
#modify phenotype to include status level
mod_pheno <- NULL
fdr       <- NULL
sig       <- NULL
mod_est   <- NULL
for(i in 1:nrow(table)){
  pheno       <- table[i,"phenotype"]
  pheno_level <- NA
  comp        <- table[i,"component"]
  name        <- unadjust_name(table[i,"name"])
  est         <- table[i,"estimate"]
  
  if(pheno %in% c("arthritis","carcinoma")){
    pheno_level <- (strsplit(comp,"_")[[1]])[2]
    mod_pheno <- c(mod_pheno,paste(pheno,pheno_level,sep="_"))
    
    if(pheno_level == "adenoma"){
      pheno_level <- "advanced_adenoma"
    }
    val <- results_wcorrection[which(results_wcorrection[,"level"] == level & 
                                       results_wcorrection[,"name"] == name &
                                       results_wcorrection[,"phenotype"] == pheno &
                                       results_wcorrection[,"phenotype_level"] == pheno_level),"fdr"]
  }else{
    mod_pheno <- c(mod_pheno,pheno)
    
    val_row <- results_wcorrection[which(results_wcorrection[,"level"] == level & 
                                           results_wcorrection[,"name"] == name &
                                           results_wcorrection[,"phenotype"] == pheno),]
    if( check_if_empty(val_row) ){
      val <- NA
    }else{
      val <- val_row[,"fdr"]
    }
  }
  
  fdr <- c(fdr,val)
  
  if(is.na(val)){
    sig <- c(sig,"NA")
  }else if(val < cutoff){
    sig <- c(sig,"*")
  }else{
    sig <- c(sig," ")
  }
  
  if(is.na(val)){
    mod_est <- c(mod_est,"NA")
  }else if(val < cutoff){
    mod_est <- c(mod_est,est)
  }else{
    mod_est <- c(mod_est,0)
  }
  
  
}
df <- data.frame(table,mod_pheno,fdr,sig,mod_est)
df$mod_est <- as.numeric(as.character(df$mod_est))
mod_est_binary <- unlist(lapply(df$mod_est,function(x) make_binary(x)))
  
df <- data.frame(df,mod_est_binary)

# select only one level for RA and CC
ra_df <- df[which(df$phenotype == "arthritis"),]
mod_ra_df <- NULL
for(name in unique(ra_df$name)){
  name_df <- ra_df[which(ra_df$name == name),]
  if(nrow(name_df) > 2){
    stop("why so many rows?!?")
  }
  if(name_df$fdr[1] < name_df$fdr[2]){
   mod_ra_df <- rbind(mod_ra_df,name_df[1,]) 
  }else{
    mod_ra_df <- rbind(mod_ra_df,name_df[2,])
  }
}
cc_df <- df[which(df$phenotype == "carcinoma"),]
mod_cc_df <- NULL
for(name in unique(cc_df$name)){
  name_df <- cc_df[which(cc_df$name == name),]
  if(nrow(name_df) > 2){
    stop("why so many rows?!?")
  }
  if(name_df$fdr[1] < name_df$fdr[2]){
    mod_cc_df <- rbind(mod_cc_df,name_df[1,]) 
  }else{
    mod_cc_df <- rbind(mod_cc_df,name_df[2,])
  }
}
other_df <- df[which(df$phenotype %in% c("cirrhosis","crohns","obesity","t2d","ulcerative_colitis")),]

mod_df <- as.data.frame(rbind(mod_ra_df,mod_cc_df,other_df))
mod_df$phenotype <- unlist(lapply(mod_df$phenotype,function(x) get_p_abbrev(x)))

coef.map        <- acast(mod_df,name ~ phenotype,value.var = "estimate")
mod_coef.map    <- acast(mod_df,name ~ phenotype,value.var = "mod_est")
binary_coef.map <- acast(mod_df,name ~ phenotype,value.var = "mod_est_binary")
fdr.map         <- acast(mod_df,name ~ phenotype,value.var = "fdr")

#make heatmap with estimate
min <- 0
for(i in 1:nrow(coef.map)){
  row_min <- min(na.omit(coef.map[i,]))
  if(row_min < min){
    min <- row_min
  }
}
max <- 0
for(i in 1:nrow(coef.map)){
  row_max <- max(na.omit(coef.map[i,]))
  if(row_max > max){
    max <- row_max
  }
}

my_palette <- colorRampPalette(c("red","white","blue"))(n = 299)
col_breaks = c(seq(min,-0.0011,length=145), 
               seq(-0.001,0.001,length=10),  
               seq(0.0011,max,length=145)) 

#pdf(paste0(outpath,level,"_coef.pdf"))
# print(heatmap.2(coef.map,
#                 col=my_palette,       
#                 breaks=col_breaks,
#                 trace="none",
#                 srtCol=45,
#                 dendrogram = "none",
#                 labRow = F))
#dev.off()

#heatmap with mod_estimate
pdf(paste0(outpath,"mod_estimate.pdf"))
print(heatmap.2(mod_coef.map,
                col=my_palette,
                #breaks=col_breaks,
                trace="none",
                labRow = F,
                dendrogram = "column",
                srtCol=360))
dev.off()

#binary estimate
pdf(paste0(outpath,"binary_mod_estimate.pdf"))
print(heatmap.2(binary_coef.map,
                col=c("red","white","blue"),
                #breaks=c(-1,0,1),
                trace="none",
                labRow = F,
                dendrogram = "column",
                srtCol=360,
                key=F))
dev.off()

#heatmap with highlighted markers
#common_markers <- shared_markers[which(shared_markers$level == "module" &
#                                         shared_markers$phenotype_count >= 4),"name"]

#markers_of_interest <- c(common_markers,"M00617","M00482","M00534","M00418","M00076","M00077","M00078","M00079",
#                         "M00363","M00233","M00350","M00056")

# markers_of_interest <- c("M00060","M000320","M00080",
#                          "M00422","M00377","M00618","M00579",
#                          "M00318","M00190","M00240","M00243","M00317","M00319",
#                          "M00072",
#                          "M00617","M00482",
#                          "M00534",
#                          "M00418",
#                          "M00076","M00077","M00078","M00079",
#                          "M00363",
#                          "M00233","M00350","M00056")
markers_of_interest <- c("M00127",
                         "M00072","M00245","M00246",
                         "M00060","M00320","M00080",
                         "M00318",
                         "M00422","M00377","M00618","M00579",
                         "M00617","M00482",
                         #"M00534",
                         "M00528","M00468",
                         "M00418",
                         "M00076","M00077","M00078","M00079",
                         #"M00122","M00123","M00573",
                         "M00567","M00470","M00092",
                         "M00363","M00175",
                         "M00233","M00350","M00056")
#markers_of_interest <- rev(markers_of_interest)

sub_mod_df <- mod_df[which(mod_df$name %in% markers_of_interest),]
sub_mod_df$name <- factor(sub_mod_df$name,levels=rev(markers_of_interest))

neg_mid <- min(sub_mod_df$mod_est)/2
pos_mid <- max(sub_mod_df$mod_est)/2
#colours <- c("#770000", "red","white", "cyan", "#007777")
colours <- c("red4","red3","red2","red","white","deepskyblue","deepskyblue2","deepskyblue3", "deepskyblue4")
breaks  <- c(min(sub_mod_df$mod_est),neg_mid,(neg_mid/2),((neg_mid/2)/2),0,((pos_mid/2)/2),(pos_mid/2),pos_mid,max(sub_mod_df$mod_est))
#labels <- c(min(sub_mod_df$mod_est),neg_mid,(neg_mid/2),0,(pos_mid/2),pos_mid,max(sub_mod_df$mod_est))
#values  <- rescale(breaks)

pdf(paste0(outpath,"common_unique_indicator.pdf"))
plot <- ggplot(sub_mod_df,aes(x=phenotype,y=name,fill=mod_est)) + 
  geom_tile(color="gray28",size=0.1) +
  geom_text(aes(label = sig)) +
  #scale_fill_gradientn(colours=colours,values=rescale(breaks),name="Estimate") +
  scale_fill_gradientn(colours=colours,values=rescale(breaks),breaks=c(-3,0,5,10,15,20),name="Estimate") +
  ylab("") + xlab("") 
print(plot)
dev.off()

# log_mod_est <- NULL
# for(i in 1:length(sub_mod_df$mod_est)){
#   val <- sub_mod_df$mod_est[i]
#   if(val == 0){
#     mod_val <- 0
#   }else if(val < 0){
#     mod_val <- log(abs(val))
#   }else{
#     mod_val <- log(val)
#   }
#   log_mod_est <- c(log_mod_est,mod_val)
# }
# sub_mod_df$log_mod_est <- log_mod_est

sub_mod_df <- merge(sub_mod_df,kegg_definitions,by="name")
sub_mod_df$definition <- unlist(lapply(sub_mod_df$definition,function(x) unlist(strsplit(x,","))[1])) 
  
sub_mod_df$definition <- factor(sub_mod_df$definition)

sub_mod_df$phenotype <- factor(sub_mod_df$phenotype)


#neg_mid <- min(sub_mod_df$log_mod_est)/2
#pos_mid <- max(sub_mod_df$log_mod_est)/2
#colours <- c("red3", "red","white", "deepskyblue", "deepskyblue4")
#breaks  <- c(min(sub_mod_df$log_mod_est),neg_mid,0,pos_mid,max(sub_mod_df$log_mod_est))
#labels <- c(min(sub_mod_df$log_mod_est),neg_mid,0,pos_mid,max(sub_mod_df$log_mod_est))

#write.table(sub_mod_df,"/Users/armourc/Desktop/sub_mod_df.txt",quote=F,sep="\t",row.names = F)
#sub_mod_df <- read.table("/Users/armourc/Desktop/sub_mod_df.txt",header=T,sep="\t")


#defs <- sub_mod_df$definition
#names(defs)<- sub_mod_df$name

#sub_mod_df$name <- factor(sub_mod_df$name,levels = markers_of_interest,labels = defs[markers_of_interest])
sub_mod_df$name <- factor(sub_mod_df$name,levels=markers_of_interest)

sub_mod_df$phenotype <- factor(sub_mod_df$phenotype,
                               levels=rev(c("RA","CC","LC","CD","OB","T2D","UC")),
                               labels=rev(c("Arthritis","Colorectal cancer","Liver cirrhosis","Crohn's disease","Obesity","Type II diabetes","Ulcerative colitis")))
  
pdf(paste0(outpath,"common_unique_indicator_wide.pdf"),width=10)
wide_plot <- ggplot(sub_mod_df,aes(x=name,y=phenotype,fill=mod_est)) + 
  geom_tile(color="gray28",size=0.1) +
  geom_text(aes(label = sig)) +
  coord_equal(ratio = 1.3) + 
  #scale_fill_gradientn(colours=colours,values=rescale(breaks),name="Estimate") +
  #scale_fill_gradientn(colours=colours,values=rescale(breaks),breaks=c(-3,-1.5,0,1.5,3),name="Estimate") +
  #scale_fill_gradientn(colours=colours,values=rescale(breaks),breaks=c(-3,0,5,10,15,20),name="Estimate") +
  scale_fill_gradientn(colours=colours,values=rescale(breaks),breaks=c(-3,-1.5,0,2,4),name="Estimate") +
  ylab("") + xlab("") + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size=12))
print(wide_plot)
dev.off()

# outfile <- paste0(outpath,level,"_coeff-heatmap.pdf")
# pdf(outfile)
#   plot <- ggplot(table,aes(x=mod_pheno,y=name,fill=estimate)) + geom_tile(color="gray28",size=0.1) + 
#     coord_equal() + geom_text(aes(label = sig)) +
#     scale_fill_gradientn(colours=colours,values=rescale(breaks)) +
#     #scale_fill_gradient2(low="#770000",mid="white",high="#007777",midpoint = 0) +
#     ggtitle(paste0(database," ",level,": coefficient")) + xlab("Phenotype") + ylab(level) +
#     theme(plot.title = element_text(size=20)) + 
#     theme(axis.title.x = element_text(size= 14)) + 
#     theme(axis.title.y = element_text(size= 14)) +
#     theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
#     scale_colour_discrete(guide=FALSE)
#   print(plot)


### MODULES

#common
outpath <- paste0(outdir,"common_markers/")
if(!dir.exists(outpath)){dir.create(outpath,recursive = T)}

common_markers <- shared_markers[which(shared_markers$level == "module" &
                                       shared_markers$phenotype_count >= 4),]

markers_of_interest <- c("M00060","M00072","M00377","M00618","M00318","M00190","M00240", "M00243", "M00317", "M00319")

abund.map <- create_abund_map(abunds,"kegg","module","copy_number")
for(marker in markers_of_interest){
  def <- mod_kegg_definitions[which(mod_kegg_definitions$name == marker),"definition"]
  
  sub_abund <- abund.map[,marker]
  df <- cbind(sample_id=names(sub_abund),relative_abundance=sub_abund)
  
  df <- as.data.frame(merge(mod_disease_to_sample,df,by="sample_id"))
  df$status <- factor(df$status,levels=c(0,1),labels=c("control","case"))
  df$relative_abundance<- as.numeric(as.character(df$relative_abundance))
  #df$phenotype <- factor(df$phenotype,levels=c("arthritis","carcinoma","cirrhosis","crohns","obesity","t2d","ulcerative_colitis"),
  #                       labels=c("Rheumatoid arthritis","Colorectal cancer","Liver cirrhosis","Crohn's disease","Obesity","Type II diabetes","Ulcerative colitis"))
  df$phenotype <- factor(df$phenotype,levels=c("arthritis","carcinoma","cirrhosis","crohns","obesity","t2d","ulcerative_colitis"),
                         labels=c("RA","CC","LC","CD","OB","T2D","UC"))
  sig_phenos <- shared_markers[which(shared_markers[,"name"] == marker),"phenotypes"]
  sig_phenos <- strsplit(sig_phenos,";")[[1]]
  
  xlabels <- NULL
  for(i in 1:length(phenotypes)){
    p <- phenotypes[i]
    pabrev <- get_p_abbrev(p)
    if(p %in% sig_phenos){
      pabrev <- paste(pabrev,"*",sep="")
    }
    xlabels <- c(xlabels,pabrev)
  }
  
  pdf(paste0(outpath,marker,".pdf"))
  plot <- ggplot(df,aes(x=phenotype,y=relative_abundance,fill=status,colour=status)) + geom_boxplot(alpha=0.6) +
    ylab("Copy number") + ggtitle(paste0(marker,": ",def)) +
    theme(#axis.text.x=element_text(angle=45,hjust=1),
      axis.title.x=element_blank()) + 
    scale_x_discrete(labels=xlabels) +
    scale_fill_manual(values=c("blue","red")) +
    scale_colour_manual(values=c("blue","red"))
  print(plot)
  dev.off()
}

#unique markers
#RA M00617

#####################
### RANDOM FOREST ###
#####################
library(randomForest)
library(pROC)

outpath <- paste0(outdir,"random_forest/")
if(!dir.exists(outpath)){dir.create(outpath,recursive = T)}

module_abunds <- abunds[which(abunds$level == "module"),]
module.df <- dcast(module_abunds,sample_id~name,value.var = "copy_number")
module.df <- merge(mod_disease_to_sample[,c("sample_id","status")],
                   module.df,
                   by="sample_id")
module.df <- unique(module.df)
rownames(module.df) <- module.df$sample_id
module.df <- module.df[,-1]
module.df$status <- as.factor(module.df$status)

module.rf <- randomForest(status ~ .,data=module.df,importance=TRUE,ntree=5000)
module.plot <- plot_importance(module.rf,shared_markers,"All Modules")
pdf(paste0(outpath,"all_module.pdf"),height = 11)
print(module.plot)
dev.off()

module1 <- shared_markers[which(shared_markers$level == "module" & shared_markers$phenotype_count >= 1),"name"]
module1.df <- module.df[,c("status",module1)]
module1.rf <- randomForest(status ~ .,data=module1.df,importance=TRUE,ntree=5000)
module1.plot <- plot_importance(module1.rf,shared_markers,"Modules for 1+ diseases")
pdf(paste0(outpath,"module1.pdf"),height = 11)
print(module1.plot)
dev.off()

module2 <- shared_markers[which(shared_markers$level == "module" & shared_markers$phenotype_count >= 2),"name"]
module2.df <- module.df[,c("status",module2)]
module2.rf <- randomForest(status ~ .,data=module2.df,importance=TRUE,ntree=5000)
module2.plot <- plot_importance(module2.rf,shared_markers,"Modules for 2+ disesases")
pdf(paste0(outpath,"module2.pdf"),height = 11)
print(module2.plot)
dev.off()

module3 <- shared_markers[which(shared_markers$level == "module" & shared_markers$phenotype_count >= 3),"name"]
module3.df <- module.df[,c("status",module3)]
module3.rf <- randomForest(status ~ .,data=module3.df,importance=TRUE,ntree=5000)
module3.plot <- plot_importance(module3.rf,shared_markers,"Modules for 3+ diseases")
pdf(paste0(outpath,"module3.pdf"),height = 11)
print(module3.plot)
dev.off()

module4 <- shared_markers[which(shared_markers$level == "module" & shared_markers$phenotype_count >= 4),"name"]
module4.df <- module.df[,c("status",module4)]
module4.rf <- randomForest(status ~ .,data=module4.df,importance=TRUE,ntree=5000)
module4.plot <- plot_importance(module4.rf,shared_markers,"Modules for 4+ diseases")
pdf(paste0(outpath,"module4.pdf"),height = 11)
print(module4.plot)
dev.off()

module.roc <- roc(module.df$status,module.rf$votes[,2])
module1.roc <- roc(module1.df$status,module1.rf$votes[,2])
module2.roc <- roc(module2.df$status,module2.rf$votes[,2])
module3.roc <- roc(module3.df$status,module3.rf$votes[,2])
module4.roc <- roc(module4.df$status,module4.rf$votes[,2])

pdf(paste0(outpath,"module_roc.pdf"))
plot(module.roc,col="red")
plot(module1.roc,col="orange",add=T)
plot(module2.roc,col="green",add=T)
plot(module3.roc,col="blue",add=T)
plot(module4.roc,col="purple",add=T)
legend(0.2,0.25,c("all","1+","2+","3+","4+"),
       col=c("red","orange","green","blue","purple"),
       lty=1)
dev.off()

#diseased vs non diseased - separated by disease
mod_module.df <- dcast(module_abunds,sample_id~name,value.var = "copy_number")
mod_module.df <- merge(mod_disease_to_sample[,c("sample_id","status","pheno.status")],
                       mod_module.df,
                       by="sample_id")
mod_module.df$pheno.status <- as.character(mod_module.df$pheno.status)
mod_module.df$pheno.status <- ifelse(mod_module.df$status == 0,0,mod_module.df$pheno.status)
mod_module.df <- unique(mod_module.df)
rownames(mod_module.df) <- mod_module.df$sample_id
mod_module.df <- mod_module.df[,-c(1,2)]
mod_module.df$pheno.status <- as.factor(mod_module.df$pheno.status)

mod_module.rf <- randomForest(pheno.status ~ .,data=mod_module.df,importance=TRUE,ntree=5000)

#by disease
p.list <- list()
colors <- c("red","orange","yellow","green","blue","purple","brown")
names(colors) <- phenotypes
#pdf(paste0(outpath,"module_by_disease_roc.pdf"))
table <- data.frame(matrix(nrow = length(phenotypes),ncol=3,dimnames = list(phenotypes,c("level","disease","AUC"))))
for(pheno in phenotypes){
  pheno_samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == pheno),"sample_id"]
  #sub_module.df <- module.df[pheno_samps,]
  sub_module.df <- module1.df[pheno_samps,] # only modules significant for 1+disease
  pheno_module.rf <- randomForest(status ~ .,data=sub_module.df,importance=TRUE,ntree=5000)
  print(pheno)
  print(pheno_module.rf)
  pheno.roc <- roc(sub_module.df$status,pheno_module.rf$votes[,2])
  if(pheno == "arthritis"){
    plot(pheno.roc,col=colors[pheno])
  }else{
    plot(pheno.roc,add=T,col=colors[pheno])
  }
  table[pheno,] <- c("module",pheno,pheno.roc$auc[1])
  p.list[[pheno]] <- pheno.roc
}
legend(0.35,0.35,phenotypes,
       col=colors,lty=1)
dev.off()

#write.table(table,paste0(outpath,"auc_values_bydisease.txt"),row.names = F,quote=F,sep="\t")

colors <-  c(rgb(128,0,128,maxColorValue = 255),rgb(128,0,0,maxColorValue = 255),rgb(192,104,22,maxColorValue = 255),
             rgb(226,168,43,maxColorValue = 255),rgb(99,130,31,maxColorValue = 255),rgb(0,124,123,maxColorValue = 255),
             rgb(126,83,0,maxColorValue = 255))
names(p.list) <- c("Rheumatoid arthritis","Colorectal cancer","Liver cirrhosis","Crohn's disease",
                   "Obesity","Type II diabetes","Ulcerative colitis")
plot <- ggroc(p.list,size=1.2) + geom_abline(intercept = 1,colour="grey") + coord_equal() + scale_color_manual(values=alpha(colors,0.7)) +
  xlab("Specificity") + ylab("Sensitivity")

pdf(paste0(outpath,"module_by_disease_roc_ggplot.pdf"))
print(plot)
dev.off()

#exclude adenoma
excl_samps <- disease_to_sample[which(disease_to_sample$phenotype == "carcinoma" & disease_to_sample$status == "advanced adenoma"),"sample_id"]
na_samps <- unique(mod_disease_to_sample[which(!(mod_disease_to_sample$sample_id %in% excl_samps)),"sample_id"])

na_module.df <- module.df[na_samps,]
na_module.df <- unique(na_module.df)

na_module.rf <- randomForest(status ~ .,data = na_module.df,importance=T,ntree=5000)

#exclude adenoma and RA
excl_samps2 <- c(excl_samps,unique(mod_disease_to_sample[which(mod_disease_to_sample$phenotype == "arthritis"),"sample_id"]))

nara_samps <- unique(mod_disease_to_sample[which(!(mod_disease_to_sample$sample_id %in% excl_samps2)),"sample_id"])

nara_module.df <- module.df[nara_samps,]

nara_module.rf <- randomForest(status ~ .,data = nara_module.df,importance=T,ntree=5000)


#lc and cd 
samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype %in% c("cirrhosis","crohns")),"sample_id"]
lc_cd_module.df <- module.df[samps,]
lc_cd_module.df <- unique(lc_cd_module.df)

lc_cd_module.rf <- randomForest(status ~ .,data = lc_cd_module.df,importance=T,ntree=5000)
lc_cd_module.roc <- roc(lc_cd_module.df$status,lc_cd_module.rf$votes[,2])

#plot
module.roc <- roc(module.df$status,module.rf$votes[,2])
na_module.roc <- roc(na_module.df$status,na_module.rf$votes[,2])
nara_module.roc <- roc(nara_module.df$status,nara_module.rf$votes[,2])
  
colors <- c("#F8766D","#00BA38","#619CFF")

plot3 <- ggroc(list("All diseases"=module.roc,"Exclude adenoma"=na_module.roc,"Exclude adenoma/RA"=nara_module.roc),size=1.2) + 
  coord_equal() + geom_abline(intercept = 1,colour="grey") + xlab("Specificity") + ylab("Sensitivity") + scale_color_manual(values=colors[c(1,2,3)]) + 
  #theme(legend.position = "none") + 
  theme(text = element_text(size=20),axis.text = element_text(size=16))

ggroc(list("all diseases" = na_module.roc),size=1.2) + 
  coord_equal() + geom_abline(intercept = 1,colour="grey") + xlab("Specificity") + ylab("Sensitivity") + scale_color_manual(values=colors[c(1,2,3)]) + 
  theme(text = element_text(size=20),axis.text = element_text(size=16))

#exclude RA and CC
alt_samps <- unique(mod_disease_to_sample[which(!(mod_disease_to_sample$phenotype %in% c("arthritis","carcinoma"))),"sample_id"])
alt_module.df <- module.df[alt_samps,]
alt_module.df <- unique(alt_module.df)

alt_module.rf <- randomForest(status ~ .,data=alt_module.df,importance=TRUE,ntree=5000)

alt_module1 <- shared_markers[which(shared_markers$level == "module" & shared_markers$phenotype_count >= 1),"name"]
alt_module1.df <- alt_module.df[,c("status",alt_module1)]
alt_module1.rf <- randomForest(status ~ .,data=alt_module1.df,importance=TRUE,ntree=5000)

alt_module2 <- shared_markers[which(shared_markers$level == "module" & shared_markers$phenotype_count >= 2),"name"]
alt_module2.df <- alt_module.df[,c("status",alt_module2)]
alt_module2.rf <- randomForest(status ~ .,data=alt_module2.df,importance=TRUE,ntree=5000)

alt_module3 <- shared_markers[which(shared_markers$level == "module" & shared_markers$phenotype_count >= 3),"name"]
alt_module3.df <- alt_module.df[,c("status",alt_module3)]
alt_module3.rf <- randomForest(status ~ .,data=alt_module3.df,importance=TRUE,ntree=5000)

alt_module4 <- shared_markers[which(shared_markers$level == "module" & shared_markers$phenotype_count >= 4),"name"]
alt_module4.df <- alt_module.df[,c("status",alt_module4)]
alt_module4.rf <- randomForest(status ~ .,data=alt_module4.df,importance=TRUE,ntree=5000)

alt_module.roc <- roc(alt_module.df$status,alt_module.rf$votes[,2])
alt_module1.roc <- roc(alt_module1.df$status,alt_module1.rf$votes[,2])
alt_module2.roc <- roc(alt_module2.df$status,alt_module2.rf$votes[,2])
alt_module3.roc <- roc(alt_module3.df$status,alt_module3.rf$votes[,2])
alt_module4.roc <- roc(alt_module4.df$status,alt_module4.rf$votes[,2])

pdf(paste0(outpath,"alt_module_roc.pdf"))
plot(alt_module.roc,col="red")
plot(alt_module1.roc,col="orange",add=T)
plot(alt_module2.roc,col="green",add=T)
plot(alt_module3.roc,col="blue",add=T)
plot(alt_module4.roc,col="purple",add=T)
legend(0.2,0.25,c("all","1+","2+","3+","4+"),
       col=c("red","orange","green","blue","purple"),
       lty=1)
dev.off()

#KOs
ko_abunds <- abunds[which(abunds$level == "gene_family"),]
ko.df <- dcast(ko_abunds,sample_id~name,value.var = "copy_number")
ko.df <- merge(mod_disease_to_sample[,c("sample_id","status")],
                   ko.df,
                   by="sample_id")
ko.df <- unique(ko.df)
rownames(ko.df) <- ko.df$sample_id
ko.df <- ko.df[,-1]
ko.df$status <- as.factor(ko.df$status)

ko.rf <- randomForest(status ~ .,data=ko.df,importance=TRUE,ntree=5000)

ko1 <- shared_markers[which(shared_markers$level == "gene" & shared_markers$phenotype_count >= 1),"name"]
ko1.df <- ko.df[,c("status",ko1)]
ko1.rf <- randomForest(status ~ .,data=ko1.df,importance=TRUE,ntree=5000)

ko2 <- shared_markers[which(shared_markers$level == "gene" & shared_markers$phenotype_count >= 2),"name"]
ko2.df <- ko.df[,c("status",ko2)]
ko2.rf <- randomForest(status ~ .,data=ko2.df,importance=TRUE,ntree=5000)

ko3 <- shared_markers[which(shared_markers$level == "gene" & shared_markers$phenotype_count >= 3),"name"]
ko3.df <- ko.df[,c("status",ko3)]
ko3.rf <- randomForest(status ~ .,data=ko3.df,importance=TRUE,ntree=5000)

ko4 <- shared_markers[which(shared_markers$level == "gene" & shared_markers$phenotype_count >= 4),"name"]
ko4.df <- ko.df[,c("status",ko4)]
ko4.rf <- randomForest(status ~ .,data=ko4.df,importance=TRUE,ntree=5000)

ko.roc <- roc(ko.df$status,ko.rf$votes[,2])
ko1.roc <- roc(ko1.df$status,ko1.rf$votes[,2])
ko2.roc <- roc(ko2.df$status,ko2.rf$votes[,2])
ko3.roc <- roc(ko3.df$status,ko3.rf$votes[,2])
ko4.roc <- roc(ko4.df$status,ko4.rf$votes[,2])

pdf(paste0(outpath,"ko_roc.pdf"))
plot(ko.roc,col="red")
plot(ko1.roc,col="orange",add=T)
plot(ko2.roc,col="green",add=T)
plot(ko3.roc,col="blue",add=T)
plot(ko4.roc,col="purple",add=T)
legend(0.2,0.25,c("all","1+","2+","3+","4+"),
       col=c("red","orange","green","blue","purple"),
       lty=1)
dev.off()


#additional metadata
module_abunds <- abunds[which(abunds$level == "module"),]
md_module.df <- dcast(module_abunds,sample_id~name,value.var = "copy_number")
md_module.df <- merge(mod_disease_to_sample[,c("sample_id","status")],
                      md_module.df,by="sample_id")
md_module.df <- merge(mod_att.map[,c("sample_id","bmi","age","sex","country")],
                      md_module.df, by="sample_id")
md_module.df <- unique(md_module.df)
rownames(md_module.df) <- md_module.df$sample_id
md_module.df <- md_module.df[,-1]
md_module.df$status <- as.factor(md_module.df$status)
md_module.df$sex <- as.factor(md_module.df$sex)
md_module.df$country <- as.factor(md_module.df$country)

md_module.rf <- randomForest(status ~ .,data=md_module.df,importance=TRUE,ntree=5000,na.action = na.omit)
md_module.plot <- plot_importance(md_module.rf,shared_markers,"All Modules with metadata")
pdf(paste0(outpath,"all_module_metadata.pdf"),height = 11)
print(md_module.plot)
dev.off()

md_module.roc <- roc(na.omit(md_module.df)[,"status"],md_module.rf$votes[,2])

# metadata no BMI
module_abunds <- abunds[which(abunds$level == "module"),]
md2_module.df <- dcast(module_abunds,sample_id~name,value.var = "copy_number")
md2_module.df <- merge(mod_disease_to_sample[,c("sample_id","status")],
                      md2_module.df,by="sample_id")
md2_module.df <- merge(mod_att.map[,c("sample_id","age","sex","country")],
                      md2_module.df, by="sample_id")
md2_module.df <- unique(md2_module.df)
rownames(md2_module.df) <- md2_module.df$sample_id
md2_module.df <- md2_module.df[,-1]
md2_module.df$status <- as.factor(md2_module.df$status)
md2_module.df$sex <- as.factor(md2_module.df$sex)
md2_module.df$country <- as.factor(md2_module.df$country)

md2_module.rf <- randomForest(status ~ .,data=md2_module.df,importance=TRUE,ntree=5000,na.action = na.omit)
#md2_module.plot <- plot_importance(md2_module.rf,shared_markers,"All Modules with metadata (no BMI)")
#pdf(paste0(outpath,"all_module_metadata2.pdf"),height = 11)
#print(md2_module.plot)
#dev.off()

md2_module.roc <- roc(na.omit(md2_module.df)[,"status"],md2_module.rf$votes[,2])

plot(md_module.roc,col="blue")
plot(md2_module.roc,col="red",add=T)
legend(0.3,0.2,c("all_metadata","no BMI"),col=c("blue","red"),lty=1)

#ROC Plot comparisons
pdf(paste0(outpath,"compare_roc.pdf"))
plot(module.roc,col="red")
plot(ko.roc,col="orange",add=T)
plot(alt_module.roc,col="green",add=T)
plot(md_module.roc,col="blue",add=T)
plot(md2_module.roc,col="purple",add=T)
legend(0.55,0.25,c("module-all diseases","ko-all diseases","module-exclude RA and CC","module-with metadata","module-with metadata no BMI"),
       col=c("red","orange","green","blue","purple"),lty=1)
dev.off()

### exclude RA and CC with metadata
alt_samps <- unique(mod_disease_to_sample[which(!(mod_disease_to_sample$phenotype %in% c("arthritis","carcinoma"))),"sample_id"])

module_abunds <- abunds[which(abunds$level == "module"),]
altmd2_module.df <- dcast(module_abunds,sample_id~name,value.var = "copy_number")
altmd2_module.df <- merge(mod_disease_to_sample[,c("sample_id","status")],
                       altmd2_module.df,by="sample_id")
altmd2_module.df <- merge(mod_att.map[,c("sample_id","age","sex","country")],
                       altmd2_module.df, by="sample_id")
altmd2_module.df <- unique(altmd2_module.df)
rownames(altmd2_module.df) <- altmd2_module.df$sample_id
altmd2_module.df <- altmd2_module.df[,-1]
altmd2_module.df$status <- as.factor(altmd2_module.df$status)
altmd2_module.df$sex <- as.factor(altmd2_module.df$sex)
altmd2_module.df$country <- as.factor(altmd2_module.df$country)
altmd2_module.df <- altmd2_module.df[alt_samps,]

altmd2_module.rf <- randomForest(status ~ .,data=altmd2_module.df,importance=TRUE,ntree=5000,na.action = na.omit)
altmd2_module.plot <- plot_importance(altmd2_module.rf,shared_markers,"Exclude CC/RA + metadata")
pdf(paste0(outpath,"all_module_metadata2.pdf"),height = 11)
print(altmd2_module.plot)
dev.off()

altmd2_module.roc <- roc(na.omit(altmd2_module.df)[,"status"],altmd2_module.rf$votes[,2])

pdf(paste0(outpath,"compare2_roc.pdf"))
plot(module.roc,col="purple")
plot(alt_module.roc,col="green",add=T)
plot(md2_module.roc,col="blue",add=T)
plot(altmd2_module.roc,col="red",add=T)
legend(0.6,0.2,c("All diseases","Exclude RA and CC","all_metadata","no CC/RA - with metadata"),
       col=c("purple","green","blue","red"),lty=1)
dev.off()

pdf(paste0(outpath,"compare2_roc2.pdf"))
plot(module.roc,col="green")
#plot(alt_module.roc,col="green",add=T)
plot(md2_module.roc,col="blue",add=T)
plot(altmd2_module.roc,col="red",add=T)
#legend(0.6,0.2,c("All diseases","Exclude RA and CC","all_metadata","no CC/RA - with metadata"),
#       col=c("purple","green","blue","red"),lty=1)
legend(0.6,0.2,c("All diseases","All diseases + metadata","Exclude CC/RA + metadata"),
       col=c("green","blue","red"),lty=1)
dev.off()

colors <- c("#F8766D","#00BA38","#619CFF")
  
plot1 <- ggroc(list("All diseases"=module.roc),size=1.2) + coord_equal() + geom_abline(intercept = 1,colour="grey") + xlab("Specificity") + ylab("Sensitivity") + scale_color_manual(values=colors[1]) + theme(legend.position = "none") + theme(text = element_text(size=20),axis.text = element_text(size=16))
plot2 <- ggroc(list("All diseases"=module.roc,"All diseases + metadata"=md2_module.roc),size=1.2) + coord_equal() + geom_abline(intercept = 1,colour="grey") + xlab("Specificity") + ylab("Sensitivity") + scale_color_manual(values=colors[c(1,2)]) + theme(legend.position = "none") + theme(text = element_text(size=20),axis.text = element_text(size=16))
plot3 <- ggroc(list("All diseases"=module.roc,"All diseases + metadata"=md2_module.roc,"Exclude CC/RA + metadata"=altmd2_module.roc),size=1.2) + coord_equal() + geom_abline(intercept = 1,colour="grey") + xlab("Specificity") + ylab("Sensitivity") + scale_color_manual(values=colors[c(1,2,3)]) + theme(legend.position = "none") + theme(text = element_text(size=20),axis.text = element_text(size=16))

pdf(paste0(outpath,"compare2_roc2_ggplot1.pdf"))
print(plot1)
dev.off()
pdf(paste0(outpath,"compare2_roc2_ggplot2.pdf"))
print(plot2)
dev.off()
pdf(paste0(outpath,"compare2_roc2_ggplot3.pdf"))
print(plot3)
dev.off()

# b_samps <- c('a187','a188','a189','a190','a191','a192','a193','a194','a195','a196','a197',
# 'a198','a199','a200','a201','a202','a203','a204','a205','a206','a207','a208','a209','a210',
# 'a211','a212','a213','a214','a215','a216','a217','a218','a219','a220','a221','a222','a223',
# 'a224','a225','a226','a228','a229','a230','a231','a232','a233','a234','a235','a236','a237',
# 'a238','a239','a240','a241','a242','a243','a244','a245','a246','a247','a248','a249','a250',
# 'a251','a252','a253','a254','a255','a256')


#table of AUC values
table <- data.frame(matrix(nrow = 4,ncol=3))
table[1,] <- c("Module","all diseases",module.roc$auc[1])
table[2,] <- c("Module","exclude RA and CC",alt_module.roc$auc[1])
table[3,] <- c("Module","with metadata (no BMI)",md2_module.roc$auc[1])
table[4,] <- c("Module","exclude RA and CC with metadata",altmd2_module.roc$auc[1])
colnames(table) <- c("level","name","AUC")
write.table(table,paste0(outpath,"auc_values.txt"),row.names = F,quote=F,sep="\t")

######################################
### BOXPLOTS (ob and t2d by study) ###
######################################
outpath <- paste0(outdir,"boxplot_by_study/")
if(!dir.exists(outpath)){dir.create(outpath)}

#module_abund.map <- create_abund_map(abunds,"kegg","module","copy_number")
module_abunds <- abunds[which(abunds$level == "module"),]

#obesity
if(!dir.exists(paste0(outpath,"obesity/log/"))){dir.create(paste0(outpath,"obesity/log/"),recursive = T)}

ob_samps <- disease_to_sample[which(disease_to_sample$phenotype == "obesity"),]
colnames(ob_samps)[2] <- "study_id"
ob_abunds <- module_abunds[which(module_abunds$sample_id %in% ob_samps$sample_id),]

ob.df <- merge(ob_samps,ob_abunds[,c("sample_id","study_id","name","copy_number")],by=c("sample_id","study_id"))

for(name in unique(ob.df$name)){
  name.df <- ob.df[which(ob.df$name == name),]
  ob_fdr <- results_wcorrection[which(results_wcorrection$name == name & 
                                      results_wcorrection$phenotype == "obesity"),"fdr"]
  ob_sig <- ifelse(ob_fdr < cutoff,"*","")
  
  # name_detail <- detail_results[which(detail_results$name == name & detail_results$phenotype == "obesity"),]
  # study_sig <- name_detail$p_value
  # names(study_sig) <- gsub("study","",name_detail$component)
  # study_sig <- ifelse(study_sig < 0.05,"*","")
  
  title <- paste0(name,ob_sig,": ",mod_kegg_definitions[name,"definition"])
  plot1 <- ggplot(name.df,aes(x=study_id,y=copy_number,
                     fill=factor(status,levels = c(0,1),labels = c("control","case")))) + geom_boxplot() +
    ggtitle(title) + xlab("Study") + ylab("Copy Number") +
    guides(fill=guide_legend(paste0("Status"))) +
    scale_fill_manual(values=c("blue","red"))
  pdf(file=paste0(outpath,"obesity/",name,".pdf"))
  print(plot1)
  dev.off()
  
  plot2 <- ggplot(name.df,aes(x=study_id,y=log(copy_number+0.01),fill=factor(status,levels = c(0,1),labels = c("control","case")))) + 
    geom_boxplot() + ggtitle(title) + xlab("Study") + ylab("log(Copy Number + 0.01)") +
    guides(fill=guide_legend(paste0("Status"))) +
    scale_fill_manual(values=c("blue","red"))
  pdf(file=paste0(outpath,"obesity/log/",name,"_log.pdf"))
  print(plot2)
  dev.off()
  
}

#type II diabetes
if(!dir.exists(paste0(outpath,"t2d/log/"))){dir.create(paste0(outpath,"t2d/log/"),recursive = T)}

t2d_samps <- disease_to_sample[which(disease_to_sample$phenotype == "t2d"),]
colnames(t2d_samps)[2] <- "study_id"
t2d_abunds <- module_abunds[which(module_abunds$sample_id %in% t2d_samps$sample_id),]

t2d.df <- merge(t2d_samps,t2d_abunds[,c("sample_id","study_id","name","copy_number")],by=c("sample_id","study_id"))
t2d.df$study_id <- ifelse(t2d.df$study_id %in% c("IGC","Richness"),"MetaHIT",t2d.df$study_id)

t2d_results <-results_wcorrection[which(results_wcorrection$level == "module" & 
                                          results_wcorrection$phenotype == "t2d"),]
t2d_sig_results <- t2d_results[which(t2d_results$fdr < cutoff),]
t2d_sig_results <- t2d_sig_results[,c("name","phenotype","fdr")]
t2d_sig_results <- merge(t2d_sig_results,shared_markers[,c("name","definition","phenotypes")],by="name")

for(name in unique(t2d.df$name)){
  name.df <- t2d.df[which(t2d.df$name == name),]
  t2d_fdr <- results_wcorrection[which(results_wcorrection$name == name & 
                                        results_wcorrection$phenotype == "t2d"),"fdr"]
  t2d_sig <- ifelse(t2d_fdr < cutoff,"*","")

  title <- paste0(name,t2d_sig,": ",mod_kegg_definitions[name,"definition"])
  plot1 <- ggplot(name.df,aes(x=study_id,y=copy_number,
                              fill=factor(status,levels = c(0,1),labels = c("control","case")))) + geom_boxplot() +
    ggtitle(title) + xlab("Study") + ylab("Copy Number") +
    guides(fill=guide_legend(paste0("Status"))) +
    scale_fill_manual(values=c("blue","red"))
  pdf(file=paste0(outpath,"t2d/",name,".pdf"))
  print(plot1)
  dev.off()
  
  plot2 <- ggplot(name.df,aes(x=study_id,y=log(copy_number+0.01),fill=factor(status,levels = c(0,1),labels = c("control","case")))) + 
    geom_boxplot() + ggtitle(title) + xlab("Study") + ylab("log(Copy Number + 0.01)") +
    guides(fill=guide_legend(paste0("Status"))) +
    scale_fill_manual(values=c("blue","red"))
  pdf(file=paste0(outpath,"t2d/log/",name,"_log.pdf"))
  print(plot2)
  dev.off()
  
}

#######################
### REPEAT SUBJECTS ###
#######################
repeat_subj <- c("31071","31267","31866","31873","530075","530258","530295","530297","530331","530373","530394","530450",
                 "530558","530600","530623","530743","530840","530867","531128","531155","531248","531274","531281","531333",
                 "531352","531361","531382","531403","531424","531469","531663","532749","532779","532796","532802","547",
                 "551","552","577","59","MH0401","MH0403","MH0404","MH0405","MH0406","MH0407","MH0408","MH0409",
                 "MH0410","MH0411","MH0412","MH0413","MH0414","MH0415","MH0416","MH0417","MH0418","MH0420","MH0421","MH0422",
                 "MH0423","MH0424","MH0425","MH0426","MH0427","MH0428","MH0429","MH0430","MH0431","MH0432","MH0433","MH0434",
                 "MH0435","MH0436","MH0437","MH0438","MH0439","MH0440","MH0441","MH0442","MH0443","MH0444","MH0445","MH0446",
                 "MH0447","MH0448","MH0450","MH0451","MH0452","MH0453","MH0454","MH0455","MH0456","MH0457","O2.UC20","O2.UC22",
                 "O2.UC23","O2.UC23","O2.UC24","O2.UC24","O2.UC24","O2.UC26","O2.UC29","O2.UC29","O2.UC30","O2.UC31","O2.UC37","O2.UC37",
                 "O2.UC37","O2.UC53","O2.UC54","O2.UC54","O2.UC55","O2.UC55","O2.UC56","O2.UC57","O2.UC58","O2.UC59","O2.UC60","SZEY-103A",
                 "SZEY-104A")

########################
### QUANTIFY OVERLAP ###
########################
library(BSDA)

all_modules <- unique(results_wcorrection[which(results_wcorrection$level == "module"),"name"])
ind_modules <- unique(results_wcorrection[which(results_wcorrection$level == "module" & 
                                                results_wcorrection$fdr < cutoff),"name"])
module_sig_results_wcorrection <- results_wcorrection[which(results_wcorrection$level == "module" &
                                                              results_wcorrection$fdr < cutoff),]
disease_modules <- list(arthritis=module_sig_results_wcorrection[which(module_sig_results_wcorrection$phenotype == "arthritis"),"name"],
                        carcinoma=module_sig_results_wcorrection[which(module_sig_results_wcorrection$phenotype == "carcinoma"),"name"],
                        cirrhosis=module_sig_results_wcorrection[which(module_sig_results_wcorrection$phenotype == "cirrhosis"),"name"],
                        crohns=module_sig_results_wcorrection[which(module_sig_results_wcorrection$phenotype == "crohns"),"name"],
                        obesity=module_sig_results_wcorrection[which(module_sig_results_wcorrection$phenotype == "obesity"),"name"],
                        t2d=module_sig_results_wcorrection[which(module_sig_results_wcorrection$phenotype == "t2d"),"name"],
                        ulcerative_colitis=module_sig_results_wcorrection[which(module_sig_results_wcorrection$phenotype == "ulcerative_colitis"),"name"])
n_disease_modules <- unlist(lapply(disease_modules,function(x) length(x)))

permutations <- 100000
table <- NULL
for(i in 1:permutations){
  random_indicators_list <- NULL
  for(phenotype in phenotypes){
    n_indicators <- length(disease_modules[[phenotype]])
    
    random_indicators <- all_modules[sample(1:length(all_modules),n_indicators)]
    
    random_indicators_list[[phenotype]] <- random_indicators
    
  }
  random_indicators_vector <- unlist(random_indicators_list)
  random_indicators_counts <- table(random_indicators_vector)
  n_2_plus <- length(which(random_indicators_counts >= 2))
  
  n_indicators <- length(unique(random_indicators_vector))
  percent_2_plus <- n_2_plus/n_indicators
  #percent_2_plus <- n_2_plus/length(all_modules)
  
  #percent_overlap <- c(percent_overlap,percent_2_plus)
  table <- rbind(table,c(n_2_plus,n_indicators))
}
table <- as.data.frame(table)
colnames(table) <- c("n_shared","n_total")
table$percent_shared <- table$n_shared/table$n_total

times_above0.78 <- nrow(table[which(table$percent_shared >= 0.78),])
pval <- times_above0.78/permutations

obs_percent_2plus <- length(shared_markers[which(shared_markers$level == "module" & shared_markers$phenotype_count >=2),"name"])/length(ind_modules)

mu <- mean(table$percent_shared)
sig <- sd(table$percent_shared)
hist(table$percent_shared)

#z <- (mu - obs_percent_2plus)/sqrt((obs_percent_2plus*(1-obs_percent_2plus))/permutations)
#CI95 <- c(mu - abs((z*sqrt((obs_percent_2plus*(1-obs_percent_2plus))/permutations))),mu + abs(z*sqrt((obs_percent_2plus*(1-obs_percent_2plus))/permutations)))

#test <- z.test(table$percent_shared,mu = obs_percent_2plus,sigma.x = sig)

# df <- as.data.frame(percent_overlap)
# ggplot(df,aes(x=percent_overlap)) + geom_histogram(binwidth = 0.005) + geom_vline(xintercept=obs_percent_2plus) + 
#   geom_vline(xintercept = mean(percent_overlap)) + xlim(c(0,1))
