#####################
### SET VARIABLES ###
#####################
#setwd("/Users/armourc/Documents/Sharpton_Lab/MetaAnalysis/new_database/")
setwd("/Users/armourc/Box Sync/Sharpton_Lab/MetaAnalysis/new_database/")

cutoff  <- 0.2
outdir  <- paste0("./stats_out/manuscript_figures/kegg_",cutoff,"/")
if( !dir.exists(outdir) ){ dir.create(outdir,recursive = TRUE) }

source("./scripts/metagenome_analysis_functions.R")

##############
### INPUTS ###
##############

### METADATA
all_attributes <- read.table("./data/all_attributes.txt",header=T,stringsAsFactors = F,sep="\t")
# convert to table with ids as rows and attibutes as columns
att.map <- convert_attibutes(all_attributes)

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
a_table <- as.data.frame(matrix(data=NA,nrow=7,ncol=3,dimnames = list(phenotypes,c("level","p_value","rsq"))))

for(phenotype in phenotypes){
  mapping_subset <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == phenotype),]
  
  ids <- mapping_subset[,"sample_id"]
  sub_abunds <- exclude_zero(abund.map[ids,])
  a_test <- adonis(sub_abunds~status,data=mapping_subset,permutations = 1000)
  a_pval <- a_test$aov.tab["status","Pr(>F)"]
  rsq <- a_test$aov.tab["status","R2"]
  
  a_table[phenotype,] <- c(level,a_pval,rsq)
}
write.table(a_table,paste0(outpath,level,"_adonis.txt"),quote=F,sep="\t")

#ob and t2d adonis with study
sub_a_table <- NULL

for(phenotype in c("obesity","t2d")){
  mapping_subset <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == phenotype),]
  
  ids <- mapping_subset[,"sample_id"]
  sub_abunds <- exclude_zero(abund.map[ids,])
  a_test <- adonis(sub_abunds~status+study,data=mapping_subset,permutations = 1000)
  stat_pval <- a_test$aov.tab["status","Pr(>F)"]
  stat_rsq <- a_test$aov.tab["status","R2"]
  stud_pval <- a_test$aov.tab["study","Pr(>F)"]
  stud_rsq <- a_test$aov.tab["study","R2"]
    
  #sub_a_table[phenotype,] <- c(level,a_pval,rsq)
  sub_a_table <- rbind(sub_a_table,c(phenotype,level,"status",stat_pval,stat_rsq))
  sub_a_table <- rbind(sub_a_table,c(phenotype,level,"study",stud_pval,stud_rsq))
  
}
colnames(sub_a_table) <- c("phenotype","level","component","p_value","rsq")
write.table(sub_a_table,paste0(outpath,level,"_ob-t2d_adonis.txt"),quote=F,sep="\t")

#NMDS
for(a in 1:length(dims)){
  for(b in 1:length(dims)){
    dim.1 <- dims[a]
    dim.2 <- dims[b]
    if(dim.1 >= dim.2){
      next()
    }
    
    plist <- list()
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
        
        a_pval <- a_table
        
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
        
        plot <- ggplot(data=df,aes(x=df[,2],y=df[,3],colour=status,shape=status)) + 
          geom_point(size=1,alpha=0.5) +
          stat_ellipse(aes(colour=status,lty=status)) +
          scale_shape_manual(values = c(20,20)) + 
          scale_color_manual(values = c("blue","red")) +
          scale_linetype_manual(values = c(1,1)) +
          xlab(colnames(df)[2]) + ylab(colnames(df)[3]) +
          ggtitle(title) +
          guides(colour=guide_legend(paste0("Status")),
                 shape=guide_legend(paste0("Status")),
                 linetype=guide_legend(paste0("Status")))
        
        pdf(file=paste0(outpath,level,"_",phenotype,"_",dim.1,"-",dim.2,".pdf"))
        print(plot)
        dev.off()
        
        plist[[i]] <<- plot
      })
      
    }
    
    pdf(file=paste0(outpath,level,"_all_dim",dim.1,"_dim",dim.2,".pdf"),width = 10,height = 13)
    all_plot <- plot_grid(plotlist = plist,nrow = 4)
    print(all_plot)
    dev.off()
  }
}



### BETA DISPERSION

level <- "gene_family"
sub_abunds <- abunds[which(abunds$level == level),]
abund.map <- create_abund_map(sub_abunds,"kegg",level,"copy_number")

results_table <- matrix(nrow=length(phenotypes),ncol=2,dimnames = list(phenotypes,c("anova_p","ptest_p")))
p_list <- list()
for(pheno in phenotypes){
  local({
    i <- pheno
    
    ids <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == pheno),"sample_id"]
    case_ids <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == pheno &
                                              mod_disease_to_sample$status == 1),"sample_id"]
    sub_abund.map <- abund.map[ids,]
    
    dist <- vegdist(sub_abund.map,method="bray")
    group <- as.factor(ifelse(row.names(sub_abund.map) %in% case_ids,"case","control"))
    
    ## Calculate multivariate dispersions
    type <- "centroid"
    mod <- betadisper(dist,group,type=type,bias.adjust=TRUE)
    
    ## Perform test
    anova_result <- anova(mod)
    a_pval <- anova_result$`Pr(>F)`[1]
    
    ## Permutation test for F
    ptest_result <- permutest(mod, pairwise = TRUE)
    p_pval <- ptest_result$tab["Groups","Pr(>F)"]
    
    ## add values to table
    results_table[pheno,] <- c(a_pval,p_pval)
    
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
    
    pdf(paste0(outpath,level,"_",pheno,"_boxplot.pdf"))
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
  })
}
pdf(file=paste0(outdir,"bias_adj-",type,"/",level,"_all.pdf"),width=4,height=12)
all_plot <- plot_grid(plotlist = p_list,nrow = 4)
print(all_plot)
dev.off()

write.table(results_table,paste0(outdir,"bias_adj-",type,"/",level,"_pvals.txt"),quote=F,sep="\t")

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
markers_of_interest <- c("M00072","M00245","M00246",
                         "M00060","M00320","M00080",
                         "M00318",
                         "M00422","M00377","M00618","M00579",
                         "M00617","M00482",
                         "M00534","M00528","M00468",
                         "M00418",
                         "M00076","M00077","M00078","M00079",
                         "M00122","M00123","M00573",
                         "M00567","M00470","M00092",
                         "M00363","M00175",
                         "M00233","M00350","M00056",
                         "M00127")

sub_mod_df <- mod_df[which(mod_df$name %in% markers_of_interest),]
sub_mod_df$name <- factor(sub_mod_df$name,levels=rev(markers_of_interest))

neg_mid <- min(sub_mod_df$mod_est)/2
pos_mid <- max(sub_mod_df$mod_est)/2
#colours <- c("#770000", "red","white", "cyan", "#007777")
colours <- c("red3", "red","white", "deepskyblue", "deepskyblue4")
breaks  <- c(min(sub_mod_df$mod_est),neg_mid,0,pos_mid,max(sub_mod_df$mod_est))
labels <- c(min(sub_mod_df$mod_est),neg_mid,0,pos_mid,max(sub_mod_df$mod_est))
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


