setwd("/Users/armourc/Box Sync/Sharpton_Lab/MetaAnalysis/new_database/")

library(randomForest)
library(pROC)
source("./scripts/metagenome_analysis_functions.R")

#####################
### SET VARIABLES ###
#####################
cutoff  <- 0.2
outdir  <- paste0("./stats_out/ms_review/randomForest/")
if( !dir.exists(outdir) ){ dir.create(outdir,recursive = TRUE) }

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
kegg_abunds <- read.table(file = "./data/abundance_function_kegg.txt",
                     header = T,
                     stringsAsFactors = F,
                     sep = "\t")
kegg_abunds <- reduce_abunds(kegg_abunds,disease_to_sample)

m2_abunds <- read.table(file="./data/abundance_taxa_metaphlan2.txt",
                        header = T,
                        stringsAsFactors = F,
                        sep = "\t")
m2_abunds <- reduce_abunds(m2_abunds,disease_to_sample)

#make sure kegg and m2 have same samples
if( length(which(!(kegg_abunds[,"sample_id"] %in% m2_abunds[,"sample_id"]))) >= 1 ){
  ksamps <- kegg_abunds[,"sample_id"]
  missing <- unique(ksamps[which(!(kegg_abunds[,"sample_id"] %in% m2_abunds[,"sample_id"]))])
  kegg_abunds <- kegg_abunds[which(!(kegg_abunds[,"sample_id"] %in% missing)),]
}
if( length(which(!(m2_abunds[,"sample_id"] %in% kegg_abunds[,"sample_id"]))) >= 1 ){
  msamps <- m2_abunds[,"sample_id"]
  missing <- unique(msamps[which(!(m2_abunds[,"sample_id"] %in% kegg_abunds[,"sample_id"]))])
  m2_abunds <- m2_abunds[which(!(m2_abunds[,"sample_id"] %in% missing)),]
}

#remove repeat samples
if( length(which(!(disease_to_sample[,"sample_id"] %in% kegg_abunds[,"sample_id"]))) >= 1 ){
  disease_to_sample <- reduce_mapping(kegg_abunds,disease_to_sample)
}
mod_disease_to_sample <- modify_mapping(disease_to_sample)

kegg_results_wcorrection <- read.table("./stats_out/manuscript_figures/kegg_0.2/kegg_output_wcorrection.txt",header=T,stringsAsFactors = F,sep="\t")
kegg_shared_markers <- read.table("./stats_out/manuscript_figures/kegg_0.2/kegg_shared_markers.txt",header=T,stringsAsFactors = F,sep="\t",quote = "")

phenotypes <- c("arthritis","carcinoma","cirrhosis","crohns",
                "obesity","t2d","ulcerative_colitis")

#see what percent of each disease's indicator modules are shared with other diseases
table <- NULL
for(pheno in phenotypes){
  pheno_results <- kegg_results_wcorrection[which(kegg_results_wcorrection$phenotype == pheno),]
  pheno_sig_module <- pheno_results[which(pheno_results$level == "module" & pheno_results$fdr < 0.2),]
  
  sub_shared_markers <- kegg_shared_markers[which(kegg_shared_markers$name %in% pheno_sig_module$name),]
  
  nrow(sub_shared_markers[which(sub_shared_markers$phenotype_count > 1),])

  info <- c(pheno,nrow(pheno_sig_module),nrow(sub_shared_markers[which(sub_shared_markers$phenotype_count > 1),]))
    
  table <- rbind(table,info)  
}
colnames(table) <- c("disease","marker_count","shared_marker_count")
table <- as.data.frame(table)
table$marker_count <- as.numeric(as.character(table$marker_count))
table$shared_marker_count <- as.numeric(as.character(table$shared_marker_count))
table$pct_shared <- table$shared_marker_count/table$marker_count


colors <-  c(rgb(128,0,128,maxColorValue = 255),rgb(128,0,0,maxColorValue = 255),rgb(192,104,22,maxColorValue = 255),
             rgb(226,168,43,maxColorValue = 255),rgb(99,130,31,maxColorValue = 255),rgb(0,124,123,maxColorValue = 255),
             rgb(126,83,0,maxColorValue = 255))
names(colors) <- phenotypes

#####################
### RANDOM FOREST ###
#####################

### METADATA ONLY
outpath <- paste0(outdir,"metadata/")
if(!dir.exists(outpath)){dir.create(outpath,recursive = T)}

mod_att.map <- merge(samp_to_subj,att.map,by="subject_id")
 
metadata.df <- merge(mod_disease_to_sample,mod_att.map[,c("subject_id","sample_id","age","bmi","sex","country")],by="sample_id")
# if excluding: 
#exclude_samps <- c(disease_to_sample[which(disease_to_sample$status == "advanced adenoma"),"sample_id"],
#                   disease_to_sample[which(disease_to_sample$phenotype == "arthritis"),"sample_id"])
#metadata.df <- metadata.df[which(!(metadata.df$sample_id %in% exclude_samps)),]
#exclude old & young
#exclude_samps <- c(mod_att.map[which(mod_att.map$age > 70),"sample_id"],mod_att.map[which(mod_att.map$age < 30),"sample_id"])
#metadata.df <- metadata.df[which(!(metadata.df$sample_id %in% exclude_samps)),]
#
metadata.df <- metadata.df[,c("sample_id","study","status","age","bmi","sex","country")]
#metadata.df <- metadata.df[,c("sample_id","study","status","bmi","sex","country")]
metadata.df <- unique(metadata.df)
rownames(metadata.df) <- metadata.df$sample_id
metadata.df <- metadata.df[,-1]
metadata.df$status <- as.factor(metadata.df$status)
metadata.df <- na.omit(metadata.df)
metadata.df$study <- as.factor(metadata.df$study)
metadata.df <- metadata.df[,which(!colnames(metadata.df) == "bmi")]

metadata.rf <- randomForest(status ~ ., data=metadata.df,importance=TRUE,ntree=1000)
metadata.roc <- roc(metadata.df$status,metadata.rf$votes[,2])
plot(metadata.roc,xlim=c(1,0),ylim=c(0,1),xaxs="i",xaxs="i",lty="twodash",main="RF with Metadata")

hist(metadata.df[which(metadata.df$status == 0),"age"],main="Control AGE",xlab = "AGE",col = "lightblue")
hist(metadata.df[which(metadata.df$status == 1),"age"],main="Case AGE",xlab = "AGE",col = "firebrick1")
#hist(metadata.df[which(metadata.df$status == 0),"bmi"],main="Control BMI",xlab = "BMI",col = "lightblue")
#hist(metadata.df[which(metadata.df$status == 1),"bmi"],main="Case BMI",xlab = "BMI",col = "firebrick1")
table(metadata.df[which(metadata.df$status == 0),"study"])
table(metadata.df[which(metadata.df$status == 1),"study"])
table(metadata.df[which(metadata.df$status == 0),"country"])
table(metadata.df[which(metadata.df$status == 1),"country"])
table(paste(metadata.df$study,metadata.df$country,metadata.df$status,sep="."))
table(metadata.df[which(metadata.df$status == 0),"sex"])
table(metadata.df[which(metadata.df$status == 1),"sex"])

#plot importance
imp.table <- importance(metadata.rf)
top_acc <- names(sort(imp.table[,"MeanDecreaseAccuracy"],decreasing = T))
top_gini <- names(sort(imp.table[,"MeanDecreaseGini"],decreasing = T))

### ACCURACY
acc.df <- imp.table[top_acc,]
acc.df <- as.data.frame(cbind(name=row.names(acc.df),acc.df))
acc.df$MeanDecreaseAccuracy <- as.numeric(as.character(acc.df$MeanDecreaseAccuracy))

acc_plot1 <- ggplot(acc.df,aes(x=MeanDecreaseAccuracy,y=factor(top_acc,levels=rev(top_acc)))) + geom_point(shape=19,size=3) + 
  ylab("") + geom_hline(yintercept = c(factor(top_acc,levels=rev(top_acc))),linetype="dashed",colour="grey") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + ggtitle("")

### GINI
gini.df <- imp.table[top_gini,]
gini.df <- as.data.frame(cbind(name=row.names(gini.df),gini.df))
gini.df$MeanDecreaseGini <- as.numeric(as.character(gini.df$MeanDecreaseGini))

gini_plot1 <- ggplot(gini.df,aes(x=MeanDecreaseGini,y=factor(top_gini,levels=rev(top_gini)))) + geom_point(shape=19,size=3) + 
  ylab("") + geom_hline(yintercept = c(factor(top_gini,levels=rev(top_gini))),linetype="dashed",colour="grey") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


plot <- plot_grid(acc_plot1,gini_plot1,nrow = 1,align = "h")
final_plot <- ggdraw(plot) + draw_figure_label("Disease vs non-diseased",position = "top",fontface = "bold")

print(final_plot)

#BY DISEASE
meta.list <- list()
colors <- c("red","orange","yellow","green","blue","purple","brown")
colors <-  c(rgb(128,0,128,maxColorValue = 255),rgb(128,0,0,maxColorValue = 255),rgb(192,104,22,maxColorValue = 255),
             rgb(226,168,43,maxColorValue = 255),rgb(99,130,31,maxColorValue = 255),rgb(0,124,123,maxColorValue = 255),
             rgb(126,83,0,maxColorValue = 255))
names(colors) <- phenotypes
table <- data.frame(matrix(nrow = length(phenotypes),ncol=3,dimnames = list(phenotypes,c("level","disease","AUC"))))
for( pheno in phenotypes){
  pheno_samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == pheno),"sample_id"]
  sub_metadata.df <- na.omit(metadata.df[pheno_samps,])
  pheno_metadata.rf <- randomForest(status ~ .,data=sub_metadata.df,importance=TRUE,ntree=5000)
  print(pheno)
  print(pheno_metadata.rf)
  pheno.roc <- roc(sub_metadata.df$status,pheno_metadata.rf$votes[,2])
  if(pheno == "arthritis"){
    plot(pheno.roc,col=colors[pheno],xlim=c(1,0),ylim=c(0,1),xaxs="i",xaxs="i",lty="twodash",main="RF with Metadata",lwd=3)
  }else{
    plot(pheno.roc,add=T,col=colors[pheno],lty="twodash",lwd=3)
  }
  table[pheno,] <- c("metadata",pheno,pheno.roc$auc[1])
  meta.list[[pheno]] <- pheno.roc
  
  # hist(sub_metadata.df[which(sub_metadata.df$status == 0),"age"],main="Control AGE",xlab = "AGE",col = "lightblue")
  # hist(sub_metadata.df[which(sub_metadata.df$status == 1),"age"],main="Case AGE",xlab = "AGE",col = "firebrick1")
  # mean(sub_metadata.df[which(sub_metadata.df$status == 0),"age"])
  # median(sub_metadata.df[which(sub_metadata.df$status == 0),"age"])
  # mean(sub_metadata.df[which(sub_metadata.df$status == 1),"age"])
  # median(sub_metadata.df[which(sub_metadata.df$status == 1),"age"])
  # #hist(sub_metadata.df[which(sub_metadata.df$status == 0),"bmi"],main="Control BMI",xlab = "BMI",col = "lightblue")
  # #hist(sub_metadata.df[which(sub_metadata.df$status == 1),"bmi"],main="Case BMI",xlab = "BMI",col = "firebrick1")
  # table(sub_metadata.df[which(sub_metadata.df$status == 0),"study"])
  # table(sub_metadata.df[which(sub_metadata.df$status == 1),"study"])
  # table(sub_metadata.df[which(sub_metadata.df$status == 0),"country"])
  # table(sub_metadata.df[which(sub_metadata.df$status == 1),"country"])
  # table(paste(sub_metadata.df$study,sub_metadata.df$country,sub_metadata.df$status,sep="."))
  # table(sub_metadata.df[which(sub_metadata.df$status == 0),"sex"])
  # table(sub_metadata.df[which(sub_metadata.df$status == 1),"sex"])
  
  
  #plot importance
  imp.table <- importance(pheno_metadata.rf)
  top_acc <- names(sort(imp.table[,"MeanDecreaseAccuracy"],decreasing = T))
  top_gini <- names(sort(imp.table[,"MeanDecreaseGini"],decreasing = T))
  
  ### ACCURACY
  acc.df <- imp.table[top_acc,]
  acc.df <- as.data.frame(cbind(name=row.names(acc.df),acc.df))
  acc.df$MeanDecreaseAccuracy <- as.numeric(as.character(acc.df$MeanDecreaseAccuracy))
  
  acc_plot1 <- ggplot(acc.df,aes(x=MeanDecreaseAccuracy,y=factor(top_acc,levels=rev(top_acc)))) + geom_point(shape=19,size=3) + 
    ylab("") + geom_hline(yintercept = c(factor(top_acc,levels=rev(top_acc))),linetype="dashed",colour="grey") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + ggtitle("")
  
  ### GINI
  gini.df <- imp.table[top_gini,]
  gini.df <- as.data.frame(cbind(name=row.names(gini.df),gini.df))
  gini.df$MeanDecreaseGini <- as.numeric(as.character(gini.df$MeanDecreaseGini))
  
  gini_plot1 <- ggplot(gini.df,aes(x=MeanDecreaseGini,y=factor(top_gini,levels=rev(top_gini)))) + geom_point(shape=19,size=3) + 
    ylab("") + geom_hline(yintercept = c(factor(top_gini,levels=rev(top_gini))),linetype="dashed",colour="grey") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  
  plot <- plot_grid(acc_plot1,gini_plot1,nrow = 1,align = "h")
  final_plot <- ggdraw(plot) + draw_figure_label(paste0(get_p_name(pheno)),position = "top",fontface = "bold")

  #pdf(paste0(outpath,pheno,"_importance.pdf"))
  #print(final_plot)
  #dev.off()
  
}
legend(0.4,0.4,phenotypes,
       col=colors,lty="longdash",lwd=5)

### METADATA + MODULE ABUNDANCE 

mod_metadata.df <- data.frame(sample_id=row.names(metadata.df),metadata.df)

module_abunds <- kegg_abunds[which(kegg_abunds$level == "module"),]
module.df <- dcast(module_abunds,sample_id~name,value.var = "copy_number")
module.df <- merge(mod_disease_to_sample[,c("sample_id","status")],
                   module.df,
                   by="sample_id")
# if excluding: 
#exclude_samps <- c(disease_to_sample[which(disease_to_sample$status == "advanced adenoma"),"sample_id"],
#                   disease_to_sample[which(disease_to_sample$phenotype == "arthritis"),"sample_id"])
#module.df <- module.df[which(!(module.df$sample_id %in% exclude_samps)),]
#
module.df <- unique(module.df)
rownames(module.df) <- module.df$sample_id
module.df$status <- as.factor(module.df$status)

data.df <- merge(mod_metadata.df,module.df,by=c("sample_id","status"))
row.names(data.df) <- data.df$sample_id
data.df <- data.df[,-1]

data.rf <- randomForest(status ~ ., data=data.df,importance=TRUE,ntree=1000)

data.roc <- roc(data.df$status,data.rf$votes[,2])
plot(data.roc,xlim=c(1,0),ylim=c(0,1),xaxs="i",xaxs="i",lty="dotted",main="RF with metadata and module abundance",col="blue")

#plot importance
imp.table <- importance(data.rf)
top_acc <- names(sort(imp.table[1:30,"MeanDecreaseAccuracy"],decreasing = T))
top_gini <- names(sort(imp.table[1:30,"MeanDecreaseGini"],decreasing = T))

### ACCURACY
acc.df <- imp.table[top_acc,]
acc.df <- as.data.frame(cbind(name=row.names(acc.df),acc.df))
acc.df$MeanDecreaseAccuracy <- as.numeric(as.character(acc.df$MeanDecreaseAccuracy))

acc_plot1 <- ggplot(acc.df,aes(x=MeanDecreaseAccuracy,y=factor(top_acc,levels=rev(top_acc)))) + geom_point(shape=19,size=3) + 
  ylab("") + geom_hline(yintercept = c(factor(top_acc,levels=rev(top_acc))),linetype="dashed",colour="grey") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + ggtitle("")

### GINI
gini.df <- imp.table[top_gini,]
gini.df <- as.data.frame(cbind(name=row.names(gini.df),gini.df))
gini.df$MeanDecreaseGini <- as.numeric(as.character(gini.df$MeanDecreaseGini))

gini_plot1 <- ggplot(gini.df,aes(x=MeanDecreaseGini,y=factor(top_gini,levels=rev(top_gini)))) + geom_point(shape=19,size=3) + 
  ylab("") + geom_hline(yintercept = c(factor(top_gini,levels=rev(top_gini))),linetype="dashed",colour="grey") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


plot <- plot_grid(acc_plot1,gini_plot1,nrow = 1,align = "h")
final_plot <- ggdraw(plot) + draw_figure_label("Disease vs non-diseased",position = "top",fontface = "bold")

print(final_plot)

### MODULE ABUNDANCE 
module_abunds <- kegg_abunds[which(kegg_abunds$level == "module"),]
module.df <- dcast(module_abunds,sample_id~name,value.var = "copy_number")
module.df <- merge(mod_disease_to_sample[,c("sample_id","status")],
                   module.df,
                   by="sample_id")
# if excluding: 
#exclude_samps <- c(disease_to_sample[which(disease_to_sample$status == "advanced adenoma"),"sample_id"],
#                   disease_to_sample[which(disease_to_sample$phenotype == "arthritis"),"sample_id"])
#module.df <- module.df[which(!(module.df$sample_id %in% exclude_samps)),]
#
module.df <- unique(module.df)
rownames(module.df) <- module.df$sample_id
module.df <- module.df[,-1]
module.df$status <- as.factor(module.df$status)

module.rf <- randomForest(status ~ .,data=module.df,importance=T,ntree=1000)

module.roc <- roc(module.df$status,module.rf$votes[,2])

plot(module.roc,xlim=c(1,0),ylim=c(0,1),xaxs="i",xaxs="i",lty="solid",main="RF with module abundance",col="red")
#plot(module.roc,xlim=c(1,0),ylim=c(0,1),xaxs="i",xaxs="i",lty="solid",main="RF with module abundance",col="orange")

#plot importance
imp.table <- importance(module.rf)
top_acc <- names(sort(imp.table[1:30,"MeanDecreaseAccuracy"],decreasing = T))
top_gini <- names(sort(imp.table[1:30,"MeanDecreaseGini"],decreasing = T))

### ACCURACY
acc.df <- imp.table[top_acc,]
acc.df <- as.data.frame(cbind(name=row.names(acc.df),acc.df))
acc.df$MeanDecreaseAccuracy <- as.numeric(as.character(acc.df$MeanDecreaseAccuracy))

acc_plot1 <- ggplot(acc.df,aes(x=MeanDecreaseAccuracy,y=factor(top_acc,levels=rev(top_acc)))) + geom_point(shape=19,size=3) + 
  ylab("") + geom_hline(yintercept = c(factor(top_acc,levels=rev(top_acc))),linetype="dashed",colour="grey") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + ggtitle("")

### GINI
gini.df <- imp.table[top_gini,]
gini.df <- as.data.frame(cbind(name=row.names(gini.df),gini.df))
gini.df$MeanDecreaseGini <- as.numeric(as.character(gini.df$MeanDecreaseGini))

gini_plot1 <- ggplot(gini.df,aes(x=MeanDecreaseGini,y=factor(top_gini,levels=rev(top_gini)))) + geom_point(shape=19,size=3) + 
  ylab("") + geom_hline(yintercept = c(factor(top_gini,levels=rev(top_gini))),linetype="dashed",colour="grey") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


plot <- plot_grid(acc_plot1,gini_plot1,nrow = 1,align = "h")
final_plot <- ggdraw(plot) + draw_figure_label("Disease vs non-diseased",position = "top",fontface = "bold")

print(final_plot)

### PARAMETERIZE MODULE ABUND RF
#set number of parameters
test_n <- c(50,100,200,300,400,500)

plot(module.roc,xlim=c(1,0),ylim=c(0,1),xaxs="i",xaxs="i",lty="solid",col="black")
cols <- c("red","orange","green","blue","purple","grey")
names(cols) <- test_n
for(n in test_n){
  imp.table <- importance(module.rf)
  top_acc <- names(sort(imp.table[1:n,"MeanDecreaseAccuracy"],decreasing = T))
  top_gini <- names(sort(imp.table[1:n,"MeanDecreaseGini"],decreasing = T))
  
  n_module.df <- module.df[,c("status",top_acc)] 
  
  n_module.rf <- randomForest(status ~ .,data=n_module.df,importance=T,ntree=1000)
  print(n)
  print(n_module.rf)
  
  n_module.roc <- roc(n_module.df$status,n_module.rf$votes[,2])
  print(n_module.roc$auc)
  
  plot(n_module.roc,xlim=c(1,0),ylim=c(0,1),xaxs="i",xaxs="i",lty="solid",add=T,col=cols[as.character(n)])
}
legend(0.4,0.38,c(test_n,521),col=c(cols,"black"),lwd=2)

#all three plots
plot(metadata.roc,xlim=c(1,0),ylim=c(0,1),xaxs="i",xaxs="i",lty="twodash")
plot(data.roc,add=T,lty="dotted",col="blue")
plot(module.roc,add=T,lty="solid",col="red")

legend(0.4,0.38,c("metadata","module","both"),lty=c("twodash","solid","dotted"),col=c("black","red","blue"),lwd=2)

### KEGG
outpath <- paste0(outdir,"kegg/")
if(!dir.exists(outpath)){dir.create(outpath,recursive = T)}

datatype <- "copy_number"
module_abunds <- kegg_abunds[which(kegg_abunds$level == "module"),]
module.df <- dcast(module_abunds,sample_id~name,value.var = datatype)
module.df <- merge(mod_disease_to_sample[,c("sample_id","status")],
                   module.df,
                   by="sample_id")
module.df <- unique(module.df)
rownames(module.df) <- module.df$sample_id
module.df <- module.df[,-1]
module.df$status <- as.factor(module.df$status)

module.rf <- randomForest(status ~ .,data=module.df,importance=TRUE,ntree=5000)
#module.plot <- plot_importance(module.rf,shared_markers,"All Modules")
#print(module.plot)

###by disease
p.list <- list()

#colors <- c("red","orange","yellow","green","blue","purple","brown")
# colors <-  c(rgb(128,0,128,maxColorValue = 255),rgb(128,0,0,maxColorValue = 255),rgb(192,104,22,maxColorValue = 255),
#              rgb(226,168,43,maxColorValue = 255),rgb(99,130,31,maxColorValue = 255),rgb(0,124,123,maxColorValue = 255),
#              rgb(126,83,0,maxColorValue = 255))
# names(colors) <- phenotypes
# 
# mod_phenotypes <- c("Rheumatoid arthritis","Colorecta carcinoma","Carcinoma + adenoma","Liver cirrhosis","Crohn's disease",
#                     "Obesity","Type II diabetes","Ulcerative colitis")
# mod_colors <- c(rgb(128,0,128,maxColorValue = 255),rgb(128,0,0,maxColorValue = 255),rgb(128,0,0,maxColorValue = 255),rgb(192,104,22,maxColorValue = 255),
#                 rgb(226,168,43,maxColorValue = 255),rgb(99,130,31,maxColorValue = 255),rgb(0,124,123,maxColorValue = 255),
#                 rgb(126,83,0,maxColorValue = 255))
# names(mod_colors) <- mod_phenotypes

mod_phenotypes <- c("arthritis","carcinoma","carcinoma2","cirrhosis","crohns","obesity","t2d","ulcerative_colitis")
mod_colors <- c(rgb(128,0,128,maxColorValue = 255),rgb(128,0,0,maxColorValue = 255),rgb(128,0,0,maxColorValue = 255),rgb(192,104,22,maxColorValue = 255),
                rgb(226,168,43,maxColorValue = 255),rgb(99,130,31,maxColorValue = 255),rgb(0,124,123,maxColorValue = 255),
                rgb(126,83,0,maxColorValue = 255))
names(mod_colors) <- mod_phenotypes

pheno_lty <- c(1,1,3,1,1,1,1,1)
names(pheno_lty) <- mod_phenotypes

#pdf(paste0(outpath,"module_by_disease_roc.pdf"))
kegg_table <- data.frame(matrix(nrow = length(phenotypes),ncol=3,dimnames = list(phenotypes,c("level","disease","AUC"))))
for(pheno in mod_phenotypes){
  if(pheno == "carcinoma2"){
    pheno_samps <- disease_to_sample[which(disease_to_sample$phenotype == "carcinoma" & disease_to_sample$status %in% c("controls","carcinoma")),"sample_id"]
  }else{
    pheno_samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == pheno),"sample_id"]
  }

  sub_module.df <- module.df[pheno_samps,]
  pheno_module.rf <- randomForest(status ~ .,data=sub_module.df,importance=TRUE,ntree=5000)
  print(pheno)
  print(pheno_module.rf)
  pheno.roc <- roc(sub_module.df$status,pheno_module.rf$votes[,2])
  if(pheno == "arthritis"){
    plot(pheno.roc,col=mod_colors[pheno],lwd=3,lty=pheno_lty[pheno])
  }else{
    plot(pheno.roc,add=T,col=mod_colors[pheno],lwd=3,lty=pheno_lty[pheno])
  }
  kegg_table[pheno,] <- c("module",pheno,pheno.roc$auc[1])
  p.list[[pheno]] <- pheno.roc
}

mod_phenotypes <- c("Rheumatoid arthritis","Carcinoma + adenoma","Carcinoma only","Liver cirrhosis","Crohn's disease",
                    "Obesity","Type II diabetes","Ulcerative colitis")
mod_colors <- c(rgb(128,0,128,maxColorValue = 255),rgb(128,0,0,maxColorValue = 255),rgb(128,0,0,maxColorValue = 255),rgb(192,104,22,maxColorValue = 255),
                rgb(226,168,43,maxColorValue = 255),rgb(99,130,31,maxColorValue = 255),rgb(0,124,123,maxColorValue = 255),
                rgb(126,83,0,maxColorValue = 255))
names(mod_colors) <- c("arthritis","carcinoma","carcinoma2","cirrhosis","crohns",
                       "obesity","T2D","ulcerative_colitis")

legend(0.45,0.4,mod_phenotypes,lwd=3,bty="n",
       col=mod_colors,lty=c(1,2,1,1,1,1,1,1))
dev.off()

mod_colors <- c(rgb(128,0,128,maxColorValue = 255),rgb(128,0,0,maxColorValue = 255),"grey30",rgb(192,104,22,maxColorValue = 255),
                rgb(226,168,43,maxColorValue = 255),rgb(99,130,31,maxColorValue = 255),rgb(0,124,123,maxColorValue = 255),
                rgb(126,83,0,maxColorValue = 255))
names(mod_colors) <- c("arthritis","carcinoma","carcinoma2","cirrhosis","crohns",
                       "obesity","T2D","ulcerative_colitis")

pdf(paste0(outpath,"module_by_disease_roc_ggplot.pdf"))
ggroc(data = p.list,size=1.2) + geom_abline(intercept = 1,colour="grey") + coord_equal() + scale_colour_manual(values=alpha(mod_colors,0.7),labels=mod_phenotypes) +
  xlab("Specificity") + ylab("Sensitivity")  
  #scale_linetype_manual(values=pheno_lty,labels=mod_phenotypes)
dev.off()

colors <-  c(rgb(128,0,128,maxColorValue = 255),rgb(128,0,0,maxColorValue = 255),rgb(192,104,22,maxColorValue = 255),
             rgb(226,168,43,maxColorValue = 255),rgb(99,130,31,maxColorValue = 255),rgb(0,124,123,maxColorValue = 255),
             rgb(126,83,0,maxColorValue = 255))
names(p.list) <- c("Rheumatoid arthritis","Colorectal cancer","Liver cirrhosis","Crohn's disease",
                   "Obesity","Type II diabetes","Ulcerative colitis")
plot <- ggroc(p.list,size=1.2) + geom_abline(intercept = 1,colour="grey") + coord_equal() + scale_color_manual(values=alpha(colors,0.7)) +
  xlab("Specificity") + ylab("Sensitivity")


#COMPARE METADATA/MODULE BY DISEASE
table <- NULL
table2 <- NULL
for( pheno in phenotypes){
  pheno_samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == pheno),"sample_id"]
  #colorectal cancer - exclude adenoma
  #pheno_samps <- disease_to_sample[which(disease_to_sample$phenotype == "carcinoma" & disease_to_sample$status %in% c("controls","carcinoma")),"sample_id"]
  
  #metadata
  sub_metadata.df <- na.omit(metadata.df[pheno_samps,])
  pheno_metadata.rf <- randomForest(status ~ .,data=sub_metadata.df,importance=TRUE,ntree=1000)
  print(paste0("metadata ",pheno))
  print(pheno_metadata.rf)
  metadata_pheno.roc <- roc(sub_metadata.df$status,pheno_metadata.rf$votes[,2])
  #metadata - no age
  sub_metadata2.df <- na.omit(metadata.df[pheno_samps,c("study","status","sex","country")])
  #sub_metadata2.df <- na.omit(metadata.df[pheno_samps,c("status","sex","age")])
  pheno_metadata2.rf <- randomForest(status ~ .,data=sub_metadata2.df,importance=TRUE,ntree=1000)
  print(paste0("metadata2 ",pheno))
  print(pheno_metadata2.rf)
  metadata_pheno2.roc <- roc(sub_metadata2.df$status,pheno_metadata2.rf$votes[,2])
  #module
  sub_module.df <- module.df[pheno_samps,]
  pheno_module.rf <- randomForest(status ~ .,data=sub_module.df,importance=TRUE,ntree=1000)
  print(paste0("module ",pheno))
  print(pheno_module.rf)
  module_pheno.roc <- roc(sub_module.df$status,pheno_module.rf$votes[,2])
  #metadata plus module
  #pheno.df <- data.frame()
  
  #pdf(paste0("./stats_out/ms_review/randomForest/metadata-module/",pheno,"_compare.pdf"))
  plot(metadata_pheno.roc,col="black",xlim=c(1,0),ylim=c(0,1),xaxs="i",xaxs="i",
       lty=1,main=paste0(get_p_name(pheno),"\nMetadata vs Module abundance"),lwd=3)
  plot(metadata_pheno2.roc,col="grey",add=T,lty=1,lwd=3)
  plot(module_pheno.roc,add=T,col=colors[pheno],lty=1,lwd=3)
  #plot(module_pheno.roc,add=T,col=colors[pheno],lty=1,lwd=3)
  legend(0.35,0.35,c("metadata","metadata-no age","module"),
         col=c("black","grey",colors[pheno]),lty=c(1,1,1),lwd=3)
  #dev.off()
  
  table <- rbind(table,
                 c(pheno,"metadata",metadata_pheno.roc$auc,pheno_metadata.rf$err.rate[nrow(pheno_metadata.rf$err.rate),"OOB"]),
                 c(pheno,"metadata-no age",metadata_pheno2.roc$auc,pheno_metadata2.rf$err.rate[nrow(pheno_metadata2.rf$err.rate),"OOB"]),
                 c(pheno,"kegg_module",module_pheno.roc$auc,pheno_module.rf$err.rate[nrow(pheno_module.rf$err.rate),"OOB"]))
  print("metadata vs module roc test")
  print(roc.test(metadata_pheno.roc,module_pheno.roc)$p.value)
  table2 <- rbind(table2,c(pheno,"metadata-module",roc.test(metadata_pheno.roc,module_pheno.roc)$p.value))
}
colnames(table) <- c("disease","model","AUC","OOB")
colnames(table2) <- c("disease","test","pvalue")
write.table(table,"/Users/armourc/Desktop/table.txt",quote=F,sep="\t",row.names = F)
write.table(table2,"/Users/armourc/Desktop/table2.txt",quote=F,sep="\t",row.names = F)

### COLORECTAL CANCER
#all subjects
pheno_samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == "carcinoma"),"sample_id"]

sub_module.df <- module.df[pheno_samps,]
pheno_module.rf <- randomForest(status ~ .,data=sub_module.df,importance=TRUE,ntree=5000)
print(pheno)
print(pheno_module.rf)
pheno.roc <- roc(sub_module.df$status,pheno_module.rf$votes[,2])
#exclude advanced adenoma subjects
pheno_samps <- disease_to_sample[which(disease_to_sample$phenotype == "carcinoma" & disease_to_sample$status %in% c("controls","carcinoma")),"sample_id"]

sub_module.df <- module.df[pheno_samps,]
pheno_module2.rf <- randomForest(status ~ .,data=sub_module.df,importance=TRUE,ntree=5000)
print(pheno)
print(pheno_module2.rf)
pheno.roc2 <- roc(sub_module.df$status,pheno_module2.rf$votes[,2])
#plot
plot(pheno.roc,col=colors["carcinoma"] )
plot(pheno.roc2,add=T,col=colors["carcinoma"],lty="dashed")
legend(0.5,0.3,c("all subjects","without adenoma"),col=colors,lty=c(1,2))

roc.test(pheno.roc,pheno.roc2)

### METAPHLAN2
datatype <- "relative_abundance"

genus_abunds <- m2_abunds[which(m2_abunds$level == "genus"),]
genus.df <- dcast(genus_abunds,sample_id~name,value.var = datatype)
genus.df <- merge(mod_disease_to_sample[,c("sample_id","status")],
                   genus.df,
                   by="sample_id")
genus.df <- unique(genus.df)
rownames(genus.df) <- genus.df$sample_id
genus.df <- genus.df[,-1]
genus.df$status <- as.factor(genus.df$status)

genus.rf <- randomForest(status ~ .,data=genus.df,importance=TRUE,ntree=5000)
genus.plot <- plot_importance(genus.rf,shared_markers,"All Genera")
print(genus.plot)


#### COMPARE
kegg_level <- "module" # gene_family, module, or pathway
m2_level <- "species"

datatype <- "copy_number"
kegg_level_abunds <- kegg_abunds[which(kegg_abunds$level == kegg_level),]
kegg.df <- dcast(kegg_level_abunds,sample_id~name,value.var = datatype)
kegg.df <- merge(mod_disease_to_sample[,c("sample_id","status")],
                   kegg.df,
                   by="sample_id")
kegg.df <- unique(kegg.df)
rownames(kegg.df) <- kegg.df$sample_id
kegg.df <- kegg.df[,-1]
kegg.df$status <- as.factor(kegg.df$status)

datatype <- "relative_abundance"
m2_level_abunds <- m2_abunds[which(m2_abunds$level == m2_level),]
m2.df <- dcast(m2_level_abunds,sample_id~name,value.var = datatype)
m2.df <- merge(mod_disease_to_sample[,c("sample_id","status")],
                  m2.df,
                  by="sample_id")
m2.df <- unique(m2.df)
rownames(m2.df) <- m2.df$sample_id
m2.df <- m2.df[,-1]
m2.df$status <- as.factor(m2.df$status)

colors <- c("red","orange","yellow","green","blue","purple","brown")
colors <-  c(rgb(128,0,128,maxColorValue = 255),rgb(128,0,0,maxColorValue = 255),rgb(192,104,22,maxColorValue = 255),
             rgb(226,168,43,maxColorValue = 255),rgb(99,130,31,maxColorValue = 255),rgb(0,124,123,maxColorValue = 255),
             rgb(126,83,0,maxColorValue = 255))
names(colors) <- phenotypes
table <- NULL
phenotypes <- c("carcinoma","cirrhosis","crohns")
colors <- c("#F8766D","#00BA38","#619CFF")
names(colors) <- phenotypes
pie(rep(1,length(colors)),col = colors,labels = phenotypes)
for(pheno in phenotypes){
  pdf(paste0(outdir,"compare/",pheno,"_",kegg_level,"-",m2_level,"_roc.pdf"))
  pheno_samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == pheno),"sample_id"]
  #for colorectal cancer to exclude advanced adenoma subjects
  if(pheno == "carcinoma"){
    pheno_samps <- disease_to_sample[which(disease_to_sample$phenotype == pheno & disease_to_sample$status %in% c("controls","carcinoma")),"sample_id"]
  }
  
  #kegg
  sub_kegg.df <- kegg.df[pheno_samps,]
  pheno_kegg.rf <- randomForest(status ~ .,data=sub_kegg.df,importance=TRUE,ntree=10000)
  print(pheno)
  print(pheno_kegg.rf)
  pheno.kegg.roc <- roc(sub_kegg.df$status,pheno_kegg.rf$votes[,2])
  
  plot(pheno.kegg.roc,col=colors[pheno],lwd=3,axes=F,xlim=c(1,0),ylim=c(0,1),xaxs="i",xaxs="i",main=get_p_name(pheno))
  
  table <- rbind(table,c(pheno,"kegg",kegg_level,pheno.kegg.roc$auc[1],pheno_kegg.rf$err.rate[nrow(pheno_kegg.rf$err.rate),1]))
  
  #metaphlan
  sub_m2.df <- m2.df[pheno_samps,]
  pheno_m2.rf <- randomForest(status ~ .,data=sub_m2.df,importance=TRUE,ntree=10000)
  print(pheno)
  print(pheno_m2.rf)
  pheno.m2.roc <- roc(sub_m2.df$status,pheno_m2.rf$votes[,2])
  
  table <- rbind(table,c(pheno,"metaphlan2",m2_level,pheno.m2.roc$auc[1],pheno_m2.rf$err.rate[nrow(pheno_m2.rf$err.rate),1]))
  
  plot(pheno.m2.roc,col=colors[pheno],add=T,lty="dashed",lwd=3)
  axis(side = 1,lwd=2)
  axis(side = 2,lwd=2)
  legend(0.45,0.3,c(paste0("kegg-",kegg_level),paste0("metaphlan2-",m2_level)),col=colors[pheno],lty=c(1,2),lwd=3)
  dev.off()
  
  #ggroc(list(pheno.kegg.roc,pheno.m2.roc),size=1.2,lty=c("dashed")) + geom_abline(intercept = 1,colour="grey") + coord_equal()  
  #  scale_color_manual(values=alpha(colors,0.7)) +
  #  xlab("Specificity") + ylab("Sensitivity")
    
#    ggroc(list(s100b=rocobj, wfns=rocobj2), linetype=2)
 
  #   ggroc(list(pheno.kegg.roc,pheno.m2.roc), aes="linetype", color="red")
    
  test <- roc.test(pheno.kegg.roc,pheno.m2.roc)
  print(paste0("roc test pvalue: ",test$p.value))
}
colnames(table) <- c("disease","database","level","AUC","OOB")
write.table(table,paste0(outdir,"compare/",kegg_level,"-",m2_level,"_values.txt"),row.names = F,quote = F,sep="\t")

#### RESULTS BY STUDY
## T2D
outpath <- paste0("./stats_out/ms_review/t2d_by_study/")
if(!dir.exists(outpath)){dir.create(outpath)}

t2d_study_results <- read.table("./model_output/t2d_study_results.txt",header=F,stringsAsFactors = F,sep="\t")
#t2d_study_results <- read.table("./model_output/",header=F,stringsAsFactors = F,sep="\t")
colnames(t2d_study_results) <- c("type","database","level","name","pheno","study",'case_n',"ctrl_n","case_mean","ctrl_mean","change_in_case","pvalue")

t2d_detail_results <- read.table("./model_output/t2d_detail_results.txt",header=T,stringsAsFactors = F,sep = "\t")

t2d_study_results <- merge(t2d_study_results,t2d_detail_results[which(t2d_detail_results$parameter == "status"),],by=c("name","pheno","study"))

t2d_studies <- unique(t2d_study_results$study)

all_results <- read.table("./model_output/kegg_output.tab",header=F,stringsAsFactors = F,sep="\t")
colnames(all_results) <- c("type","database","level","name","phenotype","phenotype_level","studies","case_n","ctrl_n","case_mean","ctrl_mean","pvalue")
all_t2d_results <- all_results[which(all_results$phenotype == "t2d" & all_results$level == "module"),]

all_detail <- detail_results <- read.table("./model_output/kegg_output_detail.txt",header=F,stringsAsFactors = F,sep="\t")
colnames(all_detail) <- c("type","database","level","name","phenotype","component","estimate","std.err","t_value","p_value")
all_t2d_detail <- all_detail[which(all_detail$phenotype == "t2d" & all_detail$level == "module"),]

all_t2d_results <- merge(all_t2d_results,all_t2d_detail,by=c("type","database","level","name","phenotype"))

all_p0.05 <- all_t2d_results[which(all_t2d_results$pvalue < 0.05),"name"]

#table <- data.frame(module=unique(t2d_study_results$name))
plots <- list()
p.cors <- NULL
for(study in t2d_studies){
  sub_study_results <- t2d_study_results[which(t2d_study_results$study == study),]
  sub_study_results.ordered <- sub_study_results[order(sub_study_results$pvalue,decreasing = F),c("name","estimate","pvalue")]
  sub_study_results.ordered <- data.frame(sub_study_results.ordered,rank=seq(1,nrow(sub_study_results.ordered),1))
  colnames(sub_study_results.ordered) <- c("module","study_estimate","study_pvalue","study_rank")
  
  #sub_all_t2d_results <- all_t2d_results[,c("name","pvalue")]
  all_t2d_results.ordered <- all_t2d_results[order(all_t2d_results$pvalue,decreasing = F),c("name","estimate","pvalue")]
  all_t2d_results.ordered <- data.frame(all_t2d_results.ordered,rank=seq(1,nrow(all_t2d_results.ordered),1))
  colnames(all_t2d_results.ordered) <- c("module","full_estimate","full_pvalue","full_rank")
  
  table <- merge(sub_study_results.ordered,all_t2d_results.ordered,by="module")
  
  cor <- cor.test(table$full_estimate,table$study_estimate,method="spearman")
  rsq <- cor$estimate
  lbl <- paste0(expression(rho)," == ",round(rsq,digits=2))
  
  plot <- ggplot(table,aes(x=full_estimate,y=study_estimate)) + geom_point(alpha=0.7) +
    ggtitle(paste0(study," slope estimates")) 
    #annotate("text",x=-20,y=median(table$study_estimate),label = lbl, parse = TRUE,size=6,col="blue") 
  
  #pdf(paste0(outpath,"study_",study,"_estimate.pdf"))
  print(plot)
  #dev.off()
  #plots[[study]] <- plot
  

  cor <- cor.test(table$full_rank,table$study_rank,method="spearman")
  rsq <- cor$estimate
  lbl <- paste0(expression(rho)," == ",round(rsq,digits=2))
  
  plot <- ggplot(table,aes(x=full_rank,y=study_rank)) + geom_point(alpha=0.7) +
    ggtitle(paste0(study," pvalue ranking")) +
    annotate("text",x=40,y=500,label = lbl, parse = TRUE,size=6,col="blue") 
  
  pdf(paste0(outpath,"study_",study,"_rank.pdf"))
  print(plot)
  dev.off()
  
  cor <- cor.test(table$full_pvalue,table$study_pvalue,method="spearman")
  rsq <- cor$estimate
  lbl <- paste0(expression(rho)," == ",round(rsq,digits=2))
  
  study_p0.05 <- sub_study_results[which(sub_study_results$pvalue < 0.05),"name"]
  count <- length(study_p0.05)
  scount <- length(study_p0.05[which(study_p0.05 %in% all_p0.05)])
  p.cors <- rbind(p.cors,c(study,rsq,count,scount))
  
  plot <- ggplot(table,aes(x=full_pvalue,y=study_pvalue)) + geom_point(alpha=0.7) +
    ggtitle(study) + xlab("Merged Dataset P-value") + ylab("Study P-value") + 
    coord_equal()
    #annotate("text",x=0.10,y=0.90,label = lbl, parse = TRUE,size=6,col="blue") 
  
  #pdf(paste0(outpath,"study_",study,"_pvalues.pdf"))
  #print(plot)
  #dev.off()
  #plots[[study]] <- plot
}

pdf(paste0(outpath,"all_study_pvalues.pdf"),width=11)
plot_grid(plotlist = plots,ncol=3)
dev.off()



# mat <- as.matrix(all_t2d_detail[which(all_t2d_detail$component == "status"),c("name","estimate")])
# colnames(mat) <- c("name","full_estimate") 
# for(study in t2d_studies){
#   sub_study_results <- t2d_study_results[which(t2d_study_results$study == study),c("name","estimate")]
#   colnames(sub_study_results) <- c("name",paste0(study,"_estimate"))
#   mat <- merge(mat,sub_study_results,by="name")
# }
# 
# rownames(mat) <- mat$name
# mat <- as.matrix((mat[,-1]))
# mat <- mapply(mat, as.numeric)
# mat <- matrix(data=mat, ncol=5, nrow=521)
# 
# heatmap.2(mat,dendrogram = "none",trace = "none",col=c("red","blue"),srtCol = 0)


## OBESITY
outpath <- paste0("./stats_out/ms_review/ob_by_study/")
if(!dir.exists(outpath)){dir.create(outpath)}

ob_study_results <- read.table("./model_output/ob_study_results.txt",header=F,stringsAsFactors = F,sep="\t")
colnames(ob_study_results) <- c("type","database","level","name","disease","study",'case_n',"ctrl_n","case_mean","ctrl_mean","change_in_case","pvalue")

ob_studies <- unique(ob_study_results$study)

all_results <- read.table("./model_output/kegg_output.tab",header=F,stringsAsFactors = F,sep="\t")
colnames(all_results) <- c("type","database","level","name","phenotype","phenotype_level","studies","case_n","ctrl_n","case_mean","ctrl_mean","pvalue")
all_ob_results <- all_results[which(all_results$phenotype == "obesity" & all_results$level == "module"),]

all_p0.05 <- all_ob_results[which(all_ob_results$pvalue < 0.05),"name"]

plots <- list()
p.cors <- NULL
for(study in ob_studies){
  sub_study_results <- ob_study_results[which(ob_study_results$study == study),]
  sub_study_results.ordered <- sub_study_results[order(sub_study_results$pvalue,decreasing = F),c("name","pvalue")]
  sub_study_results.ordered <- data.frame(sub_study_results.ordered,rank=seq(1,nrow(sub_study_results.ordered),1))
  colnames(sub_study_results.ordered) <- c("module","study_pvalue","study_rank")
  
  all_ob_results.ordered <- all_ob_results[order(all_ob_results$pvalue,decreasing = F),c("name","pvalue")]
  all_ob_results.ordered <- data.frame(all_ob_results.ordered,rank=seq(1,nrow(all_ob_results.ordered),1))
  colnames(all_ob_results.ordered) <- c("module","full_pvalue","full_rank")
  
  table <- merge(sub_study_results.ordered,all_ob_results.ordered,by="module")
  if(nrow(na.omit(table)) == 0){
    next()
  }
  cor <- cor.test(table$full_rank,table$study_rank,method="spearman")
  rsq <- cor$estimate
  lbl <- paste0(expression(rho)," == ",round(rsq,digits=2))
  
  plot <- ggplot(table,aes(x=full_rank,y=study_rank)) + geom_point(alpha=0.7) +
    ggtitle(paste0(study," pvalue ranking")) + 
    annotate("text",x=55,y=515,label = lbl, parse = TRUE,size=6,col="blue") 
  
  pdf(paste0(outpath,"study_",study,"_rank.pdf"))
  print(plot)
  dev.off()
  
  cor <- cor.test(table$full_pvalue,table$study_pvalue,method="spearman")
  rsq <- cor$estimate
  lbl <- paste0(expression(rho)," == ",round(rsq,digits=2))
  
  plot <- ggplot(table,aes(x=full_pvalue,y=study_pvalue)) + geom_point(alpha=0.7) +
    ggtitle(study) + xlab("Merged Dataset P-value") + ylab("Study P-value") + 
    coord_equal()
    #annotate("text",x=0.50,y=0.90,label = lbl, parse = TRUE,size=6,col="blue") 

  study_p0.05 <- sub_study_results[which(sub_study_results$pvalue < 0.05),"name"]
  count <- length(study_p0.05)
  scount <- length(study_p0.05[which(study_p0.05 %in% all_p0.05)])
  p.cors <- rbind(p.cors,c(study,rsq,count,scount))
  
  plots[[study]] <- plot
  
  #pdf(paste0(outpath,"study_",study,"_pvalues.pdf"))
  #print(plot)
  #dev.off()
}
pdf(paste0(outpath,"all_study_pvalues.pdf"),width=13)
plot_grid(plotlist = plots,ncol=4)
dev.off()

### RANDOM FOREST BY STUDY

# T2D
outpath <- paste0(outdir,"kegg/by_study/")

t2d_samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == "t2d"),"sample_id"]
t2d_subjs <- samp_to_subj[which(samp_to_subj$sample_id %in% t2d_samps),"subject_id"]
t2d_metadata <- att.map[which(att.map$subject_id %in% t2d_subjs),c("subject_id","age","bmi","t2d","country","sex")]
t2d_metadata <- merge(samp_to_subj,t2d_metadata,by="subject_id")
t2d_metadata <- merge(mod_disease_to_sample[which(mod_disease_to_sample$phenotype == "t2d"),],t2d_metadata,by="sample_id")
mod_t2d_metadata <- t2d_metadata
mod_t2d_metadata$study <- ifelse(mod_t2d_metadata$study %in% c("IGC","Richness"),"MetaHIT",mod_t2d_metadata$study)

t2d_studies <- unique(mod_t2d_metadata$study)
colors <- c("#F8766D","#00BA38","black","#619CFF")
names(colors) <- t2d_studies

#kegg
kegg_level <- "module"
datatype <- "copy_number"
kegg_level_abunds <- kegg_abunds[which(kegg_abunds$level == kegg_level),]
kegg.df <- dcast(kegg_level_abunds,sample_id~name,value.var = datatype)
kegg.df <- merge(mod_disease_to_sample[,c("sample_id","status")],
                 kegg.df,
                 by="sample_id")
kegg.df <- unique(kegg.df)
#with metadata
#kegg.df <- merge(t2d_metadata[,c("sample_id","status","study","age","country","sex")],kegg.df,by=c("sample_id","status"))
#kegg.df <- na.omit(kegg.df)
#kegg.df$study <- as.factor(kegg.df$study)
#
rownames(kegg.df) <- kegg.df$sample_id
kegg.df <- kegg.df[,-1]
kegg.df$status <- as.factor(kegg.df$status)

t2d_kegg.df <- kegg.df[t2d_samps,]

study <- "T2D"
study_samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == "t2d" & mod_disease_to_sample$study == study),"sample_id"]
study_kegg.df <- t2d_kegg.df[study_samps,]

#study_kegg.rf <- randomForest(status ~ .,data=t2d_kegg.df,importance=TRUE,ntree=5000,subset = study_samps)
study_kegg.rf <- randomForest(status ~ .,data=study_kegg.df,importance=TRUE,ntree=10000)

study.kegg.roc <- roc(t2d_kegg.df[study_samps,"status"],study_kegg.rf$votes[study_samps,2])

#test calculating OOB
#test <- as.numeric(colnames(study_kegg.rf$votes)[1:2][apply(study_kegg.rf$votes[study_samps,1:2],1,which.max)])
#length(which(t2d_kegg.df[study_samps,"status"] != test))/length(t2d_kegg.df[study_samps,"status"])

pdf(paste0(outdir,"t2d_by_study/train_T2D_roc.pdf"))
plot(study.kegg.roc,col=colors[study])
table <- c(study,"train",mean( predict(study_kegg.rf) != t2d_kegg.df[study_samps,"status"]),study.kegg.roc$auc)

for(study in t2d_studies){
  if(study == "T2D"){
    next()
  }
  if(study == "MetaHIT"){
    #next()
    study_samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == "t2d" & mod_disease_to_sample$study %in% c("IGC","Richness")),"sample_id"]
  }else {
    study_samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == "t2d" & mod_disease_to_sample$study == study),"sample_id"]  
  }
  
  #pred <- as.data.frame(predict(study_kegg.rf,t2d_kegg.df[study_samps,],type="response"))
  pred <- predict(study_kegg.rf,t2d_kegg.df[study_samps,], type="prob", predict.all=T, norm.votes = T)
  
  #pred$predict <- names(pred)[1:2][apply(pred[,1:2], 1, which.max)]
  #pred$observed <- t2d_kegg.df[study_samps,"status"]
  
  #roc <- roc(as.numeric(as.character(pred$observed)), as.numeric(as.character(pred$predict)))
  roc <- roc( t2d_kegg.df[study_samps,"status"], pred$aggregate[study_samps,2])
  
  #oob <- mean( predict(study_kegg.rf,t2d_kegg.df[study_samps,]) != t2d_kegg.df[study_samps,"status"])

  status_pred <- as.numeric(colnames(pred$aggregate)[1:2][apply(pred$aggregate[study_samps,1:2],1,which.max)])
  oob <- length(which(t2d_kegg.df[study_samps,"status"] != status_pred))/length(t2d_kegg.df[study_samps,"status"])    
  
  plot(roc,add=T,col=colors[study])
  
  table <- rbind(table,c(study,"test",oob,roc$auc))
}
legend(0.35,0.38,t2d_studies,col=colors,lty="solid",lwd=2)
dev.off()
colnames(table) <- c("study","datatype","OOB","AUC")
write.table(table,paste0(outdir,"t2d_by_study/values.txt"),quote=F,sep="\t",row.names = F)

#all data together
plot(study.kegg.roc,col="black",xlim=c(1,0),ylim=c(0,1),xaxs="i",xaxs="i")
study_samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == "t2d" & mod_disease_to_sample$study == "T2D"),"sample_id"]
table <- c("T2D","train",mean( predict(study_kegg.rf) != t2d_kegg.df[study_samps,"status"]),study.kegg.roc$auc)

other_samps <- rownames(t2d_kegg.df)[which(!(rownames(t2d_kegg.df) %in% study_samps))]

pred <- predict(study_kegg.rf,t2d_kegg.df[other_samps,], type="prob", predict.all=T, norm.votes = T)

roc <- roc( t2d_kegg.df[other_samps,"status"], pred$aggregate[other_samps,2])

status_pred <- as.numeric(colnames(pred$aggregate)[1:2][apply(pred$aggregate[other_samps,1:2],1,which.max)])
oob <- length(which(t2d_kegg.df[other_samps,"status"] != status_pred))/length(t2d_kegg.df[other_samps,"status"])

plot(roc,add=T,col="black",lty="dashed")
legend(0.35,0.38,c("train","test"),lty=c("solid","dashed"),lwd=2)

table <- rbind(table,c("all_other_samps","test",oob,roc$auc))


# leave one out
t2d_studies <- unique(mod_t2d_metadata$study)
colors <- c("#F8766D","#00BA38","purple","#619CFF")
names(colors) <- t2d_studies
table <- NULL
for(study in t2d_studies){
  if(study == "MetaHIT"){
    other_samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == "t2d" & !(mod_disease_to_sample$study %in% c("IGC","Richness"))),"sample_id"]
    study_samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == "t2d" & mod_disease_to_sample$study %in% c("IGC","Richness")),"sample_id"]
  }else {
    other_samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == "t2d" & !(mod_disease_to_sample$study == study)),"sample_id"] 
    study_samps <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == "t2d" & mod_disease_to_sample$study == study),"sample_id"] 
  }
  other_kegg.df <- t2d_kegg.df[other_samps,]
  
  other_kegg.rf <- randomForest(status ~ .,data=other_kegg.df,importance=TRUE,ntree=10000)
  
  other.kegg.roc <- roc(t2d_kegg.df[other_samps,"status"],other_kegg.rf$votes[other_samps,2])
  
  pdf(paste0(outdir,"t2d_by_study/l1o_",study,"_roc.pdf"))
  plot(other.kegg.roc,col="black")
  table <- rbind(table,c(paste0("other_studies(",study,")"),"train",mean( predict(other_kegg.rf) != t2d_kegg.df[other_samps,"status"]),other.kegg.roc$auc))
  
  
  pred <- predict(other_kegg.rf,t2d_kegg.df[study_samps,], type="prob", predict.all=T, norm.votes = T)
  
  roc <- roc( t2d_kegg.df[study_samps,"status"], pred$aggregate[study_samps,2])
  
  status_pred <- as.numeric(colnames(pred$aggregate)[1:2][apply(pred$aggregate[study_samps,1:2],1,which.max)])
  oob <- length(which(t2d_kegg.df[study_samps,"status"] != status_pred))/length(t2d_kegg.df[study_samps,"status"])    
  
  plot(roc,add=T,col=colors[study])
  
  
  legend(0.35,0.25,c("other studies",study),col=c("black",colors[study]),lty="solid",lwd=2)
  dev.off()
  
  table <- rbind(table,c(study,"test",oob,roc$auc))
}
colnames(table) <- c("study","datatype","OOB","AUC")




### COLORECTAL CANCER ####
### MODULE ABUNDANCE 
module_abunds <- kegg_abunds[which(kegg_abunds$level == "module"),]
module.df <- dcast(module_abunds,sample_id~name,value.var = "copy_number")
module.df <- merge(mod_disease_to_sample[,c("sample_id","status")],
                   module.df,
                   by="sample_id")
module.df <- unique(module.df)
rownames(module.df) <- module.df$sample_id
module.df <- module.df[,-1]
module.df$status <- as.factor(module.df$status)

#all colorectal cancer subjects
pheno_samps1 <- mod_disease_to_sample[which(mod_disease_to_sample$phenotype == "carcinoma"),"sample_id"]
sub_module1.df <- module.df[pheno_samps1,]
#for colorectal cancer to exclude advanced adenoma subjects
pheno_samps2 <- disease_to_sample[which(disease_to_sample$phenotype == "carcinoma" & disease_to_sample$status %in% c("controls","carcinoma")),"sample_id"]
sub_module2.df <- module.df[pheno_samps2,]

sub_module1.rf <- randomForest(status ~ .,data=sub_module1.df,importance=T,ntree=1000)
sub_module1.roc <- roc(sub_module1.df$status,sub_module1.rf$votes[,2])

sub_module2.rf <- randomForest(status ~ .,data=sub_module2.df,importance=T,ntree=1000)
sub_module2.roc <- roc(sub_module2.df$status,sub_module2.rf$votes[,2])

pdf("./stats_out/ms_review/CC_compare.pdf")
plot(sub_module1.roc,xlim=c(1,0),ylim=c(0,1),xaxs="i",xaxs="i",lty="solid",main="Colorectal cancer",col=colors["carcinoma"])
plot(sub_module2.roc,xlim=c(1,0),ylim=c(0,1),xaxs="i",xaxs="i",lty="dashed",col=colors["carcinoma"],add=T)
legend(0.45,0.2,c("all subjects","exclude adenoma"),col=colors["carcinoma"],lty=c("solid","dashed"),lwd=2)
dev.off()



### BATCH EFFECTS ###
#BETA DIVERSITY
level <- "gene_family"
outpath <- "./stats_out/ms_review/beta_diversity/"
if(!(dir.exists(outpath))){dir.create(outpath,recursive=T)}

dims <- c(1,2,3,4)

abund.map <- create_abund_map(kegg_abunds,"kegg","gene_family","copy_number")
dist <- metaMDS(abund.map,distance="bray",pc=TRUE,k=length(dims))

# ADONIS
mapping_subset <- mod_disease_to_sample[which(mod_disease_to_sample$study %in% c("T2D","MGS","Richness","LC")),]
#mapping_subset <- mod_disease_to_sample[which(mod_disease_to_sample$study %in% c("T2D","LC")),]

ids <- mapping_subset[,"sample_id"]
sub_abunds <- exclude_zero(abund.map[ids,])

metadata_subset <- mod_att.map[which(mod_att.map$sample_id %in% ids),
                               c("subject_id","sample_id","age","bmi","country","sex")]
mapping_subset2 <- merge(mapping_subset,metadata_subset,by="sample_id")
mapping_subset2 <- na.omit(mapping_subset2)
mapping_subset2$country <- as.character(mapping_subset2$country)
mapping_subset2$mod_country <- ifelse(mapping_subset2$country %in% c("Spain","Denmark"),"Europe",mapping_subset2$country)
mapping_subset2$protocol <- ifelse(mapping_subset2$study == "LC","PC","GT")
sub_abunds2 <- sub_abunds[mapping_subset2$sample_id,]

a_table <- NULL
#
a_test1 <- adonis(sub_abunds2~country,data=mapping_subset2,permutations=100)
variable <- row.names(a_test1$aov.tab)[which(!row.names(a_test1$aov.tab) %in% c("Residuals","Total"))]
p_value <- a_test1$aov.tab[variable,"Pr(>F)"]
r_squared <- a_test1$aov.tab[variable,"R2"]
a_table <- rbind(a_table,cbind("abundance~country",variable,p_value,r_squared))
#
a_test2 <- adonis(sub_abunds2~protocol,data=mapping_subset2,permutations=100)
variable <- row.names(a_test2$aov.tab)[which(!row.names(a_test2$aov.tab) %in% c("Residuals","Total"))]
p_value <- a_test2$aov.tab[variable,"Pr(>F)"]
r_squared <- a_test2$aov.tab[variable,"R2"]
a_table <- rbind(a_table,cbind("abundance~protocol",variable,p_value,r_squared))
#
a_test3 <- adonis(sub_abunds2~protocol+country,data=mapping_subset2,permutations=100)
variable <- row.names(a_test3$aov.tab)[which(!row.names(a_test3$aov.tab) %in% c("Residuals","Total"))]
p_value <- a_test3$aov.tab[variable,"Pr(>F)"]
r_squared <- a_test3$aov.tab[variable,"R2"]
a_table <- rbind(a_table,cbind("abundance~protocol+country",variable,p_value,r_squared))
#
a_test4 <- adonis(sub_abunds2~country+protocol,data=mapping_subset2,permutations=100)
variable <- row.names(a_test4$aov.tab)[which(!row.names(a_test4$aov.tab) %in% c("Residuals","Total"))]
p_value <- a_test4$aov.tab[variable,"Pr(>F)"]
r_squared <- a_test4$aov.tab[variable,"R2"]
a_table <- rbind(a_table,cbind("abundance~country+protocol",variable,p_value,r_squared))
#
a_test5 <- adonis(sub_abunds2~study,data=mapping_subset2,permutations=100)
variable <- row.names(a_test5$aov.tab)[which(!row.names(a_test5$aov.tab) %in% c("Residuals","Total"))]
p_value <- a_test5$aov.tab[variable,"Pr(>F)"]
r_squared <- a_test5$aov.tab[variable,"R2"]
a_table <- rbind(a_table,cbind("abundance~study",variable,p_value,r_squared))
#
a_test6 <- adonis(sub_abunds2~country+study,data=mapping_subset2,permutations=100)
variable <- row.names(a_test6$aov.tab)[which(!row.names(a_test6$aov.tab) %in% c("Residuals","Total"))]
p_value <- a_test6$aov.tab[variable,"Pr(>F)"]
r_squared <- a_test6$aov.tab[variable,"R2"]
a_table <- rbind(a_table,cbind("abundance~country+study",variable,p_value,r_squared))
#
a_test7 <- adonis(sub_abunds2~study+country,data=mapping_subset2,permutations=100)
variable <- row.names(a_test7$aov.tab)[which(!row.names(a_test7$aov.tab) %in% c("Residuals","Total"))]
p_value <- a_test7$aov.tab[variable,"Pr(>F)"]
r_squared <- a_test7$aov.tab[variable,"R2"]
a_table <- rbind(a_table,cbind("abundance~study+country",variable,p_value,r_squared))
#
a_test8 <- adonis(sub_abunds2~study+protocol,data=mapping_subset2,permutations=100)
variable <- row.names(a_test8$aov.tab)[which(!row.names(a_test8$aov.tab) %in% c("Residuals","Total"))]
p_value <- a_test8$aov.tab[variable,"Pr(>F)"]
r_squared <- a_test8$aov.tab[variable,"R2"]
a_table <- rbind(a_table,cbind("abundance~study+protocol",variable,p_value,r_squared))
#
a_test9 <- adonis(sub_abunds2~protocol+study,data=mapping_subset2,permutations=100)
variable <- row.names(a_test9$aov.tab)[which(!row.names(a_test9$aov.tab) %in% c("Residuals","Total"))]
p_value <- a_test9$aov.tab[variable,"Pr(>F)"]
r_squared <- a_test9$aov.tab[variable,"R2"]
a_table <- rbind(a_table,cbind("abundance~protocol+study",variable,p_value,r_squared))


#a_info <- cbind(db_level,disease,variable,p_value,r_squared)
#a_table <- rbind(a_table,pheno_info)
#colnames(a_table) <- c("level","disease","variable","p_value","r_squared")
#write.table(a_table,paste0(outpath,level,"_adonis_with_metadata.txt"),quote=F,sep="\t")
#a_table <- read.table(paste0(outpath,level,"_adonis_with_metadata.txt"),header=T,sep="\t")


#NMDS
for(a in 1:length(dims)){
  for(b in 1:length(dims)){
    dim.1 <- dims[a]
    dim.2 <- dims[b]
    if(dim.1 >= dim.2){
      next()
    }
    filename <- paste0( outpath,"batch_effects_nmds-",level,"_","dims_", dim.1,"_", dim.2,".pdf")
        
    data <- dist$points
    
    #build dataframe
    df <- data[,c(dim.1,dim.2)]
    df <- cbind(sample_id=row.names(df),df)
    df <- merge(df,mod_disease_to_sample,by="sample_id")
    df <- merge(df,cbind(sample_id=rownames(metadata.df),metadata.df),by=c("sample_id","study","status"))
    df <- merge(df, mapping_subset2,by=c("sample_id","study","phenotype","status","age","sex","country","pheno.status"))
    
    #only get controls
    df <- df[which(df$status == "0"),c("MDS1","MDS2","sample_id","study","phenotype","country","mod_country","protocol")]
    
    #fix data classes
    df[,1] <- as.numeric(as.character(df[,1]))
    df[,2] <- as.numeric(as.character(df[,2]))
    
    #diff protocols from china
    sub_df <- df[which(df$study %in% c("T2D","MGS","Richness","LC")),]
    
    country_plot <- ggplot(data=sub_df,aes(x=sub_df[,1],y=sub_df[,2],shape=protocol,fill=mod_country)) + 
      geom_point(size=2) +
      scale_shape_manual(values = c(21,24),name="Protocol") + 
      xlab(colnames(df)[1]) + ylab(colnames(df)[2]) + 
      scale_fill_manual("Region",values=c("cadetblue3","indianred3")) +
      guides(fill = guide_legend(override.aes=list(shape=21))) 
    
    study_plot <- ggplot(data=sub_df,aes(x=sub_df[,1],y=sub_df[,2],shape=protocol,fill=study)) + 
      geom_point(size=2) +
      scale_shape_manual("Protocol",values = c(21,24)) + 
      xlab(colnames(df)[1]) + ylab(colnames(df)[2]) +
      guides(fill = guide_legend(override.aes=list(shape=21),
                                 title="Study"))
    
    full_plot <- plot_grid(plotlist = list(country_plot,study_plot),nrow=1,labels = c("A","B"),align = "h")  
    pdf(filename,width=12,height=5)
    print(full_plot)
    dev.off()
  
  }
}

