library(DESeq2)
library(tidyverse)
library(tidyr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(pheatmap)

##reading in sequencing core produced normalized tpm##
counts_data<-read.table("C:/Users/kelly/Internship/clock_project/data/CD02JHU503/CD02JHU503_000_analysis/RSEM/tables/genes/genes_tpm_all_samples_norm.txt")
tpm_counts<-counts_data
colnames(tpm_counts)<-tpm_counts[1,]
tpm_counts<-tpm_counts[2:nrow(tpm_counts),]

tpm_counts<-as.data.frame(tpm_counts)
counts_data<-tpm_counts
counts_data[,2:ncol(counts_data)]<-apply(counts_data[,2:ncol(counts_data)],2,as.numeric)

names<-tpm_counts$gene_id
tpm_counts<-tpm_counts[,2:ncol(tpm_counts)]
tpm_counts<-apply(tpm_counts,2,as.numeric)
rownames(tpm_counts)<-names

##reading in sequencing core produced expected counts tpm##
setwd("C:/Users/kelly/Internship/clock_project/data/CD02JHU503/CD02JHU503_000_analysis")
raw_counts_data<- read.table("RSEM/tables/genes/genes_expected_count_all_samples.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
counts_data<-raw_counts_data

##Exploring the data##
#Violin Distribution
long_data<-counts_data[,2:ncol(counts_data)] %>%
  pivot_longer(cols = everything(), 
               names_to = "Condition", 
               values_to = "Value")
ggplot(long_data,aes(x=Condition,y=log2(Value+1),fill=Condition))+
  geom_violin(scale="width",width=0.8)+theme(axis.text.x = element_text(angle=60, hjust=1))+
  theme(legend.position="none")+labs(title="CD02JHU513 RAW Counts Distribution",x="Condition",y="Log2 Counts+1")


#MA Plot (customizable to specific conditions)
M<-log2(counts_data$expected_count_ac3stblMycMut2+1)-log2(counts_data$expected_count_2_1stblEV1+1)
A<-0.5*(log2(counts_data$expected_count_ac3stblMycMut2+1)+log2(counts_data$expected_count_2_1stblEV1+1))

plot(A,M,pch=20,col="blue",xlab="Average Expression", ylab="Log2 Fold Change",main="ac3stblMycMut2 vs. 2_1stblEv1 MA Plot")
loess_fit<-loess(M~A)
lines(A,predict(loess_fit),col="black",lwd=2)

#filter out low quality genes + loess fit
counts_log2<-counts_data
counts_log2[,2:ncol(counts_data)]<-log2(counts_data[,2:ncol(counts_data)]+1)
zero_counts<-apply(counts_log2[,2:ncol(counts_log2)],1,function(row)(sum(row==0)))
counts_log2<-counts_log2[zero_counts<=12,]
averages<-rowMeans(counts_log2[,2:ncol(counts_log2)])
vars<-rowVars(apply(counts_log2[,2:ncol(counts_log2)],2,as.numeric))
pca_loess_fit<-loess(vars~averages)
plot(averages,vars,pch=20,col="blue",xlab="Average Expression", ylab="Log2 Fold Change",main="TPM Genes MA Plot")
lines(vars,predict(pca_loess_fit),col="black",lwd=2)

#pca reduction plot
a_v<-vars/averages
highly_variable_gene_indexes<-order(a_v,decreasing=TRUE)[1:2000]
highly_variable_genes<-counts_log2$gene_id[highly_variable_gene_indexes]
pca_df<-counts_log2[highly_variable_gene_indexes,]
pca_df<-pca_df[,2:ncol(pca_df)]
rownames(pca_df)<-highly_variable_genes
pca_df<-apply(pca_df,2,as.numeric)
for(i in 1:ncol(pca_df)){
  column_mean<-mean(pca_df[,i])
  column_sd<-sd(pca_df[,i])
  pca_df[,i]<-(pca_df[,i]-column_mean)/column_sd
}
pca_df<-t(pca_df)
reduced_df<-prcomp(apply(pca_df,2,as.numeric),rank=2,scale=TRUE)
reduced_df<-as.data.frame(reduced_df$x)
cols<-c(rep("stbl_ev",2),rep("stbl_mut",2),rep("stbl_wt",2),rep("ac3_ev",2),rep("ac3_mut",2),rep("ac3_wt",2))
ggplot(reduced_df,aes(x=PC1,y=PC2,color=cols))+geom_point(size=3)+labs(title="RAW PCA Reduction of Samples",x="PC1",y="PC2")

#basic heatmap of initial top 2000 highly variable genes
pca_df<-t(pca_df)
rownames(pca_df)<-highly_variable_genes
for(i in 1:nrow(pca_df)){
  row_mean<-mean(pca_df[i,])
  row_sd<-sd(pca_df[i,])
  pca_df[i,]<-(pca_df[i,]-row_mean)/row_sd
}
intervals<-c(-2.5,-0.75,-0.5,-0.25,0,0.5,3)
colors<-c("#2E7EB2"  ,"#B2E1EE", "#D3F9FD", "#FCF1E4", "#FD9963", "#B93D13", "#2F010C")
colors<-c("#2E7EB2" ,"#4F97C1", "#70AFD0", "#91C8DF" ,"#B2E1EE", "#D3F9FD", "#FCF1E4","#FDBE8F", "#FD9963", "#FC7436" ,"#E75119", "#B93D13", "#8B2911", "#5D160E", "#2F010C")
colors<-c("#2E7EB2" ,"#91C8DF" ,"#B2E1EE", "#D3F9FD", "#FCF1E4", "#FD9963", "#B93D13", "#8B2911","#2F010C")
map<-pheatmap(pca_df,breaks=intervals,color=colors,main="Top 2000 Highly Variable Genes RAW", fontsize_row=2,fontsize_col=5,angle_col=45)

#save highly variable genes
setwd("C:/Users/kelly/Internship/clock_project/kji_DESeq_Analysis/DESeq_results_CD02JHU513")
write.csv(highly_variable_genes,"tpm_top2000_highly_variable_genes.csv")


###ACTUAL DESEQ ANALYSIS###
#designing the colData, always reread in the raw counts

samples<-colnames(raw_counts_data[c(4,5,8,9,12,13,6,7,10,11,14,15)])
condition<-c(rep("OFF",6),rep("ON",6))
colData<-as.data.frame(cbind(samples,condition))

raw_counts_data_1<-raw_counts_data
raw_counts_data<-raw_counts_data[,2:ncol(raw_counts_data)]
raw_counts_data<-apply(raw_counts_data,2,as.integer)
rownames(raw_counts_data)<-raw_counts_data_1[,1]

#SPLIT INTO GROUPS & CHECK ORDER

dds<-DESeqDataSetFromMatrix(countData=raw_counts_data[,c(3,4,7,8,11,12,5,6,9,10,13,14)],
                       colData = colData, #sample name left and a design factor on right, with rowname as condition folder,
                       design = ~condition)
dds

#remove rows with low gene counts
keep<-rowSums(counts(dds))>=10 #minimum baseline
dds<-dds[keep,]
dds


#set factor level (compare groups)
relevel(dds$condition, ref="OFF") #off or on

#RUN DESEQ
dds<-DESeq(dds)
res<-results(dds)
res_df<-as.data.frame(res)
top_genes<-res_df[order(res_df$padj,decreasing = FALSE),][1:10,]


##DESEQ Results##
#MA Plot & Volcano Plot
#ma with significant genes
alpha <- 0.05
res_df$significant <- res$padj < alpha & !is.na(res$padj)
plotMA(res,main="MA Plot of DEGs for EV vs. WT",colLine="black")
with(top_genes, text(x = baseMean, y = log2FoldChange, 
                                    labels = rownames(top_genes), cex = 0.6, pos = 3))

#new MA plot with labels (ggplot)
res_df$label<-ifelse(rownames(res_df) %in% rownames(top_genes), rownames(res_df), NA)
ggplot(res_df,aes(x=log10(baseMean),y=log2FoldChange,color=significant))+scale_color_manual(values=c("black","blue"))+geom_point(alpha=0.5)+geom_text(aes(label = label), size=2.5,vjust = -0.5, hjust = 1.5, color = "black")+labs(x = "Log10 Base Mean", y = "Log2 Fold Change", 
      title = "EV vs. WT MA Plot", color = "Significance") +
  geom_hline(yintercept = 0, linetype = "solid")

#volcano
ggplot(res_df, aes(x = log2FoldChange, y = -log(pvalue))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("black", "red")) +
  labs(
    title = "EV vs. WT Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 P-value"
  ) +
  theme(legend.position = "none")+  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)), size = 3, box.padding = 0.5,max.overlaps=20)

#organizing data by applying log2 and noramlizing by gene!
res_df_ordered<-res_df[order(res_df$padj,decreasing=FALSE),]
setwd("C:/Users/kelly/Internship/clock_project/kji_DESeq_Analysis/DESeq_results_CD02JHU503")
write.csv(res_df_ordered[,1:ncol(res_df_ordered)],"OFF_vs_ON_DEG_df_ordered.csv")
res_df_ordered<-res_df[order(abs(res_df$log2FoldChange),decreasing=TRUE),]
res_df_ordered<-res_df_ordered[which(res_df_ordered$significant==TRUE),]
variable_counts_tpm<-tpm_counts[which(rownames(tpm_counts)%in%rownames(res_df_ordered)[1:30]),c(1,2,7,8,5,6,11,12)]
variable_counts_raw<-raw_counts_data[which(rownames(raw_counts_data)%in%rownames(res_df_ordered)[1:30]),c(1,2,7,8,5,6,11,12)]
names_tpm<-rownames(variable_counts_tpm)
names_raw<-rownames(variable_counts_raw)
variable_counts_tpm<-log2(variable_counts_tpm+1)
variable_counts_tpm<-apply(variable_counts_tpm,2,as.numeric)
rownames(variable_counts_tpm)<-names_tpm

variable_counts_raw<-log2(variable_counts_raw+1)
variable_counts_raw<-apply(variable_counts_raw,2,as.numeric)
rownames(variable_counts_raw)<-names_raw

#normalizing by gene
for(i in 1:nrow(variable_counts_tpm)){
  row_mean<-mean(variable_counts_tpm[i,])
  row_sd<-sd(variable_counts_tpm[i,])
  variable_counts_tpm[i,]<-(variable_counts_tpm[i,]-row_mean)/row_sd
}
for(i in 1:nrow(variable_counts_raw)){
  row_mean<-mean(variable_counts_raw[i,])
  row_sd<-sd(variable_counts_raw[i,])
  variable_counts_raw[i,]<-(variable_counts_raw[i,]-row_mean)/row_sd
}

#produce heatmaps that are NOT CLUSTERED
pheatmap(variable_counts_raw,main="EV vs. WT RAW",angle_col=315,cluster_cols=FALSE,cluster_rows=FALSE)
pheatmap(variable_counts_tpm,main="EV vs. WT TPM",angle_col=315,cluster_cols=FALSE,cluster_rows=FALSE)

#produce clustered heatmaps
pheatmap(variable_counts_raw,main="EV vs. WT RAW",angle_col=315)
pheatmap(variable_counts_tpm,main="EV vs. WT TPM",angle_col=315)
