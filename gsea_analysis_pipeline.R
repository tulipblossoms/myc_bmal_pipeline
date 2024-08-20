library(fgsea)
library(ggplot2)
library(gridExtra)

#read in pathways
pathways<-gmtPathways("C:/Users/kelly/Internship/clock_project/scripts/mh.all.v2024.1.Mm.symbols.gmt")
#read in gene_list
setwd("C:/Users/kelly/Internship/clock_project/kji_DESeq_Analysis/DESeq_results_CD02JHU513")
deseq2_results<-read.csv("EVvs_WT_DEG_df_ordered.csv")

#sort vector
print(length(which(deseq2_results$significant==TRUE)))
gene_list<-deseq2_results$log2FoldChange
names(gene_list)<-deseq2_results$X
gene_list<-sort(gene_list,decreasing=TRUE)

#run gsea
fgsea_results<-fgsea(pathways=pathways,stats=gene_list,minSize=15,maxSize=500)

#view results
fgsea_results<-fgsea_results[order(fgsea_results$pval,decreasing=FALSE)]

#enrichment plots
enrichment_plots<-list()
for(i in 1:6){
  pathway<-fgsea_results$pathway[i]
  enrichment_plots[[i]]<-plotEnrichment(pathways[[pathway]],gene_list)+ggtitle(pathway)+theme(plot.title = element_text(size = 5))  # Adjust the size here
}

plot0<-grid.arrange(grobs=enrichment_plots, ncol=3)
ggsave("C:/Users/kelly/Internship/clock_project/kji_GSEA_Analysis/CD02JHU513_GSEA/plots/EV_vs_WT_ENRICHMENT.pdf",plot0,width=7.33,height=5.81)

#table plot (didn't create for this analyses)
#topPathwaysUp<-fgsea_results[ES>0][head(order(pval),n=10),pathway]
#topPathwaysDown<-fgsea_results[ES<0][head(order(pval),n=10),pathway]
#topPathways<-c(topPathwaysUp,rev(topPathwaysDown))
#plotGseaTable(pathways[topPathways],gene_list,fgsea_results,gseaParam=0.5)

#bar chart
topPathways<-fgsea_results[order(fgsea_results$pval,decreasing=FALSE)][1:10,]
topPathways$Color <- ifelse(topPathways$pval < 0.05, "significant", "not_significant")
setwd("C:/Users/kelly/Internship/clock_project/kji_GSEA_Analysis/CD02JHU513_GSEA/plots")
plot1<-ggplot(topPathways,aes(x=reorder(pathway,NES,decreasing=TRUE),y=NES,fill=Color))+
  geom_bar(stat="identity")+
  coord_flip()+scale_fill_manual(values = c("significant" = "#FFB6C1", "not_significant" = "#ADD8E6")) +
  geom_text(aes(label=sprintf("p=%.2g",pval)),vjust=-0.5,hjust=ifelse(topPathways$NES > 0, 1, -1))+
  labs(x="Pathway",y="Normalized ES",title="Top 10 Enriched Pathways")+
  theme(axis.text.x=element_text(angle=45,hjust=1))
plot1
ggsave("EV_vs_WT_GSEA_NES.pdf", plot = plot1, width = 8, height = 6)

plot2<-ggplot(topPathways,aes(x=reorder(pathway,ES,decreasing=TRUE),y=ES,fill=Color))+
  geom_bar(stat="identity")+
  coord_flip()+scale_fill_manual(values = c("significant" = "#FFB6C1", "not_significant" = "#ADD8E6")) +
  geom_text(aes(label=sprintf("p=%.2g",pval)),vjust=-0.5,hjust=ifelse(topPathways$ES > 0, 1, -1))+
  labs(x="Pathway",y="ES",title="Top 10 Enriched Pathways")+
  theme(axis.text.x=element_text(angle=45,hjust=1))
plot2
ggsave("EV_vs_WT_GSEA_ES.pdf", plot = plot2, width = 8, height = 6)

fgsea_results$leadingEdge <- sapply(fgsea_results$leadingEdge, paste, collapse = ", ")
write.csv(as.data.frame(fgsea_results),"C:/Users/kelly/Internship/clock_project/kji_GSEA_Analysis/CD02JHU513_GSEA/tables/EV_vs_WT_GSEA_table.csv")
