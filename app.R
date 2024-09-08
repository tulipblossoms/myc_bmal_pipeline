#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinydashboard)
library(shinythemes)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(DT)

mycER_tpm<-read.csv("MycER_Integrated_tpm.csv")
rownames(mycER_tpm)<-mycER_tpm$X
mycER_tpm<-mycER_tpm[,2:ncol(mycER_tpm)]
mycER_raw<-read.csv("MycER_Integrated_raw.csv")
rownames(mycER_raw)<-mycER_raw$X
mycER_raw<-mycER_raw[,2:ncol(mycER_raw)]
stable_tpm<-read.csv("myc_stable_tpm.csv")
rownames(stable_tpm)<-stable_tpm$X
stable_tpm<-stable_tpm[,2:ncol(stable_tpm)]
stable_raw<-read.csv("myc_stable_raw.csv")
rownames(stable_raw)<-stable_raw$X
stable_raw<-stable_raw[,2:ncol(stable_raw)]

ui <- dashboardPage(skin="purple", dashboardHeader(title="Bulk RNA-seq Myc ER/Stable Systems DEGs",
                                                   titleWidth=500
),
dashboardSidebar(
  sidebarMenu(
    menuItem("home", tabName = "Home", icon = icon("dashboard")),
    menuItem("user guide",tabName="user_guide",icon=icon("tree"))    
  )
),
dashboardBody(
  tabItems(
    tabItem("Home",fluidPage(
            selectInput("table_choice","Select a comparison:",
                        choices=list("MycER 2.1 WT ON vs. 2.1 WT OFF"="MycER_2.1_WT_ON_vs._2.1_WT_OFF",
                                     "MycER 2.1 Myc WT ON vs. 2.1 Myc EV ON"="MycER_2.1_Myc_WT_ON_vs._2.1_Myc_EV_ON",
                                     "MycER ac3 Myc WT ON vs. ac3 Myc WT OFF"="MycER_ac3_Myc_WT_ON_vs._ac3_Myc_WT_OFF",
                                     "MycER ac3 Myc WT ON vs. ac3 Myc EV ON"="MycER_ac3_Myc_WT_ON_vs._ac3_Myc_EV_ON",
                                     "MycER 2.1 Myc WT ON vs. ac3 Myc WT ON"="MycER_2.1_Myc_WT_ON_vs._ac3_Myc_WT_ON",
                                     "MycER 2.1 Myc WT OFF vs. ac3 Myc WT OFF"="MycER_2.1_Myc_WT_OFF_vs._ac3_Myc_WT_OFF",
                                     "MycER 2.1 Myc EV ON vs. ac3 Myc EV ON"="MycER_2.1_Myc_EV_ON_vs._ac3_Myc_EV_ON",
                                     "Stable 2.1 Myc WT vs. ac3 Myc WT"="Stable_2.1_Myc_WT_vs._ac3_Myc_WT",
                                     "Stable 2.1 Myc WT vs. 2.1 Myc EV"="Stable_2.1_Myc_WT_vs._2.1_Myc_EV",
                                     "Stable ac3 Myc WT vs. ac3 Myc EV"="Stable_ac3_Myc_WT_vs._ac3_Myc_EV")),
            box(width=12,solidHeader=TRUE,status="primary",title="Select Genes:",collapsible = TRUE,
                DTOutput("deg_table"),
                downloadButton("download_data", "Download Table")),
            box(width=8,solidHeader=TRUE,status="primary",title="Heatmap",
                plotOutput("heatmap_plot")),
            box(width=4,solidHeader=TRUE,status="primary",title="Customize",
                verbatimTextOutput("selected_rows"),
                selectInput("raw_expected","Select a dataset:",
                            choices=list("tpm counts"="tpm_counts","raw counts"="raw_counts")),
                selectInput("clustered","Clustered Heatmap?",
                            choices=list("cluster genes ONLY"="rows_cluster","cluster samples ONLY"="cols_cluster","clustered genes and samples"="all_cluster","no clustering"="no_cluster"))
    ))),
    tabItem("user_guide",
            verbatimTextOutput("instructions"))
  )))








# Define server logic required to draw a histogram
server <- function(input, output) {
selected_data_table<-reactive({
  file_name<-paste0(input$table_choice,".csv")
  data<-read.csv(file_name)%>%select(1:8)
  dt<-datatable(data)
  return(dt)
})

output$deg_table<-renderDT({
  selected_data_table()
})
output$download_data <- downloadHandler(
  filename=function(){
    paste0("deg_table",".csv")
  },
  content=function(file){
    write.csv(data, file)
  }
)

output$selected_rows <- renderPrint({
  selected_rows <- input$deg_table_rows_selected
  if(length(selected_rows)<=1){
    print("You must select at least 2 rows to generate a heatmap.")
  }else{
    print("Customize heatmap below!")
  }
})

output$heatmap_plot<-renderPlot({
  if(grepl("MycER",input$table_choice)){
    tpm_counts<-mycER_tpm
    raw_counts_data<-mycER_raw
    annotation_col<-data.frame(
      clock=factor(c(rep("2.1",14),rep("ac3",14))),
      myc=factor(c(rep("p",2),rep("EV",4),rep("MUT",4),rep("WT",4),rep("p",2),rep("EV",4),rep("MUT",4),rep("WT",4))),
      induction=factor(c(rep("na",2),rep("off",2),rep("on",2),rep("off",2),rep("on",2),rep('off',2),rep("on",2),rep("na",2),rep("off",2),rep("on",2),rep("off",2),rep("on",2),rep("off",2),rep("on",2)))
    )
  }else{
    tpm_counts<-stable_tpm
    raw_counts_data<-stable_raw
    annotation_col<-data.frame(
      clock=factor(c(rep("2.1",6),rep("ac3",6))),
      myc=factor(c(rep("EV",2),rep("MUT",2),rep("WT",2),rep("EV",2),rep("MUT",2),rep("WT",2)))
    )
  }
  annotation_colors <- list(
    clock = c("2.1" = "pink", "ac3" = "salmon"),
    myc = c("p" = "purple", "EV" = "lightgreen","MUT"="green","WT"="black"),
    induction=c("na"="gray","off"="darkred","on"="yellow")
  )
  selected_rows <- input$deg_table_rows_selected
  if(length(selected_rows)>1){
    res_df_ordered<-as.data.frame(selected_data_table()$x$data[selected_rows,])
    res_df_ordered<-res_df_ordered[order(abs(res_df_ordered$log2FoldChange),decreasing=TRUE),]
    res_df_ordered<-res_df_ordered[which(res_df_ordered$significant==TRUE),]
    genes<-res_df_ordered$X
    variable_counts_tpm<-tpm_counts[which(rownames(tpm_counts)%in%genes),]
    variable_counts_raw<-raw_counts_data[which(rownames(raw_counts_data)%in%genes),]
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
    if(input$raw_expected=="raw_counts"){
      if(input$clustered=="cols_cluster"){
        pheatmap(variable_counts_raw,main=paste0(input$table_choice, " Raw Counts"),angle_col=315,cluster_rows=FALSE)
      }else if(input$clustered=="all_cluster"){
        pheatmap(variable_counts_raw,main=paste0(input$table_choice, " Raw Counts"),angle_col=315)
      }else if(input$clustered=="no_cluster"){
        rownames(annotation_col) <- colnames(variable_counts_raw)
        pheatmap(variable_counts_raw,main=paste0(input$table_choice, " Raw Counts"),angle_col=315,cluster_cols=FALSE,cluster_rows=FALSE,annotation_colors=annotation_colors,annotation_col=annotation_col)
      }else{
        pheatmap(variable_counts_raw,main=paste0(input$table_choice, " Raw Counts"),angle_col=315,cluster_cols=FALSE)
      }
    }else{
      if(input$clustered=="cols_cluster"){
        pheatmap(variable_counts_tpm,main=paste0(input$table_choice, " TPM Counts"),angle_col=315,cluster_rows=FALSE)
      }else if(input$clustered=="all_cluster"){
        pheatmap(variable_counts_tpm,main=paste0(input$table_choice, " TPM Counts"),angle_col=315)
      }else if(input$clustered=="no_cluster"){
        rownames(annotation_col) <- colnames(variable_counts_tpm)
        pheatmap(variable_counts_tpm,main=paste0(input$table_choice, " TPM Counts"),angle_col=315,cluster_cols=FALSE,cluster_rows=FALSE,annotation_colors = annotation_colors,annotation_col=annotation_col)
      }else{
        pheatmap(variable_counts_tpm,main=paste0(input$table_choice, " TPM Counts"),angle_col=315,cluster_cols=FALSE)
      }
    }
  }
})

output$instructions<-renderPrint({
  lines<-c("Welcome! This R Shiny App will allow you to interactively explore differentially expressed genes for Bulk-RNA sequencing.",
  "",
  "MycER Integrated and Stable Myc Bulk-RNA seq data was used. DESeq was run on all comparisons, and the results are shown as the tables on the home page. You may download all tables onto your computer.",
  "Below the table, one should see a box to view heatmaps and customizable conditions.",
  "The 'Select a Dataset' option allows you to choose to visualize either raw gene counts or transcripts per million.",
  "The 'Clustered Heatmap' option allows you to choose to cluster just the genes, samples, both, or neither (ordered with annotations).",
  "You must select at least 2 genes from the table (by clicking on the table entries) to produce a heatmap.",
  "",
  "If you have any questions or concerns, feel free to email Kelly at kji3@jhu.edu or kxji19@gmail.com.",
  "All raw code for analyses and this dashboard can be found on github at this link: https://github.com/tulipblossoms/myc_bmal_pipeline",
  "",
  "Thank you for using Myc Bulk RNA-seq Shiny!"
  )
  cat(lines, sep = "\n")
})
}

shinyApp(ui = ui, server = server)
