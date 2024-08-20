# myc_bmal_pipeline
Below is the code used to generate DESEQ + GSEA results for the Yumm2.1/Yumm 2.1 ac3 stable and MycER cells.

The code structure is as follows:
deseq pipeline
1. Data Exploration (side-by-side violin plots, MA plots, PCA plots, and heatmaps were produced)
2. DESEQ analysis (DESEQ was run)
3. DESEQ Results + Visualization (MA plots, volcano plots, and more heatmaps generated using ggplot2 and pheatmap [for heatmaps])

gsea pipeline
1. Read in pathways (downloaded prior to analyses) and deseq gene sets
2. Ran gsea using fgsea
3. Made enrichment plots and bar charts
