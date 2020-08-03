# This script was used to do integrative analysis of transcriptome and metabolome data
# The script was modified from the manual of mixOmics package.
# This script was for paper 'Deciphering mechanisms of Oenococcus oeni combating acidic stress by integrative omics and genome-scale metabolic modeling' 

library(mixOmics) # import mixOmics package

# import transcriptome and metabolome data
transcriptome <- read.table(file.choose(), header = T) # Filename: expression_DEG.txt
metabolome <- read.table(file.choose(), header = T)    # Filename: metabolite_samples_i_filtered.txt (i=1-20)

# extract names(IDs) of each mRNA and metabolite feature.
rownames(transcriptome) <- transcriptome[,1]
rownames(metabolome) <- metabolome[,1]

# experimental design, i.e. assign each coloum into its treatment group
group_ph = c("pH4.8","pH4.8","pH4.8","pH4.8",
          "pH4.8","pH4.8","pH4.8","pH4.8",
          "pH4.8","pH3.0","pH3.0","pH3.0",
          "pH3.0","pH3.0","pH3.0")
group_time = c("0h","0h","0h","1h","1h","1h","3h","3h","3h",
               "1h","1h","1h","3h","3h","3h")

# SPLS-DA analysis
X <- list(mRNA = t(transcriptome[,2:16]),metabolite = t(metabolome[,2:16]))
MyResult.diablo <- block.splsda(X, group_ph)

# sample scatterplot from displaying the first component in transcriptome & metabolitome (upper diagonal plot) 
# and Pearson correlation between them
plotDiablo(MyResult.diablo, ncomp = 1)

# sample plot per data set
plotIndiv(MyResult.diablo, ind.names = FALSE, ellipse = TRUE, legend = TRUE, pch = c(1,1))

# Correlation Circle plot representing each type of selected features
plotVar(MyResult.diablo, cutoff=0.80, style = 'graphics', pch = c(1,1), cex = c(1.5,1.5), rad.in = 0)

# Clustered Heatmap (Euclidian distance, Complete linkage)
cimDiablo(MyResult.diablo, color.blocks = c('darkorchid', 'brown1'))

# generate cicros Plot
circosPlot(MyResult.diablo, cutoff = 0.8, line = TRUE, 
           color.blocks= c('darkorchid', 'brown1'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)

plotLoadings(MyResult.diablo, comp = 2, contrib = 'max', method = 'median')
auc.splsda = auroc(MyResult.diablo, roc.block = "mRNA", roc.comp = 2)

network(MyResult.diablo, blocks = c(1,2),
        color.node = c('darkorchid', 'lightgreen'), cutoff = 0.4)

