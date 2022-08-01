#### Isoform identification and analysis ####
#############################################
#############################################



# Install required packages 
# install.packages("tgp")
# install.packages("dplyr")
# BiocManager::install("IsoformSwitchAnalyzeR")
# BiocManager::install("cummeRbund")
# BiocManager::install("topGO")
# BiocManager::install("Mus.musculus")
# BiocManager::install("Rgraphviz")



# Load required packages 
library(IsoformSwitchAnalyzeR)
packageVersion('IsoformSwitchAnalyzeR')
library(dplyr)
library(tgp)
library(ggplot2)
library(topGO)
library(Mus.musculus)
library(Rgraphviz)
library(Repitools)


# Load data (Cufflinks/Cuffdiff)
## Input: Cufflinks->Cuffmerge->Cuffdiff quantified RNA-Seq data  
SwitchList <- readRDS("~/Desktop/Isoform_analysis_data/SwitchList.rds") 

## Load result objects
gpr_list_IF <- readRDS("~/Desktop/Isoform_analysis_data/gpr_list_IF.rds")
IsoformCountFilteredGPR <- readRDS("~/Desktop/Isoform_analysis_data/IsoformCountFilteredGPR.rds")
gpr_filtered_IF <- readRDS("~/Desktop/Isoform_analysis_data/gpr_filtered_IF.rds")
## Load results objects, filtered by Isoforms of Interest
FPKM_filtered <- readRDS("~/Desktop/Isoform_analysis_data/FPKM_filtered.rds")
IF_filtered <- readRDS("~/Desktop/Isoform_analysis_data/IF_filtered.rds")
Isoform_features_filtered <- readRDS("~/Desktop/Isoform_analysis_data/Isoform_features_filtered.rds")

# SwitchList <- importCufflinksFiles(
#   pathToGTF                      = "/Users/jonathan/Documents/Biologie_MSc/WS_21:22/Project_module_Beyer/Data/merge_out_all/merged.gtf",
#   pathToGeneDEanalysis           = "/Users/jonathan/Documents/Biologie_MSc/WS_21:22/Project_module_Beyer/Data/diff_out_all/gene_exp.diff",
#   pathToIsoformDEanalysis        = "/Users/jonathan/Documents/Biologie_MSc/WS_21:22/Project_module_Beyer/Data/diff_out_all/isoform_exp.diff",
#   pathToGeneFPKMtracking         = "/Users/jonathan/Documents/Biologie_MSc/WS_21:22/Project_module_Beyer/Data/diff_out_all/genes.fpkm_tracking",
#   pathToIsoformFPKMtracking      = "/Users/jonathan/Documents/Biologie_MSc/WS_21:22/Project_module_Beyer/Data/diff_out_all/isoforms.fpkm_tracking",
#   pathToIsoformReadGroupTracking = "/Users/jonathan/Documents/Biologie_MSc/WS_21:22/Project_module_Beyer/Data/diff_out_all/isoforms.read_group_tracking",
#   pathToSplicingAnalysis         = "/Users/jonathan/Documents/Biologie_MSc/WS_21:22/Project_module_Beyer/Data/diff_out_all/splicing.diff",
#   pathToReadGroups               = "/Users/jonathan/Documents/Biologie_MSc/WS_21:22/Project_module_Beyer/Data/diff_out_all/read_groups.info",
#   pathToRunInfo                  = "/Users/jonathan/Documents/Biologie_MSc/WS_21:22/Project_module_Beyer/Data/diff_out_all/run.info",
#   fixCufflinksAnnotationProblem  = TRUE, 
#   estimateDifferentialGeneRange  = TRUE, 
#   addIFmatrix                    = TRUE, 
#   isoformNtFasta                 = NULL, 
#   quiet                          = FALSE
# ) 



# Options to check a SwitchList
head(SwitchList)
names(SwitchList)
View(SwitchList)
View(SwitchList$conditions)
View(SwitchList$isoformFeatures)
View(SwitchList$isoformRepIF)
View(SwitchList$isoformCountMatrix)
View(SwitchList$isoformSwitchAnalysis)
View(SwitchList$isoformRepExpression)
View(SwitchList$exons)


# Pre-processing/Filtering
## Select isoforms that are present in >= 5 samples (we may use "conditions" instead?) 
### Filter via FPKM threshold (IsoformRepExpression = Replicate isoform abundance matrix (not log-transformed)) 
### Subset isoforms with expression > 1.0 in at least 5 samples:
IsoformRepExpression <- data.frame(SwitchList$isoformRepExpression[, -1], row.names = SwitchList$isoformRepExpression[, 1])
IsoformExpressionFiltered <- IsoformRepExpression[rowSums(IsoformRepExpression > 1.0) >= 5, ] 

### Alternatively: filter via NAs in Isoform Fraction (IF)
# SwitchList$isoformRepIF[SwitchList$isoformRepIF < 0.0001] = NA
# IsoformRepIFFiltered <- SwitchList$isoformRepIF[rowSums(is.na(SwitchList$isoformRepIF)) < 16, ]
# IsoformRepIFFiltered <- IsoformRepIFFiltered[-c(15,16)]


## Select isoforms with change in isoform usage over certain threshold
### The difference in isoform usage is used to measure effect size (like fold changes in gene expression analysis)
### Difference in isoform usage == difference in isoform fraction (dIF): dIF = IF2 - IF1 
### Isoform fraction (IF) = isoform_expression / gene_expression
### FYI: "Consequently, the difference in isoform usage is quantified as the difference in isoform fraction (dIF) calculated as IF2 - IF1, 
### and these dIF are used to measure the effect size (like fold changes are in gene/isoform expression analysis)."
IsoformUsageFiltered <- SwitchList$isoformFeatures[which(SwitchList$isoformFeatures$isoform_id %in% row.names(IsoformExpressionFiltered)), ]
IsoformUsageFiltered <- IsoformUsageFiltered[, c(3, 36)]
#IsoformUsageFiltered <- data.frame(IsoformUsageFiltered[, -1], row.names = IsoformUsageFiltered[, 1]) # issue: duplicate row.names 

### Select isoforms that show change in dIF threshold of > 0.1 or < -0.1 in at least one comparison of conditions, i.e. time points 
IsoformUsageFiltered <- IsoformUsageFiltered[which(IsoformUsageFiltered$dIF > 0.2 | IsoformUsageFiltered$dIF < -0.2), ]



# Gaussian Process Regression 
## Select isoform usage == isoform fractions (IF) based on filtered isoforms, i. e. present in >= 5 samples & difference in isoform usage > 10%
IsoformFractionFiltered <- SwitchList$isoformRepIF[which(SwitchList$isoformRepIF$isoform_id %in% IsoformUsageFiltered$isoform_id), ]
IsoformFractionFiltered <- data.frame(IsoformFractionFiltered[, -1], row.names = IsoformFractionFiltered[, 1])
IsoformFractionFiltered <- IsoformFractionFiltered[, c(19:21, 1:18)]
IsoformFractionFiltered <- IsoformFractionFiltered[, c(14:18, 1:13, 19:21)]
IsoformFractionFiltered <- IsoformFractionFiltered[, c(4, 5, 2, 3, 1, 6:21)] 
IsoformFractionFiltered[is.na(IsoformFractionFiltered)] <- 0 

### Prepare result matrix for GPR
IsoformFractionFilteredGPR <- matrix(NA, nrow = dim(IsoformFractionFiltered)[1], ncol=21)
row.names(IsoformFractionFilteredGPR) <- row.names(IsoformFractionFiltered)
colnames(IsoformFractionFilteredGPR) <- c("RC9_2iL_1", "RC9_2iL_2",
                                           "RC9_2i_1", "RC9_2i_2", "RC9_2h",
                                           "RC9_4h", "RC9_6h", "RC9_8h",
                                           "RC9_10h", "RC9_12h", "RC9_14h",
                                           "RC9_16h", "RC9_18h", "RC9_20h",
                                           "RC9_22h", "RC9_24h", "RC9_26h",
                                           "RC9_28h", "RC9_30h",
                                           "RC9_32h_1","RC9_32h_2")
### Include 2iL IF, in case they are needed later on
IsoformFractionFiltered <- as.matrix(IsoformFractionFiltered)

for(d in 1:dim(IsoformFractionFiltered)[1]){
  IsoformFractionFilteredGPR[d,c(1:2)] <- IsoformFractionFiltered[d,c(1:2)]
}  

### Prepare result matrices for GPR
q1_gpr <- matrix(NA, nrow = dim(IsoformFractionFiltered)[1], ncol = 19)
row.names(q1_gpr) <- row.names(IsoformFractionFiltered)
colnames(q1_gpr) <- c("RC9_2i_1", "RC9_2i_2", "RC9_2h",
                      "RC9_4h", "RC9_6h", "RC9_8h",
                      "RC9_10h", "RC9_12h", "RC9_14h",
                      "RC9_16h", "RC9_18h", "RC9_20h",
                      "RC9_22h", "RC9_24h", "RC9_26h",
                      "RC9_28h", "RC9_30h",
                      "RC9_32h_1","RC9_32h_2")

q2_gpr <- matrix(NA, nrow = dim(IsoformFractionFiltered)[1], ncol = 19)
row.names(q2_gpr) <- row.names(IsoformFractionFiltered)
colnames(q2_gpr) <- c("RC9_2i_1", "RC9_2i_2", "RC9_2h",
                      "RC9_4h", "RC9_6h", "RC9_8h",
                      "RC9_10h", "RC9_12h", "RC9_14h",
                      "RC9_16h", "RC9_18h", "RC9_20h",
                      "RC9_22h", "RC9_24h", "RC9_26h",
                      "RC9_28h", "RC9_30h",
                      "RC9_32h_1","RC9_32h_2")

X <- data.frame(time=c(0,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,32)) # time-points

gpr_models <- vector(mode="list",length=dim(IsoformFractionFilteredGPR)[1])
names(gpr_models) <- row.names(IsoformFractionFilteredGPR)

### Apply GPR on Isoform Fraction w/o 2iL
for(d in 1:dim(IsoformFractionFilteredGPR)[1]){
  if(d %% 1000 == 0){
    message(d)
  }
  
  if(all(IsoformFractionFiltered[d,c(3:21)]==0)){
    next
  }
  
  Y <- data.frame(reads = IsoformFractionFiltered[d,c(3:21)])
  set.seed(1234)
  gpr <- bgp(X,Y,verb = 0,bprior = "bmznot")
  gpr_models[[d]] <- gpr
  mean <- mean(IsoformFractionFiltered[d,c(3:21)])
  
  IsoformFractionFilteredGPR[d,c(3:21)] <- gpr[["Zp.mean"]]
  q1_gpr[d,] <- gpr[["Zp.q1"]]
  q2_gpr[d,] <- gpr[["Zp.q2"]]
  
}

IsoformFractionFilteredGPR <- IsoformFractionFilteredGPR[complete.cases(IsoformFractionFilteredGPR), ]
q1_gpr <- q1_gpr[complete.cases(q1_gpr), ]
q2_gpr <- q2_gpr[complete.cases(q2_gpr), ]
### Write results to list 
gpr_list_IF <- list(IsoformFractionFilteredGPR,q1_gpr,q2_gpr)
gpr_list_IF_mean <- gpr_list_IF[[1]]
gpr_list_IF_mean <- gpr_list_IF_mean[, c(3:21)]
### Check min and max (w/ and w/o 2iL)
min(gpr_list_IF[[1]])
max(gpr_list_IF[[1]])
min(gpr_list_IF_mean)
max(gpr_list_IF_mean)
# ==> Since min and max are still roughly around 0 and 1 we can assume that a difference of 0.2 corresponds ~ 20% change


## Select isoforms with difference in isoform usage (dIF) over threshold of 20%
gpr_filtered_IF <- data.frame("isoforms" = rownames(gpr_list_IF_mean), 
                              "min" = NA, 
                              "max" = NA, 
                              "change_min_max" = NA)

gpr_filtered_IF$min <- rowMin(gpr_list_IF_mean)
gpr_filtered_IF$max <- rowMax(gpr_list_IF_mean)
gpr_filtered_IF$change_min_max <- gpr_filtered_IF[, 3] - gpr_filtered_IF[, 2] 
#gpr_filtered_IF$fraction_min_max_ <- (gpr_filtered_IF[, 2]/gpr_filtered_IF[, 3])*100
#gpr_unfiltered_IF <- gpr_filtered_IF

### Filter for isoforms with >= 10% difference between at least 2 samples, i.e. compare max and min fraction
gpr_filtered_IF_10 <- gpr_filtered_IF[which(gpr_filtered_IF$change_min_max >= 0.1), ]
### Filter for isoforms with >= 20% difference between at least 2 samples, i.e. compare max and min fraction
gpr_filtered_IF <- gpr_filtered_IF[which(gpr_filtered_IF$change_min_max >= 0.2), ]

### Add gene names for each isoform
gpr_filtered_IF_10$gene_name <- NA
for (i in 1:nrow(gpr_filtered_IF_10)){
  if (gpr_filtered_IF_10$isoforms[i] %in% SwitchList$isoformFeatures$isoform_id){
    gpr_filtered_IF_10$gene_name[i] <- SwitchList$isoformFeatures$gene_name[which(SwitchList$isoformFeatures$isoform_id == gpr_filtered_IF_10$isoforms[i])] 
  }
}
### Add gene names for each isoform
gpr_filtered_IF$gene_name <- NA
for (i in 1:nrow(gpr_filtered_IF)){
  if (gpr_filtered_IF$isoforms[i] %in% SwitchList$isoformFeatures$isoform_id){
    gpr_filtered_IF$gene_name[i] <- SwitchList$isoformFeatures$gene_name[which(SwitchList$isoformFeatures$isoform_id == gpr_filtered_IF$isoforms[i])] 
  }
}


### 
gpr_filtered_IF_mean <- gpr_list_IF_mean[which(row.names(gpr_list_IF_mean) %in% gpr_filtered_IF$isoforms), ]
gpr_filtered_IF_mean <- as.data.frame(gpr_filtered_IF_mean)
gpr_filtered_IF_mean$gene_name <- gpr_filtered_IF$gene_name

gpr_filtered_IF$gene_name <- as.factor(gpr_filtered_IF$gene_name)



# Create histogram of maximum-dIF distribution
## Plots 
hist(gpr_filtered_IF_mean)
hist(gpr_filtered_IF$change_min_max) # Distribution of change

hist_dIF_distribution_unfiltered <- ggplot(gpr_unfiltered_IF, aes(change_min_max)) +
  geom_histogram(bins = 100) +
  labs(title = "Highest dIF between any two time points", 
       x = "dIF", 
       y = "Count of values")

hist_dIF_distribution <- ggplot(gpr_filtered_IF, aes(change_min_max)) +
  geom_histogram(bins = 100) +
  labs(title = "Highest dIF between any two time points", 
       x = "dIF", 
       y = "Count of values")

## Save plots
pdf("~/Desktop/Isoform_analysis_data/hist_dIF_distribution_unfiltered.pdf")
hist_dIF_distribution_unfiltered
dev.off()

pdf("~/Desktop/Isoform_analysis_data/hist_dIF_distribution.pdf")
hist_dIF_distribution
dev.off()


##+ I could set an additional cutoff based on the change in isoform usage? 
##+ Genes appear several times, since we have several isoforms/TCONS per gene... an issue?


# GO-term analysis w/ topGO 
## 
### Create factor vector to indicate important genes 
gene_list <- factor(as.integer(SwitchList$isoformFeatures$gene_name %in% as.character(gpr_filtered_IF$gene_name)))
names(gene_list) <- SwitchList$isoformFeatures$gene_name
names(gene_list)[which(gene_list == 1)]

### Create topGOdata object using the list of genes for Biological Processes and Molecular Function
topGO_BP <- new("topGOdata",
                ontology = "BP",
                description = "Isoform usage",
                allGenes = gene_list,
                annot = annFUN.org,
                mapping = "org.Mm.eg.db",
                ID = "symbol",
                nodeSize = 5)

topGO_MF <- new("topGOdata",
                ontology = "MF",
                description = "Isoform usage",
                allGenes = gene_list,
                annot = annFUN.org,
                mapping = "org.Mm.eg.db",
                ID = "symbol",
                nodeSize = 5)

topGO_CC <- new("topGOdata",
                ontology = "CC", 
                description = "Isoform usage", 
                allGenes = gene_list, 
                annot = annFUN.org, 
                mapping = "org.Mm.eg.db",
                ID = "symbol",
                nodeSize = 5)

### Conduct enrichment analysis by testing the over-representation of GO terms
result_BP <- runTest(object = topGO_BP, algorithm = "weight", statistic = "fisher")
result_MF <- runTest(object = topGO_MF, algorithm = "weight", statistic = "fisher")
result_CC <- runTest(object = topGO_CC, algorithm = "weight", statistic = "fisher")

### Filter for GO terms with at least 5 significant genes & less than 500 annotated genes
table_BP <- GenTable(object = topGO_BP,
                     weight = result_BP,
                     topNodes = 1000)
table_BP <- filter(.data = table_BP, Significant >= 5)
#table_BP <- filter(.data = table_BP, Annotated < 500)
table_BP$weight <- as.numeric(table_BP$weight)
table_BP <- arrange(.data = table_BP, weight)

table_MF <- GenTable(object = topGO_MF,
                     weight = result_MF,
                     topNodes = 900)
table_MF <- filter(.data = table_MF, Significant >= 5)
#table_MF <- filter(.data = table_MF, Annotated < 500)
table_MF$weight <- as.numeric(table_MF$weight)
table_MF <- arrange(.data = table_MF, weight)

table_CC <- GenTable(object = topGO_CC,
                     weight = result_CC,
                     topNodes = 900)
table_CC <- filter(.data = table_CC, Significant >= 5)
#table_MF <- filter(.data = table_MF, Annotated < 500)
table_CC$weight <- as.numeric(table_CC$weight)
table_CC <- arrange(.data = table_CC, weight)

### Prepare plot table w/ top-20 terms, convert p-value column to numeric, take log and concatenate GO ID and Term for better visualization.
### Should I use corrected p-values? 
# table_BP_plot <- table_BP[c(1:20), ]
# table_BP_plot <- mutate(.data = table_BP_plot, minus_log10_pvalue = (-1)*log10(weight))
# table_BP_plot$GO_ID_and_term <- paste(table_BP_plot$Term, table_BP_plot$GO.ID, sep = "_") 

table_BP_plot <- table_BP[c(1:20), ]
table_BP_plot <- mutate(.data = table_BP_plot, minus_log10_pvalue = (-1)*log10(weight))
table_BP_plot$GO_ID_and_term <- paste(table_BP_plot$Term)

# table_MF_plot <- table_MF[c(1:20), ]
# table_MF_plot <- mutate(.data = table_MF_plot, minus_log10_pvalue = (-1)*log10(weight))
# table_MF_plot$GO_ID_and_term <- paste(table_MF_plot$Term, table_MF_plot$GO.ID, sep = "_")

table_MF_plot <- table_MF[c(1:20), ]
table_MF_plot <- mutate(.data = table_MF_plot, minus_log10_pvalue = (-1)*log10(weight))
table_MF_plot$GO_ID_and_term <- paste(table_MF_plot$Term)

# table_CC_plot <- table_CC[c(1:20), ]
# table_CC_plot <- mutate(.data = table_CC_plot, minus_log10_pvalue = (-1)*log10(weight))
# table_CC_plot$GO_ID_and_term <- paste(table_CC_plot$Term, table_CC_plot$GO.ID, sep = "_")

table_CC_plot <- table_CC[c(1:20), ]
table_CC_plot <- mutate(.data = table_CC_plot, minus_log10_pvalue = (-1)*log10(weight))
table_CC_plot$GO_ID_and_term <- paste(table_CC_plot$Term)

### Plot 
# ggplot(table_BP_plot, aes(x = minus_log10_pvalue, y = reorder(x = GO_ID_and_term, X = minus_log10_pvalue))) + 
#   geom_bar(stat = "identity") + 
#   labs(title = "BPs of top genes w/ at least 5 significant proteins ",
#        subtitle = "341 genes; filtered by 'significant' isoform usage",
#        x = "-log10(p-value)",
#        y = "GO Term") +
#   theme(plot.title = element_text(hjust = 1),
#         plot.subtitle = element_text(hjust = 1), 
#         text = element_text(size = 18)) +
#   geom_vline(xintercept = -log10(0.05), linetype = "dotted")

g_BP <- ggplot(table_BP_plot, aes(x = minus_log10_pvalue, y = reorder(x = GO_ID_and_term, X = minus_log10_pvalue))) + 
  geom_bar(stat = "identity") + 
  labs(title = "GO enrichment: Biological process",
       x = "-log10(p-value)",
       y = "GO Term") +
  theme(plot.title = element_text(hjust = 1, size = 16),
        plot.subtitle = element_text(hjust = 1), 
        text = element_text(size = 18),
        axis.title = element_text(size = 16)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted")

tiff("~/Desktop/Isoform_analysis_data/GO_BP.tif", width = 23.35, height = 17.35, units = "cm", res = 300,
    compression = "lzw")
g_BP
dev.off()


# ggplot(table_MF_plot, aes(x = minus_log10_pvalue, y = reorder(x = GO_ID_and_term, X = minus_log10_pvalue))) + 
#   geom_bar(stat = "identity") + 
#   labs(title = "MFs of top genes w/ at least 5 significant proteins ",
#        subtitle = "341 genes; filtered by 'significant' isoform usage",
#        x = "-log10(p-value)",
#        y = "GO Term") +
#   theme(plot.title = element_text(hjust = 1),
#         plot.subtitle = element_text(hjust = 1)) +
#   geom_vline(xintercept = -log10(0.05), linetype = "dotted")

g_MF <- ggplot(table_MF_plot, aes(x = minus_log10_pvalue, y = reorder(x = GO_ID_and_term, X = minus_log10_pvalue))) + 
  geom_bar(stat = "identity") + 
  labs(title = "GO enrichment: Molecular function",
       x = "-log10(p-value)",
       y = "GO Term") +
  theme(plot.title = element_text(hjust = 1, size = 16),
        plot.subtitle = element_text(hjust = 1), 
        text = element_text(size = 18),
        axis.title = element_text(size = 16)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted")

tiff("~/Desktop/Isoform_analysis_data/GO_MF.tif", width = 23.35, height = 17.35, units = "cm", res = 300,
     compression = "lzw")
g_MF
dev.off()


# ggplot(table_CC_plot, aes(x = minus_log10_pvalue, y = reorder(x = GO_ID_and_term, X = minus_log10_pvalue))) + 
#   geom_bar(stat = "identity") + 
#   labs(title = "CCs of top genes w/ at least 5 significant proteins ",
#        subtitle = "341 genes; filtered by 'significant' isoform usage",
#        x = "-log10(p-value)",
#        y = "GO Term") +
#   theme(plot.title = element_text(hjust = 1),
#         plot.subtitle = element_text(hjust = 1)) +
#   geom_vline(xintercept = -log10(0.05), linetype = "dotted")

g_CC <- ggplot(table_CC_plot, aes(x = minus_log10_pvalue, y = reorder(x = GO_ID_and_term, X = minus_log10_pvalue))) + 
  geom_bar(stat = "identity") + 
  labs(title = "GO enrichment: Cellular Component",
       x = "-log10(p-value)",
       y = "GO Term") +
  theme(plot.title = element_text(hjust = 1, size = 16),
        plot.subtitle = element_text(hjust = 1), 
        text = element_text(size = 18),
        axis.title = element_text(size = 16)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted")

tiff("~/Desktop/Isoform_analysis_data/GO_CC.tif", width = 23.35, height = 17.35, units = "cm", res = 300,
     compression = "lzw")
g_CC
dev.off()



### Plot subgraph
showSigOfNodes(GOdata = topGO_BP, score(result_BP), firstSigNodes = 5, useInfo = "all")
showSigOfNodes(GOdata = topGO_MF, score(result_MF), firstSigNodes = 5, useInfo = "all") 
showSigOfNodes(GOdata = topGO_CC, score(result_CC), firstSigNodes = 5, useInfo = "all") 

#printGenes(object = topGO_BP, "GO:0035023", chip = "Mus.musculus")

################################################################################################################

# Plot
## Prepare data tables
FPKM_filtered <- SwitchList$isoformRepExpression[which(SwitchList$isoformRepExpression$isoform_id %in% gpr_filtered_IF$isoforms), ]
View(FPKM_filtered)

IF_filtered <- SwitchList$isoformRepIF[which(SwitchList$isoformRepIF$isoform_id %in% gpr_filtered_IF$isoforms), ]
View(IF_filtered)

Isoform_features_filtered <- SwitchList$isoformFeatures[which(SwitchList$isoformFeatures$isoform_id %in% gpr_filtered_IF$isoforms), ]
View(Gene_expression_filtered)

#Exon_filtered <- SwitchList$exons???

## Transcript plot

## FPKMxIF vs Condition (time)

## IF vs Condition (time)

## Gene expression vs Condition (time)


################################################################################################################

# Comparison to 2i vs 32h 
## Match genes that showed significant isoform switch in 2i vs 32h
genes_2ivs32h <- read.csv("/Users/jonathan/Documents/Biologie_MSc/WS_21:22/Project_module_Beyer/Data/2i_vs_32h/Results/Tables/2i_vs_32h_top_switches.csv")
## Count genes in both 2i vs 32h and complete time course
length(genes_2ivs32h$gene_name)
length(gpr_filtered_IF$gene_name)
## Add matches to time course results
gpr_filtered_IF$in_2i_vs_32h <- FALSE
for (i in 1:nrow(gpr_filtered_IF)){
  if (gpr_filtered_IF$gene_name[i] %in% genes_2ivs32h$gene_name){
    gpr_filtered_IF$in_2i_vs_32h[i] <- TRUE 
  }
}
## Count TRUE and FALSE
sum(gpr_filtered_IF$in_2i_vs_32h == TRUE)
sum(gpr_filtered_IF$in_2i_vs_32h == FALSE)

################################################################################################################

# Nonsense Mediated Decay (NMD) controlling ESC development?
## Check transcript length vs expression level to infer whether NMD acts as switch for fast differentiation control
### Add length for transcripts/isoforms of interest (+ add gene position in chromosome)
gpr_filtered_IF$length <- NA
gpr_filtered_IF$locus <- NA

for (i in 1:nrow(gpr_filtered_IF)){
  if (gpr_filtered_IF$isoforms[i] %in% Isoform_features_filtered$isoform_id){
    gpr_filtered_IF$length[i] <- Isoform_features_filtered$length[which(Isoform_features_filtered$isoform_id %in% gpr_filtered_IF$isoforms[i])]
  }
}

for (i in 1:nrow(gpr_filtered_IF)){
  if (gpr_filtered_IF$isoforms[i] %in% Isoform_features_filtered$isoform_id){
    gpr_filtered_IF$locus[i] <- Isoform_features_filtered$locus[which(Isoform_features_filtered$isoform_id %in% gpr_filtered_IF$isoforms[i])]
  }
}
### Check lengths
min(Isoform_features_filtered$length)
max(Isoform_features_filtered$length)
min(gpr_filtered_IF$length)
max(gpr_filtered_IF$length)
### Plot lengths 
TCON_length <- ggplot(gpr_filtered_IF, aes(length)) +
  geom_histogram(bins = 100) + 
  labs(title = "TCONS length distribution",
       x = "Length", 
       y = "Count")


## Calculate delta and fraction of isoform lengths
isoform_length <- gpr_filtered_IF[, c(1, 5:7)]

isoform_length$dLength <- NA
isoform_length$qLength <- NA




for (i in 1:nrow(isoform_length)){
  if (isoform_length$gene_name[i] %in% isoform_length$gene_name){
    isoform_length$dLength[i] <- isoform_length$length[i]-min(isoform_length$length[isoform_length$gene_name[i] == isoform_length$gene_name])
    isoform_length$qLength[i] <- isoform_length$dLength[i]/isoform_length$length
  }
}

### Filter, so that we keep only non-unique entries, i.e. those w/ at least 2 isoforms
isoform_length_filtered <- isoform_length[isoform_length$gene_name %in% isoform_length$gene_name[duplicated(isoform_length$gene_name)], ]
### Remove isoforms that are not annotated, i.e. have no gene-name associated
isoform_length_filtered$gene_name[isoform_length_filtered$gene_name == "-"] <- NA
isoform_length_filtered <- na.omit(isoform_length_filtered)


## Compare to proteomics data 
### Load proteomics FC 
gp_ptc_tibble_wide <- read.csv("~/Desktop/Isoform_analysis_data/gp_ptc_tibble_wide.csv")
### Subset proteomics FC by Isoforms of Interest 
proteomics_subset <- gp_ptc_tibble_wide[which(gp_ptc_tibble_wide$gene %in% gpr_filtered_IF$gene_name), ]
### Subset proteomics FC by isoforms w/ change in length
proteomics_subset <- proteomics_subset[which(proteomics_subset$gene %in% isoform_length_filtered$gene_name), ]

t_proteomics_subset <- t(proteomics_subset)
colnames(t_proteomics_subset) <- t_proteomics_subset[1,]


################################################################################################################

## Extract exon information
exons <- annoGR2DF(SwitchList$exons)
### Subset exon-info by Isoforms of Interest (instead of filtering by TCONS I could also filter by XLOC, i.e. gene_id to get all detected isoforms??)
exons_subset <- exons[which(exons$isoform_id %in% gpr_filtered_IF$isoforms), ]
### Add gene names
exons_subset$gene_name <- NA
for (i in 1:nrow(exons_subset)){
  if (exons_subset$isoform_id[i] %in% gpr_filtered_IF$isoforms){
    exons_subset$gene_name[i] <- gpr_filtered_IF$gene_name[which(gpr_filtered_IF$isoforms %in% exons_subset$isoform_id[i])]
  }
}
### Plot Exon widths
exon_width <- ggplot(exons_subset, aes(width)) +
  geom_histogram(bins = 100) + 
  labs(title = "Exon width distribution",
       x = "Exon Width", 
       y = "Count")


################################################################################################################


# Load 2i vs 32h SwitchList 
SwitchList_2ivs32h <- readRDS("~/Desktop/Isoform_analysis_data/2i_vs_32h/SwitchList_Part2.RDS")
## Extract TCON START-STOP, filtered by Isoforms of Interest



extractSwitchSummary(SwitchList_2ivs32h)
names(SwitchList_2ivs32h)
View(SwitchList_2ivs32h$ntSequence)








################################################################################################################

# Subset a whole SwitchList
## Idea: Subset SwitchList by Isoforms of Interest to take advantage of SwitchAnalysisPartII, i.e. annotation w/ external tools + visualization
test <- subsetSwitchAnalyzeRlist(
  SwitchList,
  SwitchList$isoformFeatures$gene_name == "Etv5" # can I just provide a string with all Isoforms of Interest?
)

################################################################################################################





# Export tables
write.csv(gpr_filtered_IF, "~/Downloads/gpr_filtered_IF.csv")

# Export data files  
## Save image
save.image("Isoform_detection_cufflinks.R")
## Save objects
saveRDS(gpr_list_IF, "gpr_list_IF.rds")
saveRDS(IsoformCountFilteredGPR,"IsoformCountFilteredGPR.rds")
saveRDS(gpr_filtered_IF, "~/Desktop/Isoform_analysis_data/gpr_filtered_IF.rds")
saveRDS(FPKM_filtered, "~/Desktop/Isoform_analysis_data/FPKM_filtered.rds")
saveRDS(IF_filtered, "~/Desktop/Isoform_analysis_data/IF_filtered.rds")
saveRDS(Isoform_features_filtered, "~/Desktop/Isoform_analysis_data/Isoform_features_filtered.rds")






















##################################################################################################
geneCountMatrix <- extractGeneExpression(
  SwitchList,
  extractCounts = TRUE
)


geneExpresionMatrix <- extractGeneExpression(
  SwitchList,
  extractCounts = FALSE
)


