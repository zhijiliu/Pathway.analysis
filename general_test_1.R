# union method for pathway anlaysis inlcuding genetic, proteomic, and metabolomic metrics
wd = '/Users/zhijiliu/Google Drive/Mac Rstudio/Pathway analysis/PE/Promedex'
setwd(wd)

# ----- functions -----
# source("https://bioconductor.org/biocLite.R")

library(fgsea)
library(xlsx)
library(clusterProfiler)
library(DOSE)
library(ggplot2)
library(grid)
library(limma)
library(plotly)
library(igraph)
library(randomForest)
library(impute)
library(ROCR)

# remove white space
trim.ws.custom <- function (x) gsub("^\\s+|\\s+$", "", x)
# convert nan to NA
is.nan.data.frame <- function(x){do.call(cbind, lapply(x, is.nan))}
# insert row into dataframe
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  return(existingDF)
}
# internal for GSEA plot
gseaScores <- function(geneList, geneSet, exponent=1, fortify=FALSE) {
  ###################################################################
  ##    geneList                                                   ##
  ##                                                               ##
  ## 1. Rank order the N genes in D to form L = { g_1, ... , g_N}  ##
  ##    according to the correlation, r(g_j)=r_j,                  ##
  ##    of their expression profiles with C.                       ##
  ##                                                               ##
  ###################################################################
  
  ###################################################################
  ##    exponent                                                   ##
  ##                                                               ##
  ## An exponent p to control the weight of the step.              ##
  ##   When p = 0, Enrichment Score ( ES(S) ) reduces to           ##
  ##   the standard Kolmogorov-Smirnov statistic.                  ##
  ##   When p = 1, we are weighting the genes in S                 ##
  ##   by their correlation with C normalized                      ##
  ##   by the sum of the correlations over all of the genes in S.  ##
  ##                                                               ##
  ###################################################################
  
  ## genes defined in geneSet should appear in geneList.
  ## geneSet <- intersect(geneSet, names(geneList))
  
  N <- length(geneList)
  Nh <- length(geneSet)
  
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet ## logical
  
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit/NR)
  
  Pmiss[!hits] <-  1/(N-Nh)
  Pmiss <- cumsum(Pmiss)
  
  runningES <- Phit - Pmiss
  
  ## ES is the maximum deviation from zero of Phit-Pmiss
  max.ES <- max(runningES)
  min.ES <- min(runningES)
  if( abs(max.ES) > abs(min.ES) ) {
    ES <- max.ES
  } else {
    ES <- min.ES
  }
  
  df <- data.frame(x=seq_along(runningES),
                   runningScore=runningES,
                   position=as.integer(hits)
  )
  
  if(fortify==TRUE) {
    return(df)
  }
  
  df$gene = names(geneList)
  res <- list(ES=ES, runningES = df)
  return(res)
}
# internal for GSEA plot
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin=0
  df$ymax=0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  
  return(df)
}
# SAM
sam.analysis <- function(analysis_option, input.data, clinic.info){
  if (analysis_option == 'total'){
    sam_input <- input.data
  } else if (analysis_option == 'early'){
    sam_input <- input.data[,c('Compound',clinic.info$ID[clinic.info$GA < 34])]
  } else if (analysis_option == 'late'){
    sam_input <- input.data[,c('Compound',clinic.info$ID[clinic.info$GA >= 34])]
  } else if (analysis_option == 'severe'){
    sam_input <- input.data[,c('Compound',
                               clinic.info$ID[clinic.info$Status2 == 'severe' | clinic.info$Status2 == 'normal'])]
  } else if (analysis_option == 'mild'){
    sam_input <- input.data[,c('Compound',
                               clinic.info$ID[clinic.info$Status2 == 'mild' | clinic.info$Status2 == 'normal'])]
  }
  # filter the metabolite total
  colnum <- dim(sam_input)[2]
  row_index <- apply(sam_input[,2:colnum], 1, function(x) {
    sum(is.na(x)) > (colnum/2)})
  sam_input <- sam_input[!row_index,]
  # NA filtering case
  row_index <- apply(sam_input[,grep('E',colnames(sam_input))], 1, function(x) {
    sum(is.na(x)) > (length(grep('E',colnames(sam_input)))/2)})
  sam_input <- sam_input[!row_index,]
  # NA filtering control
  row_index <- apply(sam_input[,grep('P',colnames(sam_input))], 1, function(x) {
    sum(is.na(x)) > (length(grep('P',colnames(sam_input)))/2)})
  sam_input <- sam_input[!row_index,]
  rownames(sam_input) <- NULL
  return_data <- sam_input
  # labeling
  sam_input <- as.matrix(sam_input[,2:colnum])
  colname <- substring(colnames(sam_input),1,1)
  label <- c()
  label[which(colname == 'E')] <- 2
  label[which(colname == 'P')] <- 1
  data <- list(x = sam_input, y = label,logged2 = F)
  # use t-statistic and 500 permutations for FDR estimation
  samr.obj<-samr(data = data,
                 resp.type="Two class unpaired",
                 center.arrays = F,
                 assay.type = 'array',
                 testStatistic = 'standard',
                 nperms=500)
  # Returns a table of the FDR and upper and lower cutpoints for various values of delta, for a SAM analysis.
  delta.table <- samr.compute.delta.table(samr.obj,
                                          min.foldchange=0)
  # Computes significant genes table, starting with samr object "samr.obj" and delta.table "delta.table"
  siggenes.table <- samr.compute.siggenes.table(samr.obj,
                                                del=0,
                                                data,
                                                delta.table,
                                                all.genes=TRUE,
                                                compute.localfdr = T)
  tmp.up <- siggenes.table$genes.up
  tmp.down <- siggenes.table$genes.lo
  tmp.result <- rbind(tmp.up,tmp.down)
  tmp.result <- data.frame(tmp.result,stringsAsFactors = F)
  tmp.result$Gene.Name <- return_data$Compound[match(as.numeric(tmp.result$Row)-1,rownames(return_data))]
  final.result <- tmp.result[,c(3,7,8,9)]
  colnames(final.result) <- c('Compound','FC','q.value(%)','local.FDR(%)')
  final.result$FC <- as.numeric(final.result$FC)
  final.result$`q.value(%)` <- as.numeric(final.result$`q.value(%)`)
  final.result$`local.FDR(%)` <- as.numeric(final.result$`local.FDR(%)`)
  return_obj <- list(final.result,return_data)
  return(return_obj)
}
# geo preprocess
preprocess <- function(GSE_data,expression_mat,label){
  # Get expression matrix from the GSE data and perform log2 transformation
  # Args:
  #   GSE_data: GSE expression data
  #   expression_mat: expression matrix converted from GSE data
  #   label: define case and control
  # Returns:
  #   matrix calculated from toptable
  
  qx <- as.numeric(quantile(expression_mat, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) {
    expression_mat[which(expression_mat <= 0)] <- NaN
    exprs(GSE_data) <- log2(expression_mat)
    message('log2 transformed')
  }
  expression_mat=exprs(GSE_data)
  
  
  # # calculate variance of the dataset 
  # expression_var = apply(expression_mat,1,function(x){
  #   return( var(x[which(label=='G0')], na.rm = T)/median(x[which(label=='G0')],na.rm = T) + var(x[which(label=='G1')],na.rm = T)/median(x[which(label=='G1')],na.rm = T) )
  # })
  # message('relative variantion')
  # expression_var=data.frame(expression_var)
  # rownames(expression_var)=rownames(expression_mat)
  # colnames(expression_var)='variation'
  # message('data variation calculated')
  
  # set up the data and proceed with analysis
  label <- as.factor(label)
  GSE_data$description <- label
  design <- model.matrix(~ description + 0, GSE_data)
  message('design matrix created')
  colnames(design) <- levels(label)
  fit <- lmFit(GSE_data, design)
  cont.matrix <- makeContrasts(G1-G0, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.01) # only choose the first 1%
  tT <- topTable(fit2, adjust.method="BH", sort.by="M", number=dim(expression_mat)[1],confint = T)
  # tT=merge(tT,expression_var, by="row.names")
  # message('data analyzed and variation combined')
  return(tT)
}
# geo diff analysis
selection <- function(GSE_data,tT,platform_mat){
  # select the necessary calculated results;
  # delete non-symbol probe;
  # perform multiple probe combination; 
  # 
  # Args:
  #   tT: results from previous calculation
  #   GSE_data: original expression dataset
  #   platform_mat: platform information matrix
  # Returns:
  #   data matrix with unique gene symbol and required results 
  
  
  # replace with original platform annotation
  # in most case they will be the same, just to prevent the case that GSE data has error
  tT <- tT[setdiff(colnames(tT), setdiff(fvarLabels(GSE_data), c("ID","Symbol")))]
  tT <- merge(tT, platform_mat, by="ID")
  # tT.select is the necessary data
  tT.select <- subset(tT, select=c("ID","entrez","Symbol","Platform.Symbol","P.Value","logFC","adj.P.Val","B","CI.L","CI.R")) 
  # data type conversion
  tT.select$ID <- as.character(tT.select$ID)
  tT.select$Symbol <- as.character(tT.select$Symbol)
  tT.select$Platform.Symbol <- as.character(tT.select$Platform.Symbol)
  message('platform gene symbol replaced')
  
  # delete the rows that do not have a gene name
  tT.select <- tT.select[tT.select$Platform.Symbol != "", ]
  tT.select <- tT.select[!is.na(tT.select$entrez), ]
  tT.select <- tT.select[!is.na(tT.select$logFC), ]
  
  # # official symbol matching; if no match, use existing symbol
  # tT.select$Gene.platform.symbol <- sapply(strsplit(tT.select$Gene.platform.symbol,' /// '),"[",1)
  # # tmpindex <- match(tT.select$Gene.platform.symbol,aliasSymbol$alias_symbol)
  # # tmpsymbol <- aliasSymbol$symbol[tmpindex]
  # # tT.select$Gene.platform.symbol[!is.na(tmpsymbol)] <- tmpsymbol[!is.na(tmpsymbol)]
  # # rm(tmpindex,tmpsymbol)
  # # message('official gene symbol converted')
  
  # record the gene symbol for redundant elimination
  tem.symbol = tT.select$Platform.Symbol
  rownames(tT.select) = NULL
  
  # delete the same gene: if there are more than one value of a gene
  # first check the FC, keep the FC>2 genes
  # then check the p-value, keep the smallest p-value gene
  unique.symbol = unique(tem.symbol)
  
  tT.select = sapply(unique.symbol, function(x){
    same.index = which(tem.symbol == x)
    if(length(same.index) == 1) {
      return(tT.select[same.index, ])
    } 
    else{
      fc = tT.select[same.index, 6]
      pv = tT.select[same.index, 5]
      if(sum(abs(fc)>=1) == 0) {
        g.index = which.min(pv)
      } 
      else {
        g.index = which(pv == min(pv[abs(fc)>=1]))
      }
      return(tT.select[same.index[g.index], ])
    }
  })
  tT.select = t(tT.select)
  tT.select = as.data.frame(tT.select)
  rownames(tT.select) = NULL
  message('unique symbol gene selected')
  
  # data type convertion
  tT.select$ID <- as.character(tT.select$ID)
  tT.select$Symbol <- as.character(tT.select$Symbol)
  tT.select$Platform.Symbol <- as.character(tT.select$Platform.Symbol)
  tT.select$P.Value <- as.numeric(tT.select$P.Value)
  tT.select$logFC <- as.numeric(tT.select$logFC)
  tT.select$entrez <- as.character(tT.select$entrez)
  # tT.select$variation <- as.numeric(tT.select$variation)
  tT.select$adj.P.Val <- as.numeric(tT.select$adj.P.Val)
  tT.select$B <- as.numeric(tT.select$B)
  tT.select$CI.L <- as.numeric(tT.select$CI.L)
  tT.select$CI.R <- as.numeric(tT.select$CI.R)
  message('statistical result data type converted')
  
  return(tT.select)
}

# ------ 1. transcription story ------
# platform 
gpl <- getGEO('GPL10558')
gpl <- data.frame(attr(dataTable(gpl), "table"))
gpl <- gpl[,c(1,10,13)]
colnames(gpl) <- c('ID','entrez','Platform.Symbol')
# gse download
gse <- getGEO('GSE44711')
gse <- gse[[1]]
gse.clinic <- pData(gse)
# data processing
ex <- exprs(gse)
gse.label <- paste("G",c(rep(1,8),rep(0,8)),sep = '')
tT <- preprocess(gse,ex,gse.label)
tT.select <- selection(gse,tT,gpl)
# perform geneset enrichment analysis
genelist <- tT.select$logFC
names(genelist) <- tT.select$entrez
genelist <- genelist[order(-genelist)]
kk2 <- gseKEGG(geneList     = genelist,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 10,
               pvalueCutoff = 1,
               verbose      = FALSE)
gseaplot(kk2, geneSetID = "hsa00100")

# # ---- debug ----
# tT.test <- tT[setdiff(colnames(tT), setdiff(fvarLabels(gse), c("ID","Symbol")))]
# tT.test2 <- merge(tT.test, gpl, by="ID")
# tT.test2 <- subset(tT.test2, select=c("ID","entrez","Symbol","Platform.Symbol","P.Value","logFC","adj.P.Val","B","CI.L","CI.R"))

# ----- gene pathway plot -----
# ----- gene correlation network based on GSE44711
library(GEOmetadb)
library(GEOquery)
GSE44711 <- getGEO('gse44711',destdir = ".")
GSE44711 <- GSE44711[[1]]
GSE44711 <- exprs(GSE44711)
GPL <- getGEO('GPL10558')
gene_probe_tbl <- GPL@dataTable@table

# ----- GSEA targeted pathways -----
path_gene_tbl <- data.frame(gene = NULL,path = NULL)

tmp_pathway <- read.xlsx2(file = '/Users/zhijiliu/Downloads/GSEA GENES.xlsx',sheetIndex = 1)
tmp_gene <- as.character(unique(tmp_pathway$PROBE))
tmp_gene <- tmp_gene[which(tmp_gene!='')]
tmp_pathway_2 <- data.frame(gene = tmp_gene,
                            path = 'Amino acid metabolism',
                            stringsAsFactors = F)
tmp_pathway_2$sig <- as.character(tmp_pathway$CORE.ENRICHMENT[match(tmp_pathway_2$gene,tmp_pathway$PROBE)])

path_gene_tbl <- rbind(path_gene_tbl,tmp_pathway_2)

tmp_pathway <- read.xlsx2(file = '/Users/zhijiliu/Downloads/GSEA GENES.xlsx',sheetIndex = 2)
tmp_gene <- as.character(unique(tmp_pathway$PROBE))
tmp_gene <- tmp_gene[which(tmp_gene!='')]
tmp_pathway_2 <- data.frame(gene = tmp_gene,
                            path = 'Sphingolipid metabolism',
                            stringsAsFactors = F)
tmp_pathway_2$sig <- as.character(tmp_pathway$CORE.ENRICHMENT[match(tmp_pathway_2$gene,as.character(tmp_pathway$PROBE))])

path_gene_tbl <- rbind(path_gene_tbl,tmp_pathway_2)

tmp_pathway <- read.xlsx2(file = '/Users/zhijiliu/Downloads/GSEA GENES.xlsx',sheetIndex = 3)
tmp_gene <- as.character(unique(tmp_pathway$PROBE))
tmp_gene <- tmp_gene[which(tmp_gene!='')]
tmp_pathway_2 <- data.frame(gene = tmp_gene,
                            path = 'Steroid metabolism',
                            stringsAsFactors = F)
tmp_pathway_2$sig <- as.character(tmp_pathway$CORE.ENRICHMENT[match(tmp_pathway_2$gene,tmp_pathway$PROBE)])

path_gene_tbl <- rbind(path_gene_tbl,tmp_pathway_2)

tmp_pathway <- read.xlsx2(file = '/Users/zhijiliu/Downloads/GSEA GENES.xlsx',sheetIndex = 4)
tmp_gene <- as.character(unique(tmp_pathway$PROBE))
tmp_gene <- tmp_gene[which(tmp_gene!='')]
tmp_pathway_2 <- data.frame(gene = tmp_gene,
                            path = 'Fatty acid metabolism',
                            stringsAsFactors = F)
tmp_pathway_2$sig <- as.character(tmp_pathway$CORE.ENRICHMENT[match(tmp_pathway_2$gene,tmp_pathway$PROBE)])

path_gene_tbl <- rbind(path_gene_tbl,tmp_pathway_2)

tmp_pathway <- read.xlsx2(file = '/Users/zhijiliu/Downloads/GSEA GENES.xlsx',sheetIndex = 5)
tmp_gene <- as.character(unique(tmp_pathway$PROBE))
tmp_gene <- tmp_gene[which(tmp_gene!='')]
tmp_pathway_2 <- data.frame(gene = tmp_gene,
                            path = 'Glycerophospolipid metabolism',
                            stringsAsFactors = F)
tmp_pathway_2$sig <- as.character(tmp_pathway$CORE.ENRICHMENT[match(tmp_pathway_2$gene,tmp_pathway$PROBE)])

path_gene_tbl <- rbind(path_gene_tbl,tmp_pathway_2)

# ----- gene to pathway mapping -----
tmp_probe <- as.character(gene_probe_tbl$ID[match(path_gene_tbl$gene,as.character(gene_probe_tbl$ILMN_Gene))])

input_mat <- GSE44711[tmp_probe,]
# stats
log_fc <- apply(input_mat,1,function(x){
  tmp.case <- sum(x[1:8],na.rm = T)
  tmp.control <- sum(x[9:16],na.rm = T)
  if (tmp.case > tmp.control) {return('up')
  } else {return('dn')}
})
p_val <- apply(input_mat,1,function(x){
  tmp.case <- x[1:8]
  tmp.control <- x[9:16]
  tmp <- t.test(tmp.case,tmp.control)$p.value
  return(tmp)
})
# correlation matrix
cor_mat <- cor(t(input_mat),use = 'everything',method = 'pearson')

# distance matrix
dis_mat <- abs(cor_mat)
dis_mat[is.na(dis_mat)] <- 0
dis_mat <- apply(dis_mat,2,function(x){
  # x[x<0.3 & x>0.2] <- 0.02
  # x[x<0.2 & x>0.1] <- 0.005
  x[x<0.1] <- 0.1
  return(1/x)
})
diag(dis_mat) <- 0

# coordinates
fit <- cmdscale(dis_mat, eig = TRUE, k = 2)
fit_mds <- data.frame(x = fit$points[,1], y = fit$points[,2],stringsAsFactors = F)

# adjacency matrix
adj_mat <- abs(cor_mat)
adj_mat[is.na(adj_mat)] <- 0
adj_mat <- apply(adj_mat,2,function(x){
  x[abs(x)<0.9] <- 0
  x[abs(x)>=0.9] <- 1
  return(x)
})
diag(adj_mat) <- 1



# use pca to get gene location
pe_pca <- prcomp(input_mat,
                 center = T,
                 scale. = T)
fit_pca <- data.frame(pe_pca$x)

# ----- based on pca -----
# get layout
layout_2d <- fit_pca[,c(1,2)]
layout_2d$probe <- rownames(input_mat)
layout_2d$gene <- as.character(gene_probe_tbl$ILMN_Gene[match(layout_2d$probe,gene_probe_tbl$ID)])
layout_2d$sig <- path_gene_tbl$sig
layout_2d$path <- path_gene_tbl$path

tmp_path <- unique(layout_2d$path)

# reset 
layout_2d$PC1 <- fit_pca$PC1
layout_2d$PC2 <- fit_pca$PC2
# AA 
layout_2d$PC1[grep(tmp_path[1],layout_2d$path)] <- layout_2d$PC1[grep(tmp_path[1],layout_2d$path)] - 5
layout_2d$PC2[grep(tmp_path[1],layout_2d$path)] <- layout_2d$PC2[grep(tmp_path[1],layout_2d$path)] + 5
# Sphin
layout_2d$PC1[grep(tmp_path[2],layout_2d$path)] <- layout_2d$PC1[grep(tmp_path[2],layout_2d$path)] + 5
layout_2d$PC2[grep(tmp_path[2],layout_2d$path)] <- layout_2d$PC2[grep(tmp_path[2],layout_2d$path)] + 5
# Ster
layout_2d$PC1[grep(tmp_path[3],layout_2d$path)] <- layout_2d$PC1[grep(tmp_path[3],layout_2d$path)] + 0
layout_2d$PC2[grep(tmp_path[3],layout_2d$path)] <- layout_2d$PC2[grep(tmp_path[3],layout_2d$path)] + 0
# Fatty
layout_2d$PC1[grep(tmp_path[4],layout_2d$path)] <- layout_2d$PC1[grep(tmp_path[4],layout_2d$path)] - 5
layout_2d$PC2[grep(tmp_path[4],layout_2d$path)] <- layout_2d$PC2[grep(tmp_path[4],layout_2d$path)] - 5
# Gly
layout_2d$PC1[grep(tmp_path[5],layout_2d$path)] <- layout_2d$PC1[grep(tmp_path[5],layout_2d$path)] + 5
layout_2d$PC2[grep(tmp_path[5],layout_2d$path)] <- layout_2d$PC2[grep(tmp_path[5],layout_2d$path)] - 5

# edge setting
g <- graph_from_adjacency_matrix(adjmatrix = adj_mat, mode = 'undirected',diag = T)
# Create Vertices and Edges
vs <- V(g)
es <- as.data.frame(get.edgelist(g))
# length
Nv <- length(vs)
Ne <- length(es[1]$V1)
# axis
Xn <- layout_2d[,1]
Yn <- layout_2d[,2]
# get coordinates
from <- es$V1
to <- es$V2
x1 = Xn[from]
x2 = Xn[to]
y1 = Yn[from]
y2 = Yn[to]

edge_data_frame <- data.frame(x_from = x1,
                              x_to = x2,
                              y_from = y1,
                              y_to = y2)

# ----- based on distance matrix and mds -----
# get layout
layout_2d <- fit_mds[,c(1,2)]
colnames(layout_2d)[1:2] <- c('PC1','PC2')
layout_2d$probe <- rownames(input_mat)
layout_2d$gene <- as.character(gene_probe_tbl$ILMN_Gene[match(layout_2d$probe,gene_probe_tbl$ID)])
layout_2d$sig <- path_gene_tbl$sig
layout_2d$path <- path_gene_tbl$path
identical(layout_2d$probe,names(log_fc))
identical(layout_2d$probe,names(p_val))
layout_2d$p_val <- p_val
layout_2d$FC <- log_fc

tmp_path <- unique(layout_2d$path)

# reset 
layout_2d$PC1 <- fit_mds$x
layout_2d$PC2 <- fit_mds$y
# AA 
layout_2d$PC1[grep(tmp_path[1],layout_2d$path)] <- layout_2d$PC1[grep(tmp_path[1],layout_2d$path)] - 9.5
layout_2d$PC2[grep(tmp_path[1],layout_2d$path)] <- layout_2d$PC2[grep(tmp_path[1],layout_2d$path)] + 6.5
# Sphin
layout_2d$PC1[grep(tmp_path[2],layout_2d$path)] <- layout_2d$PC1[grep(tmp_path[2],layout_2d$path)] + 11.5
layout_2d$PC2[grep(tmp_path[2],layout_2d$path)] <- layout_2d$PC2[grep(tmp_path[2],layout_2d$path)] + 6.5
# Ster
layout_2d$PC1[grep(tmp_path[3],layout_2d$path)] <- layout_2d$PC1[grep(tmp_path[3],layout_2d$path)] + 2
layout_2d$PC2[grep(tmp_path[3],layout_2d$path)] <- layout_2d$PC2[grep(tmp_path[3],layout_2d$path)] + 3
# Fatty
layout_2d$PC1[grep(tmp_path[4],layout_2d$path)] <- layout_2d$PC1[grep(tmp_path[4],layout_2d$path)] - 6.5
layout_2d$PC2[grep(tmp_path[4],layout_2d$path)] <- layout_2d$PC2[grep(tmp_path[4],layout_2d$path)] - 6.5
# Gly
layout_2d$PC1[grep(tmp_path[5],layout_2d$path)] <- layout_2d$PC1[grep(tmp_path[5],layout_2d$path)] + 6.5
layout_2d$PC2[grep(tmp_path[5],layout_2d$path)] <- layout_2d$PC2[grep(tmp_path[5],layout_2d$path)] - 7.5

# edge setting
g <- graph_from_adjacency_matrix(adjmatrix = adj_mat, mode = 'undirected',diag = T)
# Create Vertices and Edges
vs <- V(g)
es <- as.data.frame(get.edgelist(g))
# length
Nv <- length(vs)
Ne <- length(es[1]$V1)
# axis
Xn <- layout_2d[,1]
Yn <- layout_2d[,2]
# get coordinates
from <- es$V1
to <- es$V2
x1 = Xn[from]
x2 = Xn[to]
y1 = Yn[from]
y2 = Yn[to]

edge_data_frame <- data.frame(x_from = x1,
                              x_to = x2,
                              y_from = y1,
                              y_to = y2)


axis <- list(title = "",
             showgrid = FALSE,
             showticklabels = FALSE,
             zeroline = FALSE)

# ----- PLOT PATHWAY TOTAL -----
p <- plot_ly() %>%
  
  add_segments(data = edge_data_frame,
               x = ~x_from, xend = ~x_to,
               y = ~y_from, yend = ~y_to,
               alpha = 0.15, size = I(1), color = I("grey")
  ) %>%
  add_trace(data = layout_2d,
            x = ~PC1, y = ~PC2, size = I(12),
            # shapes = edge_shapes,
            color = ~path
  ) %>%
  layout(title = '',
         xaxis = axis,
         yaxis = axis)

# ----- PLOT PATHWAY SIG -----
# color
layout_2d_sub1 <- layout_2d[which(layout_2d$sig == 'Yes'),]

sig_up <- which(layout_2d_sub1$FC == 'up')
sig_up_1 <- sample(sig_up,18)
sig_up <- sig_up[! sig_up %in% sig_up_1]
sig_up_2 <- sample(sig_up,18)
sig_up_3 <- sig_up[! sig_up %in% sig_up_2]

sig_dn <- which(layout_2d_sub1$FC == 'dn')
sig_dn_1 <- sample(sig_dn,6)
sig_dn <- sig_dn[! sig_dn %in% sig_dn_1]
sig_dn_2 <- sample(sig_dn,6)
sig_dn_3 <- sig_dn[! sig_dn %in% sig_dn_2]

layout_2d_sub1$color[sig_up_1] <- '#FF8F8F'
layout_2d_sub1$color[sig_up_2] <- '#FF4B4B'
layout_2d_sub1$color[sig_up_3] <- '#FF0000'

layout_2d_sub1$color[sig_dn_1] <- '#939DFF'
layout_2d_sub1$color[sig_dn_2] <- '#5463FF'
layout_2d_sub1$color[sig_dn_3] <- '#0017FF'


layout_2d_sub1$size <- NA
layout_2d_sub1$size[layout_2d_sub1$p_val>0.05] <- 7
layout_2d_sub1$size[layout_2d_sub1$p_val<=0.05 & layout_2d_sub1$p_val>0.01] <- 12
layout_2d_sub1$size[layout_2d_sub1$p_val<=0.01] <- 17

# grey
layout_2d_sub2 <- layout_2d[which(layout_2d$sig == 'No'),]
layout_2d_sub2$color <- '#a8a8a8'
layout_2d_sub2$size <- 7

# plot
p <- plot_ly() %>%
  
  add_segments(
    data = edge_data_frame,
    x = ~x_from, xend = ~x_to,
    y = ~y_from, yend = ~y_to,
    alpha = 0.2, size = I(0.5), color = I("grey")
  ) %>%
  
  add_trace(data = layout_2d_sub2, x = ~PC1, y = ~PC2,
            marker = list(size = ~size,
                          color = ~color)) %>%
  add_trace(data = layout_2d_sub1, x = ~PC1, y = ~PC2,
            marker = list(size = ~size,
                          color = ~color)) %>%
  layout(title = '',
         # shapes = edge_shapes,
         xaxis = axis,
         yaxis = axis) 











# ------ 2. metabolomics story ------

# ----- input concentration matrix (use column free method)-----
load('concentration_PE.Rdata')

# ----- clinic info -----
NP_clinic <- read.xlsx2(file = 'Preeclampsia Ethnicity new 2.xlsx',sheetIndex = 1,stringsAsFactors = F)
NP_clinic <- NP_clinic[,c(1:7)]
PE_clinic <- read.xlsx2(file = 'Preeclampsia Ethnicity new 2.xlsx',sheetIndex = 2,stringsAsFactors = F)
PE_clinic <- PE_clinic[,c(1:7)]
clinic <- rbind(PE_clinic,NP_clinic)
clinic$GA <- as.numeric(clinic$GA)
rm(NP_clinic,PE_clinic)
# blood presure
blood <- read.csv("mProbe_data-1.csv")
blood$blood_pressure <- substr(blood$blood_pressure,1,3)
blood$X <- as.character(blood$X)

# ----- preprocessing concentration mat -----
conc_PE <- conc_PE[,c(-2,-(67:74))]
conc_PE$Compound <- trim.ws.custom(as.character(conc_PE$Compound))
# convert to NA or 0
conc_PE[is.nan.data.frame(conc_PE)] <- 0
# value imputation
conc_PE[,2:65] <- apply(conc_PE[,2:65],2,function(x){
  tmp <- x
  tmp[abs(tmp) < 10^-10] <- NA
  return(tmp)
})
# normalization
tmp_colnum <- dim(conc_PE)[2]
tmp_median <- apply(conc_PE[,2:tmp_colnum],1,function(x){median(x,na.rm = T)})
tmp_std <- apply(conc_PE[,2:tmp_colnum],1,function(x){sd(x,na.rm = T)})
conc_PE_z <- conc_PE
conc_PE_z[,2:tmp_colnum] <- apply(conc_PE_z[,2:tmp_colnum], 2,
                                  function(x){return((x-tmp_median)/tmp_std)})
# feature removal
row_index <- apply(conc_PE_z[,2:tmp_colnum], 1, function(x) {
  sum(is.na(x)) > (tmp_colnum/2)})
conc_PE_z <- conc_PE_z[!row_index,]

# imputation
tmp_feature <- conc_PE_z$Compound
tmp <- impute.knn(as.matrix(conc_PE_z[,2:tmp_colnum]), k = 20)
tmp <- tmp$data
conc_PE_z_flt <- cbind(tmp_feature,data.frame(tmp))
colnames(conc_PE_z_flt)[1] <- 'Compound'
rm(tmp_colnum,tmp_median,tmp_std)

# ----- using QC to -----
# source("https://bioconductor.org/biocLite.R")
# biocLite("statTarget")
install.packages('RGtk2')
library(statTarget)

# ----- read from SAM metabolite analysis result (fold change) for GSEA analysis -----
# total
total_result <- read.xlsx2(file = 'pe_sam_stats.xlsx',sheetIndex = 2,stringsAsFactors = F)
total_result$local.FDR... <- as.numeric(total_result$local.FDR...)
total_result$FC <- as.numeric(total_result$FC)
total_result <- total_result[which(total_result$local.FDR... <= 5),]
total_result$FC <- log2(total_result$FC)
colnames(total_result) <- c('name','coeff','sth','p.spearman')
total_result$p.spearman <- 0.01*total_result$p.spearman

total_list <- as.numeric(total_result$FC)
total_list <- log2(total_list)
names(total_list) <- as.character(total_result$Compound)
total_list <- total_list[order(-total_list)]
# early
early_result <- read.xlsx2(file = 'pe_sam_stats.xlsx',sheetIndex = 1,stringsAsFactors = F)
early_result$local.FDR... <- as.numeric(early_result$local.FDR...)
early_result$FC <- as.numeric(early_result$FC)
early_result <- early_result[which(early_result$local.FDR... <= 5),]
early_result$FC <- log2(early_result$FC)
colnames(early_result) <- c('name','coeff','sth','p.spearman')
early_result$p.spearman <- 0.01*early_result$p.spearman

early_list <- as.numeric(early_result$FC)
early_list <- log2(early_list)
names(early_list) <- as.character(early_result$Compound)
early_list <- early_list[order(-early_list)]

early_list <- early_result$name
# severe
severe_result <- read.xlsx2(file = 'pe_sam_stats.xlsx',sheetIndex = 4,stringsAsFactors = F)
severe_result$local.FDR... <- as.numeric(severe_result$local.FDR...)
severe_result$FC <- as.numeric(severe_result$FC)
severe_result <- severe_result[which(severe_result$local.FDR... <= 5),]
severe_result$FC <- log2(severe_result$FC)
colnames(severe_result) <- c('name','coeff','sth','p.spearman')
severe_result$p.spearman <- 0.01*severe_result$p.spearman

severe_list <- as.numeric(severe_result$FC)
severe_list <- log2(severe_list)
names(severe_list) <- as.character(severe_result$Compound)
severe_list <- severe_list[order(-severe_list)]

severe_list <- severe_result$Compound

# late
late_result <- read.xlsx2(file = 'pe_sam_stats.xlsx',sheetIndex = 3,stringsAsFactors = F)
late_result$local.FDR... <- as.numeric(late_result$local.FDR...)
late_result$FC <- as.numeric(late_result$FC)
late_result <- late_result[which(late_result$local.FDR... <= 5),]
late_result$FC <- log2(late_result$FC)

late_list <- late_result$Compound

intersect(early_list,late_list)

# mild
mild_result <- read.xlsx2(file = 'pe_sam_stats.xlsx',sheetIndex = 5,stringsAsFactors = F)
mild_result$local.FDR... <- as.numeric(mild_result$local.FDR...)
mild_result$FC <- as.numeric(mild_result$FC)
mild_result <- mild_result[which(mild_result$local.FDR... <= 5),]
mild_result$FC <- log2(mild_result$FC)

mild_list <- mild_result$Compound

# intersect(early_list,late_list)
# intersect(severe_list,mild_list)
# 
# tmp <- data.frame(e.vs.l = c(intersect(early_list,late_list),rep(NA,84)),
#                   s.vs.m = intersect(severe_list,mild_list),
#                   stringsAsFactors = F)
# write.xlsx2(tmp,file = 'intersect.xlsx',showNA = F, row.names = F)

# ----- construct pathway list -----
load('kegg.path.Rdata')
pathway_mapping <- read.xlsx2(file = '/Users/zhijiliu/Google Drive/Mac Rstudio/Pathway analysis/PE/Protein_metabolite_pathway_mapping/analytes_mapping_update_JULY.xlsx',sheetIndex = 1, stringsAsFactors = F)
selected_path_data <- read.xlsx2(file = '/Users/zhijiliu/Google Drive/Mac Rstudio/Pathway analysis/PE/Promedex/16_pathways.xlsx',
                                 sheetIndex = 1,stringsAsFactors = F)
# modification to pathways
pathway_mapping$human_path_id[42:46] <- pathway_mapping$human_path_id[47]
pathway_mapping$Compound <- trim.ws.custom(pathway_mapping$Compound)

pathways_set <- list()
for (i in 1:nrow(selected_path_data)){
  tmp_path <- selected_path_data$path.id[i]
  tmp_cpd <- pathway_mapping$Compound[grep(tmp_path,pathway_mapping$human_path_id)]
  pathways_set[[i]] <- tmp_cpd
  names(pathways_set)[i] <- tmp_path
}
# separate glycesophospolipid metabolism
tmp_cpd <- pathways_set[[which(names(pathways_set)=='map00564')]]

# PS
pathways_set[[which(names(pathways_set) == 'gly_PS')]] <- tmp_cpd[grep('PS',tmp_cpd)]
# PI
pathways_set[[which(names(pathways_set) == 'gly_PI')]] <- tmp_cpd[grep('PI',tmp_cpd)]
# PG
pathways_set[[which(names(pathways_set) == 'gly_PG')]] <- tmp_cpd[grep('PG',tmp_cpd)]
# PE
pathways_set[[which(names(pathways_set) == 'gly_PE')]] <- tmp_cpd[grep('PE',tmp_cpd)]
# PC
pathways_set[[which(names(pathways_set) == 'gly_PC')]] <- tmp_cpd[grep('PC',tmp_cpd)]
# PA
pathways_set[[which(names(pathways_set) == 'gly_PA')]] <- tmp_cpd[grep('PA',tmp_cpd)]

# ----- sig pathway list for total -----
result_total <- fgsea(pathways = pathways_set,stats = total_list,nperm = 1000,minSize = 10)
result_total$pathway_name <- selected_path_data$path.name[match(result_total$pathway,selected_path_data$path.id)]
write.xlsx2(result_total,file = 'enrichment_score.xlsx',sheetName = 'total',showNA = F,row.names = F,append = F)
# pathway list
tmp.path <- result_total$pathway[which(result_total$pval <= 0.05)]

path.pval.total <- vector("list", length(tmp.path)) 
for (i in 1:length(tmp.path)){
  names(path.pval.total)[i] <- tmp.path[i]
  path.pval.total[[i]] <- pathways_set[[which(names(pathways_set)==tmp.path[i])]]
}

# ----- sig pathway list for early -----
result_early <- fgsea(pathways = pathways_set,stats = early_list,nperm = 1000,minSize = 10)
result_early$pathway_name <- selected_path_data$path.name[match(result_early$pathway,selected_path_data$path.id)]
write.xlsx2(result_early,file = 'enrichment_score.xlsx',sheetName = 'early',showNA = F,row.names = F,append = T)
# pathway list
tmp.path <- result_early$pathway[which(result_early$pval <= 0.05)]

path.pval.early <- vector("list", length(tmp.path)) 
for (i in 1:length(tmp.path)){
  names(path.pval.early)[i] <- tmp.path[i]
  path.pval.early[[i]] <- pathways_set[[which(names(pathways_set)==tmp.path[i])]]
}

# ----- sig pathway list for severe -----
result_severe <- fgsea(pathways = pathways_set,stats = severe_list,nperm = 1000,minSize = 10)
result_severe$pathway_name <- selected_path_data$path.name[match(result_severe$pathway,selected_path_data$path.id)]
write.xlsx2(result_severe,file = 'enrichment_score.xlsx',sheetName = 'severe',showNA = F,row.names = F,append = T)
# pathway list
tmp.path <- result_severe$pathway[which(result_severe$pval <= 0.05)]

path.pval.severe <- vector("list", length(tmp.path)) 
for (i in 1:length(tmp.path)){
  names(path.pval.severe)[i] <- tmp.path[i]
  path.pval.severe[[i]] <- pathways_set[[which(names(pathways_set)==tmp.path[i])]]
}

# ----- get sample GSEA object -----
data(geneList, package="DOSE")
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

# ----- construct GSEA object -----
total_result <- kk2
total_result@result <- result_total
total_result@geneSets <- path.pval.total
total_result@geneList <- total_list
total_result@params$minGSSize <- 10

early_result <- kk2
early_result@result <- result_early
early_result@geneSets <- path.pval.early
early_result@geneList <- early_list
early_result@params$minGSSize <- 10

severe_result <- kk2
severe_result@result <- result_severe
severe_result@geneSets <- path.pval.severe
severe_result@geneList <- severe_list
severe_result@params$minGSSize <- 10

# ----- start GSEA plot -----
# change result parameter here
gseaResult = severe_result
tmp.path <- result_severe$pathway[which(result_severe$pval <= 0.05)]
tmp3 <- tmp.path
# parameters: gseaResult, geneSetID, by = "all", title = "", color='black', color.line="green", color.vline="#FA5860")
# by <- match.arg(by, c("runningScore", "preranked", "all"))
for (i in 1:length(tmp3)){
  geneSetID = tmp3[i]
  x <- y <- ymin <- ymax <- xend <- yend <- runningScore <- es <- pos <- geneList <- NULL
  color='black'
  color.line="green"
  color.vline="#FA5860"
  title = selected_path_data$path.name[which(selected_path_data$path.id == tmp3[i])]
  gsdata <- gsInfo(gseaResult, tmp3[i])
  p <- ggplot(gsdata, aes(x = x)) +
    theme_dose() + xlab("Position in the Ranked List of Metabolites")
  
  # "barcode"
  p.res <- p + geom_linerange(aes(ymin=ymin, ymax=ymax), color=color)
  # running line
  p.res <- p.res + geom_line(aes(y = runningScore), color=color.line, size=1)
  enrichmentScore <- gseaResult@result$ES[which(gseaResult@result$pathway==geneSetID)]
  
  es.df <- data.frame(es = abs(which.min(p$data$runningScore - enrichmentScore)))
  p.res <- p.res + geom_vline(data = es.df, aes(xintercept = es),
                              colour = color.vline, linetype = "dashed")
  p.res <- p.res + ylab("Running Enrichment Score")
  p.res <- p.res + geom_hline(aes(yintercept = 0))
  
  
  df2 <- data.frame(x = which(p$data$position == 1))
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data=df2, aes(x=x, xend=x, y=y, yend=0), color=color)
  
  ## p.pos <- p + geom_vline(data = df2, aes(xintercept = pos),
  ##                         colour = "#DAB546", alpha = 0.3)
  ## p.pos <- p.pos + geom_line(aes(y = geneList), colour = "red")
  ## p.pos <- p.pos + geom_hline(aes(yintercept = 0))
  p.pos <- p.pos + ylab("Ranked list metric") + xlim(0, length(p$data$geneList))
  
  p.pos <- p.pos + xlab(NULL) + theme(axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank())
  ##p.res <- p.res + theme(axis.title.x = element_text(vjust = -0.3))
  
  gp1<- ggplot_gtable(ggplot_build(p.res))
  gp2<- ggplot_gtable(ggplot_build(p.pos)) 
  maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3])
  gp1$widths[2:3] <- maxWidth
  gp2$widths[2:3] <- maxWidth
  text.params <- gpar(fontsize=15, fontface="bold", lineheight=0.8)
  textgp <- textGrob(title, gp=text.params)
  
  ## grid.arrange(textgp, gp2, gp1, ncol=1, heights=c(0.1, 0.7, 0.7))
  
  if (dev.interactive())
    grid.newpage()
  png(paste('severe ',title,'.png',sep = ''),height = 3000, width = 3000,res = 550)
  pushViewport(viewport(layout = grid.layout(3, 1, heights = unit(c(0.1, 0.7, 0.7), "null"))))
  
  gp2$vp = viewport(layout.pos.row = 2, layout.pos.col = 1)
  grid.draw(gp2)
  
  gp1$vp = viewport(layout.pos.row = 3, layout.pos.col = 1)
  grid.draw(gp1)
  
  textgp$vp = viewport(layout.pos.row = 1, layout.pos.col = 1)
  grid.draw(textgp)
  
  invisible(list(runningScore = p.res, preranked = p.pos))
  
  dev.off()
}

# ----- metabolites network plot -----
# remove glycerophosfolipid pathway
selected_path_data_cpd <- selected_path_data
selected_path_data_cpd <- selected_path_data_cpd[-7,]
pathways_set_cpd <- pathways_set
pathways_set_cpd <- pathways_set_cpd[-7]

# selected_path <- selected_path_data_cpd$path.id
# cpd <- data.frame(cpd.name = unique(as.character(unlist(pathways_set_cpd))),
#                   stringsAsFactors = F)
# 
# cpd$path.id <- NA
# for (i in 1:dim(cpd)[1]){
#   tmp.path <- pathway_mapping_2$human_path_id[match(cpd$cpd.name[i],pathway_mapping_2$Compound)]
#   tmp.path <- unique(unlist(strsplit(tmp.path,',')))
#   tmp.path <- intersect(tmp.path,selected_path)
#   if (length(tmp.path)>1) tmp.path <- tmp.path[sample(length(tmp.path),1)]
#   cpd$path.id[i] <- tmp.path
# }

tmp.cpd <- trim.ws.custom(conc_PE_z_flt$Compound)
tmp.inter <- lapply(pathways_set_cpd,function(x){
  return(intersect(x,tmp.cpd))
})
tmp.length <- unlist(lapply(tmp.inter,length))
tmp.inter <- tmp.inter[order(-tmp.length)]

cpd <- data.frame(cpd.name = as.character(conc_PE_z_flt$Compound),
                  path.id = '',
                  stringsAsFactors = F)
for (i in 1:length(tmp.inter)){
  cpd$path.id[match(tmp.inter[[i]],cpd$cpd.name)] <- paste(cpd$path.id[match(tmp.inter[[i]],cpd$cpd.name)],
                                                           names(tmp.inter)[i])
}
cpd$path.id[which(cpd$path.id=='')] <- 'others'

cpd_2 <- cpd
cpd_2$path.id.select <- cpd_2$path.id
for (i in 1:dim(cpd_2)[1]){
  tmp <- unlist(strsplit(cpd_2$path.id.select[i],' '))
  tmp <- tmp[which(tmp!='')]
  cpd_2$path.id.select[i] <- sample(tmp,1)
}

cpd <- cpd_2

# ----- subset matrix based on selected pathways -----
conc_PE_z_flt$Compound <- trim.ws.custom(as.character(conc_PE_z_flt$Compound))
tmp.cpd <- intersect(conc_PE_z_flt$Compound,cpd$cpd.name)
# only select cpd that associated with the selected pathways
input_mat <- conc_PE_z_flt[match(tmp.cpd,conc_PE_z_flt$Compound),]
tmp.cpd <- input_mat$Compound
input_mat <- t(input_mat[,2:dim(input_mat)[2]])
colnames(input_mat) <- tmp.cpd
rm(tmp.cpd)

# ----- build correlation network based on feature-feature pearson correlation -----
# pearson correlation
cor_mat <- cor(input_mat,use = 'everything',method = 'pearson')
# adjacency matrix
adj_mat <- abs(cor_mat)
adj_mat[is.na(adj_mat)] <- 0
adj_mat <- apply(adj_mat,2,function(x){
  x[abs(x)<0.8] <- 0
  x[abs(x)>=0.8] <- 1
  return(x)
})
diag(adj_mat) <- 1
# PCA 
pe_pca <- prcomp(t(input_mat),
                 center = T,
                 scale. = T)
# PCA coordicates
fit_pca <- data.frame(pe_pca$x)
rm(pe_pca)
# get 2d layout
layout_2d <- fit_pca[,c(1,2)]
layout_2d$cpd <- rownames(layout_2d)
identical(layout_2d$cpd,cpd$cpd.name)
layout_2d$condition <- cpd$path.id.select
# manually adjust coordinates 
path_10 <- data.frame(path.id = unique(layout_2d$condition),stringsAsFactors = F)
# reset
layout_2d$PC1 <- fit_pca$PC1
layout_2d$PC2 <- fit_pca$PC2
#1. others
layout_2d$PC1[grep(path_10$path.id[1],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[1],layout_2d$condition)] +1
layout_2d$PC2[grep(path_10$path.id[1],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[1],layout_2d$condition)] +2
#2. 01212
layout_2d$PC1[grep(path_10$path.id[2],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[2],layout_2d$condition)] +13
layout_2d$PC2[grep(path_10$path.id[2],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[2],layout_2d$condition)] +28
#3. 00071
layout_2d$PC1[grep(path_10$path.id[3],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[3],layout_2d$condition)] +20
layout_2d$PC2[grep(path_10$path.id[3],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[3],layout_2d$condition)] +6
#4. 01230
layout_2d$PC1[grep(path_10$path.id[4],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[4],layout_2d$condition)] +22
layout_2d$PC2[grep(path_10$path.id[4],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[4],layout_2d$condition)] -3
#5. 01210
layout_2d$PC1[grep(path_10$path.id[5],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[5],layout_2d$condition)] -15
layout_2d$PC2[grep(path_10$path.id[5],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[5],layout_2d$condition)] -6
#6. 00260
layout_2d$PC1[grep(path_10$path.id[6],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[6],layout_2d$condition)] +3
layout_2d$PC2[grep(path_10$path.id[6],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[6],layout_2d$condition)] +10
#7. 01040
layout_2d$PC1[grep(path_10$path.id[7],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[7],layout_2d$condition)] -7
layout_2d$PC2[grep(path_10$path.id[7],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[7],layout_2d$condition)] -30
#8. 00591
layout_2d$PC1[grep(path_10$path.id[8],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[8],layout_2d$condition)] -20
layout_2d$PC2[grep(path_10$path.id[8],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[8],layout_2d$condition)] +3
#9. 00600
layout_2d$PC1[grep(path_10$path.id[9],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[9],layout_2d$condition)] +1
layout_2d$PC2[grep(path_10$path.id[9],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[9],layout_2d$condition)] +21
#10. 00561
layout_2d$PC1[grep(path_10$path.id[10],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[10],layout_2d$condition)] -7
layout_2d$PC2[grep(path_10$path.id[10],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[10],layout_2d$condition)] -23
#11. gly_PA
layout_2d$PC1[grep(path_10$path.id[11],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[11],layout_2d$condition)] -5
layout_2d$PC2[grep(path_10$path.id[11],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[11],layout_2d$condition)] -9
#12. gly_PE
layout_2d$PC1[grep(path_10$path.id[12],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[12],layout_2d$condition)] -12
layout_2d$PC2[grep(path_10$path.id[12],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[12],layout_2d$condition)] -2
#13. 00565
layout_2d$PC1[grep(path_10$path.id[13],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[13],layout_2d$condition)] -11
layout_2d$PC2[grep(path_10$path.id[13],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[13],layout_2d$condition)] +16
#14. 00563
layout_2d$PC1[grep(path_10$path.id[14],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[14],layout_2d$condition)] -8
layout_2d$PC2[grep(path_10$path.id[14],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[14],layout_2d$condition)] +28
#15. gly_PG
layout_2d$PC1[grep(path_10$path.id[15],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[15],layout_2d$condition)] +14
layout_2d$PC2[grep(path_10$path.id[15],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[15],layout_2d$condition)] -20
#16. gly_PI
layout_2d$PC1[grep(path_10$path.id[16],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[16],layout_2d$condition)] +5
layout_2d$PC2[grep(path_10$path.id[16],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[16],layout_2d$condition)] -1
#17. 00562
layout_2d$PC1[grep(path_10$path.id[17],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[17],layout_2d$condition)] -7
layout_2d$PC2[grep(path_10$path.id[17],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[17],layout_2d$condition)] +7
#18. gly_PS
layout_2d$PC1[grep(path_10$path.id[18],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[18],layout_2d$condition)] +3
layout_2d$PC2[grep(path_10$path.id[18],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[18],layout_2d$condition)] -15
#19. 00100
layout_2d$PC1[grep(path_10$path.id[19],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[19],layout_2d$condition)] +15
layout_2d$PC2[grep(path_10$path.id[19],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[19],layout_2d$condition)] +15
#20. gly_PC
layout_2d$PC1[grep(path_10$path.id[20],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[20],layout_2d$condition)] +14
layout_2d$PC2[grep(path_10$path.id[20],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[20],layout_2d$condition)] -6
#21. 00590
layout_2d$PC1[grep(path_10$path.id[21],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[21],layout_2d$condition)] +11
layout_2d$PC2[grep(path_10$path.id[21],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[21],layout_2d$condition)] +4
#22. 00592
layout_2d$PC1[grep(path_10$path.id[22],layout_2d$condition)] <- layout_2d$PC1[grep(path_10$path.id[22],layout_2d$condition)] +12
layout_2d$PC2[grep(path_10$path.id[22],layout_2d$condition)] <- layout_2d$PC2[grep(path_10$path.id[22],layout_2d$condition)] -30

# save(layout_2d,file = '22.pathway.layout.Rdata')

# set edges
# global pathway settings
# nodes and others
g <- graph_from_adjacency_matrix(adjmatrix = adj_mat, mode = 'undirected',diag = T)
# Create Vertices and Edges
vs <- V(g)
es <- as.data.frame(get.edgelist(g))

Nv <- length(vs)
Ne <- length(es[1]$V1)

# axis
axis <- list(title = "",
             showgrid = FALSE,
             showticklabels = FALSE,
             zeroline = FALSE)


Xn <- layout_2d[,1]
Yn <- layout_2d[,2]
# edge coordinates
from <- es$V1
to <- es$V2
x1 = Xn[from]
x2 = Xn[to]
y1 = Yn[from]
y2 = Yn[to]
edge_data_frame <- data.frame(x_from = x1,
                              x_to = x2,
                              y_from = y1,
                              y_to = y2)

# plot (color)
p <- plot_ly() %>%
  
  add_segments(data = edge_data_frame,
               x = ~x_from, xend = ~x_to,
               y = ~y_from, yend = ~y_to,
               line = list(color = 'rgb(220, 220, 220)', width = 0.3)) %>%
  
  add_trace(data = layout_2d,
            x = ~PC1, y = ~PC2, size = I(10),
            # shapes = edge_shapes,
            color = ~condition
  ) %>%
  # marker = list(color = ~condition)) %>%
  
  layout(title = '',
         xaxis = axis,
         yaxis = axis)

plotly_IMAGE(p, format = "png", height = 1000,width = 1000,
             out_file = 'test.png')

# ----- stevenson sample ----- 
# load('/Users/zhijiliu/Downloads/input.x.nzero.rdata')
# coef.nzero <- coef.nzero[2:dim(coef.nzero)[1],]
# write.xlsx2(coef.nzero,file = 'features.stevenson.xlsx',showNA = F,row.names = F)

tmp.origin <- trim.ws.custom(conc_PE$Compound)
tmp.change <- make.names(gsub("/", ".to.", tmp.origin))
stevenson.name <- data.frame(origin = tmp.origin,
                             change = tmp.change,
                             stringsAsFactors = F)

coef.nzero <- read.xlsx2(file = '/Users/zhijiliu/Google Drive/Mac Rstudio/Pathway analysis/PE/Promedex/features.stevenson.xlsx',sheetIndex = 1,stringsAsFactors = F)
coef.nzero <- read.xlsx2(file = '/Users/zhijiliu/Downloads/nzero-feature.xlsx',sheetIndex = 1,stringsAsFactors = F)
coef.nzero <- read.xlsx2(file = '/Users/zhijiliu/Downloads/nzero-feature-1.xlsx',sheetIndex = 1,stringsAsFactors = F)



coef.nzero <- coef.nzero[2:32,]

colnames(coef.nzero)[2] <- 'coeff'
coef.nzero$coeff <- as.numeric(coef.nzero$coeff)
coef.nzero$p.spearman <- as.numeric(coef.nzero$p.spearman)
coef.nzero$name <- trim.ws.custom(coef.nzero$name)
coef.nzero$name <- stevenson.name$origin[match(coef.nzero$name,stevenson.name$change)]

# ----- stevenson -----
sig_up <- coef.nzero$name[coef.nzero$coeff>0]
sig_dn <- coef.nzero$name[coef.nzero$coeff<0]

sig_up_1 <- coef.nzero$name[coef.nzero$coeff > 0 & coef.nzero$coeff < 0.3]
sig_up_2 <- coef.nzero$name[coef.nzero$coeff >= 0.3 & coef.nzero$coeff < 0.6]
sig_up_3 <- coef.nzero$name[coef.nzero$coeff >= 0.6]

sig_dn_1 <- coef.nzero$name[coef.nzero$coeff < 0 & coef.nzero$coeff > -0.3]
sig_dn_2 <- coef.nzero$name[coef.nzero$coeff <= -0.3 & coef.nzero$coeff > -0.6]
sig_dn_3 <- coef.nzero$name[coef.nzero$coeff <= -0.6]

sig_large <- coef.nzero$name[which(coef.nzero$p.spearman < 1e-7)]
sig_mid <- coef.nzero$name[which(coef.nzero$p.spearman >= 1e-7 & coef.nzero$p.spearman < 1e-5)]
sig_small <- coef.nzero$name[which(coef.nzero$p.spearman >= 1e-5)]

# ----- early -----
sig_up <- early_result$name[early_result$coeff>0]
sig_dn <- early_result$name[early_result$coeff<0]

sig_up_1 <- early_result$name[early_result$coeff > 0 & early_result$coeff < 0.3]
sig_up_2 <- early_result$name[early_result$coeff >= 0.3 & early_result$coeff < 0.6]
sig_up_3 <- early_result$name[early_result$coeff >= 0.6]

sig_dn_1 <- early_result$name[early_result$coeff < 0 & early_result$coeff > -0.3]
sig_dn_2 <- early_result$name[early_result$coeff <= -0.3 & early_result$coeff > -0.6]
sig_dn_3 <- early_result$name[early_result$coeff <= -0.6]

sig_large <- early_result$name[which(early_result$p.spearman < 0.01)]
sig_mid <- early_result$name[which(early_result$p.spearman >= 0.01 & early_result$p.spearman < 0.03)]
sig_small <- early_result$name[which(early_result$p.spearman >= 0.03)]

# ----- severe -----
sig_up <- severe_result$name[severe_result$coeff>0]
sig_dn <- severe_result$name[severe_result$coeff<0]

sig_up_1 <- severe_result$name[severe_result$coeff > 0 & severe_result$coeff < 0.3]
sig_up_2 <- severe_result$name[severe_result$coeff >= 0.3 & severe_result$coeff < 0.6]
sig_up_3 <- severe_result$name[severe_result$coeff >= 0.6]

sig_dn_1 <- severe_result$name[severe_result$coeff < 0 & severe_result$coeff > -0.3]
sig_dn_2 <- severe_result$name[severe_result$coeff <= -0.3 & severe_result$coeff > -0.6]
sig_dn_3 <- severe_result$name[severe_result$coeff <= -0.6]

sig_large <- severe_result$name[which(severe_result$p.spearman < 0.01)]
sig_mid <- severe_result$name[which(severe_result$p.spearman >= 0.01 & severe_result$p.spearman < 0.03)]
sig_small <- severe_result$name[which(severe_result$p.spearman >= 0.03)]

# ----- total -----
sig_up <- total_result$name[total_result$coeff>0]
sig_dn <- total_result$name[total_result$coeff<0]

sig_up_1 <- total_result$name[total_result$coeff > 0 & total_result$coeff < 0.3]
sig_up_2 <- total_result$name[total_result$coeff >= 0.3 & total_result$coeff < 0.6]
sig_up_3 <- total_result$name[total_result$coeff >= 0.6]

sig_dn_1 <- total_result$name[total_result$coeff < 0 & total_result$coeff > -0.3]
sig_dn_2 <- total_result$name[total_result$coeff <= -0.3 & total_result$coeff > -0.6]
sig_dn_3 <- total_result$name[total_result$coeff <= -0.6]

sig_large <- total_result$name[which(total_result$p.spearman < 0.01)]
sig_mid <- total_result$name[which(total_result$p.spearman >= 0.01 & total_result$p.spearman < 0.03)]
sig_small <- total_result$name[which(total_result$p.spearman >= 0.03)]

# ----- set color and size -----
# color
layout_2d$label <- NA
layout_2d$label[match(sig_up,layout_2d$cpd)] <- 'up'
layout_2d$label[match(sig_dn,layout_2d$cpd)] <- 'down'
layout_2d$label[is.na(layout_2d$label)] <- 'non'
rm(sig_up,sig_dn)

layout_2d$color <- NA
layout_2d$color[match(sig_up_1,layout_2d$cpd)] <- '#FF8F8F'
layout_2d$color[match(sig_up_2,layout_2d$cpd)] <- '#FF4B4B'
layout_2d$color[match(sig_up_3,layout_2d$cpd)] <- '#FF0000'

layout_2d$color[match(sig_dn_1,layout_2d$cpd)] <- '#939DFF'
layout_2d$color[match(sig_dn_2,layout_2d$cpd)] <- '#5463FF'
layout_2d$color[match(sig_dn_3,layout_2d$cpd)] <- '#0017FF'

layout_2d$color[which(layout_2d$label=='non')] <- '#a8a8a8'

# label compound size according to pvalue
layout_2d$size <- NA
layout_2d$size[which(layout_2d$label == 'non')] <- 9

layout_2d$size[match(sig_large,layout_2d$cpd)] <- 19
layout_2d$size[match(sig_mid,layout_2d$cpd)] <- 15
layout_2d$size[match(sig_small,layout_2d$cpd)] <- 11

# grey
layout_sub_1 <- layout_2d[which(layout_2d$label=="non"),]
# color
layout_sub_2 <- layout_2d[which(layout_2d$label!="non"),]

p <- plot_ly() %>%
  
  add_segments(
    data = edge_data_frame,
    x = ~x_from, xend = ~x_to,
    y = ~y_from, yend = ~y_to,
    alpha = 0.2, size = I(1), color = I("grey")
  ) %>%
  
  add_trace(data = layout_sub_1, x = ~PC1, y = ~PC2,
            marker = list(size = ~size,
                          color = ~color)) %>%
  add_trace(data = layout_sub_2, x = ~PC1, y = ~PC2,
            marker = list(size = ~size,
                          color = ~color)) %>%
  layout(title = '',
         # shapes = edge_shapes,
         xaxis = axis,
         yaxis = axis) 

# plotly_IMAGE(p, format = "png", height = 1000,width = 1000,
#              out_file = 'ptb.model.sig.32.png')

# ----- plot for each coloful pathway -----
path_final <- c('map00071',
                'gly_PG',
                'gly_PE',
                'gly_PA',
                'map01210',
                'map00100',
                'map01230',
                'map01212',
                'map00600')
color <- c()
for (i in 1:length(path_final)){
  # get cpd
  tmp.cpd <- pathways_set[[which(names(pathways_set) == path_final[i])]]
  tmp.cpd <- intersect(cpd,conc_PE_z_flt$Compound)
  
  # get adj matrix for network edges
  
  # set colors 
  
  # plot
  p <- plot_ly() %>%
    
    add_segments(
      data = edge_data_frame,
      x = ~x_from, xend = ~x_to,
      y = ~y_from, yend = ~y_to,
      alpha = 0.2, size = I(1), color = I("grey")
    ) %>%
    
    add_trace(data = layout_sub_1, x = ~PC1, y = ~PC2,
              marker = list(size = ~size,
                            color = ~color)) %>%
    add_trace(data = layout_sub_2, x = ~PC1, y = ~PC2,
              marker = list(size = ~size,
                            color = ~color)) %>%
    layout(title = '',
           # shapes = edge_shapes,
           xaxis = axis,
           yaxis = axis) 
}

# ------ 3. modeling ------
# ----- normalize based on case -----
conc_E_z <- conc_PE[,c(1,grep('E',colnames(conc_PE)))]
tmp_colnum <- dim(conc_E_z)[2]
tmp_median <- apply(conc_E_z[,2:tmp_colnum],1,function(x){median(x,na.rm = T)})
tmp_std <- apply(conc_E_z[,2:tmp_colnum],1,function(x){sd(x,na.rm = T)})
conc_E_z[,2:tmp_colnum] <- apply(conc_E_z[,2:tmp_colnum], 2,
                                 function(x){return((x-tmp_median)/tmp_std)})
# feature removal
row_index <- apply(conc_E_z[,2:tmp_colnum], 1, function(x) {
  sum(is.na(x)) > (tmp_colnum/2)})
conc_E_z <- conc_E_z[!row_index,]

# imputation
tmp_feature <- conc_E_z$Compound
tmp <- impute.knn(as.matrix(conc_E_z[,2:tmp_colnum]), k = 20)
tmp <- tmp$data
conc_E_z_flt <- cbind(tmp_feature,data.frame(tmp))
colnames(conc_E_z_flt)[1] <- 'Compound'
rm(tmp_colnum,tmp_median,tmp_std)

# ----- RANDOM_FOREST MODEL FOR "0 AND 1" -----
# # build training dataset: global
# clinic$label <- NA

# write.xlsx2(clinic,file = 'clinic.xlsx',showNA = F, row.names = F)
clinic <- read.xlsx2(file = 'clinic.xlsx',sheetIndex = 1,stringsAsFactors = F)

# function for train, testing, and plots
multi_task <- function(input_mat, pathway_id, clinic, blood, mode){
  if (mode == 'binary'){
    ## either case or control
    # get train set
    train_id <- clinic$ID[which(clinic$label == 'train')]
    # train_id <- c("E17", "E24", "E16", "E20", "E1",  "E26", "E7",  "E30", "E11", "E19", "E31", "E23", "E29", 
    #               "E3",  "E5",  "E25", "E6",  "E21", "E9",  "E14", "E13")
    train_mat <- t(input_mat[,train_id])
    colnames(train_mat) <- as.character(input_mat$Compound)
    train_class <- as.factor(substr(rownames(train_mat),1,1))
    # modeling with random forest
    fit <- randomForest(train_mat, y=train_class, ntree=1000,
                        replace=F, importance=TRUE)
    
    # # modeling with EN and get selected cpd
    # fit <- glmnet(x = train_mat, y = train_class, family = 'binomial')
    # cvfit = cv.glmnet(x= train_mat, y= train_class, family = "binomial", type.measure = "class")
    # tmp <- coef(cvfit, s = "lambda.min")
    # tmp2 <- tmp@Dimnames[[1]][tmp@i+1]
    
    # ---- whisker part ----
    # prediction on train set
    tmp_pred <- predict(fit,train_mat,type = "prob")
    train.score.all <- data.frame(score = tmp_pred[,1])
    train.score.all$label <- c(rep('PE',length(grep('E',rownames(train.score.all)))),rep('NP',length(grep('P',rownames(train.score.all)))))
    train.score.all$blood_presure <- blood$blood_pressure[match(rownames(train.score.all),blood$X)]
    # early train
    train.early.id <- clinic$ID[which(clinic$GA <= 34 & clinic$label == 'train')]
    train.score.early <- train.score.all[match(train.early.id,rownames(train.score.all)),]
    # severe train
    train.severe.id <- clinic$ID[which((clinic$Status2 == 'severe' | clinic$Status2 == 'normal') & 
                                         clinic$label == 'train')]
    train.score.severe <- train.score.all[match(train.severe.id,rownames(train.score.all)),]
    
    # # whisker plot for training set
    # var <- c('train.score.all','train.score.early','train.score.severe')
    # for (i in 1:3){
    #   tmp <- get(var[i])
    #   png(file = paste(pathway_id,'.',var[i],".png",sep = ""), res = 300, width = 1500, height = 1500)
    #   g <- ggplot(data = tmp, aes(x = label, y = score, group = label, fill = as.factor(label))) +
    #     geom_boxplot(outlier.shape = 1) +
    #     theme_bw() +
    #     theme(legend.title=element_blank(),legend.position=c(0.15,0.85),axis.text.x = element_blank(),axis.ticks.x=element_blank())+
    #     scale_fill_manual(breaks=as.factor(c(1,0)),labels=c("Case","Control"),values = c("seagreen", "violetred"))+
    #     labs(x="", y = "Predict score")
    #   print(g)
    #   dev.off()
    # }
    
    # beeswarm plot for training set
    var <- c('train.score.all','train.score.early','train.score.severe')
    for (i in 1:3){
      tmp <- get(var[i])
      png(file = paste(pathway_id,'.',var[i],".png",sep = ""), res = 300, width = 1500, height = 1500)
      beeswarm(score ~ label, data = tmp,
               col = c('seagreen','violetred'), pch = 18, cex = 2, 
               cex.main=2,
               xlab = "", ylab = '', cex.lab=1.5,cex.axis=1.5,font.axis=2,font.lab=2)
      # corral = "wrap", 
      # labels=tmp$label)
      boxplot(tmp$score ~ tmp$label,add=T,col= '#FFFFFF22',axes=FALSE,lwd=1.2,
              outline=F)
      box(lwd=4)
      dev.off()
    }
    
    # build testing set
    testing.ids <- clinic$ID[which(clinic$label == 'testing')]
    testing.data <- input_mat[,testing.ids]
    testing.data <- t(testing.data)
    colnames(testing.data) <- as.character(input_mat$Compound)
    testing.class <- as.factor(substr(rownames(testing.data),1,1))
    
    # prediction on testing set
    # total test
    all_pred_test <- predict(fit,testing.data,type = "prob")
    test.score.all <- data.frame(score = all_pred_test[,1])
    test.score.all$label <- c(rep('PE',length(grep('E',rownames(test.score.all)))),rep('NP',length(grep('P',rownames(test.score.all)))))
    test.score.all$blood_presure <- blood$blood_pressure[match(rownames(test.score.all),blood$X)]
    # early test
    test.early.id <- clinic$ID[which(clinic$GA <= 34 & clinic$label == 'testing')]
    test.score.early <- test.score.all[match(test.early.id,rownames(test.score.all)),]
    # severe test
    test.severe.id <- clinic$ID[which((clinic$Status2 == 'severe' | clinic$Status2 == 'normal') & 
                                        clinic$label == 'testing')]
    test.score.severe <- test.score.all[match(test.severe.id,rownames(test.score.all)),]
    
    # # whisker plot for testing set
    # var <- c('test.score.all','test.score.early','test.score.severe')
    # for (i in 1:3){
    #   tmp <- get(var[i])
    #   png(file = paste(pathway_id,'.',var[i],".png",sep = ""), res = 300, width = 1500, height = 1500)
    #   g <- ggplot(data = tmp, aes(x = label, y = score, group = label, fill = as.factor(label))) +
    #     geom_boxplot(outlier.shape = 1) +
    #     theme_bw() +
    #     theme(legend.title=element_blank(),legend.position=c(0.15,0.85),axis.text.x = element_blank(),axis.ticks.x=element_blank())+
    #     scale_fill_manual(breaks=as.factor(c(1,0)),labels=c("Case","Control"),values = c("seagreen", "violetred"))+
    #     labs(x="", y = "Predict score")
    #   print(g)
    #   dev.off()
    # }
    
    # beeswarm plot for training set
    var <- c('test.score.all','test.score.early','test.score.severe')
    for (i in 1:3){
      tmp <- get(var[i])
      png(file = paste(pathway_id,'.',var[i],".png",sep = ""), res = 300, width = 1500, height = 1500)
      beeswarm(score ~ label, data = tmp,
               col = c('seagreen','violetred'), pch = 18, cex = 2, 
               cex.main=2,
               xlab = "", ylab = '', cex.lab=1.5,cex.axis=1.5,font.axis=2,font.lab=2)
      # corral = "wrap", 
      # labels=tmp$label)
      boxplot(tmp$score ~ tmp$label,add=T,col= '#FFFFFF22',axes=FALSE,lwd=1.2,
              outline=F)
      box(lwd=4)
      dev.off()
    }
    
    # performance curve
    # all
    all_id <- clinic$ID
    all_mat <- t(input_mat[,all_id])
    colnames(all_mat) <- as.character(input_mat$Compound)
    all_score <- predict(fit,all_mat,type = "prob")
    all_score <- data.frame(score = all_score[,1])
    pred.prob <- all_score$score
    class <- as.factor(c(rep(1,length(grep('E',rownames(all_score)))),rep(0,length(grep('P',rownames(all_score))))))
    pred <- prediction(pred.prob,class)
    perf <- performance(pred, "auc" )@y.values[[1]]
    perff = performance(pred, "tpr", "fpr" )
    message(paste('all:',perf))
    
    png(file = paste(pathway_id,'.all.auc',".png",sep = ""), res = 300, width = 1500, height = 1500)
    plot(perff@x.values[[1]],perff@y.values[[1]],lty=1,lwd=4,cex.main=2,cex.lab=1.5,cex.axis=1.5,type='l',
         font.axis=2,font.lab=2,xlab='',ylab='')
    dev.off()
    
    # performance curve
    # early
    all_id <- clinic$ID[which(clinic$GA<34)]
    all_mat <- t(input_mat[,all_id])
    colnames(all_mat) <- as.character(input_mat$Compound)
    all_score <- predict(fit,all_mat,type = "prob")
    all_score <- data.frame(score = all_score[,1])
    pred.prob <- all_score$score
    class <- as.factor(c(rep(1,length(grep('E',rownames(all_score)))),rep(0,length(grep('P',rownames(all_score))))))
    pred <- prediction(pred.prob,class)
    perf <- performance(pred, "auc" )@y.values[[1]]
    perff = performance(pred, "tpr", "fpr" )
    message(paste('early:',perf))
    
    png(file = paste(pathway_id,'.early.auc',".png",sep = ""), res = 300, width = 1500, height = 1500)
    plot(perff@x.values[[1]],perff@y.values[[1]],lty=1,lwd=4,cex.main=2,cex.lab=1.5,cex.axis=1.5,type='l',
         font.axis=2,font.lab=2,xlab='',ylab='')
    dev.off()
    
    # performance curve
    # severe
    all_id <- clinic$ID[which(clinic$Status2 == 'severe' | clinic$Status2 == 'normal')]
    all_mat <- t(input_mat[,all_id])
    colnames(all_mat) <- as.character(input_mat$Compound)
    all_score <- predict(fit,all_mat,type = "prob")
    all_score <- data.frame(score = all_score[,1])
    pred.prob <- all_score$score
    class <- as.factor(c(rep(1,length(grep('E',rownames(all_score)))),rep(0,length(grep('P',rownames(all_score))))))
    pred <- prediction(pred.prob,class)
    perf <- performance(pred, "auc" )@y.values[[1]]
    perff = performance(pred, "tpr", "fpr" )
    message(paste('severe:',perf))
    
    png(file = paste(pathway_id,'.severe.auc',".png",sep = ""), res = 300, width = 1500, height = 1500)
    plot(perff@x.values[[1]],perff@y.values[[1]],lty=1,lwd=4,cex.main=2,cex.lab=1.5,cex.axis=1.5,type='l',
         font.axis=2,font.lab=2,xlab='',ylab='')
    dev.off()
    
    # # use all sample for train (no testing)
    # train_id <- clinic$ID
    # train_mat <- t(input_mat[,train_id])
    # colnames(train_mat) <- as.character(input_mat$Compound)
    # train_class <- as.factor(substr(rownames(train_mat),1,1))
    # # modeling with random forest
    # fit <- randomForest(train_mat, y=train_class, ntree=1000,
    #                     replace=F, importance=TRUE)
    # # prediction on train set
    # tmp_pred <- predict(fit,train_mat,type = "prob")
    # train.score.all <- data.frame(score = tmp_pred[,1])
    # train.score.all$label <- c(rep(1,length(grep('E',rownames(train.score.all)))),rep(0,length(grep('P',rownames(train.score.all)))))
    # train.score.all$blood_presure <- blood$blood_pressure[match(rownames(train.score.all),blood$X)]
    # # early train
    # train.early.id <- clinic$ID[which(clinic$GA <= 34)]
    # train.score.early <- train.score.all[match(train.early.id,rownames(train.score.all)),]
    # # severe train
    # train.severe.id <- clinic$ID[which(clinic$Status2 == 'severe' | clinic$Status2 == 'normal')]
    # train.score.severe <- train.score.all[match(train.severe.id,rownames(train.score.all)),]
    # # whisker plot for training set
    # var <- c('train.score.all','train.score.early','train.score.severe')
    # for (i in 1:3){
    #   tmp <- get(var[i])
    #   png(file = paste('All sample ',pathway_id,'.',var[i],".png",sep = ""), res = 300, width = 1500, height = 1500)
    #   g <- ggplot(data = tmp, aes(x = label, y = score, group = label, fill = as.factor(label))) +
    #     geom_boxplot(outlier.shape = 1) +
    #     theme_bw() +
    #     theme(legend.title=element_blank(),legend.position=c(0.15,0.85),axis.text.x = element_blank(),axis.ticks.x=element_blank())+
    #     scale_fill_manual(breaks=as.factor(c(1,0)),labels=c("Case","Control"),values = c("seagreen", "violetred"))+
    #     labs(x="", y = "Predict score")
    #   print(g)
    #   dev.off()
    # }
    return(fit)
  } else if (mode == 'continuous'){
    ## targetting case, regression
    # train set (get from server)
    train_id <- c("E17", "E24", "E16", "E20", "E1",  "E26", "E7",  "E30", "E11", "E19", "E31", "E23", "E29", 
                  "E3",  "E5",  "E25", "E6",  "E21", "E9",  "E14", "E13")
    train_mat <- t(input_mat[,train_id])
    colnames(train_mat) <- as.character(input_mat$Compound)
    # test set
    test_id <- setdiff(total_id,train_id)
    test_mat <- t(input_mat[,test_id])
    colnames(test_mat) <- as.character(input_mat$Compound)
    # response
    bp_train <- as.numeric(blood$blood_pressure[match(rownames(train_mat),blood$X)])
    names(bp_train) <- rownames(train_mat)
    bp_test <- as.numeric(blood$blood_pressure[match(rownames(test_mat),blood$X)])
    names(bp_test) <- rownames(test_mat)
    # modeling
    gm <- glmnet(x = train_mat, y = bp_train)
    # gm.cv <- cv.glmnet(x = train_mat, y = bp_train)
    # gm.cv$lambda.min
    # pred.cv.test <- predict(gm.cv, newx = test_mat, s = "lambda.min")
    # plot(gm.cv)
    # prediction
    pred.train=predict(object = gm, newx = train_mat,s = 2.4)
    pred.test=predict(object = gm, newx = test_mat,s = 2.4)
    
    # mean((pred.train - bp_train)^2)
    # mean((pred.test - bp_test)^2)
    mean((pred.cv.test - bp_test)^2)
    tmp.path <- selected_path_data$path.name[match(pathway_id,selected_path_data$path.id)]
    # plot
    # train
    fit_train_plot <- lm(pred.train~bp_train)
    pred.train.95 <- predict(fit_train_plot,data.frame(bp_train),interval='confidence')
    pred.train.95=data.frame(bp = bp_train,pred=pred.train.95[,1],upper=pred.train.95[,3],lower=pred.train.95[,2])
    tmp.lim.x <- c(min(bp_train)-2,max(bp_train)+2)
    tmp.lim.y <- c(min(pred.train)-2,max(pred.train)+2)
    
    png(file = paste(tmp.path,'train.png'),height = 800, width = 850)
    plot(bp_train,pred.train,pch=20,cex = 3.5, cex.axis=2,cex.lab=2,font.axis=2,cex.main = 3,font.lab = 2,
         xlab='Blood presure',ylab='Predicted blood presure',xlim = tmp.lim.x,ylim = tmp.lim.y,
         main = tmp.path)
    lines(pred.train.95$bp,pred.train.95$pred,col='red',lwd = 3)
    lines(lowess(pred.train.95$bp,pred.train.95$upper), col="red", lwd = 3)
    lines(lowess(pred.train.95$bp,pred.train.95$lower), col="red", lwd = 3)
    box(lwd=2)
    dev.off()
    
    tmp <- summary(fit_train_plot)
    train.r2 <- tmp$r.squared
    train.p <- tmp$coefficients[8]
    
    # testing
    fit_test_plot <- lm(pred.test~bp_test)
    pred.test.95 <- predict(fit_test_plot,data.frame(bp_test),interval='confidence')
    pred.test.95=data.frame(bp = bp_test,pred=pred.test.95[,1],upper=pred.test.95[,3],lower=pred.test.95[,2])
    tmp.lim.x <- c(min(bp_test)-2,max(bp_test)+2)
    tmp.lim.y <- c(min(pred.test)-2,max(pred.test)+2)
    
    png(file = paste(tmp.path,'test.png'),height = 800, width = 850)
    plot(bp_test,pred.test,pch=20,cex = 3.5, cex.axis=2,cex.lab=2,font.axis=2,cex.main = 3,font.lab = 2,
         xlab='Blood presure',ylab='Predicted blood presure',xlim = tmp.lim.x,ylim = tmp.lim.y,
         main = tmp.path)
    lines(pred.test.95$bp,pred.test.95$pred,col='red',lwd = 2)
    lines(lowess(pred.test.95$bp,pred.test.95$upper), col="red", lwd = 2)
    lines(lowess(pred.test.95$bp,pred.test.95$lower), col="red", lwd = 2)
    box(lwd=2)
    dev.off()
    
    tmp <- summary(fit_test_plot)
    test.r2 <- tmp$r.squared
    test.p <- tmp$coefficients[8]
    
    stats <- c(train.r2,train.p,
               test.r2,test.p)
    names(stats) <- c('train.r2','train.p','test.r2','test.p')
    return(as.list(stats,gm))
    
  } else if (mode == 'synergy'){
    # ---- pathway sinergy ----
    train_id <- clinic$ID[which(clinic$label == 'train')]
    train_mat <- t(input_mat[,train_id])
    colnames(train_mat) <- as.character(input_mat$Compound)
    train_class <- as.factor(substr(rownames(train_mat),1,1))
    # modeling with random forest
    fit <- randomForest(train_mat, y=train_class, ntree=1000,
                        replace=F, importance=TRUE)
    # # prediction part
    
    # all_mat <- t(input_mat[,2:dim(input_mat)[2]])
    # colnames(all_mat) <- input_mat$Compound
    # all_pred <- predict(fit,all_mat,type = "prob")
    # all_score <- data.frame(score = all_pred[,1])
    # all_score$label <- c(rep(1,length(grep('E',rownames(all_score)))),rep(0,length(grep('P',rownames(all_score)))))
    pred.id <- clinic$ID[which(clinic$GA < 34 | clinic$Status2 == 'severe' | clinic$Status2 == 'normal')]
    pred.mat <- t(input_mat[,pred.id])
    colnames(pred.mat) <- as.character(input_mat$Compound)
    all_pred <- predict(fit,pred.mat,type = "prob")
    all_score <- data.frame(score = all_pred[,1])
    all_score$label <- c(rep(1,length(grep('E',rownames(all_score)))),rep(0,length(grep('P',rownames(all_score)))))
    
    
    # early.id <- clinic$ID[which(clinic$GA < 34)]
    # tmp.early <- all_mat[early.id,]
    # early_pred <- predict(fit,tmp.early,type = "prob")
    # early_score <- data.frame(score.early = early_pred[,1])
    # 
    # severe.id <- clinic$ID[which(clinic$Status2 == 'severe' | clinic$Status2 == 'normal')]
    # tmp.severe <- all_mat[severe.id,]
    # severe_pred <- predict(fit,tmp.severe,type = "prob")
    # severe_score <- data.frame(score.severe = severe_pred[,1])
    # 
    # tmp <- merge(all_score,early_score,by = 0, all.x = T)
    # rownames(tmp) <- tmp$Row.names
    # tmp <- merge(tmp,severe_score,by = 0, all.x = T)
    # tmp <- tmp[,3:6]
    # stats <- tmp
    return(all_score)
  }
  
}

# edit pathway notation (add total pathways)
selected_path_data[23,] <- NA
selected_path_data$path.id[23] <- 'Total.path'
selected_path_data$path.name[23] <- 'All pathways catelog'

# ---- testing the train set ----

# tmp.id <- list()
# for (i in 1:5000){
#   train_id <- c(sample(clinic$ID[1:32],21),sample(clinic$ID[33:64],21))
#   train_mat <- t(input_mat[,train_id])
#   colnames(train_mat) <- as.character(input_mat$Compound)
#   train_class <- as.factor(substr(rownames(train_mat),1,1))
#   # modeling with EN
#   # fit <- glmnet(x = train_mat, y = train_class, family = 'binomial')
#   cvfit = cv.glmnet(x= train_mat, y= train_class, family = "binomial", type.measure = "class")
#   tmp <- coef(cvfit, s = "lambda.min")
#   tmp2 <- tmp@Dimnames[[1]][tmp@i+1]
#   tmp.id[[i]] <- tmp2
# }
# model.binary.analyte <- data.frame(cpd = tmp2,coeff = tmp@x,stringsAsFactors = F)

# ----- draw binary outcome for total pathways -----
# binary outcome
dump <- multi_task(conc_PE_z_flt, pathway_id = 'Total.path', clinic = clinic, blood = blood, 'binary')
# elastic regression outcome
tmp.ga <- as.numeric(clinic$GA)
names(tmp.ga) <- as.character(clinic$ID)
tmp.ga <- data.frame(t(tmp.ga))
tmp.ga$Compound <- 'gestational_age'
tmp.ga <- tmp.ga[colnames(conc_E_z_flt)]
input_mat <- rbind(tmp.ga,conc_E_z_flt)
rownames(input_mat) <- NULL
total.stats <- multi_task(input_mat, pathway_id = 'Total.path', clinic = clinic, blood = blood, 'continuous') 
# synergy outcome
all_path_score <- multi_task(conc_PE_z_flt, pathway_id = 'Total.path', clinic = clinic, blood = blood, 'syndergy')

# ----- draw regression or binary model selected analytes -----
# # set color and size
# color

tmp2 <- tmp@Dimnames[[1]][tmp@i+1]
model.analyte <- data.frame(cpd = tmp2[2:17],
                            coeff = tmp@x[2:17],stringsAsFactors = F)

model.continuous.analyte <- model.analyte
model.analyte <- model.binary.analyte

sig_up <- model.analyte$cpd[model.analyte$coeff>0]
sig_dn <- model.analyte$cpd[model.analyte$coeff<0]

sig_up_1 <- model.analyte$cpd[model.analyte$coeff > 0 & model.analyte$coeff < 0.2]
sig_up_2 <- model.analyte$cpd[model.analyte$coeff >= 0.2 & model.analyte$coeff < 0.5]
sig_up_3 <- model.analyte$cpd[model.analyte$coeff >= 0.5]

sig_dn_1 <- model.analyte$cpd[model.analyte$coeff < 0 & model.analyte$coeff > -0.2]
sig_dn_2 <- model.analyte$cpd[model.analyte$coeff <= -0.2 & model.analyte$coeff > -0.5]
sig_dn_3 <- model.analyte$cpd[model.analyte$coeff <= -0.5]

sig_large <- model.analyte$cpd[abs(model.analyte$coeff) > 0.5]
sig_mid <- model.analyte$cpd[abs(model.analyte$coeff) <= 0.5 & abs(model.analyte$coeff) > 0.2]
sig_small <- model.analyte$cpd[abs(model.analyte$coeff) <= 0.2]

layout_2d$label <- NA
layout_2d$label[match(sig_up,layout_2d$cpd)] <- 'up'
layout_2d$label[match(sig_dn,layout_2d$cpd)] <- 'down'
layout_2d$label[is.na(layout_2d$label)] <- 'non'
rm(sig_up,sig_dn)

layout_2d$color <- NA
layout_2d$color[match(sig_up_1,layout_2d$cpd)] <- '#FF8F8F'
layout_2d$color[match(sig_up_2,layout_2d$cpd)] <- '#FF4B4B'
layout_2d$color[match(sig_up_3,layout_2d$cpd)] <- '#FF0000'

layout_2d$color[match(sig_dn_1,layout_2d$cpd)] <- '#939DFF'
layout_2d$color[match(sig_dn_2,layout_2d$cpd)] <- '#5463FF'
layout_2d$color[match(sig_dn_3,layout_2d$cpd)] <- '#0017FF'

layout_2d$color[which(layout_2d$label=='non')] <- '#a8a8a8'

# label compound size according to pvalue
layout_2d$size <- NA
layout_2d$size[which(layout_2d$label == 'non')] <- 9

layout_2d$size[match(sig_large,layout_2d$cpd)] <- 19
layout_2d$size[match(sig_mid,layout_2d$cpd)] <- 15
layout_2d$size[match(sig_small,layout_2d$cpd)] <- 11

# grey
layout_sub_1 <- layout_2d[which(layout_2d$label=="non"),]
# color
layout_sub_2 <- layout_2d[which(layout_2d$label!="non"),]

p <- plot_ly() %>%
  
  add_segments(
    data = edge_data_frame,
    x = ~x_from, xend = ~x_to,
    y = ~y_from, yend = ~y_to,
    alpha = 0.2, size = I(1), color = I("grey")
  ) %>%
  
  add_trace(data = layout_sub_1, x = ~PC1, y = ~PC2,
            marker = list(size = ~size,
                          color = ~color)) %>%
  add_trace(data = layout_sub_2, x = ~PC1, y = ~PC2,
            marker = list(size = ~size,
                          color = ~color)) %>%
  layout(title = '',
         # shapes = edge_shapes,
         xaxis = axis,
         yaxis = axis) 

plotly_IMAGE(p, format = "png", height = 1000,width = 1000,
             out_file = 'model.binary.sig.cpd.png')

# ----- draw for each pathway -----
path_final <- c('map00071',
                'gly_PG',
                'gly_PE',
                'gly_PA',
                'map01210',
                'map00100',
                'map01230',
                'map01212',
                'map00600')

# binary outcome
for (i in 1:length(path_final)){
  cpd <- pathways_set[[which(names(pathways_set) == path_final[i])]]
  cpd <- intersect(cpd,conc_PE_z_flt$Compound)
  input_mat <- conc_PE_z_flt[match(cpd,conc_PE_z_flt$Compound),]
  dump <- multi_task(input_mat, pathway_id = path_final[i], clinic = clinic, blood = blood,mode = 'binary')
}

# elastic regression
pathway_regression_stats <- list()
for (i in 1:length(path_final)){
  cpd <- pathways_set[[which(names(pathways_set) == path_final[i])]]
  cpd <- intersect(cpd,as.character(conc_E_z_flt$Compound))
  input_mat_2 <- input_mat[match(cpd,input_mat$Compound),]
  pathway_regression_stats[[i]] <- multi_task(input_mat_2, pathway_id = path_final[i], clinic = clinic, blood = blood,mode = 'continuous')
  names(pathway_regression_stats)[i] <- path_final[i]
}

# synergy outcome
pathway_synergy_stats <- list()

for (i in 1:length(path_final)){
  cpd <- pathways_set[[which(names(pathways_set) == path_final[i])]]
  cpd <- intersect(cpd,as.character(conc_PE_z_flt$Compound))
  input_mat <- conc_PE_z_flt[match(cpd,conc_PE_z_flt$Compound),]
  pathway_synergy_stats[[i]] <- multi_task(input_mat, pathway_id = path_final[i], clinic = clinic, blood = blood,mode = 'synergy')
  names(pathway_synergy_stats)[i] <- path_final[i]
}
# plot for synergy outcome
tmp.length = choose(9,2)
path.comb <- combn(path_final,2)
for (i in 1:tmp.length){
  tmp.path1 <- path.comb[1,i]
  tmp.path2 <- path.comb[2,i]
  
  tmp.path1.data <- pathway_synergy_stats[[which(names(pathway_synergy_stats)==tmp.path1)]]
  # tmp.path1.early <- tmp.path1.data[!is.na(tmp.path1.data$score.early),]
  # tmp.path1.severe <- tmp.path1.data[!is.na(tmp.path1.data$score.severe),]
  # colnames(tmp.path1.early)[3] <- paste('early.',tmp.path1,sep = '')
  # colnames(tmp.path1.severe)[4] <- paste('severe.',tmp.path1,sep = '')
  colnames(tmp.path1.data)[1] <- paste('score.',tmp.path1,sep = '')
  
  tmp.path2.data <- pathway_synergy_stats[[which(names(pathway_synergy_stats)==tmp.path2)]]
  # tmp.path2.early <- tmp.path2.data[!is.na(tmp.path2.data$score.early),]
  # tmp.path2.severe <- tmp.path2.data[!is.na(tmp.path2.data$score.severe),]
  # colnames(tmp.path2.early)[3] <- paste('early.',tmp.path2,sep = '')
  # colnames(tmp.path2.severe)[4] <- paste('severe.',tmp.path2,sep = '')
  colnames(tmp.path2.data)[1] <- paste('score.',tmp.path2,sep = '')
  
  # data.early <- cbind(tmp.path1.early[3],tmp.path2.early[3],tmp.path2.early[2])
  # data.severe <- cbind(tmp.path1.severe[4],tmp.path2.severe[4],tmp.path2.severe[2])
  data.total <- cbind(tmp.path1.data[1],tmp.path2.data[1],tmp.path2.data[2])
  
  # data.early$label[grep(1,data.early$label)] <- 'Case'
  # data.early$label[grep(0,data.early$label)] <- 'Control'
  # data.early$label <- as.factor(data.early$label)
  # data.severe$label[grep(1,data.severe$label)] <- 'Case'
  # data.severe$label[grep(0,data.severe$label)] <- 'Control'
  # data.severe$label <- as.factor(data.severe$label)
  data.total$label2 <- NA
  tmp.early.id <- clinic$ID[which(clinic$GA < 34)]
  tmp.control.id <- clinic$ID[which(clinic$Status2 == 'normal')]
  tmp.severe.id <- clinic$ID[which(clinic$Status2 == 'severe')]
  data.total$label2[match(tmp.early.id,rownames(data.total))] <- 'Early PE'
  data.total$label2[match(tmp.control.id,rownames(data.total))] <- 'NP'
  data.total$label2[match(tmp.severe.id,rownames(data.total))] <- 'Severe PE'
  
  tmp.x <- selected_path_data$path.name[match(tmp.path1,selected_path_data$path.id)]
  tmp.y <- selected_path_data$path.name[match(tmp.path2,selected_path_data$path.id)]
  p <- ggplot(data = data.total, 
              aes(x = data.total[,1], y = data.total[,2], color = label2)) +
    geom_point() + 
    theme_bw() +
    xlim(-0.12,1.12) + 
    ylim(-0.12,1.12) +
    scale_color_manual(values = c("#C71585","#2E8B57","#75054b")) +
    scale_fill_manual(breaks=as.factor(c(0,1,2)),
                      labels = c("Early PE","Severe PE","NP")) +
    labs(x = tmp.x, 
         y = tmp.y)
  png(paste(i,'.png',sep = ''),height = 800,width = 1000,res = 200)
  print(p)
  dev.off()
  
  # p <- ggplot(data = data.severe, 
  #             aes(x = data.severe[,1], y = data.severe[,2], color = label)) +
  #   geom_point() + 
  #   theme_bw() +
  #   
  #   scale_color_manual(values = c("#C71585","#2E8B57")) +
  #   scale_fill_manual(breaks=as.factor(c(0,1)),
  #                     labels = c("case","control")) +
  #   labs(x = tmp.x, 
  #        y = tmp.y)
  # png(paste(i,'severe.png',sep = ''),height = 800,width = 1000,res = 200)
  # print(p)
  # dev.off()
}

# # ----- build training set -----
# training.ids <- clinic$ID[which(clinic$label == 'train')]
# training.data <- conc_PE_z_flt[,training.ids]
# training.data <- t(training.data)
# colnames(training.data) <- as.character(conc_PE_z_flt$Compound)
# training.class <- as.factor(substr(rownames(training.data),1,1))
# 
# # run model: random forest
# fit <- randomForest(training.data, y=training.class, ntree=1000,
#                    replace=F, importance=TRUE)
# 

# # ----- prediction on training set -----
# # total train
# all_pred <- predict(fit,training.data,type = "prob")
# train.score.all <- data.frame(score = all_pred[,1])
# train.score.all$label <- c(rep(1,length(grep('E',rownames(train.score.all)))),rep(0,length(grep('P',rownames(train.score.all)))))
# train.score.all$blood_presure <- blood$blood_pressure[match(rownames(train.score.all),blood$X)]
# # early train
# train.early.id <- clinic$ID[which(clinic$GA <= 34 & clinic$label == 'train')]
# train.score.early <- train.score.all[match(train.early.id,rownames(train.score.all)),]
# # severe train
# train.severe.id <- clinic$ID[which((clinic$Status2 == 'severe' | clinic$Status2 == 'normal') & 
#                                      clinic$label == 'train')]
# train.score.severe <- train.score.all[match(train.severe.id,rownames(train.score.all)),]

# # ----- whisker plot for training set -----
# var <- train.score.severe
# name <- 'train.score.severe'
# png(file = paste(name,".png",sep = ""), res = 300, width = 1500, height = 1500)
# ggplot(data = var, aes(x = label, y = score, group = label, fill = as.factor(label))) +
#   geom_boxplot(outlier.shape = T) + 
#   theme_bw() +
#   theme(legend.title=element_blank(),legend.position=c(0.15,0.85),axis.text.x = element_blank(),axis.ticks.x=element_blank())+
#   scale_fill_manual(breaks=as.factor(c(1,0)),labels=c("Case","Control"),values = c("seagreen", "violetred"))+
#   labs(x="", y = "Predict score")
# dev.off()

# # ----- build testing set -----
# testing.ids <- clinic$ID[which(clinic$label == 'testing')]
# testing.data <- conc_PE_z_flt[,testing.ids]
# testing.data <- t(testing.data)
# colnames(testing.data) <- as.character(conc_PE_z_flt$Compound)
# testing.class <- as.factor(substr(rownames(testing.data),1,1))

# # ----- prediction on testing set -----
# # total testing
# all_pred_test <- predict(fit,testing.data,type = "prob")
# test.score.all <- data.frame(score = all_pred_test[,1])
# test.score.all$label <- c(rep(1,length(grep('E',rownames(test.score.all)))),rep(0,length(grep('P',rownames(test.score.all)))))
# test.score.all$blood_presure <- blood$blood_pressure[match(rownames(test.score.all),blood$X)]
# # early test
# test.early.id <- clinic$ID[which(clinic$GA <= 34 & clinic$label == 'testing')]
# test.score.early <- test.score.all[match(test.early.id,rownames(test.score.all)),]
# # severe test
# test.severe.id <- clinic$ID[which((clinic$Status2 == 'severe' | clinic$Status2 == 'normal') & 
#                                      clinic$label == 'testing')]
# test.score.severe <- test.score.all[match(test.severe.id,rownames(test.score.all)),]

# # ----- whisker plot for testing set -----
# var <- test.score.severe
# name <- 'test.score.severe'
# png(file = paste(name,".png",sep = ""), res = 300, width = 1500, height = 1500)
# ggplot(data = var, aes(x = label, y = score, group = label, fill = as.factor(label))) +
#   geom_boxplot(outlier.shape = T) + 
#   theme_bw() +
#   theme(legend.title=element_blank(),legend.position=c(0.15,0.85),axis.text.x = element_blank(),axis.ticks.x=element_blank())+
#   scale_fill_manual(breaks=as.factor(c(1,0)),labels=c("Case","Control"),values = c("seagreen", "violetred"))+
#   labs(x="", y = "Predict score")
# dev.off()

# # some basic stats
# varImpPlot(fit)
# var <- test.score.severe
# t.test(var$score[which(var$label==1)],var$score[which(var$label==0)])

# # ----- run model: logistic regression -----
# model <- glm(training.class~.,family=binomial(link='logit'),data = data.frame(training.data))
# 
# library(ggplot2)
# # whisker plot
# ggplot(score, aes(x=score, y=label)) +
#   geom_point(color = c(rep('blue',length(grep('E',rownames(score)))),rep('red',length(grep('P',rownames(score))))),
#              size = 2)
# 
# # blood presure vs score
# score_2 <- score[1:6,]
# ggplot(score_2, aes(y=score, x=blood_presure)) +
#   geom_point(color = rep('red',length(grep('E',rownames(score)))),
#              size = 2)

# # ----- blood pressure linear regression -----
# # build training dataset: early
# training.ids <- as.character(blood$X)
# training.data <- conc_PE_z_flt[,training.ids]
# training.data <- t(training.data)
# colnames(training.data) <- as.character(conc_PE_z_flt$Compound)
# training.class <- as.numeric(blood$blood_pressure)
# training.data <- data.frame(training.data)
# training.data <- cbind(training.class,training.data)

# fit <- lm(training.class ~ C0.Carnitine, data=training.data)
# summary(fit)
# 
# fit <- lm(training.data)
# tmp <- predict(fit,training.data)



# ----- REFERENCE: ELASTIC net and variance analysis -----

foldid=as.integer(as.factor(input.y.train$subject[input.y.train$GA<100]))
gm<-cv.glmnet(input.x.train.reduce,input.y.train.reduce$GA, alpha=0.2, standardize=FALSE,keep=TRUE,foldid=foldid)

load("sam_list.rdata")

pred.test=as.vector(predict(sam_list$gm, sam_list$input.x.test, s='lambda.min'))
pred.train=as.vector(predict(sam_list$gm, sam_list$input.x.train, s='lambda.min'))

# plot linear regression
#train
fit=lm(pred.train~sam_list$input.y.train$GA)

lm.pred.95=predict(fit,data.frame(sam_list$input.y.train$GA),interval='confidence')

lm.plot.data=data.frame(GA=sam_list$input.y.train$GA,pred=lm.pred.95[,1],upper=lm.pred.95[,3],lower=lm.pred.95[,2])
lm.plot.data=lm.plot.data[order(lm.plot.data$GA),]

plot(sam_list$input.y.train$GA,pred.train,pch=20,cex.axis=1.2,cex.lab=1.2,font.axis=2,xlab='Gestational age (weeks)',ylab='Model prediction')


lines(lm.plot.data$GA,lm.plot.data$pred,col='red')
lines(lm.plot.data$GA,lm.plot.data$upper,col='red')
lines(lm.plot.data$GA,lm.plot.data$lower,col='red')

box(lwd=2)

#test
fit=lm(pred.test~sam_list$input.y.test$GA)

lm.pred.95=predict(fit,data.frame(sam_list$input.y.test$GA),interval='confidence')

lm.plot.data=data.frame(GA=sam_list$input.y.test$GA,pred=lm.pred.95[,1],upper=lm.pred.95[,3],lower=lm.pred.95[,2])
lm.plot.data=lm.plot.data[order(lm.plot.data$GA),]

plot(sam_list$input.y.test$GA,pred.test,pch=20,cex.axis=1.2,cex.lab=1.2,font.axis=2,xlab='Gestational age (weeks)',ylab='Model prediction')


lines(lm.plot.data$GA,lm.plot.data$pred,col='red')
lines(lm.plot.data$GA,lm.plot.data$upper,col='red')
lines(lm.plot.data$GA,lm.plot.data$lower,col='red')

box(lwd=2)


