library(pheatmap)
library(MeDeCom)
library(RColorBrewer)
factorviz.output <- "TCGA_LUAD/FactorViz_outputs/"
K <- 7
lambda <- 0.001
cg_subset <- 1
traits <- c("gender","age")
s.id.col <- "geo_accession"
plot.path <- "plot/"
plot.name <- "heatmap_without.pdf"
max.val <- 0.6
min.val <- 0

load(paste0(factorviz.output,"medecom_set.RData"))
contris <- as.data.frame(getProportions(medecom.set,K=K,lambda=lambda,cg_subset=cg_subset))
load(paste0(factorviz.output,"ann_S.RData"))
colnames(contris) <- ann.S[,s.id.col]
if(is.null(max.val)){
  max.val <- max(app(contris))
}
if(is.null(min.val)){
  max.val <- min(rowMins(contris))
}
breaksList <- seq(min.val,max.val,by=0.01)
sel.traits <-ann.S[,traits]
row.names(sel.traits) <- ann.S[,s.id.col]
pdf(file.path(plot.path,plot.name))
pheatmap(contris,annotation_col = sel.traits,cluster_rows = F,clustering_distance_cols = "correlation",show_colnames = F,
        breaks=breaksList,
        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)))
dev.off()

#' TCGA LUAD with cancer/healthy annotated
library(pheatmap)
library(MeDeCom)
library(RColorBrewer)
factorviz.output <- "TCGA_LUAD/FactorViz_outputs/"
K <- 7
lambda <- 0.001
cg_subset <- 1
traits <- c("LUMP_estimate")
s.id.col <- "submitter_id"
plot.path <- "plot/"
plot.name <- "TCGA_contribution_plot.pdf"
max.val <- 1
min.val <- 0

load(paste0(factorviz.output,"medecom_set.RData"))
contris <- as.data.frame(getProportions(medecom.set,K=K,lambda=lambda,cg_subset=cg_subset))
load(paste0(factorviz.output,"ann_S.RData"))
colnames(contris) <- ann.S[,s.id.col]
if(is.null(max.val)){
  max.val <- max(apply(contris,1,max))
}
if(is.null(min.val)){
  min.val <- min(apply(contris,1,min))
}
breaksList <- seq(min.val,max.val,by=0.01)
sel.trait.healthy <- rep("cancer",nrow(ann.S))
sel.trait.healthy[grepl("11A",ann.S$Comment..TCGA.Barcode.)] <- "healthy"
sel.traits <- data.frame(ann.S[,traits],Cancer_Healthy=sel.trait.healthy)
row.names(sel.traits) <- ann.S[,s.id.col]
pdf(file.path(plot.path,plot.name))
pheatmap(contris,annotation_col = sel.traits,cluster_rows = F,,show_colnames = F,
        breaks=breaksList,
        color = colorRampPalette(c("white","black"))(length(breaksList)),
		cluster_distance_cols="euclidean",
		clustering_method="complete")
dev.off()

library(corrplot)
library(MeDeCom)
library(RColorBrewer)
factorviz.output <- "TCGA_LUAD/FactorViz_outputs/"
K <- 7
lambda <- 0.001
cg_subset <- 1
traits <- c("Predicted Male Probability","age","LUMP_estimate")
s.id.col <- "geo_accession"
plot.path <- "plot/"
plot.name <- "heatmap_without.pdf"
max.val <- 0.6
min.val <- 0

load(paste0(factorviz.output,"medecom_set.RData"))
contris <- as.data.frame(getProportions(medecom.set,K=K,lambda=lambda,cg_subset=cg_subset))
load(paste0(factorviz.output,"ann_S.RData"))
colnames(contris) <- ann.S[,s.id.col]
if(is.null(max.val)){
  max.val <- max(app(contris))
}
if(is.null(min.val)){
  max.val <- min(rowMins(contris))
}
breaksList <- seq(min.val,max.val,by=0.01)
sel.traits <-ann.S[,traits]
row.names(sel.traits) <- ann.S[,s.id.col]
res <- apply(contris,1,function(x){
	apply(sel.traits,2,function(y){
		na.any <- is.na(x) | is.na(y)
		cor(x[!na.any],y[!na.any])
	})
})
pdf("ICA_uncorrected_association.pdf")
corrplot(res,method="ellipse")
dev.off()


