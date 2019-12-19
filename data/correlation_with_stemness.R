library(MeDeCom)
load("FactorViz_outputs/medecom_set.RData")
load("FactorViz_outputs/ann_S.RData")
K <- 7
lambda <- 0.001
cg_subs <- 1
library(data.table)
muts <- fread("annotation/stemness_index.csv")
matchi <- match(ann.S$Normalization.Name.1,muts$TCGAlong.id)
muts <- as.data.frame(muts)
new.fr <- matrix(nrow = nrow(ann.S),ncol=ncol(muts))
for(i in 1:ncol(new.fr)){
  new.inf <- as.character(muts[matchi[!is.na(matchi)],i])
  new.inf[new.inf %in% "Not Called"] <- NA
  new.fr[!is.na(matchi),i] <- new.inf
}
new.fr <- as.data.frame(new.fr)
colnames(new.fr) <- colnames(muts)
ann.S <- data.frame(ann.S,new.fr)

ann.S$mDNAsi <- as.numeric(ann.S$mDNAsi)
ann.S$EREG.mDNAsie <- as.numeric(ann.S$EREG.mDNAsi)
ann.S$DMPsi <- as.numeric(ann.S$DMPsi)
ann.S$ENHsi <- as.numeric(ann.S$ENHsi)
sel.traits <- c("ENHsi","EREG.mDNAsie","mDNAsi","DMPsi")
props <- getProportions(medecom.set,K=7,lambda=0.001)
cors <- apply(props,1,function(x){
  sapply(sel.traits,function(trait){
    trait <- ann.S[,trait]
    na.trait <- is.na(trait)
    cor(x[!na.trait],as.numeric(trait[!na.trait]))
    #cor.test(x[!na.trait],as.numeric(trait[!na.trait]))$p.value
  })
})
corrplot(t(cors),"ellipse",addCoef.col = "black")

