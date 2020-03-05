library(MeDeCom)
load("TCGA_LUAD/FactorViz_outputs/medecom_set.RData")
load("TCGA_LUAD/FactorViz_outputs/ann_S.RData")
K <- 7
lambda <- 0
cg_subs <- 1
library(data.table)
muts <- fread("annotation/luad_tcga_pan_can_atlas_2018_clinical_data.tsv")
matchi <- match(ann.S$submitter_id,muts$'Patient ID')
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

traits <- c("Fraction.Genome.Altered","Mutation.Count")
for(x in traits){
  ann.S[,x] <- as.numeric(as.character(ann.S[,x]))
}
sel.traits <- c("X1p.Status","X5q.Status","X6q.Status","X12p.Status","X16p.Status","X16q.Status")
for(x in traits){
  foo <- ann.S[,x]
  foo[foo=="Not Called"|foo==""] <- NA
  ann.S[,x] <- foo
}

library(corrplot)
estimate <- read.table("annotation/LUAD_ESTIMATE_RNAseqV2.tab",header=T)
estimate <- as.data.frame(apply(estimate,2,function(x)as.character(x)))
matchi <- match(paste0(ann.S$submitter_id,"-01"),estimate$ID)
new.fr <- matrix(nrow = nrow(ann.S),ncol=4)
for(i in 1:ncol(new.fr)){
  new.fr[!is.na(matchi),i] <- as.character(estimate[matchi[!is.na(matchi)],i])
}
new.fr <- as.data.frame(new.fr)
colnames(new.fr) <- colnames(estimate)
ann.S <- data.frame(ann.S,new.fr)
ann.S$Immune_score <- as.numeric(ann.S$Immune_score)
ann.S$Stromal_score <- as.numeric(ann.S$Stromal_score)
ann.S$ESTIMATE_score <- as.numeric(ann.S$ESTIMATE_score)
sel.traits <- c("Mutation.Count","Fraction.Genome.Altered","LUMP_estimate","Stromal_score")
props <- getProportions(medecom.set,K=7,lambda=lambda)
cors <- apply(props,1,function(x){
  sapply(sel.traits,function(trait){
    trait <- ann.S[,trait]
    na.trait <- is.na(trait)
    #cor(x[!na.trait],as.numeric(trait[!na.trait]))
    cor.test(x[!na.trait],as.numeric(trait[!na.trait]))$p.value
  })
})
corrplot(t(cors),"ellipse",addCoef.col = "black")

eth <- ann.S$ethnicity
eth[eth == "not reported"] <- NA
ann.S$ethnicity <- eth
dis <- rep("cancer",nrow(ann.S))
dis[grepl("11A",ann.S$Comment..TCGA.Barcode.)] <- "healthy"
ann.S$class <- dis
sel.traits <- c("X1p.Status","X5q.Status","X6q.Status","X12p.Status","X16p.Status","X16q.Status")
props <- getProportions(medecom.set,K=7,lambda=0.001)
ann.S$class <- factor(ann.S$class,levels=c("healthy","cancer"))
ann.S$X1p.Status <- factor(ann.S$X1p.Status,levels=c("Gained","Lost"))
ann.S$X5q.Status <- factor(ann.S$X5q.Status,levels=c("Gained","Lost"))
ann.S$X6q.Status <- factor(ann.S$X6q.Status,levels=c("Gained","Lost"))
ann.S$X12p.Status <- factor(ann.S$X12p.Status,levels=c("Gained","Lost"))
ann.S$X16p.Status <- factor(ann.S$X16p.Status,levels=c("Gained","Lost"))
ann.S$X16q.Status <- factor(ann.S$X16q.Status,levels=c("Gained","Lost"))

madiff <- apply(props,1,function(x){
  sapply(sel.traits,function(trait){
    trait <- ann.S[,trait]
    na.trait <- is.na(trait)
    trait <- trait[!na.trait]
    x <- x[!na.trait]
#    md <- mean(x[as.character(trait)==levels(trait)[1]],na.rm=T)-mean(x[as.character(trait)==levels(trait)[2]],na.rm=T)
#	names(md) <- paste0(levels(trait)[1],"vs",levels(trait)[2])
#	md
   t.test(x[as.character(trait)==levels(trait)[1]],na.rm=T,x[as.character(trait)==levels(trait)[2]])$p.value
  })
})
to.plot <- data.frame(t(madiff),LMC=colnames(madiff))
to.plot <- melt(to.plot,id="LMC")
colnames(to.plot)[2:3] <- c("Trait","ADiffM")
mini <- min(madiff)
maxi <- max(madiff)
abs.max <- max(abs(mini),abs(maxi))
plot <- ggplot(to.plot,aes(x=Trait,y=LMC,fill=ADiffM))+geom_tile()+theme_bw()+theme(panel.grid=element_blank())+
  scale_fill_gradient2(low="#6a011f",mid="white",midpoint=0,high="#35ac35",limits=c(-abs.max,abs.max))

