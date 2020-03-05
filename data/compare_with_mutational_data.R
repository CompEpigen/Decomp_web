library(MeDeCom)
load("TCGA_LUAD/FactorViz_outputs/medecom_set.RData")
load("TCGA_LUAD/FactorViz_outputs/ann_S.RData")
K <- 7
lambda <- 0.001
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

run.trait.association(medecom.set,ann.S)

trait <- "X12p.Status"#"X1p.Status" "X16p.Status" "X16q.Status"
trait <- ann.S[,trait]
props <- getProportions(medecom.set,K=K,lambda=lambda,cg_subset=cg_subs)
to.plot <- data.frame(t(props),Trait=trait)
to.plot <- melt(to.plot,id="Trait")
colnames(to.plot)[2:3] <- c("LMC","Contribution")
plot <- ggplot(to.plot,aes(x=Trait,y=Contribution))+geom_boxplot()+theme_bw()+theme(panel.grid=element_blank())+
  facet_grid(LMC~.)

trait <- "Radiation.Therapy"
trait <- ann.S[,trait]
props <- getProportions(medecom.set,K=K,lambda=lambda,cg_subset=cg_subs)
to.plot <- data.frame(t(props),Trait=trait)
to.plot <- melt(to.plot,id="Trait")
colnames(to.plot)[2:3] <- c("LMC","Contribution")
plot <- ggplot(to.plot,aes(x=Trait,y=Contribution))+geom_boxplot()+theme_bw()+theme(panel.grid=element_blank())+
  facet_grid(LMC~.)

trait <- "Radiation.Therapy"
trait <- ann.S[,trait]
props <- getProportions(medecom.set,K=K,lambda=lambda,cg_subset=cg_subs)
wilcox.test(props[6,which(trait=="Yes")],props[6,which(trait=="No")])#is.na(trait)])

