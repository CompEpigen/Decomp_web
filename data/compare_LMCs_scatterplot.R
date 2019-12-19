library(MeDeCom)
load("FactorViz_outputs/medecom_set.RData")
load("FactorViz_outputs/ann_C.RData")
K <- 7
lambda <- 0.001
sel.lmc <- 4

lmcs <- getLMCs(medecom.set,K=K,lambda=lambda)
diffs <- lmcs[,sel.lmc]-rowMedians(lmcs[,-sel.lmc])
to.plot <- data.frame(Difference=diffs)
plot <- ggplot(to.plot,aes(x=Difference,y=..count..))+geom_histogram(fill="white",color="black",binwidth = 0.05)+
  geom_vline(xintercept = 0.75)+geom_vline(xintercept = -0.75)+xlim(-1,1)+
  theme_bw()+theme(panel.grid=element_blank())

to.plot <- data.frame(Sel.LMC=lmcs[,sel.lmc],Median.others=rowMedians(lmcs[,-sel.lmc]),Diff=diffs)
plot <- ggplot(to.plot,aes(x=Sel.LMC,y=Median.others))+geom_point(shape=ifelse(abs(to.plot$Diff)>0.75,17,16),color=ifelse(to.plot$Diff<(-0.5),"firebrick3","black"),
 size=ifelse(abs(to.plot$Diff)>0.75,2,1)
                                                                  )+
  geom_abline(slope=1,intercept=0)+xlim(0,1)+ylim(0,1)+
  theme_bw()+theme(panel.grid = element_blank())

ann.C <- ann.C[medecom.set@parameters$GROUP_LISTS[[1]],]
sel.sites <- ann.C[abs(diffs)>0.75,]
sel.sites$Difference <- diffs[abs(diffs)>0.75]
sel.sites <- sel.sites[order(abs(sel.sites$Difference),decreasing = T),]
sel.sites.gr <- makeGRangesFromDataFrame(sel.sites)
genes <- unlist(rnb.get.annotation("genes"))
dists.all <- c()
for(i in 1:length(sel.sites.gr)){
  gr <- sel.sites.gr[i]
  dists <- distance(gr,genes)
  dists.all <- rbind(dists.all,c(which.min(dists),min(dists,na.rm = T)))
}
sel.sites$Nearest.gene.symbol <- values(genes)$symbol[dists.all[,1]]
sel.sites$Nearest.gene.ENSBML <- row.names(values(genes))[dists.all[,1]]
sel.sites$Nearest.gene.distance <- dists.all[,2]
rnb.load.annotation.from.db("ensembleRegBuildBPall")
ens.build <- unlist(rnb.get.annotation("ensembleRegBuildBPall"))
op <- findOverlaps(sel.sites.gr,ens.build,select = "first")
op[!is.na(op)] <- values(ens.build)$elementType[op[!is.na(op)]]
sel.sites$EnsemblAnnotation <- op
sel.sites$Nearest.gene.ENSBML <- unlist(lapply(strsplit(sel.sites$Nearest.gene.ENSBML,"[[:punct:]]"),function(x)x[2]))
