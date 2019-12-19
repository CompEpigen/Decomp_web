library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-LUAD",
                           data.category = "Gene expression",
                           data.type = "Gene expression quantification",
                           platform = "Illumina HiSeq", 
                           experimental.strategy = "RNA-Seq",
							file.type  = "normalized_results",
                           legacy = TRUE)
GDCdownload(query, method = "api", files.per.chunk = 10)
data <- GDCprepare(query,summarizedExperiment=F)
obj <- DGEList(data)
row.names(obj$samples) <- unlist(lapply(strsplit(row.names(obj$samples),"_"),function(x)x[3]))
colnames(obj$counts) <- unlist(lapply(strsplit(colnames(obj$counts),"_"),function(x)x[3]))
row.names(obj$samples) <- substr(row.names(obj$samples),1,16)
colnames(obj$counts) <- substr(colnames(obj$counts),1,16)
cpm.obj <- cpm(obj)

load("TCGA_LUAD/FactorViz_outputs/medecom_set.RData")
props <- getProportions(medecom.set,K=7,lambda=0.001)
load("TCGA_LUAD/FactorViz_outputs/ann_S.RData")
colnames(props) <- substr(ann.S$Comment..TCGA.Barcode.,1,16)
marker.genes <- c("EPCAM","CLDN5","COL1A2","PTPRC")
in.exp <- colnames(cpm.obj) %in% colnames(props)
in.props <- colnames(props) %in% colnames(cpm.obj)
props <- props[,in.props]
cpm.obj <- cpm.obj[,in.exp]
cpm.obj <- cpm.obj[,colnames(props)]
row.names(cpm.obj) <- unlist(lapply(strsplit(row.names(cpm.obj),"[[:punct:]]"),function(x)x[1]))
cors.all <- sapply(marker.genes,function(marker){
	if(!marker %in% row.names(cpm.obj)){
		cors.gene <- NA
	}else{
		sel.exp <- cpm.obj[marker,]
		cors.gene <- apply(props,1,function(prop){
			cor(unlist(sel.exp),unlist(prop))
		})
	}
	cors.gene
})

cors.p.vals <- sapply(marker.genes,function(marker){
	if(!marker %in% row.names(cpm.obj)){
		cors.gene <- NA
	}else{
		sel.exp <- cpm.obj[marker,]
		cors.gene <- apply(props,1,function(prop){
			cor.test(unlist(sel.exp),unlist(prop))$p.value
		})
	}
	cors.gene
})
library(corrplot)
corrplot(cors.all,"ellipse")

plot.path <- "analysis/gene_expression/"
cors.all <- sapply(marker.genes,function(marker){
	if(!marker %in% row.names(cpm.obj)){
		cors.gene <- NA
	}else{
		sel.exp <- cpm.obj[marker,]
		for(j in 1:nrow(props)){
			prop <- props[j,]
			lmc <- paste0("LMC",j)
			to.plot <- data.frame(CPM=sel.exp,Proportion=prop)
			plot <- ggplot(to.plot,aes(x=Proportion,y=CPM))+geom_point(size=.1)+geom_smooth(method="lm",size=.5)+
				theme_bw()+theme(panel.grid=element_blank(),text=element_text(color="black",size=20),
								axis.ticks=element_line(size=0.5,color="black"),axis.ticks.length=unit(2,"mm"),axis.title=element_blank(),axis.text=element_blank())
			ggsave(file.path(plot.path,paste0(lmc,"_",marker,"_new.pdf")),plot,width=35,height=35,unit="mm")
		}
	}
})

