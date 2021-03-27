
# PCE cluster 2
inputGenes<-toupper(markers.filtered.names.2)
humanGene<-humanGeneMapping[humanGeneMapping$Gene.name %in% inputGenes,]
inputGenes<-humanGene$Gene
expressionData<-data[intersect(row.names(data),humanGeneMapping$Gene),]
se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),rowData = row.names(expressionData),colData = colnames(expressionData))
cellSpecificGenesExp<-teGeneRetrieval(se,expressedGeneThreshold = 1)
se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),rowData = row.names(expressionData),colData = colnames(expressionData))
cellSpecificGenesExp<-teGeneRetrieval(se,expressedGeneThreshold = 1)
print(length(inputGenes))
gs<-GeneSet(geneIds=toupper(inputGenes))
output2<-teEnrichmentCustom(gs,cellSpecificGenesExp)
enrichmentOutput<-setNames(data.frame(assay(output2[[1]]),row.names = rowData(output2[[1]])[,1]),colData(output2[[1]])[,1])
row.names(cellDetails)<-cellDetails$RName
enrichmentOutput$Tissue<- cellDetails[row.names(enrichmentOutput),"CellName"]
p2<- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 2") + ylab("-log10 p-value")
ggsave("Cluster_2.svg", dpi=700, width = 10, height = 8)
# PCE cluster 3
inputGenes<-toupper(markers.filtered.names.3)
humanGene<-humanGeneMapping[humanGeneMapping$Gene.name %in% inputGenes,]
inputGenes<-humanGene$Gene
expressionData<-data[intersect(row.names(data),humanGeneMapping$Gene),]
se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),rowData = row.names(expressionData),colData = colnames(expressionData))
cellSpecificGenesExp<-teGeneRetrieval(se,expressedGeneThreshold = 1)
se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),rowData = row.names(expressionData),colData = colnames(expressionData))
cellSpecificGenesExp<-teGeneRetrieval(se,expressedGeneThreshold = 1)
print(length(inputGenes))
gs<-GeneSet(geneIds=toupper(inputGenes))
output2<-teEnrichmentCustom(gs,cellSpecificGenesExp)
enrichmentOutput<-setNames(data.frame(assay(output2[[1]]),row.names = rowData(output2[[1]])[,1]),colData(output2[[1]])[,1])
row.names(cellDetails)<-cellDetails$RName
enrichmentOutput$Tissue<- cellDetails[row.names(enrichmentOutput),"CellName"]
p3<- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 3") + ylab("-log10 p-value")
ggsave("Cluster_3.svg", dpi=700, width = 10, height = 8)
# PCE cluster 5
inputGenes<-toupper(markers.filtered.names.5)
humanGene<-humanGeneMapping[humanGeneMapping$Gene.name %in% inputGenes,]
inputGenes<-humanGene$Gene
expressionData<-data[intersect(row.names(data),humanGeneMapping$Gene),]
se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),rowData = row.names(expressionData),colData = colnames(expressionData))
cellSpecificGenesExp<-teGeneRetrieval(se,expressedGeneThreshold = 1)
se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),rowData = row.names(expressionData),colData = colnames(expressionData))
cellSpecificGenesExp<-teGeneRetrieval(se,expressedGeneThreshold = 1)
print(length(inputGenes))
gs<-GeneSet(geneIds=toupper(inputGenes))
output2<-teEnrichmentCustom(gs,cellSpecificGenesExp)
enrichmentOutput<-setNames(data.frame(assay(output2[[1]]),row.names = rowData(output2[[1]])[,1]),colData(output2[[1]])[,1])
row.names(cellDetails)<-cellDetails$RName
enrichmentOutput$Tissue<- cellDetails[row.names(enrichmentOutput),"CellName"]
p5<- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 5") + ylab("-log10 p-value")
ggsave("Cluster_5.svg", dpi=700, width = 10, height = 8)

sct.cluster5.genes <- as.data.frame(assay(output2[[2]][["SCT"]]))[1]
sct.cluster5.genes <- humanGeneMapping[humanGeneMapping$Gene %in% as.list(sct.cluster5.genes$Gene),]
sct.cluster5.genes <- sct.cluster5.genes$Gene.name
# PCE cluster 6
inputGenes<-toupper(markers.filtered.names.6)
humanGene<-humanGeneMapping[humanGeneMapping$Gene.name %in% inputGenes,]
inputGenes<-humanGene$Gene
expressionData<-data[intersect(row.names(data),humanGeneMapping$Gene),]
se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),rowData = row.names(expressionData),colData = colnames(expressionData))
cellSpecificGenesExp<-teGeneRetrieval(se,expressedGeneThreshold = 1)
se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),rowData = row.names(expressionData),colData = colnames(expressionData))
cellSpecificGenesExp<-teGeneRetrieval(se,expressedGeneThreshold = 1)
print(length(inputGenes))
gs<-GeneSet(geneIds=toupper(inputGenes))
output2<-teEnrichmentCustom(gs,cellSpecificGenesExp)
enrichmentOutput<-setNames(data.frame(assay(output2[[1]]),row.names = rowData(output2[[1]])[,1]),colData(output2[[1]])[,1])
row.names(cellDetails)<-cellDetails$RName
enrichmentOutput$Tissue<- cellDetails[row.names(enrichmentOutput),"CellName"]
p6<- ggplot(data = enrichmentOutput, mapping = aes(x = reorder (Tissue, -Log10PValue), Log10PValue, fill= Tissue)) +
  geom_bar(stat = "identity") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 100), legend.position = "none",
        plot.title = element_text(color = "black", size=18, face="bold.italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ggtitle("Cluster 6") + ylab("-log10 p-value")
ggsave("Cluster_6.svg", dpi=700, width = 10, height = 8)
sct.cluster6.genes <- as.data.frame(assay(output2[[2]][["SCT"]]))[1]
sct.cluster6.genes <- humanGeneMapping[humanGeneMapping$Gene %in% as.list(sct.cluster6.genes$Gene),]
sct.cluster6.genes <- sct.cluster6.genes$Gene.name

fig4 = (p5 | p6) / (p2 | p3)
ggsave("Figure_4_ABCD.svg", plot = fig4, dpi=700, width = 16, height = 18)

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=9)


venn.diagram(
  x = list(sct.cluster5.genes, sct.cluster6.genes),
  category.names = c("SCT genes of Cluster 5" , "SCT genes of Cluster 6"),
  filename = 'Figure_5_A.svg',
  output=TRUE,
  
  # Output features
  imagetype="svg" ,
  height = 6 , 
  width = 6 , 
  resolution = 600,
  
  # Circles
  lwd = 2,
  lty = 'blank',
  col = c(color_list[5], color_list[6]),
  fill = c(alpha("darkturquoise",0.8), alpha('deepskyblue',0.8)),
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-2, 2),
  cat.dist = c(0.02, 0.02),
  cat.fontfamily = "sans"
)
