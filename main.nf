$HOSTNAME = ""
params.outdir = 'results'  


if (!params.Data_Path){params.Data_Path = ""} 
if (!params.min_UMI){params.min_UMI = ""} 
if (!params.max_UMI){params.max_UMI = ""} 
if (!params.mitoRatio){params.mitoRatio = ""} 
if (!params.varFeatures){params.varFeatures = ""} 
if (!params.numCells){params.numCells = ""} 

Channel.fromPath(params.Data_Path, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_12_outputFileTSV0_g_50}
Channel.value(params.min_UMI).into{g_19_text1_g_20;g_19_text1_g_47}
Channel.value(params.max_UMI).into{g_23_text2_g_20;g_23_text2_g_47}
Channel.value(params.mitoRatio).into{g_24_text3_g_20;g_24_text3_g_47}
Channel.value(params.varFeatures).into{g_25_text4_g_20;g_25_text4_g_47}
Channel.value(params.numCells).set{g_53_text1_g_41}

if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 4000
    $CPU  = 1
    $MEMORY = 400
    $QUEUE = "long"
}
process Merge_Seurat_Objects {

input:
 set val(name), file(data_path) from g_12_outputFileTSV0_g_50

output:
 file "*_seurat_obj.rds"  into g_50_rdsFile00_g_20, g_50_rdsFile00_g_47

"""
#!/usr/bin/env Rscript

# libraries
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(Seurat)

# parse file name
if(grepl(".rds","${data_path}")){
   seurat_data <- readRDS("${data_path}")
} 
if(grepl(".txt|.tsv","${data_path}")){
   seurat_data <- read.table(file = "${data_path}", header = TRUE, row.names = 1, sep = "\\t")
}
if(dir.exists("${data_path}")){
   seurat_data <- Read10X(data.dir = "${data_path}")
}
seurat_obj <- CreateSeuratObject(counts = seurat_data, project = "${name}")

# Create .RData object to load at any time
saveRDS(seurat_obj, file=paste("${name}","_seurat_obj.rds", sep=""))
"""
}


process Seurat_Rmarkdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.rmd$/) "Seurat_Rmd/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.html$/) "Seurat_Rmd_HTML/$filename"}
input:
 file seurat_obj from g_50_rdsFile00_g_47
 val minUMI from g_19_text1_g_47
 val maxUMI from g_23_text2_g_47
 val mitoRatio from g_24_text3_g_47
 val varFeatures from g_25_text4_g_47

output:
 file "*.rmd"  into g_47_rMarkdown00
 file "*.html"  into g_47_outputHTML11

shell:

'''
#!/usr/bin/env perl

my $script = <<'EOF';

---
title: "main_rmarkdown"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Libraries

```{r error=FALSE, message=FALSE, warning=FALSE, cache=FALSE, results = FALSE}
library(Seurat)
library(Matrix)
library(tidyverse)
```

## Read Data

```{r}
# Create Seurat Object for either 10x directory or any other UMI Table
seu <- readRDS("!{seurat_obj}")
```

## Filtering

```{r}
VlnPlot(seu, c("nCount_RNA", "nFeature_RNA"))
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
VlnPlot(seu, c("nCount_RNA", "nFeature_RNA","percent.mt"))
FeatureScatter(seu, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
seu <- subset(seu, subset = nCount_RNA > !{minUMI} & nCount_RNA < !{maxUMI} & percent.mt < !{mitoRatio})
```

## Normalize and Feature Selection

```{r}
seu <- NormalizeData(seu) 
seu <- FindVariableFeatures(seu, nfeatures = !{varFeatures})
top20 <- head(VariableFeatures(seu), 20)
plot_features <- VariableFeaturePlot(seu)
plot_features_label <- LabelPoints(plot = plot_features, points = top20, repel = TRUE)
plot_features_label
```

## Dimensionality Reduction

```{r}
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, verbose = FALSE, npcs = 30)
ElbowPlot(seu, ndims = 30)
DimHeatmap(seu, dims = 1:6, balanced = TRUE)
seu <- RunTSNE(seu, verbose = FALSE, dims = 1:20)
DimPlot(seu, reduction = "tsne")
seu <- RunUMAP(seu, verbose = FALSE, dims = 1:20)
DimPlot(seu, reduction = "umap")
```

## Clustering

```{r}
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.6)
DimPlot(seu, reduction = "tsne", label = T, group.by = "RNA_snn_res.0.6")
```

## Marker Analysis 

```{r}
Idents(seu) <- "RNA_snn_res.0.6"
marker_table_seu <- FindAllMarkers(seu)
marker_table_seu %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(seu, features = top10$gene) + NoLegend()
```

EOF

open OUT, ">rmark.rmd";
print OUT $script;
close OUT;

runCommand("Rscript -e 'rmarkdown::render(\\"rmark.rmd\\",\\"html_document\\", output_file = \\"report.html\\")'   ");

sub runCommand {
            my ($com) = @_;
            my $error = system($com);
            if   ($error) { die "Command failed: $error $com\\n"; }
            else          { print "Command successful: $com\\n"; }
          }

'''


}

if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 400
    $CPU  = 1
    $MEMORY = 100
    $QUEUE = "long"
}
process Filter_Seurat_Object {

input:
 file seurat_obj from g_50_rdsFile00_g_20
 val minUMI from g_19_text1_g_20
 val maxUMI from g_23_text2_g_20
 val mitoRatio from g_24_text3_g_20
 val varFeatures from g_25_text4_g_20

output:
 file "*_filtered_seurat.rds"  into g_20_rdsFile00_g_41

"""
#!/usr/bin/env Rscript

# libraries
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(Seurat)

# read data
seurat_obj <- readRDS("${seurat_obj}")

# Compute percent mito ratio
seurat_obj\$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")

# filter
filtered_seurat <- subset(x = seurat_obj, 
                          subset= (nCount_RNA > ${minUMI}) & 
                                  (nCount_RNA < ${maxUMI}) & 
                                  (mitoRatio < ${mitoRatio}))
                                  
# normalize and find most variable
filtered_seurat <- NormalizeData(filtered_seurat) 
filtered_seurat <- FindVariableFeatures(filtered_seurat, nfeatures = ${varFeatures})

# save Seurat object
saveRDS(filtered_seurat, file=paste(seurat_obj@project.name,"_filtered_seurat.rds", sep = ""))
"""
}

if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 400
    $CPU  = 1
    $MEMORY = 100
    $QUEUE = "long"
}
process Find_Clusters_Seurat {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_seurat.rds$/) "Seurat_Object/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.png$/) "Plots/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tsv$/) "Markers/$filename"}
input:
 file seurat_obj from g_20_rdsFile00_g_41
 val numCells from g_53_text1_g_41

output:
 file "*_seurat.rds"  into g_41_rdsFile00_g_52
 file "*.png" optional true  into g_41_outputFilePng11
 file "*.tsv" optional true  into g_41_txtFile22

"""
#!/usr/bin/env Rscript

# libraries
library(SingleCellExperiment)
library(tidyverse)
library(dplyr)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(Seurat)
library(svglite)

# read data
seurat_obj <- readRDS("${seurat_obj}")

# scale and run PCA
seurat_obj <- ScaleData(seurat_obj)

# terminate if ncol is lower then some number
if(ncol(seurat_obj) > ${numCells}){

	# RunPCA
	seurat_obj <- RunPCA(seurat_obj, npcs = 20)	
	
	# Determine the K-nearest neighbor graph
	seurat_obj <- FindNeighbors(object = seurat_obj, 
	                                   dims = 1:20)
	#res <- c(0.4,0.6,0.8,1.0,1.4)
	res <- c(0.4, 0.6)
	# Determine the clusters for various resolutions                                
	seurat_obj <- FindClusters(object = seurat_obj,
	                                  resolution = res)
	
	# UMAP
	seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, verbose = FALSE)
	
	# TSNE
    seurat_obj <- RunTSNE(seurat_obj, dims = 1:20, verbose = FALSE)

	# Visualize TSNE and UMAP
	for(method in c("umap", "tsne")){
		plot_list <- list()
		for(i in res){
		  plot_list[[paste0(i)]] <- DimPlot(seurat_obj, reduction = method, label = TRUE, pt.size = 1, group.by  = paste0("RNA_snn_res.",i)) + 
			labs(title = paste0(toupper(method), " Clustering (Resolution,", i,")"))
	    }
	    pdf(paste0(seurat_obj@project.name,"_",method,".pdf"), width = 5, height = 5)
	    par(mfrow = c(2,2))
	    print(plot_list)
	    dev.off()
 
		# nCount
		plot_list <- list()
		plot_list[[method]] <- FeaturePlot(seurat_obj, features="nCount_RNA", reduction = method) + 
			labs(title = paste0(toupper(method), " UMI Counts")) + 
			scale_colour_gradientn(colours = CustomPalette(low = "blue", high = "red", mid = "green", k = 100))
		pdf(paste0(seurat_obj@project.name,"_",method,"_umi.pdf"), width = 5, height = 5)
	    par(mfrow = c(2,2))
	    print(plot_list)
	    dev.off()
    }
	
	# marker tables
	markers_res_list <- list()
	for(i in res){
		Idents(seurat_obj) <- paste0("RNA_snn_res.",i)
		try({
			markers <- FindAllMarkers(seurat_obj)
			markers_res_list[[paste0(i)]] <- markers
			top_cluster_markers <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 50)
			markers_top <- FindAllMarkers(seurat_obj, logfc.threshold = 0, min.pct = 0, return.thresh = 1.1)  
			markers_top <- markers_top[markers_top\$gene %in% top_cluster_markers\$gene,]
			splitmarkers <- as_tibble(markers_top) %>% 
			  select(L2FC = avg_log2FC, gene = gene, padj = p_val_adj, cluster = cluster) %>% 
			  group_split(cluster)
			joinmarkers <- as.data.frame(splitmarkers %>% reduce(full_join, by = "gene"))
			rownames(joinmarkers) <- joinmarkers\$gene
			joinmarkers <- joinmarkers[,!grepl("cluster|gene",colnames(joinmarkers))]
			colnames(joinmarkers) <- paste("Clus",
			                         rep(0:(length(splitmarkers)-1),each=2), ".",
			                         rep(c("L2FC","padj"),length(splitmarkers)), sep = "")
			joinmarkers <- joinmarkers[order(joinmarkers\$Clus0.L2FC, decreasing = TRUE),]
			joinmarkers <- data.frame(gene = rownames(joinmarkers), joinmarkers)
			write.table(joinmarkers, file = paste0(seurat_obj@project.name,"_markers_res",i,".tsv"), 
			            quote = FALSE, sep = "\\t", row.names = FALSE)
	    })
	}
	
    plot_list <- list()
    plot_heatmap <- list()
	for(i in res){
	try({
		markers <- markers_res_list[[paste0(i)]] 
		top_cluster_markers <- as_tibble(markers) %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 12)
		plot_heatmap[[paste0(i)]] <- DoHeatmap(seurat_obj, features = top_cluster_markers\$gene, group.by = paste0("RNA_snn_res.",i)) + NoLegend()
		plot_list[[paste0(i)]] <- RidgePlot(seurat_obj,features = "nCount_RNA", group.by=paste0("RNA_snn_res.",i)) + 
		   labs(title = paste0("UMI Density of Clusters (resolution ",i,")")) + ylab(label = "Clusters")
	})
	}
	for(i in res){
	    pdf(paste0(seurat_obj@project.name,"_Heatmap.pdf"), width = 8, height = 5)
	    print(plot_heatmap)
	    dev.off()
	    pdf(paste0(seurat_obj@project.name,"_UMIdensity.pdf"), width = 8, height = 5)
	    par(mfrow = c(2,2))
	    print(plot_list)
	    dev.off()
     }
}

# saveRDS
saveRDS(seurat_obj, paste(seurat_obj@project.name,"_seurat.rds", sep = ""))
"""
}


process Create_h5ad {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.h5ad$/) "Seurat_h5ad/$filename"}
input:
 file seurat_obj from g_41_rdsFile00_g_52

output:
 file "*.h5ad"  into g_52_h5_file00

"""
#!/usr/bin/env Rscript

# libraries
library(Seurat)
library(anndata)

# read data
seurat_obj <- readRDS("${seurat_obj}")

# save h5ad file
seu_name <- gsub(".rds","","${seurat_obj}")
write_h5ad(seurat_obj, paste0(seu_name,".h5ad"))

"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
