---
library("dplyr")
library("Seurat")
library("reticulate") 
library("Matrix")
---

###'
#' Before starting, localize dataset in an accessible location. The dataset should include three files:
#'   1) barcodes.tsv
#'   2) genes.tsv | features.tsv
#'   3) matrix.mtx
#' These three files will comprise the raw data structure to be manipulated downstream. If the files
#' have an accession ID, delete this identifier. If working with multiple datasets, the parent folder
#' holding the data should be labeled accordingly to avoid confusion.
#' 
#' R package descriptions provided by https://rdrr.io/github/satijalab/seurat/
#' 
#' NOTE: This code has been written prior to Wnt/Crect Single Cell RNA Seq Experiment. This pipeline
#'       is subject to change.
###' 



###' Load data and create Seurat object ----
#' Assign dataset directory to data_dir = "PATH/TO/DATASET"
#' Check that data_dir contains the correct files with list.files(data_dir). It should return three
#' files: 2 .tsv and 1 .mtx. These three files comprise the dataset.
#' Open genes.tsv | features.tsv to check which column contains informative gene names:
#'   EX: Column 1            Column 2
#'       ENSMUSG00000051951	 Xkr4
#'       ENSMUSG00000089699	 Gm1992
#'       ENSMUSG00000102343	 Gm37381
#'       [...]               [...]
##'

data_dir = "../jacobsdc/Desktop/GSE113576/"
list.files(data_dir) 

##'
#' Load dataset and apply to NAME.data with Read10X function. Specify the column in genes.tsv | 
#' features.tsv that has informative gene names with gene.column.
##' 


Zhang.data = Read10X(data.dir = data_dir, gene.column = 2)


###' Pruning the data to fit original manuscript input ----
#' The following prunes the raw Zhang.data matrix to fit the input according to the original manuscript.
#' According to the supplemental data: 
#'   "the data were normalized to the lowest saturated sample (78.1%)
#'    leading to 101,771 reads (mean), 2,461 genes (median), and 5,513 unique molecular identifier (UMI) 
#'    counts (median) per cell."
#' The authors do not provide steps to do this, therefore the following checks were done in order to 
#' figure out which parameters need to be adjusted to fit the input matrix.
#' Using Zhang.data, create a temporary matrix ('tmp') that only keeps the genes that have non-zero
#' counts, then keep the genes expressed in 3 or more cells.
#' tmp = apply(NAME.data,                   -> Name of original data structure
#'             MARGIN     = 1 | 2 | c(1,2), -> Manipulate based on rows | columns | both 
#'             function   = f(x)            -> Some function that prunes the data based on MARGIN value
#'             )
##'

tmp  = apply(Zhang.data, 
            1, 
            function(x) sum(x>0)) 

keep = tmp>=3
tmp  = Zhang.data[keep,]
one_cell_min = apply(tmp,
                     2,
                     function(x) sum(x>0))

summary(one_cell_min)


##'
#' NOTE: summary(one_cell_min) should return 2461 genes
###' 


###' Load data and create Seurat object (continued) ---- 
#' Using NAME.data, create a Seurat object from gene expression matrix. The object will be in the
#' format of genes(rows) x cells(columns).
#' CreateSeuratObject(counts       =  NAME.data,
#'                    project      = "Project Name",
#'                    assay        = "RNA",
#'                    min.cells    =  0,   -> Include genes detected in at least this many cells
#'                    min.features =  0,   -> Include cells where at least this many genes are found
#'                    names.field  =  1,   -> See NOTE
#'                    names.delim  = "?",  -> See NOTE
#'                    )
#' NOTE: Using these functions are dependent on the format of the NAME.data matrix and the goals of
#'       the project.
#'       For names.field: if columns (== cells) are named BARCODE_CELLTYPE and the CELLTYPE identifier
#'                        is of primary interest and not the individual cells, then set names.field = 2 
#'                        to set the intial identities to CELLTYPE.
#'       For names.delim: this tells names.field how to identify element 2 by defining how column names
#'                        are labeled. For example, if they are named by BARCODE_CELLTYPE, then set 
#'                        names.delim = "_".
##'
     
Zhang  = CreateSeuratObject(counts       =  Zhang.data,
                            project      = "Zhang Paper",
                            assay        = "RNA",
                            min.cells    =  3,
                            min.features =  203)

##' 
#'  The following code can be used to check how many genes and features per cell are represented in
#'  the Seurat object (Zhang) created above:
#'  median(Zhang@meta.data[,3])
#'  median(Zhang@meta.data[,2])
#'   
#'    
#'     
#'      
#'        
###' Check the structure of the data and filter ----
head(Zhang@meta.data, 5)
plot1 <- FeatureScatter(Zhang, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Zhang <- subset(Zhang,
                subset = nFeature_RNA > 200 & nFeature_RNA < 2500)






### Normalize the data ###
Zhang <- NormalizeData(Zhang, normalization.method = "LogNormalize", scale.factor = 10000)

### Identifying unique features ###
### NOTE: finds cells that exhibit high cell-to-cell variation ###
Zhang <- FindVariableFeatures(Zhang, selection.method = "vst", nfeatures = 3000)
top15 <- head(VariableFeatures(Zhang), 15)

plot2 <- VariableFeaturePlot(Zhang)
plot3 <- LabelPoints(plot = plot2, points = top15, repel = TRUE)
CombinePlots(plots = list(plot2, plot3))

### Scale data for PCA analysis ###
all.genes <- rownames(Zhang)
Zhang <- ScaleData(Zhang, 
                   features = all.genes)

### Linear dimension reduction ###
Zhang <- RunPCA(Zhang, 
                features = VariableFeatures(object = Zhang)
### Check the PCA in text form ###
print(Zhang[["pca"]], dims = 1:5, nfeatures = 5)
### Or graphically ###
DimPlot(Zhang, reduction = "pca")

### RunJackStraw to determine how many PCAs to include in the clustering analysis ###
### NOTE: Check PCA 1 through 20 ###
Zhang <- JackStraw(Zhang,
                    num.replicate = 100)
Zhang <- ScoreJackStraw(Zhang, dims = 1:20)
JackStrawPlot(Zhang, dims = 1:15)

ElbowPlot(Zhang) # Shows that after 15 PCAs most of the variance is captured 

### Clustering cells based on PCA ###
Zhang <- FindNeighbors(Zhang, dims = 1:15)
Zhang <- FindClusters(Zhang, resolution = 0.5)
Zhang <- RunUMAP(Zhang, dims = 1:15)

DimPlot(Zhang, reduction = "umap")

### Find the markers ###
cluster1.markers <- FindMarkers(Zhang, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

cluster5.markers <- FindMarkers(Zhang, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

Zhang.markers <- FindAllMarkers(Zhang, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster1.markers <- FindMarkers(Zhang, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(Zhang, features = c("Kiss1", "Kiss1r"))
VlnPlot(Zhang, features = c("Kiss1", "Kiss1r"), slot = "counts", log = TRUE)
FeaturePlot(Zhang, 
            features = c("Kiss1", "Kiss1r"))