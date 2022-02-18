#function to get "not in" from %in% in dplyr
'%!in%' <- function(x,y)!('%in%'(x,y))

#read and analyze
read_loom_and_analyze<-function(loom_file="./loom/s_fca_biohub_body_wall_10x.loom", plot=T) {

  require(assertthat)

  assert_that(
    see_if(is.readable(loom_file))
  )


  print(paste0("running sample ",loom_file))

  #read loom
  loom_object <- connect(filename = loom_file, mode = "r",skip.validate = T)

  #extract genes
  gene.names <- loom_object[["row_attrs/Gene"]][]

  #extract reads
  reads_matrix<-as.matrix(loom_object[["matrix"]][,])
  colnames(reads_matrix)<-gene.names
  reads_matrix<-t(reads_matrix[,colSums(reads_matrix[,])>30])

  # Pull three bits of metadata from the column attributes
  attrs <- c("n_counts", "n_genes", "age","ClusterID","CellID","tissue","sex","sample_id","percent_mito","fly_genetics","dissection_lab","batch_id","annotation")
  attr.df <- loom_object$get.attribute.df(MARGIN = 2, attributes =  attrs)

  #close loom
  loom_object$close_all()
  rm(loom_object)

  #read h5
  Convert(gsub(loom_file, pattern = ".loom",replacement = ".h5ad"), dest = "h5seurat", overwrite = TRUE)
  loom_object_h5<- LoadH5Seurat(gsub(loom_file, pattern = ".loom",replacement = ".h5seurat"))

  #attach names to reads matrix
  colnames(reads_matrix)<-colnames(loom_object_h5@assays$RNA@counts)

  #attach meta data
  loom_object_h5@meta.data=cbind(loom_object_h5@meta.data,attr.df)

  #create name for file and plots
  nn=gsub(pattern = "s_fca_biohub_",replacement = "",basename(loom_file))
  nn=gsub(pattern = "_10x.loom",replacement = "",nn)

  #create plots with basic analysis
  if (plot){
    pdf(paste0("General_",nn,"_DimPlot_fromloom.pdf"))
    p=DimPlot(loom_object_h5,repel = T,label = T,reduction = "umap",group.by = "annotation")+ NoLegend() + labs(title =  paste0("Clusters from loom ",nn))
    print(p)
    p=DimPlot(loom_object_h5,reduction = "umap",group.by = "tissue")+ labs(title = paste0("tissue ",nn))
    print(p)
    p=DimPlot(loom_object_h5,reduction = "umap",group.by = "sex")+ labs(title =  paste0("sex ",nn))
    print(p)
    p=DimPlot(loom_object_h5,reduction = "umap",group.by = "age")+ labs(title =  paste0("age ",nn))
    print(p)
    dev.off()
  }


  #To be able to see any gene you are interested in you have to re-do the basic analysis from screatch. This is because the provided file only has most variable genes

  print("entering new seurat analysis")
  #create new seurat
  new_seurat_obj <- CreateSeuratObject(counts = reads_matrix)
  new_seurat_obj <- NormalizeData(object = new_seurat_obj)
  new_seurat_obj <- FindVariableFeatures(object = new_seurat_obj)
  new_seurat_obj <- ScaleData(object = new_seurat_obj)
  new_seurat_obj <- RunPCA(object = new_seurat_obj)
  new_seurat_obj <- FindNeighbors(object = new_seurat_obj,dims = 1:30)
  new_seurat_obj <- FindClusters(object = new_seurat_obj)
  new_seurat_obj <- RunTSNE(object = new_seurat_obj,dims = 1:30)
  new_seurat_obj <- RunUMAP(object = new_seurat_obj,dims = 1:30)

  #save old idents and transfere them
  new_seurat_obj[["old.ident"]] <- Idents(object = new_seurat_obj)
  new_seurat_obj@meta.data=cbind(new_seurat_obj@meta.data,attr.df)

  #Now we can see our new analysis
  if (plot){
    pdf(paste0(nn,"_DimPlot_newanalysis.pdf"))

    p=DimPlot(new_seurat_obj,repel = T,label = T,reduction = "tsne")+ labs(title = nn)
    print(p)
    p=DimPlot(new_seurat_obj,repel = T,label = T,reduction = "tsne", group.by = "annotation")+ NoLegend() + labs(title = nn)
    print(p)
    p=DimPlot(new_seurat_obj,repel = T,label = T,reduction = "umap")+ labs(title = nn)
    print(p)
    p=DimPlot(new_seurat_obj,repel = T,label = T,reduction = "umap", group.by = "annotation")+ NoLegend() + labs(title = nn)

    print(p)

    dev.off()
  }

  print("calculting markers")
  #get markers for new ID
  new_seurat_obj@active.ident<-as.factor(new_seurat_obj$old.ident)
  new_seurat_obj<-RenameCells(new_seurat_obj,new.names =new_seurat_obj$CellID)
  all.markers.new<-FindAllMarkers(new_seurat_obj,only.pos = T)
  #get markers for paper ID
  new_seurat_obj@active.ident<-as.factor(gsub(pattern = " ",replacement = "",new_seurat_obj$annotation))
  all.markers.paper<-FindAllMarkers(new_seurat_obj,only.pos = T)
  #extract top10
  top10markers <- all.markers.paper %>% group_by(cluster) %>%top_n(n = 10, wt = avg_log2FC)

  if(plot){
    pdf(paste0(nn,"_markers_plots",".pdf"),width = 15,height = 30)

    cluster.averages <- AverageExpression(new_seurat_obj, return.seurat = TRUE, add.ident = "sex")
    p=DoHeatmap(cluster.averages, features = top10markers$gene) + NoLegend()
    print(p)

    dev.off()
  }

  #now we rename the objects and save the data
  print("saving data")
  assign(x = paste0(nn,"_h5"),value = loom_object_h5)
  assign(x = paste0(nn,"_new"),value = new_seurat_obj)
  assign(x = paste0(nn,"_new_markers"),value = all.markers.new)
  assign(x = paste0(nn,"_paper_markers"),value = all.markers.paper)
  assign(x = paste0(nn,"_matadata"),value = attr.df)

  #this will save the R objects
  save(list = ls(pattern = paste0("^",nn)),file = paste0("data_",nn,".RData"))

  #rm(list =c("loom_object_h5","new_seurat_obj","cluster.averages","all.markers.new","all.markers.paper","newdata","x","reads_matrix","loom_object","attr.df") )

}

loom_too_Seurat<-function(loom_file="./loom/s_fca_biohub_body_wall_10x.loom", plot=T) {

  require(assertthat)

  assert_that(
    see_if(is.readable(loom_file))
  )


  print(paste0("running sample ",loom_file))

  #read loom
  loom_object <- connect(filename = loom_file, mode = "r",skip.validate = T)

  #extract genes
  gene.names <- loom_object[["row_attrs/Gene"]][]

  #extract reads
  reads_matrix<-as.matrix(loom_object[["matrix"]][,])
  colnames(reads_matrix)<-gene.names
  reads_matrix<-t(reads_matrix[,colSums(reads_matrix[,])>30])

  # Pull three bits of metadata from the column attributes
  attrs <- c("n_counts", "n_genes", "age","ClusterID","CellID","tissue","sex","sample_id","percent_mito","fly_genetics","dissection_lab","batch_id","annotation")
  attr.df <- loom_object$get.attribute.df(MARGIN = 2, attributes =  attrs)

  #close loom
  loom_object$close_all()
  rm(loom_object)

  #read h5
  Convert(gsub(loom_file, pattern = ".loom",replacement = ".h5ad"), dest = "h5seurat", overwrite = TRUE)
  loom_object_h5<- LoadH5Seurat(gsub(loom_file, pattern = ".loom",replacement = ".h5seurat"))

  #attach names to reads matrix
  colnames(reads_matrix)<-colnames(loom_object_h5@assays$RNA@counts)

  #attach meta data
  loom_object_h5@meta.data=cbind(loom_object_h5@meta.data,attr.df)

  #create name for file and plots
  nn=gsub(pattern = "s_fca_biohub_",replacement = "",basename(loom_file))
  nn=gsub(pattern = "_10x.loom",replacement = "",nn)

  return(loom_object_h5)

  #create plots with basic analysis
  if (plot){
    pdf(paste0("General_",nn,"_DimPlot_fromloom.pdf"))
    p=DimPlot(loom_object_h5,repel = T,label = T,reduction = "umap",group.by = "annotation")+ NoLegend() + labs(title =  paste0("Clusters from loom ",nn))
    print(p)
    p=DimPlot(loom_object_h5,reduction = "umap",group.by = "tissue")+ labs(title = paste0("tissue ",nn))
    print(p)
    p=DimPlot(loom_object_h5,reduction = "umap",group.by = "sex")+ labs(title =  paste0("sex ",nn))
    print(p)
    p=DimPlot(loom_object_h5,reduction = "umap",group.by = "age")+ labs(title =  paste0("age ",nn))
    print(p)
    dev.off()
  }
}

loom_too_Seurat_brandNew<-function(loom_file="./loom/s_fca_biohub_body_wall_10x.loom", plot=T,preprocesing=T) {

  require(assertthat)

  assert_that(
    see_if(is.readable(loom_file))
  )


  print(paste0("running sample ",loom_file))

  #read loom
  loom_object <- connect(filename = loom_file, mode = "r",skip.validate = T)

  #extract genes
  gene.names <- loom_object[["row_attrs/Gene"]][]

  #extract reads
  reads_matrix<-as.matrix(loom_object[["matrix"]][,])
  colnames(reads_matrix)<-gene.names
  reads_matrix<-t(reads_matrix[,colSums(reads_matrix[,])>30])

  # Pull three bits of metadata from the column attributes
  attrs <- c("n_counts", "n_genes", "age","ClusterID","CellID","tissue","sex","sample_id","percent_mito","fly_genetics","dissection_lab","batch_id","annotation")
  attr.df <- loom_object$get.attribute.df(MARGIN = 2, attributes =  attrs)

  #close loom
  loom_object$close_all()
  rm(loom_object)

  #read h5
  Convert(gsub(loom_file, pattern = ".loom",replacement = ".h5ad"), dest = "h5seurat", overwrite = TRUE)
  loom_object_h5<- LoadH5Seurat(gsub(loom_file, pattern = ".loom",replacement = ".h5seurat"))

  #attach names to reads matrix
  colnames(reads_matrix)<-colnames(loom_object_h5@assays$RNA@counts)

  #attach meta data
  loom_object_h5@meta.data=cbind(loom_object_h5@meta.data,attr.df)

  #create name for file and plots
  nn=gsub(pattern = "s_fca_biohub_",replacement = "",basename(loom_file))
  nn=gsub(pattern = "_10x.loom",replacement = "",nn)

  return(loom_object_h5)

  #create plots with basic analysis
  if (plot){
    pdf(paste0("General_",nn,"_DimPlot_fromloom.pdf"))
    p=DimPlot(loom_object_h5,repel = T,label = T,reduction = "umap",group.by = "annotation")+ NoLegend() + labs(title =  paste0("Clusters from loom ",nn))
    print(p)
    p=DimPlot(loom_object_h5,reduction = "umap",group.by = "tissue")+ labs(title = paste0("tissue ",nn))
    print(p)
    p=DimPlot(loom_object_h5,reduction = "umap",group.by = "sex")+ labs(title =  paste0("sex ",nn))
    print(p)
    p=DimPlot(loom_object_h5,reduction = "umap",group.by = "age")+ labs(title =  paste0("age ",nn))
    print(p)
    dev.off()
  }

  #To be able to see any gene you are interested in you have to re-do the basic analysis from screatch. This is because the provided file only has most variable genes

  print("entering new seurat analysis")
  #create new seurat
  new_seurat_obj <- CreateSeuratObject(counts = reads_matrix)

  if(preprocesing){
    new_seurat_obj <- NormalizeData(object = new_seurat_obj)
    new_seurat_obj <- FindVariableFeatures(object = new_seurat_obj)
    new_seurat_obj <- ScaleData(object = new_seurat_obj)
    new_seurat_obj <- RunPCA(object = new_seurat_obj)
    new_seurat_obj <- FindNeighbors(object = new_seurat_obj,dims = 1:30)
    new_seurat_obj <- FindClusters(object = new_seurat_obj)
    new_seurat_obj <- RunTSNE(object = new_seurat_obj,dims = 1:30)
    new_seurat_obj <- RunUMAP(object = new_seurat_obj,dims = 1:30)

    #save old idents and transfere them
    new_seurat_obj[["old.ident"]] <- Idents(object = new_seurat_obj)
    new_seurat_obj@meta.data=cbind(new_seurat_obj@meta.data,attr.df)

    #Now we can see our new analysis
    if (plot){
      pdf(paste0(nn,"_DimPlot_newanalysis.pdf"))

      p=DimPlot(new_seurat_obj,repel = T,label = T,reduction = "tsne")+ labs(title = nn)
      print(p)
      p=DimPlot(new_seurat_obj,repel = T,label = T,reduction = "tsne", group.by = "annotation")+ NoLegend() + labs(title = nn)
      print(p)
      p=DimPlot(new_seurat_obj,repel = T,label = T,reduction = "umap")+ labs(title = nn)
      print(p)
      p=DimPlot(new_seurat_obj,repel = T,label = T,reduction = "umap", group.by = "annotation")+ NoLegend() + labs(title = nn)

      print(p)

      dev.off()
    }

  }

  return(new_seurat_obj)


}


