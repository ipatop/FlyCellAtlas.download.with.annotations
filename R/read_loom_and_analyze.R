#function to get "not in" from %in% in dplyr
'%!in%' <- function(x,y)!('%in%'(x,y))

#' Read and process loom and H5AD files from FlyCellAtlas
#'
#'This function requires loom and h5ad files from FlyCellAtlas and process them into Seurat objects. It then also re-analyze them to get the full gene-expression matrix and not only the extract marker genes and also get further information.
#'
#' @param loom_file character with Name and directory of loom file -h5da has to be in same directory-
#' @param plot Boolean indicating if plots should be generated
#' @param dir Directory to save RData and plots in
#'
#' @return It saves old and new Seurat onjects along with marker genes and metadata.
#'   *"_h5" is seurat object from the paper, this only has most variable genes
#'
#'   *"_new" is the new seurat object
#'
#'   *"_new_markers" is all.marker genes from the new seurat object with all genes
#'
#'   *"_paper_markers" is all.marker genes from the paper
#'
#'   *,"_matadata" is the metadata
#' @export
#'
#' @examples
#' There are examples in the folder names "./loom". To run it:
#' read_loom_and_analyze(loom_file = "~/Documents/scripts/FlyCellAtlas_download_with_annotations/loom/s_fca_biohub_body_wall_10x.loom")
#'
read_loom_and_analyze<-function(loom_file="./loom/s_fca_biohub_body_wall_10x.loom", plot=T, dir="./") {

  assertthat::assert_that(
    assertthat::see_if(assertthat::is.readable(loom_file))
  )

  print(paste0("running sample ",loom_file))

  #read loom
  loom_object <- loomR::connect(filename = loom_file, mode = "r+",skip.validate = T)

  #extract genes
  gene.names <- loom_object[["row_attrs/Gene"]][]

  #extract reads
  reads_matrix<-as.matrix(loom_object[["matrix"]][,])
  colnames(reads_matrix)<-gene.names
  reads_matrix<-t(reads_matrix[,colSums(reads_matrix[,])>30])

  # Pull three bits of metadata from the column attributes
  attrs <- c("n_counts", "n_genes", "age","ClusterID","CellID","tissue","sex","sample_id","percent_mito","fly_genetics","dissection_lab","batch_id","annotation")

  #add cell_name atribute for it to run
  loom_object$add.attribute(MARGIN = 2,list(cell_names = loom_object[["col_attrs/CellID"]][]))

  #extract metadata
  attr.df <- loom_object$get.attribute.df(MARGIN = 2, attribute.names  =  attrs)

  #close loom
  loom_object$close_all()

  #read h5
  SeuratDisk::Convert(gsub(loom_file, pattern = "x.loom",replacement = "x.h5ad"), dest = "h5seurat", overwrite = TRUE)
  loom_object_h5<- SeuratDisk::LoadH5Seurat(gsub(loom_file, pattern = "x.loom",replacement = "x.h5seurat"))

  #attach names to reads matrix
  colnames(reads_matrix)<-colnames(loom_object_h5@assays$RNA@counts)

  #attach meta data
  loom_object_h5@meta.data=cbind(loom_object_h5@meta.data,attr.df)

  #create name for file and plots
  nn=gsub(pattern = "s_fca_biohub_",replacement = "",basename(loom_file))
  nn=gsub(pattern = "_10x.loom",replacement = "",nn)

  #create plots with basic analysis
  if (plot){
   pdf(paste0(dir,"/General_",nn,"_DimPlot_fromloom.pdf"))
    p=Seurat::DimPlot(loom_object_h5,repel = T,label = T,reduction = "umap",group.by = "annotation")+ NoLegend() + labs(title =  paste0("Clusters from loom ",nn))
    print(p)
    p=Seurat::DimPlot(loom_object_h5,reduction = "umap",group.by = "tissue")+ labs(title = paste0("tissue ",nn))
    print(p)
    p=Seurat::DimPlot(loom_object_h5,reduction = "umap",group.by = "sex")+ labs(title =  paste0("sex ",nn))
    print(p)
    p=Seurat::DimPlot(loom_object_h5,reduction = "umap",group.by = "age")+ labs(title =  paste0("age ",nn))
    print(p)
    dev.off()
  }

  #To be able to see any gene you are interested in you have to re-do the basic analysis from scratch. This is because the provided file only has most variable genes

  print("entering new seurat analysis")
  #create new Seurat
  new_seurat_obj <- SeuratObject::CreateSeuratObject(counts = reads_matrix)
  new_seurat_obj <- Seurat::NormalizeData(object = new_seurat_obj)
  new_seurat_obj <- Seurat::FindVariableFeatures(object = new_seurat_obj)
  new_seurat_obj <- Seurat::ScaleData(object = new_seurat_obj)
  new_seurat_obj <- Seurat::RunPCA(object = new_seurat_obj)
  new_seurat_obj <- Seurat::FindNeighbors(object = new_seurat_obj,dims = 1:30)
  new_seurat_obj <- Seurat::FindClusters(object = new_seurat_obj)
  new_seurat_obj <- Seurat::RunTSNE(object = new_seurat_obj,dims = 1:30)
  new_seurat_obj <- Seurat::RunUMAP(object = new_seurat_obj,dims = 1:30)

  #save old idents and transfere them
  new_seurat_obj[["old.ident"]] <- Seurat::Idents(object = new_seurat_obj)
  new_seurat_obj@meta.data=cbind(new_seurat_obj@meta.data,attr.df)

  #Now we can see our new analysis
  if (plot){

    grDevices::pdf(paste0(dir,"/DimPlot_",nn,"_newanalysis.pdf"))

    p=Seurat::DimPlot(new_seurat_obj,repel = T,label = T,reduction = "tsne")+ labs(title = nn)
    print(p)
    p=Seurat::DimPlot(new_seurat_obj,repel = T,label = T,reduction = "tsne", group.by = "annotation")+ NoLegend() + labs(title = nn)
    print(p)
    p=Seurat::DimPlot(new_seurat_obj,repel = T,label = T,reduction = "umap")+ labs(title = nn)
    print(p)
    p=Seurat::DimPlot(new_seurat_obj,repel = T,label = T,reduction = "umap", group.by = "annotation")+ NoLegend() + labs(title = nn)

    print(p)

    dev.off()
  }

  print("calculting markers")
  #get markers for new ID
  new_seurat_obj@active.ident<-as.factor(new_seurat_obj$old.ident)
  new_seurat_obj<-Seurat::RenameCells(new_seurat_obj,new.names =new_seurat_obj$CellID)
  all.markers.new<-Seurat::FindAllMarkers(new_seurat_obj,only.pos = T)
  #get markers for paper ID
  new_seurat_obj@active.ident<-as.factor(gsub(pattern = " ",replacement = "",new_seurat_obj$annotation))
  all.markers.paper<-Seurat::FindAllMarkers(new_seurat_obj,only.pos = T)
  #extract top10
  top10markers <- all.markers.paper %>% group_by(cluster) %>%top_n(n = 10, wt = avg_log2FC)

  if(plot){
    pdf(paste0(dir,"/markers_plots_",nn,".pdf"),width = 15,height = 30)

    cluster.averages <- Seurat::AverageExpression(new_seurat_obj, return.seurat = TRUE, add.ident = "sex")
    p=Seurat::DoHeatmap(cluster.averages, features = top10markers$gene) + NoLegend()
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
  save(list = ls(pattern = paste0("^",nn)),file = paste0(dir,"/data_",nn,".RData"))

  #rm(list =c("loom_object_h5","new_seurat_obj","cluster.averages","all.markers.new","all.markers.paper","newdata","x","reads_matrix","loom_object","attr.df") )

}



#' Read and process loom and H5AD files from FlyCellAtlas
#'
#'This function requires loom and h5ad files from FlyCellAtlas and process them into Seurat objects.
#'
#' @param loom_file character with Name and directory of loom file -h5da has to be in same directory-
#' @param plot Boolean indicating if plots should be generated
#' @param dir Directory to save RData and plots in
#'
#' @return Seurat object from paper without any further processing
#' @export
#'
#' @examples
#' There are examples in the folder names "./loom". To run it:
#' loom_to_seurat(loom_file = "~/Documents/scripts/FlyCellAtlas_download_with_annotations/loom/s_fca_biohub_body_wall_10x.loom")
#'
loom_too_Seurat<-function(loom_file="./loom/s_fca_biohub_body_wall_10x.loom", plot=T) {


  assertthat::assert_that(
    assertthat::see_if(assertthat::is.readable(loom_file))
  )

  print(paste0("running sample ",loom_file))

  #read loom
  loom_object <- loomR::connect(filename = loom_file, mode = "r+",skip.validate = T)

  #extract genes
  gene.names <- loom_object[["row_attrs/Gene"]][]

  #extract reads
  reads_matrix<-as.matrix(loom_object[["matrix"]][,])
  colnames(reads_matrix)<-gene.names
  reads_matrix<-t(reads_matrix[,colSums(reads_matrix[,])>30])

  # Pull three bits of metadata from the column attributes
  attrs <- c("n_counts", "n_genes", "age","ClusterID","CellID","tissue","sex","sample_id","percent_mito","fly_genetics","dissection_lab","batch_id","annotation")

  #add cell_name atribute for it to run
  loom_object$add.attribute(MARGIN = 2,list(cell_names = loom_object[["col_attrs/CellID"]][]))

  #extract metadata
  attr.df <- loom_object$get.attribute.df(MARGIN = 2, attribute.names  =  attrs)

  #close loom
  loom_object$close_all()

  #read h5
  SeuratDisk::Convert(gsub(loom_file, pattern = "x.loom",replacement = "x.h5ad"), dest = "h5seurat", overwrite = TRUE)
  loom_object_h5<- SeuratDisk::LoadH5Seurat(gsub(loom_file, pattern = "x.loom",replacement = "x.h5seurat"))

  #attach names to reads matrix
  colnames(reads_matrix)<-colnames(loom_object_h5@assays$RNA@counts)

  #attach meta data
  loom_object_h5@meta.data=cbind(loom_object_h5@meta.data,attr.df)



  #create plots with basic analysis
  if (plot){
    pdf(paste0(dir,"/General_",nn,"_DimPlot_fromloom.pdf"))
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
  return(loom_object_h5)
}


#' Read and process loom and H5AD files from FlyCellAtlas
#'
#'This function requires loom and h5ad files from FlyCellAtlas and process them into Seurat objects. It creates a new object form scratch and tranfere celltype assignment
#'
#' @param loom_file character with Name and directory of loom file -h5da has to be in same directory-
#' @param plot Boolean indicating if plots should be generated
#' @param dir Directory to save plots in
#' @param preprocesing Boolean indicating if preprocessing should be done (ie. findVariableFeatues, Scale Data, etc) or not
#'
#' @return It saves old and new Seurat onjects along with marker genes and metadata.
#'   *"_h5" is seurat object from the paper, this only has most variable genes
#'
#'   *"_new" is the new seurat object
#'
#'   *"_new_markers" is all.marker genes from the new seurat object with all genes
#'
#'   *"_paper_markers" is all.marker genes from the paper
#'
#'   *,"_matadata" is the metadata
#' @export
#'
#' @examples
#' There are examples in the folder names "./loom". To run it:
#' loom_too_Seurat_brandNew(loom_file = "~/Documents/scripts/FlyCellAtlas_download_with_annotations/loom/s_fca_biohub_body_wall_10x.loom")
#'
loom_too_Seurat_brandNew<-function(loom_file="./loom/s_fca_biohub_body_wall_10x.loom", plot=T,preprocesing=T,dir="./") {
  assertthat::assert_that(
    assertthat::see_if(assertthat::is.readable(loom_file))
  )

  print(paste0("running sample ",loom_file))

  #read loom
  loom_object <- loomR::connect(filename = loom_file, mode = "r+",skip.validate = T)

  #extract genes
  gene.names <- loom_object[["row_attrs/Gene"]][]

  #extract reads
  reads_matrix<-as.matrix(loom_object[["matrix"]][,])
  colnames(reads_matrix)<-gene.names
  reads_matrix<-t(reads_matrix[,colSums(reads_matrix[,])>30])

  # Pull three bits of metadata from the column attributes
  attrs <- c("n_counts", "n_genes", "age","ClusterID","CellID","tissue","sex","sample_id","percent_mito","fly_genetics","dissection_lab","batch_id","annotation")

  #add cell_name atribute for it to run
  loom_object$add.attribute(MARGIN = 2,list(cell_names = loom_object[["col_attrs/CellID"]][]))

  #extract metadata
  attr.df <- loom_object$get.attribute.df(MARGIN = 2, attribute.names  =  attrs)

  #close loom
  loom_object$close_all()

  #read h5
  SeuratDisk::Convert(gsub(loom_file, pattern = "x.loom",replacement = "x.h5ad"), dest = "h5seurat", overwrite = TRUE)
  loom_object_h5<- SeuratDisk::LoadH5Seurat(gsub(loom_file, pattern = "x.loom",replacement = "x.h5seurat"))

  #attach names to reads matrix
  colnames(reads_matrix)<-colnames(loom_object_h5@assays$RNA@counts)

  #attach meta data
  loom_object_h5@meta.data=cbind(loom_object_h5@meta.data,attr.df)

  #create name for file and plots
  nn=gsub(pattern = "s_fca_biohub_",replacement = "",basename(loom_file))
  nn=gsub(pattern = "_10x.loom",replacement = "",nn)

  #create plots with basic analysis
  if (plot){
   pdf(paste0(dir,"/General_",nn,"_DimPlot_fromloom.pdf"))
    p=Seurat::DimPlot(loom_object_h5,repel = T,label = T,reduction = "umap",group.by = "annotation")+ NoLegend() + labs(title =  paste0("Clusters from loom ",nn))
    print(p)
    p=Seurat::DimPlot(loom_object_h5,reduction = "umap",group.by = "tissue")+ labs(title = paste0("tissue ",nn))
    print(p)
    p=Seurat::DimPlot(loom_object_h5,reduction = "umap",group.by = "sex")+ labs(title =  paste0("sex ",nn))
    print(p)
    p=Seurat::DimPlot(loom_object_h5,reduction = "umap",group.by = "age")+ labs(title =  paste0("age ",nn))
    print(p)
    dev.off()
  }


  #To be able to see any gene you are interested in you have to re-do the basic analysis from screatch. This is because the provided file only has most variable genes

  print("entering new seurat analysis")
  #create new seurat
  new_seurat_obj <- SeuratObject::CreateSeuratObject(counts = reads_matrix)

  if(preprocesing){
    new_seurat_obj <- Seurat::NormalizeData(object = new_seurat_obj)
    new_seurat_obj <- Seurat::FindVariableFeatures(object = new_seurat_obj)
    new_seurat_obj <- Seurat::ScaleData(object = new_seurat_obj)
    new_seurat_obj <- Seurat::RunPCA(object = new_seurat_obj)
    new_seurat_obj <- Seurat::FindNeighbors(object = new_seurat_obj,dims = 1:30)
    new_seurat_obj <- Seurat::FindClusters(object = new_seurat_obj)
    new_seurat_obj <- Seurat::RunTSNE(object = new_seurat_obj,dims = 1:30)
    new_seurat_obj <- Seurat::RunUMAP(object = new_seurat_obj,dims = 1:30)

    #save old idents and transfere them
    new_seurat_obj[["old.ident"]] <- Seurat::Idents(object = new_seurat_obj)
    new_seurat_obj@meta.data=cbind(new_seurat_obj@meta.data,attr.df)

    #Now we can see our new analysis
    if (plot){
      pdf(paste0(nn,"_DimPlot_newanalysis.pdf"))

      p=Seurat::DimPlot(new_seurat_obj,repel = T,label = T,reduction = "tsne")+ labs(title = nn)
      print(p)
      p=Seurat::DimPlot(new_seurat_obj,repel = T,label = T,reduction = "tsne", group.by = "annotation")+ NoLegend() + labs(title = nn)
      print(p)
      p=Seurat::DimPlot(new_seurat_obj,repel = T,label = T,reduction = "umap")+ labs(title = nn)
      print(p)
      p=Seurat::DimPlot(new_seurat_obj,repel = T,label = T,reduction = "umap", group.by = "annotation")+ NoLegend() + labs(title = nn)

      print(p)

      dev.off()
    }

}

    new_seurat_obj@meta.data=cbind(new_seurat_obj@meta.data,attr.df)

    #Now we can see our new analysis
    if (plot){
      grDevices::pdf(paste0(dir,"/DimPlot_",nn,"_newanalysis.pdf"))

      p=Seurat::DimPlot(new_seurat_obj,repel = T,label = T,reduction = "tsne")+ ggplot2::labs(title = nn)
      print(p)
      p=Seurat::DimPlot(new_seurat_obj,repel = T,label = T,reduction = "tsne", group.by = "annotation")+ Seurat::NoLegend() + ggplot2::labs(title = nn)
      print(p)
      p=Seurat::DimPlot(new_seurat_obj,repel = T,label = T,reduction = "umap")+ ggplot2::labs(title = nn)
      print(p)
      p=Seurat::DimPlot(new_seurat_obj,repel = T,label = T,reduction = "umap", group.by = "annotation")+ Seurat::NoLegend() + ggplot2::labs(title = nn)

      print(p)

      grDevices::dev.off()
    }
  return(new_seurat_obj)

}

