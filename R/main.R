

# ### Tools used in DR.SC -------------------------------------------------

# One sample function
sce2seurat <- function(sce, verbose= TRUE){

  ## Transfer SingleCellExperiment object to a Seurat object for preparation for DR.SC model fitting.
  if(verbose){
    message("Transfer SingleCellExperiment object to a Seurat object")
    message("preparation for  model fitting")
  }

  require(SingleCellExperiment)
  count <- counts(sce)
  meta_data <- as.data.frame(colData(sce))
  require(Seurat)
  seu <- CreateSeuratObject(counts=count, meta.data = meta_data)

  return(seu)
}

spe2seurat <- function(spe, verbose= TRUE){

  ## Transfer SpatialExperiment object to a Seurat object for preparation for DR.SC model fitting.
  if(verbose){
    message("Transfer SpatialExperiment object to a Seurat object")
    message("preparation for  model fitting")
  }
  require(Seurat)
  suppressPackageStartupMessages(require(SpatialExperiment))
  if(is.null(colnames(spe))) stop("spe2seurat: Check argument spe! spe must have colnames!")
  if(is.null(row.names(spe))) stop("spe2seurat: Check argument spe! spe must have row.names!")

  col_meta_data <- as.data.frame(colData(spe))
  meta.data <- cbind.data.frame(col_meta_data,data.frame(row=spatialCoords(spe)[,1],
                                           col=spatialCoords(spe)[,2]))
  ret <- CreateSeuratObject(
    counts=assays(spe)$counts,
    meta.data=meta.data
  )

  return(ret)
}


## spe2seuratList(spe)
# Tools used in PRECAST ---------------------------------------------------


# Multiple sample

spe2seuratList <- function(spe, batch=NULL, verbose=TRUE){

  ## Transfer one SpatialExperiment object by its samples(batches)  to a Seurat list object for preparation for PRECAST model fitting.

  if(verbose){
    message("Transfer one SpatialExperiment object by its samples(batches)  to a Seurat list object for preparation for PRECAST model fitting.")
  }
   require(purrr)
  if(is.null(batch) && (!is.list(spe))) stop("If batch=NULL, spe must be a list with each component to be SpatialExperiment object!")

  if(!is.list(spe) && (!is.null(batch))){
    uniq_sample_id <- unique(colData(spe)[,batch])
    # Create a seurate object for each unique sample_id
    seuList <- map(uniq_sample_id,
        .f = function(smp_id, spe){
          # browser()
          ret_spe <- spe[, colData(spe)[,batch] == smp_id]
          ret_seurat <- spe2seurat(ret_spe)

          return(ret_seurat)
        },
        spe = spe)


  }

  if(is.list(spe)){
    seuList <- pbapply::pblapply(spe, spe2seurat, verbose=verbose)
  }

  return(seuList)
}


seu2seuList <- function(seu, batch, verbose=TRUE){


  ## Transfer one Seurat object by its samples(batches)  to a Seurat list object for preparation for PRECAST model fitting.
  if(verbose){
    message("Transfer one Seurat object by its samples(batches)  to a Seurat list object for preparation for PRECAST model fitting.")
  }
    require(purrr)
    meta_data <- seu@meta.data
    uniq_sample_id <- unique(meta_data[,batch])
    # Create a seurate object for each unique sample_id
    map(uniq_sample_id,
        .f = function(smp_id, seu){
          # browser()
          ret_seurat <- seu[, meta_data[,batch] == smp_id]
          return(ret_seurat)
        },
        seu = seu)

}


AddCoord2Reduc <- function(seu, spatial_coords= c("row", "col"), embed_name= "position"){

   require(PRECAST)
   embeds <- as.matrix(seu@meta.data[,spatial_coords])
   #embed_name <- "position"
   seu <- Add_embed(embed = embeds, seu=seu, embed_name = embed_name, assay= DefaultAssay(seu))
   return(seu)
}

# Generate simulated data using splatter package-------------------------------------------

generate_count <- function(seu, annotated_label, NumSpatialDomain=7, NumBatches = 3,
                           sim_seed=1,J = 2000, batch_facLoc=0.1, batch_facScale = 0.1){

  library(SingleCellExperiment)
  library(splatter)
  ## read position and annotation label from real data dlpfc

  #pos = seu@meta.data[,spatial_coords]
  #y <- as.numeric(seu@meta.data[,c("layer_guess_reordered")])
  y <- as.numeric(seu@meta.data[,annotated_label])
  assay <- DefaultAssay(seu)
  cnts = as.matrix(seu[[assay]]@counts)
  n_spots <- ncol(cnts)
  init_params <- splatEstimate(cnts)


  #batch_facLoc = 0.1 ## Batch effects in location
  #batch_facScale = 0.1
  #C = length(unique(y)) # the number of spatial Domains

  C <- NumSpatialDomain
  I = NULL
  N = n_spots*2
  L = NumBatches ## number of NumBatches

  debug = FALSE

  # 1.simulate count data
  noBatch <- ifelse(batch_facLoc == 0, TRUE, FALSE)
  group_prob <- as.vector(table(y)/length(y))


  params <- setParams(
    init_params,
    batchCells = rep(N, L), # 3N here represents a large number such that
    # we have sufficient cells of each type to be
    # allocated to the spatial transcriptomics data
    batch.rmEffect = noBatch,
    batch.facLoc = batch_facLoc,
    batch.facScale = batch_facScale,
    nGenes = J,
    group.prob = group_prob,
    seed = sim_seed)

  sim_groups <- splatSimulate(
    params = params,
    method = "groups",
    verbose = FALSE)
  return(sim_groups)

}

simu_seuList <- function(seu, annotated_label, spatial_coords= c("row", "col"),
                             NumBatches = 3, NumSpotsRatio=0.9, ngenes = 2000, sim_seed=1,
                             batch_facLoc=0.1, batch_facScale = 0.1){


  # annotated_label= "annotation"; spatial_coords= c("row", "col");
  # NumBatches = 3; NumSpotsRatio=0.2; ngenes = 20; sim_seed=1;
  # batch_facLoc=0.1; batch_facScale = 0.1

  if(NumSpotsRatio<0.1) stop("simu_seuList:: Check argument: NumSpotsRatio! this argument must be greater than 0.1!")
  require(Seurat)
  ### Spatial domain annotations
  y <- as.numeric(seu@meta.data[,annotated_label])
  NumSpatialDomain <- length(unique(y))
  sim_groups <- generate_count(seu,annotated_label, NumSpatialDomain ,sim_seed = sim_seed, J= ngenes,
                               batch_facLoc=batch_facLoc, batch_facScale = batch_facScale)
  ## reorder the sample to  be the same as y
  meta_data_all <- colData(sim_groups)


  Groupy = paste0("Group", y)
  nbatch <- NumBatches

  pos <- as.matrix(seu@meta.data[,spatial_coords])
  gen_pos_id <- function(iseed, pos){
    set.seed(iseed)
    quantx <- runif(1, NumSpotsRatio-0.1 ,NumSpotsRatio+0.1)
    quanty <- runif(1, NumSpotsRatio-0.1 ,NumSpotsRatio+0.1)
    pos_id <- which(pos[,1]< quantile(pos[,1], quantx) & pos[,2]< quantile(pos[,2], quanty)  )
    return(pos_id)
  }

  idxList_pos <- list()
  for(i_batch in 1: nbatch){

    # list(1:length(y), which(pos[,1]< quantile(pos[,1], 0.9)), which(pos[,2]< quantile(pos[,2], 0.9)))
    idxList_pos[[i_batch]] <- gen_pos_id(i_batch, pos)
  }
  posList <- lapply(idxList_pos, function(idx) pos[idx, ])


  yList <- lapply(idxList_pos, function(idx) y[idx])

  idxList_sim <- list()
  for(r in 1:nbatch){
    message("r = ", r)
    y_tmp <- yList[[r]]
    num_each_celltype = table(y_tmp)
    i_vec <- as.numeric(names(num_each_celltype))
    idx1 = rep(0, length(y_tmp))
    for (i in i_vec){
      idx1[y_tmp==i] = which(meta_data_all$Group == paste0("Group",i) & meta_data_all$Batch==paste0("Batch", r))[1:num_each_celltype[i]]
    }
    idxList_sim[[r]] <- idx1 ## align with posList

  }

  sceList_sim <- lapply(idxList_sim, function(idx) sim_groups[,idx])
  ## Add spatial coordinates
  sceList_sim <- lapply(1: nbatch, function(r){
    sce <- sceList_sim[[r]]
    colData(sce)$row <- posList[[r]][,1]
    colData(sce)$col <- posList[[r]][,2]
    return(sce)
  })

  seuList_use <- pbapply::pblapply(sceList_sim, sce2seurat)
  return(seuList_use)
}


