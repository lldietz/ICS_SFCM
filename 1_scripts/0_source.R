# LOAD LIBRARIES
library(data.table)
library(formulaic)
library(flowCore)
library(flowViz)
library(PeacoQC)
library(CATALYST)
library(SingleCellExperiment)
library(uwot)
library(knitr)
library(stringr)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(glue)
library(cowplot)
library(scater)
library(diffcyt)
library(patchwork)
library(ComplexHeatmap)
library(flowWorkspace)
library(glue)
library(scales)
library(ggcyto)
library(CytoExploreR)
library(cyCombine)
library(outliers)
library(clustree)
library(ggrepel)
library(cytoqc)
library(gridExtra)
library(grid)
library(tibble)
library(tidyr)
library(ggpointdensity)
library(ggeasy)

# SET AND WRITE DIRECTORIES
dir_project = dirname(getwd())
scriptname = str_replace(basename(rstudioapi::getSourceEditorContext()$path),".Rmd|.R","")
dir_pipeline = glue("{dir_project}/2_pipeline/{scriptname}")
for (folder in c("tmp","store","out")) dir.create(glue("{dir_pipeline}/{folder}"), showWarnings = FALSE, recursive = TRUE)
dir_tmp = glue("{dir_pipeline}/tmp")
dir_store = glue("{dir_pipeline}/store")
dir_out = glue("{dir_pipeline}/out")
rm(list = c("folder","scriptname"))
folder_transform = grep("transform",list.dirs(glue("{dir_project}/2_pipeline"), full.names = FALSE, recursive = FALSE), value=TRUE,ignore.case = T)
folder_pregate =  grep("pregate",list.dirs(glue("{dir_project}/2_pipeline"), full.names = FALSE, recursive = FALSE), value=TRUE,ignore.case = T)
folder_QC = grep("QC|quality",list.dirs(glue("{dir_project}/2_pipeline"), full.names = FALSE, recursive = FALSE), value=TRUE,ignore.case = T)
folder_norm = grep("norm",list.dirs(glue("{dir_project}/2_pipeline"), full.names = FALSE, recursive = FALSE), value=TRUE,ignore.case = T)
folder_celltype =  grep("celltype",list.dirs(glue("{dir_project}/2_pipeline"), full.names = FALSE, recursive = FALSE), value=TRUE,ignore.case = T)
folder_marker = grep("marker",list.dirs(glue("{dir_project}/2_pipeline"), full.names = FALSE, recursive = FALSE), value=TRUE,ignore.case = T)
folder_clustering = grep("clustering",list.dirs(glue("{dir_project}/2_pipeline"), full.names = FALSE, recursive = FALSE), value=TRUE,ignore.case = T)

# LOAD FUNCTIONS

q_DefinePanel = function(x,list_type,list_state){
  
  # get column names
  fcs_colname = colnames(x)
  
  # get antigen from description
  antigen = pData(parameters(x[[1]]))$desc
  antigen = sapply(str_split(antigen, " - ",2),"[",1)
  antigen[c(grep(" : ",antigen,value=FALSE,invert=TRUE))]=NA
  antigen = sapply(str_split(antigen, " : ",2),"[",1)
  
  # table of panel
  panel = data.frame(fcs_colname, antigen, row.names = NULL)
  
  # create dict
  values_type = rep("type", length(list_type))
  values_state = rep("state", length(list_state))
  keys = c(list_type,list_state)
  values = c(values_type,values_state)
  names(values) = keys
  
  # add column of marker_class
  panel = panel %>% 
    mutate(marker_class = recode(antigen, !!!values, .default="none")) %>% 
    replace_na(list(marker_class = "none"))
  
  return(panel)
}

q_PrefixSuffixSampleNames = function(fs, add = c("prefix","suffix"), to_add){
  if(add == "prefix"){
    filenames_new = glue("{to_add}{sampleNames(fs)}.fcs")
  }
  
  if(add == "suffix"){
    filenames = str_replace(sampleNames(fs),".fcs","")
    filenames_new = glue("{filenames}{to_add}.fcs")
  }
  
  sampleNames(fs) = filenames_new
  return(fs)
}

q_SplitExtractString = function(string, split_by, number_element){
  
  extracted_string = sapply(strsplit(string, 
                                     split = split_by),
                            "[[",
                            number_element)
  
  return(extracted_string)
}

q_AddExperimentInfo = function(fs){
  
  # find filenames 
  filenames = sampleNames(fs)
  
  # find well id
  well_ids = fsApply(fs, keyword, "$WELLID")
  well_ids = matrix(unlist(well_ids), ncol=1, byrow=FALSE)
  
  # find dates
  dates = fsApply(fs, keyword, "$DATE")
  filedates = matrix(unlist(dates), ncol=1, byrow=FALSE)
  filedates = lubridate::dmy(filedates)
  # sorted dates
  filedates_sorted = unique(sort(filedates))
  
  # define batch
  batches = factor(filedates, levels = filedates_sorted)
  levels(batches) = LETTERS[1:length(unique(batches))]
  
  # make df
  ei = data.frame(name = filenames, 
                  well_id = well_ids,
                  batch = batches, 
                  date = filedates)
  
  rownames(ei) = ei$name
  
  pData(fs) = ei
  
  pData(fs)$sample_id = str_replace(q_SplitExtractString(q_SplitExtractString(pData(fs)$name, "_WLSM",1),(glue("{pData(fs)$well_id}_")),2),"_Fullstain","")
  
  pData(fs) = pData(fs)[,c("date","batch","sample_id","well_id")]
  
  return(fs)
}

q_densityplot = function(parameter, data_ff_fs, ...){
  
  if (is(data_ff_fs, "flowFrame")){
    
    minRange_maxRange_df = range(data_ff_fs, type = "data")
    minRange_maxRange_df_t = data.frame(t(minRange_maxRange_df))
    rownames(minRange_maxRange_df_t) = colnames(minRange_maxRange_df)
    colnames(minRange_maxRange_df_t) = c("minRange","maxRange")
    
    data_ff_fs@parameters@data$minRange = minRange_maxRange_df_t$minRange
    data_ff_fs@parameters@data$maxRange = minRange_maxRange_df_t$maxRange
    
  }
  
  if (is(data_ff_fs, "flowSet")) {
    
    n_flowframes = nrow(pData(data_ff_fs))
    
    for (flowframe in (1:n_flowframes)){
      
      minRange_maxRange_df = range(data_ff_fs[[flowframe]], type = "data")
      minRange_maxRange_df_t = data.frame(t(minRange_maxRange_df))
      rownames(minRange_maxRange_df_t) = colnames(minRange_maxRange_df)
      colnames(minRange_maxRange_df_t) = c("minRange","maxRange")
      
      data_ff_fs[[flowframe]]@parameters@data$minRange = minRange_maxRange_df_t$minRange
      data_ff_fs[[flowframe]]@parameters@data$maxRange = minRange_maxRange_df_t$maxRange
      
    }
    
  }
  
  return(flowViz::densityplot(parameter, data_ff_fs, ...))
  
}

q_plotDR = function(x, color_by, facet_by=NULL, ncol=NULL, labelsize=3, highlight_clust = NULL, background_color = "grey90", label_repel = TRUE,...){
  
  
  # create plot and df
  plt = plotDR(x,"UMAP",color_by = color_by,facet_by=facet_by, ...)
  df = plt$data
  
  # without facet_by  
  if (is.null(facet_by)){
    df = aggregate(x = df[c("x", "y")], 
                   by = list(c = get(color_by,df)), 
                   FUN = median)
    
    no_clusters = length(unique(df$c))
    palette = plt[["plot_env"]][["k_pal"]][1:no_clusters]
    
    if(!is.null(highlight_clust)){
      
      all_clust = unique(df$c)
      grey_clusters = which(!all_clust %in% highlight_clust)
      palette[grey_clusters] = background_color
      
      df$c = ifelse(df$c %in% highlight_clust, as.character(df$c), NA)
      
    }
    
    if(label_repel == TRUE){
      plt + 
        geom_label_repel(aes(x, y, label = c),
                         size=labelsize,
                         col=palette, 
                         data = df, 
                         show.legend = FALSE,
                         inherit.aes = FALSE,
                         seed = 1234) + 
        scale_colour_manual(values=palette)
    } else {
      plt + 
        geom_label(aes(x, y, label = c),
                         size=labelsize,
                         col=palette, 
                         data = df, 
                         show.legend = FALSE,
                         inherit.aes = FALSE) + 
        scale_colour_manual(values=palette)
    }
    
    
    
    # with facet_by    
  } else {
    df = aggregate(x = df[c("x", "y")], 
                   by = list(c = get(color_by,df),f=get(facet_by,df)), 
                   FUN = median)
    
    no_clusters = length(unique(df$c))
    no_facets = length(unique(df$f))
    palette = plt[["plot_env"]][["k_pal"]][1:no_clusters]
    
    if(!is.null(highlight_clust)){
      
      all_clust = unique(df$c)
      grey_clusters = which(!all_clust %in% highlight_clust)
      palette[grey_clusters] = background_color
      
      df$c = ifelse(df$c %in% highlight_clust, as.character(df$c), NA)
      
    }
    
    df = df %>% rename(!!facet_by := "f")
    
    if(label_repel == TRUE){
      
      plt + 
        geom_label_repel(data = df,
                         mapping = aes(x = x, y = y, label = c),
                         size=labelsize,
                         col=rep(palette,no_facets),
                         show.legend = FALSE,
                         inherit.aes = FALSE,
                         seed = 1234) + 
        scale_colour_manual(values=palette)
    } else {
      
      plt + 
        geom_label(data = df,
                         mapping = aes(x = x, y = y, label = c),
                         size=labelsize,
                         col=rep(palette,no_facets),
                         show.legend = FALSE,
                         inherit.aes = FALSE) + 
        scale_colour_manual(values=palette)
    }
    
    
  }
  
  
}

q_plotDR_shuffled = function(sce, ...) {
  
  sce = sce[, sample(ncol(sce))]
  
  plotDR(sce, dr = "UMAP", ...)
  
}

# wrapper function inspired by FlowSOM::GetFlowJoLabels
q_Gating_matrix_aggregates = function(flowSet, gatingTemplate,...){
  
  gates = GatingSet(flowSet)
  cyto_gatingTemplate_apply(gates, gatingTemplate)
  
  # parameter inputs for function below:
  cellTypes = sapply(str_split(gs_get_leaf_nodes(gs),"/"),tail,1)
  
  # MODIFIED VERSION OF FUNCTION FlowSOM::GetFlowJoLabels
  files_in_wsp <- flowWorkspace::sampleNames(gates)
  # get counts manually (these are normally as suffix when using the original function)
  counts2 = fsApply(flowSet, keyword, "$TOT")
  # check the order of files:
  stopifnot(names(counts2) == files_in_wsp)
  counts = as.numeric(unlist(counts2, use.names=FALSE))
  names(counts) = names(counts2)
  files_in_wsp <- gsub("_[0-9]*$", "", files_in_wsp)
  result <- list()
  for(file in files_in_wsp){
    gate_names <- flowWorkspace::gs_get_pop_paths(gates[[file]],
                                                  path = "auto")
    
    gatingMatrix <- matrix(NA,
                           nrow = counts[file],
                           ncol = length(gate_names),
                           dimnames = list(NULL,
                                           gate_names))
    for(gate in gate_names){
      gatingMatrix[, gate] <-
        flowWorkspace::gh_pop_get_indices(gates[[file]], gate)
    }
    
    if(is.null(cellTypes)){
      cellTypes_tmp <- flowWorkspace::gs_get_leaf_nodes(gates[[file]],
                                                        path = "auto")
    } else {
      cellTypes_tmp <- cellTypes
    }
    
    manual <- FlowSOM::ManualVector(gatingMatrix, cellTypes_tmp)
    
    result[[file]] <- list("matrix" = gatingMatrix,
                           "manual" = manual)
    
  }
  gating = result
  
  # get boolean matrix
  matrices = sapply(gating, "[[", 1)
  combined_matrix = do.call(rbind,matrices)
  
  # create aggregate file
  agg = FlowSOM::AggregateFlowFrames(fs_clean,
                                     cTotal =  max(counts)*length(counts),
                                     writeOutput = F,
                                     keepOrder = T,
                                     sampleWithReplacement=F)
  # get aggregate labels
  # original ID col name
  if(length(colnames(agg@exprs)[grep("Original_ID",colnames(agg@exprs))] > 1)){
    
    new_id_col = colnames(agg@exprs)[grep("Original_ID",colnames(agg@exprs))][-1]
    
  } else {
    
    new_id_col = colnames(agg@exprs)[grep("Original_ID",colnames(agg@exprs))][1]
    
  }
  
  cell_types_of_interest <- setdiff(levels(gating[[1]][["manual"]]),"Unlabeled")
  aggregate_labels <- c()
  for (file in unique(exprs(agg)[, "File"])) {
    aggregate_labels <- c(aggregate_labels,
                          as.character(FlowSOM::ManualVector(gating[[file]][["matrix"]],
                                                             cell_types_of_interest)
                                       [exprs(agg)[, new_id_col]
                                         [exprs(agg)[, "File"] == file]]))
  }
  
  stopifnot(nrow(combined_matrix)==sum(counts))
  stopifnot(length(aggregate_labels) ==sum(counts))
  
  result_list = list()
  result_list[["bool_matrix"]] = combined_matrix
  result_list[["aggregate_labels"]] = aggregate_labels
  
  return(result_list)
  
}

# TODO: find better column names for the dataframe created and maybe better name for function? (to make more intuitive)
q_ICSsfcm_GatingMarkers_BackgroundSub_BulkClust = function(sce,k,condition_col,background_value){
  
  n_cells_sce = nrow(colData(sce))
  colnames_marker_pos = setdiff(names(colData(sce))[grep(pattern = "*_pos",names(colData(sce)))],c("roo_pos","root_pos"))
  n_markers_gated = length(colnames_marker_pos)
  n_conditions = length(unique(colData(sce)[[condition_col]]))
  metaname = k
  n_metaclusters = length(unique(colData(sce)[[metaname]]))
  
  # add meta25 to sce
  colData(sce)[metaname] = cluster_ids(sce, metaname)
  
  # make dataframe from sce coldata (rows = n cells in sce)
  df = data.frame(colData(sce)[,(names(colData(sce)) %in% c(condition_col,metaname,colnames_marker_pos))])
  
  #check lenght of df
  stopifnot(nrow(df) == n_cells_sce)
  
  # make the df long format (rows = n cells in sce * 18 gated markers = 2,472,912)
  df_long = df %>% 
    pivot_longer(cols = -c(metaname, condition_col), 
                 names_to = "marker",
                 values_to = "pos")
  # check lenght
  stopifnot(nrow(df_long) == nrow(colData(sce))*n_markers_gated)
  
  # get number of positive cells
  counts_pos = df_long %>% 
    group_by(!!! rlang::syms(metaname),!!! rlang::syms(condition_col),marker) %>% 
    summarise(
      n_pos = sum(pos))
  
  # get cluster sizes
  cluster_sizes = df %>% 
    group_by(!!! rlang::syms(metaname),!!! rlang::syms(condition_col)) %>% 
    tally()
  
  # check the sum is equal to n cells in sce
  stopifnot(sum(cluster_sizes$n) == n_cells_sce)
  
  # combine the counts of positive cells with the cluster sizes 
  df_counts = counts_pos %>% right_join(cluster_sizes, by=c(metaname,condition_col))
  
  # check lenght 
  stopifnot(nrow(df_counts) == n_conditions*n_markers_gated*n_metaclusters)
  
  # calculate pers_pos
  df_counts$perc_pos = (df_counts$n_pos/df_counts$n)*100
  
  # subtract negative value per cluster and marker 
  df_final = df_counts %>%
    group_by(!!! rlang::sym(metaname)) %>%
    mutate(perc_pos_backgroundsubtracted = perc_pos - perc_pos[get({{condition_col}}) == background_value])
  
  # check all negative are 0
  stopifnot(sum(df_final[df_final[condition_col] == background_value,"perc_pos_backgroundsubtracted"])==0)
  
  # remove all stimulation == NEG
  df_final = df_final[df_final[condition_col]!=background_value,]
  
  # recalculate number of positive cells by multiplying per_pos with n cells (--> gets a result relative to cluster size)
  df_final$n_pos_backgroundsubtracted = (df_final$perc_pos_backgroundsubtracted/100)*df_final$n
  
  # set negative values to zero (if not, the pie chart looks weird with blank space)
  df_final$n_pos_backgroundsubtracted = pmax(df_final$n_pos_backgroundsubtracted,0)
  
  # # calculate the cluster distribution for each marker and stimulation
  # dataframe = df_final %>% 
  #   group_by(!!! rlang::syms(condition_col), marker) %>%
  #   mutate(percentage = n_pos_backgroundsubtracted / sum(n_pos_backgroundsubtracted) * 100)
  #  return(dataframe)
  
  return(df_final)
  
}
