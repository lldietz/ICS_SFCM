#### Load source ####
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("0_source.R")

#### Load data #### 

# Panel
panel = read.csv(file = glue("{dir_project}/2_pipeline/{folder_QC}/out/Panel.csv"))[,-1]
features_state = panel[panel$marker_class=="state","antigen"]

# Clean flowSet
fs_clean = readRDS(file = glue("{dir_project}/2_pipeline/{folder_QC}/out/flowSet_clean.rds"))


#### Create gatingSet #### 
gs = GatingSet(fs_clean)

# add metadata
cyto_details(gs)$stimulation = c("GAG","NEG","SEB")
cyto_details(gs)$stimulation = factor(cyto_details(gs)$stimulation,levels = c("NEG","GAG","SEB"))
cyto_details(gs)

# apply gatingTemplate
my_gatingTemplate = "gatingTemplate_marker+.csv"
cyto_gatingTemplate_apply(gs, my_gatingTemplate)


######## DRAW GATES ######## 
# cyto_gate_draw(gs,
#                parent = "root",
#                channels = c("Ki67","SSC-A"),
#                alias = "Ki67+",
#                axes_limits = "data",
#                type = "rectangle",
#                gatingTemplate = my_gatingTemplate)

######## PLOT GATES ######## 
pdf(file = glue("{dir_store}/gates_markers.pdf"), width=10,height=4)
for(feature in features_state){
  channel = panel[(panel$marker_class=="state"&panel$antigen==feature),"fcs_colname"]
  p = cyto_plot(gs,
                parent = "root",
                channels = c(channel,"SSC-A"),
                alias = "",
                axes_limits = "data",
                group_by = "stimulation")
  print(p)

}
dev.off()



png(filename = glue("{dir_project}/3_results/gatesMarkers_IFN.png"),width = 10, height = 3, unit = "in", res = 300)
p = cyto_plot(gs,
              parent = "root",
              channels = c("BV421-A","SSC-A"),
              alias = "",
              axes_limits = "data",
              group_by = "stimulation")
print(p)

dev.off()


######## Save matrix ######## 

gating_result = q_Gating_matrix_aggregates(fs_clean,my_gatingTemplate)
boolean_matrix = gating_result$bool_matrix
write.csv(boolean_matrix, file = glue("{dir_out}/boolean_matrix_marker+.csv"))

