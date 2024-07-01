#### Load source ####
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("0_source.R")

#### Load data #### 

## Panel
panel = read.csv(file = glue("{dir_project}/2_pipeline/{folder_QC}/out/Panel.csv"))[,-1]

## Clean flowSet
fs_clean = readRDS(file = glue("{dir_project}/2_pipeline/{folder_QC}/out/FlowSet_clean.rds"))
pData(fs_clean)


#### Create gatingSet #### 
gs = GatingSet(fs_clean)

# apply gatingTemplate
my_gatingTemplate = "gatingTemplate_celltype.csv"
cyto_gatingTemplate_apply(gs, my_gatingTemplate)

######## DRAW GATES ######## 

# ##### /CD3+/CD4-CD8-
# ##### /CD3+/CD4+CD8+
# ##### /CD3+/CD4 T cells
# ##### /CD3+/CD8 T cells
# cyto_gate_draw(gs,
#                parent = "root",
#                channels = c("BUV496-A","SparkBlue-574-A"),
#                alias = c("CD4-CD8-","CD4+CD8+","CD4 T cells", "CD8 T cells"),
#                type = c("rectangle","rectangle","rectangle","rectangle"),
#                contour_lines = 15,
#                gatingTemplate = my_gatingTemplate,
#                axes_limits = "data"
# )
#
# ##### /CD3+/CD4+/Treg
# ##### /CD3+/CD4+/FoxP3-
# cyto_gate_draw(gs,
#                parent = "CD4 T cells",
#                channels = c("SparkNIR-685-A","SSC-A"),
#                alias = c("Treg","FoxP3-"),
#                type = c("rectangle", "rectangle"),
#                contour_lines = 15,
#                gatingTemplate = my_gatingTemplate,
#                axes_limits = "data")
#
##### /CD3+/CD4+/FoxP3-/EM
##### /CD3+/CD4+/FoxP3-/CM
##### /CD3+/CD4+/FoxP3-/CD45RA+CCR7+
##### /CD3+/CD4+/FoxP3-/TEMRA
# cyto_gate_draw(gs,
#                parent = "FoxP3-",
#                channels = c("[RB780]-A","BV785-A"),
#                alias = c("CD4 Effector memory", "CD4 Central memory", "CD4 CD45RA+CCR7+","CD4 TEMRA"),
#                type = c("rectangle","rectangle","rectangle","rectangle"),
#                contour_lines = 15,
#                gatingTemplate = my_gatingTemplate,
#                axes_limits = "data"
#                )
# 
# ##### /CD3+/CD8+/EM
# ##### /CD3+/CD8+/CM
# ##### /CD3+/CD8+/CD45RA+CCR7+
# ##### /CD3+/CD8+/TEMRA
# cyto_gate_draw(gs,
#                parent = "CD8 T cells",
#                channels = c("[RB780]-A","BV785-A"),
#                alias = c("CD8 Effector memory", "CD8 Central memory", "CD8 CD45RA+CCR7+","CD8 TEMRA"),
#                type = c("rectangle","rectangle","rectangle","rectangle"),
#                contour_lines = 15,
#                gatingTemplate = my_gatingTemplate,
#                axes_limits = "data"
# )
# 
##### /CD3+/CD4+/FoxP3-/CD45RA+CCR7+/Naïve
##### /CD3+/CD4+/FoxP3-/CD45RA+CCR7+/SCM
# cyto_gate_draw(gs,
#                parent = "CD4 CD45RA+CCR7+",
#                channels = c("BV650-A","PE-Fire640-A"),
#                alias = c("CD4 Naïve","CD4 SCM"),
#                type = c("rectangle","rectangle"),
#                contour_lines = 15,
#                gatingTemplate = my_gatingTemplate,
#                axes_limits = "data"
# )
# 
# ##### /CD3+/CD8+/CD45RA+CCR7+/Naïve
# ##### /CD3+/CD8+/CD45RA+CCR7+/SCM
# cyto_gate_draw(gs,
#                parent = "CD8 CD45RA+CCR7+",
#                channels = c("BV650-A","PE-Fire640-A"),
#                alias = c("CD8 Naïve","CD8 SCM"),
#                type = c("rectangle","rectangle"),
#                contour_lines = 15,
#                gatingTemplate = my_gatingTemplate,
#                axes_limits = "data"
# )


######## PLOT GATES ######## 

pdf(file = glue("{dir_store}/gating_scheme.pdf"),height = 4, width = 18)
cyto_plot_gating_scheme(gs,
                        layout = c(1,6))
dev.off()

png(file = glue("{dir_store}/gatingtree.png"),height = 4, width = 5,unit = "in",res = 300)
plot(gs)
dev.off()


######## Get annotation vector ######## 

gating_results = q_Gating_matrix_aggregates(fs_clean,my_gatingTemplate)

aggregate_labels = gating_results$aggregate_labels

saveRDS(aggregate_labels, file = glue("{dir_out}/aggregate_labels.rds"))

