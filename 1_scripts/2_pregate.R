#### Load source ####
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("0_source.R")

#### Load data #### 
fs_transform = readRDS(file = glue("{dir_project}/2_pipeline/{folder_transform}/out/FlowSet_transform.rds"))
pData(fs_transform)$name = rownames(pData(fs_transform))

#### Create gatingSet #### 
gs = GatingSet(fs_transform)
cyto_details(gs)$name = rownames(cyto_details(gs)) # makes sure there are no discrepancies between names

# apply gatingTemplate
my_gatingTemplate = "gatingTemplate_pregating.csv"
cyto_gatingTemplate_apply(gs, my_gatingTemplate)

######## DRAW GATES ######## 

######## Lymphocytes ########

cyto_gate_draw(gs,
               parent = "root",
               channels = c("FSC-A", "SSC-A"),
               alias = c("Lymphocytes"),
               type = c("polygon"),
               contour_lines = 15,
               gatingTemplate = my_gatingTemplate,
               xlim = c(0,800000),
               ylim = c(0,2000000)
               )

cyto_plot(gs[1:3],
          parent = "root",
          channels = c("FSC-A", "SSC-A"),
          alias = "",
          xlim = c(0,800000),
          ylim = c(0,2000000),
          layout = c(1,3))

######## Singlets ########

cyto_gate_draw(gs,
               parent = "Lymphocytes",
               channels = c("FSC-A", "FSC-H"),
               alias = c("Singlets"),
               type = c("polygon"),
               contour_lines = 15,
               gatingTemplate = my_gatingTemplate,
               axes_limits = "data"
)

cyto_plot(gs[1:3],
          parent = "Lymphocytes",
          channels = c("FSC-A", "FSC-H"),
          alias = "",
          axes_limits = "auto",
          layout = c(1,3))


######## Live ########

cyto_gate_draw(gs,
               parent = "Singlets",
               channels = c("LiveDeadFixableBlue-A", "SSC-A"),
               alias = c("Live"),
               type = c("rectangle"),
               contour_lines = 15,
               gatingTemplate = my_gatingTemplate,
               axes_limits = "auto"
)

cyto_plot(gs[1:3],
          parent = "Singlets",
          channels = c("LiveDeadFixableBlue-A", "SSC-A"),
          alias = "",
          axes_limits = "auto",
          layout = c(1,3))

######## Dump ########

cyto_gate_draw(gs,
               parent = "Live",
               channels = c("BUV563-A", "SSC-A"),
               alias = c("Dump-"),
               type = c("rectangle"),
               contour_lines = 15,
               gatingTemplate = my_gatingTemplate,
               axes_limits = "auto"
)

cyto_plot(gs[1:3],
          parent = "Live",
          channels = c("BUV563-A", "SSC-A"),
          alias = "",
          axes_limits = "auto",
          layout = c(1,3))


######## CD3+CD19- ########

cyto_gate_draw(gs,
               parent = "Dump-",
               channels = c("BUV615-A", "BV570-A"),
               alias = c("CD3+CD19-"),
               type = c("rectangle"),
               contour_lines = 15,
               gatingTemplate = my_gatingTemplate,
               axes_limits = "auto"
)

cyto_plot(gs[1:3],
          parent = "Dump-",
          channels = c("BUV615-A", "BV570-A"),
          alias = "",
          axes_limits = "auto",
          layout = c(1,3))

######## PLOT GATES ########

height = 4
width = 10

png(filename = glue("{dir_store}/Gatingscheme.png"), height = height, width = width, unit = "in", res = 400)
cyto_plot_gating_scheme(gs,
                        layout = c(1,5),
                        axes_limits = "instrument")
dev.off()

png(filename = glue("{dir_store}/gatingtree.png"), height = height, width = width, unit = "in", res = 400)
plot(gs)
dev.off()

png(filename = glue("{dir_store}/Gating1_lymphocytes.png"), height = height, width = width, unit = "in", res = 400)
cyto_plot(gs[1:3],
          parent = "root",
          channels = c("FSC-A", "SSC-A"),
          alias = "",
          xlim = c(0,800000),
          ylim = c(0,2000000),
          layout = c(1,3),
          display = 1)
dev.off()



png(filename = glue("{dir_store}/Gating2_singlets.png"), height = height, width = width, unit = "in", res = 400)
cyto_plot(gs[1:3],
          parent = "Lymphocytes",
          channels = c("FSC-A", "FSC-H"),
          alias = "",
          axes_limits = "auto",
          layout = c(1,3),
          display = 1)

dev.off()

png(filename = glue("{dir_store}/Gating3_live.png"), height = height, width = width, unit = "in", res = 400)
cyto_plot(gs[1:3],
          parent = "Singlets",
          channels = c("LiveDeadFixableBlue-A", "SSC-A"),
          alias = "",
          axes_limits = "auto",
          layout = c(1,3))
dev.off()

png(filename = glue("{dir_store}/Gating4_Dump-.png"), height = height, width = width, unit = "in", res = 400)
cyto_plot(gs[1:3],
          parent = "Live",
          channels = c("BUV563-A", "SSC-A"),
          alias = "",
          axes_limits = "auto",
          layout = c(1,3))
dev.off()

png(filename = glue("{dir_store}/Gating5_CD3+CD19-.png"), height = height, width = width, unit = "in", res = 400)
cyto_plot(gs[1:3],
          parent = "Dump-",
          channels = c("BUV615-A", "BV570-A"),
          alias = "",
          axes_limits = "auto",
          layout = c(1,3))
dev.off()


######## SAVE ########

fs_pregated = cytoset_to_flowSet(gs_pop_get_data(gs, "CD3+CD19-"))
fs_pregated = q_PrefixSuffixSampleNames(fs_pregated, add = "suffix", to_add = "_pregated")
pData(fs_pregated)$name = rownames(pData(fs_pregated))

dir_write = glue("{dir_store}/FCS_pregated")
dir.create(dir_write, showWarnings = FALSE)
write.flowSet(fs_pregated, outdir = dir_write)

saveRDS(fs_pregated, glue("{dir_out}/FlowSet_pregated.rds"))



