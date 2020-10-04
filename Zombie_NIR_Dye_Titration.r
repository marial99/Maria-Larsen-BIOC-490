library(stringr)
library(CytoExploreR)
##read in fcs files
gs <- cyto_setup("050820NIRDYE BMDC titration/",
                 gatingTemplate = "Activation-gatingTemplate.csv")

##annotate experiment details

##logicle transform
cyto_transform(gs)
cyto_details_edit(gs)
cyto_gatingTemplate_apply(gs)
##draw FSC/SSC gate (using all cells-"root")
cyto_gate_draw(gs,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"))
##draw Single cell gate from drilled down parent (in this case "cells")
cyto_gate_draw(gs,
               parent = "Cells",
               alias = "Single Cells",
               channels = c("FSC-A","FSC-H"))
cyto_details(gs)
cyto_plot_save("Rplots/LDtitration.png")
cyto_plot(gs,
          parent = "Single Cells",
          channels= c("R_780/60-A"), group_by = c("Stain"),
          density_modal = FALSE, 
          legend = TRUE, 
          density_fill_alpha = 0.5,
          density_stack = 1.8, popup = TRUE)
