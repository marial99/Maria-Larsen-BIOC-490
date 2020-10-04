setwd("/Volumes/conormcguinness/Desktop/PhD/flow_cytometry/bd fortessa/220720 BMDC fresh vs frozen/")
library(stringr)
library(CytoExploreR)
##read in fcs files
gs <- cyto_setup("fcsfiles",
                 gatingTemplate = "Activation-gatingTemplate.csv")

##annotate experiment details
cyto_compensate(gs)
cyto_details_edit(gs)

##logicle transform
cyto_transform(gs)
cyto_gatingTemplate_apply(gs)
##draw FSC/SSC gate (using all cells-"root")
cyto_details_edit(gs)
cyto_gate_edit(gs,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"))
##draw Single cell gate from drilled down parent (in this case "cells")
cyto_gate_draw(gs,
               parent = "Cells",
               alias = "Single Cells",
               channels = c("FSC-A","FSC-H"))

##define some fmos to draw gates with
cyto_names(gs)
LDFMO<-cyto_extract(gs, parent = "Single Cells")[[12]]
cyto_gate_draw(gs, parent = "Single Cells",
               alias = c("Dead cells", "Live cells"),
               channels= c("R_780/60-A"), type = "threshold",
               overlay = LDFMO,
               density_fill_alpha =0.5, negate = TRUE)
cyto_plot_save("Rplots/LDplots.png")
cyto_plot(gs, parent = "Single Cells",
               alias = c("Dead cells", "Live cells"),
               channels= c("R_780/60-A"),
               density_fill_alpha =1,
               density_stack = 1,
               density_modal= FALSE, popup=TRUE, legend = TRUE)
##fresh vs frozen
# cyto_plot(gs, parent = "Cells", channels = c("R_780/60-A"), 
#           #legend = TRUE, legend_text = c("Frozen", "Fresh", "Unstained live dead"),
#           # group_by = "Fresh_Frozen",
#          # select = list(FMO=c("N","LD")), 
#           title = "", xlab = "Zombie NIRdye L/D stain Absorbance", density_fill_alpha = 0.5)
cyto_details(cyto_select(gs, FMO = c("N", "LD")))
cyto_names(gs)
CD11FMO<-cyto_extract(gs, parent = "Single Cells")[[9]]
cyto_gate_draw(gs, parent = "Live cells",
               alias = c("CD11c+"),
               overlay = CD11FMO,
               channels= c("BV421-A"), type = "threshold",
               select = list(FMO = c("CD11", "F480")))

cyto_names(gs)
plotorder<-rev(c(15, 9:12, 13, 14, 3:4, 1:2, 7:8, 5:6))
legendtext<-cyto_details(gs)$Sample[plotorder]
legendtextCd11c<-legendtext[-c(2,4,6,8,10,11,12,13)]
cyto_plot_save("Rplots/LivecellsCD11c.png")
cyto_plot(gs[plotorder],
          parent = "Live cells",
          channels= c("BV421-A"), 
          #group_by = c("Timepoint", "Treatment"), 
          alias = c("CD11c+"),
          select = list(Fresh_Frozen = "fresh", FMO = c("N", "CD11", "Unstained")),
          display = 10000,
          density_modal = FALSE, 
          legend = TRUE, 
          legend_text = legendtextCd11c,
          density_fill_alpha = 0.5,
          density_stack = 1.52, popup = TRUE)

CD80FMO<-cyto_extract(gs, parent = "Single Cells")[[10]]
cyto_gate_edit(gs, parent = "Live cells",
               alias = c("CD80+"),
               overlay = CD80FMO,
               channels= c("B_530/30-A"), type = "threshold")
Cd80legendtext<-legendtext[-c(2,4,6,8,10,11,12,14)]
cyto_plot_save("Rplots/LiveCellsCD80.png")
cyto_plot(gs[plotorder], parent = "Live cells",
          alias = c("Cd80+"),
          channels = c("B_530/30-A"), type = "threshold",
          #group_by = c("Timepoint", "Treatment"),
          select = list(Fresh_Frozen = "fresh", FMO = c("N", "CD80", "Unstained")), 
          density_modal = FALSE, 
          legend = TRUE,
          legend_text = Cd80legendtext,
          density_stack = 1.3, density_fill_alpha = 0.5,
          popup = TRUE)
cyto_plot_save("Rplots/LiveCellsCD80_nogate.png")
cyto_plot(gs[plotorder], parent = "Live cells",
          # alias = c("Cd80+"),
          channels = c("B_530/30-A"), type = "threshold",
          #group_by = c("Timepoint", "Treatment"),
          select = list(Fresh_Frozen = "fresh", FMO = c("N", "CD80", "Unstained")), 
          density_modal = FALSE, 
          legend = TRUE,
          legend_text = Cd80legendtext,
          density_stack = 1.3, density_fill_alpha = 0.5,
          popup = TRUE)

cyto_gate_draw(gs, parent = "CD11c+",
               alias = c("CD80+_CD11c+"),
               overlay = CD80FMO,
               channels= c("B_530/30-A"),
               type = "threshold")
cyto_gatingTemplate_edit(gs)
cyto_gate(rename)
cyto_plot_save("Rplots/CD11c_CD80.png")
cyto_plot(gs[plotorder], parent = "CD11c+",
          alias = c("Cd80+"),
          channels = c("B_530/30-A"), type = "threshold",
          #group_by = c("Timepoint", "Treatment"),
          select = list(Fresh_Frozen = "fresh", FMO = c("N", "CD80")), 
          density_modal = FALSE, 
          legend = TRUE,
          legend_text = Cd80legendtext,
          density_stack = 1.25, density_fill_alpha = 0.5,
          popup = TRUE)
cyto_plot_save("Rplots/CD11c_CD80_nogate.png")
cyto_plot(gs[plotorder], parent = "CD11c+",
          # alias = c("CD80+"),
          channels = c("B_530/30-A"), type = "threshold",
          #group_by = c("Timepoint", "Treatment"),
          select = list(Fresh_Frozen = "fresh", FMO = c("N", "CD80")), 
          density_modal = FALSE, 
          legend = TRUE,
          legend_text = Cd80legendtext,
          density_stack = 1.25, density_fill_alpha = 0.5,
          popup = TRUE)
cyto_details_edit(gs)
cyto_names(gs)
