
##own analysis
setwd("./fcsfiles/")
library(CytoExploreR)
gs <- cyto_setup(".",
                 gatingTemplate = "Activation-gatingTemplate.csv")
cyto_details_edit(gs)
##compensates using values generated on facsdiva
cyto_compensate(gs)
##transforms the fluorescent gates to logicle format
cyto_transform(gs)
cyto_gatingTemplate_apply(gs)
##draw the FSC-H and SSC-H plot
cyto_gate_edit(gs,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"))

cyto_gate_edit(gs,
               parent = "Cells",
               alias = "Single Cells",
               channels = c("FSC-A","FSC-H"))
##list the names of the samples
cyto_names(gs)
##cyto_extract "pulls out" one sample
LDFMO<-cyto_extract(gs, parent = "Single Cells")[[3]]
##LD gate draw
cyto_gate_edit(gs,
               parent = "Single Cells",
               alias = "Live cells",
               channels = c("R_780/60-A"),
               overlay = LDFMO,
               type = "interval")
##live dead plot (stacked)

library(forcats)
##order for this plot

#define the order to plot the samples-plots from the bottom up (hence rev())
plotorder<-rev(c(10,1,2,3,9,11,7,5,6,4,12,8))
legendtext<-cyto_details(gs)$Sample[plotorder]
legendtextLD<-legendtext[-c(11,10)]

cyto_plot_save("../Rplots/LDplot.png")
cyto_plot(gs[plotorder],
               parent = "Single Cells",
               alias = "Live cells",
               channels = c("R_780/60-A"),
          select = list(FMO=c("No", "LD", "Unstained")),
          density_modal = FALSE, 
          label_position = (-1000),
          legend = TRUE,
          legend_text =  legendtextLD,
          density_fill_alpha = 1,
          xlab = "Zombie NIRdye Live Dead Stain R_780/60-A",
          density_stack = 1.3, popup = TRUE)
CD8FMO<-cyto_extract(gs, parent = "Live cells")[[1]]
##plot gate for CD8 positive cells
cyto_gate_edit(gs, parent = "Live cells",
               alias = "CD8+",
               channels = c("R_670/14-A"),
               overlay = CD8FMO,
               type = "threshold")

legendtextCD8<-legendtext[-c(9,10)]
cyto_plot_save("../Rplots/CD8plot.png")
cyto_plot(gs[plotorder], parent = "Live cells", alias = "",
          channels=c("R_670/14-A"), 
          #group_by = c("Sample"),
          select = list(FMO=c("No", "CD8", "Unstained")),
          density_modal = FALSE, 
          density_stack = 1.5,
          density_fill_alpha = 0.5,
          legend = TRUE,
          legend_text = legendtextCD8,
          popup = TRUE)

PHA<-cyto_extract(gs, parent = "Live cells")[[9]]
cyto_gate_edit(gs, parent = "Live cells",
               alias = "Dividing cells",
               channels = c("B_530/30-A"),
               overlay = PHA,
               density_fill_alpha = 0.5)
cyto_gate_edit(gs, parent = "CD8+",
               alias = "CD8 Dividing cells",
               channels = c("B_530/30-A"),
               overlay = PHA)
legendtextLivCFSE<-legendtext[-c(9)]
cyto_plot_save("../Rplots/LivecellsCFSEplot_with_CD8_FMO.png")
cyto_plot(gs[plotorder], parent = "Live cells", 
          alias = "Dividing cells",
          channels=c("B_530/30-A"), 
          #group_by = c("Sample"),
          select = list(FMO=c("No", "CFSE", "Unstained", "CD8")),
          density_modal = FALSE, density_stack = 1.3,
          density_fill_alpha = 1, popup = TRUE, 
          label_text = rep("", 10),
          legend = TRUE,
          legend_text = legendtextLivCFSE)
cyto_plot_save("../Rplots/CD8cellsCFSEplot.png")
plotorder
cyto_names(gs)[plotorder]
cyto_plot(gs[plotorder], parent = "CD8+", 
          alias = "",
          channels=c("B_530/30-A"), 
          #group_by = c("Sample"),
          select = list(FMO=c("No", "CFSE", "Unstained")),
          density_modal = FALSE, density_stack = 1.3,
          density_fill_alpha = 1, popup = TRUE, 
          label_text = rep("", 10),
          legend = TRUE,
          legend_text = legendtextLivCFSE)

cyto_gate_edit(gs,
               parent = "Live cells",
               alias = c("CD8+ dividing", "CD8+ cytostatic",
                         "CD8- cytostatic", "CD8- dividing"),
               channels=c("B_530/30-A", "R_670/14-A"),
               type = "quadrant")
legendtext2Dplot<-legendtext
legendtext2Dplot[9]<-"UNSTIM"
cyto_plot_save("../Rplots/2Dplotallgates.png", width = 15)
cyto_plot(gs[rev(plotorder)], parent = "Live cells",
         alias = "",
         channels=c("B_530/30-A", "R_670/14-A"), 
         #group_by = c("Sample"),
         select = list(FMO=c("No", "CFSE", "Unstained", "LD", "CD8")),
         xlab ="",
         ylab = "",
         #axes_text = c("TRUE", "FALSE"),
         title = rev(legendtext2Dplot),
         layout = c(3,4),
         label_position = "manual",
         popup = TRUE)
cyto_plot_save("../Rplots/2Dplotnogates.png", width = 15)
cyto_plot(gs[rev(plotorder)], parent = "Live cells",
          channels=c("B_530/30-A", "R_670/14-A"), 
          #group_by = c("Sample"),
          select = list(FMO=c("No", "CFSE", "Unstained")),
          xlab ="",
          ylab = "",
          axes_label_text_size = 1,
          axes_text = c("TRUE", "FALSE"),
          title = rev(legendtextLivCFSE),
          layout = c(2,5),
          popup = TRUE)
cyto_plot_save("../Rplots/2Dplotdivgatesonly.png", width = 15)
cyto_plot(gs[rev(plotorder)], parent = "Live cells",
          alias = (c("CD8- dividing", "CD8+ dividing")),
          channels=c("B_530/30-A", "R_670/14-A"), 
          #group_by = c("Sample"),
          select = list(FMO=c("No", "CFSE", "Unstained")),
          xlab ="",
          ylab = "",
          #axes_text = c("TRUE", "FALSE"),
          title = rev(legendtextLivCFSE),
          layout = c(2,5),
          popup = TRUE)
plotorder
par("mar")
par(mar = c(1,1,1,1))
