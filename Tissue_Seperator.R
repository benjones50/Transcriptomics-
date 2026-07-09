library(gatepoints)
#now seperating tissues 




output_dir <- "/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/out"
data_dir <- "/projects/nagy_lab_projects/projects_benjones/Filtering Michail Data/data"
raw_data_dir <- "/projects/nagy_lab_projects/projects_benjones/bladder_cancer_spatial_expression/data/spaceranger_count_manual_align_outs/"



#tools to manually seperate multiple tissues on one assay

#type of tissue label determined via earlier plots using gene markers
#used this knowledge to select each tissue, adding tag to data table coordinates

#if old tissue seperation loadable
#loads in gate selections defining tissue boundaries
gates <- readRDS(file = file.path(data_dir, sample_tissue,"tissue_gates.rds"))

left_gate <- gates$left_gate
right_gate <- gates$right_gate



#creates column specifying the coordinates 
object$tissue <- "Unassigned"

object$tissue[as.numeric(left_gate)] <- "left_Tissue"
object$tissue[as.numeric(right_gate)] <- "right_Tissue"




#otherwise, process to make selections

#plots coordinates where there tissue
coords <- GetTissueCoordinates(object)


plot(
  coords$x,
  coords$y,
  pch = 16,
  cex = 0.2,
  asp = 1
)


#creates a data frame that has all of the coordinates of bins 
xy <- data.frame(
  x = coords$x,
  y = coords$y
)



#if first time, or making adjustment:

if (FALSE){ #stops from autorunning
  
  #creates column specifying the coordinates 
  object$tissue <- "Unassigned"
  
  
  left_gate <- fhs(xy)
  object$tissue[as.numeric(left_gate)] <- "left"
  
  
  right_gate <- fhs(xy)
  object$tissue[as.numeric(right_gate)] <- "right"
  

  
  
  #incase of mistakes:
  # To reset left
  object$tissue[object$tissue == "left"] <- "Unassigned"
  
  # To reset right
  object$tissue[object$tissue == "right"] <- "Unassigned"

  
}



#pdf output of tissue seperation 

pdf(file = file.path(sample_output_dir, "tissue_visualization.pdf"), width = 8, height = 8)

#visulizes tissue seperation
tissue_factor <- factor(object$tissue)

plot(
  xy$x,
  xy$y,
  col = tissue_factor,
  pch = 16,
  cex = 0.2,
  asp = 1
)

legend(
  "topright",
  legend = levels(tissue_factor),
  col = 1:length(levels(tissue_factor)),
  pch = 16,
  title = "Tissue",
  bty = "n"
)

dev.off()


#allows us to save the gates to reload in the future, 
#dont have to draw every time

saveRDS(
  list(
    left_gate = left_gate,
    right_gate = right_gate,
  ),
  file = file.path(data_dir, sample_tissue,"tissue_gates.rds")
)