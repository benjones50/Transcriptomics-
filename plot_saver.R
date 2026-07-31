#helper function to save plots easily bc its very repetitive

save_plot <- function(
    plot_object,           #plot
    file_stub,             #name of plot
    save_dir = output_dir, #where saving
    width = 16,            #dimensions
    height = 6,
    dpi = 600,              #dpi for png
    save_pdf = FALSE,
    limitsize = FALSE
) {
  
  #makes directories for where saving
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (save_pdf){
  #saves pdf
  message("saving pdf: ", file_stub)
  ggsave(
    filename = file.path(save_dir, paste0(file_stub, ".pdf")),
    plot = plot_object,
    width = width,
    height = height,
    limitsize = limitsize
  )
  }
  
  #saves png
  message("saving png: ", file_stub)
  ggsave(
    filename = file.path(save_dir, paste0(file_stub, ".png")),
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    limitsize = limitsize,
  )
  message("saved: ", file_stub)
}


#use example:

# save_plot(
#   qc_plot, #plot
#   "QC_scatter_plot", #name of plot
#   save_dir = sample_output_dir, #where saving
#   width = 10, 
#   height = 5,
#   dpi = 600 
# ) 









