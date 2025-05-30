###################################Fig.2d
data <- read.xlsx('D:/CMML3/ICA2/application/spot_deconvolution/evaluate/Diff_pcc_srcc.xlsx')
first_plot_order <- c(
  "RCTD", "Tangram", "Cell2location", "spatialDWLS", "DSTG", 
  "Seurat", "SPOTlight", "DestVI", "Stereoscope"
)
data$X1 <- factor(data$X1, levels=rev(first_plot_order))
color_map <- c(
  "RCTD" = "#984EA3",         
  "Tangram" = "#FFFF33",        
  "Cell2location" = "#E41A1C",  
  "spatialDWLS" = "#4DAF4A",   
  "DSTG" = "#999999" ,           
  "Seurat" = "#FF7F00",         
  "SPOTlight" = "#A65628",      
  "DestVI" = "#377EB8",         
  "Stereoscope" = "#F781BF"   
)
#data,aes(y = reorder(X1, -differ_SRCC), x = differ_SRCC,fill = X1)
ggplot(data,aes(y = X1, x = differ_SRCC,fill = X1))+
  geom_bar(stat = "identity", position = position_dodge())+
  scale_fill_manual(values = color_map) +
  theme_classic() +
  theme(legend.position = "bottom")  +      
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 10,color = 'black'),  
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = 'bottom'
  )




