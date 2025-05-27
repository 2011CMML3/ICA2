#####################Fig.2b
data <- read.xlsx('D:/CMML3/ICA2/application/spot_deconvolution/evaluate/all_spot.csv') #Average value per spot for all deconvolution methods obtained from pcc_all_spot.csv, ssim_all_spot.csv, rmse_all_spot.csv,jsd_all_spot.csv
data$Tool <- factor(data$Tool, levels=c('Cell2location', 'DestVI','spatialDWLS','RCTD','Seurat','Tangram','SPOTlight','Stereoscope','DSTG'))
data_filtered <- subset(data, Metrics %in% c("PCC")) #SRCC, RMSE, JS

# Visualization style
color_map <- c(
  "Cell2location" = "#E41A1C",
  "DestVI" = "#377EB8",         
  "spatialDWLS" = "#4DAF4A",    
  "RCTD" = "#984EA3",          
  "Seurat" = "#FF7F00",        
  "Tangram" = "#FFFF33",        
  "SPOTlight" = "#A65628",      
  "Stereoscope" = "#F781BF",    
  "DSTG" = "#999999"            
)

shape_map <- c(
  "Cell2location" = 16,  
  "DestVI" = 15,         
  "spatialDWLS" = 17,    
  "RCTD" = 18,           
  "Seurat" = 19,         
  "Tangram" = 15,        
  "SPOTlight" = 17,      
  "Stereoscope" = 18,    
  "DSTG" = 19            
)

linetype_map <- c(
  "Cell2location" = "solid",
  "DestVI" = "dashed",
  "spatialDWLS" = "dotted",
  "RCTD" = "dotdash",
  "Seurat" = "longdash",
  "Tangram" = "twodash",
  "SPOTlight" = "solid",
  "Stereoscope" = "dashed",
  "DSTG" = "dotted"
)
ggplot(data_filtered, aes(x = n_cell, y = AverageScore, 
                          color = Tool,     
                          shape = Tool,      
                          linetype = Tool,   
                          group = Tool)) +   
  geom_line(size = 0.5) +                    
  geom_point(size = 2) +                  
  scale_color_manual(values = color_map) +  
  scale_shape_manual(values = shape_map) + 
  scale_linetype_manual(values = linetype_map) +  
  labs(x = "Number of Cells (n_cell)", y = "JS Score") +
  theme_classic() +
  theme(legend.position = "bottom")  +      
  theme(
    plot.title = element_text(hjust = 0.5),
    #panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),  
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 10,color = 'black'),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = 'bottom'
    )

 # scale_y_continuous(limits = c(0, 0.9))  


