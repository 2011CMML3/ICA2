############################Fig. 2e
data <- read.xlsx('D:/2025春夏学期/CMML3/ICA2/R/application/spot_deconvolution/evaluate/AS.xlsx',colNames = F)
data$X1 <- factor(data$X1, levels=c('Cell2location', 'DestVI','spatialDWLS','RCTD','Seurat','Tangram','SPOTlight','Stereoscope','DSTG'))
ggplot(data,aes(y = reorder(X1, X2), x = X2,fill = X1)) +
  stat_boxplot(geom = "errorbar",
               width=0.3)+
  geom_boxplot(outlier.shape = NA, outlier.size = 0) +
  #facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  #scale_x_continuous(breaks = seq(0,0.6, by = 0.2), limits = c(0,0.6))+
  scale_fill_manual(values = color_map)+
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 10,color = 'black'),  
    axis.text.y = element_text(size = 10, color = "black")
  )
