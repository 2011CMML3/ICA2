#############Fig.2c
library(ggplot2)
library(openxlsx)

#PCC GEV

data <- read.xlsx('D:/CMML3/ICA2/R/zinb_scdesign/PCC.xlsx')
ggplot(data, aes(x = method, y = pcc, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +         
  #geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +    
  labs(
    #title = "PCC Distribution by Method Group",
    x = "Method",
    y = "PCC GEV",
    fill = "Method"
  ) +
  scale_fill_brewer(palette = "Set3") +                  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),   
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),  
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black") 
  )


# MAE

data2 <- read.xlsx('D:/CMML3/ICA2/R/zinb_scdesign/MAE.xlsx')
ggplot(data2, aes(x = method, y = mae, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +          
  #geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +    
  labs(
   # title = "MAE Distribution by Method Group",
    x = "Method",
    y = "MAE GEV",
    fill = "Method"
  ) +
  scale_fill_brewer(palette = "Set3") +                   
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),    
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(colour = "black")  
  )


    legend.text.align = 1
  )

