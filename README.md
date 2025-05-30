# SscCube for spot deconvolution analysis
## Software Requirements
1. scCube 2.0.1
1. R 4.3.3
1. scDesign 1.6.0
1. ZINB-WaVE 1.30.0
1. RCTD 2.0
1. Cell2location 0.1.4
2. SpatialDWLS 0.0.1
3. Stereoscope 0.2.0
4. Tangram 1.0.0
5. DestVI 0.20.0
6. Seurat 5.00
7. SPOTlight 0.1.0
8. DSTG 0.1.0
## Reference-based SRT data simulation & Bechmark ##
Code for simulating SRT data using scCube, scDesign3 and ZINB-WaVE is provided in the **SRT_Benchmark**. Two datasetes are provided as references, which are Mouse hippocampus MERFISH data and human DLPFC 10X Visium data. The code for running scDesign3 and ZINB-WaVE are stored at **scDesign3**, **ZINB-WaVE**, respectively. And the 

## Generate SRT data using scRNA-seq data ## 
1. The single-cell RNA-seq data of mouse brain Tabula_Muris_TM_facs_Brain_Non_Myeloid_adata.h5ad as well as the pretrained model Tabula_Muris_TM_facs_Brain_Non_Myeloid_epoch10000.pth were downloaded from https://github.com/ZJUFanLab/scCube/blob/main/tutorial/statistics.md. 

2. The generation pipeline of spatial transcriptomic data with reslutions of n = 5, 10, 20, 30, 50, 100 is stored at **Resolution/Data_Generation/tutorial_resolution.ipynb** which uses the single cell RNA-seq data to simulate gene expression profiles as well as random spatial patterns. 
   
## Spot deconvolution benchmark ##
1. To evaluate the spot deconvolution capability of tools including RCTD, Cell2location, SpatialDWLS, Tangram, DestVI, Seurat, SPOTlight and DSTG, parameters to evaluate including the RMSE, SRCC, JS, PCC, and AS.
2. The performance of the deconvolution tools were analyzed based on the result file of spot deconvolution, **cell_proportions.csv**. Code for benchmarking is deposited at **Resolution/Evaluation/Evaluation.ipynb**.

## Visualization ##
1. Code for visulazation of the evaluation result for SRT data simulation is stored at **R/SRT_evaluation_draw.R**, which generate boxplots for the PCC GEV and MAE.
2. Code for visualization of the evaluation result for deconvolution benchmark is stored at **R/Deconvolution_evaluation.R**, **R/Difference_bar_plot.R**, and **Average_AS.R**, which generate line plot, bar plot and boxplot.
