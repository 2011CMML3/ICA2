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
Code for simulating SRT data using scCube, scDesign3 and ZINB-WaVE is provided in the **SRT_Benchmark**. Two datasetes are provided as references, which are Mouse hippocampus MERFISH data and human DLPFC 10X Visium data.

## Generate SRT data using scRNA-seq data ## 
The single-cell RNA-seq data of mouse brain Tabula_Muris_TM_facs_Brain_Non_Myeloid_adata.h5ad as well as the pretrained model Tabula_Muris_TM_facs_Brain_Non_Myeloid_epoch10000.pth were downloaded from https://github.com/ZJUFanLab/scCube/blob/main/tutorial/statistics.md. 
The generation pipeline of spatial transcriptomic data with reslutions of n = 5, 10, 20, 30, 50, 100 is stored at **Resolution/Data_Generation/tutorial_resolution.ipynb** which uses the single cell RNA-seq data to simulate gene expression profiles as well as random spatial patterns. 
   
## Spot deconvolution benchmark ##
To evaluate the spot deconvolution capability of the above  9 tools, parameters including the RMSE, SRCC, JS, PCC, and AS values were analyzed based on the spot deconvolution file, 'cell_type_proportions.csv'. Code for benchmarking the nine deconvolution methods is deposited at **Resolution/Evaluation/Evaluation.ipynb**.

## Visualization ##
The visulazation of the evaluation result for both SRT data simulation and the deconvolution benchmark are stored at **R**, which generate boxplot, bar plot, line plot generation.
