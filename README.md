# Thesis upplimentary---Codes

This repository contains 2 folders for imaging processing and scRNAseq analysis respectively. 

For scRNAseq analysis, the 7 files (written in R) in the folder cover the full pipeline of scRNAseq analysis involved in the thesis. File 1~3 processes the scRNAseq database, which covers quality control/cell filtering (file 1), integration (file 2) and clustering (file 3). Files 4~7 aim at analysing differentially expressed genes (DEG) based on pseudo bulk RNAseq databases. File 4 creates pseudo bulk RNAseq count matrices and analyses expression fold changes. Based on the output of file 4, file 5 compares the expression of cluster marker genes between WT and mutant; file 6 performs GO analysis and KEGG pathway analysis; file 7 analyses pan-habenular gene expression in WT and mutant 

The pipeline was built based on published algorithms/pipeline of Seurat for general scRNAseq processing and pseudo bulk count matrix generation, DEseq2 for pseudo bulk DEG analysis, and Bioconductor zebrafish genome-wide annotation database for GO and KEGG analysis. Other algorithms were incorporated for specific purposes. These algorithms include miQC, DoubletFinder, Scater, DAseq, Harmony, and tascCODA. The files are self-explanatory. Regarding how the whole pipeline works, please read the annotations in each file and the documentation for published algorithms.

The raw input should be 4 Seurat objects corresponding to WT and mutant scRNAseq count matrices (2 bio replicates for each genotype). A published zebrafish larval telencephalon database is required for file 7. Files need to run in order 1 to 7 since the output of the previous file is required for the following files. 




For image processing, 5 files are included for different purposes:

“Integrated signal analysis” is built to measure gene expression in the left and right habenulae.  This algorithm is applied to measure asymmetric expression of pou4f1/pou4f2 in WT, and expression of tac3a and kctd12.2 in different genotypes (to study synergistic interaction between pou4f1 and pou4f2). This algorithm takes a 3D confocal image of whole-mount fluorescent staining. The algorithm creates a mask for the fluorescent signal by thresholding the signal intensity based on the average brightness of the image. The output mask is then imposed onto the original confocal image to remove autofluorescence and noise. The output file is then used to calculate integrated fluorescent intensity. 

To run the code, 4 folders need to be built. The first folder should contain 3D (only one channel of fluorescent signal is permitted) microscopy images written in TIF format. The second folder should be named “masks”. This folder should be empty before running the code, and the intermediate output (masks for fluorescent signal) will be written in this folder.  The third folder should be empty as well. The masked fluorescent signal will be written to this folder. Integrated fluorescent intensity will be automatically summarised in csv file. The left and right halves of the image and the total image will be calculated for each sample. It is recommended to manually crop the final output image to remove autofluorescence and noise at the edges of the image, then store the cropped image in the 4th folder named “cropped masked signal”.

“Habenula cell number measurement” is built to count cell numbers of a 2D image of a DAPI-stained sample. This algorithm is adapted from the first half of the standard watershed cell segmentation workflow. What the algorithm does is: 1) masks the cells, 2) seeds the cells based on masks. 3) Use the manually drawn border as a mask to exclude detected cells that are not habenula. 4) The seeding result is overlaid onto to DAPI image, then manually corrected. 5) The number of seeds (therefore cell number) is counted. To run the algorithms, DAPI staining of each sample needs to be prepared in one folder (named “slices”) in a TIF file format. For this algorithm, 3 planes were manually selected for each sample. Thus, a TIF file for each sample was designed to be a 3D array with the shape of 3 X image width X image height. To mark the border of the habenula, masks of the habenula were manually drawn and stored as TIF files with the same shape in the folder “Hb mask”. The following empty folder needs to be built to store intermediate output files before running the code: a. “cell mask” to store an image of cell masks; b. “detected cells” to store automatically detected cells (seedings) on the entire slice, including non-habenular cells; c. “seed_in_Hb” to store seeding results of only habenular cells. 

The 3 files named as “Hb06 measurement 1/2/3_***” are built to measure the population size of Hb06. The idea is to mask the area co-stained by foxa1-gng8 or pou3f1-pnoca. The first file will generate masks of cell membrane based on DAPI-negative area. In file 2, the fluorescent signal of Hb06 marker gene is masked by thresholding based on auto-fluorescent intensity (i.e., signal strength at cell membrane) for each Z plane of the 3D image. Masks of 2 marker genes co-stained for the same sample are multiplied to generate masks of areas that are co-stained by both markers. This area reflects where Hb06 is located, and the pixel number of this area was calculated. To run the pipeline, one folder needs to be created for each sample. For each sample, a multi-channel image was taken (DAPI channel, channel for marker gene 1, channel for marker gene 3). Split the 3 channels and store them as three 3D TIF files in a folder for each sample. Intermediate outputs will also be stored in the same folder. The code needs to be run one by one from 1 to 3. 












