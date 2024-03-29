# fast BD-RPCA
 
This MATLAB package is a collection of scripts allowing to generate figures in the paper [1], but for a **simuation** dataset. This paper introduces a computationally efficient technique for estimating high-resolution Doppler blood flow from an ultrafast ultrasound image sequence. More precisely, it consists in a new fast alternating minimization algorithm that implements a blind deconvolution method based on robust principal component analysis. Numerical investigation carried out on in vivo data shows the efficiency of the proposed approach in comparison with state-of-the-art methods


## Instructions
1. Download the package in .zip file (click green Code above) and then unzip it. Note that the **name** of the unzipped folder should be **fast_BDRPCA**.  
2. Set **Current Folder** of MATLAB to this unzipped folder, i.e. **fast_BDRPCA**.  
3. Download all **simulation** data from the following link: 
https://cloud.irit.fr/index.php/s/846gUKURnYbehVl and then put them into the folder **Data**.
4. Run each file **Fig?.m** corresponding to each figure (from Fig. 1a to Fig. 1d) in [1], but for a **simulation** dataset. 
5. The Matlab code for a fast version of DRPCA (non-blind) was also provided in this package whose *running file* is **Code_fast_DRPCA.m**. 
6. To print nice pdf figures, the **export_fig** package was used, which required a supporting software downloaded from the following link (there is also a portable version of this software): https://www.ghostscript.com/download/gpcldnld.html. In the codes, change **FigFeatures.print= 1** if you want to print the .pdf figure using this package. 


[1] D.-H. Pham, A. Basarab, JP. Remenieras, P. Rodriguez and D. Kouame," Fast High Resolution Blood Flow Estimation and Clutter Rejection via an Alternating Optimization Problem," in *ISBI 2021*, Nice, France. Available: https://arxiv.org/pdf/2011.01811.pdf.

## COPYRIGHT

RobustPCA_Doppler.m comes from https://github.com/dlaptev/RobustPCA.

fast BD-RPCA is copyright reserved. If any issue related this package appears, please contact Duong Hung PHAM at duong-hung.pham@irit.fr.
