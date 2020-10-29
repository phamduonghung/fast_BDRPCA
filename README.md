# fast BD-RPCA
 

This MATLAB package is a collection of scripts allowing to generate figures in the paper [1]. This paper introduces a computationally efficient technique for estimating high-resolution Doppler blood flow from an ultrafast ultrasound image sequence. More precisely, it consists in a new fast alternating minimization algorithm that implements a blind deconvolution method based on robust principal component analysis. Numerical investigation carried out on in vivo data shows the efficiency of the proposed approach in comparison with state-of-the-art methods



## Instructions: 
1. Download the package in .zip file (click green Code above) and then unzip it. Note that the **name** of the unzipped folder should be **fast_BDRPCA**.  
2. Set **Current Folder** of MATLAB to this unzipped folder, i.e. **fast_BDRPCA**.  
3. Download all "simulation" data from the following link: 
https://cloud.irit.fr/index.php/s/PVcA9S1OyiCRcZW and then put them into the folder **Data**
4. Run each file **Fig?.m** corresponding to each figure (from Fig. 1 to Fig. 2a-2e) in [1]. 
5. To print nice pdf figures, the **export_fig** package was used, which required a software support downloaded from the following link (there is also a portable version of this software): https://www.ghostscript.com/download/gpcldnld.html. In the codes, change **FigFeatures.print= 1** if you want to print the .pdf figure when using this package. 


[1] D.-H. Pham, A. Basarab, JP. Remenieras, P. Rodriguez and D. Kouame," Fast High Resolution Blood Flow Estimation and Clutter Rejection via an Alternating Optimization Problem," *Submitted to ISBI 2021*, Nice, France.

COPYRIGHT

fast BD-RPCA is copyright reserved. If any issue related this package appears, please contact Duong Hung PHAM at duong-hung.pham@irit.fr.
