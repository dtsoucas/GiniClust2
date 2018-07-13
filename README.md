## GiniClust2

GiniClust2 is a clustering algorithm for the simultaneous detection of common and rare cell types from single-cell gene expression data.  It uses a novel cluster-aware weighted consensus clustering algorithm to combine GiniClust and Fano-based k-means clustering results, by maximizing the strengths of these individual clustering methods in detecting rare and common clusters, respectively.  The full manuscript is available at: https://doi.org/10.1186/s13059-018-1431-3.

GiniClust2 is written in the R programming language.  To use, please download the **GiniClust2\_download** folder and follow the instructions provided in the **Reference\_Manual.pdf**.

The following additional folders are included, which contain GiniClust2 applications to several data sets: 

   - **Main**: This folder contains the following main scripts for running GiniClust2 on four simulated and real data sets:
   
      - **Main_Simulation.R**: Data simulated from intestinal cells [1], to create both common and rare cell types.
      - **Main_10X_subsample.R**: Data subsampled at various proportions from 10X PBMC data [2], to create 140 datasets with cell types of varying rarity.
      - **Main_inDrop_day4.R**: inDrop data from mouse embryonic stem cells 4 days post-LIF [3].
      - **Main_10X_full.R**: The full 68k-cell 10X PBMC data set [2].
      
   - **Proj**: This folder contains four subfolders for each of the data sets mentioned above.  Each of these contain data, results, and figures from the GiniClust2 analyses.
   
   - **Rfunction**: This folder contains the R scripts for GiniClust2, called in the **Main** files.


###### Data sources:<br />
<sub> [1] Grün, D., Lyubimova, A., Kester, L., Wiebrands, K., Basak, O., Sasaki, N., et al. (2015). Single-cell messenger RNA sequencing reveals rare intestinal cell types. Nature, 525(7568), 251-255. <br />
[2] Zheng, G. X., Terry, J. M., Belgrader, P., Ryvkin, P., Bent, Z. W., Wilson, R., et	al. (2017). Massively parallel digital transcriptional profiling of single	cells. Nat Commun, 8, 14049. <br />
[3] Klein, A. M., Mazutis, L., Akartuna, I., Tallapragada, N., Veres, A., Li, V., et al. (2015). Droplet barcoding for single-cell transcriptomics applied to embryonic stem cells. Cell, 161(5), 1187-1201.
