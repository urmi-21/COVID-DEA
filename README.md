# Differential expression of COVID-19-related genes in European Americans and African Americans
This file will provide instructions on how to reproduce the results described in the pre-print ["Differential expression of COVID-19-related genes in European Americans and African Americans"](https://www.biorxiv.org/content/10.1101/2020.06.09.143271v3).
All analysis were performed using [MetaOmGraph](https://github.com/urmi-21/MetaOmGraph) version 1.8.1. Violin plots were made using ggplot in R.

## Setting up tools and data
* Download MetaOmGraph from [here](http://metnetweb.gdcb.iastate.edu/MetNet_MetaOmGraph.htm). This will download a zip file. Unzip the file to get mog1.8.1.jar file.
* Download the Human cancer RNA-Seq mog project from [here](http://metnetweb.gdcb.iastate.edu/MetNet_MetaOmGraph.htm). This will download a zip file. Unzip the file to get a folder containing three files: Data file, Metadata file and MOG project file (.mog)
* MetaOmGraph user guide is available [here](https://github.com/urmi-21/MetaOmGraph/tree/master/manual)

## Open the project with MetaOmGraph
1. Double click on the .jar file to start MetaOmGraph.
2. Click on open new project and locate the .mog file for the Human cancer RNA-Seq mog project
3. MetaOmGraph will open and display the project.

## Perform Differential Expression Analysis with MetaOmGraph
For more detailed explanation, please go through section 8 of the [MetaOmGraph user manual](https://github.com/urmi-21/MetaOmGraph/tree/master/manual).
1. In the top menubar, go to `Tools --> Differential Expresion Analysis`
2. In the Differential Expression Analysis window, search the groups and perform the analysis.


## Execute R via MetaOmGraph
1. Select the features (genes) from the MetaOmGraph's `Feature Metadata` tab.
2. Go to `Plot --> Selected Rows --> Using R`
3. Browse to the R script in `rscripts/makeplots.R`
4. Enter output directory name. NOTE: the output directory is relative to the project directory.

*NOTE: Please have these packages installed in R:*
* readr
* dplyr
* diptest
* plyr
* scales
* data.table
* ggplot2
* ggthemes
* ggpubr


## Perform Correlation Analysis with MetaOmGraph
1. Apply log_2 transformation: In the main menubar, to `Edit --> Transform data --> Log_2`
2. Select the required row in the MetaOmGraph's `Feature Metadata` tab.
3. Go to `Statistical Analysis --> Correlation` or choose other appropriate option.
4. Click the green button next to `Statistical Analysis` to save the `Feature Metadata` table containing the correlation values.


### References
1. Singh, Urminder, Manhoi Hur, Karin Dorman, and Eve Syrkin Wurtele. "MetaOmGraph: a workbench for interactive exploratory data analysis of large expression datasets." Nucleic acids research 48, no. 4 (2020): e23-e23.

## R Session Info
```
R version 3.6.0 (2019-04-26)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19041)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252 
[2] LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggpubr_0.3.0      ggthemes_4.2.0    ggplot2_3.3.1     data.table_1.12.2
[5] scales_1.0.0      plyr_1.8.4        diptest_0.75-7    dplyr_0.8.4      
[9] readr_1.3.1      

loaded via a namespace (and not attached):
 [1] zip_2.0.3        Rcpp_1.0.2       cellranger_1.1.0 pillar_1.4.4    
 [5] compiler_3.6.0   forcats_0.4.0    tools_3.6.0      digest_0.6.20   
 [9] lifecycle_0.2.0  tibble_3.0.1     gtable_0.3.0     nlme_3.1-139    
[13] lattice_0.20-38  pkgconfig_2.0.2  rlang_0.4.6      openxlsx_4.1.0.1
[17] cli_1.1.0        curl_3.3         haven_2.1.0      rio_0.5.16      
[21] withr_2.1.2      stringr_1.4.0    generics_0.0.2   vctrs_0.3.0     
[25] hms_0.4.2        grid_3.6.0       tidyselect_1.1.0 glue_1.4.1      
[29] R6_2.4.0         rstatix_0.5.0    readxl_1.3.1     foreign_0.8-71  
[33] carData_3.0-2    car_3.0-3        purrr_0.3.2      tidyr_1.0.2     
[37] magrittr_1.5     backports_1.1.4  ellipsis_0.2.0.1 abind_1.4-5     
[41] assertthat_0.2.1 colorspace_1.4-1 ggsignif_0.5.0   stringi_1.4.3   
[45] munsell_0.5.0    broom_0.5.2      crayon_1.3.4    

```
