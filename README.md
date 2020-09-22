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


## Execute R via MetaOmGraph (violin plots)
1. Select the features (genes) from the MetaOmGraph's `Feature Metadata` tab.
2. Go to `Plot --> Selected Rows --> Using R`
3. Browse to the R script in `rscripts/violinPlots.R`
4. Enter output directory name. NOTE: the output directory is relative to the project directory.

*NOTE: Please have these packages installed in R:* Please see full [sessionInfo](https://github.com/urmi-21/COVID-DEA/blob/master/rscripts/sessionInfo.txt)
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


## List of other analysis

Here is a list of all the analysis performed along with code to reproduce results.

* Perform dip-test for all genes: `rscripts/diptest_allgenes.R`
* Perform GSEA for DE lists from pooled GTEx and pooled TCGA samples and generate cnet/ridge plots: `rscripts/gsea_gtex_tcga.R`
* Perform GSEA tissue/tumorwise and generate cnet/ridge plots: `rscripts/gsea_tissuewise.R`
* Perform limma DE analysis in tissue/tumor-wise manner: `rscripts/limma_tissuewise.R`
* Perform limma DE analysis for BRCA samples including molecular-subtyes in the model (uses TCGABiolinks): `rscripts/limma_tissuewise.R`
* Compare Mann-Whitney and limma DE results; performs GSEA and compare enriched terms: `rscripts/limma_MW_compare_gsea.R`


### References
1. Singh, Urminder, Manhoi Hur, Karin Dorman, and Eve Syrkin Wurtele. "MetaOmGraph: a workbench for interactive exploratory data analysis of large expression datasets." Nucleic acids research 48, no. 4 (2020): e23-e23.

## R Session Info

Please see [sessionInfo](https://github.com/urmi-21/COVID-DEA/blob/master/rscripts/sessionInfo.txt)
