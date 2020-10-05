## Scripts used for data analysis

Data used in these scripts is available [here](http://metnetweb.gdcb.iastate.edu/MetNet_MetaOmGraph.htm). Note that the GTEx race assignments are missing from the metadata as those are controlled access and can be requested from GTEx. The DE genes used in GSEA analysis are present under supplementary data folder in this repository.

Here is a list of all the analysis performed along with code to reproduce results.

* Perform dip-test for all genes: `diptest_allgenes.R`
* Perform GSEA for DE lists from pooled GTEx and pooled TCGA samples and generate cnet/ridge plots: `gsea_gtex_tcga.R`
* Perform limma DE analysis in tissue/tumor-wise manner: `limma_tissuewise.R`
* Perform limma DE analysis for BRCA samples including molecular-subtyes in the model (uses TCGABiolinks): `limma_tissuewise.R`
* Compare Mann-Whitney and limma DE results; performs GSEA using MW and limma DE genes and compare enriched terms and DE genes between methods: `limma_MW_compare_gsea.R`
* Generate violin plots via MetaOmGraph: `violinPlots.R` **Note:** To reproduce figures in the paper, make sure the samples from PRAD, STAD, LIHC, ESCA, KICH and UCS are filtered out in MetaOmGraph
