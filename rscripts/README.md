## Scripts used for data analysis
Here is a list of all the analysis performed along with code to reproduce results.

* Perform dip-test for all genes: `diptest_allgenes.R`
* Perform GSEA for DE lists from pooled GTEx and pooled TCGA samples and generate cnet/ridge plots: `gsea_gtex_tcga.R`
* Perform GSEA tissue/tumorwise and generate cnet/ridge plots: `gsea_tissuewise.R`
* Perform limma DE analysis in tissue/tumor-wise manner: `limma_tissuewise.R`
* Perform limma DE analysis for BRCA samples including molecular-subtyes in the model (uses TCGABiolinks): `limma_tissuewise.R`
* Compare Mann-Whitney and limma DE results; performs GSEA and compare enriched terms: `limma_MW_compare_gsea.R`
* Plot violin plot via MetaOmGraph: `violinPlots.R`