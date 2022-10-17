# Cell type annotation

To predict the cell labels I decided to use singleR (v 1.8.1) and as reference the dataset from [Triana *et al* 2021](https://www.nature.com/articles/s41590-021-01059-0), precisely the healthy transcriptome containing  13165 annotated cells using Ab and proteomics.

The queries datasets are:
- Collection of all MM datasets after QC, raw counts
- Collection of HEalthy datasets after QC, raw counts

After annotating the cells I converted the files in anndata format using sceasy (v 0.0.5) available from the container [batch_correct](https://singularityhub.github.io/singularityhub-archive/containers/Sarah145-batch_correct-latest/) build by [Sarah Ennis](https://github.com/Sarah145).

The conversion was done to:
- To run the batch correction with scVI available [here](../Merge_batch-correct/).
- Score the set of markers to identify those NK that show exhausted phenotypes (set of markers is [here](https://github.com/VergaJU/ScoreMarkers/blob/main/data/nk_exhaustion.csv) retrieved from the paper: [Foroutan *et al* 2021](https://doi.org/10.1158/2326-6066.CIR-21-0137))

The report on batch correction and cell type annotation is available [here](../Reports/Batch_annotation.ipynb) and [here](../Reports/cell-proportions.ipynb)

The report on NK cells proportion is [here](../Reports/nk_cells.ipynb)



