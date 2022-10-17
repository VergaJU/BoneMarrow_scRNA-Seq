# Cell type annotation

To predict the cell labels I decided to use singleR (v 1.8.1) and as reference the dataset from [Triana *et al* 2021](https://www.nature.com/articles/s41590-021-01059-0), precisely the healthy transcriptome containing  13165 annotated cells using Ab and proteomics.

The queries datasets are:
- Collection of all MM datasets after QC, raw counts
- Collection of HEalthy datasets after QC, raw counts

After annotating the cells I converted the files in anndata format using sceasy (v 0.0.5) available from the container [batch_correct](https://singularityhub.github.io/singularityhub-archive/containers/Sarah145-batch_correct-latest/) build by [Sarah Ennis](https://github.com/Sarah145).

To evaluate the annotation quality and proportions I created some jupyter notebooks available [here](./report/)