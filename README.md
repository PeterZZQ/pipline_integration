## Preprocessing count matrix of scATAC-Seq and scRNA-Seq for integration

### Description

This is the preprocessing pipeline used in `scDART` and `scCFRM`, and these applications show that the pipeline works well in preserving cellular indentity in the datasets while significantly reduce the dimensionality.



### Before using

Before using the pipeline, make sure:

* Download the **reference genome**, link can be found in `./reference_genome/link`, or you can download it with dropbox [link](https://www.dropbox.com/sh/ss03bj9rr3ea0n9/AADfR5l2uueObCTqMVVvUrSha?dl=0).

* Save all count matrices into the `./raw/` folder, the data should include

  * `barcodes_rna.txt`: barcodes in scRNA-Seq data, one line correspond to one barcode, no header.
  * `genes.txt`: genes in scRNA-Seq data, one line correspond to one gene, no header.
  * `G.mtx`: the scRNA-Seq count matrix, of the shape **(ngenes, ncells)**, matrix market format.
  * `meta_rna.csv`: the meta info for cells in scRNA-Seq (e.g. ground truth cluster idenity, etc). A data frame, should have index matching `barcodes_rna.txt`. The cluster annotation is stored in the column `celltype`. (Optional, note that validation cannot be performed without the file).

  

  * `barcodes_atac.txt`: barcodes in scATAC-Seq data, one line correspond to one barcode, no header.
  * `regions.txt`: regions in scATAC-Seq data, one line correspond to one region, no header. The region should follow `chrX_start_end` format.
  * `R.mtx`: the scRNA-Seq count matrix, of the shape **(nregions, ncells)**, matrix market format.
  * `meta_atac.csv`: the meta info for cells in scATAC-Seq (e.g. ground truth cluster idenity, etc). A data frame, should have index matching `barcodes_atac.txt`. The cluster annotation is stored in the column `celltype`. (Optional, note that validation cannot be performed without the file).



### Usage

Check `preprocessing.ipynb`.



### Dependency

**Python:**

* `anndata`, ver. 0.7.5 
* `scanpy`, ver. 1.6.0

**R**:  

* [`Seurat` ver 3.2.3](https://satijalab.org/seurat/articles/install.html) (or alternatively `GenomicRanges`, `BiocGenerics`, `rtracklayer`, `GenomeInfoDb`, `Signac`, `GenomicRanges`

