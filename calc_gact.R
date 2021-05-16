rm(list = ls())
gc()

library(Seurat)
library(Matrix)

find_geneact <- function(peak.matrix, annotation.file, seq.levels, upstream = 2000, downstream = 0, verbose = FALSE){
    # peak.matrix is the region by cell matrix, peak.df is all regions
    peak.df <- rownames(x = peak.matrix)
    
    # reformualte peak.df of the form "chromosome", "start", "end"
    peak.df <- do.call(what = rbind, args = strsplit(x = peak.df, split = "_"))
    peak.df <- as.data.frame(x = peak.df)
    colnames(x = peak.df) <- c("chromosome", "start", "end")
    
    # peak.df -> peaks.gr
    peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)
    BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1
    
    # gtf stores the annotation (reference genome)
    gtf <- rtracklayer::import(con = annotation.file)
    gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = "coarse")
    if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
        GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
    }
    # gtf.genes stores the genes 
    gtf.genes <- gtf[gtf$type == "gene"]
    
    # update the regions correspond to each gtf.genes, gtf.body_prom
    gtf.body_prom <- Signac::Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
    
    # assign peaks.gr to nearest gene region
    gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
    # only leave the ones(regions) overlap with the gene regions(distance = 0)
    keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 
                                        0]
    peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
    gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]
    gene.ids$gene_name[is.na(gene.ids$gene_name)] <- gene.ids$gene_id[is.na(gene.ids$gene_name)]
    peak.ids$gene.name <- gene.ids$gene_name
    peak.ids <- as.data.frame(x = peak.ids)
    peak.ids$peak <- rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)]
    # new annotations should include peaks and corresponding gene.name
    annotations <- peak.ids[, c("peak", "gene.name")]
    
    return(annotations)
}



# hyper-parameters
species <- "Mouse"
# upstream region size (base-pair)
upstream <- 2000
# downstream region size (base-pair)
downstream <- 0

path <- "./raw/"

# counts.atac, mm, of the shape regions by barcodes
counts.atac <- readMM(paste0(path, "./RxC1.mtx"))
# regions chrX_start_end
row.names(counts.atac) <- read.table(file = paste0(path, "./regions.txt"), sep = ",", header = FALSE)[[1]]
# barcodes
colnames(counts.atac) <- read.table(file = paste0(path, "./barcodes_atac.txt"), sep = ",", header = FALSE)[[1]]

if(species == "Mouse"){
    A = find_geneact(peak.matrix = counts.atac, annotation.file = "./reference_genome/Mus_musculus.GRCm38.84.gtf", 
                           seq.levels = c(1:19, "X", "Y"), upstream = upstream, downstream = downstream, verbose = TRUE)
} else if(species == "Human"){
    A = find_geneact(peak.matrix = counts.atac, annotation.file = "./reference_genome/Homo_sapiens.GRCh37.82.gtf", 
                           seq.levels = c(1:22, "X", "Y"), upstream = upstream, downstream = downstream, verbose = TRUE)
} else{
    stop("species can only be Human or Mouse")
}

# output gene activity matrix
write.table(A, file = paste0(path, "gact.csv"), sep = ",")
