# Synthetic target, with no experimental gene, sgRNA or sample identifiers.
synthetic_reference <- read_fasta_base("example/synthetic_reference.fa")[[1L]]$seq
synthetic_guide <- "ATGCTACTCGTACGATTCGA"
stopifnot(length(find_guide_hits(synthetic_reference, synthetic_guide, "NGG")$hits) == 1L)
