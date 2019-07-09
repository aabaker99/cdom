library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(argparse)

# TODO actually cant use this package
# try to look up NF1:
# rv = genes(edb, filter = ~ gene_name=="NF1")  
# in return value you see  31094927-31382116 for chr17
# however, you cant just query it back
#> gnm = GRanges("17:31094927-31382115")                                                                                
#> gnm_prt <- genomeToTranscript(gnm, edb)                  
#Warning message:                                           
#1 genomic region(s) could not be mapped to a transcript; hint: see ?seqlevelsStyle if you used UCSC chromosome names   
#> gnm_prt                                                                                                              
#IRangesList of length 1                                    
#[[1]]                        
#IRanges object with 1 range and 7 metadata columns:                                                                    
#          start       end     width |       tx_id     exon_id exon_rank                                                
#      <integer> <integer> <integer> | <character> <character> <integer>
#  [1]        -1        -1         1 |        <NA>        <NA>      <NA>                                                
#      seq_start   seq_end    seq_name  seq_strand                                                                      
#      <integer> <integer> <character> <character>                                                                      
#  [1]  31094927  31382115          17           *     

tmp = function() {
	edb = EnsDb.Hsapiens.v86 # GRCh38
	# mcols(SYP) contains metadata associated with the columns argument
	SYP = proteins(edb, filter = ~ genename == "SYP", columns=c("protein_id", "tx_id", listColumns(edb, "protein_domain")), return.type = "AAStringSet")
	SYP = SYP[names(SYP) == "ENSP00000263233"]
	SYP_rng = IRanges(start = mcols(SYP)$prot_dom_start, end=mcols(SYP)$prot_dom_end)

	# map protein domains to genomic coordinates
	# return a list of GRanges (GenomicRanges)
	SYP_gnm = proteinToGenome(SYP_rng, edb, id = "protein_id")

	# Get chromosome, start, end for each protein domain in genomic coordinates
	levels(runValue(seqnames(SYP_gnm[[1]][1])))
	start(ranges(SYP_gnm[[1]][1]))
	end(ranges(SYP_gnm[[1]][1]))
}

main = function() {
  
}

main()
