# XCI_analysis

XCI_main_script.R contains code to reproduce all main and supplemental figures in the manuscript. 
The code depends on files and data that are compiled in the XCI_data/data directory available in the figshare repository.

Link: []

Code to regenerate those datasets are provided in the other scripts in this directory:

curate_25mers_from_variants.R curates the 7k-ish k-mers that discriminate between reference and alternate alleles on the X chromosome --
	output is the all_kmers_chrX.rds file provided in XCI_paper/data/kmer_data
	functions necessary for this script are found in masterKmers_source.R
scripts to collate & process k-mer counts (the counts themselves are generated from msBWT: https://github.com/holtjma/msbwt)
	takes counts from chrX_alt_counts_sp1.csv; chrX_ref_counts_sp1.csv; all files in sp2_counts directory
	runs JAGS models and perfoms MCMC
	requires ase_summary_source.R and jags_kmers_source.R for necessary functions
	