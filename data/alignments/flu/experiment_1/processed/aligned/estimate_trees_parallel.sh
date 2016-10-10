 ls *.fasta | parallel --nice 10 --max-procs 14  Rscript phyml_R.r 
