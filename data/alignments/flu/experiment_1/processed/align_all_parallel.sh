 ls *.fasta | parallel --nice 10 --max-procs 20  ./align.sh
