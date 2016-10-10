file=$1
name=$(echo $file)
parts=${name//./}
size=${#parts[@]}
stem0=${parts[size-1]}
stem=${stem0//fasta/}
mafft --localpair --maxiterate 1000 $file > aligned/aln_$stem.fasta
#echo "$file is done"

