#!/bin/bash
## Usage: enter a .fasta file name ($1), a template ($2) and a prefix ($3)
name=$(echo $1)
parts=${name//./}
size=${#parts[@]}
stem0=${parts[size-1]}
stem=${stem0//fasta/}
beastgen -D "chain_length=100000000,log_every=10000" -date_order -1 -date_prefix \_ $2 $1 $3_$stem.xml
