#!/bin/bash
for file in *.fasta
do
./convert_to_xml.sh $file beastgen.3part.HKYG.UCLN.CS.template hky_ucln
./convert_to_xml.sh $file beastgen.3part.HKYG.STRICT.CS.template hky_strict
done

