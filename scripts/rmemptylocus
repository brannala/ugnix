#!/bin/bash
./gsum -l $1 > temp1 &&
cat temp1 | grep  "no_alleles: 0" | awk '{print $2}' > temp2 &&
grep -v -f temp2 $1    
# echo $LN
