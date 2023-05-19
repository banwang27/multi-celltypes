#!/bin/bash

for f in *.maf
do
  python callSNPsFromMAF.py $f $f.coords
  echo "Finished processing file: $f"
done
