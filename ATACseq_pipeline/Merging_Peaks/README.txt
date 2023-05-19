The pipeline_merge.sh script included here is if there are multiple conditions, cell types
etc. and you want to create a unified peaklist.  It works very similarly to the other pipeline.sh
and relies on some of the scripts and data in there.  Read pipeline_merge.sh for more details.
It tries to preserve information about which cell types peaks were called in.
The input files are the final output bed files of the Peak_Calling pipeline.sh.

deal_fast.py merges overlapping peaks found by repeated bedtools intersect and
adds back peaks that are found in only one cell type.

merge_all.py and merge_all2.py merge overlapping peaks and change their names

find_partners_all.py finds orthologous peaks between species