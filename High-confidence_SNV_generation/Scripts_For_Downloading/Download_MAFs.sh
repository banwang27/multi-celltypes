#!/bin/bash

#SBATCH --output=MAF_Download.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

#Note, I used UCSC mafs when available and ensembl mafs if there was none from UCSC
#The goal was to ensure the greatest degree of compatibility with the UCSC fasta and gtf files
#which include much much better gene annotations for many NHPs.

#Chimp using UCSC version and hg38/PanTro6
wget -P chimp ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPanTro6/hg38.panTro6.synNet.maf.gz
gunzip chimp/hg38.panTro6.synNet.maf.gz
#Gorilla using hg38/gorgor6
wget -P gorilla ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsGorGor6/hg38.gorGor6.synNet.maf.gz
gunzip gorilla/hg38.panTro6.synNet.maf.gz
#Rhesus using hg38/rheMac10
wget -P macaque ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsRheMac10/hg38.rheMac10.synNet.maf.gz
gunzip macaque/hg38.rheMac10.synNet.maf.gz
#Mouse using hg38/mm39
wget -P mouse ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsMm39/hg38.mm39.synNet.maf.gz
gunzip mouse/hg38.mm39.synNet.maf.gz
#Mouse Lemur using hg38/micMur2
wget -P mouse_lemur ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsMicMur2/hg38.micMur2.synNet.maf.gz
gunzip mouse_lemur/hg38.micMur2.synNet.maf.gz
#Bonobo using hg38/panPan3
wget -P bonobo ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPanPan3/hg38.panPan3.synNet.maf.gz
gunzip bonobo/hg38.panPan3.synNet.maf.gz
#Orangutan using hg38/ponAbe3
wget -P orangutan ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPonAbe3/hg38.ponAbe3.synNet.maf.gz
gunzip orangutan/hg38.ponAbe3.synNet.maf.gz
#Gibbon using hg38/nleu_3
wget -P gibbon ftp://ftp.ebi.ac.uk/ensemblgenomes/pub/metazoa/current/maf/ensembl-compara/pairwise_alignments/hsap_grch38.v.nleu_nleu_3.0.lastz_net.tar.gz
gunzip hsap_grch38.v.nleu_nleu_3.0.lastz_net.tar.gz
tar -xf hsap_grch38.v.nleu_nleu_3.0.lastz_net.tar
#Marmoset using hg38/calJac4
wget -P marmoset ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsCalJac4/hg38.calJac4.synNet.maf.gz
gunzip marmoset/hg38.calJac4.synNet.maf.gz
#Crab eating Macaque using hg38/mfas_6
wget -P macaque_fasciularis ftp://ftp.ebi.ac.uk/ensemblgenomes/pub/metazoa/current/maf/ensembl-compara/pairwise_alignments/hsap_grch38.v.mfas_macaca_fascicularis_6.0.lastz_net.tar.gz
gunzip macaque_fascicularis/hsap_grch38.v.mfas_macaca_fascicularis_6.0.lastz_net.tar.gz
tar -xf macaque_fascicularis/hsap_grch38.v.mfas_macaca_fascicularis_6.0.lastz_net.tar
#Anubis Baboon using hg38/PapAnu4
wget -P anubis_baboon ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPapAnu4/hg38.papAnu4.synNet.maf.gz
gunzip anubis_baboon/hg38.papAnu4.synNet.maf.gz