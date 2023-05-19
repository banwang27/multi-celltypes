#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G

#Create .coords files from MAF.  Relies on callSNPsFromMAF.py
for f in *.maf
do
  python callSNPsFromMAF.py $f z$f.coords
  echo "Finished processing file: $f"
done

#Concatenate all the files by chromosome.
for i in {X,Y,MT}
do
  cat Compara.H.sap_P.tro_lastz_net_on_H.sap.$i.*.coords > Compara.H.sap_P.tro_lastz_net_on_H.sap.$i.ALL.maf.coords
  echo "Finished processing chromosome $i."
done

for i in {1..22}
do
  cat Compara.H.sap_P.tro_lastz_net_on_H.sap.$i.*.coords > Compara.H.sap_P.tro_lastz_net_on_H.sap.$i.ALL.maf.coords
  echo "Finished processing chromosome $i."
don

#Grab only the ones with differences
for f in *ALL*.coords
do
  awk '{ if ($4 != $5) {print}}' $f > $f.differ.txt
  echo "Finished filtering file: $f"
done

#Merge all the SNPs into one big file
cat *differ.txt > all_snps.txt
echo "Merged SNPs"

#Split the file for memory concerns
split all_snps.txt -l 5000000
echo "Split SNPs"

#Reformat for downstream filtering.  It is not really a bed format.  
for f in xa*
do
  python make_bed.py $f
  echo "Finished making bedfile for: $f"
done
echo "Made bed files"

#Merge again
cat xa*chimp* > chimp_referenced_chp_hum_snps.bed
cat xa*human* > human_referenced_chp_hum_snps.bed
echo "Merged bed files"

#Liftover if desired
#sed 's/^/chr/' chimp_referenced_chp_hum_snps.bed > chimp_referenced_chp_hum_snps.lift
#echo "Lifting over"
#./liftOver minMatch=0.9 chimp_referenced_chp_hum_snps.lift panTro4ToPanTro5.over.chain.gz chimp_referenced_chp_hum_snps_lifted.bed chimp_referenced_chp_hum_snps_lifted.err
#sed 's/...//' chimp_referenced_chp_hum_snps_lifted.bed > chimp_referenced_chp_hum_snps_over.bed
#echo "Done with liftover"

#Switch to refseq chromosome names.  Need to change these in make_bed.py and here for other species.
#Need to add _over to the end of chimp referenced if you uncomment out the liftover part.
python switchChrom.py human_referenced_chp_hum_snps.bed human_referenced_chp_hum_snps_NC.bed humanChroms.txt refseq
python switchChrom.py chimp_referenced_chp_hum_snps.bed chimp_referenced_chp_hum_snps_NC.bed chimpChroms.txt refseq
echo "Chromsosome names swapped, ready for filtering"