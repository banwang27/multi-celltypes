#!/bin/bash
#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4G

for f in *.maf
do
  python callINDELsFromMAF.py $f $f.indels
  echo "Finished processing file: $f"
done
