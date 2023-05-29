folder_name="FABP7"
ml R/4.2.0
cd /scratch/users/banwang3/references/outputs/${folder_name}

Rscript /scratch/users/banwang3/references/scripts/summariseresult.R --generator ${folder_name}
Rscript /scratch/users/banwang3/references/scripts/summariseresult_chp.R --generator ${folder_name}
Rscript /scratch/users/banwang3/references/scripts/sei_pred_boxplot.R --generator ${folder_name}
Rscript /scratch/users/banwang3/references/scripts/sei_pred_boxplot_chp.R --generator ${folder_name}
