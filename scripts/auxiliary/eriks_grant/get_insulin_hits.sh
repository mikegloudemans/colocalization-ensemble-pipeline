cat rsid_list.txt <(sed s//\\n/g /users/mgloud/projects/gwas/data/from_bruna/GWAS_SumStat/WHRadjBMI_GIANT_Europeans_Hits.txt | tail -n +2 | cut -f1)  <(sed s//\\n/g /users/mgloud/projects/gwas/data/from_bruna/GWAS_SumStat/TG_GCLC_Mixed_Hits.txt | tail -n +2 | cut -f2) <(sed s//\\n/g /users/mgloud/projects/gwas/data/from_bruna/GWAS_SumStat/FastInsu_adjBMI_MAGIC_Europeans_Hits.txt | tail -n +2 | cut -f2) <(sed s//\\n/g /users/mgloud/projects/gwas/data/from_bruna/GWAS_SumStat/FastGlu_MAGIC_Europeans_Hits.txt | tail -n +2 | cut -f2) | sort | uniq > insulin_rsids.txt
