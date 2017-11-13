# Prep splicing QTLs

sqtl_dir=/users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/sqtls
header="rsid\tchr\tsnp_pos\tref\talt\tfeature\tcluster\tfeature_distance\tpvalue\tbeta\tse"
echo -e $header > $sqtl_dir/hcasmc_sqtls.txt

join -1 3 -2 3 \
	<(zcat /srv/persistent/bliu2/HCASMC_eQTL/data/joint3/asvcf/phased_and_imputed.chr*.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz | grep -v "\#" | cut -f1,2,3,4,5 | sort -k3,3) \
	<(zcat /srv/persistent/bliu2/HCASMC_eQTL/processed_data/sqtl/fastQTL/nominal/all.nominal.txt.gz | awk '{if ($2 != ".") print $0}' | sed s/:/\\t/g | awk '{print $1 "_" $2 "_" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' | sort -k3,3) | sed s/chr//g | grep -v nan | sort -k2,2 -k3,3n | sed s/\ /\\t/g >> $sqtl_dir/hcasmc_sqtls.txt

bgzip -f $sqtl_dir/hcasmc_sqtls.txt
tabix -S 1 -s 2 -b 3 -e 3 $sqtl_dir/hcasmc_sqtls.txt.gz
