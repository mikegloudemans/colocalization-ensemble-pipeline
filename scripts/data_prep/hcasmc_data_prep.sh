#!/bin/bash
# Author: Mike Gloudemans
# Updated: 11/13/2017

# Format and index CAD data
hcasmc_gwas_dir=/users/mgloud/projects/brain_gwas/data/gwas/hcasmc

# Note that order of alt and ref is swapped in this example
header="Markername\tsnptestid\tchr\tsnp_pos\talt\tref\talt_freq\tlog_or\tse\tpvalue\tn_samples\texome\tinfo_ukbb"

# Copy file from Bosh, sort it, and relabel columns to make it compatible with coloc pipeline
cat <(echo -e $header) \
	<(zcat /srv/persistent/bliu2/HCASMC_eQTL/data/gwas/ukbb/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz | awk '{if (($3 != "NA") && ($4 != "NA")) print $0}' | tail -n +2 | awk 'BEGIN {OFS="\t"}; {$4 = sprintf("%d", $4); print}' | sort -k3,3 -k4,4n) > $hcasmc_gwas_dir/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt

bgzip -f $hcasmc_gwas_dir/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt
tabix -s 3 -b 4 -e 4 -S 1  $hcasmc_gwas_dir/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz

################################################################################################

# Format and index CARDIoGRAM data

# Het p-value = heterogeneity p-value, not what we're interested in
# Again, alt and ref order is swapped in this file.
header="markername\tchr\tsnp_pos\talt\tref\talt_freq\tmedian_info\tmodel\tbeta\tse\tpvalue\thet_pvalue\tn_studies"

cat <(echo -e $header) \
	<(cat /srv/persistent/bliu2/HCASMC_eQTL/data/gwas/CARDIoGRAMplusC4D/cad.add.160614.website.txt | awk '{if (($3 != "NA") && ($4 != "NA")) print $0}' | tail -n +2 | awk 'BEGIN {OFS="\t"}; {$3 = sprintf("%d", $3); print}' | sort -k2,2 -k3,3n) > $hcasmc_gwas_dir/CARDIoGRAM_cad.add.160614.website.txt

bgzip -f $hcasmc_gwas_dir/CARDIoGRAM_cad.add.160614.website.txt
tabix -s 2 -b 3 -e 3 -S 1 $hcasmc_gwas_dir/CARDIoGRAM_cad.add.160614.website.txt.gz

################################################################################################

# Compile all eQTL data into a single file, sort it, compress it, and index it
rm -f /users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/hcasmc_eqtls_presorting.txt
for f in `find /srv/persistent/bliu2/HCASMC_eQTL/processed_data/rasqual/output_pval/ -name *pval.txt`;
do
	tail -n +2 $f >> /users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/hcasmc_eqtls_presorting.txt
done

find /srv/persistent/bliu2/HCASMC_eQTL/processed_data/rasqual/output_pval/ -name *pval.txt | head -n 1 | xargs cat | head -n 1 | sed s/pos/snp_pos/g | sed s/fid/gene/g  > /users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/hcasmc_eqtls.txt
# Remove SNPs where ref and alt allele are the same. I don't know how this is even happening, but we don't want to mess with these.
sort -k3,3 -k4,4g /users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/hcasmc_eqtls_presorting.txt | sed s/chr//g | awk '{if ($5 != $6) print $0}' >> /users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/hcasmc_eqtls.txt 
bgzip -f /users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/hcasmc_eqtls.txt
tabix -s 3 -b 4 -e 4 -S 1 /users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/hcasmc_eqtls.txt.gz

################################################################################################

# Prep splicing QTLs

sqtl_dir=/users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/sqtls
header="rsid\tchr\tsnp_pos\tref\talt\tfeature\tcluster\tfeature_distance\tpvalue\tbeta\tse"
echo -e $header > $sqtl_dir/hcasmc_sqtls.txt

join -1 3 -2 3 \
	        <(zcat /srv/persistent/bliu2/HCASMC_eQTL/data/joint3/asvcf/phased_and_imputed.chr*.rename.dr2.hwe.indellt51.rnasample.hg19.vcf.new.gz | grep -v "\#" | cut -f1,2,3,4,5 | sort -k3,3) \
		        <(zcat /srv/persistent/bliu2/HCASMC_eQTL/processed_data/sqtl/fastQTL/nominal/all.nominal.txt.gz | awk '{if ($2 != ".") print $0}' | sed s/:/\\t/g | awk '{print $1 "_" $2 "_" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' | sort -k3,3) | sed s/chr//g | grep -v nan | sort -k2,2 -k3,3n | sed s/\ /\\t/g >> $sqtl_dir/hcasmc_sqtls.txt

bgzip -f $sqtl_dir/hcasmc_sqtls.txt
tabix -S 1 -s 2 -b 3 -e 3 $sqtl_dir/hcasmc_sqtls.txt.gz

