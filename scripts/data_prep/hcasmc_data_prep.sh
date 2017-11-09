#!/bin/bash

# TODO: Reset relative paths to appropriate absolute paths, now that I've moved this to the scripts folder.

# Format and index CAD data

header="Markername\tsnptestid\tchr\tsnp_pos\ta1\ta2\teffect_allele_freq\tlogOR\tse_gc\tpvalue\tn_samples\texome\tinfo_ukbb"

# Copy file from Bosh, sort it, and relabel columns to make it compatible with coloc pipeline
cat <(echo -e $header) \
	<(zcat /srv/persistent/bliu2/HCASMC_eQTL/data/gwas/ukbb/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz | awk '{if (($3 != "NA") && ($4 != "NA")) print $0}' | tail -n +2 | sort -k3,3 -k4,4n) > UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt

bgzip -f UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt
tabix -s 3 -b 4 -e 4 -S 1  UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz

################################################################################################

# Format and index CARDIoGRAM data

# Het p-value = heterogeneity p-value, not what we're interested in
header="markername\tchr\tsnp_pos\ta1\ta2\teffect_allele_freq\tmedian_info\tmodel\tbetax\tse_dgc\tpvalue\thet_pvalue\tn_studies"

cat <(echo -e $header) \
	        <(cat /srv/persistent/bliu2/HCASMC_eQTL/data/gwas/CARDIoGRAMplusC4D/cad.add.160614.website.txt | awk '{if (($3 != "NA") && ($4 != "NA")) print $0}' | tail -n +2 | sort -k2,2 -k3,3n) > CARDIoGRAM_cad.add.160614.website.txt

bgzip -f CARDIoGRAM_cad.add.160614.website.txt
tabix -s 2 -b 3 -e 3 -S 1 CARDIoGRAM_cad.add.160614.website.txt.gz

################################################################################################

# Compile all eQTL data into a single file, sort it, compress it, and index it

#for f in `find /srv/persistent/bliu2/HCASMC_eQTL/processed_data/rasqual/output_pval/ -name *pval.txt`;
#do
#	tail -n +2 $f >> /users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/hcasmc_eqtls_presorting.txt
#done

find /srv/persistent/bliu2/HCASMC_eQTL/processed_data/rasqual/output_pval/ -name *pval.txt | head -n 1 | xargs cat | head -n 1 | sed s/pos/snp_pos/g | sed s/fid/gene/g  > /users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/hcasmc_eqtls.txt
sort -k3,3 -k4,4n /users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/hcasmc_eqtls_presorting.txt | sed s/chr//g >> /users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/hcasmc_eqtls.txt 
bgzip -f /users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/hcasmc_eqtls.txt
tabix -s 3 -b 4 -e 4 -S 1 /users/mgloud/projects/brain_gwas/data/eqtls/hcasmc/hcasmc_eqtls.txt.gz

