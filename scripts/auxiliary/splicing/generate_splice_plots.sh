# Index and prepare GTF
grep exon /users/mgloud/projects/lines/data/genes.hg19.gtf | sort -k1,1 -k4,4g | bgzip > /users/mgloud/projects/brain_gwas/data/spliceplot/hg19.exons.gtf.gz
tabix -s 1 -b 4 -e 5 /users/mgloud/projects/brain_gwas/data/spliceplot/hg19.exons.gtf.gz
# ^ Important to include -b and -e because we want anything overlapping the full exon, not just part of it.

# Run SplicePlot
# Note: With this tool, there must be one position in common that is shared between all the sites.
#python /users/mgloud/software/SplicePlot/initialize_data.py chr12:56115778 chr12:56115278-56117670,chr12:56115278-56115473 --vcf /srv/persistent/bliu2/rpe/data/genotype/asvcf/glucose_nodup/rpe.imputed.chr12.all_filters.vcf.new.gz --gtf /users/mgloud/projects/brain_gwas/data/spliceplot/hg19.exons.gtf.gz --mf /users/mgloud/projects/brain_gwas/data/spliceplot/rpe_glucose_map_file.txt


