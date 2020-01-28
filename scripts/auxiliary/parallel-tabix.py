# Sort, zip, tabix edQTL files
# Now in parallel!

from multiprocessing import Pool
import subprocess
import glob
import sys
import traceback

max_cores = 15

def main():
    files = glob.glob("/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sqtl/GTEx_Analysis_v8_sQTL_all_associations/*")
    pool = Pool(max_cores)
    for f in files:
        tissue = f.split("/")[-1].split(".v8.")[0]
        pool.apply_async(tabix_tissue_wrapper, args=(tissue,))
    pool.close()
    pool.join()

def tabix_tissue_wrapper(tissue):
    try:
        tabix_tissue(tissue)
    except Exception:
        traceback.print_exc(file=sys.stdout)


def tabix_tissue(tissue):
    with open("/srv/scratch/mgloud/brain_gwas/data/sqtls/gtex_v8/{0}.sQTLs.txt".format(tissue), "w") as w:
        w.write("feature\tfeature2\tchr\tsnp_pos\tref\talt\tbuild\ttss_distance\tma_samples\tma_count\tmaf\tpvalue\tbeta\tse\n") 
            
    subprocess.check_call("zcat /mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/sqtl/GTEx_Analysis_v8_sQTL_all_associations/{0}.v8.sqtl_allpairs.txt.gz | sed s/_/\\\\t/g | sed s/chr//g | tail -n +2 | sort -k3,3 -k4,4n >> /srv/scratch/mgloud/brain_gwas/data/sqtls/gtex_v8/{0}.sQTLs.txt".format(tissue), shell=True)
    subprocess.check_call("bgzip -f /srv/scratch/mgloud/brain_gwas/data/sqtls/gtex_v8/{0}.sQTLs.txt".format(tissue), shell=True)
    subprocess.check_call("tabix -f -S 1 -s 3 -b 4 -e 4 /srv/scratch/mgloud/brain_gwas/data/sqtls/gtex_v8/{0}.sQTLs.txt.gz".format(tissue), shell=True)

if __name__ == "__main__":
    main()
