import subprocess
import sys
#from scipy import stats
#from shutil import copyfile

in_file = sys.argv[1]
N_trait1 = sys.argv[2]
s_trait1 = sys.argv[3]
type_trait1 = sys.argv[4]
N_trait2 = sys.argv[5]
s_trait2 = sys.argv[6]
type_trait2 = sys.argv[7]
out_file = sys.argv[8]
snakefile = sys.argv[9]

def main():

	# Run COLOC with R helper script and get results
	coloc_status = subprocess.run(f"Rscript {snakefile}/../scripts/methods/coloc/run_coloc.R {in_file} {N_trait1} {s_trait1} {type_trait1} {N_trait2} {s_trait2} {type_trait2}".split(), capture_output=True)
	coloc_status = coloc_status.stdout.decode("utf-8")

	with open(out_file, "a") as a:
		a.write(coloc_status + "\n")


if __name__ == "__main__":
	main()
