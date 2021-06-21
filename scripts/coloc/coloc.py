import subprocess
import sys
#from scipy import stats
#from shutil import copyfile

in_file = sys.argv[1]
N_source = sys.argv[2]
s_source = sys.argv[3]
type_source = sys.argv[4]
N_lookup = sys.argv[5]
s_lookup = sys.argv[6]
type_lookup = sys.argv[7]
out_file = sys.argv[8]

def main():

	# Run COLOC with R helper script and get results
	coloc_status = subprocess.run(f"Rscript ../scripts/coloc/run_coloc.R {in_file} {N_source} {s_source} {type_source} {N_lookup} {s_lookup} {type_lookup}".split(), capture_output=True)
	print(coloc_status.stderr.decode("utf-8"))
	coloc_status = coloc_status.stdout.decode("utf-8")

	print(coloc_status)
	with open(out_file, "a") as a:
		a.write(coloc_status + "\n")


if __name__ == "__main__":
	main()
