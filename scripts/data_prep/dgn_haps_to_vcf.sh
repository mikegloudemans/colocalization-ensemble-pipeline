#!/bin/bash
# Author: Mike Gloudemans
# Date: 9/22/2017

for i in {1..22}; do

	cp /mnt/lab_data/montgomery/shared/datasets/dgn/genotypes_from_kim/phasing/chr$i.haps.gz .
	cp /mnt/lab_data/montgomery/shared/datasets/dgn/genotypes_from_kim/phasing/chr$i.sample.gz .
	gunzip chr$i.haps.gz
	gunzip chr$i.sample.gz

	/srv/persistent/bliu2/tools/shapeit.v2.r790.RHELS_5.4.dynamic/shapeit -convert --input-haps chr$i \
		--output-vcf chr$i.haps.phased.vcf

	rm -f chr$i.haps
	rm -f chr$i.sample
done
