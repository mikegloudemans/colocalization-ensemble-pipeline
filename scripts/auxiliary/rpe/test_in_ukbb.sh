echo "Eye problems/disorders: Macular degeneration"
cat <(zcat /mnt/lab_data/montgomery/shared/datasets/ukbb/gwas/6148_5.assoc.sorted.tsv.gz | head -n 1) <(zcat /mnt/lab_data/montgomery/shared/datasets/ukbb/gwas/6148_5.assoc.sorted.tsv.gz | egrep "rs3138142|rs1532277|rs12899503|rs45563938|rs35064805|rs7229687|rs28661406" | column -t)
echo ""

echo "Non-cancer illness code, self-reported: macular degeneration"
cat <(zcat /mnt/lab_data/montgomery/shared/datasets/ukbb/gwas/6148_5.assoc.sorted.tsv.gz | head -n 1) <(zcat /mnt/lab_data/montgomery/shared/datasets/ukbb/gwas/20002_1528.assoc.sorted.tsv.gz | egrep "rs3138142|rs1532277|rs12899503|rs45563938|rs35064805|rs7229687|rs28661406" | column -t)
echo ""

echo "Wears glasses or contact lenses"
cat <(zcat /mnt/lab_data/montgomery/shared/datasets/ukbb/gwas/6148_5.assoc.sorted.tsv.gz | head -n 1) <(zcat /mnt/lab_data/montgomery/shared/datasets/ukbb/gwas/2207.assoc.sorted.tsv.gz | egrep "rs3138142|rs1532277|rs12899503|rs45563938|rs35064805|rs7229687|rs28661406" | column -t)
echo ""

echo "Reason for glasses/contact lenses: For short-sightedness, i.e. only or mainly for distance viewing such as driving, cinema etc (called 'myopia')"
cat <(zcat /mnt/lab_data/montgomery/shared/datasets/ukbb/gwas/6148_5.assoc.sorted.tsv.gz | head -n 1) <(zcat /mnt/lab_data/montgomery/shared/datasets/ukbb/gwas/6147_1.assoc.sorted.tsv.gz | egrep "rs3138142|rs1532277|rs12899503|rs45563938|rs35064805|rs7229687|rs28661406" | column -t)
echo ""

echo "Reason for glasses/contact lenses: For long-sightedness, i.e. for distance and near, but particularly for near tasks like reading (called 'hypermetropia')"
cat <(zcat /mnt/lab_data/montgomery/shared/datasets/ukbb/gwas/6148_5.assoc.sorted.tsv.gz | head -n 1) <(zcat /mnt/lab_data/montgomery/shared/datasets/ukbb/gwas/6147_2.assoc.sorted.tsv.gz | egrep "rs3138142|rs1532277|rs12899503|rs45563938|rs35064805|rs7229687|rs28661406" | column -t)
echo ""
