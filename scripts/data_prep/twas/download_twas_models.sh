files=(http://gusevlab.org/projects/fusion/weights/GTEx.Adipose_Subcutaneous.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Adipose_Visceral_Omentum.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Adrenal_Gland.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Artery_Aorta.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Artery_Coronary.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Artery_Tibial.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Brain_Amygdala.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Brain_Anterior_cingulate_cortex_BA24.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Brain_Caudate_basal_ganglia.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Brain_Cerebellar_Hemisphere.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Brain_Cerebellum.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Brain_Cortex.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Brain_Frontal_Cortex_BA9.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Brain_Hippocampus.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Brain_Hypothalamus.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Brain_Nucleus_accumbens_basal_ganglia.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Brain_Putamen_basal_ganglia.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Brain_Spinal_cord_cervical_c-1.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Brain_Substantia_nigra.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Breast_Mammary_Tissue.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Cells_EBV-transformed_lymphocytes.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Cells_Transformed_fibroblasts.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Colon_Sigmoid.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Colon_Transverse.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Esophagus_Gastroesophageal_Junction.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Esophagus_Mucosa.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Esophagus_Muscularis.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Heart_Atrial_Appendage.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Heart_Left_Ventricle.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Liver.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Lung.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Minor_Salivary_Gland.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Muscle_Skeletal.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Nerve_Tibial.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Ovary.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Pancreas.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Pituitary.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Prostate.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Skin_Not_Sun_Exposed_Suprapubic.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Skin_Sun_Exposed_Lower_leg.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Small_Intestine_Terminal_Ileum.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Spleen.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Stomach.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Testis.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Thyroid.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Uterus.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Vagina.P01.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Whole_Blood.P01.tar.bz2)

#for file in ${files[@]};
#do
#	#wget $file -P /users/mgloud/projects/brain_gwas/data/twas_models
#	echo $file
#	file=`echo $file | sed s/http:\\/\\/gusevlab.org\\/projects\\/fusion\\/weights//g`
#	echo $file
#	tar -xvf /users/mgloud/projects/brain_gwas/data/twas_models/$file
#done

for file in `ls /users/mgloud/projects/brain_gwas/data/twas_models/`; 
do
	tar -xvf /users/mgloud/projects/brain_gwas/data/twas_models/$file
done
