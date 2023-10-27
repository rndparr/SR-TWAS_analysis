#!/bin/bash

# load tabix
module load tabix/0.2.6

# without window:
# 19:1040101-1065572

tabix_region=19:40101-2065572

# GENOTYPE DATA

# # ROSMAP
# ROSMAP_in=/mnt/YangFSS/data2/rparrish/ROSMAP_WGS/genotype/CHR19_ROSMAP_WGS_b38_lifted.txt.gz
# ROSMAP_out=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/genotype/ROSMAP_ABCA7_raw.vcf

# ## header row to out file
# zgrep -m1 -E '^#CHROM' ${ROSMAP_in} > ${ROSMAP_out}


# # bgzip -d ${ROSMAP_in}
# # bgzip /mnt/YangFSS/data2/rparrish/ROSMAP_WGS/genotype/CHR19_ROSMAP_WGS_b38_lifted.txt
# # tabix -fp vcf ${ROSMAP_in}

# ## data lines to out file
# tabix ${ROSMAP_in} ${tabix_region} | \
# awk -c -v fmt="GT" 'BEGIN { FS = OFS = "\t" }{
# 	# split row
# 	N_cols = split($0, row, "\t"); 

# 	# save alt alleles to an array; save number of alt alleles
# 	N_alt = split(row[5], ALT, ",");

# 	# split FORMAT field, get position of DS/GT from FORMAT field,
# 	# save value for use with data columns
# 	split(row[9], FORMAT, ":");
# 	for (j in FORMAT) {
# 		if (FORMAT[j] == fmt) {
# 			fmt_pos = j
# 		}
# 	}
# 	row[9] = fmt

# 	# iterate through the sampleid columns, get fmt
# 	for (i = 10; i <= NF; i++) {
# 		# save only DS/GT value of data fields
# 		split(row[i], DATA, ":");
# 		row[i]=DATA[fmt_pos];
# 	}

# 	# iterate through alt alleles
# 	for (j = 1; j <= N_alt; j++){
# 		# columns
# 		printf row[1]"\t"row[2]"\t.\t"row[4]"\t"ALT[j]"\t.\t.\t.\t"row[9]"\t";

# 		# iterate through the sampleid values, make biallelic
# 		for (i = 10; i <= NF; i++) {
# 			# prevent gsub from changing row[i]
# 			val=row[i];

# 			# replace all references to alt alleles != j
# 			for (k = 1; k <= N_alt; k++) {
# 				if (k != j) {
# 					gsub(k, ".", val)
# 				}
# 			}
# 			# set jth alt allele designation to 1 (has to be last)
# 			gsub(j, "1", val);

# 			# add the data value to the line
# 			printf val (i==N_cols?"\n":"\t");
# 		}
# 	}
# }' \
# >> ${ROSMAP_out}

# ## bgzip & tabix
# bgzip -f ${ROSMAP_out}
# tabix -f -p vcf ${ROSMAP_out}.gz

ROSMAP_in=/mnt/YangFSS/data2/rparrish/ROSMAP_WGS/genotype/CHR19_ROSMAP_WGS_b38_lifted_TIGAR.txt.gz
ROSMAP_out=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/genotype/ROSMAP_ABCA7_raw.vcf

zgrep -m1 -E '^#CHROM' ${ROSMAP_in} > ${ROSMAP_out}
tabix ${ROSMAP_in} ${tabix_region} >> ${ROSMAP_out}
bgzip -f ${ROSMAP_out}

tabix -f -p vcf ${ROSMAP_out}.gz


# GTEx 
GTEx_in=/mnt/YangFSS/data/rparrish/GTEx_V8/GenotypeFiles/CHR19_GTEx_WGS.vcf.gz
GTEx_out=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/genotype/GTEx_ABCA7_raw.vcf

zgrep -m 1 -E '^#CHROM' ${GTEx_in} > ${GTEx_out}
tabix ${GTEx_in} ${tabix_region} >> ${GTEx_out}
bgzip -f ${GTEx_out}

tabix -f -p vcf ${GTEx_out}.gz


# unload
module unload tabix/0.2.6



####

# Headers:
dataset=GTEx
dataset=ROSMAP

geno=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/genotype/${dataset}_ABCA7_raw.vcf.gz
sampleid_out=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/sampleid/${dataset}_genotype_sampleid.txt

zgrep -m 1 -E '^#CHROM' ${geno} | \
	awk 'BEGIN { FS = OFS = "\t" } {for (i = 1; i <= NF; i++) { if (i > 9) { printf $i"\n" }}}' \
	> ${sampleid_out}

wc -l ${sampleid_out}
# GTEx: 715
# ROSMAP: 465


####

# only european samples for gtex
geno=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/genotype/GTEx_ABCA7_raw.dosage.gz
sampleid_out=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/sampleid/GTEx_genotype_sampleid.txt
sampleid_out=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/sampleid/GTEx_sampleid.txt

zgrep -m 1 -E '^CHROM' ${geno} | \
	awk 'BEGIN { FS = OFS = "\t" } {for (i = 1; i <= NF; i++) { if (i > 5) { printf $i"\n" }}}' \
	> ${sampleid_out}


####

key=/mnt/YangFSS/data2/AMP-AD/ROSMAP/WGS_JointCall/AMP-AD_rosmap_WGS_id_key.csv

sed 's/\r/\n/g' $key | awk -F"," '{print $2}' | wc -l


awk -F"^M" '{print $1}' $key | less -S







# zgrep -m 1 -E '^CHROM' ${geno} | less -S

# zgrep -m 1 -E '^CHROM' ${geno} | awk 'BEGIN { FS = OFS = "\t" } {for (i = 1; i <= NF; i++) { if (i > 5) { printf $i"\n" }}}' | less -S

# zgrep -m 1 -E '^CHROM' ${geno} | awk 'BEGIN { FS = OFS = "\t" } {for (i = 1; i <= NF; i++) { if (i > 9) { printf $i"\n" }}}' | less -S

# zgrep -m 1 -E '^#CHROM' ${geno} | awk 'BEGIN { FS = OFS = "\t" } {for (i = 1; i <= NF; i++) { if (i > 9) { printf $i"\n" }}}' > 

# zgrep -m1 -E '^#CHROM' ${in} | awk 'BEGIN { FS = OFS = "\t" } {for (i = 1; i <= NF; i++) { if (i > 9) { printf $i"\n" }}}' \
#   > /mnt/YangFSS/data2/rparrish/ROSMAP_WGS/sampleid/ROSMAP_WGS_sampleid_genotype.txt

######

# wc -l /mnt/YangFSS/data/rparrish/thesis/sim/sampleid/ROSMAP_gt_sampleid.txt

