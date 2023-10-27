#!/bin/bash
#
#$ -q b.q
#$ -wd /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/logs
#$ -pe smp 1
#$ -t 1-4
#$ -j y

# chmod 755 /home/rparrish/YangFSSdata/SR_TWAS/sim/type1_error/submit_split_pred.sh

# qsub /home/rparrish/YangFSSdata/SR_TWAS/sim/type1_error/submit_split_pred.sh

i=${SGE_TASK_ID:-1}

dir=/home/rparrish/YangFSSdata/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/pred
out_dir=${dir}/split_files

datasets=(Naive_EN SR_EN GTEx_EN ROSMAP)
filenames=(Naive_pred_EN SR_pred_EN GTEx_EN_pred ROSMAP_pred)

# split into files by sample
awk -v dataset=${datasets[$i-1]} -v filename=${filenames[$i-1]} -v out_dir=${out_dir} 'BEGIN { FS = OFS = "\t" }{
	# print header
	if ($0 ~ /^CHROM/){
		split($0, h, "\t");
		# print files list 
		for (j = 6; j <= NF; j++) {
			print out_dir"/"dataset"_"h[j]".txt" > out_dir"/"dataset"_paths.txt"
		}
		next
	}
	# fields
	targetid=$4;
	for (j = 6; j <= NF; j++) {
		print targetid"\t"$j > out_dir"/"dataset"_"h[j]".txt"
	}
}' ${dir}/${filenames[$i-1]}.txt

# sort output
for file in `cat ${out_dir}/${datasets[$i-1]}_paths.txt`; do 
	sort -n -o ${file} ${file} 
done




