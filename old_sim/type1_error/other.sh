





dir=/home/rparrish/YangFSSdata/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/pred/
out_dir=${dir}/split_files
datasets=(Naive_EN SR_EN GTEx_EN ROSMAP)
filenames=(Naive_pred_EN SR_pred_EN GTEx_EN_pred ROSMAP_pred)

for i in `seq 0 3`; do 

	awk -v dataset=${datasets[$i]} -v filename=${filenames[$i]} -v out_dir=${out_dir} 'BEGIN { FS = OFS = "\t" }{
		# print header
		if ($0 ~ /^CHROM/){
			split($0, h, "\t");
			# print "TargetID\t"dataset;
			next
		}
		# fields
		targetid=$4;
		for (j = 6; j <= NF; j++) {
			print targetid"\t"$j > out_dir"/"dataset"_"h[j]".txt"
		}
	}' ${dir}/${filenames[$i]}.txt

done



# ############
# ## melt pred job output

# dir=/home/rparrish/YangFSSdata/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/pred/
# datasets=(Naive SR)
# GTEx_model=_EN

# for dataset in ${datasets[@]}; do
# 	awk -v dataset=${dataset} 'BEGIN { FS = OFS = "\t" }{
# 		# print header
# 		if ($0 ~ /^CHROM/){
# 			split($0, a, "\t");
# 			print "TargetID\tsample_id\t"dataset;
# 			next
# 		}
# 		# fields
# 		targetid=$4;
# 		for (i = 6; i <= NF; i++) {
# 			print targetid"\t"a[i]"\t"$i
# 		}
# 	}' ${dir}/${dataset}_pred${GTEx_model}.txt > ${dir}/${dataset}_pred${GTEx_model}_melt.txt

# done

# ############

# dir=/home/rparrish/YangFSSdata/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/pred/

# datasets=(Naive SR)

# head -n1 ${dir}/Naive_pred_EN.txt | \
# awk 'BEGIN { FS = OFS = "\t" } {print "dataset\t"$0'} \
# > ${dir}/all_pred_results.txt

# # combine prediction results
# for dataset in ${datasets[@]}; do

# 	path=${dir}/${dataset}_pred_EN.txt

# 	tail -n+2 ${path} | \
# 		awk -v dataset=${dataset} 'BEGIN { FS = OFS = "\t" } {print dataset"\t"$0}' \
# 		>> ${dir}/all_pred_results.txt

# done

#################
## get info file for base models

dir=/mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/train

echo -e 'CHROM\tGeneStart\tGeneEnd\tTargetID\tGeneName' > ${dir}/GTEx_info_EN.txt

for i in `seq 1 1000000`; do
	echo -e "19\t1040101\t1065572\t"${i}"\tABCA7" >> ${dir}/GTEx_info_EN.txt
done

ln -s ${dir}/GTEx_info_EN.txt ${dir}/ROSMAP_info.txt


#########
# ## get base model header
# echo -e "CHROM\tPOS\tID\tREF\tALT\tTargetID\tES" > /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/train//sims/GTEx_EN_weights_header.txt
# bgzip -f /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/train//sims/GTEx_EN_weights_header.txt


# echo -e "CHROM\tPOS\tID\tREF\tALT\tTargetID\tES" > /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/train//sims/ROSMAP_weights_header.txt
# bgzip -f /mnt/YangFSS/data2/rparrish/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/train//sims/ROSMAP_weights_header.txt

#########

dir=/home/rparrish/YangFSSdata/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/pred/

echo -e "TargetID\tsample_id\tNaive\tSR" > ${dir}/all_pred.txt


awk 'BEGIN { FS = OFS = "\t" } {
	if(FNR==NR){ a[$1"\t"$2]=$3; next }} {
	if($1"\t"$2 in a) {print $1"\t"$2"\t"a[$1"\t"$2]"\t"$3 } }' \
	${dir}/Naive_pred_EN_melt.txt  ${dir}/SR_pred_EN_melt.txt  \
>> ${dir}/all_pred.txt


a[$1"\t"$2]=$3
if($1"\t"$2 in a) {print $1"\t"$2"\t"a[$1"\t"$2]"\t"$3 }


${dir}/Naive_pred_EN_melt.txt 
${dir}/all_pred.txt


TargetID        sample_id       Naive
TargetID        sample_id       SR



###### 


dir=/home/rparrish/YangFSSdata/SR_TWAS/sim/type1_error/overlap_0.1_0.1_0.5_0.1_0.1_508/pred/
datasets=(Naive SR)
GTEx_model=_EN


less -S ${dir}/Naive_pred_EN.txt


ABCA7

dataset=Naive


# REMOVE UNNECCESSARY COLUMNS
head -n1 ${dir}/${dataset}_pred${GTEx_model}.txt | sed 's/CHROM\tGeneStart\tGeneEnd\t/dataset\t/g; s/\tGeneName//g' > ${dir}/${dataset}_pred${GTEx_model}_2.txt
tail -n+2 ${dir}/${dataset}_pred${GTEx_model}.txt | sed 's/19\t1040101\t1065572\t/'${dataset}'\t/g; s/\tABCA7//g' >> ${dir}/${dataset}_pred${GTEx_model}_2.txt



header=`head -n1 ${dir}/Naive_pred_EN.txt`
echo -e "dataset\t"${header} > ${dir}/all_pred.txt


awk 'BEGIN { FS = OFS = "\t" }{

}'

echo -e "dataset\t" > ${dir}/all_pred.txt

for dataset in ${datasets[@]}; do
	awk -v dataset=${dataset} 'BEGIN { FS = OFS = "\t" }{
		# print header
		if ($0 ~ /^CHROM/){
			next
		} else {
			print dataset"\t"$0
		}
	}' ${dir}/${dataset}_pred${GTEx_model}.txt >> ${dir}/${dataset}_pred${GTEx_model}_2.txt

done


grep -m1 "13"$'\t'"1010"$'\t' ${dir}/Naive_pred_EN_melt.txt

grep -m1 "1000000"$'\t'"SM-CTENG"$'\t' ${dir}/Naive_pred_EN_melt.txt



