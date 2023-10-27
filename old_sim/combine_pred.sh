

sim_dir=/home/rparrish/YangFSSdata/SR_TWAS/sim


# ###### 
# dir=${sim_dir}/power

# cp ${dir}/all_pred_results_EN_old.txt ${dir}/all_pred_results_EN.txt
# tail -n+2 ${dir}/all_pred_results_EN3.txt >> ${dir}/all_pred_results_EN.txt

# cp ${dir}/all_power_results_EN_old.txt ${dir}/all_power_results_EN.txt
# tail -n+2 ${dir}/all_power_results_EN3.txt >> ${dir}/all_power_results_EN.txt



###### 
dir=${sim_dir}/pred


# ln -s ${dir}/GTEx_EN_pred.txt ${dir}/GTEx_pred_EN.txt
# ln -s ${dir}/GTEx_EN_pred_120621.txt ${dir}/GTEx_pred_EN_120621.txt
# ln -s ${dir}/GTEx_EN_pred_121721.txt ${dir}/GTEx_pred_EN_121721.txt
# ln -s ${dir}/GTEx_EN_pred_overlap.txt ${dir}/GTEx_pred_EN_overlap.txt

datasets=(GTEx ROSMAP Naive SR)

# add header
head -n1 ${dir}/${datasets[0]}_pred.txt | \
awk 'BEGIN { FS = OFS = "\t" } {print "dataset\tsuffix\t"$0'} \
> ${dir}/all_pred_results.txt

# combine prediction results
for dataset in ${datasets[@]}; do

	for suffix in '' _120621 _121721 _overlap _EN _EN_123121 _EN_120621 _EN_121721 _EN_overlap; do
		path=${dir}/${dataset}_pred${suffix}.txt
		if [[ "$suffix"x == ""x ]]; then
			suffix2=NA
		else
			suffix2=${suffix}
		fi

		tail -n+2 ${path} | \
		awk -v dataset=${dataset} -v suffix=${suffix2} 'BEGIN { FS = OFS = "\t" } {print dataset"\t"suffix"\t"$0}' \
		>> ${dir}/all_pred_results.txt

	done

done


# _EN_123121




###### 
# dir=${sim_dir}/power
# # cp ${dir}/all_power_results_ENoverlap.txt ${dir}/all_power_results_EN.txt
# tail -n+2 ${dir}/all_power_results_ENoverlap.txt >> ${dir}/all_power_results_EN.txt


# dir=${sim_dir}/power
# # cp ${dir}/all_pred_results_EN1.txt ${dir}/all_pred_results_EN.txt
# tail -n+2 ${dir}/all_pred_results_ENoverlap.txt >> ${dir}/all_pred_results_EN.txt



# dir=${sim_dir}/power
# cp ${dir}/all_power_results_EN1.txt ${dir}/all_power_results_EN.txt
# tail -n+2 ${dir}/all_power_results_EN2.txt >> ${dir}/all_power_results_EN.txt


# dir=${sim_dir}/power
# cp ${dir}/all_pred_results_EN1.txt ${dir}/all_pred_results_EN.txt
# tail -n+2 ${dir}/all_pred_results_EN2.txt >> ${dir}/all_pred_results_EN.txt


# ln -s ${dir}/GTEx_EN_pred.txt ${dir}/GTEx_pred_EN.txt

# add header
# cp ${dir}/all_power_results_EN1.txt ${dir}/all_power_results_EN.txt
# tail -n+2 ${dir}/all_power_results_EN2.txt >> ${dir}/all_power_results_EN.txt



# for dataset in GTEx ROSMAP Naive SR; do
# 	for dataset in ${datasets[@]}; do
# 		for suffix in '' _120621 _121721; do
# 			tail -n+2 ${dir}/${dataset}_pred${suffix}.txt | \
# 			awk -v dataset=${dataset} suffix=${suffix} 'BEGIN { FS = OFS = "\t" } {print dataset"\t"suffix"\t"$0}' \
# 			>> ${dir}/all_pred_results.txt
# 		done
# 	done

# 	for suffix in '' _120621 _121721; do


# 		tail -n+2 ${dir}/${dataset}_pred${suffix}.txt \
# 		>> ${dir}/${dataset}_pred_all.txt
# 	done

# done



# 	for suffix in '' _120621 _121721; do
# 		echo $suffix
# 	done

######

# pred
for dataset in GTEx ROSMAP Naive SR; do
	cp ${dir}/${dataset}_pred.txt ${dir}/${dataset}_pred_all.txt

	for suffix in _120621 _121721; do
		tail -n+2 ${dir}/${dataset}_pred${suffix}.txt \
		>> ${dir}/${dataset}_pred_all.txt
	done

done