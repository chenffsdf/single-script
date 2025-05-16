#!/bin/bash -e
# ================================================================
# 脚本名称: 05.suppa.sh
# 功能描述: 
#   使用 SUPPA2 工具进行可变剪接分析，生成事件注释(ioe)和事件量化(ioi)文件
# 主要步骤:
#   1. 根据输入的GTF文件生成可变剪接事件注释
#   2. 创建输出目录结构
# 适用环境:
#   - 需要在SLURM集群的CU分区运行
#   - 需要预先配置SUPPA2的conda环境
# 
# 输入参数:
#   $1 : 输入GTF文件路径 (e.g. /data/genome.gtf)
#   $2 : 输出目录路径 (e.g. /results/suppa_out)
#
# 输出文件:
#   - {output_dir}/ioe_result/ : 存放事件注释文件(.ioe)
#   - {output_dir}/ioi_result/ : 存放事件量化文件(.ioi)
#
# 资源需求:
#   - 16 CPU核心
#   - 5GB 内存
#
# 作者: guowenbin & chenjiawei
# 创建日期: 20250516
# 版本: 1.0
# 
# 使用示例:
#   sbatch 05.suppa.sh genome.gtf ./suppa_results
#
# 注意事项:
#   1. 需要先加载conda环境
#   2. 输出目录会自动创建
#   3. 运行失败会发送邮件通知(END/FAIL状态)
# ================================================================

#SBATCH --partition=cu
#SBATCH --export=all
#SBATCH --job-name=suppa_cjw
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=chenjiawei@higentec.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o %x_%j.out
#SBATCH --mem=5G
source /home/chenjw/miniconda3/etc/profile.d/conda.sh
conda activate suppa

echo "This script is running on $HOSTNAME"
echo "Start time: $(date)"
start=$SECONDS

gtf_file=$1
# tpm=$2
output_dir=$2
# format="ioe"

mkdir -p ${output_dir}/ioe_result
mkdir -p ${output_dir}/ioi_result
##########################################################################
  
echo $output_dir

# case "$format" in
# "ioe")
echo Suppa for local events $format

suppa.py generateEvents \
-i $gtf_file \
-o ${output_dir}/ioe_result/events \
-f ioe \
-e SE SS MX RI FL
#;;

# "ioi")
echo Suppa for transcript events $format

suppa.py generateEvents \
-i $gtf_file \
-o ${output_dir}/ioi_result/events \
-f ioi

# esac

# tpm_dir=$(dirname $tpm)

# head -1 $tpm|awk -F "\t" -v OFS="\t" '{print $1,$2,$3}' > ${tpm_dir}/control.tpm
# cat $tpm|tail -n +2|awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$4}' >> ${tpm_dir}/control.tpm
# head -1 $tpm|awk -F "\t" -v OFS="\t" '{print $4,$5,$6}' > ${tpm_dir}/treat.tpm
# cat $tpm|tail -n +2|awk -F "\t" -v OFS="\t" '{print $1,$5,$6,$7}' >> ${tpm_dir}/treat.tpm

# for event in ioi SE A3 A5 MX RI AF AL; do
# 	mkdir -p ${output_dir}/psi_result/${event}
# 	mkdir -p ${output_dir}/diff_result/${event}
# 	if [[ "$event" == "ioi" ]]; then
# 		suppa.py psiPerEvent \
# 		-i ${output_dir}/ioi_result/events.ioi \
# 		-e $tpm \
# 		-o ${output_dir}/psi_result/${event}/${event}
# 	else
# 		suppa.py psiPerEvent \
# 		-i ${output_dir}/ioe_result/events_${event}_strict.ioe \
# 		-e $tpm \
# 		-o ${output_dir}/psi_result/${event}/${event}
# 	fi

# 	head -1 ${output_dir}/psi_result/${event}/${event}.psi|awk -F "\t" -v OFS="\t" '{print $1,$2,$3}' > ${output_dir}/psi_result/${event}/control.psi
# 	cat ${output_dir}/psi_result/${event}/${event}.psi|tail -n +2|awk -F "\t" -v OFS="\t" '{print $1,$2,$3,$4}' >> ${output_dir}/psi_result/${event}/control.psi
# 	head -1 ${output_dir}/psi_result/${event}/${event}.psi|awk -F "\t" -v OFS="\t" '{print $4,$5,$6}' > ${output_dir}/psi_result/${event}/treat.psi
# 	cat ${output_dir}/psi_result/${event}/${event}.psi|tail -n +2|awk -F "\t" -v OFS="\t" '{print $1,$5,$6,$7}' >> ${output_dir}/psi_result/${event}/treat.psi

# 	if [[ "$event" == "ioi" ]]; then
# 		suppa.py diffSplice \
# 		-m empirical -gc \
# 		-i ${output_dir}/ioi_result/events.ioi \
# 		--save_tpm_events \
# 		-p ${output_dir}/psi_result/${event}/treat.psi ${output_dir}/psi_result/${event}/control.psi \
# 		-e ${tpm_dir}/treat.tpm ${tpm_dir}/control.tpm \
# 		-o ${output_dir}/diff_result/${event}/diff_${event}
# 	else
# 		suppa.py diffSplice \
# 		-m empirical -gc \
# 		-i ${output_dir}/ioe_result/events_${event}_strict.ioe \
# 		--save_tpm_events \
# 		-p ${output_dir}/psi_result/${event}/treat.psi ${output_dir}/psi_result/${event}/control.psi \
# 		-e ${tpm_dir}/treat.tpm ${tpm_dir}/control.tpm \
# 		-o ${output_dir}/diff_result/${event}/diff_${event}
# 	fi
# done


###############################################################################################
###--->Time of the analysis
echo "done"
end=$SECONDS
taken=$(($end-$start))
echo "Time taken: $(($taken/3600)) hours $((($taken/60)%60)) minutes $(($taken%60)) seconds"
echo "End time: $(date)"
