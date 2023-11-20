db=/public/home/2022122/xugang/project/tair10/tair10.bowtie2
bed=/public/home/2022122/xugang/project/tair10/tair10.gene.bed
bed2=/public/home/2022122/xugang/project/tair10/tair10.mRNA.bed
bed3=/public/home/2022122/xugang/project/tair10/tair10.ncRNA.bed
bed4=/public/home/2022122/xugang/project/tair10/tair10.lncRNA.bed
bed5=/public/home/2022122/xugang/project/tair10/tair10.miRNA.bed
bed6=/public/home/2022122/xugang/project/tair10/tair10.rRNA.bed
bed7=/public/home/2022122/xugang/project/tair10/tair10.snoRNA.bed
bed8=/public/home/2022122/xugang/project/tair10/tair10.snRNA.bed
bed9=/public/home/2022122/xugang/project/tair10/tair10.tRNA.bed
datapath2=`pwd`/rawdata
output=`pwd`/output
thread=50
partnum=6
#node=Fnode2
#node=Fnode1
node=Cnode
node=Cnode2
#node=Gnode
counter=0
#################

bowtief(){
#map to tair10 genome
echo '' 
}

samtools(){
#convert sam to bam
#sort bam
#build index
echo '' 
}
macscallpeak(){
#call peaks
log=$output/a3-callpeak/log
[[ -d $log ]] || mkdir -p  $log

#remove background file
((counter++))
name1=$1
name2=$2

echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 3${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc
macs2 callpeak  -B -t $output/a2-bam/${name1}.bam -c $output/a2-bam/${name2}.bam -n cond1 --nomodel --extsize 120 --outdir $output/a3-callpeak/ --name ${name1}

">a3.callpeak.$counter.$name1.$name2.sh


}
macspeakdiff(){
log=$output/a8-macspeakdiff/log
[[ -d $log ]] || mkdir -p  $log

#remove background file
((counter++))
name1=$1
name2=$2

echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 8${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc
macs2 bdgdiff --t1 $output/a3-callpeak/${name1}_treat_pileup.bdg --c1 $output/a3-callpeak/${name1}_control_lambda.bdg --t2 $output/a3-callpeak/${name2}_treat_pileup.bdg --c2 $output/a3-callpeak/${name2}_control_lambda.bdg -g 60 -l 120 --outdir $output/a8-macspeakdiff/ --o-prefix ${name1}.vs.${name2}
">a8.peakdiff.$counter.$name1.sh


}
bamtobw(){
log=$output/a4-bw/log
[[ -d $log ]] || mkdir -p  $log

#remove background file
((counter++))
name1=$1
name2=$2

echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 4${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc
conda run -n deeptool bamCoverage -b ${output}/a2-bam/${name1}.bam -of bigwig -o ${output}/a4-bw/${name2}.bw -p $((thread)) --ignoreDuplicates --binSize 1000 --normalizeUsing RPKM

">a5.$counter.$name1.sh
}


bigcomparef(){
log=$output/a4-bw/log
[[ -d $log ]] || mkdir -p  $log

#remove background file
((counter++))
name1=$1
name2=$2

echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 4${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc
conda run -n deeptool bamCompare -b1 ${output}/a2-bam/${name1}.bam -b2  $output/a2-bam/${name2}.bam -o ${output}/a4-bw/${name1}.log2.bw -p $((thread*2)) --ignoreDuplicates --binSize 1000 
#bigWigToBedGraph ${output}/a4-bw/${name1}.log2.bw ${output}/a4-bw/${name1}.log2.bed

">a4.$counter.$name1.sh
}
bwcomparef(){
log=$output/a4-bw/log
[[ -d $log ]] || mkdir -p  $log

#remove background file
((counter++))
name1=$1
name2=$2

echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 4${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc
conda run -n deeptool bigwigCompare -b1 ${output}/a4-bw/${name1}.log2.bw -b2  $output/a4-bw/${name2}.log2.bw -o ${output}/a4-bw/${name1}.vs.${name2}.log2.bw -p $((thread*2)) --binSize 1000


">a5.$counter.$name1.sh





}

computmatrix(){
#calculate the matrix to calculate
log=$output/a5-matrix/log
[[ -d $log ]] || mkdir -p  $log

#remove background file
((counter++))
name1=$1
echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 5${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc

conda run -n deeptool computeMatrix scale-regions -R ${bed2} ${bed3} ${bed4} -S  ${output}/a4-bw/${name1}.log2.bw --smartLabels -p $((thread)) --binSize 10 -b 3000 -a 3000 --regionBodyLength 5000 --sortRegions keep -o $output/a5-matrix/${name1}.gz --outFileSortedRegions $output/a5-matrix/computeMatrix_${name1}.bed --outFileNameMatrix $output/a5-matrix/matrix_${name1}.tab
">a5.computematrix.$counter.$name1.sh
}

computmatrixpoint(){
#calculate the matrix to calculate
log=$output/a5-matrix/log
[[ -d $log ]] || mkdir -p  $log

#remove background file
((counter++))
name1=$1
echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 5${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc

conda run -n deeptool computeMatrix reference-point  --referencePoint TSS -R ${bed} -S  ${output}/a4-bw/${name1}.bw --smartLabels -p $((thread)) -b 3000 -a 3000  --skipZeros -o $output/a5-matrix/${name1}.gz --outFileSortedRegions $output/a5-matrix/computeMatrix_${name1}.bed
">a5.computematrix.$counter.$name1.sh
}

plotprofile(){
#use deeptools plot data profiles

log=$output/a6-profile/log
[[ -d $log ]] || mkdir -p  $log
#remove background file
((counter++))
name1=$1
echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 5${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc

conda run -n deeptool plotProfile -m  $output/a5-matrix/${name1}.gz -out $output/a6-profile/Profile_${name1}.eps --outFileNameData $output/a6-profile/plotProfile_${name1}.tab
conda run -n deeptool plotProfile -m  $output/a5-matrix/${name1}.gz -out $output/a6-profile/Profile_${name1}.pdf --outFileNameData $output/a6-profile/plotProfile_${name1}.tab
">a6.plotprofile.$counter.$name1.sh 
}

plotheatmap(){
#plot peak heatmap.
log=$output/a7-heatmap/log
[[ -d $log ]] || mkdir -p  $log

#remove background file
((counter++))
name1=$1


echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 5${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc

conda run -n deeptool plotHeatmap --heatmapWidth 12 --heatmapHeight 50 --zMax 10 --colorList \" #4393C3,white,#A50026 \" --missingDataColor white -m  $output/a5-matrix/${name1}.gz -out $output/a7-heatmap/${name1}_Heatmap.eps --boxAroundHeatmaps no
conda run -n deeptool plotHeatmap --heatmapWidth 12 --heatmapHeight 50 --zMax 10 --colorList \" #4393C3,white,#A50026 \" --missingDataColor white -m  $output/a5-matrix/${name1}.gz -out $output/a7-heatmap/${name1}_Heatmap.pdf --boxAroundHeatmaps no
">a7.heatmap.$counter.$name1.sh

}
#macscallpeak WT_2_1_IP.unique WT_2_1_input.unique
#macscallpeak WT_1_1_IP.unique WT_1_1_input.unique
#macscallpeak B_1_1_IP.unique B_1_1_input.unique
#macspeakdiff WT_2_1_IP.unique  WT_1_1_IP.unique
#macspeakdiff col_1_IP B_1_1_IP.unique
#bigcomparef col_1_IP col_1_Input 
#bigcomparef col_2_IP col_2_Input
#bigcomparef M_2_1_IP.unique M_2_1_input.unique
#bigcomparef M_1_1_IP.unique M_1_1_input.unique
#bigcomparef WT_2_1_IP.unique WT_2_1_input.unique
#bigcomparef WT_1_1_IP.unique WT_1_1_input.unique
#bamtobw WT_2_1_IP.unique wt_smxl7_h3k27me2
#bamtobw M_2_1_IP.unique mut_smxl7_h3k27me2

bamtobw  ATAC_mut_smxl7_2.unique  ATAC_mut_smxl7_2 
bamtobw  ATAC_oe_smxl7_1.unique ATAC_oe_smxl7_1
#bamtobw  ATAC_oe_smxl7_2.unique ATAC_oe_smxl7_2
#bamtobw  ATAC_WT_1.unique ATAC_WT_1
#bamtobw  ATAC_WT_2.unique ATAC_WT_2

#bwcomparef M_2_1_IP.unique M_1_1_IP.unique
#bwcomparef WT_2_1_IP.unique WT_1_1_IP.unique

#computmatrix col_1_IP
#computmatrix col_2_IP
#computmatrix M_2_1_IP.unique.vs.M_1_1_IP.unique
#computmatrix WT_2_1_IP.unique.vs.WT_1_1_IP.unique

#computmatrixpoint wt_smxl7_h3k27me2
#computmatrixpoint mut_smxl7_h3k27me2
#computmatrixpoint ATAC_mut_smxl7_2
#computmatrixpoint ATAC_oe_smxl7_1
#computmatrixpoint ATAC_oe_smxl7_2
#computmatrixpoint ATAC_WT_1
#computmatrixpoint ATAC_WT_2

#plotprofile col_1_IP
#plotprofile col_2_IP
#plotprofile M_2_1_IP.unique.vs.M_1_1_IP.unique
#plotprofile WT_2_1_IP.unique.vs.WT_1_1_IP.unique
#plotprofile wt_smxl7_h3k27me2
#plotprofile mut_smxl7_h3k27me2 


plotprofile ATAC_mut_smxl7_2
plotprofile ATAC_oe_smxl7_1
plotprofile ATAC_oe_smxl7_2
plotprofile  ATAC_WT_1
plotprofile  ATAC_WT_2


#plotheatmap col_1_IP
#plotheatmap col_2_IP
#plotheatmap WT_2_1_IP.unique.vs.WT_1_1_IP.unique
#plotheatmap M_2_1_IP.unique.vs.M_1_1_IP.unique
#plotheatmap wt_smxl7_h3k27me2
#plotheatmap mut_smxl7_h3k27me2 

plotheatmap ATAC_mut_smxl7_2
plotheatmap ATAC_oe_smxl7_1
plotheatmap ATAC_oe_smxl7_2
plotheatmap ATAC_WT_1
plotheatmap ATAC_WT_2



peakoverlap(){
cd /public/home/2022122/xugang/project/yaoruifeng/output/a3-callpeak  
bedtools window -a col_1_IP_summits.bed -b col_2_IP_summits.bed -w 120 |wc -l
bedtools window -a col_1_IP_summits.bed -b col_2_IP_summits.bed -w 120  -v |wc -l
bedtools window -a col_2_IP_summits.bed -b col_1_IP_summits.bed -w 120  -v |wc -l
}

#na1=col_1_IP_summits.bed
#wc -l ${na1}
#bedtools intersect -a ${na1} -b /public/home/2022122/xugang/project/tair10/tair10.exon.bed -wa |sort -u | wc -l
#bedtools intersect -a ${na1} -b /public/home/2022122/xugang/project/tair10/tair10.exon.bed -wa -v > ${na1}.r1.bed








