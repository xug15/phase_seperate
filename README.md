# phase seperate project

## Content
* [Sample information](#Sample-information)  
### CHIP-seq, ATAC-seq
* [Set up environment variables](#1-set-up-environment-variables)   
### mRNA-seq
* [Calculate the gene count.](#calculate-the-gene-count)  
* [Differential expression gene](script/gene.df.r)   
### Hi-C seq
* [Hi-C seq information](hic.md)
* [Hi-C process](#hi-c-library-preparation-sequencing-and-analysis)
### Bisulfite-seq
* [Bisulfite-seq analysis](#bisulfite-seq-analysis)


## Sample information 
M-1-1    突变体 H3 ChIP-seq (SMX7)  
M-2-1    突变体 H3K27me3 ChIP-seq (SMX7)  
WT_1_1    野生型 H3 ChIP-seq (SMX7)  
WT_2_1    野生型 H3K27me3 ChIP-seq (SMX7)  
O_1    过表达 smxl7 H3  
O_2    过表达 smxl7 H3k27me3  
O-3-1    smax7 -(转录因子) GFP ChIP-seq  
O-4-1    smax7   -(转录因子) GFP ChIP-seq 重复  
B_1_1    空白对照 GFP ChIP-seq 

col_1 DAP-seq Chip-seq SMXL7
col_2 DAP-seq Chip-seq SMXL7 

 
ATAC_M_1   突变体  (SMX7)  
ATAC_M_2   突变体重复 (SMX7)  
ATAC_O_1   过表达 (SMX7)  
ATAC_O_2   过表达重复 (SMX7)  
ATAC_WT_1 野生型 (SMX7)  
ATAC_WT_2 野生型重复 (SMX7)  

## 0. set conda, pip, mamba
```sh

 pip install --upgrade --force-reinstall mamba
 

```
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")

```

## 1. Set up environment variables.
```sh
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

```

##   
```sh
fastqcf(){
((counter++))
name=$1
path=$2
file='a0-fastqc'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 1${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date
source /public/home/2022122/xugang/bashrc
fastqc -t ${thread} ${path}/${name} -o $output/${file}

" > a0.fastqc.$counter.${name}.sh
}

fastqcf mod.A-1_1.fq.gz rawdata/modification/
fastqcf mod.A-1_2.fq.gz rawdata/modification/
fastqcf mod.A-2_1.fq.gz rawdata/modification/
fastqcf mod.A-2_2.fq.gz rawdata/modification/
fastqcf mod.A-3_1.fq.gz rawdata/modification/
fastqcf mod.A-3_2.fq.gz rawdata/modification/
fastqcf mod.Mu-1_1.fq.gz rawdata/modification/
fastqcf mod.Mu-1_2.fq.gz rawdata/modification/
fastqcf mod.Mu-2_1.fq.gz rawdata/modification/
fastqcf mod.Mu-2_2.fq.gz rawdata/modification/
fastqcf mod.Mu-3_1.fq.gz rawdata/modification/
fastqcf mod.Mu-3_2.fq.gz rawdata/modification/
fastqcf mod.OE-1_1.fq.gz rawdata/modification/
fastqcf mod.OE-1_2.fq.gz rawdata/modification/
fastqcf mod.OE-2_1.fq.gz rawdata/modification/
fastqcf mod.OE-2_2.fq.gz rawdata/modification/
fastqcf mod.OE-3_1.fq.gz rawdata/modification/
fastqcf mod.OE-3_2.fq.gz rawdata/modification/
fastqcf mod.WT-1_1.fq.gz rawdata/modification/
fastqcf mod.WT-1_2.fq.gz rawdata/modification/
fastqcf mod.WT-2_1.fq.gz rawdata/modification/
fastqcf mod.WT-2_2.fq.gz rawdata/modification/
fastqcf mod.WT-3_1.fq.gz rawdata/modification/
fastqcf mod.WT-3_2.fq.gz rawdata/modification/
#
fastqcf hic_mut_2_1.fq.gz rawdata/hic-seq
fastqcf hic_mut_2_2.fq.gz rawdata/hic-seq
fastqcf hic_oe_smxl7_2_1.fq.gz rawdata/hic-seq
fastqcf hic_oe_smxl7_2_2.fq.gz rawdata/hic-seq
fastqcf hic_WT_2_1.fq.gz rawdata/hic-seq
fastqcf hic_WT_2_2.fq.gz rawdata/hic-seq
```



## Using the bowtie. Implemete the read into genome sequence.
Fistly using software align the reads into genome, the filter the result and transfor formate into bam files.
sort the bam.

```sh

bowtief(){
name=$1
((counter++))
file='a1-map'
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 1${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date
source /public/home/2022122/xugang/bashrc
bowtie2 -p  ${thread} --end-to-end --sensitive-local -x ${db} -1 $output/a0-clean/$name.1.fq.gz -2 $output/a0-clean/$name.2.fq.gz  --un-conc-gz  $output/$file/${name}.decontaminate.fq.gz -S $output/$file/${name}.sam >  $output/$file/${name}.txt 2>  $output/$file/${name}.txt

" > a1.mapf.$counter.${name}.sh
}
```

```sh
samtobam(){
#convert sam to bam
#sort bam
#build index
name=$1
name2=$2
file='a2-bam'
#remove background file
((counter++))
[[ -d $output/${file}/log ]] || mkdir -p $output/${file}/log
echo -e "#!/bin/bash
#SBATCH -o ${output}/$file/log/${name}.%j.out
#SBATCH -e ${output}/$file/log/${name}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 2${name}
#SBATCH -N 1
#SBATCH -n ${thread}
echo date
source /public/home/2022122/xugang/bashrc

samtools view -q 30 --threads ${thread} -bhS $output/a1-bowtie2/${name}.sam > $output/$file/${name}.bam
samtools sort --threads ${thread} $output/$file/${name}.bam > $output/$file/${name}.sort.bam 
rm  $output/$file/${name}.bam
mv $output/$file/${name}.sort.bam $output/$file/${name2}.unique.bam
samtools index $output/$file/${name2}.unique.bam
" > a2.samtobam.$counter.${name}.sh 
}
samtobam col_1_Input dap-smxl7-r1-input
samtobam col_1_IP dap-smxl7-r1-ip
samtobam col_2_Input dap-smxl7-r2-input
samtobam col_2_IP dap-smxl7-r2-ip
```

## Using the macs2 software to call peaks.
the output files contain the peak locations.

```sh
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


macscallpeak WT_2_1_IP.unique WT_2_1_input.unique
macscallpeak WT_1_1_IP.unique WT_1_1_input.unique
macscallpeak B_1_1_IP.unique B_1_1_input.unique
```
## calculate the difference of peaks between the different conditions with macs software.

```sh
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
macspeakdiff WT_2_1_IP.unique WT_1_1_IP.unique
macspeakdiff col_1_IP B_1_1_IP.unique
```

## Implemete the file into bigwig files.
Just tranform the bam to bigwig.

```sh
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

bamtobw ATAC_mut_smxl7_1.unique ATAC_mut_SMXL7_1
bamtobw ATAC_mut_smxl7_2.unique ATAC_mut_SMXL7_2
bamtobw ATAC_oe_smxl7_1.unique ATAC_oe_SMXL7_1
bamtobw ATAC_oe_smxl7_2.unique ATAC_oe_SMXL7_2
bamtobw ATAC_WT_1.unique ATAC_WT_1
bamtobw ATAC_WT_2.unique ATAC_WT_2

bamtobw mRNA.M_1.sorted mRNA.mut_1
bamtobw mRNA.M_2.sorted mRNA.mut_2
bamtobw mRNA.M_3.sorted mRNA.mut_3
bamtobw mRNA.O_1.sorted mRNA.oe_1
bamtobw mRNA.O_2.sorted mRNA.oe_2
bamtobw mRNA.O_3.sorted mRNA.oe_3
bamtobw mRNA.WT_1.sorted mRNA.wt_1
bamtobw mRNA.WT_2.sorted mRNA.wt_2
bamtobw mRNA.WT_3.sorted mRNA.wt_3
}


```
## remove background signal, and generate the output file with bigwig file format.  
* Compare the signal between two conditions by deeptools bamCompare.  



```sh
bigcomparef(){
log=$output/a4-bw/log
[[ -d $log ]] || mkdir -p  $log
#remove background file
((counter++))
name1=$1
name2=$2
name3=$3 
echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 4${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc
conda run -n deeptool bamCompare -b1 ${output}/a2-bam/${name1}.bam -b2  $output/a2-bam/${name2}.bam -o ${output}/a4-bw/${name3}.rmbg.bw -p $((thread)) --ignoreDuplicates --binSize 1000 
bigWigToWig ${output}/a4-bw/${name3}.rmbg.bw ${output}/a4-bw/${name3}.rmbg.wig

">a4.$counter.$name1.sh
}
bigcomparef dap-smxl7-r1-ip.unique dap-smxl7-r1-input.unique dap-smxl7-r1 
bigcomparef dap-smxl7-r2-ip.unique dap-smxl7-r2-input.unique dap-smxl7-r2
bigcomparef M_h3_1_IP.unique M_h3_1_input.unique M_h3
bigcomparef M_h3k27me3_1_IP.unique M_h3k27me3_1_input.unique M_h3k27me3
bigcomparef wt_h3_IP.unique wt_h3_input.unique wt_h3 
bigcomparef wt_h3k27me3_IP.unique wt_h3k27me3_input.unique wt_h3k27me3 
bigcomparef oe_h3k27me3_IP.unique oe_h3k27me3_input.unique oe_h3k27me3
bigcomparef oe_h3_IP.unique oe_h3_input.unique oe_h3
bigcomparef wt_h3k27me3_IP.unique wt_h3_IP.unique wt_h3k27me3_h3_withoutbackgorund
```


## Calculate the H3K27me3/H3 signal.
* Compare the bigwig files by deeptools bigwigCompare funciton.

```sh
bwcomparef(){
log=$output/a4-bw/log
[[ -d $log ]] || mkdir -p  $log
#remove background file
((counter++))
name1=$1
name2=$2
name3=$3 
echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 4${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc
conda run -n deeptool bigwigCompare -b1 ${output}/a4-bw/${name1}.rmbg.bw -b2  $output/a4-bw/${name2}.rmbg.bw -o ${output}/a4-bw/${name3}.bw -p $((thread)) --binSize 1000
">a5.$counter.$name1.sh
}
bwcomparef M_h3k27me3 M_h3 M_h3k27me3.h3
bwcomparef oe_h3k27me3 oe_h3 oe_h3k27me3.h3
bwcomparef wt_h3k27me3 wt_h3 wt_h3k27me3.h3 

```

## Compute the matrix of gene body and get the matrix.

```sh
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
computmatrix col_1_IP
computmatrix col_2_IP
computmatrix M_2_1_IP.unique.vs.M_1_1_IP.unique
computmatrix WT_2_1_IP.unique.vs.WT_1_1_IP.unique
```
## multiple input files (scale-regions mode)

```sh
computmatrix_multi(){
#calculate the matrix to calculate
log=$output/a5-matrix/log
[[ -d $log ]] || mkdir -p  $log

#remove background file
((counter++))
name1=$1
name2=$2
name3=$3
name4=$4
echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 5${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc

conda run -n deeptool computeMatrix scale-regions -R ${bed2} ${bed3} -S ${output}/a4-bw/${name1}.bw ${output}/a4-bw/${name2}.bw ${output}/a4-bw/${name3}.bw --smartLabels -p $((thread)) --binSize 10 -b 3000 -a 3000 --regionBodyLength 5000 --sortRegions keep -o $output/a5-matrix/${name4}.gz --outFileSortedRegions $output/a5-matrix/computeMatrix_${name4}.bed --outFileNameMatrix $output/a5-matrix/matrix_${name4}.tab
">a5.computematrix.multi.$counter.$name1.sh
}
computmatrix_multi wt_h3k27me3.h3 M_h3k27me3.h3 oe_h3k27me3.h3 wt.m.oe.h3k27.h3
```
## Compute the signal density cover reference point(TSS) and get the matrix of reference points.

```sh
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
computmatrixpoint wt_smxl7_h3k27me2
computmatrixpoint mut_smxl7_h3k27me2
```
## Compute the signal density cover reference point(TSS) with multiple files and get the matrix of reference points.

```sh
computmatrixpoint_multi(){
#calculate the matrix to calculate
log=$output/a5-matrix/log
[[ -d $log ]] || mkdir -p  $log

#remove background file
((counter++))
name1=$1
name2=$2
name3=$3
name4=$4
echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 5${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc

conda run -n deeptool computeMatrix reference-point  --referencePoint TSS -R ${bed} -S ${output}/a4-bw/${name1}.bw ${output}/a4-bw/${name2}.bw ${output}/a4-bw/${name3}.bw --smartLabels -p $((thread)) -b 3000 -a 3000  --skipZeros -o $output/a5-matrix/${name4}.point.gz --outFileSortedRegions $output/a5-matrix/computeMatrix.${name4}.point.bed
">a5.computematrix.point.$counter.${name4}.sh
}
computmatrixpoint_multi wt_h3k27me3.h3 M_h3k27me3.h3 oe_h3k27me3.h3 wt.m.oe.h3k27.h3
```


## Plot the signal distribution of gene body or reference points such as TSS, TTS.

```sh
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
plotprofile wt.m.oe.h3k27.h3
plotprofile wt.m.oe.h3k27.h3.point

plotprofile col_1_IP
plotprofile col_2_IP
plotprofile M_2_1_IP.unique.vs.M_1_1_IP.unique
plotprofile WT_2_1_IP.unique.vs.WT_1_1_IP.unique
plotprofile wt_smxl7_h3k27me2
plotprofile mut_smxl7_h3k27me2 
plotprofile col_1_IP.M_2_1_IP.unique.vs.M_1_1_IP.unique.WT_2_1_IP.unique.vs.WT_1_1_IP.unique

```

## Plot the heatmap of Chip-seq signal coverage in the gene body or reference points.

```sh
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

conda run -n deeptool plotHeatmap --heatmapWidth 12 --heatmapHeight 50 --zMax 2 --colorList \" #4393C3,white,#A50026 \" --missingDataColor white -m  $output/a5-matrix/${name1}.gz -out $output/a7-heatmap/${name1}_Heatmap.eps --boxAroundHeatmaps no
conda run -n deeptool plotHeatmap --heatmapWidth 12 --heatmapHeight 50 --zMax 2 --colorList \" #4393C3,white,#A50026 \" --missingDataColor white -m  $output/a5-matrix/${name1}.gz -out $output/a7-heatmap/${name1}_Heatmap.pdf --boxAroundHeatmaps no
">a7.heatmap.$counter.$name1.sh
}
plotheatmap wt.m.oe.h3k27.h3
plotheatmap wt.m.oe.h3k27.h3.point

plotheatmap col_1_IP
plotheatmap col_2_IP
plotheatmap WT_2_1_IP.unique.vs.WT_1_1_IP.unique
plotheatmap M_2_1_IP.unique.vs.M_1_1_IP.unique
plotheatmap wt_smxl7_h3k27me2
plotheatmap mut_smxl7_h3k27me2 
plotheatmap col_1_IP.M_2_1_IP.unique.vs.M_1_1_IP.unique.WT_2_1_IP.unique.vs.WT_1_1_IP.unique
```

## Impletement the two replicates peaks.

```sh
peakoverlap(){
cd /public/home/2022122/xugang/project/yaoruifeng/output/a3-callpeak  
bedtools window -a col_1_IP_summits.bed -b col_2_IP_summits.bed -w 120 |wc -l
bedtools window -a col_1_IP_summits.bed -b col_2_IP_summits.bed -w 120  -v |wc -l
bedtools window -a col_2_IP_summits.bed -b col_1_IP_summits.bed -w 120  -v |wc -l
}
```

##

```sh
na1=col_1_IP_summits.bed
wc -l ${na1}
bedtools intersect -a ${na1} -b /public/home/2022122/xugang/project/tair10/tair10.exon.bed -wa |sort -u | wc -l
bedtools intersect -a ${na1} -b /public/home/2022122/xugang/project/tair10/tair10.exon.bed -wa -v > ${na1}.r1.bed
```

## Calculate the gene count.

```sh
genecount(){
gtf=/public/home/2022122/xugang/project/tair10/tair10.gtf
name=$1
thread=2
#plot peak heatmap.
log=$output/a9-genecount/log
[[ -d $log ]] || mkdir -p  $log
#remove background file
((counter++))
name1=$1

[[ -d $output/a9-genecount/ ]] || mkdir -p $output/a9-genecount/  

echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 5${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc


samtools view $output/a2-bam/${name}.sorted.bam | gfold count -ann ${gtf} -tag stdin -o $output/a9-genecount/${name}.read_cnt  >$output/a9-genecount/${name}.err

">a9.genecount.$counter.$name.sh
}
genecount mRNA.M_1
genecount mRNA.M_2
genecount mRNA.M_3
genecount mRNA.O_1
genecount mRNA.O_2
genecount mRNA.O_3
genecount mRNA.WT_1
genecount mRNA.WT_2
genecount mRNA.WT_3
```
```sh
mergegf(){
for i in `ls $output/a9-genecount/ |grep .read_cnt$`;
do
file=`echo $i| cut -f 2 -d '.'`;
echo -e "gene\t$file">$output/a9-genecount/$i.txt ;
cut -f 1,3 $output/a9-genecount/$i >> $output/a9-genecount/$i.txt;
done

cd $output/a9-genecount
file=(`ls |grep .read_cnt.txt$`);
begin1=${file[0]};
echo "join $begin1"
cp $begin1 tmp;
echo $begin1;
file2=("${file[@]:1}");
echo ${file2[0]};
for i in ${file2[@]};
do echo $i;
join tmp $i > tmp2
mv tmp2 tmp
done
mv tmp merge.gene.tsv
sed -i 's/ \+/\t/g' merge.gene.tsv
rm *.read_cnt.txt
cd -
}
mergegf
```

## Hi-C library preparation, sequencing and analysis
* Sequencing reads were mapped to the TAIR10 reference genome using the HiC-Pro (v.2.11.1) pipeline90. 
* The bam files (bwt2merged.bam) generated by HiC-Pro containing mapped reads were used as input files for FAN-C (v.0.9.8). 
* The module ‘fanc auto’ was applied to generate 500, 100, 50, 10 and 1 kb contact matrices (hic files). 
* The resultant hic files with 100-kb resolution were directed to the ‘fanc expected’ module to calculate the expected interaction probability against genomic distance for intrachromosomal interaction. 
* For matrix and score comparisons, the default comparison method of fold-change was used with the ‘fanc compare’ command. 
* The outputs (hic object) were transferred to text files by ‘fanc dump’ and were visualized as heatmaps in R using ggplot2. 
* To explore whether higher contacts observed in p35S::H2B.8–eGFP depend on H2B.8 incorporation, the genome was binned into 1-kb windows. 
* H2B.8 signals (log2(IP/input)) of each window were generated and sorted into 20 quantiles by strength. 
* The interaction frequency differences (values generated by FAN-C at 1-kb resolution) of each quantile pair for either short-range interactions or interactions between pericentromeric regions and chromosome arms were averaged and plotted as heatmaps. 
* The resolution of our Hi-C data was estimated as previously reported. 
* Our Hi-C data was deemed to achieve 1-kb resolution as 80% of genomic bins (1 kb) had >1,000 contacts. Our WT data were compared to published contact matrices.

### Prepare environment.

```sh
https://github.com/nservant/HiC-Pro
wget https://github.com/nservant/HiC-Pro/archive/refs/heads/master.zip
mv master.zip hicpro.zip
unzip hicpro.zip
cd HiC-Pro-master/
conda env create -f environment.yml -p hicpro
conda activate hicpro
vi config-install.txt
#PREFIX = /public/home/2022122/xugang/app/HiC-Pro-master
make configure
make install
#export path
export PATH="/public/home/2022122/xugang/app/HiC-Pro-master/bin:"$PATH

```
### The raw file fq.gz have unregnoize characters.
The below script is save the origin data.

```sh
#!/bin/bash
#BATCH -o /public/home/2022122/xugang/project/yaoruifeng/output/a0-fastqc/log/mod.A-1_1.fq.gz.%j.out
#SBATCH -e /public/home/2022122/xugang/project/yaoruifeng/output/a0-fastqc/log/mod.A-1_1.fq.gz.%j.error
#SBATCH --partition=Cnode2
#SBATCH -J 1mod.A-1_1.fq.gz
#SBATCH -N 1
#SBATCH -n 5

echo date
zcat hic_mut_1_1.fq.gz > hic_mut_1_1.fq.gz.zcat
zcat hic_mut_2_1.fq.gz > hic_mut_2_1.fq.gz.zcat
zcat hic_mut_2_2.fq.gz > hic_mut_2_2.fq.gz.zcat
zcat hic_oe_smxl7_1_1.fq.gz > hic_oe_smxl7_1_1.fq.gz.zcat
zcat hic_oe_smxl7_2_1.fq.gz > hic_oe_smxl7_2_1.fq.gz.zcat
zcat hic_oe_smxl7_2_2.fq.gz > hic_oe_smxl7_2_2.fq.gz.zcat
zcat hic_WT_2_1.fq.gz > hic_WT_2_1.fq.gz.zcat
zcat hic_WT_2_2.fq.gz > hic_WT_2_2.fq.gz.zcat

wc -l *zcat
#   269319907 hic_mut_1_1.fq.gz.zcat
#   342152479 hic_mut_2_1.fq.gz.zcat
#   278716231 hic_mut_2_2.fq.gz.zcat
#   256988755 hic_oe_smxl7_1_1.fq.gz.zcat
#   449189785 hic_oe_smxl7_2_1.fq.gz.zcat
#   279289099 hic_oe_smxl7_2_2.fq.gz.zcat
#   286123953 hic_WT_2_1.fq.gz.zcat
#   255456739 hic_WT_2_2.fq.gz.zcat
mkdir log
mv hic*gz log

head -n 269319900 hic_mut_1_1.fq.gz.zcat > hic_mut_1_1.fq
head -n 278716200 hic_mut_2_1.fq.gz.zcat > hic_mut_2_1.fq
head -n 278716200 hic_mut_2_2.fq.gz.zcat > hic_mut_2_2.fq
head -n 256988700 hic_oe_smxl7_1_1.fq.gz.zcat > hic_oe_smxl7_1_1.fq
head -n 279289000 hic_oe_smxl7_2_1.fq.gz.zcat > hic_oe_smxl7_2_1.fq
head -n 279289000 hic_oe_smxl7_2_2.fq.gz.zcat > hic_oe_smxl7_2_2.fq
head -n 255456700 hic_WT_2_1.fq.gz.zcat > hic_WT_2_1.fq
head -n 255456700 hic_WT_2_2.fq.gz.zcat > hic_WT_2_2.fq

gzip -f hic_mut_1_1.fq
gzip -f hic_mut_2_1.fq
gzip -f hic_mut_2_2.fq
gzip -f hic_oe_smxl7_1_1.fq
gzip -f hic_oe_smxl7_2_1.fq
gzip -f hic_oe_smxl7_2_2.fq
gzip -f hic_WT_2_1.fq
gzip -f hic_WT_2_2.fq

cp  hic_*gz /public/home/2022122/xugang/project/yaoruifeng/rawdata/hic-seq
cd /public/home/2022122/xugang/project/yaoruifeng/rawdata/hic-seq
mv hic_mut_2_1.fq.gz hic_mut_R1.fastq.gz
mv hic_mut_2_2.fq.gz hic_mut_R2.fastq.gz
mv hic_oe_smxl7_2_1.fq.gz hic_oe_R1.fastq.gz
mv hic_oe_smxl7_2_2.fq.gz hic_oe_R2.fastq.gz
mv hic_WT_2_1.fq.gz hic_wt_R1.fastq.gz 
mv hic_WT_2_2.fq.gz hic_wt_R2.fastq.gz 
mkdir -p mydata/mut
mkdir -p mydata/oe
mkdir -p mydata/wt
cp hic_mut* mydata/mut
cp hic_oe* mydata/oe
cp hic_wt* mydata/wt

```



### Reads were mapped tair genome with HiC-pro software.

### How can I generate the list of restriction fragments after genome digestion ?
```sh
# site https://github.com/nservant/HiC-Pro
# https://nservant.github.io/HiC-Pro/UTILS.html#utils
/public/home/2022122/xugang/app/HiC-Pro/HiC-Pro-master/bin/utils/digest_genome.py -r hindiii dpnii -o /public/home/2022122/xugang/project/tair10/tair10.hindiii_dpnii.bed /public/home/2022122/xugang/project/tair10/tair10.fa 
```
### Generate the chromosome's file
```sh
cut -f 1,2 /public/home/2022122/xugang/project/tair10/tair10.fa.fai > /public/home/2022122/xugang/project/tair10/tair10.chromesize
# The output file
/public/home/2022122/xugang/project/tair10/tair10.chromesize
```

### Run HiC-Pro
```sh
mkdir -p /public/home/2022122/xugang/project/yaoruifeng/output/a10.hicseq/log

```
### Write the below content into the /public/home/2022122/xugang/project/yaoruifeng/hic.txt
```sh
# Please change the variable settings below if necessary

#########################################################################
## Paths and Settings  - Do not edit !
#########################################################################

TMP_DIR = tmp
LOGS_DIR = logs
BOWTIE2_OUTPUT_DIR = bowtie_results
MAPC_OUTPUT = hic_results
RAW_DIR = rawdata

#######################################################################
## SYSTEM AND SCHEDULER - Start Editing Here !!
#######################################################################
N_CPU = 50
SORT_RAM = 50000M
LOGFILE = hicpro.log

JOB_NAME =
JOB_MEM =
JOB_WALLTIME =
JOB_QUEUE =
JOB_MAIL =

#########################################################################
## Data
#########################################################################
PAIR1_EXT = _R1
PAIR2_EXT = _R2

#######################################################################
## Alignment options
#######################################################################

MIN_MAPQ = 10

BOWTIE2_IDX_PATH = /public/home/2022122/xugang/project/tair10/
BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder
BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder

#######################################################################
## Annotation files
#######################################################################

REFERENCE_GENOME = tair10.bowtie2
GENOME_SIZE = /public/home/2022122/xugang/project/tair10/tair10.chromesize

#######################################################################
## Allele specific analysis
#######################################################################

ALLELE_SPECIFIC_SNP =

#######################################################################
## Capture Hi-C analysis
#######################################################################

CAPTURE_TARGET =
REPORT_CAPTURE_REPORTER = 1

#######################################################################
## Digestion Hi-C
#######################################################################

GENOME_FRAGMENT = /public/home/2022122/xugang/project/tair10/tair10.hindiii_dpnii.bed
LIGATION_SITE = AAGCTAGCTT
MIN_FRAG_SIZE =
MAX_FRAG_SIZE =
MIN_INSERT_SIZE =
MAX_INSERT_SIZE =

#######################################################################
## Hi-C processing
#######################################################################

MIN_CIS_DIST =
GET_ALL_INTERACTION_CLASSES = 1
GET_PROCESS_SAM = 0
RM_SINGLETON = 1
RM_MULTI = 1
RM_DUP = 1

#######################################################################
## Contact Maps
#######################################################################

BIN_SIZE = 20000 40000 150000 500000 1000000
MATRIX_FORMAT = upper

#######################################################################
## Normalization
#######################################################################
MAX_ITER = 100
FILTER_LOW_COUNT_PERC = 0.02
FILTER_HIGH_COUNT_PERC = 0
EPS = 0.1
```

```sh
#conda activate /public/home/2022122/xugang/app/HiC-Pro-master/hicpro
#plot peak heatmap.
hicf(){
log=$output/a10.hicseq/log
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
rm -rf /public/home/2022122/xugang/project/yaoruifeng/output/a10.hicseq/Hi-C
/public/home/2022122/xugang/app/HiC-Pro-master/bin/HiC-Pro -i rawdata/hic-seq/mydata -o ${output}/a10.hicseq/Hi-C -c hic.txt
">a10.hic.$counter.$name1.sh

}
hicf

```

###  FAN-C (v.0.9.8).
visit https://support.hdfgroup.org/ftp/HDF5/current/src/   
download hdf5-1.10.5.tar.gz
```sh
# make a new directory
mkdir hdf5-build
cd hdf5-build
# replace xx with current version number
wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.xx.tar.gz
# unpack
tar xzf hdf5-1.8.xx.tar.gz
cd hdf5-1.8.xx/
# use --prefix to set the folder in which HDF5 should be installed
# alternatively, you can omit --prefix=... here and run
# sudo make install to install globally (requires admin rights)
./configure --prefix=/public/home/2022122/xugang/app/hdf5-1.8.21
make
make install
export HDF5_DIR=/public/home/2022122/xugang/app/hdf5-1.8.21
```

### FAN-C  

* The bam files (bwt2merged.bam) generated by HiC-Pro containing mapped reads were used as input files for FAN-C (v.0.9.8). 
* The module ‘fanc auto’ was applied to generate 500, 100, 50, 10 and 1 kb contact matrices (hic files). 
* The resultant hic files with 100-kb resolution were directed to the ‘fanc expected’ module to calculate the expected interaction probability against genomic distance for intrachromosomal interaction. 
* For matrix and score comparisons, the default comparison method of fold-change was used with the ‘fanc compare’ command. 
* The outputs (hic object) were transferred to text files by ‘fanc dump’ and were visualized as heatmaps in R using ggplot2. 
* To explore whether higher contacts observed in p35S::H2B.8–eGFP depend on H2B.8 incorporation, the genome was binned into 1-kb windows. 
* H2B.8 signals (log2(IP/input)) of each window were generated and sorted into 20 quantiles by strength. 
* The interaction frequency differences (values generated by FAN-C at 1-kb resolution) of each quantile pair for either short-range interactions or interactions between pericentromeric regions and chromosome arms were averaged and plotted as heatmaps. 
* The resolution of our Hi-C data was estimated as previously reported. 
* Our Hi-C data was deemed to achieve 1-kb resolution as 80% of genomic bins (1 kb) had >1,000 contacts. Our WT data were compared to published contact matrices.

We strongly recommend installing FAN-C in a fresh virtual environment to prevent dependency issues. You can do this with tools like pyenv or manually using venv, for example:  

```sh
#下载
#1.1 创建新环境
conda create -n FAN-C python=3.8 -y
conda activate FAN-C
#1.2 pip 下载
conda install libgcc=5.2.0
pip install fanc
pip install --upgrade --force-reinstall fanc
#note，如果出现报错可以尝试以下办法
pip uninstall fanc
pip uninstall pysam
pip install Cython
pip install pysam
pip install fanc
# 2.测试
# 2.1 下载测试数据
wget -O examples.zip "https://keeper.mpdl.mpg.de/d/147906745b634c779ed3/files/?p=/examples.zip&dl=1"
unzip examples.zip
cd examples
#
source /public/home/2022122/xugang/bashrc
cd /public/home/2022122/xugang/app/fanc_examples

conda run -n FAN-C fanc auto SRR4271982_chr18_19_1.fastq.gzip SRR4271982_chr18_19_2.fastq.gzip output/ \
          -g hg19_chr18_19.fa -i hg19_chr18_19/hg19_chr18_19 -n fanc_example -t 4 -r HindIII \
          --split-ligation-junction -q 30

#2.2 测试软件
fanc auto SRR4271982_chr18_19_1.fastq.gzip SRR4271982_chr18_19_2.fastq.gzip output/ \
          -g hg19_chr18_19.fa -i hg19_chr18_19/hg19_chr18_19 -n fanc_example -t 4 -r HindIII \
          --split-ligation-junction -q 30 --run-with test

2023-12-12 14:51:25,903 INFO FAN-C version: 0.9.27
2023-12-12 14:51:25,994 INFO Output folder: output/
2023-12-12 14:51:25,994 INFO Input files: SRR4271982_chr18_19_1.fastq.gzip, SRR4271982_chr18_19_2.fastq.gzip
2023-12-12 14:51:25,994 INFO Input file types: fastq, fastq
2023-12-12 14:51:25,994 INFO Final basename: fanc_example (you can change this with the -n option!)
0 (threads: 4, depends on: ): fanc map -m 25 -s 3 -t 4 --no-iterative --restriction-enzyme HindIII SRR4271982_chr18_19_1.fastq.gzip hg19_chr18_19/hg19_chr18_19 output/sam/SRR4271982_chr18_19_1.bam
1 (threads: 4, depends on: ): fanc map -m 25 -s 3 -t 4 --no-iterative --restriction-enzyme HindIII SRR4271982_chr18_19_2.fastq.gzip hg19_chr18_19/hg19_chr18_19 output/sam/SRR4271982_chr18_19_2.bam
2 (threads: 1, depends on: 0, 1): fanc sort_sam -t 4 --no-sambamba output/sam/SRR4271982_chr18_19_1.bam
3 (threads: 1, depends on: 0, 1): fanc sort_sam -t 4 --no-sambamba output/sam/SRR4271982_chr18_19_2.bam
4 (threads: 4, depends on: 2, 3): fanc pairs -f -g hg19_chr18_19.fa -t 4 -us -r HindIII -q 30.0 -S output/sam/SRR4271982_chr18_19_1.bam output/sam/SRR4271982_chr18_19_2.bam output/pairs/fanc_example.pairs
5 (threads: 1, depends on: 4): fanc pairs -d 10000 -l -p 2 --statistics-plot output/plots/stats/fanc_example.pairs.stats.pdf output/pairs/fanc_example.pairs
6 (threads: 1, depends on: 5): fanc pairs --ligation-error-plot output/plots/stats/fanc_example.pairs.ligation_error.pdf output/pairs/fanc_example.pairs
7 (threads: 1, depends on: 5): fanc pairs --re-dist-plot output/plots/stats/fanc_example.pairs.re_dist.pdf output/pairs/fanc_example.pairs
8 (threads: 1, depends on: 5): fanc hic -f output/pairs/fanc_example.pairs output/hic/fanc_example.hic
9 (threads: 4, depends on: 8): fanc hic -f -b 5000000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_5mb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_5mb.hic
10 (threads: 4, depends on: 8): fanc hic -f -b 2000000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_2mb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_2mb.hic
11 (threads: 4, depends on: 8): fanc hic -f -b 1000000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_1mb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_1mb.hic
12 (threads: 4, depends on: 8): fanc hic -f -b 500000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_500kb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_500kb.hic
13 (threads: 4, depends on: 8): fanc hic -f -b 250000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_250kb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_250kb.hic
14 (threads: 4, depends on: 8): fanc hic -f -b 100000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_100kb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_100kb.hic
15 (threads: 4, depends on: 8): fanc hic -f -b 50000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_50kb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_50kb.hic
16 (threads: 4, depends on: 8): fanc hic -f -b 25000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_25kb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_25kb.hic
17 (threads: 4, depends on: 8): fanc hic -f -b 10000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_10kb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_10kb.hic
18 (threads: 4, depends on: 8): fanc hic -f -b 5000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_5kb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_5kb.hic


```

```sh
usage: fanc <command> [options]

-- Matrix generation --
auto              Automatically process an entire Hi-C data set
map               Map reads in a FASTQ file to a reference genome
pairs             Process and filter read pairs
hic               Process, filter, and correct Hic files

-- Matrix analysis --
cis-trans         Calculate cis/trans ratio of this Hi-C object
expected          Calculate Hi-C expected values (distance decay)
pca               Do a PCA on multiple Hi-C objects
compartments      Calculate AB compartment matrix
insulation        Calculate insulation scores for Hic object
directionality    Calculate directionality index for Hic object
boundaries        Determine domain boundaries
compare           Create pairwise comparisons of Hi-C comparison maps
loops             Call loops in a Hic object using FAN-C implementation of HICCUPS
aggregate         Make aggregate plots with FAN-C

-- Other helpers --
fragments         In-silico genome digestion
sort-sam          Convenience function to sort a SAM file by name
from-juicer       Import a Hi-C object from juicer (Aiden lab)
from-txt          Import a Hi-C object from a sparse matrix txt format
to-cooler         Convert a Hic file into cooler format
to-juicer         Convert a ReadPairs file to Juicer 
dump              Dump Hic file to txt file(s)
overlap-peaks     Overlap peaks from multiple samples
subset            Create a new Hic object by subsetting
stats             Get statistics on number of reads used at each step of a pipeline
write-config      Write default config file to specified location
downsample        Downsample contacts from a Hic object
upgrade           Upgrade objects from old FAN-C versions

#Automatically process an entire Hi-C data set
fanc auto SRR4271982_chr18_19_1.fastq.gzip SRR4271982_chr18_19_2.fastq.gzip output/ \
          -g hg19_chr18_19.fa -i hg19_chr18_19/hg19_chr18_19 -n fanc_example -t 4 -r HindIII \
          --split-ligation-junction -q 30 --run-with test
# Map reads in a FASTQ file to a reference genome
fanc map -m 25 -s 3 -t 4 --no-iterative --restriction-enzyme HindIII SRR4271982_chr18_19_1.fastq.gzip hg19_chr18_19/hg19_chr18_19 output/sam/SRR4271982_chr18_19_1.bam
fanc map -m 25 -s 3 -t 4 --no-iterative --restriction-enzyme HindIII SRR4271982_chr18_19_2.fastq.gzip hg19_chr18_19/hg19_chr18_19 output/sam/SRR4271982_chr18_19_2.bam
#Convenience function to sort a SAM file by name
fanc sort_sam -t 4 --no-sambamba output/sam/SRR4271982_chr18_19_1.bam
fanc sort_sam -t 4 --no-sambamba output/sam/SRR4271982_chr18_19_2.bam
#Process and filter read pairs
fanc pairs -f -g hg19_chr18_19.fa -t 4 -us -r HindIII -q 30.0 -S output/sam/SRR4271982_chr18_19_1.bam output/sam/SRR4271982_chr18_19_2.bam output/pairs/fanc_example.pairs
fanc pairs -d 10000 -l -p 2 --statistics-plot output/plots/stats/fanc_example.pairs.stats.pdf output/pairs/fanc_example.pairs
fanc pairs --ligation-error-plot output/plots/stats/fanc_example.pairs.ligation_error.pdf output/pairs/fanc_example.pairs
fanc pairs --re-dist-plot output/plots/stats/fanc_example.pairs.re_dist.pdf output/pairs/fanc_example.pairs
#Process, filter, and correct Hic files
fanc hic -f output/pairs/fanc_example.pairs output/hic/fanc_example.hic
fanc hic -f -b 5000000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_5mb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_5mb.hic
fanc hic -f -b 2000000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_2mb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_2mb.hic
fanc hic -f -b 1000000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_1mb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_1mb.hic
fanc hic -f -b 500000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_500kb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_500kb.hic
fanc hic -f -b 250000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_250kb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_250kb.hic
fanc hic -f -b 100000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_100kb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_100kb.hic
fanc hic -f -b 50000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_50kb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_50kb.hic
fanc hic -f -b 25000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_25kb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_25kb.hic
fanc hic -f -b 10000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_10kb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_10kb.hic
fanc hic -f -b 5000 -r 0.1 -t 4 --statistics-plot output/plots/stats/fanc_example_5kb.stats.pdf -n --norm-method kr output/hic/fanc_example.hic output/hic/binned/fanc_example_5kb.hic
#
fanc expected -p output/expected/fanc_example_100kb_expected.png  -c chr19  output/hic/binned/fanc_example_100kb.hic  output/expected/fanc_example_100kb_expected.txt
#fanc expected -l "HindIII 100k" "HindIII 5M" "MboI 100k" "MboI 1M" "MboI 50k"  -c chr19 -p architecture/expected/expected_multi.png  architecture/other-hic/lowc_hindiii_100k_1mb.hic  architecture/other-hic/lowc_hindiii_5M_1mb.hic architecture/other-hic/lowc_mboi_100k_1mb.hic  architecture/other-hic/lowc_mboi_1M_1mb.hic  architecture/other-hic/lowc_mboi_50k_1mb.hic architecture/expected/expected_multi.txt
fancplot -o output/expected/fanc_example_100kb_chr18_oe.png  chr18:1-78mb -p triangular -e output/hic/binned/fanc_example_100kb.hic  -vmin -2 -vmax 2

```

```sh
#!/bin/bash
#SBATCH -o /public/home/2022122/xugang/app/fanc_examples/log/log.%j.out
#SBATCH -e /public/home/2022122/xugang/app/fanc_examples/log/log.%j.error
#SBATCH --partition=Cnode2
#SBATCH -J 5fanc
#SBATCH -N 1
#SBATCH -n 50
source /public/home/2022122/xugang/bashrc
#run test
# Map reads in a FASTQ file to a reference genome
echo "conda run -n FAN-C fanc map -m 25 -s 3 -t 50 --no-iterative --restriction-enzyme HindIII SRR4271982_chr18_19_1.fastq.gzip hg19_chr18_19/hg19_chr18_19 output2/sam/SRR4271982_chr18_19_1.bam"
conda run -n FAN-C fanc map -m 25 -s 3 -t 50 --no-iterative --restriction-enzyme HindIII SRR4271982_chr18_19_1.fastq.gzip hg19_chr18_19/hg19_chr18_19 output2/sam/SRR4271982_chr18_19_1.bam
echo "conda run -n FAN-C fanc map -m 25 -s 3 -t 50 --no-iterative --restriction-enzyme HindIII SRR4271982_chr18_19_2.fastq.gzip hg19_chr18_19/hg19_chr18_19 output2/sam/SRR4271982_chr18_19_2.bam"
conda run -n FAN-C fanc map -m 25 -s 3 -t 50 --no-iterative --restriction-enzyme HindIII SRR4271982_chr18_19_2.fastq.gzip hg19_chr18_19/hg19_chr18_19 output2/sam/SRR4271982_chr18_19_2.bam

#Convenience function to sort a SAM file by name
echo "conda run -n FAN-C fanc sort_sam -t 50 --no-sambamba output2/sam/SRR4271982_chr18_19_1.bam"
conda run -n FAN-C fanc sort_sam -t 50 --no-sambamba output2/sam/SRR4271982_chr18_19_1.bam
echo "conda run -n FAN-C fanc sort_sam -t 50 --no-sambamba output2/sam/SRR4271982_chr18_19_2.bam"
conda run -n FAN-C fanc sort_sam -t 50 --no-sambamba output2/sam/SRR4271982_chr18_19_2.bam
#Process and filter read pairs
echo "conda run -n FAN-C fanc pairs -f -g hg19_chr18_19.fa -t 50 -us -r HindIII -q 30.0 -S output2/sam/SRR4271982_chr18_19_1.bam output2/sam/SRR4271982_chr18_19_2.bam output2/pairs/fanc_example.pairs"
conda run -n FAN-C fanc pairs -f -g hg19_chr18_19.fa -t 50 -us -r HindIII -q 30.0 -S output2/sam/SRR4271982_chr18_19_1.bam output2/sam/SRR4271982_chr18_19_2.bam output2/pairs/fanc_example.pairs

echo "conda run -n FAN-C fanc pairs -d 10000 -l -p 2 --statistics-plot output2/plots/stats/fanc_example.pairs.stats.pdf output2/pairs/fanc_example.pairs"
conda run -n FAN-C fanc pairs -d 10000 -l -p 2 --statistics-plot output2/plots/stats/fanc_example.pairs.stats.pdf output2/pairs/fanc_example.pairs

echo "conda run -n FAN-C fanc pairs --ligation-error-plot output2/plots/stats/fanc_example.pairs.ligation_error.pdf output2/pairs/fanc_example.pairs"
conda run -n FAN-C fanc pairs --ligation-error-plot output2/plots/stats/fanc_example.pairs.ligation_error.pdf output2/pairs/fanc_example.pairs

echo "conda run -n FAN-C fanc pairs --re-dist-plot output2/plots/stats/fanc_example.pairs.re_dist.pdf output2/pairs/fanc_example.pairs"
conda run -n FAN-C fanc pairs --re-dist-plot output2/plots/stats/fanc_example.pairs.re_dist.pdf output2/pairs/fanc_example.pairs

#Process, filter, and correct Hic files
echo "conda run -n FAN-C fanc hic -f output2/pairs/fanc_example.pairs output2/hic/fanc_example.hic"
conda run -n FAN-C fanc hic -f output2/pairs/fanc_example.pairs output2/hic/fanc_example.hic
#
echo "conda run -n FAN-C fanc hic -f -b 5000000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_5mb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_5mb.hic"
conda run -n FAN-C fanc hic -f -b 5000000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_5mb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_5mb.hic
#
echo "conda run -n FAN-C fanc hic -f -b 2000000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_2mb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_2mb.hic"
conda run -n FAN-C fanc hic -f -b 2000000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_2mb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_2mb.hic

echo "conda run -n FAN-C fanc hic -f -b 1000000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_1mb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_1mb.hic"
conda run -n FAN-C fanc hic -f -b 1000000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_1mb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_1mb.hic

echo "conda run -n FAN-C fanc hic -f -b 500000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_500kb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_500kb.hic"
conda run -n FAN-C fanc hic -f -b 500000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_500kb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_500kb.hic

echo "conda run -n FAN-C fanc hic -f -b 250000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_250kb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_250kb.hic"
conda run -n FAN-C fanc hic -f -b 250000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_250kb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_250kb.hic

echo "conda run -n FAN-C fanc hic -f -b 100000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_100kb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_100kb.hic"
conda run -n FAN-C fanc hic -f -b 100000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_100kb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_100kb.hic

echo "conda run -n FAN-C fanc hic -f -b 50000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_50kb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_50kb.hic"
conda run -n FAN-C fanc hic -f -b 50000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_50kb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_50kb.hic

echo "conda run -n FAN-C fanc hic -f -b 25000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_25kb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_25kb.hic"
conda run -n FAN-C fanc hic -f -b 25000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_25kb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_25kb.hic

echo "conda run -n FAN-C fanc hic -f -b 10000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_10kb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_10kb.hic"
conda run -n FAN-C fanc hic -f -b 10000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_10kb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_10kb.hic

echo "conda run -n FAN-C fanc hic -f -b 5000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_5kb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_5kb.hic"
conda run -n FAN-C fanc hic -f -b 5000 -r 0.1 -t 50 --statistics-plot output2/plots/stats/fanc_example_5kb.stats.pdf -n --norm-method kr output2/hic/fanc_example.hic output2/hic/binned/fanc_example_5kb.hic

```
### Run the steps with FAN-C
```sh
fancf(){
log=$output/a11.fanc/log
[[ -d $log ]] || mkdir -p $log
[[ -d $output/a11.fanc/bam ]] || mkdir -p $output/a11.fanc/bam
[[ -d $output/a11.fanc/pairs ]] ||  mkdir -p $output/a11.fanc/pairs
[[ -d $output/a11.fanc/plots/stats/ ]] || mkdir -p $output/a11.fanc/plots/stats/
[[ -d $output/a11.fanc/hic ]] || mkdir -p $output/a11.fanc/hic
[[ -d $output/a11.fanc/hic/binned ]] || mkdir -p $output/a11.fanc/hic/binned

#remove background file
((counter++))
name=$1
name1=$2
name2=$3


echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 5${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc

cp ${output}/a10.hicseq/Hi-C/bowtie_results/bwt2/${name}/${name1} $output/a11.fanc/bam/${name}.R1.bam
cp ${output}/a10.hicseq/Hi-C/bowtie_results/bwt2/${name}/${name2} $output/a11.fanc/bam/${name}.R2.bam
#Convenience function to sort a SAM file by name
#sort
conda run -n FAN-C fanc sort_sam -t ${thread} --no-sambamba $output/a11.fanc/bam/${name}.R1.bam $output/a11.fanc/bam/${name}.sort.R1.bam
conda run -n FAN-C fanc sort_sam -t ${thread} --no-sambamba $output/a11.fanc/bam/${name}.R2.bam $output/a11.fanc/bam/${name}.sort.R2.bam
#pairs
conda run -n FAN-C fanc pairs -f -g /public/home/2022122/xugang/project/tair10/tair10.fa  -t ${thread} -us -r HindIII,MboI -q 30.0 -S $output/a11.fanc/bam/${name}.sort.R1.bam $output/a11.fanc/bam/${name}.sort.R2.bam $output/a11.fanc/pairs/${name}.pairs
#hic
conda run -n FAN-C fanc hic -f $output/a11.fanc/pairs/${name}.pairs $output/a11.fanc/hic/${name}.hic
conda run -n FAN-C fanc hic -f -b 100000 -r 0.1 -t ${thread} --statistics-plot $output/a11.fanc/plots/stats/${name}_100kb.stats.pdf -n --norm-method kr $output/a11.fanc/hic/${name}.hic $output/a11.fanc/hic/binned/fa${name}_100kb.hic
#plot
conda run -n FAN-C fanc pairs -d 10000 -l -p 2 --statistics-plot $output/a11.fanc/plots/stats/${name}.pairs.stats.pdf $output/a11.fanc/pairs/${name}.pairs

conda run -n FAN-C fanc pairs --ligation-error-plot $output/a11.fanc/plots/stats/${name}.pairs.ligation_error.pdf $output/a11.fanc/pairs/${name}.pairs

conda run -n FAN-C fanc pairs --re-dist-plot $output/a11.fanc/plots/stats/${name}.pairs.re_dist.pdf $output/a11.fanc/pairs/${name}.pairs


">a11.fanc.$counter.$name1.sh
}
fancf mut hic_mut_R1_tair10.bowtie2.bwt2merged.bam hic_mut_R2_tair10.bowtie2.bwt2merged.bam
fancf oe  hic_oe_R1_tair10.bowtie2.bwt2merged.bam  hic_oe_R2_tair10.bowtie2.bwt2merged.bam
fancf wt  hic_wt_R1_tair10.bowtie2.bwt2merged.bam  hic_wt_R2_tair10.bowtie2.bwt2merged.bam


fancdumpf(){
log=$output/a11.fanc/log
#remove background file
((counter++))
name=$1
[[ -d $output/a11.fanc/dump ]] || mkdir -p $output/a11.fanc/dump
[[ -d $output/a11.fanc/expected ]] || mkdir -p $output/a11.fanc/expected

echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 5${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc

conda run -n FAN-C fanc dump $output/a11.fanc/hic/binned/fa${name}_100kb.hic $output/a11.fanc/dump/${name}_100kb.dum.tsv
conda run -n FAN-C fanc expected -p $output/a11.fanc/expected/${name}_100kb_expected.png  $output/a11.fanc/hic/binned/fa${name}_100kb.hic  $output/a11.fanc/expected/${name}_100kb_expected.txt
">a12.fanc.dump.$counter.$name.sh
}
fancdumpf mut
fancdumpf oe
fancdumpf wt


fanccomparef(){
log=$output/a11.fanc/log
#remove background file
((counter++))
name1=$1
name2=$2
[[ -d $output/a11.fanc/compare ]] || mkdir -p $output/a11.fanc/compare

echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 5${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc

conda run -n FAN-C fanc compare $output/a11.fanc/hic/binned/fa${name1}_100kb.hic $output/a11.fanc/hic/binned/fa${name2}_100kb.hic   $output/a11.fanc/compare/${name1}.vs.${name2}.comparsion.hic
conda run -n FAN-C fanc dump $output/a11.fanc/compare/${name1}.vs.${name2}.comparsion.hic $output/a11.fanc/compare/${name1}.vs.${name2}.comparsion.tsv
">a13.fanc.compare.$counter.$name1.sh
}
fanccomparef mut wt
fanccomparef oe wt

```


#### FANC help.
fanc auto
```sh
# fanc auto -h
2023-12-13 11:12:46,677 INFO FAN-C version: 0.9.27
usage: fanc auto [-h] [-g GENOME] [-r RESTRICTION_ENZYME] [-i GENOME_INDEX] [-n BASENAME] [-s STEP_SIZE] [-b BIN_SIZES [BIN_SIZES ...]] [-t THREADS]
                 [--max-restriction-site-distance MAX_RESTRICTION_SITE_DISTANCE] [--fanc-parallel] [--split-fastq] [--memory-map] [--ice] [--norm-method NORM_METHOD] [-q QUALITY_CUTOFF]
                 [--iterative-quality-cutoff ITERATIVE_QUALITY_CUTOFF] [--le-inward-cutoff INWARD_CUTOFF] [--le-outward-cutoff OUTWARD_CUTOFF] [--auto-le-cutoff] [-tmp] [--iterative] [--no-sam-sort]
                 [--restore-coverage] [--split-ligation-junction] [--no-filter-pairs] [--no-hic] [--run-with RUN_WITH] [--job-prefix JOB_PREFIX] [--grid-startup-commands GRID_STARTUP_COMMANDS]
                 [--grid-cleanup-commands GRID_CLEANUP_COMMANDS] [-f] [--no-sambamba USE_SAMBAMBA]
                 input [input ...] output_folder

Automatically process an entire Hi-C data set.

positional arguments:
  input                 Input files. fanc will try to guess the file type by its extension.
  output_folder         Output folder. All output files and folders will be generated under this directory.

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        Genome for the Hi-C object.Path to region-based file (BED, GFF, ...) containing the non-overlapping regions to be used for Hi-C object construction. Typically restriction-enzyme
                        fragments. Alternatively: Path to genome file (FASTA, folder with FASTA, hdf5 file), which will be used in conjunction with the type of restriction enzyme (-r) to calculate fragments
                        directly.
  -r RESTRICTION_ENZYME, --restriction-enzyme RESTRICTION_ENZYME
                        Restriction enzyme name. Used for in silico digestion of genomic sequences and splitting of reads at Hi-C ligation junctions. (e.g. HindIII, case-sensitive). Separate multiple enzymes
                        with ','. Restriction names can be any supported by Biopython, which obtains data from REBASE (http://rebase.neb.com/rebase/rebase.html).
  -i GENOME_INDEX, --genome-index GENOME_INDEX
                        Bowtie 2 or BWA genome index. Only required when passing FASTQ files as input.
  -n BASENAME, --basename BASENAME
                        Basename for output files. If not provided, will be guessed based on input file names.
  -s STEP_SIZE, --step-size STEP_SIZE
                        Step size for iterative mapping. Default: 3
  -b BIN_SIZES [BIN_SIZES ...], --bin-sizes BIN_SIZES [BIN_SIZES ...]
                        Bin sizes for Hi-C matrix generation. Default: 5mb, 2mb, 1mb, 500kb, 250kb, 100kb, 50kb, 25kb, 10kb, 5kb.
  -t THREADS, --threads THREADS
                        Maximum number of threads. The number provided here will not be exceeded by analysis steps running alone or in parallel. Default: 1
  --max-restriction-site-distance MAX_RESTRICTION_SITE_DISTANCE
                        Insert / ligation fragment sizes are inferred from the sum of distances of both reads to the nearest restriction sites. Fragment larger than this value are filtered from the Pairs
                        object. The default value of 10000 only removes extremely large fragments, and is thus considered conservative.Default: 10000
  --fanc-parallel       Use FAN-C parallelisation, which launches multiple mapper jobs. This may be faster in some cases than relying on the internal parallelisation of the mapper, but has potentially high
                        disk I/O and memory usage.
  --split-fastq         Split fastq files into chunks of 10M reads. Reads will be merged again on the SAM level. Splitting and merging bypasses the -tmp flag. This option reduces disk usage in tmp, in case
                        the system has a small tmp partition.
  --memory-map          Map Bowtie2 index to memory. Recommended if running on medium-memory systems and using many parallel threads).
  --ice                 DEPRECATED. Correct Hi-C matrices using ICE instead of Knight-Ruiz matrix balancing. Slower, but much more memory-friendly.
  --norm-method NORM_METHOD
                        Normalisation method. Options are: KR (default) = Knight-Ruiz matrix balancing (Fast, accurate, but memory-intensive); ICE = ICE matrix balancing (more CPU-intensive, but also more
                        memory-efficient); VC = vanilla coverage (a single round of ICE balancing); VC-SQRT = vanilla coverage square root (reduces overcorrection compared to VC)
  -q QUALITY_CUTOFF, --quality-cutoff QUALITY_CUTOFF
                        Cutoff for the minimum mapping quality of a read. For numbers larger than 1, will filter on MAPQ. If a number between 0 and 1 is provided, will filter on the AS tag instead of mapping
                        quality (only BWA). The quality cutoff is then interpreted as the fraction of bases that have to be matched for any given read. Only applies to SAM/BAM input!. Default is not to
                        filter on mapping quality.
  --iterative-quality-cutoff ITERATIVE_QUALITY_CUTOFF
                        MAPQ cutoff for mapped reads. Only applies when iterative mapping is enabled: if a mapped read has MAPQ below this cutoff,it will be sent to another iteration in an attempt to find a
                        higher quality alignment. Default is 3 for BWA and 30 for Bowtie2.
  --le-inward-cutoff INWARD_CUTOFF
                        Ligation error inward cutoff. Default: no ligation error filtering.
  --le-outward-cutoff OUTWARD_CUTOFF
                        Ligation error outward cutoff. Default: no ligation error filtering.
  --auto-le-cutoff      Automatically determine ligation error cutoffs. Use with caution, this setting has a tendency to choose large cutoffs, removing many pairs close to the diagonal.
  -tmp, --work-in-tmp   Work in temporary directory. Copies input files and generates output files in a temporary directory. Files will be moved to their intended destination once an analysis step finishes.
                        Reduces network I/O if using remote file systems.
  --iterative           Map reads iteratively. Can improve mappability, especially with low-quality reads. Reads are initially trimmed to 25bp and mapped to the reference genome. If no unique mapping
                        location is found, the read is extended by 3bp and the process is repeated until the full length of the read is reached or a unique mapping location is found.
  --no-sam-sort         Do not sort SAM/BAM files. Sorted files are required for the pair generating step. Only omit this if you are supplying presorted (by read name) SAM/BAM files.
  --restore-coverage    Restore coverage to the original total number of reads. Otherwise matrix entries will be contact probabilities. Only available for KR matrix balancing.
  --split-ligation-junction
                        Split reads at predicted ligation junction before mapping. Requires the -r argument.
  --no-filter-pairs     Do not filter read pairs. By default, the following filters are applied: self-ligations, PCR duplicates,restriction distance (>10kb)
  --no-hic              Do not process pairs into Hi-C maps (stop after read pairing step).
  --run-with RUN_WITH   Choose how to run the commands in fanc auto. Options: 'parallel' (default): Run fanc commands on local machine, use multiprocessing parallelisation. 'sge': Submit fanc commands to a
                        Sun/Oracle Grid Engine cluster. 'slurm': Submit fanc commands to a Slurm cluster. 'test': Do not run fanc commands but print all commands and their dependencies to stdout for review.
  --job-prefix JOB_PREFIX
                        Job Prefix for SGE and Slurm. Works with '--run-with sge' and --run-with slurm. Default: 'fanc_<6 random letters>_'
  --grid-startup-commands GRID_STARTUP_COMMANDS
                        Path to a file with BASH commands that are executed before every FAN-C command that is run on a grid engine / cluster. This could, for example, include environment-specific settings,
                        such as activation of a Python virtualenv.
  --grid-cleanup-commands GRID_CLEANUP_COMMANDS
                        Path to a file with BASH commands that are executed after every FAN-C command that is run on a grid engine / cluster. Use this to clean the file system or environment set up with
                        --grid-startup-commands.
  -f, --force-overwrite
                        Force overwriting of existing files. Otherwise you will be prompted before files are overwritten.
  --no-sambamba USE_SAMBAMBA
                        Do not use sambamba for sorting SAM file, even when it is available. Use pysam instead.
(FAN-C) [2022122@admin1 fanc_examples]$
```

fanc map
```sh
# fanc map
usage: fanc map [-h] [-m MIN_SIZE] [-s STEP_SIZE] [--trim-front] [-t THREADS] [-q QUALITY] [-r RESTRICTION_ENZYME] [-k MAX_ALIGNMENTS] [-a] [-b BATCH_SIZE] [--fanc-parallel] [--split-fastq] [--memory-map]
                [--no-iterative] [--mapper-type MAPPER_TYPE] [-tmp]
                input [input ...] index output

Map reads in a FASTQ file to a reference genome.

positional arguments:
  input                 File name of the input FASTQ file (or gzipped FASTQ)
  index                 Bowtie 2 or BWA genome index base. Index type will be determined automatically.
  output                Output file or folder. When providing multiple input files, this must be the path to an output folder.

optional arguments:
  -h, --help            show this help message and exit
  -m MIN_SIZE, --min-size MIN_SIZE
                        Minimum length of read before extension. Default 25.
  -s STEP_SIZE, --step-size STEP_SIZE
                        Number of base pairs to extend at each round of mapping. Default is 10.
  --trim-front          Trim reads from front instead of back.
  -t THREADS, --threads THREADS
                        Number of threads used for mapping. Default: 1
  -q QUALITY, --quality QUALITY
                        Mapping quality cutoff. Alignments with a quality score lower than this will be sent to another mapping iteration. Default: 3 (BWA), 30 (Bowtie2)
  -r RESTRICTION_ENZYME, --restriction-enzyme RESTRICTION_ENZYME
                        Name (case sensitive) of restriction enzyme used in Hi-C experiment. Will be used to split reads by predicted ligation junction before mapping. You can omit this if you do not want to
                        split your reads by ligation junction. Restriction names can be any supported by Biopython, which obtains data from REBASE (http://rebase.neb.com/rebase/rebase.html). For restriction
                        enzyme cocktails, separate enzyme names with ","
  -k MAX_ALIGNMENTS, --max-alignments MAX_ALIGNMENTS
                        Maximum number of alignments per read to be reported.
  -a, --all-alignments  Report all valid alignments of a read Warning: very slow!.
  -b BATCH_SIZE, --batch-size BATCH_SIZE
                        Number of reads processed (mapped and merged) in one go per worker. The default 100000 works well for large indexes (e.g. human, mouse). Smaller indexes (e.g. yeast) will finish
                        individual bowtie2 processes very quickly - set this number higher to spawn new processes less frequently.
  --fanc-parallel       Use FAN-C parallelisation, which launches multiple mapper jobs. This may be faster in some cases than relying on the internal paralellisation of the mapper, but has potentially high
                        disk I/O and memory usage.
  --split-fastq         Split FASTQ file into 10M chunks before mapping. Easier on tmp partitions.
  --memory-map          Map Bowtie2 index to memory. Only enable if your system has enough memory to hold the entire Bowtie2 index.
  --no-iterative        Do not use iterative mapping strategy. (much faster, less sensitive).
  --mapper-type MAPPER_TYPE
                        Manually set mapper type. Currently supported: bowtie2, bwa, and bwa-mem2. Not that this is generally auto-detected from the index path.
  -tmp, --work-in-tmp   Copy original file to temporary directory.Reduces network I/O.
```
fanc sort_sam
```sh
#fanc sort_sam
usage: fanc sort_sam [-h] [-t THREADS] [-S] [-tmp] sam [output]

Convenience function to sort a SAM file by name. Exactly the same as 'samtools sort -n', but potentiallyfaster if sambamba is available.

positional arguments:
  sam                   Input SAM/BAM
  output                Output SAM/BAM. If not provided, will replace input file with sorted version after sorting.

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of sorting threads (only when sambamba is available). Default: 1
  -S, --no-sambamba     Do not use sambamba, even when it is available. Use pysam instead.
  -tmp, --work-in-tmp   Work in temporary directory
```
fanc pairs
```sh
# fanc pairs
usage: fanc pairs [-h] [-g GENOME] [-r RESTRICTION_ENZYME] [-m] [-u] [-us] [-q QUALITY] [-c CONTAMINANT] [-i INWARD] [-o OUTWARD] [--filter-ligation-auto] [-d REDIST] [-l] [-p DUP_THRESH] [-s STATS]
                  [--reset-filters] [--statistics-plot STATS_PLOT] [--re-dist-plot RE_DIST_PLOT] [--ligation-error-plot LIGATION_ERROR_PLOT] [-t THREADS] [-b BATCH_SIZE] [-S] [-f] [--bwa] [-tmp]
                  input [input ...]

Process and filter read pairs

positional arguments:
  input                 IMPORTANT: The last positional argument will be the output file, unless only a single Pairs object is provided. In that case, filtering and correcting will be done in place. Possible
                        inputs are: two SAM/BAM files (paired-end reads, sorted by read name using "samtools sort -n" or equivalent) and an output file; a HiC-Pro pairs file (format:
                        name<tab>chr1<tab>pos1<tab>strand1<tab>chr2<tab>pos2<tab>strand2) and an output file; a pairs file in 4D Nucleome format (https://github.com/4dn-
                        dcic/pairix/blob/master/pairs_format_specification.md) and an output file, or an existing fanc Pairs object. In case of SAM/BAM, HiC-Pro, or 4D Nucleome you must also provide the
                        --genome argument, and if --genome is not a file with restriction fragments (or Hi-C bins), you must also provide the --restriction-enzyme argument.

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        Path to region-based file (BED, GFF, ...) containing the non-overlapping regions to be used for Hi-C binning. Typically restriction-enzyme fragments. Alternatively: Path to genome
                        file (FASTA, folder with FASTA, HDF5 file), which will be used in conjunction with the type of restriction enzyme to calculate fragments directly.
  -r RESTRICTION_ENZYME, --restriction_enzyme RESTRICTION_ENZYME
                        Name of the restriction enzyme used in the experiment, e.g. HindIII, or MboI. Case-sensitive, only necessary when --genome is provided as FASTA. Restriction names can be any supported
                        by Biopython, which obtains data from REBASE (http://rebase.neb.com/rebase/rebase.html). Separate multiple restriction enzymes with ","
  -m, --filter-unmappable
                        Filter read pairs where one or both halves are unmappable. Only applies to SAM/BAM input!
  -u, --filter-multimapping
                        Filter reads that map multiple times. If the other mapping locations have a lower score than the best one, the best read is kept. Only applies to SAM/BAM input!
  -us, --filter-multimapping-strict
                        Strictly filter reads that map multiple times. Only applies to SAM/BAM input!
  -q QUALITY, --filter-quality QUALITY
                        Cutoff for the minimum mapping quality of a read. For numbers larger than 1, will filter on MAPQ. If a number between 0 and 1 is provided, will filter on the AS tag instead of mapping
                        quality (only BWA). The quality cutoff is then interpreted as the fraction of bases that have to be matched for any given read. Only applies to SAM/BAM input! Default: no mapping
                        quality filter.
  -c CONTAMINANT, --filter-contaminant CONTAMINANT
                        Filter contaminating reads from other organism. Path to mapped SAM/BAM file. Will filter out reads with the same name. Only applies to SAM/BAM input! Default: no contaminant filter
  -i INWARD, --filter-inward INWARD
                        Minimum distance for inward-facing read pairs. Default: no inward ligation error filter
  -o OUTWARD, --filter-outward OUTWARD
                        Minimum distance for outward-facing read pairs. Default: no outward ligation error filter
  --filter-ligation-auto
                        Auto-guess settings for inward/outward read pair filters. Overrides --filter-outward and --filter-inward if set. This is highly experimental and known to overshoot in some cases. It
                        is generally recommended to specify cutoffs manually.
  -d REDIST, --filter-re-distance REDIST
                        Maximum distance for a read to the nearest restriction site. Default: no RE distance filter
  -l, --filter-self-ligations
                        Remove read pairs representing self-ligated fragments.Default: no self-ligation filter.
  -p DUP_THRESH, --filter-pcr-duplicates DUP_THRESH
                        If specified, filter read pairs for PCR duplicates. Parameter determines distance between alignment starts below which they are considered starting at same position. Sensible values
                        are between 1 and 5. Default: no PCR duplicates filter
  -s STATS, --statistics STATS
                        Path for saving filter statistics
  --reset-filters       Remove all filters from the ReadPairs object.
  --statistics-plot STATS_PLOT
                        Path for saving filter statistics plot (PDF)
  --re-dist-plot RE_DIST_PLOT
                        Plot the distribution of restriction site distances of all read pairs (sum left and right read).
  --ligation-error-plot LIGATION_ERROR_PLOT
                        Plot the relative orientation of read pairs mapped to the reference genome as a fraction of reads oriented in the same direction. Allows the identification of ligation errors as a
                        function of genomic distance.
  -t THREADS, --threads THREADS
                        Number of threads to use for extracting fragment information. Default: 1
  -b BATCH_SIZE, --batch-size BATCH_SIZE
                        Batch size for read pairs to be submitted to individual processes. Default: 1000000
  -S, --no-check-sorted
                        Assume SAM files are sorted and do not check if that is actually the case
  -f, --force-overwrite
                        If the specified output file exists, it will be overwritten without warning.
  --bwa                 Use filters appropriate for BWA and not Bowtie2. This will typically be identified automatically from the SAM/BAM header. Set this flag if you are having problems during filtering
                        (typically 0 reads pass the filtering threshold).
  -tmp, --work-in-tmp   Work in temporary directory
```
fanc hic
```sh
# fanc hic
usage: fanc hic [-h] [-b BIN_SIZE] [-l FILTER_LOW_COVERAGE] [-r FILTER_LOW_COVERAGE_RELATIVE] [-a] [-d FILTER_DIAGONAL] [--marginals-plot MARGINALS_PLOT] [--reset-filters] [--downsample DOWNSAMPLE]
                [--subset SUBSET] [-i] [-k] [-n] [-m NORM_METHOD] [-w] [-c] [--only-inter] [-s STATS] [--statistics-plot STATS_PLOT] [--chromosomes CHROMOSOMES [CHROMOSOMES ...]] [-f] [-t THREADS]
                [--deepcopy] [-tmp]
                input [input ...]

Process, filter, and correct Hic files

positional arguments:
  input                 IMPORTANT: The last positional argument will be the output file, unless only a single Hic object is provided. In that case, binning, filtering and correcting will be done in place.
                        Input files. If these are FAN-C Pairs objects (see "fanc pairs"), they will be turned into Hic objects. Hic objects (also the ones converted from Pairs) will first be merged and the
                        merged object will be binned, filtered and corrected as specified in the remaining parameters.

optional arguments:
  -h, --help            show this help message and exit
  -b BIN_SIZE, --bin-size BIN_SIZE
                        Bin size in base pairs. You can use human-readable formats,such as 10k, or 1mb. If omitted, the command will end after the merging step.
  -l FILTER_LOW_COVERAGE, --filter-low-coverage FILTER_LOW_COVERAGE
                        Filter bins with low coverage (lower than specified absolute number of contacts)
  -r FILTER_LOW_COVERAGE_RELATIVE, --filter-low-coverage-relative FILTER_LOW_COVERAGE_RELATIVE
                        Filter bins using a relative low coverage threshold (lower than the specified fraction of the median contact count)
  -a, --low-coverage-auto
                        Filter bins with "low coverage" (under 10% of median coverage for all non-zero bins)
  -d FILTER_DIAGONAL, --diagonal FILTER_DIAGONAL
                        Filter bins along the diagonal up to this specified distance. Use 0 for only filtering the diagonal.
  --marginals-plot MARGINALS_PLOT
                        Plot Hi-C marginals to determine low coverage thresholds.
  --reset-filters       Remove all filters from the Hic object.
  --downsample DOWNSAMPLE
                        Downsample a binned Hi-C object before filtering and correcting. Sample size or reference Hi-C object. If sample size is < 1,will be interpreted as a fraction of valid pairs.
  --subset SUBSET       Comma-separated list of regions that will be used in the output Hic object. All contacts between these regions will be in the output object. For example, "chr1,chr3" will result in a
                        Hic object with all regions in chromosomes 1 and 3, plus all contacts within chromosome 1, all contacts within chromosome 3, and all contacts between chromosome 1 and 3. "chr1" will
                        only contain regions and contacts within chromosome 1.
  -i, --ice-correct     DEPRECATED. Use ICE iterative correction on the binned Hic matrix
  -k, --kr-correct      DEPRECATED. Use Knight-Ruiz matrix balancing to correct the binned Hic matrix
  -n, --normalise       Normalise Hi-C matrix according to --norm-method
  -m NORM_METHOD, --norm-method NORM_METHOD
                        Normalisation method used for -n. Options are: KR (default) = Knight-Ruiz matrix balancing (Fast, accurate, but memory-intensive normalisation); ICE = ICE matrix balancing (less
                        accurate, but more memory-efficient); VC = vanilla coverage (a single round of ICE balancing);VC-SQRT = vanilla coverage square root (reduces overcorrection compared to VC)
  -w, --whole-matrix    Correct the whole matrix at once, rather than individual chromosomes.
  -c, --restore-coverage
                        Restore coverage to the original total number of reads. Otherwise matrix entries will be contact probabilities.
  --only-inter          Calculate bias vector only on inter-chromosomal contacts. Ignores all intra-chromosomal contacts. Always uses whole-matrix balancing, i.e. implicitly sets -w
  -s STATS, --statistics STATS
                        Path for saving filter statistics
  --statistics-plot STATS_PLOT
                        Path for saving filter statistics plot (PDF)
  --chromosomes CHROMOSOMES [CHROMOSOMES ...]
                        Limit output Hic object to these chromosomes. Only available in conjunction with "-b" option.
  -f, --force-overwrite
                        If the specified output file exists, it will be overwritten without warning.
  -t THREADS, --threads THREADS
                        Number of threads (currently used for binning only)
  --deepcopy            Deep copy Hi-C file. Copies a Hi-C file to FAN-C format by duplicating individual bins, pixels, and bias information. Can be used to upgrade an existing FAN-C file with an older
                        version or to convert Cooler or Juicer files to FAN-C format.
  -tmp, --work-in-tmp   Work in temporary directory

```
fancplot  
```sh
#fancplot -h
usage: fancplot [<fancplot global parameters>] <region> [<region> ...]
            --plot <plot type> [<plot parameters>] <plot data file(s)> [...]

            Run fancplot --plot <plot type> -h for help on a specific subplot.

Plot types:

-- Matrix --
triangular    Triangular Hi-C plot
square        Square Hi-C plot
split         Matrix vs matrix plot
mirror        "Mirrored" matrix comparison plot

-- Region --
scores        Region scores plot with parameter dependency
line          Line plot
bar           Bar plot for region scores
layer         Layered feature plot
gene          Gene plot

fancplot plotting tool for fanc

positional arguments:
  regions               List of region selectors (<chr>:<start>-<end>) or files with region information (BED, GTF, ...).

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Suppresses interactive plotting window and redirects plot to file. Specify path to file when plotting a single region, and path to a folder for plotting multiple regions.
  -s SCRIPT, --script SCRIPT
                        Use a script file to define plot.
  -p PLOT, --plot PLOT  New plot, type will be chosen automatically by file type, unless "-t" is provided.
  -n NAME, --name NAME  Plot name to be used as prefix when plotting multiple regions. Is ignored for single region and interactive plot.
  --width WIDTH         Width of the figure in inches. Default: 4
  -w WINDOW_SIZE, --window-size WINDOW_SIZE
                        Plotting region size in base pairs. If provided, the actual size of the given region is ignored and instead a region <chromosome>:<region center - window size/2> - <region center +
                        window size/2> will be plotted.
  --invert-x            Invert x-axis for this plot
  --tick-locations TICK_LOCATIONS [TICK_LOCATIONS ...]
                        Manually define the locations of the tick labels on the genome axis.
  --verbose, -v         Set verbosity level: Can be chained like "-vvv" to increase verbosity. Default is to show errors, warnings, and info messages (same as "-vv"). "-v" shows only errors and warnings,
                        "-vvv" shows errors, warnings, info, and debug messages.
  --silent              Do not print log messages to command line.
  -V, --version         Print version information
  --pdf-text-as-font    When saving a plot to PDF, save text as a font instead of a path. This will increase the file size, sometimes by a lot, but it makes the text in plots editable in vector graphics
                        programs such as Inkscape or Illustrator.
```
fanc expected
```sh
# expect fanc expected -h
2023-12-16 13:03:59,612 INFO FAN-C version: 0.9.27
usage: fanc expected [-h] [-p PLOT_FILE] [-l LABELS [LABELS ...]] [-c CHROMOSOME] [-tmp] [--recalculate] [-N] input [input ...] output

Calculate Hi-C expected values (distance decay)

positional arguments:
  input                 Input matrix (Hi-C, fold-change map, ...)
  output                Output expected contacts (tsv).

optional arguments:
  -h, --help            show this help message and exit
  -p PLOT_FILE, --plot PLOT_FILE
                        Output file for distance decay plot (pdf).
  -l LABELS [LABELS ...], --labels LABELS [LABELS ...]
                        Labels for input objects.
  -c CHROMOSOME, --chromosome CHROMOSOME
                        Specific chromosome to calculate expected values for.
  -tmp, --work-in-tmp   Work in temporary directory
  --recalculate         Recalculate expected values regardless of whether they are already stored in the matrix object.
  -N, --no-norm         Calculate expected values on unnormalised data.

```
fanc compare
```sh
#fanc compare -h
2023-12-14 15:53:46,109 INFO FAN-C version: 0.9.27
usage: fanc compare [-h] [-c COMPARISON] [-o OUTPUT_FORMAT] [-S] [-l] [--log-matrix] [-Z] [-I] [-e] [-u] [-tmp]
                    input input output

Create pairwise comparisons of Hi-C comparison maps

positional arguments:
  input                 Input matrix (e.g. Hic) files.
  output                Output ComparisonMatrix file.

optional arguments:
  -h, --help            show this help message and exit
  -c COMPARISON, --comparison COMPARISON
                        Type of comparison. Default: fold-change, other options are: difference
  -o OUTPUT_FORMAT, --output-format OUTPUT_FORMAT
                        Output format for region-based comparisons. Only relevant when using BED, GFF, or another
                        region-based format as input.
  -S, --no-scale        Do not scale input matrices to the same number of valid pairs
  -l, --log             Log2-convert comparison values (AFTER the comparison)
  --log-matrix          Log2-convert matrices (BEFORE the comparison)
  -Z, --ignore-zero     Do not consider pixels where one matrix entry is zero
  -I, --ignore-infinite
                        Do not consider pixels where the comparison yields "inf"
  -e, --observed-expected
                        O/E transform matrix values before comparison. Only has an effect on matrix comparisons.
  -u, --uncorrected     Compare uncorrected matrices. Only has an effect on matrix comparisons.
  -tmp, --work-in-tmp   Work in temporary directory
(FAN-C) [2022122@admin1 fanc_examples]$

```
fanc  dump -h
```sh
# fanc  dump -h
2023-12-14 16:07:32,209 INFO FAN-C version: 0.9.27
usage: fanc dump [-h] [-s SUBSET] [-S] [--only-intra] [-e] [-l] [-u] [-tmp] hic [matrix] [regions]

Dump Hic file to txt file(s).

positional arguments:
  hic                   Hic file
  matrix                Output file for matrix entries. If not provided, will write to stdout.
  regions               Output file for Hic regions. If not provided, will write regions into matrix file.

optional arguments:
  -h, --help            show this help message and exit
  -s SUBSET, --subset SUBSET
                        Only output this matrix subset. Format: <chr>[:<start>-<end>][--<chr>[:<start><end>]], e.g.:
                        "chr1--chr1" to extract only the chromosome 1 submatrix; "chr2:3400000-4200000" to extract
                        contacts of this region on chromosome 2 to all other regions in the genome;
  -S, --no-sparse       Store full, square matrix instead of sparse format.
  --only-intra          Only dump intra-chromosomal data. Dumps everything by default.
  -e, --observed-expected
                        O/E transform matrix values.
  -l, --log2            Log2-transform matrix values. Useful for O/E matrices (-e option)
  -u, --uncorrected     Output uncorrected (not normalised) matrix values).
  -tmp, --work-in-tmp   Work in temporary directory
(FAN-C) [2022122@admin1 fanc_examples]$
```

```sh
The module ‘fanc auto’ was applied to generate 500, 100, 50, 10 and 1 kb contact matrices (hic files). 
* The resultant hic files with 100-kb resolution were directed to the ‘fanc expected’ module to calculate the expected interaction probability against genomic distance for intrachromosomal interaction. 
* For matrix and score comparisons, the default comparison method of fold-change was used with the ‘fanc compare’ command. 
* The outputs (hic object) were transferred to text files by ‘fanc dump’ and were visualized as heatmaps in R using ggplot2. 

```



## Bisulfite-seq analysis

* 1. Downloaded sequencing reads were processed using TrimGalore (v.0.4.1) (https://github.com/FelixKrueger/TrimGalore) with default parameters. 
* 2. Reads were mapped to TAIR10 using Bismark (v.0.22.2), and methylation was called using MethylDackel (v.0.5.2) (https://github.com/dpryan79/MethylDackel), selecting --CHG and --CHH options. 
* 3. CG methylation data were used in PCA.

```sh
# install reads processed.
conda install trim-galore
# Check that cutadapt is installed
cutadapt --version
# Check that FastQC is installed
fastqc -v
# Install Trim Galore
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
# Run Trim Galore
/public/home/2022122/xugang/app/TrimGalore-0.6.10/trim_galore

```
### Install Bismark.
```sh
# install Bismark.
wget https://github.com/FelixKrueger/Bismark/archive/refs/tags/v0.24.2.tar.gz

tar -xvzf  TrimGalore-0.6.10.tar.gz

vi ~/xugang/bashrc
#/public/home/2022122/xugang/app/Bismark-0.24.2:\

```
### Prepare the genome file
```sh
bismark_genome_preparation /public/home/2022122/xugang/project/tair10/bismark

```

### Run Bismark.
```sh
bismarkf(){
#plot peak heatmap.
log=$output/a11-bismark/log
[[ -d $log ]] || mkdir -p  $log
#remove background file
((counter++))
path=$1
name1=$2
name2=$3
name=$4

[[ -d $output/a9-genecount/ ]] || mkdir -p $output/a9-genecount/  

echo "#!/bin/bash
#SBATCH -o $log/${name1}.%j.out
#SBATCH -e $log/${name1}.%j.error
#SBATCH --partition=${node}
#SBATCH -J 5${1}
#SBATCH -N 1
#SBATCH -n ${thread}
source /public/home/2022122/xugang/bashrc


bismark --parallel ${thread} --temp_dir $log -o $output/a11-bismark/${name} /public/home/2022122/xugang/project/tair10/bismark -1 $path/$name1 -2 $path/$name2

">a11.bismark.$counter.$name.sh
}

bismarkf rawdata/modification mod.A-1_1.fq.gz mod.A-1_2.fq.gz mod.A-1 
bismarkf rawdata/modification mod.A-2_1.fq.gz mod.A-2_2.fq.gz mod.A-2
bismarkf rawdata/modification mod.A-3_1.fq.gz mod.A-3_2.fq.gz mod.A-3

bismarkf rawdata/modification mod.Mu-1_1.fq.gz mod.Mu-1_2.fq.gz mod.Mu-1
bismarkf rawdata/modification mod.Mu-2_1.fq.gz mod.Mu-2_2.fq.gz mod.Mu-2
bismarkf rawdata/modification mod.Mu-3_1.fq.gz mod.Mu-3_2.fq.gz mod.Mu-3

bismarkf rawdata/modification mod.OE-1_1.fq.gz mod.OE-1_2.fq.gz mod.OE-1
bismarkf rawdata/modification mod.OE-2_1.fq.gz mod.OE-2_2.fq.gz mod.OE-2
bismarkf rawdata/modification mod.OE-3_1.fq.gz mod.OE-3_2.fq.gz mod.OE-3

bismarkf rawdata/modification mod.WT-1_1.fq.gz mod.WT-1_2.fq.gz mod.WT-1
bismarkf rawdata/modification mod.WT-2_1.fq.gz mod.WT-2_2.fq.gz mod.WT-2
bismarkf rawdata/modification mod.WT-3_1.fq.gz mod.WT-3_2.fq.gz mod.WT-3


```


```sh

# install MethylDackel.
wget https://github.com/dpryan79/MethylDackel/archive/refs/tags/0.6.1.tar.gz


```


### Module
Usage  
To get a list of all currently loaded modules:
```sh


$ module list
Currently Loaded Modulefiles:
1) DEVELOP                3) intelmpi/2017.4
2) intel/16.0             4) likwid/system-default
detailing that from the category DEVELOP, the Intel Compiler ICC in version 16, the Intel implementation of MPI in Version 2017.4 and the system default version of Likwid are currently loaded and usable.

Calling

$ module avail
lists all available (loadable) modules and module groups. With the information of these two commands, one can:
```
load a specific module
unload a specific module
```sh

$ module load x

$ module unload x
```
load a specific module  
unload a specific module
```sh

$ module load x

$ module unload x
```
swap a specific module for another one (especially useful to switch between different versions of the same program)
```sh

$ module switch x y
```



















































