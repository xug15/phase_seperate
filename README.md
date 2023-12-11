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


bismark --parallel ${thread} -o $output/a11-bismark/${name} /public/home/2022122/xugang/project/tair10/bismark -1 $path/$name1 -2 $path/$name2

">a9.genecount.$counter.$name.sh
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


























































