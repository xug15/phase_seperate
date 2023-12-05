## 一、技术原理
Hi-C是染色质区域捕获（Chromosome conformation capture）与高通量测序（High-throughput sequencing）相结合而产生的一种新技术。

以整个细胞核为研究对象，利用高通量测序技术，结合生物信息分析方法，研究全基因组范围内整个染色质DNA在空间位置上的关系，获得高分辨率的染色质调控元件相互作用图谱。

通过甲醛处理将DNA与蛋白质交联在一起从而固定DNA的构象，并进行酶切、生物素引入、酶连和提取，然后对处理好的核酸进行文库构建和高通量测序，最终通过对测序数据进行分析即可揭示染色体片段间的交互信息。

该技术仅利用单个个体就可以将基因组序列定位到染色体，从而进行染色体水平的研究。Hi-C可以与RNA-Seq、ChIP-Seq等数据进行联合分析，从基因调控网络和表观遗传网络来阐述生物体性状形成的相关机制。

强术语版：

Hi-C技术主要将空间结构临近的DNA片段进行交联，并将交联的DNA片段富集，然后进行高通量测序，对测序数据进行分析即可揭示全基因组范围内的染色体片段间的交互作用。利用Hi-C技术可以揭示基因组的一般结构特征，包括从隔室（动物中A/B Compartments，植物中为CSD/LSD）到拓扑相关结构域domain（动物中TAD，植物中TAD-like）,最后再到环（loop）的染色质的这种层级结构。

染色体的三维（3D）结构深刻影响DNA复制，转录及DNA损伤修复。最近研究揭示了人基因组中存在数百万的潜在的顺式调控元件，其中大量位于基因间区和远离其靶基因的启动子区。这些主要由增强子组成的远端调控元件，在生物发育过程中通过基因组中形成的Loop结构影响目标基因的转录，从而实现对目标基因的调控．另外通过结合BS-seq、ChIP-seq和RNA-seq数据分析发现揭示染色体空间结构与基因组表观遗传修饰和转录活性关系，共同来阐释核基因组的包装模式。



1、用甲醛对细胞进行固定，使DNA与蛋白，蛋白与蛋白之间进行交联；  
2、进行酶切（如Hind III等限制性内切酶），使交联两侧产生粘性末端；  
3、末端修复，引入生物素标记，连接；  
4、解交联，使DNA和蛋白、蛋白和蛋白分开，提取DNA，打断，捕获带有生物素标记片段，进行建库；  
5、测序。


## 二、技术优势
1. 通过Scaffold间的交互频率大小，可以对已组装的基因组序列进行纠错（基因组更准确）。

2. 基因信息不再仅仅是contig片段，而是被划分至染色体上，成为染色体水平。

3. 无需辛苦的构建群体，单一一个体就能实现染色体定位。

4. 相比遗传图谱，标记密度更大，序列定位更完整（能把更多的contig挂至染色体上!信息更全面！）

5. 染色体重排等结构变异研究可以开展啦~(研究可以更深入！）

6. QTL、GWAS可以定位区间到某个染色体啦~（追踪变异！）

7. 该物种的三维基因结构、染色体互作及动态变化可以解析啦~（从基因到表观！全方位解析）

8. 周期短。

总结：快（周期短）、准（准确性高）、省（单个体）

## 三、应用：
### 1. 染色体跨度的单倍型图谱构建

为疾病风险预测提供思路

为肿瘤形成机制提供依据

为农业动植物经济性状连锁标记及基因组进化奠定基础

### 2. 探索基因组的3D结构

有助于了解基因组折叠对基因表达的影响

有助于了解基因组三维结构对细胞发育，分化及细胞命运的决定

### 3. 开发调控基因的DNA元件

揭示基因组远程调控元件介导的分子网络

通过分析调控元件空间互作热点，开发潜在的药物作用靶点

### 4. Hi-C辅助基因组组装

辅助动植物基因组组装

宏基因组组装及菌群精确分类

#### 4.1 辅助组装应用

1）PacBio+Hi-C

利用了PacBio三代分子测序和北京百迈客生物技术有限公司的Hi-C技术绘制了三个高质量棉花基因图谱（题目：Genome sequence ofGossypium herbaceumand genome updates ofGossypium arboreumandGossypium hirsutumprovide insights into cotton A-genome evolution  期刊：Nature genetics ）

2）SMRT+Hi-C

利用单分子实时测序SMRT和北京百迈客生物技术有限公司的Hi-C技术组装了一个高质量、染色体体水平的珙桐基因组，研究发现苞片中光合作用相关基因几乎缺失或表达减少，而抗菌、抗冷、抗水等抗逆相关基因在苞片中高度表达，突出了苞片在保护花和吸引授粉者中的重要作用。（题目：Genomic analyses of a “living fossil”: The endangered dove-tree 期刊：Molecular Ecology Resources ）

3）Illumina+PacBio组装+Hi-C技术将95.15%的序列挂载

该研究利用Illumina+PacBio组装了1.12Gb油桐基因组，并通过北京百迈客生物技术有限公司的Hi-C技术将95.15%的序列挂载到11条染色体上。基于比较基因组等系列相关研究解析了油桐基因组进化和桐油生物合成的分子机制。（题目：Tung Tree (Vernicia fordii) Genome Provides A Resource for Understanding Genome Evolution and Improved Oil Production 期刊：Genomics Proteomics & Bioinformatics ）

4）Nanopore测序+Hi-C染色体构象捕获技术

构建了高质量的枇杷基因组。（题目：Chromosome-level genome assembly and annotation of the loquat (Eriobotrya japonica) genome 期刊：Gigascience）

## 四、Hi-C分析常用工具
### 数据标准化：

1. HiCNorm

http://www.people.fas.harvard.edu/~junliu/HiCNorm/

2. ICE

https://mirnylab.bitbucket.io/hiclib/index.html

3. HiC-Pro

https://github.com/nservant/HiC-Pro

### TAD鉴定

1. HiCseg:Models the uncertainty in Hi‐C data

https://cran.r-project.org/web/packages/HiCseg/index.html

2. TADbit

https://github.com/3DGenomes/TADbit

3. DomainCaller

http://chromosome.sdsc.edu/mouse/hi-c/download.html

4. InsulationScore:Robust to different sequencing depth;

can detect dynamics of TAD boundaries

https://github.com/dekkerlab/crane-nature-2015

5. Arrowhead:High computational efficiency with

### dynamic programming

https://github.com/theaidenlab/juicer/wiki/Download

6. TADtree

compbio.cs.brown.edu/projects/tadtree/

7. Armatus:TAD calling robust in different resolutions

https://github.com/kingsfordgroup/armatus

8. Topdom

http://zhoulab.usc.edu/TopDom/

### 交互片段鉴定（interaction）  

1. Fit-Hi-C: Accurate background model using
non-parametric spline

http://noble.gs.washington.edu/proj/fit-hi-c

2. GOTHiC: Models contact-frequency uncertainty as binomial distribution

http://bioconductor.org/packages/release/bioc/html/GOTHiC.html

3. HOMER

homer.ucsd.edu/homer/download.html

4. HIPPIE

wanglab.pcbi.upenn.edu/hippie

5. diffHic

https://bioconductor.org/packages/release/bioc/html/diffHic.html

6. HiCCUPS:Designed for high-resolution Hi‐C data

https://github.com/theaidenlab/juicer/wiki/Download

### 3D构象：  

1. 3D-GNOME

https://bitbucket.org/3dome/3dome_mmc

2. Tadbit

https://github.com/3DGenomes/tadbit

### 可视化：

1. HiCPlotter

https://github.com/kcakdemir/HiCPlotter

