# Library
library(ggplot2)

setwd('D:/a-document/github/phase_seperate/data/hic-data')
data=read.table('mut_100kb.dum.tsv',header = F,sep="\t")


my_data=data
head(my_data)
my_data2=my_data
my_data2$V4=my_data$V1
my_data2$V5=my_data$V2
my_data2$V6=my_data$V3
my_data2$V1=my_data$V4
my_data2$V2=my_data$V5
my_data2$V3=my_data$V6
my_data3=rbind(my_data,my_data2)
my_data=my_data3

my_data$zscore=(my_data$V7 - mean(my_data$V7 )) / sd(my_data$V7 )
my_data$x=as.numeric(my_data$V1)*1000000000+my_data$V2
my_data$y=as.numeric(my_data$V4)*1000000000+my_data$V6

my_data$labx=paste0('chr',my_data$V1,':',my_data$V2)
my_data$laby=paste0('chr',my_data$V4,':',my_data$V6)
#set the order.How can I order a ggplot heatmap based on the value of a column?
my_data$labx <- reorder(my_data$labx,my_data$x)
my_data$laby = reorder(my_data$laby,my_data$y)

head(my_data)
head(my_data)
max(my_data$V7)
min(my_data$V7)
max(my_data$zscore)
min(my_data$zscore)
my_data$groups=cut(my_data$zscore,
                  breaks=c(-0.1, 0, 0.1, 0.2, 0.3, 100))
head(my_data)
#ggplot(my_data, aes(labx, laby, fill= V7)) + 
#  geom_tile() 
#color1

pdf('mut.color1.pdf')
ggplot(my_data, aes(labx, laby, fill = groups)) +          # Specify colors manually
  geom_tile() +
  scale_fill_manual(breaks = levels(my_data$groups),
                    values = c("#ffffe5", "#fee391", "#fe9929", "#cc4c02", "#662506"))+
  scale_x_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))+
  scale_y_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))
dev.off()





#color2
pdf('mut.color2.pdf')
ggplot(my_data, aes(labx, laby, fill = groups)) +          # Specify colors manually
  geom_tile() +
  scale_fill_manual(breaks = levels(my_data$groups),
                    values = c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D", "#67000D"))+
  scale_x_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))+
  scale_y_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))
dev.off()
#color3
pdf('mut.color3.pdf')
ggplot(my_data, aes(labx, laby, fill = groups)) +          # Specify colors manually
  geom_tile() +
  scale_fill_manual(breaks = levels(my_data$groups),
                    values = c("#3C52A1", "#69ACDE", "#CFEAFA", "#F4AA73", "#CD473E"))+
  scale_x_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))+
  scale_y_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))
dev.off()

########################################
data=read.table('wt_100kb.dum.tsv',header = F,sep="\t")

head(data)

my_data=data
head(my_data)
my_data2=my_data
my_data2$V4=my_data$V1
my_data2$V5=my_data$V2
my_data2$V6=my_data$V3
my_data2$V1=my_data$V4
my_data2$V2=my_data$V5
my_data2$V3=my_data$V6
my_data3=rbind(my_data,my_data2)
my_data=my_data3


my_data$zscore=(my_data$V7 - mean(my_data$V7 )) / sd(my_data$V7 )
my_data$x=as.numeric(my_data$V1)*1000000000+my_data$V2
my_data$y=as.numeric(my_data$V4)*1000000000+my_data$V6

my_data$labx=paste0('chr',my_data$V1,':',my_data$V2)
my_data$laby=paste0('chr',my_data$V4,':',my_data$V6)
#set the order.How can I order a ggplot heatmap based on the value of a column?
my_data$labx <- reorder(my_data$labx,my_data$x)
my_data$laby = reorder(my_data$laby,my_data$y)

head(my_data)
head(my_data)
max(my_data$V7)
min(my_data$V7)
max(my_data$zscore)
min(my_data$zscore)
my_data$groups=cut(my_data$zscore,
                   breaks=c(-0.1, 0, 0.1, 0.2, 0.3, 100))
head(my_data)
#ggplot(my_data, aes(labx, laby, fill= V7)) + 
#  geom_tile() 
#color1
pdf('wt.color1.pdf')
ggplot(my_data, aes(labx, laby, fill = groups)) +          # Specify colors manually
  geom_tile() +
  scale_fill_manual(breaks = levels(my_data$groups),
                    values = c("#ffffe5", "#fee391", "#fe9929", "#cc4c02", "#662506"))+
  scale_x_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))+
  scale_y_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))
dev.off()
#color2
pdf('wt.color2.pdf')
ggplot(my_data, aes(labx, laby, fill = groups)) +          # Specify colors manually
  geom_tile() +
  scale_fill_manual(breaks = levels(my_data$groups),
                    values = c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D", "#67000D"))+
  scale_x_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))+
  scale_y_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))
dev.off()
#color3
pdf('wt.color3.pdf')
ggplot(my_data, aes(labx, laby, fill = groups)) +          # Specify colors manually
  geom_tile() +
  scale_fill_manual(breaks = levels(my_data$groups),
                    values = c("#3C52A1", "#69ACDE", "#CFEAFA", "#F4AA73", "#CD473E"))+
  scale_x_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))+
  scale_y_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))
dev.off()
##############################################################
data=read.table('oe_100kb.dum.tsv',header = F,sep="\t")

my_data=data
head(my_data)
my_data2=my_data
my_data2$V4=my_data$V1
my_data2$V5=my_data$V2
my_data2$V6=my_data$V3
my_data2$V1=my_data$V4
my_data2$V2=my_data$V5
my_data2$V3=my_data$V6
my_data3=rbind(my_data,my_data2)
my_data=my_data3


my_data$zscore=(my_data$V7 - mean(my_data$V7 )) / sd(my_data$V7 )
my_data$x=as.numeric(my_data$V1)*1000000000+my_data$V2
my_data$y=as.numeric(my_data$V4)*1000000000+my_data$V6

my_data$labx=paste0('chr',my_data$V1,':',my_data$V2)
my_data$laby=paste0('chr',my_data$V4,':',my_data$V6)
#set the order.How can I order a ggplot heatmap based on the value of a column?
my_data$labx <- reorder(my_data$labx,my_data$x)
my_data$laby = reorder(my_data$laby,my_data$y)

head(my_data)
head(my_data)
max(my_data$V7)
min(my_data$V7)
max(my_data$zscore)
min(my_data$zscore)
my_data$groups=cut(my_data$zscore,
                   breaks=c(-0.1, 0, 0.1, 0.2, 0.3, 100))
head(my_data)
#ggplot(my_data, aes(labx, laby, fill= V7)) + 
#  geom_tile() 
#color1
pdf('oe.color1.pdf')
ggplot(my_data, aes(labx, laby, fill = groups)) +          # Specify colors manually
  geom_tile() +
  scale_fill_manual(breaks = levels(my_data$groups),
                    values = c("#ffffe5", "#fee391", "#fe9929", "#cc4c02", "#662506"))+
  scale_x_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))+
  scale_y_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))
dev.off()
#color2
pdf('oe.color2.pdf')
ggplot(my_data, aes(labx, laby, fill = groups)) +          # Specify colors manually
  geom_tile() +
  scale_fill_manual(breaks = levels(my_data$groups),
                    values = c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D", "#67000D"))+
  scale_x_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))+
  scale_y_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))
dev.off()
#color3
pdf('oe.color3.pdf')
ggplot(my_data, aes(labx, laby, fill = groups)) +          # Specify colors manually
  geom_tile() +
  scale_fill_manual(breaks = levels(my_data$groups),
                    values = c("#3C52A1", "#69ACDE", "#CFEAFA", "#F4AA73", "#CD473E"))+
  scale_x_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))+
  scale_y_discrete(breaks=c('chr1:1','chr2:1','chr3:1','chr4:1','chr5:1'))
dev.off()

