#### relative abundance
ASV.relative <- decostand(ASV.table, "total")
rowSums(ASV.relative)

ASV.MG <- ASV.relative[,rownames(taxa.MG)]

taxa.MT <- taxa.MT[which(taxa.MT$Family!="Methylacidiphilaceae"),]
ASV.MT <- ASV.relative[,rownames(taxa.MT)]

sample.data <- data.frame(sample_data(physeq2))
sample.data$Locality <- gsub("Control","Intact",sample.data$Locality) 
sample.data$Locality[88:108] <- "Intermediate"

sample.data$group <- paste(sample.data$year,sample.data$Locality,sep = "-")

## pairwise Kruskal-Wallis
library(PMCMR)
posthoc.kruskal.dunn.test(MT ~ factor(Locality), data=sample.data[which(sample.data$year==2019),], p.adjust="fdr")
posthoc.kruskal.dunn.test(MT ~ factor(Locality), data=sample.data[which(sample.data$year==2021),], p.adjust="fdr")
posthoc.kruskal.dunn.test(MG ~ factor(year), data=sample.data, p.adjust="fdr")

## overall abundance
# MG
ASV.MG$total <- rowSums(ASV.MG)
sample.data$MG <- ASV.MG$total



library(plyr)
MG.site <- ddply(sample.data, c("year","Locality"), summarise,
                          N    = length(MG),
                          mean = mean(MG*100),
                          sd   = sd(MG*100),
                          se   = sd / sqrt(N)
)


MG.depth <-  ddply(sample.data, c("year","Depth_proxy"), summarise,
                   N    = length(MG),
                   mean = mean(MG*100),
                   sd   = sd(MG*100),
                   se   = sd / sqrt(N)
)

MG.site$xvalue <- paste(MG.site$Locality,MG.site$N,sep = ", n =")
MG.depth$xvalue <- paste(MG.depth$Depth_proxy,MG.depth$N,sep = ", n =")

addline_format <- function(x,...){
  gsub(',','\n',x)
}

MG.site$xvalue <- addline_format(MG.site$xvalue)
MG.depth$xvalue <- addline_format(MG.depth$xvalue)
MG.depth$xvalue <- factor(MG.depth$xvalue,levels = c("0-15cm\n n =17", "0-15cm\n n =24",  "16-30cm\n n =11",  "16-30cm\n n =3",   "31-45cm\n n =7",   "31-45cm\n n =8",
                                                     "46-60cm\n n =13",  "46-60cm\n n =6",   "61-75cm\n n =1" ,  "61-75cm\n n =4",  
                           "76-90cm\n n =1",   "76-90cm\n n =3",   "91-105cm\n n =2" , "91-105cm\n n =4", "106-120cm\n n =4"))
MG.depth$xvalue <- factor(MG.depth$xvalue,levels =rev(levels(MG.depth$xvalue)))


library(ggplot2)
p1 <- ggplot(MG.site[1:3,]) +
  geom_col(aes(x=xvalue, y=mean),color="black",fill="grey50")+
  geom_errorbar(aes(x=xvalue, ymin=mean, ymax=mean+se))+
  xlab("")+ ylab("Relative abundance (%)")+
  scale_y_continuous(limits = c(0,0.055),expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90, size = 15, color = "black",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.y = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),strip.text.x = element_text(size = 17, color = "black",face = "bold"),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())+facet_wrap(year~.,scales = "free",ncol = 1)
p1


p2 <- ggplot(MG.site[4:6,]) +
  geom_col(aes(x=xvalue, y=mean),color="black",fill="grey50")+
  geom_errorbar(aes(x=xvalue, ymin=mean, ymax=mean+se))+
  xlab("")+ ylab("Relative abundance (%)")+
  scale_y_continuous(limits = c(0,0.055),expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90, size = 15, color = "black",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.y = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),strip.text.x = element_text(size = 17, color = "black",face = "bold"),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())+facet_wrap(year~.,scales = "free",ncol = 1)
p2

p3 <- ggplot(MG.depth[1:8,]) +
  geom_col(aes(x=xvalue, y=mean),color="black",fill="grey50")+
  geom_errorbar(aes(x=xvalue, ymin=mean, ymax=mean+se))+
  xlab("")+ ylab("Relative abundance (%)")+
  scale_y_continuous(limits = c(0,0.055),expand = c(0,0))+
  theme(axis.text.y = element_text(size = 15, color = "black",face = "bold"),
        axis.text.x = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.x = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),strip.text.x = element_text(size = 17, color = "black",face = "bold"),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())+
  facet_wrap(.~year,scales = "free",ncol = 1)+
  coord_flip()
p3

p4 <- ggplot(MG.depth[9:15,]) +
  geom_col(aes(x=xvalue, y=mean),color="black",fill="grey50")+
  geom_errorbar(aes(x=xvalue, ymin=mean, ymax=mean+se))+
  xlab("")+ ylab("Relative abundance (%)")+
  scale_y_continuous(limits = c(0,0.055),expand = c(0,0))+
  theme(axis.text.y = element_text(size = 15, color = "black",face = "bold"),
        axis.text.x = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.x = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),strip.text.x = element_text(size = 17, color = "black",face = "bold"),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())+
  facet_wrap(.~year,scales = "free",ncol = 1)+
  coord_flip()
p4


library(cowplot)
plot_grid(p1,p3,p2,p4,ncol = 2,align = "vh",axis = "tblr")


# MT
ASV.MT$total <- rowSums(ASV.MT)
sample.data$MT <- ASV.MT$total



library(plyr)
MT.site <- ddply(sample.data, c("year","Locality"), summarise,
                 N    = length(MT),
                 mean = mean(MT*100),
                 sd   = sd(MT*100),
                 se   = sd / sqrt(N)
)


MT.depth <-  ddply(sample.data, c("year","Depth_proxy"), summarise,
                   N    = length(MT),
                   mean = mean(MT*100),
                   sd   = sd(MT*100),
                   se   = sd / sqrt(N)
)

MT.site$xvalue <- paste(MT.site$Locality,MT.site$N,sep = ", n =")
MT.depth$xvalue <- paste(MT.depth$Depth_proxy,MT.depth$N,sep = ", n =")

addline_format <- function(x,...){
  gsub(',','\n',x)
}

MT.site$xvalue <- addline_format(MT.site$xvalue)
MT.depth$xvalue <- addline_format(MT.depth$xvalue)
MT.depth$xvalue <- factor(MT.depth$xvalue,levels = c("0-15cm\n n =17", "0-15cm\n n =24",  "16-30cm\n n =11",  "16-30cm\n n =3",   "31-45cm\n n =7",   "31-45cm\n n =8",
                                                     "46-60cm\n n =13",  "46-60cm\n n =6",   "61-75cm\n n =1" ,  "61-75cm\n n =4",  
                                                     "76-90cm\n n =1",   "76-90cm\n n =3",   "91-105cm\n n =2" , "91-105cm\n n =4", "106-120cm\n n =4"))
MT.depth$xvalue <- factor(MT.depth$xvalue,levels =rev(levels(MT.depth$xvalue)))


library(ggplot2)
p5 <- ggplot(MT.site[1:3,]) +
  geom_col(aes(x=xvalue, y=mean),color="black",fill="grey50")+
  geom_errorbar(aes(x=xvalue, ymin=mean, ymax=mean+se))+
  xlab("")+ ylab("Relative abundance (%)")+
  scale_y_continuous(limits = c(0,0.25),expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90, size = 15, color = "black",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.y = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),strip.text.x = element_text(size = 17, color = "black",face = "bold"),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())+facet_wrap(year~.,scales = "free",ncol = 1)
p5


p6 <- ggplot(MT.site[4:6,]) +
  geom_col(aes(x=xvalue, y=mean),color="black",fill="grey50")+
  geom_errorbar(aes(x=xvalue, ymin=mean, ymax=mean+se))+
  xlab("")+ ylab("Relative abundance (%)")+
  scale_y_continuous(limits = c(0,0.25),expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90, size = 15, color = "black",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.y = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),strip.text.x = element_text(size = 17, color = "black",face = "bold"),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())+facet_wrap(year~.,scales = "free",ncol = 1)
p6

p7 <- ggplot(MT.depth[1:8,]) +
  geom_col(aes(x=xvalue, y=mean),color="black",fill="grey50")+
  geom_errorbar(aes(x=xvalue, ymin=mean, ymax=mean+se))+
  xlab("")+ ylab("Relative abundance (%)")+
  scale_y_continuous(limits = c(0,0.35),expand = c(0,0))+
  theme(axis.text.y = element_text(size = 15, color = "black",face = "bold"),
        axis.text.x = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.x = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),strip.text.x = element_text(size = 17, color = "black",face = "bold"),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())+
  facet_wrap(.~year,scales = "free",ncol = 1)+
  coord_flip()
p7

p8 <- ggplot(MT.depth[9:15,]) +
  geom_col(aes(x=xvalue, y=mean),color="black",fill="grey50")+
  geom_errorbar(aes(x=xvalue, ymin=mean, ymax=mean+se))+
  xlab("")+ ylab("Relative abundance (%)")+
  scale_y_continuous(limits = c(0,0.35),expand = c(0,0))+
  theme(axis.text.y = element_text(size = 15, color = "black",face = "bold"),
        axis.text.x = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.x = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),strip.text.x = element_text(size = 17, color = "black",face = "bold"),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())+
  facet_wrap(.~year,scales = "free",ncol = 1)+
  coord_flip()
p8


library(cowplot)
plot_grid(p5,p7,p6,p8,ncol = 2,align = "vh",axis = "tblr")


#### bubble plots
# MG
sample.data$group2 <- paste(sample.data$group,sample.data$Depth_proxy,sep = "-")
MG.bubble <- aggregate(ASV.MG[,1:5],list(sample.data$group2),mean)
MG.bubble <- reshape2::melt(MG.bubble,id.vars="Group.1")

library(stringr)
MG.bubble <- cbind(MG.bubble,str_split_fixed(MG.bubble$Group.1, "-", 3))

colnames(MG.bubble) <- c("group","ASV","value","Year","Site","Depth")

taxa.MG2 <- read.csv("../Erik_thesis2/revision/tax_silva138.2_MG_ASV.csv",row.names = 1)

MG.bubble$family <- MG.bubble$ASV
for (i in 1:5) {
  MG.bubble$family <- gsub(rownames(taxa.MG2)[i],taxa.MG2$Family[i],MG.bubble$family)
}

MG.bubble$order<- MG.bubble$ASV
for (i in 1:5) {
  MG.bubble$order <- gsub(rownames(taxa.MG2)[i],taxa.MG2$Order[i],MG.bubble$order)
}


MG.bubble$genus<- MG.bubble$ASV
for (i in 1:5) {
  MG.bubble$genus <- gsub(rownames(taxa.MG2)[i],taxa.MG2$Genus2[i],MG.bubble$genus)
}

# MG.bubble$group2 <- paste(MG.bubble$order,MG.bubble$ASV,sep = "-")

p9 <- ggplot(MG.bubble[which(MG.bubble$Year=="2019"),],aes(x=factor(Depth,levels = c("0-15cm",   "16-30cm"  , "31-45cm"   ,"46-60cm" , "61-75cm",   "76-90cm" ,  "91-105cm", "106-120cm" )),
                               y=genus,size=ifelse(value==0, NA, value)*100))+
  geom_point(aes(fill=order),color="black",shape=21,alpha=0.8)+scale_y_discrete(position = "right",limits=rev(taxa.MG2$Genus2))+
  scale_size(name="Relative abundance (%)",range = c(0.5,15),limits = c(0,0.25),breaks = c(0.05,0.1,0.15,0.2))+
  guides(fill=guide_legend(title="",override.aes = list(size=3.5)),size="none")+
  scale_fill_brewer(palette = "Dark2")+
  facet_grid(cols=vars(Site), scales = "free",space="free")+
  theme(axis.title = element_blank(),axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12,angle=90,face = "bold",hjust = 1,vjust = 0.5),legend.position = "left",
        strip.text = element_text(size = 12,face = "bold"),legend.background = element_rect(fill = NA),legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 10),legend.title = element_text(size = 12,face = "bold"),panel.grid.major = element_line(linetype = "dotted", colour = "grey90"),
        panel.background = element_rect(fill = NA), panel.border = element_rect(linetype = "solid", colour = "black",fill=NA,size = 0.5))
  
p9


p10 <- ggplot(MG.bubble[which(MG.bubble$Year=="2021"),],aes(x=factor(Depth,levels = c("0-15cm",   "16-30cm"  , "31-45cm"   ,"46-60cm" , "61-75cm",   "76-90cm" ,  "91-105cm", "106-120cm" )),
                                                           y=genus,size=ifelse(value==0, NA, value)*100))+
  geom_point(aes(fill=order),color="black",shape=21,alpha=0.8)+scale_y_discrete(position = "right",limits=rev(taxa.MG2$Genus2))+
  scale_size(name="Relative abundance (%)",range = c(0.5,15),limits = c(0,0.25),breaks = c(0.05,0.1,0.15,0.2))+
  guides(fill="none",size="none")+
  scale_fill_brewer(palette = "Dark2")+
  facet_grid(cols=vars(Site), scales = "free",space="free")+
  theme(axis.title = element_blank(),axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12,angle=90,face = "bold",hjust = 1,vjust = 0.5),legend.position = "left",
        strip.text = element_text(size = 12,face = "bold"),legend.background = element_rect(fill = NA),legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 10),legend.title = element_text(size = 12,face = "bold"),panel.grid.major = element_line(linetype = "dotted", colour = "grey90"),
        panel.background = element_rect(fill = NA), panel.border = element_rect(linetype = "solid", colour = "black",fill=NA,size = 0.5))

p10




# MT
MT.bubble <- aggregate(ASV.MT[,1:26],list(sample.data$group2),mean)
MT.bubble <- reshape2::melt(MT.bubble,id.vars="Group.1")

MT.bubble <- cbind(MT.bubble,str_split_fixed(MT.bubble$Group.1, "-", 3))

colnames(MT.bubble) <- c("group","ASV","value","Year","Site","Depth")

taxa.MT2 <- read.csv("../Erik_thesis2/revision/tax_silva138.2_MT_ASV.csv",row.names = 1)

MT.bubble$family <- MT.bubble$ASV
for (i in 1:26) {
  MT.bubble$family <- gsub(rownames(taxa.MT2)[i],taxa.MT2$Family[i],MT.bubble$family)
}

MT.bubble$genus<- MT.bubble$ASV
for (i in 1:26) {
  MT.bubble$genus <- gsub(rownames(taxa.MT2)[i],taxa.MT2$Genus2[i],MT.bubble$genus)
}

MG.bubble$genus<- MG.bubble$ASV
for (i in 1:5) {
  MG.bubble$genus <- gsub(rownames(taxa.MG2)[i],taxa.MG2$Genus2[i],MG.bubble$genus)
}

# MT.bubble$group2 <- paste(MT.bubble$family,MT.bubble$genus,sep = ":")
# MT.bubble$group2 <- paste(MT.bubble$group2,MT.bubble$ASV,sep = "-")
# 
# MT.bubble$group2 <- gsub("Candidatus","Ca.", MT.bubble$group2)
# MT.bubble$group2 <- gsub("NA-","Methylobacter-",MT.bubble$group2)
# MT.bubble$genus[is.na(MT.bubble$genus)] <- "Methylobacter"

p11 <- ggplot(MT.bubble[which(MT.bubble$Year=="2019"),],aes(x=factor(Depth,levels = c("0-15cm",   "16-30cm"  , "31-45cm"   ,"46-60cm" , "61-75cm",   "76-90cm" ,  "91-105cm", "106-120cm" )),
                                                           y=genus,size=ifelse(value==0, NA, value)*100))+
  geom_point(aes(fill=family),color="black",shape=21,alpha=0.8)+scale_y_discrete(position = "right",limits=rev(taxa.MT2$Genus2))+
  scale_size(name="Relative abundance (%)",range = c(0.5,15),limits = c(0,0.25),breaks = c(0.05,0.1,0.15,0.2))+
  guides(fill=guide_legend(title="",override.aes = list(size=3.5)),size=guide_legend(override.aes = list(shape=19)))+
  scale_fill_brewer(palette = "Dark2")+
  facet_grid(cols=vars(Site), scales = "free",space="free")+
  theme(axis.title = element_blank(),axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12,angle=90,face = "bold",hjust = 1,vjust = 0.5),legend.position = "left",
        strip.text = element_text(size = 12,face = "bold"),legend.background = element_rect(fill = NA),legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 10),legend.title = element_text(size = 12,face = "bold"),panel.grid.major = element_line(linetype = "dotted", colour = "grey90"),
        panel.background = element_rect(fill = NA), panel.border = element_rect(linetype = "solid", colour = "black",fill=NA,size = 0.5))

p11


p12 <- ggplot(MT.bubble[which(MT.bubble$Year=="2021"),],aes(x=factor(Depth,levels = c("0-15cm",   "16-30cm"  , "31-45cm"   ,"46-60cm" , "61-75cm",   "76-90cm" ,  "91-105cm", "106-120cm" )),
                                                            y=genus,size=ifelse(value==0, NA, value)*100))+
  geom_point(aes(fill=family),color="black",shape=21,alpha=0.8)+scale_y_discrete(position = "right",limits=rev(taxa.MT2$Genus2))+
  scale_size(name="Relative abundance (%)",range = c(0.5,15),limits = c(0,0.25),breaks = c(0.05,0.1,0.15,0.2))+
  guides(fill="none",size="none")+
  scale_fill_brewer(palette = "Dark2")+
  facet_grid(cols=vars(Site), scales = "free",space="free")+
  theme(axis.title = element_blank(),axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12,angle=90,face = "bold",hjust = 1,vjust = 0.5),legend.position = "left",
        strip.text = element_text(size = 12,face = "bold"),legend.background = element_rect(fill = NA),legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 10),legend.title = element_text(size = 12,face = "bold"),panel.grid.major = element_line(linetype = "dotted", colour = "grey90"),
        panel.background = element_rect(fill = NA), panel.border = element_rect(linetype = "solid", colour = "black",fill=NA,size = 0.5))

p12

plot_grid(p9,p10,p11,p12,ncol = 2,align = "vh",axis = "tblr",rel_heights = c(1,2.8))


write.csv(sample.data,"sample.data.csv")


############ correlations between MT & parameters
cor.2019 <- read.csv("sample.data.rev.2019.csv",header = T, stringsAsFactors = F, row.names = 1)
cor.2021 <- read.csv("sample.data.rev.2021.csv",header = T, stringsAsFactors = F, row.names = 1)

taxa.MT.typeI <- taxa.MT[which(taxa.MT$Family=="Methylomonadaceae"),]
ASV.MT.typeI <- ASV.relative[,rownames(taxa.MT.typeI )]

ASV.MT.typeI $total <- rowSums(ASV.MT.typeI)
sample.data$MT.typeI <- ASV.MT.typeI$total


taxa.MT.typeII <- taxa.MT[which(taxa.MT$Family=="Beijerinckiaceae"),]
ASV.MT.typeII <- ASV.relative[,rownames(taxa.MT.typeII)]

ASV.MT.typeII$total <- rowSums(ASV.MT.typeII)
sample.data$MT.typeII <- ASV.MT.typeII$total

cor.2019$MT.typeI <- sample.data$MT.typeI[53:108]
cor.2019$MT.typeII <- sample.data$MT.typeII[53:108]

cor.2021 <- cor.2021[rownames(sample.data)[1:52],]
cor.2021$MT.typeI <- sample.data$MT.typeI[1:52]
cor.2021$MT.typeII <- sample.data$MT.typeII[1:52]

#### 2019
cor.coeffecints.2019 <- cor.2019[1:3,10:19]
rownames(cor.coeffecints.2019) <- c("Methanogens","Type I MOB","Type II MOB")

cor.p.values.2019 <- cor.coeffecints.2019

cor.2019.sort <- cor.2019[,c(7,20,21,10:19)]

library(vegan)

for (j in 1:10) {
  for (i in 1:3) {
    cor.coeffecints.2019[i,j] <- cor.test(cor.2019.sort[,i],cor.2019.sort[,j+3],method = "spearman")$estimate
    cor.p.values.2019[i,j] <- cor.test(cor.2019.sort[,i],cor.2019.sort[,j+3],method = "spearman")$p.value
  }
}


cor.p.values.adjusted.2019 <- matrix(p.adjust(as.vector(as.matrix(cor.p.values.2019)), method='fdr'),ncol=10)

cor.p.values.adjusted.2019 <- data.frame(cor.p.values.adjusted.2019)


colnames(cor.coeffecints.2019) <- c("Depth","Moisture","pH","SOC","TN","C/N","BS","Temperature","DOC","DTN")


cor.test(cor.2019.sort[,1],cor.2019.sort[,3],method = "pearson")


library(reshape2)

cor.coeffecints.2019$clade <- rownames(cor.coeffecints.2019)
cor.coeffecints.melt.2019 <- melt(cor.coeffecints.2019,id.vars = "clade")

rownames(cor.p.values.adjusted.2019) <- rownames(cor.coeffecints.2019)
colnames(cor.p.values.adjusted.2019) <- colnames(cor.p.values.adjusted.2019)
cor.p.values.adjusted.2019$clade <- rownames(cor.coeffecints.2019)
cor.p.values.adjusted.melt.2019 <- melt(cor.p.values.adjusted.2019,id.vars = "clade")

cor.coeffecints.melt.2019$p.adjust <- cut(cor.p.values.adjusted.melt.2019$value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

library(ggplot2)

p13 <- ggplot(cor.coeffecints.melt.2019, aes(variable, factor(clade), fill= value)) +
  geom_tile(color="black") + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C",midpoint = 0,limits = c(-0.6,0.6))+
  scale_y_discrete(expand = c(0,0))+scale_x_discrete(expand = c(0,0))+
  labs(x=NULL,y=NULL,fill="Coefficient")+geom_text(aes(label=p.adjust),color="black",size=6)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text.y = element_text(size=12,face = "bold"),axis.text.x = element_text(size=12,face = "bold",angle = 30,vjust = 0.7),
                   legend.text = element_text(size = 10,face = "bold"),legend.title = element_text(size = 12,face = "bold"),strip.text = element_blank())


#### 2021
cor.coeffecints.2021 <- cor.2021[1:3,10:22]
rownames(cor.coeffecints.2021) <- c("Methanogens","Type I MOB","Type II MOB")

cor.p.values.2021 <- cor.coeffecints.2021

cor.2021.sort <- cor.2021[,c(7,23,24,10:22)]


for (j in 1:13) {
  for (i in 1:3) {
    cor.coeffecints.2021[i,j] <- cor.test(cor.2021.sort[,i],cor.2021.sort[,j+3],method = "spearman")$estimate
    cor.p.values.2021[i,j] <- cor.test(cor.2021.sort[,i],cor.2021.sort[,j+3],method = "spearman")$p.value
  }
}


cor.p.values.adjusted.2021 <- matrix(p.adjust(as.vector(as.matrix(cor.p.values.2021)), method='fdr'),ncol=13)

cor.p.values.adjusted.2021 <- data.frame(cor.p.values.adjusted.2021)


colnames(cor.coeffecints.2021) <- c("Depth","Moisture","pH","SOC","TN","C/N","BS","MiA%","MiA.OC","MiA.N","MaA%","MaA.OC","MaA.N")


cor.test(cor.2021.sort[,1],cor.2021.sort[,3],method = "spearman")

cor.coeffecints.2021$clade <- rownames(cor.coeffecints.2021)
cor.coeffecints.melt.2021 <- melt(cor.coeffecints.2021,id.vars = "clade")

rownames(cor.p.values.adjusted.2021) <- rownames(cor.coeffecints.2021)
colnames(cor.p.values.adjusted.2021) <- colnames(cor.p.values.adjusted.2021)
cor.p.values.adjusted.2021$clade <- rownames(cor.coeffecints.2021)
cor.p.values.adjusted.melt.2021 <- melt(cor.p.values.adjusted.2021,id.vars = "clade")

cor.coeffecints.melt.2021$p.adjust <- cut(cor.p.values.adjusted.melt.2021$value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))


p14 <- ggplot(cor.coeffecints.melt.2021, aes(variable, factor(clade), fill= value)) +
  geom_tile(color="black") + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C",midpoint = 0,limits = c(-0.6,0.6))+
  scale_y_discrete(expand = c(0,0))+scale_x_discrete(expand = c(0,0))+
  labs(x=NULL,y=NULL,fill="Coefficient")+geom_text(aes(label=p.adjust),color="black",size=6)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text.y = element_text(size=12,face = "bold"),axis.text.x = element_text(size=12,face = "bold",angle = 30,vjust = 0.7),
                   legend.text = element_text(size = 10,face = "bold"),legend.title = element_text(size = 12,face = "bold"),strip.text = element_blank())

library(cowplot)
plot_grid(p13,p14,ncol = 2,rel_widths = c(10,12.1),align = "vh")



##### metatranscriptome 2022 & qPCR 2019
metatrans.2022 <- read.csv("../RNA_samples/Metatrans_AK/data_kegg_abs.csv")
metatrans.2022 <- metatrans.2022[grep("K10944",metatrans.2022$F5),]
metatrans.2022 <- data.frame(t(metatrans.2022[,1:9]))
metatrans.2022$site <- c("Dry","Dry","Dry","Wet","Wet","Wet","Intact","Intact","Intact")
library(plyr)
metatrans.site <- ddply(metatrans.2022, c("site"), summarise,
                 N    = length(X1273),
                 mean = mean(X1273),
                 sd   = sd(X1273),
                 se   = sd / sqrt(N)
)
library(ggplot2)
ggplot(metatrans.site) +
  geom_col(aes(x=site, y=mean),color="black",fill="grey50")+
  geom_errorbar(aes(x=site, ymin=mean, ymax=mean+se))+
  xlab("")+ ylab("amoA/pmoA transcripts per g dry soil")+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90, size = 15, color = "black",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.y = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),strip.text.x = element_text(size = 17, color = "black",face = "bold"),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())
library(PMCMR)
PMCMRplus::kwAllPairsDunnTest(X1273 ~ factor(site), data=metatrans.2022, p.adjust="fdr")



sample.data.2019 <- read.csv("sample.data.rev.2019.csv",header = T, stringsAsFactors = F, row.names = 1)
# mcrA total
library(plyr)
mcrA.total.site <- ddply(sample.data.2019, c("Locality"), summarise,
                 N    = length(mcrA_Total),
                 mean = mean(mcrA_Total),
                 sd   = sd(mcrA_Total),
                 se   = sd / sqrt(N)
)


mcrA.total.depth <-  ddply(sample.data.2019, c("Depth_proxy"), summarise,
                   N    = length(mcrA_Total),
                   mean = mean(mcrA_Total),
                   sd   = sd(mcrA_Total),
                   se   = sd / sqrt(N)
)

mcrA.total.site$xvalue <- paste(mcrA.total.site$Locality,mcrA.total.site$N,sep = ", n =")
mcrA.total.depth$xvalue <- paste(mcrA.total.depth$Depth_proxy,mcrA.total.depth$N,sep = ", n =")

addline_format <- function(x,...){
  gsub(',','\n',x)
}

mcrA.total.site$xvalue <- addline_format(mcrA.total.site$xvalue)
mcrA.total.depth$xvalue <- addline_format(mcrA.total.depth$xvalue)
mcrA.total.depth$xvalue <- factor(mcrA.total.depth$xvalue,levels = c("0-15cm\n n =17",  "16-30cm\n n =11",  "31-45cm\n n =7", 
                                                      "46-60cm\n n =6",   "61-75cm\n n =4",  
                                                       "76-90cm\n n =3",    "91-105cm\n n =4", "106-120cm\n n =4"))
mcrA.total.depth$xvalue <- factor(mcrA.total.depth$xvalue,levels =rev(levels(mcrA.total.depth$xvalue)))

p15 <- ggplot(mcrA.total.site) +
  geom_col(aes(x=xvalue, y=mean),color="black",fill="grey50")+
  geom_errorbar(aes(x=xvalue, ymin=mean, ymax=mean+se))+
  xlab("")+ ylab("mcrA copies per g dry soil")+
  scale_y_continuous(limits = c(0,8.4e+8),breaks=c(0,2e+8,4e+8,6e+8,8e+8),labels= c(0,2,4,6,8),expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90, size = 15, color = "black",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.y = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),strip.text.x = element_text(size = 17, color = "black",face = "bold"),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())
p15


p16 <- ggplot(mcrA.total.depth) +
  geom_col(aes(x=xvalue, y=mean),color="black",fill="grey50")+
  geom_errorbar(aes(x=xvalue, ymin=mean, ymax=mean+se))+
  xlab("")+ ylab("mcrA copies per g dry soil")+
  scale_y_continuous(limits = c(0,8.4e+8),breaks=c(0,2e+8,4e+8,6e+8,8e+8),labels= c(0,2,4,6,8),expand = c(0,0))+
  theme(axis.text.y = element_text(size = 15, color = "black",face = "bold"),
        axis.text.x = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.x = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),strip.text.x = element_text(size = 17, color = "black",face = "bold"),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())+coord_flip()
p16

# pmoA total
library(plyr)
pmoA.total.site <- ddply(sample.data.2019, c("Locality"), summarise,
                         N    = length(pmoA_Total),
                         mean = mean(pmoA_Total),
                         sd   = sd(pmoA_Total),
                         se   = sd / sqrt(N)
)


pmoA.total.depth <-  ddply(sample.data.2019, c("Depth_proxy"), summarise,
                           N    = length(pmoA_Total),
                           mean = mean(pmoA_Total),
                           sd   = sd(pmoA_Total),
                           se   = sd / sqrt(N)
)

pmoA.total.site$xvalue <- paste(pmoA.total.site$Locality,pmoA.total.site$N,sep = ", n =")
pmoA.total.depth$xvalue <- paste(pmoA.total.depth$Depth_proxy,pmoA.total.depth$N,sep = ", n =")


pmoA.total.site$xvalue <- addline_format(pmoA.total.site$xvalue)
pmoA.total.depth$xvalue <- addline_format(pmoA.total.depth$xvalue)
pmoA.total.depth$xvalue <- factor(pmoA.total.depth$xvalue,levels = c("0-15cm\n n =17",  "16-30cm\n n =11",  "31-45cm\n n =7", 
                                                                     "46-60cm\n n =6",   "61-75cm\n n =4",  
                                                                     "76-90cm\n n =3",    "91-105cm\n n =4", "106-120cm\n n =4"))
pmoA.total.depth$xvalue <- factor(pmoA.total.depth$xvalue,levels =rev(levels(pmoA.total.depth$xvalue)))

p17 <- ggplot(pmoA.total.site) +
  geom_col(aes(x=xvalue, y=mean),color="black",fill="grey50")+
  geom_errorbar(aes(x=xvalue, ymin=mean, ymax=mean+se))+
  xlab("")+ ylab("pmoA copies per g dry soil")+
  scale_y_continuous(limits = c(0,1.7e+6),breaks=c(0,0.5e+6,1.0e+6,1.5e+6),labels= c(0,0.5,1,1.5),expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90, size = 15, color = "black",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.y = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),strip.text.x = element_text(size = 17, color = "black",face = "bold"),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())
p17


p18 <- ggplot(pmoA.total.depth) +
  geom_col(aes(x=xvalue, y=mean),color="black",fill="grey50")+
  geom_errorbar(aes(x=xvalue, ymin=mean, ymax=mean+se))+
  xlab("")+ ylab("pmoA copies per g dry soil")+
  scale_y_continuous(limits = c(0,1.7e+6),breaks=c(0,0.5e+6,1.0e+6,1.5e+6),labels= c(0,0.5,1,1.5),expand = c(0,0))+
  theme(axis.text.y = element_text(size = 15, color = "black",face = "bold"),
        axis.text.x = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.x = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),strip.text.x = element_text(size = 17, color = "black",face = "bold"),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())+coord_flip()
p18


# pmoA active
library(plyr)
pmoA.active.site <- ddply(sample.data.2019, c("Locality"), summarise,
                         N    = length(pmoA_Active),
                         mean = mean(pmoA_Active),
                         sd   = sd(pmoA_Active),
                         se   = sd / sqrt(N)
)


pmoA.active.depth <-  ddply(sample.data.2019, c("Depth_proxy"), summarise,
                           N    = length(pmoA_Active),
                           mean = mean(pmoA_Active),
                           sd   = sd(pmoA_Active),
                           se   = sd / sqrt(N)
)

pmoA.active.site$xvalue <- paste(pmoA.active.site$Locality,pmoA.active.site$N,sep = ", n =")
pmoA.active.depth$xvalue <- paste(pmoA.active.depth$Depth_proxy,pmoA.active.depth$N,sep = ", n =")


pmoA.active.site$xvalue <- addline_format(pmoA.active.site$xvalue)
pmoA.active.depth$xvalue <- addline_format(pmoA.active.depth$xvalue)
pmoA.active.depth$xvalue <- factor(pmoA.active.depth$xvalue,levels = c("0-15cm\n n =17",  "16-30cm\n n =11",  "31-45cm\n n =7", 
                                                                     "46-60cm\n n =6",   "61-75cm\n n =4",  
                                                                     "76-90cm\n n =3",    "91-105cm\n n =4", "106-120cm\n n =4"))
pmoA.active.depth$xvalue <- factor(pmoA.active.depth$xvalue,levels =rev(levels(pmoA.active.depth$xvalue)))

p19 <- ggplot(pmoA.active.site) +
  geom_col(aes(x=xvalue, y=mean),color="black",fill="grey50")+
  geom_errorbar(aes(x=xvalue, ymin=mean, ymax=mean+se))+
  xlab("")+ ylab("pmoA copies per g dry soil")+
  scale_y_continuous(limits = c(0,5.5e+6),breaks=c(0,1e+6,2e+6,3e+6,4e+6,5e+6),labels= c(0,1,2,3,4,5),expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90, size = 15, color = "black",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.y = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),strip.text.x = element_text(size = 17, color = "black",face = "bold"),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())
p19


p20 <- ggplot(pmoA.active.depth) +
  geom_col(aes(x=xvalue, y=mean),color="black",fill="grey50")+
  geom_errorbar(aes(x=xvalue, ymin=mean, ymax=mean+se))+
  xlab("")+ ylab("pmoA copies per g dry soil")+
  scale_y_continuous(limits = c(0,5.5e+6),breaks=c(0,1e+6,2e+6,3e+6,4e+6,5e+6),labels= c(0,1,2,3,4,5),expand = c(0,0))+
  theme(axis.text.y = element_text(size = 15, color = "black",face = "bold"),
        axis.text.x = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.x = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),strip.text.x = element_text(size = 17, color = "black",face = "bold"),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())+coord_flip()
p20

library(cowplot)
plot_grid(p15,p16,p17,p18,p19,p20,ncol = 2,align = "vh",axis = "tblr")

## pairwise Kruskal-Wallis
library(PMCMR)
PMCMRplus::kwAllPairsDunnTest(pmoA_Active ~ factor(Depth_proxy), data=sample.data.2019, p.adjust="fdr")


