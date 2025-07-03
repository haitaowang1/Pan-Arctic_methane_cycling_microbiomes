#creates PCoA plots of the Methanogens and Methanotrophs communities in the samples, each plot is saved once uncoloured, once coloured by site and once coloured by soil horizon. Statistical results of PERMANOVA can be found in the created data.frames 'adonis_allMethanogens' and 'adonis_allMethanotrophs'.


library(vegan)
library(tidyverse)
library(writexl)
library(car)

load("data_frames_spec")
load("data_frames_sitecat")
load("OTU_table_Methanogens")
load("OTU_table_Methanotrophs")
load("Species_table_Methanogens")
load("Species_table_Methanotrophs")
load("Site_table")

##### venn diagram
load("Site_table")
load("OTU_table_Methanogens")
load("OTU_table_Methanotrophs")

OTU_table_Methanogens <- OTU_table_Methanogens[,rownames(Site_table)]
OTU_table_Methanogens <- data.frame(t(OTU_table_Methanogens))

OTU_table_Methanogens <- OTU_table_Methanogens[which(Site_table$site=="Cherskiy"),]  # switch site names
Site_table <- Site_table[which(Site_table$site=="Cherskiy"),]

OTU_table_Methanogens.organic <- OTU_table_Methanogens[which(Site_table$category=="organic layer"),]
OTU_table_Methanogens.organic <- OTU_table_Methanogens.organic[,colSums(OTU_table_Methanogens.organic) > 0]
organic <- colnames(OTU_table_Methanogens.organic)

OTU_table_Methanogens.topsoil <- OTU_table_Methanogens[which(Site_table$category=="topsoil" | Site_table$category=="subsoil"),]
OTU_table_Methanogens.topsoil <- OTU_table_Methanogens.topsoil[,colSums(OTU_table_Methanogens.topsoil) > 0]
topsoil <- colnames(OTU_table_Methanogens.topsoil)

# OTU_table_Methanogens.subsoil <- OTU_table_Methanogens[which(Site_table$category=="subsoil"),]
# OTU_table_Methanogens.subsoil <- OTU_table_Methanogens.subsoil[,colSums(OTU_table_Methanogens.subsoil) > 0]
# subsoil <- colnames(OTU_table_Methanogens.subsoil)

OTU_table_Methanogens.cryoOM <- OTU_table_Methanogens[which(Site_table$category=="cryoOM"),]
OTU_table_Methanogens.cryoOM <- OTU_table_Methanogens.cryoOM[,colSums(OTU_table_Methanogens.cryoOM) > 0]
cryoOM <- colnames(OTU_table_Methanogens.cryoOM)

OTU_table_Methanogens.permafrost <- OTU_table_Methanogens[which(Site_table$category=="permafrost"),]
OTU_table_Methanogens.permafrost <- OTU_table_Methanogens.permafrost[,colSums(OTU_table_Methanogens.permafrost) > 0]
permafrost <- colnames(OTU_table_Methanogens.permafrost)


# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("yanlinlin82/ggvenn")

library(ggvenn)

x <- list(A=organic, B=topsoil, C=cryoOM, D=permafrost)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 1, set_name_size = 6
)

###### methanotrophs
load("Site_table")
load("OTU_table_Methanogens")
load("OTU_table_Methanotrophs")

OTU_table_Methanotrophs <- OTU_table_Methanotrophs[,rownames(Site_table)]
OTU_table_Methanotrophs <- data.frame(t(OTU_table_Methanotrophs))

OTU_table_Methanotrophs <- OTU_table_Methanotrophs[which(Site_table$site=="Herschel"),]  # switch site names
Site_table <- Site_table[which(Site_table$site=="Herschel"),]

OTU_table_Methanotrophs.organic <- OTU_table_Methanotrophs[which(Site_table$category=="organic layer"),]
OTU_table_Methanotrophs.organic <- OTU_table_Methanotrophs.organic[,colSums(OTU_table_Methanotrophs.organic) > 0]
organic <- colnames(OTU_table_Methanotrophs.organic)

OTU_table_Methanotrophs.topsoil <- OTU_table_Methanotrophs[which(Site_table$category=="topsoil" | Site_table$category=="subsoil"),]
OTU_table_Methanotrophs.topsoil <- OTU_table_Methanotrophs.topsoil[,colSums(OTU_table_Methanotrophs.topsoil) > 0]
topsoil <- colnames(OTU_table_Methanotrophs.topsoil)

# OTU_table_Methanotrophs.subsoil <- OTU_table_Methanotrophs[which(Site_table$category=="subsoil"),]
# OTU_table_Methanotrophs.subsoil <- OTU_table_Methanotrophs.subsoil[,colSums(OTU_table_Methanotrophs.subsoil) > 0]
# subsoil <- colnames(OTU_table_Methanotrophs.subsoil)

OTU_table_Methanotrophs.cryoOM <- OTU_table_Methanotrophs[which(Site_table$category=="cryoOM"),]
OTU_table_Methanotrophs.cryoOM <- OTU_table_Methanotrophs.cryoOM[,colSums(OTU_table_Methanotrophs.cryoOM) > 0]
cryoOM <- colnames(OTU_table_Methanotrophs.cryoOM)

OTU_table_Methanotrophs.permafrost <- OTU_table_Methanotrophs[which(Site_table$category=="permafrost"),]
OTU_table_Methanotrophs.permafrost <- OTU_table_Methanotrophs.permafrost[,colSums(OTU_table_Methanotrophs.permafrost) > 0]
permafrost <- colnames(OTU_table_Methanotrophs.permafrost)


library(ggvenn)

x <- list(A=organic, B=topsoil, C=cryoOM, D=permafrost)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 1, set_name_size = 6
)


#all Methanogens
load("OTU_table_Methanogens")
load("OTU_table_Methanotrophs")

OTU_table_Methanogens <- OTU_table_Methanogens[,rownames(Site_table)]

OTU_table_Methanogens <- OTU_table_Methanogens[rowSums(OTU_table_Methanogens)>0,]
OTU_table_Methanogens <- OTU_table_Methanogens[,colSums(OTU_table_Methanogens)>0]

OTU_table_Methanogens <- data.frame(t(OTU_table_Methanogens))

OTU_table_Methanogens <- decostand(OTU_table_Methanogens, method = "total", MARGIN = 1)

distance <- vegdist(OTU_table_Methanogens, distance="bray")
dist <- distance
distance <- as.matrix(distance)

library(ape)
pcoa1 <- ape::pcoa(distance)
pcoa1$values
data <- as.data.frame(pcoa1$vectors[,1:2])

#data <- as.data.frame(cmdscale(distance))

data$type <- Site_table[rownames(OTU_table_Methanogens),]$category
data$site <- Site_table[rownames(OTU_table_Methanogens),]$site

colnames(data)[1:2] <- c("V1","V2")

data$group <- paste(data$type,data$site,sep = "_")

library(plyr)  #caculate mean, SD and SE for each axis
mean.V1 <- ddply(data, .(group), summarize,  mean=mean(as.numeric(V1)), SD=sd(as.numeric(V1)), SE = sd(as.numeric(V1))/sqrt(length(V1)))
mean.V2 <- ddply(data, .(group), summarize,  mean=mean(as.numeric(V2)), SD=sd(as.numeric(V2)), SE = sd(as.numeric(V2))/sqrt(length(V2)))
names(mean.V1)<-c("group1","Mean1","SD1","SE1")
names(mean.V2)<-c("group2","Mean2","SD2","SE2")
ordinate.nmds <- cbind(mean.V1,mean.V2)

x=data.frame(t(data.frame(strsplit(ordinate.nmds$group1,"_"))))

ordinate.nmds <- cbind(ordinate.nmds,x)

# set.seed(2)
# test <- betadisper(dist,data$site,type="median")
# anova(test)
# set.seed(2)
# test <- betadisper(dist,data$type,type="median")
# anova(test)
# adonis_allMethanogens <- adonis(distance~site*type, data)$aov.tab

# library(RColorBrewer) #randomly choose colors
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

centroids1 <- aggregate(data[,1:2],list(data$site), mean)
colnames(centroids1) <- c("site","mean1","mean2")
data <- merge(data,centroids1,by="site")

centroids2 <- aggregate(data[,2:3],list(data$type), mean)
colnames(centroids2) <- c("type","mean3","mean4")
data <- merge(data,centroids2,by="type")

# col1 <- sample(col_vector, 8)

# ggplot()+
#   geom_hline(yintercept = 0,linetype="dotted",color="grey")+geom_vline(xintercept = 0,linetype="dotted",color="grey")+
#   geom_point(data=data,aes(V1,V2,color=site,fill=site),stroke=1,size=5,alpha=0.6,shape=21)+
#   geom_point(data=data,aes(x=mean1,y=mean2,fill=site),size=7,shape=23)+
#   scale_fill_manual(values = col1)+scale_color_manual(values = col1)+
#   theme(axis.title=element_text(size=12,face = "bold"),axis.text = element_text(size=10),
#         panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=1),
#         panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
#         legend.text = element_text(size = 15), legend.key.size = unit(1, "cm"), legend.key = element_rect(fill = NA),legend.title = element_text(size = 20))+
#   labs(x="PCoA1 (46.9%)",y="PCoA2 (22.2%)")
# 
# 
# col2 <- sample(col_vector, 5)
# 
# ggplot()+
#   geom_hline(yintercept = 0,linetype="dotted",color="grey")+geom_vline(xintercept = 0,linetype="dotted",color="grey")+
#   geom_point(data=data,aes(V1,V2,color=type,fill=type),stroke=1,size=5,alpha=0.6,shape=21)+
#   geom_point(data=data,aes(x=mean3,y=mean4,fill=type),size=7,shape=23)+
#   scale_fill_manual(values = col2)+scale_color_manual(values = col2)+
#   theme(axis.title=element_text(size=12,face = "bold"),axis.text = element_text(size=10),
#         panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=1),
#         panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
#         legend.text = element_text(size = 15), legend.key.size = unit(1, "cm"), legend.key = element_rect(fill = NA),legend.title = element_text(size = 20))+
#   labs(x="PCoA1 (46.9%)",y="PCoA2 (22.2%)")

adonis2(distance~type*site,data = data)


###
# col1 <- sample(col_vector, 8)
p1 <- ggplot(ordinate.nmds,aes(Mean1,Mean2))+
  geom_hline(yintercept = 0,linetype="dotted",color="black")+geom_vline(xintercept = 0,linetype="dotted",color="black")+
  geom_errorbar(aes(ymin=Mean2-SE2, ymax=Mean2+SE2),color="grey",size=1)+geom_errorbarh(aes(xmin=Mean1-SE1, xmax=Mean1+SE1),color="grey",size=1)+
  geom_point(aes(Mean1,Mean2,shape=factor(X1,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost")),
                 fill=factor(X2,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy"))),stroke=1,size=7,color="black",alpha=1)+
  scale_fill_brewer(palette = "Dark2")+scale_shape_manual(values=c(21:25))+
  guides(fill = guide_legend(title = "Site",override.aes = list(shape=21,size=4,color="white")), shape=guide_legend(title = "Layer",override.aes = list(size=4)))+
  theme(axis.title=element_text(size=12,face = "bold"),axis.text = element_text(size=10,face="bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 10), legend.key.size = unit(0.5, "cm"), legend.key = element_rect(fill = NA),legend.title = element_text(size = 12,face="bold"))+
  labs(x="PCoA1 (46.9%)",y="PCoA2 (22.2%)")

p2 <- ggplot(data)+geom_boxplot(aes(x=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")),y=V1),lwd=1,width=0.8)+
  geom_jitter(aes(x=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")),y=V1,
                              fill=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy"))),width=0.2,shape=21,size=4,alpha=0.9)+
  scale_fill_brewer(palette = "Dark2")+labs(x="",y="PCoA1")+scale_y_continuous(expand=c(0.1,0.1))+
  theme(axis.title.y=element_text(size=12,face = "bold"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,face="bold"),axis.text.x = element_text(size=12,face="bold",angle = 30,vjust = 0.6),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.position = "none")

p3 <- ggplot(data)+geom_boxplot(aes(x=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost")),y=V1),lwd=1,width=0.5)+
  geom_jitter(aes(x=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost")),y=V1,
                              shape=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))),width=0.2*5/8,size=4,alpha=0.9,fill="grey50")+
  scale_shape_manual(values=c(21:25))+
  labs(x="",y="PCoA1")+  scale_y_continuous(expand=c(0.1,0.1))+
  theme(axis.title=element_text(size=12,face = "bold"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,face="bold"),axis.text.x = element_text(size=12,face="bold",angle = 30,vjust = 0.6),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.position = "none")


p4 <- ggplot(data)+geom_boxplot(aes(x=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")),y=V2),lwd=1,width=0.8)+
  geom_jitter(aes(x=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")),y=V2,
                  fill=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy"))),width=0.2,shape=21,size=4,alpha=0.9)+
  scale_fill_brewer(palette = "Dark2")+labs(x="",y="PCoA2")+scale_y_continuous(expand=c(0.1,0.1))+
  theme(axis.title.y=element_text(size=12,face = "bold"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,face="bold"),axis.text.x = element_text(size=12,face="bold",angle = 30,vjust = 0.6),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.position = "none")

p5 <- ggplot(data)+geom_boxplot(aes(x=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost")),y=V2),lwd=1,width=0.5)+
  geom_jitter(aes(x=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost")),y=V2,
                  shape=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))),width=0.2*5/8,size=4,alpha=0.9,fill="grey50")+
  scale_shape_manual(values=c(21:25))+
  labs(x="",y="PCoA2")+  scale_y_continuous(expand=c(0.1,0.1))+
  theme(axis.title=element_text(size=12,face = "bold"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,face="bold"),axis.text.x = element_text(size=12,face="bold",angle = 30,vjust = 0.6),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.position = "none")

library(cowplot)
middle_col <- plot_grid(p2,p3,ncol = 1,align = "v",axis = "l",rel_heights = c(1,1))
right_col <- plot_grid(p4,p5,ncol = 1,align = "v",axis = "l",rel_heights = c(1,1))
plot_grid(p1,middle_col,right_col,ncol = 3,align = "v",axis = "tblr",rel_widths = c(1,1))

library(PMCMR)
posthoc.kruskal.dunn.test(V2 ~ factor(site),data = data, p.adjust = "fdr")

#####################
# col2 <- sample(col_vector, 5)
ggplot(ordinate.nmds,aes(Mean1,Mean2))+
  geom_hline(yintercept = 0,linetype="dotted",color="black")+geom_vline(xintercept = 0,linetype="dotted",color="black")+
  geom_errorbar(aes(ymin=Mean2-SE2, ymax=Mean2+SE2),color="grey",size=1)+geom_errorbarh(aes(xmin=Mean1-SE1, xmax=Mean1+SE1),color="grey",size=1)+
  geom_point(aes(Mean1,Mean2,shape=X1,fill=X1),stroke=1,size=7,color="black",alpha=1)+
  scale_fill_manual(values=col2)+scale_shape_manual(values=c(21:25))+
  guides(fill = guide_legend(title = "group",override.aes = list(shape=c(21:25),color=col2)), shape=guide_legend(title = "group"))+
  theme(axis.title=element_text(size=12,face = "bold"),axis.text = element_text(size=10),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15), legend.key.size = unit(1, "cm"), legend.key = element_rect(fill = NA),legend.title = element_text(size = 20))+
  labs(x="PCoA1 (46.9%)",y="PCoA2 (22.2%)")

#all Methanotrophs

OTU_table_Methanotrophs <- OTU_table_Methanotrophs[,rownames(Site_table)]

OTU_table_Methanotrophs <- OTU_table_Methanotrophs[rowSums(OTU_table_Methanotrophs)>0,]
OTU_table_Methanotrophs <- OTU_table_Methanotrophs[,colSums(OTU_table_Methanotrophs)>0]

OTU_table_Methanotrophs <- data.frame(t(OTU_table_Methanotrophs))

OTU_table_Methanotrophs <- decostand(OTU_table_Methanotrophs, method = "total", MARGIN = 1)


distance2 <- vegdist(OTU_table_Methanotrophs, distance="bray")
distance2 <- as.matrix(distance2)

pcoa2 <- ape::pcoa(distance2)
pcoa2$values
data2 <- as.data.frame(pcoa2$vectors[,1:2])

data2$site <- Site_table[rownames(OTU_table_Methanotrophs),]$site
data2$type <- Site_table[rownames(OTU_table_Methanotrophs),]$category

# set.seed(2)
# test <- betadisper(dist,data$site,type="median")
# anova(test)
# set.seed(2)
# test <- betadisper(dist,data$type,type="median")
# anova(test)
# adonis_allMethanotrophs <- adonis(distance~site*type, data)$aov.tab
colnames(data2)[1:2] <- c("V1","V2")
data2$group <- paste(data2$type,data2$site,sep = "_")

library(plyr)  #caculate mean, SD and SE for each axis
mean.V1 <- ddply(data2, .(group), summarize,  mean=mean(as.numeric(V1)), SD=sd(as.numeric(V1)), SE = sd(as.numeric(V1))/sqrt(length(V1)))
mean.V2 <- ddply(data2, .(group), summarize,  mean=mean(as.numeric(V2)), SD=sd(as.numeric(V2)), SE = sd(as.numeric(V2))/sqrt(length(V2)))
names(mean.V1)<-c("group1","Mean1","SD1","SE1")
names(mean.V2)<-c("group2","Mean2","SD2","SE2")
ordinate.nmds2 <- cbind(mean.V1,mean.V2)

x=data.frame(t(data.frame(strsplit(ordinate.nmds2$group1,"_"))))

ordinate.nmds2 <- cbind(ordinate.nmds2,x)


adonis2(distance2~type*site,data = data2)


###
# col1 <- sample(col_vector, 8)
p6 <- ggplot(ordinate.nmds2,aes(Mean1,Mean2))+
  geom_hline(yintercept = 0,linetype="dotted",color="black")+geom_vline(xintercept = 0,linetype="dotted",color="black")+
  geom_errorbar(aes(ymin=Mean2-SE2, ymax=Mean2+SE2),color="grey",size=1)+geom_errorbarh(aes(xmin=Mean1-SE1, xmax=Mean1+SE1),color="grey",size=1)+
  geom_point(aes(Mean1,Mean2,shape=factor(X1,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost")),
                 fill=factor(X2,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy"))),stroke=1,size=7,color="black",alpha=1)+
  scale_fill_brewer(palette = "Dark2")+scale_shape_manual(values=c(21:25))+
  guides(fill = guide_legend(title = "Site",override.aes = list(shape=21,size=4,color="white")), shape=guide_legend(title = "Layer",override.aes = list(size=4)))+
  theme(axis.title=element_text(size=12,face = "bold"),axis.text = element_text(size=10,face="bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 10), legend.key.size = unit(0.5, "cm"), legend.key = element_rect(fill = NA),legend.title = element_text(size = 12,face="bold"))+
  labs(x="PCoA1 (23.3%)",y="PCoA2 (18.5%)")

p7 <- ggplot(data2)+ geom_boxplot(aes(x=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")),y=V1),lwd=1,width=0.8)+
  geom_jitter(aes(x=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")),y=V1,
                                   fill=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy"))),width=0.2,shape=21,size=4,alpha=0.9)+
  scale_fill_brewer(palette = "Dark2")+labs(x="",y="PCoA1")+scale_y_continuous(expand=c(0.1,0.1))+
  theme(axis.title.y=element_text(size=12,face = "bold"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,face="bold"),axis.text.x = element_text(size=12,face="bold",angle = 30,vjust = 0.6),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.position = "none")

p8 <- ggplot(data2)+  geom_boxplot(aes(x=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost")),y=V1),lwd=1,width=0.5)+
  geom_jitter(aes(x=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost")),y=V1,
                                   shape=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))),width=0.2*5/8,size=4,alpha=0.9,fill="grey50")+
  scale_shape_manual(values=c(21:25))+
  labs(x="",y="PCoA1")+  scale_y_continuous(expand=c(0.1,0.1))+
  theme(axis.title=element_text(size=12,face = "bold"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,face="bold"),axis.text.x = element_text(size=12,face="bold",angle = 30,vjust = 0.6),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.position = "none")

p9 <- ggplot(data2)+ geom_boxplot(aes(x=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")),y=V2),lwd=1,width=0.8)+
  geom_jitter(aes(x=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")),y=V2,
                  fill=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy"))),width=0.2,shape=21,size=4,alpha=0.9)+
  scale_fill_brewer(palette = "Dark2")+labs(x="",y="PCoA2")+scale_y_continuous(expand=c(0.1,0.1))+
  theme(axis.title.y=element_text(size=12,face = "bold"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,face="bold"),axis.text.x = element_text(size=12,face="bold",angle = 30,vjust = 0.6),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.position = "none")

p10 <- ggplot(data2)+  geom_boxplot(aes(x=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost")),y=V2),lwd=1,width=0.5)+
  geom_jitter(aes(x=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost")),y=V2,
                  shape=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))),width=0.2*5/8,size=4,alpha=0.9,fill="grey50")+
  scale_shape_manual(values=c(21:25))+
  labs(x="",y="PCoA2")+  scale_y_continuous(expand=c(0.1,0.1))+
  theme(axis.title=element_text(size=12,face = "bold"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,face="bold"),axis.text.x = element_text(size=12,face="bold",angle = 30,vjust = 0.6),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.position = "none")

aggregate(data2$V1,list(data2$type),mean)

library(cowplot)
middle_col2 <- plot_grid(p7,p8,ncol = 1,align = "v",axis = "l",rel_heights = c(1,1))
right_col2 <- plot_grid(p9,p10,ncol = 1,align = "v",axis = "l",rel_heights = c(1,1))
plot_grid(p6,middle_col2,right_col2,ncol = 3,align = "v",axis = "tblr")

library(PMCMR)
posthoc.kruskal.dunn.test(V2 ~ factor(type),data = data2, p.adjust = "fdr")


col2 <- sample(col_vector, 5)
ggplot(ordinate.nmds2,aes(Mean1,Mean2))+
  geom_hline(yintercept = 0,linetype="dotted",color="black")+geom_vline(xintercept = 0,linetype="dotted",color="black")+
  geom_errorbar(aes(ymin=Mean2-SE2, ymax=Mean2+SE2),color="grey",size=1)+geom_errorbarh(aes(xmin=Mean1-SE1, xmax=Mean1+SE1),color="grey",size=1)+
  geom_point(aes(Mean1,Mean2,shape=X1,fill=X1),stroke=1,size=7,color="black",alpha=1)+scale_fill_manual(values=col2)+scale_shape_manual(values=c(21:25))+
  guides(fill = guide_legend(title = "group",override.aes = list(shape=c(21:25),color=col2)), shape=guide_legend(title = "group"))+
  theme(axis.title=element_text(size=12,face = "bold"),axis.text = element_text(size=10),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank(),
        legend.text = element_text(size = 15), legend.key.size = unit(1, "cm"), legend.key = element_rect(fill = NA),legend.title = element_text(size = 20))+
  labs(x="PCoA1 (23.3%)",y="PCoA2 (18.5%)")


# mantel test
load("env_data")
env_data1 <- env_data[!is.na(env_data$Watercontent),]
env_data1 <- env_data1[!is.na(env_data1$pH),]
intersect(rownames(env_data1),rownames(OTU_table_Methanogens))

env_data2 <- env_data1[intersect(rownames(env_data1),rownames(OTU_table_Methanogens)),]
OTU_table_Methanogens2 <- OTU_table_Methanogens[intersect(rownames(env_data1),rownames(OTU_table_Methanogens)),]

distance.MG <-vegdist(OTU_table_Methanogens2, distance="bray")
distance.water = dist(env_data2$Watercontent, method = "euclidean")
distance.pH = dist(env_data2$pH, method = "euclidean")

mantel(distance.MG, distance.pH, method = "spearman", permutations = 999, na.rm = TRUE)

intersect(rownames(env_data1),rownames(OTU_table_Methanotrophs))
env_data3 <- env_data1[intersect(rownames(env_data1),rownames(OTU_table_Methanotrophs)),]
OTU_table_Methanotrophs2 <- OTU_table_Methanotrophs[intersect(rownames(env_data1),rownames(OTU_table_Methanotrophs)),]

distance.MT <-vegdist(OTU_table_Methanotrophs2, distance="bray")
distance.water = dist(env_data3$Watercontent, method = "euclidean")
distance.pH = dist(env_data3$pH, method = "euclidean")

mantel(distance.MT, distance.water, method = "spearman", permutations = 999, na.rm = TRUE)


