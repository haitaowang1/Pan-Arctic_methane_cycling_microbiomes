library(ggplot2)

load("../Erik_thesis/OTU_table_Methanogens_cat")
load("Site_table_new")
load("../Erik_thesis/col_taxa")

OTU_table_Methanogens_cat <- data.frame(t(OTU_table_Methanogens_cat))
OTU_table_Methanogens_cat <- OTU_table_Methanogens_cat[rownames(Site_table),]

OTU_table_Methanogens_site <- aggregate(OTU_table_Methanogens_cat,list(Site_table$site),mean)
OTU_table_Methanogens_site <- reshape2::melt(OTU_table_Methanogens_site,id.vars="Group.1")

library(RColorBrewer) #randomly choose colors
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# col_taxa=sample(col_vector, 7)
# save(col_taxa,file="col_taxa")

ggplot(OTU_table_Methanogens_site,aes(x=factor(Group.1,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")),y=value*100,fill=variable)) + 
  geom_bar(stat="identity",position="stack",color="black")+
  guides(fill=guide_legend(nrow = 21,title="Taxonomy"))+ylab("Relative abundance (%)")+xlab("")+
  scale_fill_manual(values=col_taxa)+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.title = element_text(size = 20,face = "bold"),axis.text.y = element_text(size = 15),axis.text.x = element_text(size = 15,angle=90,face = "bold"),
        legend.text = element_text(size = 10),legend.title = element_text(size = 12),
        panel.background = element_rect(fill = NA),panel.border = element_rect(linetype = "solid", colour = "black",fill=NA))

OTU_table_Methanogens_horizon <- aggregate(OTU_table_Methanogens_cat,list(Site_table$category),mean)
OTU_table_Methanogens_horizon <- reshape2::melt(OTU_table_Methanogens_horizon,id.vars="Group.1")


ggplot(OTU_table_Methanogens_horizon,aes(x=factor(Group.1,levels = c("permafrost","cryoOM","subsoil","topsoil","organic layer")),y=value*100,fill=variable)) + 
  geom_bar(stat="identity",position="stack",color="black")+
  guides(fill=guide_legend(nrow = 21,title="Taxonomy"))+ylab("Relative abundance (%)")+xlab("")+
  scale_fill_manual(values=col_taxa)+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.title = element_text(size = 20,face = "bold"),axis.text.y = element_text(size = 15),axis.text.x = element_text(size = 15,angle=90,face = "bold"),
        legend.text = element_text(size = 10),legend.title = element_text(size = 12),
        panel.background = element_rect(fill = NA),panel.border = element_rect(linetype = "solid", colour = "black",fill=NA))+coord_flip()

################################
load("../Erik_thesis/OTU_table_Methanotrophs_cat2")

OTU_table_Methanotroph_cat <- data.frame(t(OTU_table_Methanotrophs_cat2))
OTU_table_Methanotroph_cat <- OTU_table_Methanotroph_cat[rownames(Site_table),]

OTU_table_Methanotroph_site <- aggregate(OTU_table_Methanotroph_cat,list(Site_table$site),mean)
OTU_table_Methanotroph_site <- reshape2::melt(OTU_table_Methanotroph_site,id.vars="Group.1")

unique(OTU_table_Methanotroph_site$variable)
col_taxa2=sample(col_vector, 14)

ggplot(OTU_table_Methanotroph_site,aes(x=factor(Group.1,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")),y=value*100,fill=variable)) + 
  geom_bar(stat="identity",position="stack",color="black")+
  guides(fill=guide_legend(nrow = 21,title="Taxonomy"))+ylab("Relative abundance (%)")+xlab("")+
  scale_fill_manual(values=col_taxa2,labels=rownames(OTU_table_Methanotrophs_cat2))+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.title = element_text(size = 20,face = "bold"),axis.text.y = element_text(size = 15),axis.text.x = element_text(size = 15,angle=90,face = "bold"),
        legend.text = element_text(size = 10),legend.title = element_text(size = 12),
        panel.background = element_rect(fill = NA),panel.border = element_rect(linetype = "solid", colour = "black",fill=NA))

OTU_table_Methanotroph_horizon <- aggregate(OTU_table_Methanotroph_cat,list(Site_table$category),mean)
OTU_table_Methanotroph_horizon <- reshape2::melt(OTU_table_Methanotroph_horizon,id.vars="Group.1")


ggplot(OTU_table_Methanotroph_horizon,aes(x=factor(Group.1,levels = c("permafrost","cryoOM","subsoil","topsoil","organic layer")),y=value*100,fill=variable)) + 
  geom_bar(stat="identity",position="stack",color="black")+
  guides(fill=guide_legend(nrow = 21,title="Taxonomy"))+ylab("Relative abundance (%)")+xlab("")+
  scale_fill_manual(values=col_taxa2,labels=rownames(OTU_table_Methanotrophs_cat2))+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.title = element_text(size = 20,face = "bold"),axis.text.y = element_text(size = 15),axis.text.x = element_text(size = 15,angle=90,face = "bold"),
        legend.text = element_text(size = 10),legend.title = element_text(size = 12),
        panel.background = element_rect(fill = NA),panel.border = element_rect(linetype = "solid", colour = "black",fill=NA))+coord_flip()

######################## specialist and generalist ASVs
load("Site_table_new")
load("../Erik_thesis/OTU_table_Methanogens")
OTU_table_Methanogens_site2 <- data.frame(t(OTU_table_Methanogens))
OTU_table_Methanogens_site2 <- OTU_table_Methanogens_site2[rownames(Site_table),]
OTU_table_Methanogens_site2 <-  aggregate(OTU_table_Methanogens_site2,list(Site_table$site),mean)
OTU_table_Methanogens_horizon2 <- data.frame(t(OTU_table_Methanogens))
OTU_table_Methanogens_horizon2 <- OTU_table_Methanogens_horizon2[rownames(Site_table),]
OTU_table_Methanogens_horizon2 <-  aggregate(OTU_table_Methanogens_horizon2,list(Site_table$category),mean)

OTU_table_Methanogens_bubble <- rbind(OTU_table_Methanogens_site2,OTU_table_Methanogens_horizon2)
OTU_table_Methanogens_bubble <- reshape2::melt(OTU_table_Methanogens_bubble,id.vars="Group.1")

OTU_table_Methanogens_bubble$type <- "Site"
OTU_table_Methanogens_bubble$type[grep("subsoil|topsoil|cryoOM|organic|permafrost",OTU_table_Methanogens_bubble$Group.1)] <- "Horizon"

load("../Erik_thesis/Species_table_Methanogens")
Species_table_Methanogens <- read.csv("revision/tax_silva138.2_MG_zotu.csv",row.names = 1)
Species_table_Methanogens$group <- paste(Species_table_Methanogens$Family,Species_table_Methanogens$Genus2,sep = ":")

OTU_table_Methanogens_bubble$taxa <- OTU_table_Methanogens_bubble$variable
for (i in 1:22) {
  OTU_table_Methanogens_bubble$taxa <- gsub(rownames(Species_table_Methanogens)[i],Species_table_Methanogens$group[i],OTU_table_Methanogens_bubble$taxa)
}

OTU_table_Methanogens_bubble$taxa3 <- OTU_table_Methanogens_bubble$variable
for (i in 1:22) {
  OTU_table_Methanogens_bubble$taxa3 <- gsub(rownames(Species_table_Methanogens)[i],Species_table_Methanogens$Family[i],  OTU_table_Methanogens_bubble$taxa3)
}

OTU_table_Methanogens_bubble$taxa4 <- OTU_table_Methanogens_bubble$variable
for (i in 1:22) {
  OTU_table_Methanogens_bubble$taxa4 <- gsub(rownames(Species_table_Methanogens)[i],Species_table_Methanogens$Genus2[i],  OTU_table_Methanogens_bubble$taxa4)
}

OTU_table_Methanogens_bubble$variable <- gsub("Zotu","ZOTU",OTU_table_Methanogens_bubble$variable)
OTU_table_Methanogens_bubble$taxa2 <- paste(OTU_table_Methanogens_bubble$taxa4,OTU_table_Methanogens_bubble$variable,sep = "-")

## fix some taxonomy after BLAST
#OTU_table_Methanogens_bubble$taxa3 <- gsub("Rice Cluster II","Ca. Methanoflorentaceae",OTU_table_Methanogens_bubble$taxa3)
#OTU_table_Methanogens_bubble$taxa2 <- gsub("Rice Cluster II","Ca. Methanoflorentaceae",OTU_table_Methanogens_bubble$taxa2)

x_order2 <- data.frame(sort(unique(OTU_table_Methanogens_bubble$taxa2)),stringsAsFactors = F)
x_order2$ranks <- rownames(x_order2)
x_order2 <- x_order2[c(21,22,6:1,7,16,11:8,18:17,15:12,20:19),]

split_x_order2 <- do.call(rbind, strsplit(as.character(x_order2$sort.unique.OTU_table_Methanogens_bubble.taxa2..), "-"))
x_order2$taxa <- split_x_order2[,1]



p1 <- ggplot(OTU_table_Methanogens_bubble,aes(x=factor(Group.1,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy","organic layer","topsoil","subsoil","cryoOM","permafrost")),
                                        y=taxa2,size=ifelse(value==0, NA, value)*100))+geom_point(aes(fill=taxa3),color="black",shape=21,stroke=0.7)+
  scale_y_discrete(position = "right",limits=x_order2[,1],labels=x_order2$taxa)+scale_fill_manual(values = col_taxa)+
  scale_size(name="RA (%)",range = c(1,15),limits=c(0,0.55),breaks = c(0.01,0.1,0.5))+
  guides(fill=guide_legend(title="Family",override.aes = list(size=3.5)),size=guide_legend(override.aes = list(shape=19)))+
  theme(axis.title = element_blank(),axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12,angle=90,face = "bold",hjust = 1,vjust = 0.5),legend.position = "left",
        strip.text = element_text(size = 12,face = "bold"),legend.background = element_rect(fill = NA),legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 10),legend.title = element_text(size = 12,face = "bold"),panel.grid.major = element_line(linetype = "dotted", colour = "grey90"),
        panel.background = element_rect(fill = NA), panel.border = element_rect(linetype = "solid", colour = "black",fill=NA,size = 0.5))+
  facet_grid(.~factor(type,levels = c("Site","Horizon")),scales = "free",space = "free")
p1
#############################################################
load("../Erik_thesis/OTU_table_Methanotrophs")
OTU_table_Methanotrophs_site2 <- data.frame(t(OTU_table_Methanotrophs))
OTU_table_Methanotrophs_site2 <- OTU_table_Methanotrophs_site2[rownames(Site_table),]
OTU_table_Methanotrophs_site2 <-  aggregate(OTU_table_Methanotrophs_site2,list(Site_table$site),mean)
OTU_table_Methanotrophs_horizon2 <- data.frame(t(OTU_table_Methanotrophs))
OTU_table_Methanotrophs_horizon2 <- OTU_table_Methanotrophs_horizon2[rownames(Site_table),]
OTU_table_Methanotrophs_horizon2 <-  aggregate(OTU_table_Methanotrophs_horizon2,list(Site_table$category),mean)

## to calculate how much the generalists account for
OTU_table_Methanotrophs_site3 <- OTU_table_Methanotrophs_site2[,-1]
rownames(OTU_table_Methanotrophs_site3) <- OTU_table_Methanotrophs_site2[,1]
OTU_table_Methanotrophs_site3 <- data.frame(t(OTU_table_Methanotrophs_site3))
OTU_table_Methanotrophs_site3 <- rbind(OTU_table_Methanotrophs_site3,colSums(OTU_table_Methanotrophs_site3))
OTU_table_Methanotrophs_site3 <- rbind(OTU_table_Methanotrophs_site3,colSums(OTU_table_Methanotrophs_site3[c("Zotu366","Zotu1727","Zotu1519","Zotu1239","Zotu1199","Zotu827"),]))
OTU_table_Methanotrophs_site3 <- rbind(OTU_table_Methanotrophs_site3,OTU_table_Methanotrophs_site3[28,]/OTU_table_Methanotrophs_site3[27,]*100)
OTU_table_Methanotrophs_site3[29,] 
rowMeans(OTU_table_Methanotrophs_site3[29,])
##############################################

OTU_table_Methanotrophs_bubble <- rbind(OTU_table_Methanotrophs_site2,OTU_table_Methanotrophs_horizon2)
OTU_table_Methanotrophs_bubble <- reshape2::melt(OTU_table_Methanotrophs_bubble,id.vars="Group.1")

OTU_table_Methanotrophs_bubble$type <- "Site"
OTU_table_Methanotrophs_bubble$type[grep("subsoil|topsoil|cryoOM|organic|permafrost",OTU_table_Methanotrophs_bubble$Group.1)] <- "Horizon"

#load("../Erik_thesis/Species_table_Methanotrophs")
Species_table_Methanotrophs <- read.csv("revision/tax_silva138.2_MT_zotu.csv",row.names = 1)
Species_table_Methanotrophs$group <- paste(Species_table_Methanotrophs$Family,Species_table_Methanotrophs$Genus2,sep = ":")

OTU_table_Methanotrophs_bubble$taxa <- OTU_table_Methanotrophs_bubble$variable
for (i in 1:26) {
  OTU_table_Methanotrophs_bubble$taxa <- gsub(rownames(Species_table_Methanotrophs)[i],Species_table_Methanotrophs$group[i],OTU_table_Methanotrophs_bubble$taxa)
}

OTU_table_Methanotrophs_bubble$taxa3 <- OTU_table_Methanotrophs_bubble$variable
for (i in 1:26) {
  OTU_table_Methanotrophs_bubble$taxa3 <- gsub(rownames(Species_table_Methanotrophs)[i],Species_table_Methanotrophs$Family[i],  OTU_table_Methanotrophs_bubble$taxa3)
}

OTU_table_Methanotrophs_bubble$taxa4 <- OTU_table_Methanotrophs_bubble$variable
for (i in 1:26) {
  OTU_table_Methanotrophs_bubble$taxa4 <- gsub(rownames(Species_table_Methanotrophs)[i],Species_table_Methanotrophs$Genus2[i],  OTU_table_Methanotrophs_bubble$taxa4)
}

OTU_table_Methanotrophs_bubble$variable <- gsub("Zotu","ZOTU",OTU_table_Methanotrophs_bubble$variable)
OTU_table_Methanotrophs_bubble$taxa2 <- paste(OTU_table_Methanotrophs_bubble$taxa4,OTU_table_Methanotrophs_bubble$variable,sep = "-")

## fix some taxonomy after BLAST
# OTU_table_Methanotrophs_bubble$taxa3 <- gsub("Other Methylococcales","Methylococcaceae",OTU_table_Methanotrophs_bubble$taxa3)
# OTU_table_Methanotrophs_bubble$taxa2 <- gsub("NA-ZOTU1253","Methylotuvimicrobium-ZOTU1253",OTU_table_Methanotrophs_bubble$taxa2)
# OTU_table_Methanotrophs_bubble$taxa2 <- gsub("NA-ZOTU5478","Methylotuvimicrobium-ZOTU5478",OTU_table_Methanotrophs_bubble$taxa2)
# OTU_table_Methanotrophs_bubble$taxa2 <- gsub("NA","Methylobacter",OTU_table_Methanotrophs_bubble$taxa2)
# OTU_table_Methanotrophs_bubble$taxa2 <- gsub("Methylobacter-ZOTU13655","NA-ZOTU13655",OTU_table_Methanotrophs_bubble$taxa2)
# 
# OTU_table_Methanotrophs_bubble$taxa3 <- gsub("pLW-20","Methylococcaceae",OTU_table_Methanotrophs_bubble$taxa3)

x_order <- data.frame(sort(unique(OTU_table_Methanotrophs_bubble$taxa2)),stringsAsFactors = F)
x_order$ranks <- rownames(x_order)
x_order <- x_order[c(3:1,14:13,26:22,12:4,21:15),]

split_x_order <- do.call(rbind, strsplit(as.character(x_order$sort.unique.OTU_table_Methanotrophs_bubble.taxa2..), "-"))
x_order$taxa <- split_x_order[,1]

p2 <- ggplot(OTU_table_Methanotrophs_bubble,aes(x=factor(Group.1,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy","organic layer","topsoil","subsoil","cryoOM","permafrost")),
                                        y=taxa2,size=ifelse(value==0, NA, value)*100))+geom_point(aes(fill=taxa3),color="black",shape=21,stroke=0.7)+
  scale_y_discrete(position = "right",limits=x_order[,1],labels=x_order$taxa)+scale_fill_manual(values = col_taxa)+
  scale_size(name="RA (%)",range = c(1,15),limits=c(0,0.1),breaks = c(0.01,0.05,0.1))+
  guides(fill=guide_legend(title="Family",order=2,override.aes = list(size=3.5)),size=guide_legend(order=1,override.aes = list(shape=19)))+
  theme(axis.title = element_blank(),axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12,angle=90,face = "bold",hjust = 1,vjust = 0.5),legend.position = "left",
        strip.text = element_text(size = 12,face = "bold"),legend.background = element_rect(fill = NA),legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 10),legend.title = element_text(size = 12,face = "bold"),panel.grid.major = element_line(linetype = "dotted", colour = "grey90"),
        panel.background = element_rect(fill = NA), panel.border = element_rect(linetype = "solid", colour = "black",fill=NA,size = 0.5))+
  facet_grid(.~factor(type,levels = c("Site","Horizon")),scales = "free",space = "free")
p2
library(cowplot)
plot_grid(p1,p2,ncol = 1,align = "v",axis = "tblr")

### caculating proportions
Species_table_Methanotrophs <- Species_table_Methanotrophs[rownames(OTU_table_Methanotrophs),]
OTU_table_Methanotrophs_cat3 <- aggregate(OTU_table_Methanotrophs, list(Species_table_Methanotrophs$Genus),sum)
rownames(OTU_table_Methanotrophs_cat3) <- OTU_table_Methanotrophs_cat3$Group.1
OTU_table_Methanotrophs_cat3 <- data.frame(t(OTU_table_Methanotrophs_cat3[,-1]))
OTU_table_Methanotrophs_cat3 <- aggregate(OTU_table_Methanotrophs_cat3, list(Site_table$site),mean)
OTU_table_Methanotrophs_cat3$Bacter.ratio <- OTU_table_Methanotrophs_cat3$Methylobacter/rowSums(OTU_table_Methanotrophs_cat3[,-1])
mean(OTU_table_Methanotrophs_cat3$Bacter.ratio)

