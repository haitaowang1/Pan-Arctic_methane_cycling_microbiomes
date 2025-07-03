#correlates Methanogens and Methanotrophs abundances with each other and with pH and Water content data, on a per-sample basis.


library(vegan)
library(tidyverse)
library(writexl)
library(car)

load("../Erik_thesis/OTU_table_Methanogens")
load("../Erik_thesis/OTU_table_Methanotrophs")
load("../Erik_thesis/OTU_table_Methanogens_cat")
load("../Erik_thesis/OTU_table_Methanotrophs_cat")
load("Site_table_new")
load("env_data")
Site_table$rownames <- rownames(Site_table)
Site_table$category <- as.character(Site_table$category)
Site_table$site <- as.character(Site_table$site)

OTU_table_Methanogens <- as.data.frame(t(OTU_table_Methanogens))
OTU_table_Methanogens <- OTU_table_Methanogens[rownames(Site_table),]
OTU_table_Methanotrophs <- as.data.frame(t(OTU_table_Methanotrophs))
OTU_table_Methanotrophs <- OTU_table_Methanotrophs[rownames(Site_table),]
OTU_table_Methanogens_cat <- decostand(OTU_table_Methanogens_cat, method = "total", MARGIN = 2)
OTU_table_Methanogens_cat <- as.data.frame(t(OTU_table_Methanogens_cat))
OTU_table_Methanogens_cat <- OTU_table_Methanogens_cat[rownames(Site_table),]
OTU_table_Methanotrophs_cat <- decostand(OTU_table_Methanotrophs_cat, method = "total", MARGIN = 2)
OTU_table_Methanotrophs_cat <- as.data.frame(t(OTU_table_Methanotrophs_cat))
OTU_table_Methanotrophs_cat <- OTU_table_Methanotrophs_cat[rownames(Site_table),]

OTU_table_Methanogens$sum <- rowSums(OTU_table_Methanogens)
OTU_table_Methanotrophs$sum <- rowSums(OTU_table_Methanotrophs)

Data_frame_corr <- data.frame(sum_Methanogens = OTU_table_Methanogens$sum, sum_Methanotrophs = OTU_table_Methanotrophs$sum, pH = rep(0), WC = rep(0), site = rep(0), type = rep(0))
Data_frame_corr <- cbind(Data_frame_corr, OTU_table_Methanogens_cat, OTU_table_Methanotrophs_cat)

rownames(Data_frame_corr) <- rownames(OTU_table_Methanogens)

for (i in c(1:nrow(Data_frame_corr))) {
  True_False <- rownames(Site_table) == rownames(Data_frame_corr)[i]
  if(any(True_False) == TRUE){
    Data_frame_corr$site[i] <- Site_table$site[True_False == TRUE]
  }
}
for (i in c(1:nrow(Data_frame_corr))) {
  True_False <- rownames(Site_table) == rownames(Data_frame_corr)[i]
  if(any(True_False) == TRUE){
    Data_frame_corr$type[i] <- Site_table$category[True_False == TRUE]
  }
}
for (i in c(1:nrow(Data_frame_corr))) {
  True_False <- rownames(env_data) == rownames(Data_frame_corr)[i]
  if(any(True_False) == TRUE){
    Data_frame_corr$pH[i] <- env_data$pH[True_False == TRUE]
  }
}
for (i in c(1:nrow(Data_frame_corr))) {
  True_False <- rownames(env_data) == rownames(Data_frame_corr)[i]
  if(any(True_False) == TRUE){
    Data_frame_corr$WC[i] <- env_data$Watercontent[True_False == TRUE]
  }
}

# between WC and MG/MT
library(car)
library(ggplot2)

Data_frame_corr2 <- Data_frame_corr[!is.na(Data_frame_corr$WC),]
which(Data_frame_corr2$WC==0)
which(Data_frame_corr2$sum_Methanogens==0)
Data_frame_corr2 <- Data_frame_corr2[which(Data_frame_corr2$WC>0),]
#Data_frame_corr2 <- Data_frame_corr2[which(Data_frame_corr2$sum_Methanogens>0),]

shapiro.test(log10(Data_frame_corr2$sum_Methanogens+0.001))
shapiro.test(Data_frame_corr2$WC)

# pred <- predict(lm(WC ~ log10(sum_Methanogens+0.001), Data_frame_corr2),
#                 se.fit = TRUE, interval = "confidence")
# limits <- as.data.frame(pred$fit)
# 
# p1 <- ggplot(Data_frame_corr2,aes(x=log10(sum_Methanogens+0.001), y=WC))+
#   geom_point(aes(color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))),size=2,stroke=1.5,shape=19)+
#   geom_smooth(method = "lm",color="black",level=0.95,size=2,se=F)+
#   geom_line(aes(x = log10(sum_Methanogens+0.001), y = limits$lwr), linetype = 2,size=1.5) +
#   geom_line(aes(x = log10(sum_Methanogens+0.001), y = limits$upr), linetype = 2,size=1.5) +
#   scale_color_brewer(palette="Dark2")+
#   guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)))+
#   theme(axis.title=element_text(size=20,face = "bold"),axis.text= element_text(size=12,face = "bold"),
#         panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
#         axis.line.y.left   = element_line(color = 'black',size=1),
#         panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
#   labs(x="log10(MG RA)",y="Water content (%)")

p1 <- ggplot(Data_frame_corr2,aes(x=log10(sum_Methanogens+0.001), y=WC,
                                  color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost")),
                                  fill=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))))+
  geom_point(size=2,stroke=1.5,shape=19)+
  geom_smooth(method = "lm",size=2,se=T,level=0.95,alpha=0.2)+
  scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+
  guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)),fill="none")+
  theme(axis.title=element_text(size=20,face = "bold"),axis.text= element_text(size=12,face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
  labs(x="log10(MG RA)",y="Water content (%)")

p1

# cor.test(Data_frame_corr2$WC, log10(Data_frame_corr2$sum_Methanogens+0.001),method = "pearson")
# cor.test(Data_frame_corr2$WC, log10(Data_frame_corr2$sum_Methanogens+0.001),method = "spearman")

cor.test(Data_frame_corr2$WC[Data_frame_corr2$type=="organic layer"], log10(Data_frame_corr2$sum_Methanogens[Data_frame_corr2$type=="organic layer"]+0.001),method = "pearson")
cor.test(Data_frame_corr2$WC[Data_frame_corr2$type=="topsoil"], log10(Data_frame_corr2$sum_Methanogens[Data_frame_corr2$type=="topsoil"]+0.001),method = "pearson")
cor.test(Data_frame_corr2$WC[Data_frame_corr2$type=="subsoil"], log10(Data_frame_corr2$sum_Methanogens[Data_frame_corr2$type=="subsoil"]+0.001),method = "pearson")
cor.test(Data_frame_corr2$WC[Data_frame_corr2$type=="cryoOM"], log10(Data_frame_corr2$sum_Methanogens[Data_frame_corr2$type=="cryoOM"]+0.001),method = "pearson")
cor.test(Data_frame_corr2$WC[Data_frame_corr2$type=="permafrost"], log10(Data_frame_corr2$sum_Methanogens[Data_frame_corr2$type=="permafrost"]+0.001),method = "pearson")

cor.test(Data_frame_corr2$WC[Data_frame_corr2$type=="organic layer"], log10(Data_frame_corr2$sum_Methanogens[Data_frame_corr2$type=="organic layer"]+0.001),method = "spearman")
cor.test(Data_frame_corr2$WC[Data_frame_corr2$type=="topsoil"], log10(Data_frame_corr2$sum_Methanogens[Data_frame_corr2$type=="topsoil"]+0.001),method = "spearman")
cor.test(Data_frame_corr2$WC[Data_frame_corr2$type=="subsoil"], log10(Data_frame_corr2$sum_Methanogens[Data_frame_corr2$type=="subsoil"]+0.001),method = "spearman")
cor.test(Data_frame_corr2$WC[Data_frame_corr2$type=="cryoOM"], log10(Data_frame_corr2$sum_Methanogens[Data_frame_corr2$type=="cryoOM"]+0.001),method = "spearman")
cor.test(Data_frame_corr2$WC[Data_frame_corr2$type=="permafrost"], log10(Data_frame_corr2$sum_Methanogens[Data_frame_corr2$type=="permafrost"]+0.001),method = "spearman")

# Data_frame_corr2 <- Data_frame_corr2[which(Data_frame_corr2$site != "Tazovskiy"),]
# 
# p1 <- ggplot(Data_frame_corr2,aes(x=log10(sum_Methanogens*100), y=WC))+
#   geom_point(aes(color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))),size=5,stroke=1,shape=19,alpha=0.9)+
#   geom_smooth(aes(group=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy"))),method = "lm",level=0.95,size=1.5,se=T,color="grey20")+
#   scale_color_brewer(palette="Dark2")+
#   guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)))+
#   theme(axis.title=element_text(size=15,face = "bold"),axis.text= element_text(size=12,face = "bold"),
#         panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
#         axis.line.y.left   = element_line(color = 'black',size=1),strip.text = element_text(size = 12,face = "bold"),
#         panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
#   labs(x="log10(Methanogen RA)",y="Water content (%)")+facet_grid(.~factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")))
# 
# cor.test(log10(Data_frame_corr2$sum_Methanogens[which(Data_frame_corr2$site=="Herschel")]), Data_frame_corr2$WC[which(Data_frame_corr2$site=="Herschel")],method = "spearman")
# 
# cor.test(Data_frame_corr2$sum_Methanogens[which(Data_frame_corr2$site=="Beaufort coast")], Data_frame_corr2$WC[which(Data_frame_corr2$site=="Beaufort coast")],method = "spearman")
# cor.test(log10(Data_frame_corr2$sum_Methanogens[which(Data_frame_corr2$type=="permafrost")]), Data_frame_corr2$WC[which(Data_frame_corr2$type=="permafrost")],method = "pearson")


Data_frame_corr3 <- Data_frame_corr[!is.na(Data_frame_corr$WC),]
which(Data_frame_corr3$sum_Methanotrophs==0)
Data_frame_corr3 <- Data_frame_corr3[which(Data_frame_corr3$WC>0),]
#Data_frame_corr3 <- Data_frame_corr3[which(Data_frame_corr3$sum_Methanotrophs>0),]

shapiro.test(log10(Data_frame_corr3$sum_Methanotrophs+0.001))

# pred <- predict(lm(WC ~ log10(sum_Methanotrophs+0.001), Data_frame_corr3),
#                 se.fit = TRUE, interval = "confidence")
# limits2 <- as.data.frame(pred$fit)
# 
# p2 <- ggplot(Data_frame_corr3,aes(x=log10(sum_Methanotrophs+0.001), y=WC))+
#   geom_point(aes(color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))),size=2,stroke=1.5,shape=19)+
#   geom_smooth(method = "lm",color="black",level=0.95,size=2,se=F)+
#   geom_line(aes(x = log10(sum_Methanotrophs+0.001), y = limits2$lwr), linetype = 2,size=1.5) +
#   geom_line(aes(x = log10(sum_Methanotrophs+0.001), y = limits2$upr), linetype = 2,size=1.5) +
#   scale_color_brewer(palette="Dark2")+
#   guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)))+
#   theme(axis.title=element_text(size=20,face = "bold"),axis.text= element_text(size=12,face = "bold"),
#         panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
#         axis.line.y.left   = element_line(color = 'black',size=1),
#         panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
#   labs(x="log10(MT RA)",y="Water content (%)")

p2 <- ggplot(Data_frame_corr3,aes(x=log10(sum_Methanotrophs+0.001), y=WC,
                                  color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost")),
                                  fill=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))))+
  geom_point(size=2,stroke=1.5,shape=19)+
  geom_smooth(method = "lm",size=2,se=T,level=0.95,alpha=0.2)+
  scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+
  guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)),fill="none")+
  theme(axis.title=element_text(size=20,face = "bold"),axis.text= element_text(size=12,face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
  labs(x="log10(MT RA)",y="Water content (%)")

p2


cor.test(Data_frame_corr3$WC[Data_frame_corr3$type=="organic layer"], log10(Data_frame_corr3$sum_Methanotrophs[Data_frame_corr3$type=="organic layer"]+0.001),method = "pearson")
cor.test(Data_frame_corr3$WC[Data_frame_corr3$type=="topsoil"], log10(Data_frame_corr3$sum_Methanotrophs[Data_frame_corr3$type=="topsoil"]+0.001),method = "pearson")
cor.test(Data_frame_corr3$WC[Data_frame_corr3$type=="subsoil"], log10(Data_frame_corr3$sum_Methanotrophs[Data_frame_corr3$type=="subsoil"]+0.001),method = "pearson")
cor.test(Data_frame_corr3$WC[Data_frame_corr3$type=="cryoOM"], log10(Data_frame_corr3$sum_Methanotrophs[Data_frame_corr3$type=="cryoOM"]+0.001),method = "pearson")
cor.test(Data_frame_corr3$WC[Data_frame_corr3$type=="permafrost"], log10(Data_frame_corr3$sum_Methanotrophs[Data_frame_corr3$type=="permafrost"]+0.001),method = "pearson")

cor.test(Data_frame_corr3$WC[Data_frame_corr3$type=="organic layer"], log10(Data_frame_corr3$sum_Methanotrophs[Data_frame_corr3$type=="organic layer"]+0.001),method = "spearman")
cor.test(Data_frame_corr3$WC[Data_frame_corr3$type=="topsoil"], log10(Data_frame_corr3$sum_Methanotrophs[Data_frame_corr3$type=="topsoil"]+0.001),method = "spearman")
cor.test(Data_frame_corr3$WC[Data_frame_corr3$type=="subsoil"], log10(Data_frame_corr3$sum_Methanotrophs[Data_frame_corr3$type=="subsoil"]+0.001),method = "spearman")
cor.test(Data_frame_corr3$WC[Data_frame_corr3$type=="cryoOM"], log10(Data_frame_corr3$sum_Methanotrophs[Data_frame_corr3$type=="cryoOM"]+0.001),method = "spearman")
cor.test(Data_frame_corr3$WC[Data_frame_corr3$type=="permafrost"], log10(Data_frame_corr3$sum_Methanotrophs[Data_frame_corr3$type=="permafrost"]+0.001),method = "spearman")
# Data_frame_corr3 <- Data_frame_corr3[which(Data_frame_corr3$site != "Tazovskiy"),]
# p2 <- ggplot(Data_frame_corr3,aes(x=log10(sum_Methanotrophs*100), y=WC))+
#   geom_point(aes(color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))),size=5,stroke=1,shape=19,alpha=0.9)+
#   geom_smooth(aes(group=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy"))),method = "lm",level=0.95,size=1.5,se=T,color="grey20")+
#   scale_color_brewer(palette="Dark2")+
#   guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)))+
#   theme(axis.title=element_text(size=15,face = "bold"),axis.text= element_text(size=12,face = "bold"),
#         panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
#         axis.line.y.left   = element_line(color = 'black',size=1),strip.text = element_text(size = 12,face = "bold"),legend.position = "none",
#         panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
#   labs(x="log10(Methanotroph RA)",y="Water content (%)")+facet_grid(.~factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")))
# 
# cor.test(log10(Data_frame_corr3$sum_Methanotrophs[which(Data_frame_corr3$site=="Cherskiy")]), Data_frame_corr3$WC[which(Data_frame_corr3$site=="Cherskiy")],method = "pearson")

# intersect(rownames(Data_frame_corr2),rownames(Data_frame_corr3))
# 
# Data_frame_corr6 <- Data_frame_corr[intersect(rownames(Data_frame_corr2),rownames(Data_frame_corr3)),]
# Data_frame_corr6$MG_MT <- log10(Data_frame_corr6$sum_Methanogens)-log10(Data_frame_corr6$sum_Methanotrophs)
# Data_frame_corr6$MG2MT <- log10(Data_frame_corr6$sum_Methanogens/Data_frame_corr6$sum_Methanotrophs)
# which(Data_frame_corr6$MG_MT==0)
# 
# pred <- predict(lm(WC ~ MG_MT, Data_frame_corr6),
#                 se.fit = TRUE, interval = "confidence")
# limits5 <- as.data.frame(pred$fit)
# 
# p3 <- ggplot(Data_frame_corr6,aes(x=MG_MT, y=WC))+
#   geom_point(aes(color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))),size=2,stroke=1.5,shape=19)+
#   geom_smooth(method = "lm",color="black",level=0.95,size=2,se=F)+
#   geom_line(aes(x = MG_MT, y = limits5$lwr), linetype = 2,size=1.5) +
#   geom_line(aes(x = MG_MT, y = limits5$upr), linetype = 2,size=1.5) +
#   scale_color_brewer(palette="Dark2")+
#   guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)))+
#   theme(axis.title=element_text(size=20,face = "bold"),axis.text= element_text(size=12,face = "bold"),
#         panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
#         axis.line.y.left   = element_line(color = 'black',size=1),
#         panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
#   labs(x="log10(MG RA)-log10(MT RA)",y="Water content (%)")
# p3
# cor.test(Data_frame_corr6$WC, Data_frame_corr6$MG2MT,method = "pearson")

# ggplot(Data_frame_corr6,aes(x=MG_MT, y=WC))+
#   geom_point(aes(color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))),size=5,stroke=1,shape=19,alpha=0.9)+
#   geom_smooth(aes(group=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy"))),method = "lm",level=0.95,size=1.5,se=T,color="grey20")+
#   scale_color_brewer(palette="Dark2")+
#   guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)))+
#   theme(axis.title=element_text(size=15,face = "bold"),axis.text= element_text(size=12,face = "bold"),
#         panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
#         axis.line.y.left   = element_line(color = 'black',size=1),strip.text = element_text(size = 12,face = "bold"),
#         panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
#   labs(x="log10(Methanogens/Methanotrophs)",y="Water content (%)")+facet_grid(.~factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")))
# 
# cor.test(Data_frame_corr6$WC[which(Data_frame_corr6$site=="Herschel")], Data_frame_corr6$MG_MT[which(Data_frame_corr6$site=="Herschel")],method = "pearson")
# 
# cor.test(Data_frame_corr6$WC[which(Data_frame_corr6$type=="permafrost")], Data_frame_corr6$MG_MT[which(Data_frame_corr6$type=="permafrost")],method = "pearson")

# between pH and MG/MT
Data_frame_corr4 <- Data_frame_corr[!is.na(Data_frame_corr$pH),]
which(Data_frame_corr4$pH==0)
which(Data_frame_corr4$sum_Methanogens==0)
Data_frame_corr4 <- Data_frame_corr4[which(Data_frame_corr4$pH>0),]
#Data_frame_corr4 <- Data_frame_corr4[which(Data_frame_corr4$sum_Methanogens>0),]

shapiro.test(log10(Data_frame_corr4$sum_Methanogens+0.001))
shapiro.test(log10(Data_frame_corr4$pH))


p4 <- ggplot(Data_frame_corr4,aes(x=log10(sum_Methanogens+0.001), y=pH,
                                   color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost")),
                                  fill=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))))+
  geom_point(size=2,stroke=1.5,shape=19)+
  geom_smooth(method = "lm",size=2,se=T,level=0.95,alpha=0.2)+
  scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+
  guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)),fill="none")+
  theme(axis.title=element_text(size=20,face = "bold"),axis.text= element_text(size=12,face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
  labs(x="log10(MG RA)",y="pH")
p4


cor.test(Data_frame_corr4$pH[Data_frame_corr4$type=="organic layer"], log10(Data_frame_corr4$sum_Methanogens[Data_frame_corr4$type=="organic layer"]+0.001),method = "pearson")
cor.test(Data_frame_corr4$pH[Data_frame_corr4$type=="topsoil"], log10(Data_frame_corr4$sum_Methanogens[Data_frame_corr4$type=="topsoil"]+0.001),method = "pearson")
cor.test(Data_frame_corr4$pH[Data_frame_corr4$type=="subsoil"], log10(Data_frame_corr4$sum_Methanogens[Data_frame_corr4$type=="subsoil"]+0.001),method = "pearson")
cor.test(Data_frame_corr4$pH[Data_frame_corr4$type=="cryoOM"], log10(Data_frame_corr4$sum_Methanogens[Data_frame_corr4$type=="cryoOM"]+0.001),method = "pearson")
cor.test(Data_frame_corr4$pH[Data_frame_corr4$type=="permafrost"], log10(Data_frame_corr4$sum_Methanogens[Data_frame_corr4$type=="permafrost"]+0.001),method = "pearson")

cor.test(Data_frame_corr4$pH[Data_frame_corr4$type=="organic layer"], log10(Data_frame_corr4$sum_Methanogens[Data_frame_corr4$type=="organic layer"]+0.001),method = "spearman")
cor.test(Data_frame_corr4$pH[Data_frame_corr4$type=="topsoil"], log10(Data_frame_corr4$sum_Methanogens[Data_frame_corr4$type=="topsoil"]+0.001),method = "spearman")
cor.test(Data_frame_corr4$pH[Data_frame_corr4$type=="subsoil"], log10(Data_frame_corr4$sum_Methanogens[Data_frame_corr4$type=="subsoil"]+0.001),method = "spearman")
cor.test(Data_frame_corr4$pH[Data_frame_corr4$type=="cryoOM"], log10(Data_frame_corr4$sum_Methanogens[Data_frame_corr4$type=="cryoOM"]+0.001),method = "spearman")
cor.test(Data_frame_corr4$pH[Data_frame_corr4$type=="permafrost"], log10(Data_frame_corr4$sum_Methanogens[Data_frame_corr4$type=="permafrost"]+0.001),method = "spearman")

# Data_frame_corr4 <- Data_frame_corr4[which(Data_frame_corr4$site != "Tazovskiy"),]
# p3 <- ggplot(Data_frame_corr4,aes(x=log10(sum_Methanogens*100), y=pH))+
#   geom_point(aes(color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))),size=5,stroke=1,shape=19,alpha=0.9)+
#   geom_smooth(aes(group=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy"))),method = "lm",level=0.95,size=1.5,se=T,color="grey20")+
#   scale_color_brewer(palette="Dark2")+
#   guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)))+
#   theme(axis.title=element_text(size=15,face = "bold"),axis.text= element_text(size=12,face = "bold"),
#         panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
#         axis.line.y.left   = element_line(color = 'black',size=1),strip.text = element_text(size = 12,face = "bold"),legend.position = "none",
#         panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
#   labs(x="log10(Methanogen RA)",y="pH")+facet_grid(.~factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")))
# 
# cor.test(log10(Data_frame_corr4$sum_Methanogens[which(Data_frame_corr4$site=="Cherskiy")]), Data_frame_corr4$pH[which(Data_frame_corr4$site=="Cherskiy")],method = "pearson")



Data_frame_corr5 <- Data_frame_corr[!is.na(Data_frame_corr$pH),]
Data_frame_corr5 <- Data_frame_corr5[which(Data_frame_corr5$pH>0),]
#Data_frame_corr5 <- Data_frame_corr5[which(Data_frame_corr5$sum_Methanotrophs>0),]

# pred <- predict(lm(pH ~ log10(sum_Methanotrophs+0.001), Data_frame_corr5),
#                 se.fit = TRUE, interval = "confidence")
# limits4 <- as.data.frame(pred$fit)

p5 <- ggplot(Data_frame_corr5,aes(x=log10(sum_Methanotrophs+0.001), y=pH,
                                  color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost")),
                                  fill=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))))+
  geom_point(size=2,stroke=1.5,shape=19)+
  geom_smooth(method = "lm",size=2,se=T,level=0.95,alpha=0.2)+
  scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+
  guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)),fill="none")+
  theme(axis.title=element_text(size=20,face = "bold"),axis.text= element_text(size=12,face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
  labs(x="log10(MT RA)",y="pH")
p5

cor.test(Data_frame_corr5$pH[Data_frame_corr5$type=="organic layer"], log10(Data_frame_corr5$sum_Methanotrophs[Data_frame_corr5$type=="organic layer"]+0.001),method = "pearson")
cor.test(Data_frame_corr5$pH[Data_frame_corr5$type=="topsoil"], log10(Data_frame_corr5$sum_Methanotrophs[Data_frame_corr5$type=="topsoil"]+0.001),method = "pearson")
cor.test(Data_frame_corr5$pH[Data_frame_corr5$type=="subsoil"], log10(Data_frame_corr5$sum_Methanotrophs[Data_frame_corr5$type=="subsoil"]+0.001),method = "pearson")
cor.test(Data_frame_corr5$pH[Data_frame_corr5$type=="cryoOM"], log10(Data_frame_corr5$sum_Methanotrophs[Data_frame_corr5$type=="cryoOM"]+0.001),method = "pearson")
cor.test(Data_frame_corr5$pH[Data_frame_corr5$type=="permafrost"], log10(Data_frame_corr5$sum_Methanotrophs[Data_frame_corr5$type=="permafrost"]+0.001),method = "pearson")

cor.test(Data_frame_corr5$pH[Data_frame_corr5$type=="organic layer"], log10(Data_frame_corr5$sum_Methanotrophs[Data_frame_corr5$type=="organic layer"]+0.001),method = "spearman")
cor.test(Data_frame_corr5$pH[Data_frame_corr5$type=="topsoil"], log10(Data_frame_corr5$sum_Methanotrophs[Data_frame_corr5$type=="topsoil"]+0.001),method = "spearman")
cor.test(Data_frame_corr5$pH[Data_frame_corr5$type=="subsoil"], log10(Data_frame_corr5$sum_Methanotrophs[Data_frame_corr5$type=="subsoil"]+0.001),method = "spearman")
cor.test(Data_frame_corr5$pH[Data_frame_corr5$type=="cryoOM"], log10(Data_frame_corr5$sum_Methanotrophs[Data_frame_corr5$type=="cryoOM"]+0.001),method = "spearman")
cor.test(Data_frame_corr5$pH[Data_frame_corr5$type=="permafrost"], log10(Data_frame_corr5$sum_Methanotrophs[Data_frame_corr5$type=="permafrost"]+0.001),method = "spearman")
# Data_frame_corr5 <- Data_frame_corr5[which(Data_frame_corr5$site != "Tazovskiy"),]
# p4 <- ggplot(Data_frame_corr5,aes(x=log10(sum_Methanotrophs*100), y=pH))+
#   geom_point(aes(color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))),size=5,stroke=1,shape=19,alpha=0.9)+
#   geom_smooth(aes(group=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy"))),method = "lm",level=0.95,size=1.5,se=T,color="grey20")+
#   scale_color_brewer(palette="Dark2")+
#   guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)))+
#   theme(axis.title=element_text(size=15,face = "bold"),axis.text= element_text(size=12,face = "bold"),
#         panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
#         axis.line.y.left   = element_line(color = 'black',size=1),strip.text = element_text(size = 12,face = "bold"),legend.position = "none",
#         panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
#   labs(x="log10(Methanotroph RA)",y="pH")+facet_grid(.~factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")))
# 
# cor.test(log10(Data_frame_corr5$sum_Methanotrophs[which(Data_frame_corr5$site=="Herschel")]), Data_frame_corr5$pH[which(Data_frame_corr5$site=="Herschel")],method = "pearson")


# 
# intersect(rownames(Data_frame_corr4),rownames(Data_frame_corr5))
# 
# Data_frame_corr7 <- Data_frame_corr[intersect(rownames(Data_frame_corr4),rownames(Data_frame_corr5)),]
# Data_frame_corr7$MG_MT <- log10(Data_frame_corr7$sum_Methanogens)-log10(Data_frame_corr7$sum_Methanotrophs)
# Data_frame_corr7$MG2MT <- log10(Data_frame_corr7$sum_Methanogens/Data_frame_corr7$sum_Methanotrophs)
# which(Data_frame_corr7$MG_MT==0)
# 
# pred <- predict(lm(pH ~ MG_MT, Data_frame_corr7),
#                 se.fit = TRUE, interval = "confidence")
# limits6 <- as.data.frame(pred$fit)
# 
# p6 <- ggplot(Data_frame_corr7,aes(x=MG_MT, y=pH))+geom_point(aes(color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))),size=2,stroke=1.5,shape=19)+
#   geom_smooth(method = "lm",color="black",level=0.95,size=2,se=F)+
#   geom_line(aes(x = MG_MT, y = limits6$lwr), linetype = 2,size=1.5) +
#   geom_line(aes(x = MG_MT, y = limits6$upr), linetype = 2,size=1.5) +
#   scale_color_brewer(palette="Dark2")+
#   guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)))+
#   theme(axis.title=element_text(size=20,face = "bold"),axis.text= element_text(size=12,face = "bold"),
#         panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
#         axis.line.y.left   = element_line(color = 'black',size=1),
#         panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
#   labs(x="log10(MG RA)-log10(MT RA)",y="pH")
# p6
# cor.test(Data_frame_corr7$pH, Data_frame_corr7$MG_MT,method = "pearson")

# ggplot(Data_frame_corr7,aes(x=MG_MT, y=pH))+
#   geom_point(aes(color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))),size=5,stroke=1,shape=19,alpha=0.9)+
#   geom_smooth(aes(group=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy"))),method = "lm",level=0.95,size=1.5,se=T,color="grey20")+
#   scale_color_brewer(palette="Dark2")+
#   guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)))+
#   theme(axis.title=element_text(size=15,face = "bold"),axis.text= element_text(size=12,face = "bold"),
#         panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
#         axis.line.y.left   = element_line(color = 'black',size=1),strip.text = element_text(size = 12,face = "bold"),
#         panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
#   labs(x="log10(Methanogens/Methanotrophs)",y="pH")+facet_grid(.~factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")))
# 
# cor.test(Data_frame_corr7$pH[which(Data_frame_corr7$site=="Cherskiy")], Data_frame_corr7$MG_MT[which(Data_frame_corr7$site=="Cherskiy")],method = "pearson")

# between MG and MT
#Data_frame_corr8 <- Data_frame_corr[which(Data_frame_corr$sum_Methanogens>0 & Data_frame_corr$sum_Methanotrophs>0),]
Data_frame_corr8 <- Data_frame_corr

# pred <- predict(lm(log10(sum_Methanotrophs+0.001) ~ log10(sum_Methanogens+0.001), Data_frame_corr8),
#                 se.fit = TRUE, interval = "confidence")
# limits7 <- as.data.frame(pred$fit)

p7 <- ggplot(Data_frame_corr8,aes(x=log10(sum_Methanogens+0.001), y=log10(sum_Methanotrophs+0.001),
                                  color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost")),
                                  fill=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))))+
  geom_point(size=2,stroke=1.5,shape=19)+
  geom_smooth(method = "lm",size=2,se=T,level=0.95,alpha=0.2)+
  scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+
  guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)),fill="none")+
  theme(axis.title=element_text(size=20,face = "bold"),axis.text= element_text(size=12,face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
  labs(x="log10(MG RA)",y="log10(MT RA)")
p7
# cor.test(log10(Data_frame_corr8$sum_Methanogens+0.001), log10(Data_frame_corr8$sum_Methanotrophs+0.001),method = "pearson")
# cor.test(log10(Data_frame_corr8$sum_Methanogens+0.001), log10(Data_frame_corr8$sum_Methanotrophs+0.001),method = "spearman")


cor.test(log10(Data_frame_corr8$sum_Methanogens[Data_frame_corr8$type=="organic layer"]+0.001), log10(Data_frame_corr8$sum_Methanotrophs[Data_frame_corr8$type=="organic layer"]+0.001),method = "pearson")
cor.test(log10(Data_frame_corr8$sum_Methanogens[Data_frame_corr8$type=="topsoil"]+0.001), log10(Data_frame_corr8$sum_Methanotrophs[Data_frame_corr8$type=="topsoil"]+0.001),method = "pearson")
cor.test(log10(Data_frame_corr8$sum_Methanogens[Data_frame_corr8$type=="subsoil"]+0.001), log10(Data_frame_corr8$sum_Methanotrophs[Data_frame_corr8$type=="subsoil"]+0.001),method = "pearson")
cor.test(log10(Data_frame_corr8$sum_Methanogens[Data_frame_corr8$type=="cryoOM"]+0.001), log10(Data_frame_corr8$sum_Methanotrophs[Data_frame_corr8$type=="cryoOM"]+0.001),method = "pearson")
cor.test(log10(Data_frame_corr8$sum_Methanogens[Data_frame_corr8$type=="permafrost"]+0.001), log10(Data_frame_corr8$sum_Methanotrophs[Data_frame_corr8$type=="permafrost"]+0.001),method = "pearson")

cor.test(log10(Data_frame_corr8$sum_Methanogens[Data_frame_corr8$type=="organic layer"]+0.001), log10(Data_frame_corr8$sum_Methanotrophs[Data_frame_corr8$type=="organic layer"]+0.001),method = "spearman")
cor.test(log10(Data_frame_corr8$sum_Methanogens[Data_frame_corr8$type=="topsoil"]+0.001), log10(Data_frame_corr8$sum_Methanotrophs[Data_frame_corr8$type=="topsoil"]+0.001),method = "spearman")
cor.test(log10(Data_frame_corr8$sum_Methanogens[Data_frame_corr8$type=="subsoil"]+0.001), log10(Data_frame_corr8$sum_Methanotrophs[Data_frame_corr8$type=="subsoil"]+0.001),method = "spearman")
cor.test(log10(Data_frame_corr8$sum_Methanogens[Data_frame_corr8$type=="cryoOM"]+0.001), log10(Data_frame_corr8$sum_Methanotrophs[Data_frame_corr8$type=="cryoOM"]+0.001),method = "spearman")
cor.test(log10(Data_frame_corr8$sum_Methanogens[Data_frame_corr8$type=="permafrost"]+0.001), log10(Data_frame_corr8$sum_Methanotrophs[Data_frame_corr8$type=="permafrost"]+0.001),method = "spearman")


# Data_frame_corr8 <- Data_frame_corr8[which(Data_frame_corr8$site != "Tazovskiy"),]
# 
# 
# p5 <- ggplot(Data_frame_corr8,aes(x=log10(sum_Methanogens*100), y=log10(sum_Methanotrophs*100)))+
#   geom_point(aes(color=factor(type,levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))),size=5,stroke=1,shape=19,alpha=0.9)+
#   geom_smooth(aes(group=factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy"))),method = "lm",level=0.95,size=1.5,se=T,color="grey20")+
#   scale_color_brewer(palette="Dark2")+
#   guides(color=guide_legend(title="Horizon",title.theme = element_text(size = 15,face = "bold"),label.theme = element_text(size=15),override.ae=list(size=8)))+
#   theme(axis.title=element_text(size=15,face = "bold"),axis.text= element_text(size=12,face = "bold"),
#         panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
#         axis.line.y.left   = element_line(color = 'black',size=1),strip.text = element_text(size = 12,face = "bold"),legend.position = "none",
#         panel.background = element_rect(fill = NA),panel.grid.major = element_blank()) +
#   labs(x="log10(Methanogen RA)",y="log10(Methanotroph RA)")+facet_grid(.~factor(site,levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy")))
# 
# cor.test(log10(Data_frame_corr8$sum_Methanogens[which(Data_frame_corr8$site=="Cherskiy")]), log10(Data_frame_corr8$sum_Methanotrophs[which(Data_frame_corr8$site=="Cherskiy")]),method = "spearman")

library(cowplot)
# plot_grid(p1,p2,p3,p4,p5,ncol = 1,align = "hv",axis = "tblr")
#plot_grid(p1,p2,p3,p4,p5,p6,p7, ncol = 3,align = "hv",axis = "tblr")

plot_grid(p1,p2,p4,p5,p7, ncol = 2,align = "hv",axis = "tblr")
