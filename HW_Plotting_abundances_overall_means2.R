#creates plots showing the mean relative abundnaces of Methanogens and Methanotrophs within sites and soil horizons with error bars.

library(vegan)
library(tidyverse)
library(writexl)
library(ggplot2)

load("data_frames_sitecat")
load("data_frames_cat")
load("data_frames_spec")

### statistics
list_colsums_stat <- list()
for (i in names(data_frames_cat)) {
  table <- data_frames_cat[[i]]
  list_colsums_stat[[i]] <- colSums(table)
}
list_colsums_stat[[1]]

methangen_stat <- data.frame()
for (i in seq(1,15,by=2)) {
  table <- data.frame(list_colsums_stat[[i]])
  methangen_stat <- rbind(methangen_stat,table)
}

load("Site_table_new")
length(intersect(rownames(Site_table),rownames(methangen_stat)))

Site_table <- Site_table[rownames(methangen_stat),]
methangen_stat <- cbind(methangen_stat,Site_table)
colnames(methangen_stat)[1] <- "MG"

methangen_stat2 <- methangen_stat[which(methangen_stat$MG!=0),]

plyr::count(methangen_stat2$site)


library(PMCMR)
posthoc.kruskal.dunn.test(MG ~ factor(site),data = methangen_stat, p.adjust = "fdr")
posthoc.kruskal.dunn.test(MG ~ factor(category),data = methangen_stat, p.adjust = "fdr")
pairwise.t.test(log10(methangen_stat$MG+0.001),methangen_stat$category,p.adjust.method = "fdr")
pairwise.t.test(log10(methangen_stat$MG+0.001),methangen_stat$site,p.adjust.method = "fdr")


methanotroph_stat <- data.frame()
for (i in seq(2,16,by=2)) {
  table <- data.frame(list_colsums_stat[[i]])
  methanotroph_stat <- rbind(methanotroph_stat,table)
}

Site_table <- Site_table[rownames(methanotroph_stat),]
methanotroph_stat <- cbind(methanotroph_stat,Site_table)
colnames(methanotroph_stat)[1] <- "MT"

posthoc.kruskal.dunn.test(MT ~ factor(site),data = methanotroph_stat, p.adjust = "fdr")
posthoc.kruskal.dunn.test(MT ~ factor(category),data = methanotroph_stat, p.adjust = "fdr")
pairwise.t.test(log10(methanotroph_stat$MT+0.01),methanotroph_stat$category,p.adjust.method = "fdr")
pairwise.t.test(log10(methanotroph_stat$MT+0.01),methanotroph_stat$site,p.adjust.method = "fdr")

#############
list_colsums <- list()
list_SD <- list()
list_n <- list()
for (i in names(data_frames_cat)) {
  table <- data_frames_cat[[i]]
  list_colsums[[i]] <- mean(colSums(table))
  list_SD[[i]] <- sd(colSums(table))
  list_n[[i]] <- ncol(table)
}

data_frames_true_false <- grepl("Methanogens", names(list_colsums))
names_to_use <- names(list_colsums)[data_frames_true_false == TRUE]
names_to_use_ii <- names(list_colsums)[data_frames_true_false == FALSE]

Methanogens <- list()
Methanogens_SD <- list()
Methanogens_n <- list()
for (i in names_to_use) {
  result <- list_colsums[[i]]
  Methanogens[[i]] <- result
  result <- list_SD[[i]]
  Methanogens_SD[[i]] <- result
  result <- list_n[[i]]
  Methanogens_n[[i]] <- sqrt(result)
}

Methanotrophs <- list()
Methanotrophs_SD <- list()
Methanotrophs_n <- list()
for (i in names_to_use_ii) {
  result <- list_colsums[[i]]
  Methanotrophs[[i]] <- result
  result <- list_SD[[i]]
  Methanotrophs_SD[[i]] <- result
  result <- list_n[[i]]
  Methanotrophs_n[[i]] <- sqrt(result)
}

addline_format <- function(x,...){
  gsub(',','\n',x)
}

#Methanogens sites
xValue <- c("AriMas","Beaufort coast","Cherskiy","Disko","Herschel","Logata","Tazovskiy","Zackenberg")
yValue <- 100*as.numeric(Methanogens[1:8])
SEValue <- 100*as.numeric(as.numeric(Methanogens_SD[1:8])/as.numeric(Methanogens_n[1:8]))
xValue <- factor(xValue, levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy"))
list_n <- list()
for (i in names(data_frames_cat)[c(1,3,5,7,9,11,13,15)]) {
  table <- data_frames_cat[[i]]
  list_n[[i]] <- ncol(table)
}
levels <- list()
for (n in levels(xValue)) {
  table <- data_frames_cat[[names(data_frames_cat)[names(data_frames_cat) == paste(n, " Methanogens (cat.)", sep = "")]]]
  levels[[n]] <- paste(n, ncol(table), sep = ", n = ")
}
xValue <- factor(paste(xValue, list_n, sep = ", n = "), levels = as.character(levels))
data <- data.frame(xValue,yValue,SEValue)
data <- arrange(data, xValue)
p1 <- ggplot(data) +
  geom_col(aes(x=xValue, y=yValue),color="black",fill="grey50")+
  geom_errorbar(aes(x=xValue, ymin=yValue, ymax=yValue+SEValue))+
  scale_x_discrete(labels=addline_format(data$xValue))+
  xlab("")+ ylab("Relative abundance (%)")+
  scale_y_continuous(limits = c(0,1.1),expand = c(0,0),breaks = c(0,0.3,0.6,0.9))+
  theme(axis.text.x = element_text(angle = 90, size = 15, color = "black",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.y = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())
p1

#Methanogens layers
xValue <- c("cryoOM","organic layer","permafrost","subsoil","topsoil")
yValue <- 100*as.numeric(Methanogens[9:13])
SEValue <- 100*as.numeric(as.numeric(Methanogens_SD[9:13])/as.numeric(Methanogens_n[9:13]))
xValue <- factor(xValue, levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))
list_n <- list()
for (i in names(data_frames_cat)[c(17,19,21,23,25)]) {
  table <- data_frames_cat[[i]]
  list_n[[i]] <- ncol(table)
}
levels <- list()
for (n in levels(xValue)) {
  table <- data_frames_cat[[names(data_frames_cat)[names(data_frames_cat) == paste(n, " Methanogens (cat.)", sep = "")]]]
  levels[[n]] <- paste(n, ncol(table), sep = ", n = ")
}
xValue <- factor(paste(xValue, list_n, sep = ", n = "), levels = as.character(levels))
data <- data.frame(xValue,yValue,SEValue)
data <- arrange(data, xValue)
p2 <- ggplot(data) +
  geom_col(aes(x=fct_rev(xValue), y=yValue),color="black",fill="grey50")+
  geom_errorbar(aes(x=fct_rev(xValue), ymin=yValue, ymax=yValue+SEValue))+
  scale_x_discrete(labels=rev(addline_format(data$xValue)))+
  coord_flip()+
  xlab("")+ylab("Relative abundance (%)")+
  scale_y_continuous(limits = c(0,0.52),expand = c(0,0))+
  theme(axis.text.y = element_text(size = 15, color = "black",face = "bold"),
        axis.text.x = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.x = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())
p2


#Methanotrophs sites
xValue <- c("AriMas","Beaufort coast","Cherskiy","Disko","Herschel","Logata","Tazovskiy","Zackenberg")
yValue <- 100*as.numeric(Methanotrophs[1:8])
SEValue <- 100*as.numeric(as.numeric(Methanotrophs_SD[1:8])/as.numeric(Methanotrophs_n[1:8]))
xValue <- factor(xValue, levels = c("Herschel","Beaufort coast","Disko","Zackenberg","Tazovskiy","AriMas","Logata","Cherskiy"))
list_n <- list()
for (i in names(data_frames_cat)[c(2,4,6,8,10,12,14,16)]) {
  table <- data_frames_cat[[i]]
  list_n[[i]] <- ncol(table)
}
levels <- list()
for (n in levels(xValue)) {
  table <- data_frames_cat[[names(data_frames_cat)[names(data_frames_cat) == paste(n, " Methanotrophs (cat.)", sep = "")]]]
  levels[[n]] <- paste(n, ncol(table), sep = ", n = ")
}
xValue <- factor(paste(xValue, list_n, sep = ", n = "), levels = as.character(levels))
data <- data.frame(xValue,yValue,SEValue)
data <- arrange(data, xValue)
p3 <- ggplot(data, aes(x=xValue, y=yValue)) +
  geom_col(aes(x=xValue, y=yValue),color="black",fill="grey50")+
  geom_errorbar(aes(x=xValue, ymin=yValue, ymax=yValue+SEValue))+
  scale_x_discrete(labels=addline_format(data$xValue))+
  xlab("")+ ylab("Relative abundance (%)")+
  scale_y_continuous(limits = c(0,0.52),expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90, size = 15, color = "black",face = "bold",hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.y = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())
p3

#Methanotrophs layers
xValue <- c("cryoOM","organic layer","permafrost","subsoil","topsoil")
yValue <- 100*as.numeric(Methanotrophs[9:13])
SEValue <- 100*as.numeric(as.numeric(Methanotrophs_SD[9:13])/as.numeric(Methanotrophs_n[9:13]))
xValue <- factor(xValue, levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))
list_n <- list()
for (i in names(data_frames_cat)[c(18,20,22,24,26)]) {
  table <- data_frames_cat[[i]]
  list_n[[i]] <- ncol(table)
}
levels <- list()
for (n in levels(xValue)) {
  table <- data_frames_cat[[names(data_frames_cat)[names(data_frames_cat) == paste(n, " Methanotrophs (cat.)", sep = "")]]]
  levels[[n]] <- paste(n, ncol(table), sep = ", n = ")
}
xValue <- factor(paste(xValue, list_n, sep = ", n = "), levels = as.character(levels))
data <- data.frame(xValue,yValue,SEValue)
data <- arrange(data, xValue)
p4 <- ggplot(data) +
  geom_col(aes(x=fct_rev(xValue), y=yValue),color="black",fill="grey50")+
  geom_errorbar(aes(x=fct_rev(xValue), ymin=yValue, ymax=yValue+SEValue))+
  scale_x_discrete(labels=rev(addline_format(data$xValue)))+
  coord_flip()+
  xlab("")+ylab("Relative abundance (%)")+
  scale_y_continuous(limits = c(0,0.32),expand = c(0,0))+
  theme(axis.text.y = element_text(size = 15, color = "black",face = "bold"),
        axis.text.x = element_text(size = 13.5, color = "black"),legend.position = "none",
        axis.title.x = element_text(size = 17, color = "black",face = "bold"),
        panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
        axis.line.y.left   = element_line(color = 'black',size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_blank())
p4
library(cowplot)
plot_grid(p1,p3,p2,p4,ncol = 2,align = "v",axis = "tblr",rel_heights = c(2,1.3))









