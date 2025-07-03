#creates plots showing the community distribution and the abundances of Methanogens and Methanotrophs, averaged for soil horizons within sites and for sites within soil horizons.

library(vegan)
library(tidyverse)
library(writexl)
library(ggplot2)

load("data_frames_sitecat")
load("data_frames_cat")
load("data_frames_spec")
load("Site_table_new")
Site_table$rownames <- rownames(Site_table)
Site_table$category <- as.character(Site_table$category)

addline_format <- function(x,...){
  gsub(',','\n',x)
}

#site Methanogens
for (i in names(data_frames_spec)[c(1,3,5,7,9,11,13,15)]) {
  table <- data_frames_spec[[i]]
  table <- as.data.frame(t(table))
  
  for (n in c(1:nrow(table))) {
    True_False <- rownames(Site_table) == rownames(table)[n]
    if(any(True_False) == TRUE){
      table$type[n] <- Site_table$category[True_False == TRUE]
    }
  }
  
  for (n in c(1:(ncol(table)-1))) {
    table[,n] <- as.numeric(table[,n])
  }
  table$colsums <- rowSums(table[,1:(ncol(table)-1)])
  table_ii <- table
  
  table <- data.frame()
  for(n in unique(table_ii$type)){
    subtable <- table_ii[table_ii$type == n,]
    subtable[,ncol(subtable)] <- as.numeric(subtable[,ncol(subtable)])
    mean_colsums <- mean(subtable$colsums)
    SE <- sd(subtable$colsums)/sqrt(length(subtable$colsums))
    number <- nrow(subtable)
    row <- data.frame(mean_colsums, SE, number, row.names = n)
    table <- rbind(table, row)
  }
  table <- as.data.frame(t(table))
  
  xValue <- rep(colnames(table[1:(nrow(table)-1),]), each = 1)
  xValue <- factor(xValue, levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))
  levels <- list()
  for (n in levels(xValue)) {
    levels[[n]] <- paste(n, table[nrow(table),n], sep = ", n = ")
  }
  xValue <- factor(paste(xValue, table[nrow(table),], sep = ", n = "), levels = as.character(levels))
  table <- as.data.frame(t(table))
  yValue <- 100*table$mean_colsums
  SEValue <- 100*table$SE
  data <- data.frame(xValue,yValue,SEValue)
  data <- arrange(data, xValue)
  ggplot(data) +
    geom_col(aes(x=fct_rev(xValue), y=yValue))+
    geom_errorbar(aes(x=fct_rev(xValue), ymin=yValue, ymax=yValue+SEValue))+
    coord_flip()+
    scale_x_discrete(labels=addline_format(rev(data$xValue)))+
    xlab("")+
    ylab("Relative Abundance (%)")+
    scale_y_continuous(limits = c(0,2.2),expand = c(0,0))+
    theme(axis.text.y = element_text(size = 15, color = "black",face = "bold"),
          axis.text.x = element_text(size = 13.5, color = "black"),legend.position = "none",
          axis.title.x = element_text(size = 17, color = "black",face = "bold"),
          panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
          axis.line.y.left   = element_line(color = 'black',size=1),
          panel.background = element_rect(fill = NA),panel.grid.major = element_blank())
  ggsave(filename = paste(i, " colsums.png"), device = png(), width = 20, height = 12, units = "cm")
}
graphics.off()

#site Methanotrophs
for (i in names(data_frames_spec)[c(2,4,6,8,10,12,14,16)]) {
  table <- data_frames_spec[[i]]
  table <- as.data.frame(t(table))
  
  for (n in c(1:nrow(table))) {
    True_False <- rownames(Site_table) == rownames(table)[n]
    if(any(True_False) == TRUE){
      table$type[n] <- Site_table$category[True_False == TRUE]
    }
  }
  
  for (n in c(1:(ncol(table)-1))) {
    table[,n] <- as.numeric(table[,n])
  }
  table$colsums <- rowSums(table[,1:(ncol(table)-1)])
  table_ii <- table
  
  table <- data.frame()
  for(n in unique(table_ii$type)){
    subtable <- table_ii[table_ii$type == n,]
    subtable[,ncol(subtable)] <- as.numeric(subtable[,ncol(subtable)])
    mean_colsums <- mean(subtable$colsums)
    SE <- sd(subtable$colsums)/sqrt(length(subtable$colsums))
    number <- nrow(subtable)
    row <- data.frame(mean_colsums, SE, number, row.names = n)
    table <- rbind(table, row)
  }
  table <- as.data.frame(t(table))
  
  xValue <- rep(colnames(table[1:(nrow(table)-1),]), each = 1)
  xValue <- factor(xValue, levels = c("organic layer","topsoil","subsoil","cryoOM","permafrost"))
  levels <- list()
  for (n in levels(xValue)) {
    levels[[n]] <- paste(n, table[nrow(table),n], sep = ", n = ")
  }
  xValue <- factor(paste(xValue, table[nrow(table),], sep = ", n = "), levels = as.character(levels))
  table <- as.data.frame(t(table))
  yValue <- 100*table$mean_colsums
  SEValue <- 100*table$SE
  data <- data.frame(xValue,yValue,SEValue)
  data <- arrange(data, xValue)
  ggplot(data) +
    geom_col(aes(x=fct_rev(xValue), y=yValue))+
    geom_errorbar(aes(x=fct_rev(xValue), ymin=yValue, ymax=yValue+SEValue))+
    coord_flip()+
    scale_x_discrete(labels=addline_format(rev(data$xValue)))+
    xlab("")+
    ylab("Relative Abundance (%)")+
    scale_y_continuous(limits = c(0,1.4),expand = c(0,0))+
    theme(axis.text.y = element_text(size = 15, color = "black",face = "bold"),
          axis.text.x = element_text(size = 13.5, color = "black"),legend.position = "none",
          axis.title.x = element_text(size = 17, color = "black",face = "bold"),
          panel.border = element_blank(),axis.line.x.bottom = element_line(color = 'black',size=1),
          axis.line.y.left   = element_line(color = 'black',size=1),
          panel.background = element_rect(fill = NA),panel.grid.major = element_blank())
  ggsave(filename = paste(i, " colsums.png"), device = png(), width = 20, height = 12, units = "cm")
}
graphics.off()

###### I stopped here


