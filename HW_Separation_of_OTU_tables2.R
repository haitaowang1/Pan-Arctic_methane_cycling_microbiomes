#splits up the original OTU tables, taxonomically categorised OTU tables and functionally categorised OTU tables by site and by soil horizon, creating subtables for each of these categories. Saves these as lists of data.frames. Saves both the tables containing the individual samples (no filename addition) as well as versions where abundnaces for sites within the soil horizon subtables and soil horizons within the site subtables are averaged (addition 'sitecat').


library(vegan)
library(tidyverse)
library(writexl)

load("../Erik_thesis/Site_table")  # this is the wrong table
load("../Erik_thesis/OTU_table_Methanogens")
load("../Erik_thesis/OTU_table_Methanotrophs")
load("../Erik_thesis/OTU_table_Methanogens_cat")
load("../Erik_thesis/OTU_table_Methanotrophs_cat")
load("../Erik_thesis/OTU_table_Methanogens_types")
load("../Erik_thesis/OTU_table_Methanotrophs_types")


Site_table <- read.csv("bacteria_mapping_final.csv",row.names = 1)

Site_table <- Site_table[colnames(OTU_table_Methanogens),]

data_frames <- list()
data_frames_sitecat <- list()
data_frames_Methanogens_types <- list()
data_frames_Methanogens_types_sitecat <- list()
data_frames_Methanotrophs_types <- list()
data_frames_Methanotrophs_types_sitecat <- list()
data_total_Methanogens <- data.frame()
data_total_Methanotrophs <- data.frame()
for(site in c("AriMas","Beaufort coast","Cherskiy","Disko","Herschel","Logata","Tazovskiy","Zackenberg")) {
  Site_table_red <- Site_table[which(Site_table$site == site),]
  list <- rownames(Site_table_red)
  table_Methanogens <- select(OTU_table_Methanogens, list)
  data_frames[[paste(site, " Methanogens", sep = "")]] <- table_Methanogens
  table_Methanogens_cat <- select(OTU_table_Methanogens_cat, list)
  data_frames[[paste(site, " Methanogens (cat.)", sep = "")]] <- table_Methanogens_cat
  table_Methanotrophs <- select(OTU_table_Methanotrophs, list)
  data_frames[[paste(site, " Methanotrophs", sep = "")]] <- table_Methanotrophs
  table_Methanotrophs_cat <- select(OTU_table_Methanotrophs_cat, list)
  data_frames[[paste(site, " Methanotrophs (cat.)", sep = "")]] <- table_Methanotrophs_cat
  
  types_table <- data.frame()
  list_number <- list()
  for (y in unique(Site_table_red$category)) {
    list_type <- rownames(Site_table_red[Site_table_red$category == y,])
    type <- as.data.frame(t(rowMeans(table_Methanogens_cat[list_type])), row.names = y)
    list_number[[y]] <- nrow(Site_table_red[Site_table_red$category == y,])
    types_table <- rbind(types_table, type)
  }
  types_table$n <- as.character(list_number)
  types_table <- as.data.frame(t(types_table))
  data_frames_sitecat[[paste(site, " Methanogens (sitecat.)", sep = "")]] <- types_table
  
  types_table <- data.frame()
  list_number <- list()
  for (y in unique(Site_table_red$category)) {
    list_type <- rownames(Site_table_red[Site_table_red$category == y,])
    type <- as.data.frame(t(rowMeans(table_Methanotrophs_cat[list_type])), row.names = y)
    list_number[[y]] <- nrow(Site_table_red[Site_table_red$category == y,])
    types_table <- rbind(types_table, type)
  }
  types_table$n <- as.character(list_number)
  types_table <- as.data.frame(t(types_table))
  data_frames_sitecat[[paste(site, " Methanotrophs (sitecat.)", sep = "")]] <- types_table
  
  table_Methanogens_types <- select(OTU_table_Methanogens_types, list)
  data_frames_Methanogens_types[[paste(site, " Methanogens types", sep = "")]] <- table_Methanogens_types
  
  table_Methanotrophs_types <- select(OTU_table_Methanotrophs_types, list)
  data_frames_Methanotrophs_types[[paste(site, " Methanotrophs types", sep = "")]] <- table_Methanotrophs_types
  
  types_table_Methanogens_types <- data.frame()
  list_number <- list()
  for (y in unique(Site_table_red$category)) {
    list_type <- rownames(Site_table_red[Site_table_red$category == y,])
    type <- as.data.frame(t(rowMeans(table_Methanogens_types[list_type])), row.names = y)
    list_number[[y]] <- nrow(Site_table_red[Site_table_red$category == y,])
    types_table_Methanogens_types <- rbind(types_table_Methanogens_types, type)
  }
  types_table_Methanogens_types$n <- as.character(list_number)
  types_table_Methanogens_types <- as.data.frame(t(types_table_Methanogens_types))
  data_frames_Methanogens_types_sitecat[[paste(site, " Methanogens types (sitecat.)", sep = "")]] <- types_table_Methanogens_types
  
  types_table_Methanotroph_types <- data.frame()
  list_number <- list()
  for (y in unique(Site_table_red$category)) {
    list_type <- rownames(Site_table_red[Site_table_red$category == y,])
    type <- as.data.frame(t(rowMeans(table_Methanotrophs_types[list_type])), row.names = y)
    list_number[[y]] <- nrow(Site_table_red[Site_table_red$category == y,])
    types_table_Methanotroph_types <- rbind(types_table_Methanotroph_types, type)
  }
  types_table_Methanotroph_types$n <- as.character(list_number)
  types_table_Methanotroph_types <- as.data.frame(t(types_table_Methanotroph_types))
  data_frames_Methanotrophs_types_sitecat[[paste(site, " Methanotrophs types (sitecat.)", sep = "")]] <- types_table_Methanotroph_types
  
  row_total_Methanogens <- as.data.frame(t(rowMeans(table_Methanogens_cat)), row.names = site)
  row_total_Methanogens$n <- as.character(ncol(table_Methanogens_cat))
  data_total_Methanogens <- rbind(data_total_Methanogens, row_total_Methanogens)
  row_total_Methanotrophs <- as.data.frame(t(rowMeans(table_Methanotrophs_cat)), row.names = site)
  row_total_Methanotrophs$n <- as.character(ncol(table_Methanotrophs_cat))
  data_total_Methanotrophs <- rbind(data_total_Methanotrophs, row_total_Methanotrophs)
}

for(type in c("cryoOM","organic layer","permafrost","subsoil","topsoil")) {
  Site_table_red <- Site_table[which(Site_table$category == type),]
  list <- rownames(Site_table_red)
  table_Methanogens <- select(OTU_table_Methanogens, list)
  data_frames[[paste(type, " Methanogens", sep = "")]] <- table_Methanogens
  table_Methanogens_cat <- select(OTU_table_Methanogens_cat, list)
  data_frames[[paste(type, " Methanogens (cat.)", sep = "")]] <- table_Methanogens_cat
  table_Methanotrophs <- select(OTU_table_Methanotrophs, list)
  data_frames[[paste(type, " Methanotrophs", sep = "")]] <- table_Methanotrophs
  table_Methanotrophs_cat <- select(OTU_table_Methanotrophs_cat, list)
  data_frames[[paste(type, " Methanotrophs (cat.)", sep = "")]] <- table_Methanotrophs_cat
  
  sites_table <- data.frame()
  list_number <- list()
  for (y in unique(Site_table_red$site)) {
    list_site <- rownames(Site_table_red[Site_table_red$site == y,])
    site <- as.data.frame(t(rowMeans(table_Methanogens_cat[list_site])), row.names = y)
    list_number[[y]] <- nrow(Site_table_red[Site_table_red$site == y,])
    sites_table <- rbind(sites_table, site)
  }
  sites_table$n <- as.character(list_number)
  sites_table <- as.data.frame(t(sites_table))
  data_frames_sitecat[[paste(type, " Methanogens (sitecat.)", sep = "")]] <- sites_table
  
  sites_table <- data.frame()
  list_number <- list()
  for (y in unique(Site_table_red$site)) {
    list_site <- rownames(Site_table_red[Site_table_red$site == y,])
    site <- as.data.frame(t(rowMeans(table_Methanotrophs_cat[list_site])), row.names = y)
    list_number[[y]] <- nrow(Site_table_red[Site_table_red$site == y,])
    sites_table <- rbind(sites_table, site)
  }
  sites_table$n <- as.character(list_number)
  sites_table <- as.data.frame(t(sites_table))
  data_frames_sitecat[[paste(type, " Methanotrophs (sitecat.)", sep = "")]] <- sites_table
  
  table_Methanogens_types <- select(OTU_table_Methanogens_types, list)
  data_frames_Methanogens_types[[paste(type, " Methanogens types", sep = "")]] <- table_Methanogens_types
  
  table_Methanotrophs_types <- select(OTU_table_Methanotrophs_types, list)
  data_frames_Methanotrophs_types[[paste(type, " Methanotrophs types", sep = "")]] <- table_Methanotrophs_types
  
  sites_table_Methanogens_types <- data.frame()
  list_number <- list()
  for (y in unique(Site_table_red$site)) {
    list_site <- rownames(Site_table_red[Site_table_red$site == y,])
    site <- as.data.frame(t(rowMeans(table_Methanogens_types[list_site])), row.names = y)
    list_number[[y]] <- nrow(Site_table_red[Site_table_red$site == y,])
    sites_table_Methanogens_types <- rbind(sites_table_Methanogens_types, site)
  }
  sites_table_Methanogens_types$n <- as.character(list_number)
  sites_table_Methanogens_types <- as.data.frame(t(sites_table_Methanogens_types))
  data_frames_Methanogens_types_sitecat[[paste(type, " Methanogens types (sitecat.)", sep = "")]] <- sites_table_Methanogens_types
  
  sites_table_Methanotrophs_types <- data.frame()
  list_number <- list()
  for (y in unique(Site_table_red$site)) {
    list_site <- rownames(Site_table_red[Site_table_red$site == y,])
    site <- as.data.frame(t(rowMeans(table_Methanotrophs_types[list_site])), row.names = y)
    list_number[[y]] <- nrow(Site_table_red[Site_table_red$site == y,])
    sites_table_Methanotrophs_types <- rbind(sites_table_Methanotrophs_types, site)
  }
  sites_table_Methanotrophs_types$n <- as.character(list_number)
  sites_table_Methanotrophs_types <- as.data.frame(t(sites_table_Methanotrophs_types))
  data_frames_Methanotrophs_types_sitecat[[paste(type, " Methanotrophs types (sitecat.)", sep = "")]] <- sites_table_Methanotrophs_types
  
  row_total_Methanogens <- as.data.frame(t(rowMeans(table_Methanogens_cat)), row.names = type)
  row_total_Methanogens$n <- as.character(ncol(table_Methanogens_cat))
  data_total_Methanogens <- rbind(data_total_Methanogens, row_total_Methanogens)
  row_total_Methanotrophs <- as.data.frame(t(rowMeans(table_Methanotrophs_cat)), row.names = type)
  row_total_Methanotrophs$n <- as.character(ncol(table_Methanotrophs_cat))
  data_total_Methanotrophs <- rbind(data_total_Methanotrophs, row_total_Methanotrophs)
}
data_total_Methanogens <- as.data.frame(t(data_total_Methanogens))
data_total_Methanotrophs <- as.data.frame(t(data_total_Methanotrophs))

save(data_frames, file="data_frames")
save(data_frames_sitecat, file = "data_frames_sitecat")
save(data_frames_Methanogens_types, file = "data_frames_Methanogens_types")
save(data_frames_Methanogens_types_sitecat, file = "data_frames_Methanogens_types_sitecat")
save(data_frames_Methanotrophs_types, file = "data_frames_Methanotrophs_types")
save(data_frames_Methanotrophs_types_sitecat, file = "data_frames_Methanotrophs_types_sitecat")
save(data_total_Methanogens, file = "data_total_Methanogens")
save(data_total_Methanotrophs, file = "data_total_Methanotrophs")

data_frames_true_false <- grepl("cat", names(data_frames))
names_to_use <- names(data_frames)[data_frames_true_false == TRUE]
names_to_use_ii <- names(data_frames)[data_frames_true_false == FALSE]

data_frames_cat <- list()
for (i in names_to_use) {
  result <- data_frames[[i]]
  data_frames_cat[[i]] <- result
}
save(data_frames_cat, file = "data_frames_cat")

data_frames_spec <- list()
for (i in names_to_use_ii) {
  result <- data_frames[[i]]
  data_frames_spec[[i]] <- result
}
save(data_frames_spec, file = "data_frames_spec")

save(Site_table,file = "Site_table_new")
