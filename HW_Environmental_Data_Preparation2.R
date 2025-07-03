#reads in the three files of environmental data and combines them into one. Creates a table with pH and water content for each individual sample ('env_data') as well as a list of 'sitecat' tables with soil horizon means within sites and site means within soil horizons ('data_frames_env_sitecat') and a table with overall site and soil horizon means ('table_means')

library(vegan)
library(tidyverse)
library(writexl)

load("../Erik_thesis/OTU_table_rel")
load("Site_table_new")
Site_table$rownames <- rownames(Site_table)
Site_table$category <- as.character(Site_table$category)
Site_table$site <- as.character(Site_table$site)

environmental_data <- read.csv("../Erik_thesis/environmental_data_cryocarb.csv", header = TRUE,fileEncoding="latin1",stringsAsFactors = F)
environmental_data_II <- read.csv("../Erik_thesis/Herschel_dataset_Physicochemical_MV.csv", header = TRUE,fileEncoding="latin1",stringsAsFactors = F)
environmental_data_III <- read.csv("../Erik_thesis/disko-parameters_arranged_SampleIDs.csv", header = TRUE,fileEncoding="latin1",stringsAsFactors = F)

for(i in c(7:nrow(environmental_data))){
  environmental_data$Full.Code[i] <- paste(substr(environmental_data$Full.Code[i],1,2),substr(environmental_data$Full.Code[i],4,20),sep = "")
}

# substr(environmental_data$Full.Code[22],1,2)
# substr(environmental_data$Full.Code[22],4,20)

env_data <- data.frame(environmental_data$pH..H2O.[7:nrow(environmental_data)], environmental_data$Watercontent[7:nrow(environmental_data)],stringsAsFactors = F)
rownames(env_data) <- environmental_data$Full.Code[7:nrow(environmental_data)]
colnames(env_data) <- c("pH","Watercontent")

env_data_II <- data.frame(environmental_data_II$pH..H2O., environmental_data_II$Moisture....,stringsAsFactors = F)
rownames(env_data_II) <- environmental_data_II$Sample.ID.1
colnames(env_data_II) <- c("pH","Watercontent")

env_data_III <- data.frame(environmental_data_III$pH,stringsAsFactors = F)
rownames(env_data_III) <- environmental_data_III$Sample.Ids
env_data_III$WC <- rep(0)
colnames(env_data_III) <- c("pH","Watercontent")


env_data <- rbind(env_data, env_data_II, env_data_III)

pH <- env_data$pH
pH <- gsub(",",".",pH)
env_data$pH <- pH

Watercontent <- env_data$Watercontent
Watercontent <- gsub(",",".",Watercontent)
env_data$Watercontent <- Watercontent

for(i in c(1:ncol(env_data))){
  env_data[,i] <- as.numeric(env_data[,i])
}

for(i in c(1:ncol(env_data))){
  for (y in c(1:nrow(env_data))) {
    if (is.na(env_data[y,i])){
      env_data[y,i] <- 0
    }
  }
}
env_data <- env_data[which(rowSums(env_data)>0),]
for(i in c(1:ncol(env_data))){
  for (y in c(1:nrow(env_data))) {
    if (env_data[y,i] == 0){
      env_data[y,i] <- NA
    }
  }
}
env_data$site <- rep(0)
env_data$type <- rep(0)
env_data$rownames <- rownames(env_data)
Site_table$rownames <- rownames(Site_table)

for (i in c(1:nrow(env_data))) {
  True_False <- Site_table$rownames == env_data$rownames[i]
  if(any(True_False) == TRUE){
    env_data$site[i] <- Site_table$site[True_False == TRUE]
  }
}
for (i in c(1:nrow(env_data))) {
  True_False <- Site_table$rownames == env_data$rownames[i]
  if(any(True_False) == TRUE){
    env_data$type[i] <- Site_table$category[True_False == TRUE]
  }
}
env_data <- env_data[which(env_data$site != 0),]
env_data <- env_data[,1:4]

data_frames_env_sitecat <- list()
table_means <- data.frame()
for(i in c("AriMas","Beaufort coast","Cherskiy","Herschel","Logata","Tazovskiy","Zackenberg")){
  table <- data.frame()
  env_subset <- env_data[which(env_data$site == i),]
  
  env_subset <- as.data.frame(t(env_subset))
  env_subset$means <- rep(0)
  for (n in c(1:2)) {
    row <- env_subset[n,1:ncol(env_subset)-1]
    env_subset$means[n] <- mean(as.numeric(row[is.na(row) == FALSE]))
  }
  env_subset <- as.data.frame(t(env_subset))
  row <- env_subset[nrow(env_subset),1:2]
  rownames(row) <- i
  table_means <- rbind(table_means, row)

  if(i == "Tazovskiy"){
    for(y in c("cryoOM","organic layer","subsoil","topsoil")){
      env_subset_type <- env_subset[which(env_subset$type == y),]
      env_subset_type <- as.data.frame(t(env_subset_type))
      env_subset_type$means <- rep(0)
      for (n in c(1:2)) {
        row <- env_subset_type[n,1:ncol(env_subset_type)-1]
        env_subset_type$means[n] <- mean(as.numeric(row[is.na(row) == FALSE]))
      }
      env_subset_type <- as.data.frame(t(env_subset_type))
      row <- env_subset_type[nrow(env_subset_type),1:2]
      rownames(row) <- c(y)
      table <- rbind(table, row)
    }
  } else {
    for(y in c("cryoOM","organic layer","permafrost","subsoil","topsoil")){
      env_subset_type <- env_subset[which(env_subset$type == y),]
      env_subset_type <- as.data.frame(t(env_subset_type))
      env_subset_type$means <- rep(0)
      for (n in c(1:2)) {
        row <- env_subset_type[n,1:ncol(env_subset_type)-1]
        env_subset_type$means[n] <- mean(as.numeric(row[is.na(row) == FALSE]))
      }
      env_subset_type <- as.data.frame(t(env_subset_type))
      row <- env_subset_type[nrow(env_subset_type),1:2]
      rownames(row) <- c(y)
      table <- rbind(table, row)
    }
  }
  data_frames_env_sitecat[[i]] <- table 
}

for(i in c("cryoOM","organic layer","permafrost","subsoil","topsoil")){
  table <- data.frame()
  env_subset <- env_data[which(env_data$type == i),]
  
  env_subset <- as.data.frame(t(env_subset))
  env_subset$means <- rep(0)
  for (n in c(1:2)) {
    row <- env_subset[n,1:ncol(env_subset)-1]
    env_subset$means[n] <- mean(as.numeric(row[is.na(row) == FALSE]))
  }
  env_subset <- as.data.frame(t(env_subset))
  row <- env_subset[nrow(env_subset),1:2]
  rownames(row) <- i
  table_means <- rbind(table_means, row)
  
  if(i == "permafrost"){
    for(y in c("AriMas","Cherskiy","Logata","Zackenberg")){
      env_subset_site <- env_subset[which(env_subset$site == y),]
      env_subset_site <- as.data.frame(t(env_subset_site))
      env_subset_site$means <- rep(0)
      for (n in c(1:2)) {
        row <- env_subset_site[n,1:ncol(env_subset_site)-1]
        env_subset_site$means[n] <- mean(as.numeric(row[is.na(row) == FALSE]))
      }
      env_subset_site <- as.data.frame(t(env_subset_site))
      row <- env_subset_site[nrow(env_subset_site),1:2]
      rownames(row) <- c(y)
      table <- rbind(table, row)
    }
  } else {
    for(y in c("AriMas","Cherskiy","Logata","Tazovskiy","Zackenberg")){
      env_subset_site <- env_subset[which(env_subset$site == y),]
      env_subset_site <- as.data.frame(t(env_subset_site))
      env_subset_site$means <- rep(0)
      for (n in c(1:2)) {
        row <- env_subset_site[n,1:ncol(env_subset_site)-1]
        env_subset_site$means[n] <- mean(as.numeric(row[is.na(row) == FALSE]))
      }
      env_subset_site <- as.data.frame(t(env_subset_site))
      row <- env_subset_site[nrow(env_subset_site),1:2]
      rownames(row) <- c(y)
      table <- rbind(table, row)
    }
  }
  data_frames_env_sitecat[[i]] <- table 
}

save(environmental_data, file = "environmental_data")
save(env_data, file = "env_data")
save(data_frames_env_sitecat, file = "data_frames_env_sitecat")
save(table_means, file = "table_means")
