#dada2
library(dada2); packageVersion("dada2")
library(ShortRead)
library(ggplot2)

path <- "fastq"
fns <- list.files(path)
fns

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#Examine quality profiles of forward and reverse reads
plotQualityProfile(fnFs[[1]])
plotQualityProfile(fnFs[[90]])

plotQualityProfile(fnRs[[1]])
plotQualityProfile(fnRs[[90]])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# use trimleft to remove adapter and primer
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,150),trimLeft = c(25,26),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


#plot errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]] 
dadaRs[[1]]


#merge
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

rm(derepFs);rm(derepRs)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#making sequence table
seqtab <- makeSequenceTable(mergers)

dim(seqtab)


#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)

#saveRDS(seqtab.nochim,"seqtab.nochim.rds")

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.nochim)))

seqtab2 <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% seq(200,300)]

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab2))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#taxonomy

taxa <- assignTaxonomy(seqtab2, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")


colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")
unname(head(taxa))


#representative seqs
uniquesToFasta(seqtab2, fout="uniqueSeqs_Alaska.fasta", 
               ids = paste0("ASV", seq(length(getSequences(seqtab2)))))


seqtab.nochim.new <- seqtab2
colnames(seqtab.nochim.new) <- paste0("ASV", seq(length(getSequences(seqtab2))))

rownames(taxa) <- paste0("ASV", seq(length(getSequences(seqtab2))))

rm(mergers)


sample.data <- read.csv("Mapping_Alaska2019_2021.csv",stringsAsFactors = F,row.names = 1)
ASV.table <- seqtab.nochim.new[rownames(sample.data),]
# merging the resequenced samples
ASV.table <- aggregate(ASV.table,list(sample.data$SampleID),sum)
rownames(ASV.table) <- ASV.table$Group.1
ASV.table <- ASV.table[,-1]

head(colSums(ASV.table))
head(colSums(seqtab.nochim.new))
rowSums(ASV.table)
rowSums(seqtab.nochim.new)

sample.data <- sample.data[1:108,]
rownames(sample.data) <- sample.data$SampleID

taxa <- data.frame(taxa,stringsAsFactors = F)

library(phyloseq)

# Construct phyloseq object (straightforward from dada2 outputs)
physeq <- phyloseq(otu_table(ASV.table, taxa_are_rows=FALSE), 
                   sample_data(sample.data),
                   tax_table(as.matrix(taxa)))

#pre-processing
# physeq <- prune_taxa(taxa_sums(physeq)>1,physeq)

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

physeq.chlo <- subset_taxa(physeq, Order=="Chloroplast")

physeq2 <- pop_taxa(physeq,taxa_names(physeq.chlo))

physeq.mito <- subset_taxa(physeq, Family=="Mitochondria")

physeq2 <- pop_taxa(physeq2,taxa_names(physeq.mito))

sum(sample_sums(physeq2))/sum(sample_sums(physeq))

sample_sums(physeq2)

physeq3 <- prune_samples(sample_sums(physeq2)>1000,physeq2)

## to look at overall microbiome composition
ASV.table2 <- data.frame(otu_table(physeq3))
sample.data2 <- data.frame(sample_data(physeq3))

library(vegan)
ASV.table2.nor <- decostand(ASV.table2, method="hellinger")

distance.bray <- vegdist(ASV.table2.nor, method = "bray")
nmds.bray <- metaMDS(distance.bray, k = 2, trymax = 100)


nmds_df_sites <- as.data.frame(nmds.bray$points)
nmds_df_sites <- cbind(sample.data2,nmds_df_sites)

library(ggplot2)
ggplot()+geom_point(data=nmds_df_sites,aes(MDS1,MDS2,fill=Locality,shape=factor(year)))+
  scale_shape_manual(values = c(21,24))+
  theme_bw()


ggplot()+geom_point(data=nmds_df_sites,aes(MDS1,MDS2,fill=Locality,shape=factor(year)),
                    stroke=1.5,
                    size=8)+
  scale_shape_manual(values = c(21,24))+
  scale_fill_brewer(palette ="Dark2")+
  guides(fill=guide_legend(title="Site",override.aes = list(shape = 21)),shape=guide_legend(title="Year"))+
  theme(axis.title=element_text(size=20),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey",linetype = "dotted"),
        legend.text = element_text(size = 15), legend.key.size = unit(1, "cm"), legend.key = element_rect(fill = NA),legend.title = element_text(size = 20))


ggplot()+geom_point(data=nmds_df_sites,aes(MDS1,MDS2,fill=factor(Depth_proxy,levels = c("0-15cm","16-30cm","31-45cm","46-60cm","61-75cm","76-90cm","91-105cm","106-120cm")),
                                           shape=factor(Locality)),
                    stroke=1.5,
                    size=8)+
  scale_shape_manual(values = c(21,24,22,23))+
  scale_fill_brewer(palette ="Dark2")+
  guides(fill=guide_legend(title="Depth",override.aes = list(shape = 21)),shape=guide_legend(title="Site"))+
  theme(axis.title=element_text(size=20),
        panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size=1),
        panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey",linetype = "dotted"),
        legend.text = element_text(size = 15), legend.key.size = unit(1, "cm"), legend.key = element_rect(fill = NA),legend.title = element_text(size = 20))

##### look at MG & MT
ASV.table3 <- data.frame(otu_table(physeq2))
taxa2 <- data.frame(tax_table(physeq2),stringsAsFactors = F)

MG.names <- c("Methanobacteriales","Methanococcales","Methanopyrales","Methanofastidiosales","Methanocellales","Methanomicrobiales","Methanonatronaechaeales","Methanosarciniales",
              "Methanomassiliicoccales")

taxa.MG <- taxa2[taxa2$Order %in% MG.names,]

MT.names <- c("Methylocapsa","Methylocella","Methylocystis","Methyloferula","Methylosinus","Methylococcales","Methylacidiphilaceae","Candidatus Methylomirabilis")

taxa.MT <- rbind(taxa2[taxa2$Genus %in% MT.names,], taxa2[taxa2$Order %in% MT.names,], taxa2[taxa2$Family %in% MT.names,])

### to save the fastq files
seqtab2 <- data.frame(seqtab2)
seqtab.nochim.new <- data.frame(seqtab.nochim.new)

fasta.seqs <- data.frame(cbind(colnames(seqtab.nochim.new),colnames(seqtab2)))
rownames(fasta.seqs) <- fasta.seqs$X1

MG.fasta.seqs <- fasta.seqs[rownames(taxa.MG),]
seqtab.MG <- seqtab2[,as.character(MG.fasta.seqs$X2)]

uniquesToFasta(as.matrix(seqtab.MG), fout="MG_Alaska.fasta", 
               ids = MG.fasta.seqs$X1)

MT.fasta.seqs <- fasta.seqs[rownames(taxa.MT),]
seqtab.MT <- seqtab2[,as.character(MT.fasta.seqs$X2)]

uniquesToFasta(as.matrix(seqtab.MT), fout="MT_Alaska.fasta", 
               ids = MT.fasta.seqs$X1)

