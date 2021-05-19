

library(ape)
library(Biostrings)
library(ggplot2)
library(stats)
library(ade4)
library(ape)
library(adegenet)
library(phangorn)
library(stringr)
library(seqinr)

# Obtain SARS Cov-2 Complete Genomes from India
df <- read.csv("NA_Seqrecord_2021_3_6.csv")    # read csv file
df1<-subset(df,df$Nuc_Completeness=="complete")
write.csv(df1, "NA_complete.csv", row.names=FALSE)
df_India <- df1[df1$Geo_Location %like% "India",]
df2 <-df_India[!grepl("USA", df_India$Geo_Location),]
write.csv(df2, "NA_India.csv", row.names = FALSE)
df <- read.csv("NA_India.csv", header = TRUE, sep = ",")
accession <- df$Accession
seq <- read.GenBank(accession, as.character = TRUE)
write.dna(seq, file= paste("NA_India","fasta",sep="."),format = "fasta", append =FALSE,nbcol = 6, colw = 10)

# Data cleaning

fastaFile <- readDNAStringSet("NA_India.fasta", format="fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)
write.csv(df, "df_na_india.csv", row.names = FALSE)
df <- read.csv("df_na_india.csv")
df1 <- df[df$sequence %like% "N",] 
write.csv(df1, "df_na_dirtyseq_india.csv", row.names = FALSE)
fasta <- read.csv("df_na_dirtyseq_india.csv")
fasta.fas<-paste0(">",fasta$seq_name, "|", 
                  fasta$Year_Month, "|",
                  "\n",fasta$sequence,"\n")
writeLines(fasta.fas,"df_na_dirtyseq_india.fasta")
df <- read.csv("df_na_dirtyseq_india.csv")
ID <- df$seq_name
length(ID)
df1 <- read.csv("df_na_india.csv")
rows <- which(df1$seq_name %in% ID)
df2 <- df1[-rows, ]
write.csv(df2, "df_na_cleanseq_india.csv", row.names = FALSE)
fasta <- read.csv("df_na_cleanseq_india.csv")
fasta.fas<-paste0(">",fasta$seq_name, "|", 
                  fasta$Year_Month, "|",
                  "\n",fasta$sequence,"\n")
writeLines(fasta.fas,"df_na_cleanseq_india.fasta")
df <- read.csv("df_na_cleanseq_india.csv")
df1 <- df[!duplicated(df$sequence),]
nrow(df1)
write.csv(df1, "df_na_cleanseq_unique_india.csv", row.names = FALSE)

fasta <- read.csv("df_na_cleanseq_unique_india.csv")
fasta.fas<-paste0(">",fasta$SN,
                  "\n",fasta$sequence,"\n")
writeLines(fasta.fas,"df_na_cleanseq_unique_india.fasta")

#Unique seq individual CDS

setwd()

#remove redundant sequences
filenames <- list.files(pattern="*_df_seq.csv") 
genename <- gsub("_df_seq.csv", "", filenames) 
for (i in 1:length(filenames)) {
  tmp <- read.csv(filenames[i])
  row <- which(tmp$seqname == 'Ref')
  tmp <- tmp[c(row, (1:nrow(tmp))[-row]),] 
  tmp <- tmp[!duplicated(tmp$seq),]   
  out_file <- paste0(genename[i], "_df_UniqueSeq", ".csv")
  write.csv(tmp, out_file, row.names=FALSE)
}

filenames <- list.files(pattern="*_df_UniqueSeq.csv") 
for (i in 1:length(filenames)) {
  tmp <- read.csv(filenames[i])
  tmp$seq <- substr(tmp$seq, 1, nchar(tmp$seq)-3)
  out_file <- paste0(genename[i], "_df_UniqueSeq_nostop", ".csv")
  write.csv(tmp, out_file, row.names=FALSE)
}

filenames <- list.files(pattern="*_df_UniqueSeq_nostop.csv") 
df <- data.frame(matrix(ncol = 4, nrow = 12))
x = c("Gene", "Unique_Seq", "Gene_Len", "Norm_Variation")
colnames(df) <- x
df[1:12, 1] = genename

for (i in 1:length(filenames)) {
  tmp <- read.csv(filenames[i])
  x = nrow(tmp)
  y = nchar(tmp [1, 2])
  z = (x/(y/1000))
  df[i, 2] = x    
  df[i, 3] = y    
  df[i, 4] = z    
}
write.csv(df, "unique_genes.csv", row.names = FALSE)
df <- read.csv("unique_genes.csv")
Gene <- df$Gene 
Unique_Seq <- df$Unique_Seq
Gene_Length <- df$Gene_Len

Normalized_Variation <- df$Norm_Variation
ggplot(df, mapping=aes(x=Gene_Length , y = Normalized_Variation, color= Gene))+ 
  geom_point()+
  geom_smooth (method ="loess")+
  ggtitle("Normalized unique sequences in different coding regions in 517 SARS-CoV-2 Genomes")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle=90, size = 8, hjust=1))+
  facet_wrap(facets = vars(Gene))

filenames <- list.files(pattern="*_df_UniqueSeq_nostop.csv") 
genename <- gsub("_df_UniqueSeq_nostop.csv", "", filenames) 

for (i in 1:length(filenames)) {
  tmp <- read.csv(filenames[i])
  tmp.fas<-paste0(">",tmp$seqname,
                  "\n",tmp$seq,"\n")
  out_file <- paste0(genename[i], "_Uni_NoStop_aln", ".fasta")
  writeLines(tmp.fas, out_file)
}

#selection dat
setwd()
filenames <- list.files(pattern="*_df_seq.csv") 
genename <- gsub("_df_seq.csv", "", filenames) 

for (i in 1:length(filenames)) {
  tmp <- read.csv(filenames[i])
  tmp$seq = str_replace_all(tmp$seq, "-", "n")  
  row <- which(tmp$seqname == 'Ref')
  tmp <- tmp[c(row, (1:nrow(tmp))[-row]),] 
  tmp$seq <- substr(tmp$seq, 1, nchar(tmp$seq)-3) 
  out_file <- paste0(genename[i], "_df_seq_ref1_n", ".csv")
  write.csv(tmp, out_file, row.names=FALSE)
}

filenames <- list.files(pattern="*_df_seq_ref1_n.csv") 
for (i in 1:length(filenames)) {
  tmp <- read.csv(filenames[i])
  tmp.fas<-paste0(">",tmp$seqname,
                  "\n",tmp$seq,"\n")
  out_file <- paste0(genename[i], ".fasta")
  writeLines(tmp.fas, out_file)
}

tmp <- read.csv("ORF1ab_df_seq_ref1_n.csv")
tmp$seq <- substr(tmp$seq, 1, nchar(tmp$seq)-2)
write.csv(tmp, "ORF1ab_mo3.csv", row.names=FALSE)
df <- read.csv("ORF1ab_mo3.csv")
df.fas<-paste0(">",df$seqname,
               "\n",df$seq,"\n")
writeLines(df.fas, "ORF1ab_mo3.fasta")

aln <- read.alignment(file = "Fas_Muscle_timewise.fasta", format = "fasta")
df <- as.matrix(aln)
df1 <- as.data.frame(t(df))

nsp1_10 <- df1[266:13441,]
write.csv(nsp1_10, "nsp1_10_matrix.csv", row.names = FALSE)

RdRp_nsp16 <-  df1[13468:21552,]
write.csv(RdRp_nsp16, "RdRp_nsp16_matrix.csv", row.names = FALSE)

filenames <- list.files(pattern="*.matrix.csv") 
genename <- gsub("_matrix.csv", "", filenames) 
for (i in 1:length(filenames)) {
  tmp <- read.csv(filenames[i])
  tmp <- apply(tmp, 2, paste, collapse ="") 
  out_file <- paste0(genename[i], "_df", ".csv")
  write.csv(tmp, out_file, row.names=TRUE)
}

filenames <- list.files(pattern="*_df.csv") 
for (i in 1:length(filenames)) {
  tmp <- read.csv(filenames[i])
  colnames(tmp)[1] <- "seqname"
  colnames(tmp)[2] <- "seq"
  row <- which(tmp$seqname == 'Ref')
  tmp <- tmp[c(row, (1:nrow(tmp))[-row]),] 
  out_file <- paste0(genename[i], "_df_seq", ".csv")
  write.csv(tmp, out_file, row.names=FALSE)
}

filenames <- list.files(pattern="*_df_seq.csv") 
genename <- gsub("_df_seq.csv", "", filenames)
for (i in 1:length(filenames)) {
  tmp <- read.csv(filenames[i])
  tmp.fas<-paste0(">",tmp$seqname,
                  "\n",tmp$seq,"\n")
  out_file <- paste0(genename[i], "_aln", ".fasta")
  writeLines(tmp.fas, out_file)
}

