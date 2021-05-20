
library(msa)
library(seqinr)
library(stats)
library(ape)
library(phangorn)
library(dplyr)
library(readr)

fas <- readAAStringSet("df_na_cleanseq_unique_india.fasta")
fas_Muscle <- msa(fas, "Muscle", verbose = TRUE) 
fas_Muscle
print(fas_Muscle, show="complete")

aln_to_fas <- function(alignment, filename) {
  sink(filename)
  
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  
  }
  
  sink(NULL)
}

aln_to_fas (fas_Muscle, 'fas_Muscle_timewise.fasta') # to convert the file into fasta

#Make a new folder to store files for cordinate-specific alignments
setwd()
list.files()
aln <- read.alignment(file = "Fas_Muscle_timewise.fasta", format = "fasta")
df <- as.matrix(aln)
df1 <- as.data.frame(t(df))
ncol(df1)   #517 columns
nrow(df1)   #rows 29871
write.csv(df1, "final_aln_matrix.csv")

# matrix subset coordinate wise
df <- read.csv("final_aln_matrix.csv")
ncol(df)  #518 including the first column of coordinate
nrow(df)  #29871 characters

#Subset aln for ORF1ab (266 to 21555)
ORF1ab <- df[266:21555,]
ncol(ORF1ab)
nrow(ORF1ab)
write.csv(ORF1ab, "ORF1ab_matrix.csv", row.names = FALSE)

#Subset aln for spike (21563 to 25384)
spike <- df[21563:25384,]
ncol(spike)
nrow(spike)
write.csv(spike, "spike_matrix.csv", row.names = FALSE)

#Subset aln for ORF3a (25393 to 26220)
ORF3a <- df[25393:26220,]
ncol(ORF3a)
nrow(ORF3a)
write.csv(ORF3a, "ORF3a_matrix.csv", row.names = FALSE)

#Subset aln for envelope (26245 to 26472)
envelope <- df[26245:26472,]
ncol(envelope)
nrow(envelope)
write.csv(envelope, "envelope_matrix.csv", row.names = FALSE)

#Subset aln for mem_glycoprotein (26523 to 27191)
mem_glycoprotein <- df[26523:27191,]
ncol(mem_glycoprotein)
nrow(mem_glycoprotein)
write.csv(mem_glycoprotein, "mem_glycoprotein_matrix.csv", row.names = FALSE)

#Subset aln for ORF6 (27202:27387)
ORF6 <- df[27202:27387,]
ncol(ORF6)
nrow(ORF6)
write.csv(ORF6, "ORF6_matrix.csv", row.names = FALSE)

#Subset aln for ORF7a (27394:27759)
ORF7a <- df[27394:27759,]
ncol(ORF7a)
nrow(ORF7a)
write.csv(ORF7a, "ORF7a_matrix.csv", row.names = FALSE)

#Subset aln for ORF7b (27756:27887)
ORF7b <- df[27756:27887,]
ncol(ORF7b)
nrow(ORF7b)
write.csv(ORF7b, "ORF7b_matrix.csv", row.names = FALSE)

#Subset aln for ORF8 (27894:28259)
ORF8 <- df[27894:28259,]
ncol(ORF8)
nrow(ORF8)
write.csv(ORF8, "ORF8_matrix.csv", row.names = FALSE)

#Subset aln for nucleocapsid (28274:29533)
nucleocapsid <- df[28274:29533,]
ncol(nucleocapsid)
nrow(nucleocapsid)
write.csv(nucleocapsid, "nucleocapsid_matrix.csv", row.names = FALSE)

#Subset aln for ORF10 (29558:29674)
ORF10 <- df[29558:29674,]
ncol(ORF10)
nrow(ORF10)
write.csv(ORF10, "ORF10_matrix.csv", row.names = FALSE)

#Copy all the files in a new folder 
setwd()
list.files()
filenames <- list.files(pattern="*.csv") 
genename <- gsub("_matrix.csv", "", filenames) 
for (i in 1:length(filenames)) {
  tmp <- read.csv(filenames[i])
  tmp <- subset(tmp, select = -X)
  tmp <- apply(tmp, 2, paste, collapse ="") 
  out_file <- paste0(genename[i], "_df", ".csv")
  write.csv(tmp, out_file, row.names=TRUE)
}

filenames <- list.files(pattern="*_df.csv") 
genename <- gsub("_matrix.csv", "", filenames) 
for (i in 1:length(filenames)) {
  tmp <- read.csv(filenames[i])
  colnames(tmp)[1] <- "seqname"
  colnames(tmp)[2] <- "seq"
  out_file <- paste0(genename[i], "_df_seq", ".csv")
  write.csv(tmp, out_file, row.names=FALSE)
}

filenames <- list.files(pattern="*_df_seq.csv") 
for (i in 1:length(filenames)) {
  tmp <- read.csv(filenames[i])
  tmp.fas<-paste0(">",tmp$seqname,
                  "\n",tmp$seq,"\n")
  out_file <- paste0(genename[i], "_aln", ".fasta")
  writeLines(tmp.fas, out_file)
}


filenames <- list.files(pattern="*.fasta") 
for (i in 1:length(filenames)) {
  tmp <- read.dna(filenames[i], format = "fasta")
  D <- dist.dna(tmp, model = "TN93")
  distance_matrix <- as.matrix(D)[1:517, "Ref", drop=FALSE] 
  out_file <- paste0(genename[i], "_dis_matrix", ".csv")
  write.csv(distance_matrix, out_file, row.names=TRUE)
}

setwd()
list.files()
df <- list.files(full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

write.csv(df, "All_dist_matrix.csv", row.names = FALSE)

Time <- df$Time 
Gene <- df$Gene

ggplot(df, mapping=aes(x=Time , y = Distance, color= Time))+ 
  geom_point()+
  geom_smooth (method ="loess")
ggtitle("Genetic distance of different regions of SARS-CoV-2")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle=90, size = 8, hjust=1))+
  facet_wrap(facets = vars(Gene))

