
library(adegenet)
library(ape)


setwd() 

ent_genome <- read.dna("Fas_Muscle_timewise.fasta", format="fasta")    
ent_genome_gi <- DNAbin2genind(ent_genome, polyThres=0)     

envelope <- read.dna("envelope_Uni_NoStop_aln.fasta", format="fasta")    
envelope_gi <- DNAbin2genind(envelope, polyThres=0)     

mem_glycoprotein <- read.dna("mem_glycoprotein_Uni_NoStop_aln.fasta", format="fasta")    
mem_glycoprotein_gi <- DNAbin2genind(mem_glycoprotein, polyThres=0)     

nucleocapsid <- read.dna("nucleocapsid_Uni_NoStop_aln.fasta", format="fasta")    
nucleocapsid_gi <- DNAbin2genind(nucleocapsid, polyThres=0)     

ORF10 <- read.dna("ORF10_Uni_NoStop_aln.fasta", format="fasta")    
ORF10_gi <- DNAbin2genind(ORF10, polyThres=0)     

ORF1ab <- read.dna("ORF1ab_Uni_NoStop_aln.fasta", format="fasta")    
ORF1ab_gi <- DNAbin2genind(ORF1ab, polyThres=0)     

ORF3a <- read.dna("ORF3a_Uni_NoStop_aln.fasta", format="fasta")    
ORF3a_gi <- DNAbin2genind(ORF3a, polyThres=0)     

ORF6 <- read.dna("ORF6_Uni_NoStop_aln.fasta", format="fasta")    
ORF6_gi <- DNAbin2genind(ORF6, polyThres=0)     

ORF7a <- read.dna("ORF7a_Uni_NoStop_aln.fasta", format="fasta")    
ORF7a_gi <- DNAbin2genind(ORF7a, polyThres=0)     

ORF7b <- read.dna("ORF7b_Uni_NoStop_aln.fasta", format="fasta")    
ORF7b_gi <- DNAbin2genind(ORF7b, polyThres=0)     

ORF8 <- read.dna("ORF8_Uni_NoStop_aln.fasta", format="fasta")    
ORF8_gi <- DNAbin2genind(ORF8, polyThres=0)     

spike <- read.dna("spike_Uni_NoStop_aln.fasta", format="fasta")    
spike_gi <- DNAbin2genind(spike, polyThres=0)     


ent_genome_gi
envelope_gi
mem_glycoprotein_gi
nucleocapsid_gi
ORF10_gi
ORF1ab_gi
ORF3a_gi
ORF6_gi
ORF7a_gi
ORF7b_gi
ORF8_gi
spike_gi

length(ent_genome_gi$loc.n.all)
which(ent_genome_gi$loc.n.all == 2) 
length(which(ent_genome_gi$loc.n.all == 2)) 
length(which(ent_genome_gi$loc.n.all == 3))  
length(envelope_gi$loc.n.all)
which(envelope_gi$loc.n.all == 2)   
length(which(envelope_gi$loc.n.all == 2))  
length(which(envelope_gi$loc.n.all == 3))  
which(mem_glycoprotein_gi$loc.n.all == 2)   
length(which(mem_glycoprotein_gi$loc.n.all == 2))  
length(which(mem_glycoprotein_gi$loc.n.all == 3))  
which(nucleocapsid_gi$loc.n.all == 2)   
length(which(nucleocapsid_gi$loc.n.all == 2))  
length(which(nucleocapsid_gi$loc.n.all == 3))  
which(ORF10_gi$loc.n.all == 2)   
length(which(ORF10_gi$loc.n.all == 2))  
length(which(ORF10_gi$loc.n.all == 3))  
which(ORF1ab_gi$loc.n.all == 2)   
length(which(ORF1ab_gi$loc.n.all == 2))  
length(which(ORF1ab_gi$loc.n.all == 3))  
which(ORF3a_gi$loc.n.all == 2)   
length(which(ORF3a_gi$loc.n.all == 2))  
length(which(ORF3a_gi$loc.n.all == 3))  
which(ORF6_gi$loc.n.all == 2)   
length(which(ORF6_gi$loc.n.all == 2))  
length(which(ORF6_gi$loc.n.all == 3))  
which(ORF7a_gi$loc.n.all == 2)   
length(which(ORF7a_gi$loc.n.all == 2))  
length(which(ORF7a_gi$loc.n.all == 3))  
which(ORF7b_gi$loc.n.all == 2)   
length(which(ORF7b_gi$loc.n.all == 2))  
length(which(ORF7b_gi$loc.n.all == 3))  
which(ORF8_gi$loc.n.all == 2)   
length(which(ORF8_gi$loc.n.all == 2))  
length(which(ORF8_gi$loc.n.all == 3))  
which(spike_gi$loc.n.all == 2)   
length(which(spike_gi$loc.n.all == 2))  
length(which(spike_gi$loc.n.all == 3))  

# % polymorphic sites
100*mean(unlist((ent_genome_gi)$loc.n.all)>1)
100*mean(unlist((envelope_gi)$loc.n.all)>1)
100*mean(unlist((mem_glycoprotein_gi)$loc.n.all)>1)
100*mean(unlist((nucleocapsid_gi)$loc.n.all)>1)
100*mean(unlist((ORF10_gi)$loc.n.all)>1)
100*mean(unlist((ORF1ab_gi)$loc.n.all)>1)
100*mean(unlist((ORF3a_gi)$loc.n.all)>1)
100*mean(unlist((ORF6_gi)$loc.n.all)>1)
100*mean(unlist((ORF7a_gi)$loc.n.all)>1)
100*mean(unlist((ORF7b_gi)$loc.n.all)>1)
100*mean(unlist((ORF8_gi)$loc.n.all)>1)
100*mean(unlist((spike_gi)$loc.n.all)>1)

#tables

df <- ent_genome_gi$tab
write.csv(df, "ent_genome_gi_allele.csv")
df <- envelope_gi$tab
write.csv(df, "envelope_allele.csv")
df <- mem_glycoprotein_gi$tab
write.csv(df, "mem_glycoprotein_allele.csv")
df <- nucleocapsid_gi$tab
write.csv(df, "nucleocapsid_allele.csv")
df <-ORF10_gi$tab
write.csv(df, "ORF10_gi_allele.csv")
df <- ORF1ab_gi$tab
write.csv(df, "ORF1ab_allele.csv")
df <- ORF3a_gi$tab
write.csv(df, "ORF3a_allele.csv")
df <- ORF6_gi$tab
write.csv(df, "ORF6_allele.csv")
df <- ORF7a_gi$tab
write.csv(df, "ORF7a_allele.csv")
df <- ORF7b_gi$tab
write.csv(df, "ORF7b_allele.csv")
df <- ORF8_gi$tab
write.csv(df, "ORF8_allele.csv")
df <- spike_gi$tab
write.csv(df, "spike_allele.csv")


# plot number of alleles per locus
par(mar = c(3,4,4,4)) 
par(mfrow=c(3,4)) #plots in 3 rows & 4 columns
temp <- table(unlist((ent_genome_gi)$loc.n.all))
barplot(temp, main="ent_genome",
        xlab="Number of alleles", ylab="Number of sites", col=heat.colors(4))
temp <- table(unlist((envelope_gi)$loc.n.all))
barplot(temp, main="envelope",
        xlab="Number of alleles", ylab="Number of sites", col=heat.colors(4))
temp <- table(unlist((mem_glycoprotein_gi)$loc.n.all))
barplot(temp, main="mem_glycoprotein",
        xlab="Number of alleles", ylab="Number of sites", col=heat.colors(4))
temp <- table(unlist((nucleocapsid_gi)$loc.n.all))
barplot(temp, main="nucleocapsid",
        xlab="Number of alleles", ylab="Number of sites", col=heat.colors(4))
temp <- table(unlist((ORF10_gi)$loc.n.all))
barplot(temp, main="ORF10",
        xlab="Number of alleles", ylab="Number of sites", col=heat.colors(4))
temp <- table(unlist((ORF1ab_gi)$loc.n.all))
barplot(temp, main="ORF1ab",
        xlab="Number of alleles", ylab="Number of sites", col=heat.colors(4))
temp <- table(unlist((ORF3a_gi)$loc.n.all))
barplot(temp, main="ORF3a",
        xlab="Number of alleles", ylab="Number of sites", col=heat.colors(4))
temp <- table(unlist((ORF6_gi)$loc.n.all))
barplot(temp, main="ORF6",
        xlab="Number of alleles", ylab="Number of sites", col=heat.colors(4))
temp <- table(unlist((ORF7a_gi)$loc.n.all))
barplot(temp, main="ORF7a",
        xlab="Number of alleles", ylab="Number of sites", col=heat.colors(4))
temp <- table(unlist((ORF7b_gi)$loc.n.all))
barplot(temp, main="ORF7b",
        xlab="Number of alleles", ylab="Number of sites", col=heat.colors(4))
temp <- table(unlist((ORF8_gi)$loc.n.all))
barplot(temp, main="ORF8",
        xlab="Number of alleles", ylab="Number of sites", col=heat.colors(4))
temp <- table(unlist((spike_gi)$loc.n.all))
barplot(temp, main="spike",
        xlab="Number of alleles", ylab="Number of sites", col=heat.colors(4))


#distribution
ent_genome     
envelope     
mem_glycoprotein     
nucleocapsid     
ORF10     
ORF1ab     
ORF3a     
ORF6     
ORF7a     
ORF7b     
ORF8     
spike     

#Function
snpposi.plot <- function(...){
  UseMethod("snpposi.plot")
}

## METHOD FOR INTEGER - BASIC METHOD
#' @method snpposi.plot integer
#' @export
snpposi.plot.integer <- function(x, genome.size, smooth=0.1, col="royalblue", alpha=.2,
                                 codon=TRUE, start.at=1, ...){
  ## IF WE REPRESENT DENSITY PER CODON POSITION ##
  if(codon){
    ## define base positions (1/2/3) ##
    codon.posi <- ((2 + x) %% 3) + 1
    fac <- factor(codon.posi, levels=1:3)
    
    ## make ggplot output ##
    out <- ggplot(data.frame(x=x, codon=fac), aes(x=x)) + xlim(0, genome.size)
    theme_update(plot.title = element_text(hjust = 0.5))
    out <- out + geom_density(adjust=smooth, aes(fill=codon, colour=codon),alpha=I(alpha)) + geom_rug(aes(colour=codon),alpha=.7)
    out <- out + labs(x="Nucleotide position", title="Distribution of SNPs in spike")
    out <- out + guides(fill=guide_legend(title="Codon position"), colour=guide_legend(title="Codon position"))
  } else {
    ## OTHERWISE, JUST ONE DENSITY ##
    ## make ggplot output ##
    out <- ggplot(data.frame(x=x), aes(x=x)) + xlim(0, genome.size)
    out <- out + geom_density(adjust=smooth, fill=transp(col,alpha=alpha), colour=col) + geom_rug(colour=col,alpha=.7)
    out <- out + labs(x="Nucleotide position", title="Distribution of SNPs in ent_genome")
  }
  
  ## return ##
  return(out)
} # end snpposi.plot.integer

#' @method snpposi.plot numeric
#' @export
snpposi.plot.numeric <- function(x, ...){
  out <- snpposi.plot(as.integer(x), ...)
  return(out)
}

#' @method snpposi.plot DNAbin
#' @export
snpposi.plot.DNAbin <- function(x, ...){
  out <- snpposi.plot(x=as.integer(seg.sites(x)),
                      genome.size=ncol(x), ...)
  return(out)
} # end snpposi.plot.DNAbin


snpposi.plot(ent_genome, codon=FALSE) #As its aln of whole genome so we can't define the codon position
snpposi.test(ent_genome)
snpposi.plot(envelope, codon=TRUE) #As its aln of CDS so we can define the codon position
snpposi.test(envelope)
snpposi.plot(mem_glycoprotein, codon=TRUE) #As its aln of CDS so we can define the codon position
snpposi.test(mem_glycoprotein)
snpposi.plot(nucleocapsid, codon=TRUE) #As its aln of CDS so we can define the codon position
snpposi.test(nucleocapsid)
snpposi.plot(ORF10, codon=TRUE) #As its aln of CDS so we can define the codon position
snpposi.test(ORF10)
snpposi.plot(ORF1ab, codon=TRUE) #As its aln of CDS so we can define the codon position
snpposi.test(ORF1ab)
snpposi.plot(ORF3a, codon=TRUE) #As its aln of CDS so we can define the codon position
snpposi.test(ORF3a)
snpposi.plot(ORF6, codon=TRUE) #As its aln of CDS so we can define the codon position
snpposi.test(ORF6)
snpposi.plot(ORF7a, codon=TRUE) #As its aln of CDS so we can define the codon position
snpposi.test(ORF7a)
snpposi.plot(ORF7b, codon=TRUE) #As its aln of CDS so we can define the codon position
snpposi.test(ORF7b)
snpposi.plot(ORF8, codon=TRUE) #As its aln of CDS so we can define the codon position
snpposi.test(ORF8)
snpposi.plot(spike, codon=TRUE) #As its aln of CDS so we can define the codon position
snpposi.test(spike)

#Mutations
ent_genome_mut <- findMutations(ent_genome, from=1)
df <- as.matrix(ent_genome_mut)
write.table(df, "ent_genome_mut.txt")
write.csv(df, "ent_genome_mut.csv")
g <- graphMutations(ent_genome, from=1)
g <- gengraph(ent_genome, cutoff=5)
plot(g$graph)

envelope_mut <- findMutations(envelope, from=1)
df <- as.matrix(envelope_mut)
write.table(df, "envelope_mut.txt")
write.csv(df, "envelope_mut.csv")
g <- graphMutations(envelope, from=1)
g <- gengraph(envelope, cutoff=5)
plot(g$graph)

mem_glycoprotein_mut <- findMutations(mem_glycoprotein, from=1)
df <- as.matrix(mem_glycoprotein_mut)
write.table(df, "mem_glycoprotein_mut.txt")
write.csv(df, "mem_glycoprotein_mut.csv")
g <- graphMutations(mem_glycoprotein, from=1)
g <- gengraph(mem_glycoprotein, cutoff=5)
plot(g$graph)

nucleocapsid_mut <- findMutations(nucleocapsid, from=1)
df <- as.matrix(nucleocapsid_mut)
write.table(df, "nucleocapsid_mut.txt")
write.csv(df, "nucleocapsid_mut.csv")
g <- graphMutations(nucleocapsid, from=1)
g <- gengraph(nucleocapsid, cutoff=5)
plot(g$graph)

ORF1ab_mut <- findMutations(ORF1ab, from=1)
df <- as.matrix(ORF1ab_mut)
write.table(df, "ORF1ab_mut.txt")
write.csv(df, "ORF1ab_mut.csv")
g <- graphMutations(ORF1ab, from=1)
g <- gengraph(ORF1ab, cutoff=5)
plot(g$graph)

ORF3a_mut <- findMutations(ORF3a, from=1)
df <- as.matrix(ORF3a_mut)
write.table(df, "ORF3a_mut.txt")
write.csv(df, "ORF3a_mut.csv")
g <- graphMutations(ORF3a, from=1)
g <- gengraph(ORF3a, cutoff=5)
plot(g$graph)

ORF6_mut <- findMutations(ORF6, from=1)
df <- as.matrix(ORF6_mut)
write.table(df, "ORF6_mut.txt")
write.csv(df, "ORF6_mut.csv")
g <- graphMutations(ORF6, from=1)
g <- gengraph(ORF6, cutoff=5)
plot(g$graph)

ORF7a_mut <- findMutations(ORF7a, from=1)
df <- as.matrix(ORF7a_mut)
write.table(df, "ORF7a_mut.txt")
write.csv(df, "ORF7a_mut.csv")
g <- graphMutations(ORF7a, from=1)
g <- gengraph(ORF7a, cutoff=5)
plot(g$graph)

ORF7b_mut <- findMutations(ORF7b, from=1)
df <- as.matrix(ORF7b_mut)
write.table(df, "ORF7b_mut.txt")
write.csv(df, "ORF7b_mut.csv")
g <- graphMutations(ORF7b, from=1)
g <- gengraph(ORF7b, cutoff=5)
plot(g$graph)

ORF8_mut <- findMutations(ORF8, from=1)
df <- as.matrix(ORF8_mut)
write.table(df, "ORF8_mut.txt")
write.csv(df, "ORF8_mut.csv")
g <- graphMutations(ORF8, from=1)
g <- gengraph(ORF8, cutoff=5)
plot(g$graph)


ORF10_mut <- findMutations(ORF10, from=1)
df <- as.matrix(ORF10_mut)
write.table(df, "ORF10_mut.txt")
write.csv(df, "ORF10_mut.csv")
g <- graphMutations(ORF10, from=1)
g <- gengraph(ORF10, cutoff=5)
plot(g$graph)

spike_mut <- findMutations(spike, from=1)
df <- as.matrix(spike_mut)
write.table(df, "spike_mut.txt")
write.csv(df, "spike_mut.csv")

g <- graphMutations(spike, from=1)
g <- gengraph(spike, cutoff=5)
plot(g$graph)


# oci with SNP
loci_envelope <- envelope_gi$loc.n.all
write.csv(loci_envelope, "loci_envelope.csv")
loci_mem_glycoprotein <- mem_glycoprotein_gi$loc.n.all
write.csv(loci_mem_glycoprotein, "loci_mem_glycoprotein.csv")
loci_nucleocapsid <- nucleocapsid_gi$loc.n.all
write.csv(loci_nucleocapsid, "loci_nucleocapsid.csv")
loci_ORF1ab <- ORF1ab_gi$loc.n.all
write.csv(loci_ORF1ab, "loci_ORF1ab.csv")
loci_ORF10 <- ORF10_gi$loc.n.all
write.csv(loci_ORF10, "loci_ORF10.csv")
loci_ORF3a <- ORF3a_gi$loc.n.all
write.csv(loci_ORF3a, "loci_ORF3a.csv")
loci_ORF6 <- ORF6_gi$loc.n.all
write.csv(loci_ORF6, "loci_ORF6.csv")
loci_ORF7a <- ORF7a_gi$loc.n.all
write.csv(loci_ORF7a, "loci_ORF7a.csv")
loci_ORF7b <- ORF7b_gi$loc.n.all
write.csv(loci_ORF7b, "loci_ORF7b.csv")
loci_ORF8 <- ORF8_gi$loc.n.all
write.csv(loci_ORF8, "loci_ORF8.csv")
loci_spike <- spike_gi$loc.n.all
write.csv(loci_spike, "loci_spike.csv")

#PCA 

sum(is.na(ent_genome_gi$tab))
X <- scaleGen(ent_genome_gi, NA.method="mean")
class(X)
dim(X)
par(mar = c(3,3,3,3)) 
par(mfrow=c(1,1))
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3) 
barplot(pca1$eig,main="PCA eigenvalues", col=heat.colors(50))
pca1$eig
pca1$li
pca1$c1
s.label(pca1$li)
title("PCA of ent_genome dataset")
add.scatter.eig(pca1$eig[1:20], 3,1,2)