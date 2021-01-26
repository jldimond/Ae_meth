setwd("~/Documents/Projects/Anthopleura/Methylation-splicing/Ae_meth/analyses/")

library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

#import methylation files produced by nanopolish with threshold of 3 CpG counts per locus
lst <- list.files(path="~/Documents/Projects/Anthopleura/Methylation-splicing/Ae_meth/analyses/", 
                  pattern = "thresh.tsv", full.names=TRUE)
filelist <- lapply(lst, read.table)
names(filelist) <- list.files(path="~/Documents/Projects/Anthopleura/Methylation-splicing/Ae_meth/analyses/", 
                              pattern = "thresh.tsv", full.names=FALSE)

#merge/join dataframes in list according to first three columns (contig, start, end)
mergedData <- Reduce(function(...) merge(..., by = c("V1","V2","V3")), filelist)

#just get columns of interest and provide headers
mergedData2 <- mergedData[,c(1:7,14:17,24:27,9:12,19:22,29:32)]
colnames(mergedData2) <- c("chrom", "start", "end",	"num_cpgs_bar1", "called_bar1",	"called_meth_bar1",	"meth_freq_bar1",
                           "num_cpgs_bar3", "called_bar3",	"called_meth_bar3",	"meth_freq_bar3",
                           "num_cpgs_bar5", "called_bar5",	"called_meth_bar5",	"meth_freq_bar5",
                           "num_cpgs_bar2", "called_bar2",	"called_meth_bar2",	"meth_freq_bar2",
                           "num_cpgs_bar4", "called_bar4",	"called_meth_bar4",	"meth_freq_bar4",
                           "num_cpgs_bar6", "called_bar6",	"called_meth_bar6",	"meth_freq_bar6")

#calculate and plot means
meth_means <- colMeans(subset(mergedData2, select = c(meth_freq_bar1,meth_freq_bar3,meth_freq_bar5,
                                                      meth_freq_bar2,meth_freq_bar4,meth_freq_bar6)))
plot(meth_means)

asym_mean <- mean(meth_means[1:3])
sym_mean <- mean(meth_means[4:6])

asym_sd <- sd(meth_means[1:3])
sym_sd <- sd(meth_means[4:6])

st.err <- function(x) {
  sd(x)/sqrt(length(x))
}

asym_se <- st.err(meth_means[1:3])
sym_se <- st.err(meth_means[4:6])


summary_meth <- as.data.frame(cbind(c("Aposymbiotic", "Symbiotic"), c(asym_mean, sym_mean), c(asym_se, sym_se)))
summary_meth$V2 <- as.numeric(as.character(summary_meth$V2))
summary_meth$V3 <- as.numeric(as.character(summary_meth$V3))


dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = summary_meth$V2 + summary_meth$V3,
              ymin = summary_meth$V2 - summary_meth$V3)

p <- ggplot(data = summary_meth, aes(x = V1, y = V2, fill = V1)) +
  geom_bar(stat = "identity", position = dodge, fill = c("#DE77AE", "#7FBC41")) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1)))+
  theme_bw() +
  labs(x="",y="Methylation frequency") +
  theme(legend.position = "none")  +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 12, color = "black")) +
  theme(axis.ticks = element_line(color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(size = 1)) 
  theme(text = element_text(size = 20)) 


#Overall variance homogeneity test and ttest
var.test(meth_means[1:3], meth_means[4:6])
t.test(meth_means[1:3], meth_means[4:6], paired = FALSE)

#get only meth freq columns
meth_cols <- mergedData2[,c(7,11,15,19,23,27)]

#ttest cannot deal with all zeros, so find these instances and remove them
indx <- !rowSums(meth_cols[,1]==meth_cols)==ncol(meth_cols)
meth_cols2 <- meth_cols[indx,]

#get rows with "perfect data" (instances of identical values in each of two sets of replicates); ttest cannot handle these
perfect1 <- as.integer(as.logical(apply(meth_cols2[,1:3], 1, function(i) length(unique(i)) > 1)))
perfect2 <- as.integer(as.logical(apply(meth_cols2[,4:6], 1, function(i) length(unique(i)) > 1)))
perfect <- cbind(perfect1, perfect2)
indx2 <- which(rowSums(perfect) == 0)
indx3 <- row.names(rowSums(perfect) == 0)
#drop these rows
methcols3 <- meth_cols2[!seq_len(nrow(meth_cols2)) %in% indx2, ]

#perform ttest on all rows
methcols3$stat <- sapply(1:nrow(methcols3), function(i) t.test((methcols3[i,c(1,2,3)]), ##note values changed
                                                               (methcols3[i,c(4,5,6)]))["p.value"])

methcols3$stat <- as.numeric(methcols3$stat)

#FDR (Benjamini-Hochberg) correction of P-values
methcols3$padjust <- p.adjust(methcols3$stat, method = "BH", n = length(methcols3$stat))

#number of p values <0.05
sigp <- sum(methcols3$padjust <= 0.05)
  
#no values less than 0.05. But, get "perfect data"...
sigmethindx <- c("4249", "11202", "11203", "19933", "66251", "164020")

sig_meth_data <- mergedData2[sigmethindx,]

#**ctg1212:11000-25000, AIPGENE15872	Swiss-Prot	sp|Q61493|DPOLZ_MOUSE, DNA polymerase zeta catalytic subunit-like [Exaiptasia pallida], 4e-06, 54.9% ident, 
#DNA polymerase zeta catalytic subunit-like isoform [Acropora millepora], 5e-138, 71.43% ident 

#PCA
pca <- prcomp(meth_cols, scale=T)
scores <- data.frame(pca$rotation)
scores2 <- cbind(c("Apo", "Apo", "Apo", "Sym", "Sym", "Sym"), scores)
colnames(scores2)[1] <- "SymState"

#component loadings
(pca$sdev)^2 / sum(pca$sdev^2) 
# Or cumulative 
cumsum((pca$sdev)^2) / sum(pca$sdev^2) 

pcaplot <- ggplot(scores2, aes(PC1, PC2, colour = SymState)) + 
  theme_bw() +
  geom_point(size = 5) +
  scale_colour_manual(values = c("#DE77AE", "#7FBC41")) +
  xlab("PC1 62.2%") +
  ylab("PC2 9.4%")+
  theme(legend.position="none") +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 12, color = "black")) +
  theme(axis.ticks = element_line(color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(size = 1)) 
  theme(text = element_text(size = 20)) 
  
pcaplot

#plot barplot and PCA together
library(cowplot)
plot_grid(p, pcaplot, labels=c("A", "B"), ncol = 2, nrow = 1)


#ctg1814:1-2500 no exonerate hits, no mapped ONT seq, 1 mapped Aele seq
#ctg2526:400-2000 no exonerate hits, no mapped ONT seq, 7 mapped Aele seq
#ctg1857:8000-11000 one mapped ONT seq AIPGENE7907	TrEMBL	tr|A0A015K5V9|A0A015K5V9_9GLOM	Pif1p, AIPGENE22229	Swiss-Prot	sp|Q9UUA2|PIF1_SCHPO, AIPGENE26387	Swiss-Prot	sp|Q5AXT5|PIF1_EMENI, AIPGENE23631	TrEMBL	tr|A7SKZ0|A7SKZ0_NEMVE	Predicted protein, AIPGENE25249	Swiss-Prot	sp|Q6AZB8|HARB1_DANRE
#ctg2890:14000-15500 no mapped ONT seq, no exonerate hits, no mapped Aele seq
#*ctg107:105500-107500 several mapped seq, protein SFI1 homolog [Exaptaisia pallida] 9e-12, 35.85% ident
#**ctg1212:11000-25000 no mapped ONT reads, AIPGENE15872	Swiss-Prot	sp|Q61493|DPOLZ_MOUSE, DNA polymerase zeta catalytic subunit-like [Exaiptasia pallida], 4e-06, 54.9% ident, DNA polymerase zeta catalytic subunit-like isoform [Acropora millepora], 5e-138, 71.43% ident 
#*ctg139:129000-131000 AIPGENE21864	Swiss-Prot	sp|Q6AYY8|ACATN_RAT, acetyl-coenzyme A transporter 1 [Exaiptasia pallida] 9e-47, 41.96% ident
#ctg882:23000-24000 ignored, no exonerate hits
#ctg662:114000-116000 no exonerate hits, no ONT mapped seq, transmembrane protease serine 9 [Exaiptasia pallida], 6e-05, 95.83 ident; last sequence: splicing factor U2AF 65 kDa subunit [Exaiptasia pallida], 2e-98, 86.5& ident 
#ctg559:49000-52000 AIPGENE20105	Swiss-Prot	sp|O43299|AP5Z1_HUMAN, AP-5 complex subunit seta-1 {Exaiptasia pallida}, 1e-41, 36.22% ident
#ctg485:44000-45000 many aip genes without hits
#*ctg4024:3000-8000 no exonerate hits, hits with retrotransposon elements
#ctg1695:39800-40500 no exonerate hits

