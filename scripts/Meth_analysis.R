setwd("~/Documents/Projects/Anthopleura/Methylation-splicing/Ae_meth/analyses/")

library(ggplot2)
library(factoextra)

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

meth_means2 <- as.data.frame(cbind(meth_means, c("Aposymbiotic", "Aposymbiotic", "Aposymbiotic", "Symbiotic",
           "Symbiotic", "Symbiotic")))
meth_means2$meth_means <- as.numeric(as.character((meth_means2$meth_means)))


#boxplot of overall methylation frequency
p <- ggplot(meth_means2, aes(x = V2, y = meth_means, fill = meth_means)) +
  geom_boxplot(fill = c("#DE77AE", "#7FBC41")) +
  geom_text(aes(fontface=3), x=1, y=0.06, label="p = 0.073") +
  theme_bw() +
  labs(x="",y="Methylation frequency") +
  theme(legend.position = "none")  +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 12, color = "black")) +
  theme(axis.ticks = element_line(color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(size = 1)) +
  theme(text = element_text(size = 20)) 


#Overall variance homogeneity test and ttest
var.test(meth_means[1:3], meth_means[4:6])
t.test(meth_means[1:3], meth_means[4:6], paired = FALSE)

#get only meth freq columns
meth_cols <- mergedData2[,c(7,11,15,19,23,27)]
row.names(meth_cols) <- mergedData2$chrom

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

#scree plot
eig <- fviz_eig(pca, geom = "bar", barfill = "white", barcolor = "black",
                ggtheme = theme_classic(), main= "", ylab= "% Variance")


#component loadings
(pca$sdev)^2 / sum(pca$sdev^2) 
# Or cumulative 
cumsum((pca$sdev)^2) / sum(pca$sdev^2) 

#plot the pca
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
  theme(panel.border = element_rect(size = 1)) +
  theme(text = element_text(size = 20)) 
  
pcaplot

#add scree plot
pcaplot2 <- pcaplot +
  annotation_custom(ggplotGrob(eig), 
                    ymin = 0.30, ymax= 0.9, xmin=0.405, xmax=0.425)

#plot barplot and PCA together
library(cowplot)
plot_grid(p, pcaplot2, labels=c("A", "B"), ncol = 2, nrow = 1)


