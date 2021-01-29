#This script compares transcriptome CpG O/E to nanopore-derived methylation of reads
#mapped to the transcriptome

setwd("~/Documents/Projects/Anthopleura/Methylation-splicing/Ae_meth/analyses/")

library(ggplot2)
library(ggpubr)

#open nanopore methylation data
Aele_meth_freq <-read.delim("methylation_frequency_thresh10.tsv", sep = "\t", header=T)

#get average meth frequency
Aele_meth_avg <- aggregate(methylated_frequency ~ chromosome, Aele_meth_freq, mean)

#get sum of CpG motifs
Aele_CpG_motifs <- aggregate(num_motifs_in_group ~ chromosome, Aele_meth_freq, sum)

#merge meth frequency and sum motifs
Meth_avg_sum_motifs <-merge(Aele_meth_avg,Aele_CpG_motifs,by="chromosome")

#open file with seq length and CpG counts (length is V4, CpG count is V5)
Aele_comb <-read.delim("comb", header=F)

#just get seq length and CpG counts
Aele_comb2 <- Aele_comb[,c(1,4,5)]

#open CpG O/E data
Aele_cpgoe <-read.delim("ID_CpG", header=F)

#remove extra space in seq id
Aele_cpgoe2 <- as.data.frame(apply(Aele_cpgoe,2,function(x)gsub('\\s+', '',x)))

#change CpG data from factor to numeric
Aele_cpgoe2[,2] <-as.numeric(as.character(Aele_cpgoe2[,2]))

#merge nanopore data with seq length data
temp_merge <-merge(Meth_avg_sum_motifs,Aele_comb2,by.x="chromosome",by.y="V1")

#merge above data with CpG O/E
All_data <-merge(temp_merge,Aele_cpgoe2,by.x="chromosome",by.y="V1")

#rename columns
colnames(All_data)[4] <- "seq_length"
colnames(All_data)[5] <- "num_CpGs"
colnames(All_data)[6] <- "CpG_OE"

#now determine CpG coverage of nanopore data
cpg_cov <- All_data[,3]/All_data[,5]

#bind cpg_cov to dataframe
All_data <- cbind(All_data, cpg_cov)

#now set threshold of 20% CpG coverage to data. This ensures that there is a reasonable amount
#of nanopore data for each transcriptome contig.
All_data_thresh <- All_data[All_data[,7] >= 0.2,]

#now set threshold of 2 for CpG O/E 
All_data_thresh <- All_data_thresh[All_data_thresh[,6] <= 2,]

#scatter plot of nanopore meth data vs. CpG O/E data
plot(All_data_thresh[,2], All_data_thresh[,6])

#linear regression of nanopore meth data vs. CpG O/E data
model <- lm(All_data_thresh[,2] ~ All_data_thresh[,6])
abline(model, col = "green")

#compare fully methylated to fully unmethylated

CpGOE_0 <- as.numeric(All_data_thresh$CpG_OE[All_data_thresh[,2] == 0])
CpGOE_1 <- as.numeric(All_data_thresh$CpG_OE[All_data_thresh[,2] == 1])

var.test(CpGOE_0, CpGOE_1)
t.test(CpGOE_0, CpGOE_1, paired = FALSE, var.equal = FALSE)

CpGOE_0 <- cbind(CpGOE_0, rep("Fully unmethylated", length(CpGOE_0)))
CpGOE_1 <- cbind(CpGOE_1, rep("Fully methylated", length(CpGOE_1)))

CpG_01 <- as.data.frame(rbind(CpGOE_1, CpGOE_0))
CpG_01[,1] <- as.numeric(as.character(CpG_01[,1]))
CpG_01$V2 = factor(CpG_01$V2,c("Fully unmethylated","Fully methylated"))

var.test(meth_means[1:3], meth_means[4:6])
t.test(meth_means[1:3], meth_means[4:6], paired = FALSE)

#plotting

# Scatter plot using ggplot and ggpubr, plus boxplot
sp <- ggscatter(All_data_thresh, x = "methylated_frequency", y = "CpG_OE",
                color = "#868686FF", size = 0.4, alpha = 0.1,
                add = "loess", add.params = list(color = "#b35806"),
                cor.coef = TRUE, 
                cor.coeff.args = list(method = "pearson", 
                                      label.x.npc = "center", label.y.npc = "top"),
                ylab = "CpG O/E", xlab = "Methylated frequency")+
  border()                                         
#Violinplot of fully meth and umeth by CpG O/E
viol <- ggplot(CpG_01, aes(x=V2, y=CpGOE_1, color =V2),xlab = "", ylab = "CpG O/E", show.legend = FALSE) + 
  geom_violin(size =1) +
  scale_color_manual(values=c("#f1a340", "#998ec3")) +
  labs(x="", y = "CpG O/E")+
  annotate(geom = "text", label = "italic(p) < 2.2e-16", x=1.5, y=1.5, size = 4, parse = TRUE) +
  theme_bw() +
  theme(legend.position = "none")  +
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 12, color = "black")) +
  theme(axis.ticks = element_line(color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(size = 1)) 


# Plot both
library(cowplot)
plot_grid(sp, viol, labels=c("A", "B"), ncol = 2, nrow = 1)

