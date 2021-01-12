setwd("~/Documents/Projects/Anthopleura/Methylation&splicing/Anthopleura_genome")

#open nanopore methylation data
Aele_meth <-read.delim("methylation_avg.tsv", sep = "", header=F)

#set min and max thresholds for data
Aele_meth2 <- Aele_meth[Aele_meth$V2 >= 0.001 & Aele_meth$V2 <= 50,]

#change directory 
setwd("~/Documents/Projects/Anthopleura/Methylation&splicing/Aele/")

#open CpG O/E data
Aele_cpg <-read.delim("ID_CpG", header=F)

#set min and max thresholds
Aele_cpg2 <- Aele_cpg[Aele_cpg$V2 >= 0.001 & Aele_cpg$V2 <= 2,]

#remove extra space in seq id
Aele_cpg2 <- as.data.frame(apply(Aele_cpg,2,function(x)gsub('\\s+', '',x)))

#open sequence length data
Aele_seq_length <-read.delim("Aele_seq_length.tsv", sep = "\t", header=F) 

#remove extra space in seq id
Aele_seq_length2 <- as.data.frame(apply(Aele_seq_length,2,function(x)gsub('\\s+', '',x)))

#change seq length data from factor to numeric
Aele_seq_length2[,2] <-as.numeric(as.character(Aele_seq_length2[,2]))

#merge nanopore data with seq length data
Meth_length <-merge(Aele_meth2,Aele_seq_length2,by="V1")

#change CpG data from factor to numeric
Aele_cpg2[,2] <-as.numeric(as.character(Aele_cpg2[,2]))

#merge CpG O/E data with the other data
Ae_merged <- na.omit(merge(Meth_length,Aele_cpg2,by="V1"))

#smooth scatter plot of nanopore meth data vs. CpG O/E data
smoothScatter(log(Ae_merged[,2]), Ae_merged[,4], ylim = c(0,2.5))

#linear regression of nanopore meth data vs. CpG O/E data
model <- lm(log(Ae_merged[,2]) ~ Ae_merged[,4])
abline(model)

#check sequence length distribution
hist(Aele_seq_length[,2], xlim = c(0,5000), breaks = 500)

#set threshold of seq length = 300 for transcriptome seq length
Ae_merged2 <- Ae_merged[Ae_merged$V2.y >= 300,]

smoothScatter(log(Ae_merged2[,2]), Ae_merged2[,4], ylim = c(0,2.5))

#######################################################################

#run same analysis with nanopore data set at called site threshold of 10


setwd("~/Documents/Projects/Anthopleura/Methylation&splicing/Anthopleura_genome")

#open nanopore methylation data
Aele_meth <-read.delim("methylation_avg_thresh10.tsv", sep = "", header=F)

#set min and max thresholds for data
Aele_meth2 <- Aele_meth[Aele_meth$V2 >= 0.001 & Aele_meth$V2 <= 10,]

#change directory 
setwd("~/Documents/Projects/Anthopleura/Methylation&splicing/Aele/")

#open CpG O/E data
Aele_cpg <-read.delim("ID_CpG", header=F)

#set min and max thresholds
Aele_cpg2 <- Aele_cpg[Aele_cpg$V2 >= 0.001 & Aele_cpg$V2 <= 2,]

#remove extra space in seq id
Aele_cpg2 <- as.data.frame(apply(Aele_cpg2,2,function(x)gsub('\\s+', '',x)))

#open sequence length data
Aele_seq_length <-read.delim("Aele_seq_length.tsv", sep = "\t", header=F) 

#remove extra space in seq id
Aele_seq_length2 <- as.data.frame(apply(Aele_seq_length,2,function(x)gsub('\\s+', '',x)))

#change seq length data from factor to numeric
Aele_seq_length2[,2] <-as.numeric(as.character(Aele_seq_length2[,2]))

#merge nanopore data with seq length data
Meth_length <-merge(Aele_meth2,Aele_seq_length2,by="V1")

#change CpG data from factor to numeric
Aele_cpg2[,2] <-as.numeric(as.character(Aele_cpg2[,2]))

#merge CpG O/E data with the other data
Ae_merged <- na.omit(merge(Meth_length,Aele_cpg2,by="V1"))

#smooth scatter plot of nanopore meth data vs. CpG O/E data
smoothScatter(log(Ae_merged[,2]), Ae_merged[,4], ylim = c(0,2.5))

#linear regression of nanopore meth data vs. CpG O/E data
model <- lm(log(Ae_merged[,2]) ~ Ae_merged[,4])
abline(model)

#check sequence length distribution
hist(Aele_seq_length[,2], xlim = c(0,10000), breaks = 500)

#set threshold of seq length = 300 for transcriptome seq length
Ae_merged2 <- Ae_merged[Ae_merged$V2.y >= 300,]

#smooth scatter of thresholded data
smoothScatter(sqrt(Ae_merged2[,2]), Ae_merged2[,4], ylim = c(0,1.5), xlab = "sqrt(Nanopore methylation level)", ylab = "CpG O/E")

#linear regression of nanopore meth data vs. CpG O/E data (with 300 seq length threshold)
model <- lm(sqrt(Ae_merged2[,2]) ~ Ae_merged2[,4])
abline(model)

#srt transform
srt_meth <- sqrt(Ae_merged2[,2])

Ae_merged3 <- cbind(Ae_merged2, srt_meth, model$fitted.values)

# ggpolot

library(ggplot2)
library(ggpubr)

ggscatterhist(Ae_merged3, x = "srt_meth", y = "V2",
          size = 0.2, alpha = 0.1,
          margin.params = list(fill = "lightgray",
          add = "model$fitted.values"))+
          border()

#see http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

# Are there better ways of transforming data? Square root? Log? Arcsin? inverse?
# Are these thresholds on base calls and seq length necessary?
# Are there any other thresholds that we should be using instead?
# Are there better ways of visualizing this data, e.g. side by side density plots? (e.g. https://github.com/daattali/ggExtra/issues/61)
# Does mixture model analysis of nanopore methylation data suggest two distributions as with CpG O/E?
# Should we be thresholding out extra long sequences?


dens <- density(Ae_merged2[,2])
plot(dens)

################################################################

setwd("~/Documents/Projects/Anthopleura/Methylation&splicing/Anthopleura_genome")

#open nanopore methylation data
Aele_meth_freq <-read.delim("methylation_frequency_thresh10.tsv", sep = "\t", header=T)

Aele_meth_avg <- aggregate(methylated_frequency ~ chromosome, Aele_meth_freq, mean)

#change directory 
setwd("~/Documents/Projects/Anthopleura/Methylation&splicing/Aele/")

#open CpG O/E data
Aele_cpg <-read.delim("ID_CpG", header=F)

#set min and max thresholds
Aele_cpg2 <- Aele_cpg[Aele_cpg$V2 >= 0.001 & Aele_cpg$V2 <= 2,]

#remove extra space in seq id
Aele_cpg2 <- as.data.frame(apply(Aele_cpg2,2,function(x)gsub('\\s+', '',x)))

#open sequence length data
Aele_seq_length <-read.delim("Aele_seq_length.tsv", sep = "\t", header=F) 

#remove extra space in seq id
Aele_seq_length2 <- as.data.frame(apply(Aele_seq_length,2,function(x)gsub('\\s+', '',x)))

#change seq length data from factor to numeric
Aele_seq_length2[,2] <-as.numeric(as.character(Aele_seq_length2[,2]))

#merge nanopore data with seq length data
Meth_length <-merge(Aele_meth_avg,Aele_seq_length2,by.x="chromosome", by.y = "V1")

#change CpG data from factor to numeric
Aele_cpg2[,2] <-as.numeric(as.character(Aele_cpg2[,2]))

#merge CpG O/E data with the other data
Ae_merged <- na.omit(merge(Meth_length,Aele_cpg2,by.x="chromosome", by.y = "V1"))

#smooth scatter plot of nanopore meth data vs. CpG O/E data
plot(Ae_merged[,2], Ae_merged[,4])

#linear regression of nanopore meth data vs. CpG O/E data
model <- lm(Ae_merged[,2] ~ Ae_merged[,4])
abline(model, col = "green")

#check sequence length distribution
hist(Aele_seq_length[,2], xlim = c(0,5000), breaks = 500)

#set threshold of seq length = 300 for transcriptome seq length
Ae_merged2 <- Ae_merged[Ae_merged$V2.x >= 300,]


#log+1 transform
log_meth <- log(Ae_merged2[,2]+1)

#bind transformed data to dataframe
Ae_merged3 <- cbind(Ae_merged2, log_meth)

#plot log+1 transform
smoothScatter(Ae_merged3[,5], Ae_merged3[,4], ylim = c(0,2))

#linear regression of log transformed nanopore meth data vs. CpG O/E data
model <- lm(Ae_merged3[,5] ~ Ae_merged3[,4])
abline(model)




# ggpolot

library(ggplot2)
library(ggpubr)

ggscatterhist(Ae_merged3, x = "log_meth", y = "V2.y",
              size = 0.2, alpha = 0.1,
              margin.params = list(fill = "lightgray"))+
  border()

##################################################################################

setwd("~/Documents/Projects/Anthopleura/Methylation&splicing/Anthopleura_genome")

#open nanopore methylation data
Aele_meth_freq <-read.delim("methylation_frequency_thresh10.tsv", sep = "\t", header=T)

#get average meth frequency
Aele_meth_avg <- aggregate(methylated_frequency ~ chromosome, Aele_meth_freq, mean)

#get sum of CpG motifs
Aele_CpG_motifs <- aggregate(num_motifs_in_group ~ chromosome, Aele_meth_freq, sum)

#merge meth frequency and sum motifs
Meth_avg_sum_motifs <-merge(Aele_meth_avg,Aele_CpG_motifs,by="chromosome")

#change directory 
setwd("~/Documents/Projects/Anthopleura/Methylation&splicing/Aele/")

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

#now set threshold of 20% CpG coverage to data
All_data_thresh <- All_data[All_data[,7] >= 0.2,]

#now set threshold of 2 for CpG O/E 
All_data_thresh <- All_data_thresh[All_data_thresh[,6] <= 2,]

#smooth scatter plot of nanopore meth data vs. CpG O/E data
plot(All_data_thresh[,2], All_data_thresh[,6])

#linear regression of nanopore meth data vs. CpG O/E data
model <- lm(All_data_thresh[,2] ~ All_data_thresh[,6])
abline(model, col = "green")


# load ggpolot and ggpubr packages

library(ggplot2)
library(ggpubr)

# Scatter plot with density panels using ggplot and ggpubr
sp <- ggscatter(All_data_thresh, x = "methylated_frequency", y = "CpG_OE",
                color = "#0073C2FF", size = 0.4, alpha = 0.1,
                add = "loess", add.params = list(color = "#EFC000FF"),
                cor.coef = TRUE, 
                cor.coeff.args = list(method = "pearson", 
                                      label.x.npc = "right", label.y.npc = "top"),
                ylab = "CpG O/E", xlab = "Methylated frequency")+
  border()                                         
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(All_data_thresh, "methylated_frequency")
yplot <- ggdensity(All_data_thresh, "CpG_OE")+
  rotate()
# Cleaning the plots
yplot <- yplot + clean_theme() 
xplot <- xplot + clean_theme()
# Arranging the plot
ggarrange(xplot, NULL, sp, yplot, 
          ncol = 2, nrow = 3,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)

#compare fully methylated to fully unmethylated

CpGOE_0 <- as.numeric(All_data_thresh$CpG_OE[All_data_thresh[,2] == 0])
CpGOE_1 <- as.numeric(All_data_thresh$CpG_OE[All_data_thresh[,2] == 1])
CpGOE_0 <- cbind(CpGOE_0, rep("Fully unmethylated", length(CpGOE_0)))
CpGOE_1 <- cbind(CpGOE_1, rep("Fully methylated", length(CpGOE_1)))

CpG_01 <- as.data.frame(rbind(CpGOE_1, CpGOE_0))
CpG_01[,1] <- as.numeric(as.character(CpG_01[,1]))
CpG_01$V2 = factor(CpG_01$V2,c("Fully unmethylated","Fully methylated"))

#boxplot
bxp <- ggboxplot(CpG_01, x = "V2", y = "CpGOE_1",
                 color = "V2", palette = c("#0073C2FF", "#EFC000FF"),
                 xlab = "", ylab = "CpG O/E", show.legend = FALSE)
bxp


# Scatter plot with density panels using ggplot and ggpubr, plus boxplot
sp <- ggscatter(All_data_thresh, x = "methylated_frequency", y = "CpG_OE",
                color = "#868686FF", size = 0.4, alpha = 0.1,
                add = "loess", add.params = list(color = "#CD534CFF"),
                cor.coef = TRUE, 
                cor.coeff.args = list(method = "pearson", 
                                      label.x.npc = "right", label.y.npc = "top"),
                ylab = "CpG O/E", xlab = "Methylated frequency")+
  border()                                         
#Boxplot of fully meth and umeth by CpG O/E
bxp <- ggboxplot(CpG_01, x = "V2", y = "CpGOE_1",
                 color = "V2", palette = c("#0073C2FF", "#EFC000FF"),
                 xlab = "", ylab = "CpG O/E", show.legend = FALSE)
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(All_data_thresh, "methylated_frequency")
yplot <- ggdensity(All_data_thresh, "CpG_OE")+
  rotate()
# Cleaning the plots
yplot <- yplot + clean_theme() 
xplot <- xplot + clean_theme()
# Arranging the plots
ggarrange(xplot, NULL, sp, yplot, bxp, NULL, 
          ncol = 2, nrow = 3,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)

