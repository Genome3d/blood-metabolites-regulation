setwd("/Users/tayazafadason/Downloads")
met.data <- read.delim("significant_eqtls.txt", header = TRUE)
met.data$SNP_Gene <- paste0(met.data$SNP, '_', met.data$Gene_Name)
library(ggplot2)
library(plyr)
s <- unique(met.data[,c(1,4)])
s <- count(s$SNP)
ggplot(s, aes(s$x, s$freq))+
  geom_point()+
  coord_flip()
  theme(axis.text.x = element_text(angle=90))
  
# Plot eQTL effect size
met.data <- met.data[order(abs(met.data$Effect_Size), decreasing = T), ]

ggplot(met.data, 
       aes(x=factor(SNP_Gene, levels=unique(SNP_Gene)), # Order by effect size
           y = Effect_Size))+
  geom_point(aes(color=cis_SNP.gene_interaction), alpha=0.5, size=1)+
  theme_classic()+
  scale_y_continuous(limits = c(-2,2), breaks=seq(-2, 2, 0.5),
                     expand = c(0,0))+
  scale_x_discrete(expand=c(0.02,0))+
  scale_color_manual(values=c("red", "grey80"), labels=c("Trans", "Cis"))+
  labs(y="eQTL efect size")+
  guides(color=guide_legend(NULL))+
  theme(axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())+
  geom_hline(yintercept = 0)
