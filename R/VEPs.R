library(tidyverse)
library(ggstatsplot)
library(scales)
library(ggpubr)
library(pROC)
library(patchwork)
library(predtools)
library(yardstick)
library(PRROC)
library(scales)
library(ggrepel)


path <- "/home/snouto/projects/genesis/second/newestfigures"

theme_Publication <- function(base_size=14, base_family="sans",...) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.direction = "horizontal",
            legend.box = "vetical",
            legend.key.size= unit(0.5, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold"),...
    ))
  
}

nucleotide_level_veps <- read_csv("/home/snouto/projects/simple/nucleotides.level.csv")
veps.types <- read.csv("/home/snouto/projects/genesis/second/veps/veps.types.csv")
bens.latest <- read_csv("/home/snouto/projects/genesis/second/veps.latest.classification.csv")
eval <- read_csv("/home/snouto/projects/anomaly/predictors.auc.csv") %>% 
  inner_join(bens.latest,by="vep")

eval %>% write_csv("/home/snouto/projects/genesis/second/out/predictors.AUC.csv")
## Load Metrics Irrespective of protein structural region
evals <- read_csv("/home/snouto/projects/genesis/second/metrics/metrics.all.csv") %>% 
  filter((vep %in% bens.latest$vep)) %>% filter(!(vep %in% c("MOIpred_recesive","sequence_unet",
                                                            "MutationTaster",
                                                            "BLOSUM62",
                                                            "Grantham",
                                                            "AlphScore",
                                                            nucleotide_level_veps$VEP))) %>%
  left_join(bens.latest,by='vep')

evals %>% write_csv("/home/snouto/projects/genesis/second/out/veps.evaluation.csv")

paper_veps <- data.frame(vep=unique(evals$vep))

write_csv(paper_veps,"/home/snouto/projects/genesis/second/paper.veps.csv")

(p1 <- evals  %>% drop_na(type) %>%
  ggplot(aes(x=specificity,y=recall,label=vep)) + 
  geom_point(aes(colour=type),size=5) + geom_text_repel(max.overlaps = 500,
                                                        size = 5) +
  theme_bw() + 
  labs(x="",y="Sensitivity") +
  guides(size="none",colour=guide_legend(title="Type")) + 
  #scale_colour_manual(values=c("red","green","blue","black")) +
  scale_y_continuous(breaks = seq(0,1,0.1)) + 
  scale_x_continuous(breaks = seq(0,1,0.1))+
  theme(legend.position = "top",
        legend.text = element_text(size=25),
        axis.text = element_text(size=20),
        axis.title = element_text(size=30),
        strip.text = element_text(size=20),
        legend.spacing.x = unit(0.5,"cm")
  ) + 
  guides(colour=guide_legend(title=""))+
  facet_wrap(~region,labeller = labeller(region=c(disordered="Disordered",intermediate="Intermediate",
                                                  ordered="Ordered")))
)

evals2 <- read_csv("/home/snouto/projects/genesis/second/metrics/regions.metrics.csv") %>% 
  filter((vep %in% veps.types$vep)) %>% filter(!(vep %in% c("MOIpred_recesive","sequence_unet",
                                                            "MutationTaster",
                                                            "BLOSUM62",
                                                            "Grantham",
                                                            nucleotide_level_veps$VEP))) %>%
  left_join(bens.latest,by='vep') %>% filter(vep %in% paper_veps$vep)

##### Quantify the gain in sensitivity when using region-specific optimal threshold


tbl1 <- evals %>% group_by(type,region) %>% summarise(
  avg_recall = mean(recall,na.rm = T),
  avg_specificity = mean(specificity,na.rm=T),
  sem_recall = sd(recall,na.rm=T)/sqrt(n()),
  sem_specificity = sd(specificity,na.rm=T) / sqrt(n()),
  recall_low = avg_recall - sem_recall,
  recall_high = avg_recall + sem_recall,
  specificity_low = avg_specificity - sem_specificity,
  specificity_high = avg_specificity + sem_specificity
  
) %>% mutate(condition="Non specific")

tbl2 <- evals2 %>% group_by(type,region) %>% summarise(
  avg_recall = mean(recall,na.rm = T),
  avg_specificity = mean(specificity,na.rm=T),
  sem_recall = sd(recall,na.rm=T)/sqrt(n()),
  sem_specificity = sd(specificity,na.rm=T) / sqrt(n()),
  recall_low = avg_recall - sem_recall,
  recall_high = avg_recall + sem_recall,
  specificity_low = avg_specificity -  sem_specificity,
  specificity_high = avg_specificity + sem_specificity
  
) %>% mutate(condition="Region specific")


evals_df <- rbind(tbl1,tbl2)


evals_df2 <- rbind(
  evals %>% select(region,type,vep,recall) %>% mutate(condition="Non specific"),
  evals2 %>% select(region,type,vep,recall) %>% mutate(condition="Region specific")
) %>% filter(region == "disordered")

gain_sensitivity <- evals_df2 %>% group_by(vep) %>%
  mutate(recall_diff = recall[condition == "Region specific"] - recall[condition == "Non specific"]) %>%
  ggplot(aes(x=reorder(vep,recall_diff),y=recall_diff,fill=type))+
  geom_col(position=position_dodge(width=0.9))+
  theme_pubr()+
  theme(
    axis.text = element_text(size=15,angle=90),
    axis.title = element_text(size=20),
    legend.text = element_text(size=20)
  )+
  labs(x="",y="Difference in sensitivity",fill="")
  

ggsave(
  paste0(path,"/","S3.png"),
  plot=gain_sensitivity,
  width=21.22,
  height = 10.573,
  units = "in", dpi = 300
)

################################################################################
#Can you make a new version of Figure S3, but expanded to 6 panels.
#Left side – sensitivity difference (as is); right side – specificity difference
#Top – disordered, Middle – intermediate; Bottom – ordered
S3 <- rbind(
  evals %>% select(region,type,vep,recall,specificity) %>% mutate(condition="Non specific"),
  evals2 %>% select(region,type,vep,recall,specificity) %>% mutate(condition="Region specific")
)

(TPR_disorder <- S3 %>% filter(region == "disordered") %>% group_by(vep) %>%
  mutate(recall_diff = recall[condition == "Region specific"] - recall[condition == "Non specific"]) %>%
  ggplot(aes(x=reorder(vep,recall_diff),y=recall_diff,fill=type))+
  scale_y_continuous(breaks=seq(-0.20,0.20,0.1),limits=c(-0.20,0.20))+
  geom_col(position=position_dodge(width=0.9))+
  theme_pubr()+
  theme(
    axis.text = element_text(size=25,angle=90,face = "bold"),
    axis.title = element_text(size=30),
    legend.text = element_text(size=20)
  )+
  labs(x="",y="Difference in sensitivity",fill="")
)

ggsave(
  paste0(path,"/","TPR_disorder.png"),
  plot=TPR_disorder,
  width=21.22,
  height = 10.573,
  units = "in", dpi = 300
)

(SP_disorder <- S3 %>% filter(region == "disordered") %>% group_by(vep) %>%
    mutate(SP_diff = specificity[condition == "Region specific"] - specificity[condition == "Non specific"]) %>%
    ggplot(aes(x=reorder(vep,SP_diff),y=SP_diff,fill=type))+
    geom_col(position=position_dodge(width=0.9))+
    scale_y_continuous(breaks=seq(-0.20,0.20,0.1),limits=c(-0.20,0.20))+
    theme_pubr()+
    theme(
      axis.text = element_text(size=25,angle=90,face="bold"),
      axis.title = element_text(size=30),
      legend.text = element_text(size=20)
    )+
    labs(x="",y="Difference in specificity",fill="")
)

ggsave(
  paste0(path,"/","SP_disorder.png"),
  plot=SP_disorder,
  width=21.22,
  height = 10.573,
  units = "in", dpi = 300
)

# Intermediate
(TPR_intermediate <- S3 %>% filter(region == "intermediate") %>% group_by(vep) %>%
    mutate(recall_diff = recall[condition == "Region specific"] - recall[condition == "Non specific"]) %>%
    ggplot(aes(x=reorder(vep,recall_diff),y=recall_diff,fill=type))+
    geom_col(position=position_dodge(width=0.9))+
    scale_y_continuous(breaks=seq(-0.20,0.20,0.1),limits=c(-0.20,0.20))+
    theme_pubr()+
    theme(
      axis.text = element_text(size=25,angle=90,face="bold"),
      axis.title = element_text(size=30),
      legend.text = element_text(size=20)
    )+
    labs(x="",y="Difference in sensitivity",fill="")
)

ggsave(
  paste0(path,"/","TPR_intermediate.png"),
  plot=TPR_intermediate,
  width=21.22,
  height = 10.573,
  units = "in", dpi = 300
)

(SP_intermediate <- S3 %>% filter(region == "intermediate") %>% group_by(vep) %>%
    mutate(SP_diff = specificity[condition == "Region specific"] - specificity[condition == "Non specific"]) %>%
    ggplot(aes(x=reorder(vep,SP_diff),y=SP_diff,fill=type))+
    geom_col(position=position_dodge(width=0.9))+
    scale_y_continuous(breaks=seq(-0.20,0.20,0.1),limits=c(-0.20,0.20))+
    theme_pubr()+
    theme(
      axis.text = element_text(size=25,angle=90,face="bold"),
      axis.title = element_text(size=30),
      legend.text = element_text(size=20)
    )+
    labs(x="",y="Difference in specificity",fill="")
)
ggsave(
  paste0(path,"/","SP_intermediate.png"),
  plot=SP_intermediate,
  width=21.22,
  height = 10.573,
  units = "in", dpi = 300
)

# Ordered
(TPR_ordered <- S3 %>% filter(region == "ordered") %>% group_by(vep) %>%
    mutate(recall_diff = recall[condition == "Region specific"] - recall[condition == "Non specific"]) %>%
    ggplot(aes(x=reorder(vep,recall_diff),y=recall_diff,fill=type))+
    geom_col(position=position_dodge(width=0.9))+
    scale_y_continuous(breaks=seq(-0.20,0.20,0.1),limits=c(-0.20,0.20))+
    theme_pubr()+
    theme(
      axis.text = element_text(size=25,angle=90,face="bold"),
      axis.title = element_text(size=30),
      legend.text = element_text(size=20)
    )+
    labs(x="",y="Difference in sensitivity",fill="")
)
ggsave(
  paste0(path,"/","TPR_ordered.png"),
  plot=TPR_ordered,
  width=21.22,
  height = 10.573,
  units = "in", dpi = 300
)

(SP_ordered <- S3 %>% filter(region == "ordered") %>% group_by(vep) %>%
    mutate(SP_diff = specificity[condition == "Region specific"] - specificity[condition == "Non specific"]) %>%
    ggplot(aes(x=reorder(vep,SP_diff),y=SP_diff,fill=type))+
    geom_col(position=position_dodge(width=0.9))+
    scale_y_continuous(breaks=seq(-0.20,0.20,0.1),limits=c(-0.20,0.20))+
    theme_pubr()+
    theme(
      axis.text = element_text(size=25,angle=90,face="bold"),
      axis.title = element_text(size=30),
      legend.text = element_text(size=20)
    )+
    labs(x="",y="Difference in specificity",fill="")
)
ggsave(
  paste0(path,"/","SP_ordered.png"),
  plot=SP_ordered,
  width=21.22,
  height = 10.573,
  units = "in", dpi = 300
)

ggarrange(TPR_disorder,SP_disorder,TPR_intermediate,SP_intermediate,
         TPR_ordered,SP_ordered,ncol=2,nrow=3,legend="top",common.legend=TRUE,
         labels=c("Disordered","Disordered","Intermediate","Intermediate",
                  "Ordered","Ordered"))

################################################################################

##### Can you make boxplots of different VEP types across different structural regions
nonspecific <- evals %>% select(region,vep,recall,specificity,f1_score,type) %>% mutate(
  condition = "Non specific"
)
specific <- evals2 %>% select(region,vep,recall,specificity,f1_score,type) %>% mutate(
  condition = "Region specific"
)

evals_df <- rbind(nonspecific,specific)

(p2 <- evals2  %>% drop_na(type) %>%
    ggplot(aes(x=specificity,y=recall,label=vep)) + 
    geom_point(aes(colour=type),size=5) + geom_text_repel(max.overlaps = 500,
                                                          size=5) +
    theme_bw() + 
    labs(x="Specificity",y="Sensitivity") +
    guides(size="none",colour=guide_legend(title="Type")) + 
    #scale_colour_manual(values=c("red","green","blue","black")) +
    scale_y_continuous(breaks = seq(0,1,0.1)) + 
    scale_x_continuous(breaks = seq(0,1,0.1))+
    theme(legend.position = "top",
          legend.text = element_text(size=25),
          axis.text = element_text(size=20),
          axis.title = element_text(size=30),
          strip.text = element_text(size=20),
          legend.spacing.x = unit(0.5,"cm")
    ) + 
    guides(colour=guide_legend(title=""))+
    facet_wrap(~region,labeller = labeller(region=c(disordered="Disordered",intermediate="Intermediate",
                                                    ordered="Ordered")))
)

(vepsfig <- ggarrange(p1,p2,common.legend = TRUE,nrow=2,ncol=1))

## Figure 5
ggsave(
  paste0(path,"/","f4.png"),
  plot=vepsfig,
  width=25.22,
  height = 20.573,
  units = "in", dpi = 300
)

## Specificity versus Recall
evals  %>% drop_na(type) %>%
  ggplot(aes(x=recall,y=specificity,label=vep)) + 
  geom_point(aes(colour=type),size=5) + geom_text_repel(max.overlaps = 500,label.size=0.05,label.r=0.15) +
  theme_bw() + 
  labs(x="Recall",y="Specificity") +
  guides(size="none",colour=guide_legend(title="Type")) + 
  #scale_colour_manual(values=c("red","green","blue","black")) +
  scale_y_continuous(breaks = seq(0,1,0.1)) + 
  scale_x_continuous(breaks = seq(0,1,0.1))+
  theme(legend.position = "top",
        legend.text = element_text(size=15),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20),
        strip.text = element_text(size=15)
  ) + 
  guides(colour=guide_legend(title=""))+
  facet_wrap(~region,labeller = labeller(region=c(disordered="Disordered",intermediate="Intermediate",
                                                  ordered="Ordered")))

######################################################################################################
#-What if some or all of the clinical trained VEPs are 
#using more appropriate properties that are more predictive in the context of disorder 
#e.g. binding partners, 3D structure and the surrounding residue context.

# Can you check this in disordered regions

veps <- read_excel("/home/snouto/projects/genesis/second/veps.info.xlsx")

veps.classification <- veps %>% 
  mutate(training_type = ifelse(str_detect(str_to_lower(Features),"structur"),"Structure","Other")) %>%
  select(VEP_name,training_type) %>% rename(vep=VEP_name)

evals.augmented <- evals %>% left_join(veps.classification,by="vep")

evals.augmented %>% drop_na(type) %>%
  ggplot(aes(x=1 - specificity,y=recall,label=vep)) + 
  geom_point(aes(colour=training_type),size=5) + geom_text_repel(max.overlaps = 500) +
  theme_bw() + 
  labs(x="False positive rate",y="") +
  guides(size="none",colour=guide_legend(title="Type")) + 
  #scale_colour_manual(values=c("red","green","blue","black")) +
  scale_y_continuous(breaks = seq(0,1,0.1)) + 
  scale_x_continuous(breaks = seq(0,1,0.1))+
  theme(legend.position = "top",
        legend.text = element_text(size=15),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20),
        strip.text = element_text(size=15)
  ) + 
  guides(colour=guide_legend(title=""))+
  facet_wrap(~region,labeller = labeller(region=c(disordered="Disordered",intermediate="Intermediate",
                                                  ordered="Ordered")))

(ea <- evals.augmented %>% group_by(region,training_type,type) %>%
  summarise(
    avg_recall = mean(recall,na.rm = T),
    sem_recall = sd(recall,na.rm=T) / sqrt(n()),
    recall_low = avg_recall - sem_recall,
    recall_high = avg_recall + sem_recall
  ) %>% filter(type != "Population-tuned") %>% ggplot(aes(x=type,y=avg_recall,fill=training_type)) +
  geom_col(position=position_dodge(width=0.9))+
  geom_errorbar(aes(ymin=recall_low,ymax=recall_high),width=0.2,color='black',
                position=position_dodge(width=0.9))+
  scale_fill_manual(values = c("#690B22","#E07A5F"))+
  theme_pubr()+
  theme(
    axis.text = element_text(size=20),
    axis.title = element_text(size=20),
    legend.text = element_text(size=20),
    strip.text = element_text(size=20)
  )+
  labs(x="",y="Average Sensitivity",fill="")+
  guides(colour=guide_legend(title=""))+
  facet_wrap(~region,labeller = labeller(region=c(disordered="Disordered",intermediate="Intermediate",
                                                  ordered="Ordered")))
)


## Figure 5
ggsave(
  paste0(path,"/","S4.png"),
  plot=ea,
  width=25.22,
  height = 20.573,
  units = "in", dpi = 300
)

######################################################################################################

#### Saving optimal thresholds for publication
ot <- read_csv("/home/snouto/projects/genesis/second/newestfigures/metrics/veps.optimal.thresholds.csv")
ot.regions <- read_csv("/home/snouto/projects/genesis/second/newestfigures/metrics/optimal.thresholds.regions.csv")


ot.final <- ot %>% filter(vep %in% paper_veps$vep)
ot.regions.final <- ot.regions %>% filter(vep %in% paper_veps$vep)

ot.final %>% write_csv("/home/snouto/projects/genesis/second/newestfigures/metrics/optimal.thresholds.final.csv",
                       )

ot.regions.final %>% write_csv("/home/snouto/projects/genesis/second/newestfigures/metrics/optimal.thresholds.per.regions.csv")
