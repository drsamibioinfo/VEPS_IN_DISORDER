library(tidyverse)
library(ggpubr)
library(scales)
library(patchwork)
library(ggsignif)
library(rstatix)
library(ggstatsplot)
library(readxl)
library(ggrepel)
library(envalysis)
library(reshape2)
####################################################### Functions ###############################################################################
path <- "/home/snouto/projects/genesis/second/newestfigures"
### function to calculate the log odds and confidence interval
calc_odds_ci <- function(groups_tbl,subjects = NULL ,grp_name="groups"){
  groups <- groups_tbl %>% as_tibble(.) %>% distinct(!!sym(grp_name))
  grp_names <- c()
  odds <- c()
  lower_ci <- c()
  upper_ci <- c()
  lors <- c()
  lower_ci_l <- c()
  upper_ci_l <- c()
  first_grp <- subjects[1]
  second_grp <- subjects[2]
  for (group in groups){
    first_freq <- groups_tbl[first_grp,group]
    second_freq <- groups_tbl[second_grp,group]
    first_total <- sum(groups_tbl[first_grp,])
    second_total <- sum(groups_tbl[second_grp,])
    # Odds
    odds_first <- (first_freq / first_total) / (1 - (first_freq / first_total))
    odds_second <- (second_freq / second_total) / (1 - (second_freq / second_total))
    or <-  odds_first / odds_second
    lor <- log(or)
    lci <- lor - (1.96 * sqrt((1/first_freq) + (1/second_freq)))
    uci <-  lor + (1.96 * sqrt((1/first_freq) + (1/second_freq)))
    
    
    grp_names <- c(grp_names,group)
    odds <- c(odds,or)
    lower_ci <- c(lower_ci,exp(lci))
    upper_ci <- c(upper_ci,exp(uci))
    lower_ci_l <- c(lower_ci_l,lci)
    upper_ci_l <- c(upper_ci_l,uci)
    lors <- c(lors,lor)
  }
  df <- data.frame(groups=grp_names,odds=odds,lower_ci=lower_ci,upper_ci = upper_ci,lors=lors,lower_ci_l=lower_ci_l,
                   upper_ci_l = upper_ci_l)
  return(df)
}

# Define a custom labelling function
custom_labeller <- function(variable,value) {
  # Define the titles for each panel based on the category values
  if (value == "LOF") {
    return("Loss of function")
  } else if (value == "GOF") {
    return("Gain of function")
  } else {
    return("Dominant negative")
  }
}

location_labeller <- function(variable){
  return (str_to_title(variable))
}

get_binom_conf <- function(conf_matrix,category,target){
  df <- conf_matrix %>% as_tibble(.)
  cat_u <- unlist(unique(df[category])[category],use.names = FALSE)
  target_u <- unlist(unique(df[target])[target],use.names = FALSE)
  low <- c()
  high <- c()
  sign <- c()
  cats <- c()
  tars <- c()
  for (cat in cat_u){
    for(t in target_u){
      test <- prop_test(
        x = conf_matrix[cat,t],
        n = sum(conf_matrix[,t]),
        correct = FALSE,
        detailed = TRUE
      )
      cats <- c(cats,cat)
      tars <- c(tars,t)
      low <- c(low,test$conf.low)
      high <- c(high,test$conf.high)
      sign <- c(sign,test$p.signif)
      
    }
  }
  df <- setNames(data.frame(cats, tars, low, high, sign), c(category, target, "low", "high", "sign"))
  return (df)
}

#Add a panel showing how the classes are defined, with cartoon representations of disorder, intermediate and ordered
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
# This function will calculate the 95% confidence interval for the differences in proportions between two groups
num_supscript <- function(w, digits=0) {
  sprintf(paste0("%.", digits, "f*x*10^%d"), w/10^floor(log10(abs(w))), floor(log10(abs(w))))
}


#############################################################################################################

## Final paper plots

#############################################################################################################
# Pick Collagen-containing proteins to exclude them from the analysis
comparisons <- list(c("disordered","intermediate","ordered"))

# Extracting all proteins containing the collagen pFAM domain
pfam <- read_tsv("/home/snouto/projects/secondpaper/human_pfam.tsv") %>% filter(domain == "PF01391")

excluded_collagen_proteins <- unique(pfam$uniprot_id)
general <- read_csv("/home/snouto/projects/genesis/second/general.stats.withoutcollagen.csv")
#general %>% View()
total_proteome <- as.numeric(general[general$item == "total_proteome","value"])
disordered <- as.numeric(general[general$item == "disordered","value"])
ordered <- as.numeric(general[general$item == "ordered","value"])
intermediate <- as.numeric(general[general$item == "intermediate","value"])
pcts <- c(disordered/total_proteome,intermediate/total_proteome,ordered/total_proteome)
pcts.names <- c("Disordered","Intermediate","Ordered")
intervals.disordered <- prop_test(x=disordered,n=total_proteome,correct = FALSE,detailed = TRUE)
intervals.intermediate <- prop_test(x=intermediate,n=total_proteome,correct = FALSE,detailed=TRUE)
intervals.ordered <- prop_test(x=ordered,n=total_proteome,correct=FALSE,detailed = TRUE)
general.df <- data.frame(region=pcts.names,pct = pcts,
                         low = c(intervals.disordered$conf.low,intervals.intermediate$conf.low,
                                 intervals.ordered$conf.low),
                         high = c(intervals.disordered$conf.high,intervals.intermediate$conf.high,
                                  intervals.ordered$conf.high),
                         sign = c(
                           intervals.disordered$p.signif,intervals.intermediate$p.signif,intervals.ordered$p.signif
                         )
)
variants <- read_csv("/home/snouto/projects/genesis/second/cdata.missense.gz",lazy=T
                     ,col_select = c("uniprot_id","gene","ref","alt","pos","significance","location","wt_fcr","wt_ncpr",
                                     "source")) %>% 
  filter(pos != 1) %>% filter(!(uniprot_id %in% excluded_collagen_proteins)) %>% 
  distinct(uniprot_id,ref,pos,alt,significance,.keep_all = T)

variants <- variants %>% mutate(label=ifelse(significance == 0,"Putatively benign","Pathogenic"))

section1 <- variants %>% drop_na(location)

## Figure 1

whole.idr <- read_csv("/home/snouto/projects/genesis/second/wholeidr.stats.csv")

whole.idr <- whole.idr %>% mutate(opct = 1 - (pct + ipct))

# whole.idr %>% mutate(error = sqrt((pct - 0.3)^2 + (ipct - 0.2)^2)) %>% 
#   filter(str_starts(protein,"P")) %>% arrange(desc(size),error) %>% View()



(p0 <- general.df %>%
    mutate(pct = pct * 100,
           region = factor(region, levels = c("Disordered", "Intermediate", "Ordered"))
    ) %>%
    ggplot(aes(x = "", y = pct, fill = reorder(region, region))) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 20,hjust = 0.5),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 15, hjust = 0.5),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()
    ) +
    #scale_fill_brewer(palette = "Set3") +
    scale_fill_manual(values=c("red","purple","darkblue"))+
    labs(
      title = "Human Proteome",
      subtitle = "",
      fill = ""
    ) +
    geom_text(aes(label = paste0(round(pct, 1), "%")),
              position = position_stack(vjust = 0.5),
              size = 9, color = "white") +
    guides(fill = guide_legend(reverse = FALSE))
)


section1 <- section1 %>% mutate(
  source = case_when(
    source == "Clinvar_Benign" ~ "ClinVar benign",
    source == "Clinvar_Pathogenic" ~ "ClinVar pathogenic",
    source == "gnomAD_Benign" ~ "gnomAD putatively benign"
  )
)
# this script continues "all.in.all.R" script
tbl.counts <- with(section1,table(location,source))

tbl.counts.props <- prop.table(tbl.counts,margin = 2)

(p2 <- tbl.counts.props %>% as_tibble(.) %>% ggplot(aes(x = location, y = n * 100, fill = source)) +
    geom_col(position = position_dodge(width=0.9)) +
    geom_label(aes(label = paste0(round(n * 100, 2), "%")), position=position_dodge(width=0.9),
               size=4,show.legend = FALSE,
               colour = "white",fontface="bold") +
    # geom_errorbar(aes(ymin=low*100,ymax=high*100),
    #               position = position_dodge(width=0.9),
    #               color="black",width=0.3)+
    labs(x = "Structural Location(s)", y = "% of variants") +
    scale_y_continuous(breaks = seq(0, 100, 10)) +
    theme_publish() +
    theme(legend.position = "bottom", 
          axis.text = element_text(size = 20),
          axis.title = element_text(size=15),
          legend.text = element_text(size=20,margin=margin(r=8)),
          legend.title = element_text(size=20)) +
    scale_fill_manual(values=c("red","purple","darkblue"))+
    guides(fill=guide_legend(title=""))
)

##
## Panel C: I think it might look better grouping by variant type rather than structural region. 
# So have ClinVar pathogenic, and show bars for disordered, intermediate and ordered 
# (using the same colouring as for the other panels), 
# then ClinVar benign, then gnomAD putatively benign.
##
ordering <- data.frame(source=c("ClinVar benign","ClinVar pathogenic",
                                "gnomAD putatively benign"),ordering=c(2,1,3))

section1 <- section1 %>% left_join(ordering,by='source')

panelC_tbl <- with(section1,table(source,location))

panelC.conf <- as.data.frame(as.table(panelC_tbl)) |> 
  mutate(total = sum(Freq), .by = source) |> 
  mutate(conf_int = map2(Freq, total, ~ binom.test(.x, .y)$conf.int),
         p_value = map2(Freq, total, ~ binom.test(.x, .y)$p.value),
         lower = map_dbl(conf_int, 1),
         upper = map_dbl(conf_int, 2)) |>
  select(-conf_int)

(p2 <- prop.table(panelC_tbl,margin=1) %>% as_tibble(.) %>% 
    left_join(ordering,by='source') %>% mutate(n=n*100) %>% left_join(panelC.conf,by=c("source","location")) %>%
    mutate(lower = lower * 100 , upper = upper * 100) %>%
    ggplot(aes(x=reorder(source,ordering),y=n,fill=location))+
    geom_col(position=position_dodge(width=0.9))+
    # geom_label(aes(label=paste0(round(n,2)," %")),
    #            size=4,show.legend = FALSE,
    #            colour = "white",fontface="bold",
    #            position=position_dodge(width=0.9))+
    geom_errorbar(aes(ymin=lower,ymax=upper),position=position_dodge(width=0.9),width=0.2,color='black')+
    geom_text(aes(label=comma(Freq)),vjust=-1.4,
              position=position_dodge(width=0.9),
              size=6,
              fontface = "bold")+
    # geom_text(aes(label = paste0(round(n,2), " %")), 
    #            position=position_dodge(width=0.9),
    #            show.legend = FALSE,
    #            color="white",
    #            fontface="bold",
    #            vjust = 1.3,
    #            size=7
    #            )+
    labs(x = "", y = "% of variants") +
    scale_y_continuous(breaks = seq(0, 100, 10)) +
    theme_publish() +
    theme(legend.position = "bottom", 
          axis.text = element_text(size = 20),
          axis.title = element_text(size=15),
          legend.text = element_text(size=20,margin=margin(r=8)),
          legend.title = element_text(size=20)) +
    scale_fill_manual(values=c("red","purple","darkblue"))+
    guides(fill=guide_legend(title=""))
)

(fig1 <- ggarrange(p0,p2,widths = c(1,2),
                   common.legend = TRUE,legend = "bottom")+
    theme(legend.title = element_text(margin = margin(r = 10),size=10)))



ggsave(
  paste0(path,"/","f1bc.png"),
  plot=fig1,
  width=21.22,
  height = 10.573,
  units = "in", dpi = 300
)

##### this section will print the hexcode for the purple color used in these plots
color <- "purple"

# Get RGB values for the color
rgb_value <- col2rgb(color)

# Convert RGB values to hexadecimal code
hex_code <- rgb(rgb_value[1], rgb_value[2], rgb_value[3], maxColorValue = 255)

# Print the hex code
print(hex_code)
####################################################################################

#########################################################################################################################
# Supplementary Figure
stats <- read_csv("/home/snouto/projects/genesis/second/wholeidr.stats.csv") %>% mutate(
  disorder = as.numeric(size * pct))

bins <- seq(0, 1.1, by = 0.1)
bin_labels <- as.character(bins[-length(bins)])
stats$bin <- cut(stats$pct, breaks = bins, labels = bin_labels, right = FALSE,include.lowest = TRUE)
percentage_data <- data.frame(table(stats$bin) / nrow(stats) * 100)
(percentage_data <- percentage_data %>% mutate(
  x = as.numeric(levels(Var1))* 100
)
)
dummy_labels <- data.frame(
  x = percentage_data$x,
  lbl = c("0","0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100")
)

percentage_data <- percentage_data %>% left_join(dummy_labels,by="x")
percentage_data$lbl <- as.character(percentage_data$lbl)

# Plot a bar chart using the calculated percentages
(S1 <- ggplot(percentage_data, aes(x = lbl, y = Freq)) +
  geom_bar(stat = "identity", fill = "#69b3a2") +
  geom_label(aes(label=round(Freq,2)),size=8)+
  labs(
       y = "Percentage of proteins (%)",
       x = "Percentage of disordered residues per protein") +
  theme_pubr()+
    theme(
      axis.text = element_text(size=20),
      axis.title = element_text(size=20)
    )
)

ggsave(
  paste0(path,"/","S1.png"),
  plot=S1,
  width=21.22,
  height = 10.573,
  units = "in", dpi = 300
)

ggplot(stats, aes(x = pct*100)) +
  geom_histogram(binwidth = 10, color = "black", fill = "#69b3a2") +
  scale_x_continuous(breaks = seq(0,100,5))+
  labs(title = "Distribution of Disorder Percentage in Proteins",
       x = "Percentage of Disorder",
       y = "Frequency") +
  theme_minimal()
##########################################################################################################################


#### Second figure 
inheritance <- read_tsv("/home/snouto/projects/secondpaper/data/omim.latest.2023.tsv") %>% select(-gene) %>% 
  filter(inheritance %in% c("Autosomal dominant","Autosomal recessive"))

section2 <- variants %>% left_join(inheritance,by="uniprot_id",relationship = "many-to-many") %>% 
  drop_na(location,inheritance) %>% filter(label == "Pathogenic")

inheritance_tbl <- with(section2,table(inheritance,location))

inheritance.confint <- as.data.frame(as.table(inheritance_tbl)) |> 
  mutate(total = sum(Freq), .by = inheritance) |> 
  mutate(conf_int = map2(Freq, total, ~ binom.test(.x, .y)$conf.int),
         p_value = map2(Freq, total, ~ binom.test(.x, .y)$p.value),
         lower = map_dbl(conf_int, 1),
         upper = map_dbl(conf_int, 2)) |>
  select(-conf_int)

inheritance.confint$p_value <- as.numeric(inheritance.confint$p_value)

fisher_test(inheritance_tbl,workspace=2^20)

####################
## Doing fisher test again
fisher_data <- section2 %>% mutate(fisher_location=ifelse(location == "disordered","disordered","non-disordered"))
fisher_tbl <- with(fisher_data,table(inheritance,fisher_location))

(fisher_test(fisher_tbl,workspace = 2^20))
##############################
fisher_intermediate <- section2 %>% mutate(fisher_location=ifelse(location == "intermediate","intermediate","non-intermediate"))
fisher_inttbl <- with(fisher_intermediate,table(inheritance,fisher_location))
fisher_test(fisher_inttbl,workspace=2^20)
### Generating 95% confidence interval

(inh1 <- prop.table(inheritance_tbl, margin = 1) %>%
    as_tibble(.) %>% rename(pct="n") %>% left_join(inheritance.confint,by = c("inheritance","location")) %>%
    ggplot(aes(x = location, y = pct*100,fill=inheritance)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = Freq), show.legend = FALSE ,
              position = position_dodge(width = 1),
              vjust=-1.5,
              color = "black", fontface = "bold", size = 6) +
    #geom_text(aes(label = paste0("P value = ",sprintf("%.2e",p_value))),parse = FALSE,x=1.5,y=72,size=6)+ 
    geom_errorbar(aes(ymin=lower*100,ymax=upper*100),width=0.2,color="black",
                  position=position_dodge(width=0.9))+
    scale_y_continuous(breaks=seq(0,100,10))+
    theme_pubclean() +
    theme(
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 20),
      strip.text = element_text(size = 20),
      legend.position = "top",
      legend.text = element_text(size = 20)
    ) +
    scale_fill_manual(values = c("#FF6500", "#1E3E62")) +
    labs(x = "", y = "% ") +
    guides(fill = guide_legend(title = ""))
)


#################################################################################################################

mihaly <- read_csv("/home/snouto/projects/genesis/second/plots_comments/mihaly.final.classification.csv") %>%
  mutate(
    class = str_to_upper(class)
  )

mihaly$gene <- NULL

#variants %>% filter(label == "Pathogenic")
section3 <-  section2 %>%
  filter( inheritance == "Autosomal dominant") %>% left_join(mihaly,by='uniprot_id')

section3 <- section3 %>% distinct(uniprot_id,ref,alt,pos,.keep_all = T)

mm.tbl <- with(section3,table(class,location))

mm.confint <- as.data.frame(as.table(mm.tbl)) |> 
  mutate(total = sum(Freq), .by = class) |> 
  mutate(conf_int = map2(Freq, total, ~ binom.test(.x, .y)$conf.int),
         p_value = map2(Freq, total, ~ binom.test(.x, .y)$p.value),
         lower = map_dbl(conf_int, 1),
         upper = map_dbl(conf_int, 2)) |>
  select(-conf_int)

mm.final <- prop.table(mm.tbl,margin=1) %>% as_tibble() %>% rename(pct="n") %>%
  left_join(mm.confint,by=c("class","location"))

(mmfig <- 
    mm.final %>%
    ggplot(aes(x = location, y = pct*100,fill=class)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = Freq), show.legend = FALSE ,
              position = position_dodge(width = 1),
              vjust=-1.5,
              color = "black", fontface = "bold", size = 6) +
    #geom_text(aes(label = paste0("P value = ",sprintf("%.2e",p_value))),parse = FALSE,x=2,y=72,size=6)+ 
    geom_errorbar(aes(ymin=lower*100,ymax=upper*100),width=0.2,color="black",position=position_dodge(width=0.9))+
    theme_pubclean() +
    theme(
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 20),
      strip.text = element_text(size = 20),
      legend.position = "top",
      legend.text = element_text(size = 20)
    ) +
    scale_y_continuous(breaks=seq(0,100,10))+
    scale_fill_manual(values = c("#4C4B16", "#E6C767","#F87A53")) +
    labs(x = "", y = "% ") +
    guides(fill = guide_legend(title = ""))
  
)

### Calculate fisher exact test between DN and "GOF/LOF" together
(mm2.tbl <- with(section3 %>% mutate(class = case_when(
  class %in% c("LOF","GOF") ~ "Other",
  TRUE ~ class
)),table(class,location)))

mmt <- fisher_test(mm2.tbl,workspace = 2^20)
print(mmt$p)

########################################################################################################################
# Variant Effect Predictors ##################################################
nucleotide_level_veps <- read_csv("/home/snouto/projects/simple/nucleotides.level.csv")
veps.types <- read.csv("/home/snouto/projects/genesis/second/veps/veps.types.csv")
bens.latest <- read_csv("/home/snouto/projects/genesis/second/veps.latest.classification.csv")
eval <- read_csv("/home/snouto/projects/anomaly/predictors.auc.csv") %>% 
  inner_join(bens.latest,by="vep")


## Load Metrics Irrespective of protein structural region
evals <- read_csv("/home/snouto/projects/genesis/second/metrics/metrics.all.csv") %>% 
  filter((vep %in% bens.latest$vep)) %>% filter(!(vep %in% c("MOIpred_recesive","sequence_unet",
                                                             "MutationTaster",
                                                             "BLOSUM62",
                                                             "Grantham",
                                                             "AlphScore",
                                                             nucleotide_level_veps$VEP))) %>%
  left_join(bens.latest,by='vep')

paper_veps <- data.frame(vep=unique(evals$vep))

## Figure 3
nucleotide_level_veps <- read_csv("/home/snouto/projects/simple/nucleotides.level.csv")
bens.latest <- read_csv("/home/snouto/projects/genesis/second/veps.latest.classification.csv")
eval <- read_csv("/home/snouto/Dropbox (PhD)/Genesis/second/metrics2/metrics.all.csv") %>% 
  filter(!(vep %in% nucleotide_level_veps$VEP))

eval <- eval %>% inner_join(bens.latest,by='vep')
eval <- eval %>% filter(!(vep %in% c("MutationTaster")))

eval <- eval %>% filter(vep %in% paper_veps$vep)

regions_aucs <- read_csv("/home/snouto/Dropbox (PhD)/Genesis/second/metrics2/regions.metrics.csv",
                         col_select = c("region","vep","auc"))

eval <- eval %>% left_join(regions_aucs,by=c("region","vep"))

#Update to the three group classifications of Ben's updated preprint
#https://www.biorxiv.org/content/10.1101/2024.05.12.593741v2. Maybe it's ok to include population-tuned and 
#population-free in the same panel, as there aren't many population-tuned in your analysis 
#(just AlphaMissense and LIST-S2 I think)
#-Make sure the colour of the disorder, intermediate, ordered groups matches across all the figures
veps.ranks <- eval %>% filter(region == "disordered") %>%
  mutate(
    auc_rank = rank(auc)
  ) %>% select(vep,auc_rank)
eval$auc_rank <- NULL
eval <- eval %>% left_join(veps.ranks,by="vep")
(poptuned <- eval %>% filter(type %in% c("Population-tuned","Population-free")) %>%
    mutate(auc_adjusted = auc - 0.5) %>%
    arrange(desc(auc_rank)) %>%
    ggplot(aes(x=reorder(vep,auc_rank),y=auc_adjusted,fill=region))+
    geom_col(position=position_dodge(width=0.9))+
    scale_y_continuous(labels = function(x) x + 0.5) +
    scale_fill_manual(values=c("red","purple","darkblue"))+
    theme_pubclean()+
    theme(
      legend.text = element_text(size=20),
      axis.title = element_text(size=20),
      axis.text  = element_text(size=16),
      plot.title = element_text(size=20),
      legend.spacing.x = unit(0.9,"cm")
    )+
    labs(x="Variant effect predictors",y="AUROC",title="Population-tuned and population-free VEPs",fill="")+
    coord_flip()
  
)

(popclinical <- eval %>% filter(type == "Clinical-trained") %>%
    mutate(auc_adjusted = auc - 0.5) %>%
    ggplot(aes(x=reorder(vep,auc_rank),y=auc_adjusted,fill=region))+
    geom_col(position=position_dodge(width=0.9))+
    scale_y_continuous(labels = function(x) x + 0.5) +
    scale_fill_manual(values=c("red","purple","darkblue"))+
    theme_pubclean()+
    theme(
      legend.text = element_text(size=20),
      axis.title = element_text(size=20),
      axis.text  = element_text(size=16),
      plot.title = element_text(size=20),
      legend.spacing.x = unit(0.9,"cm")
    )+
    labs(x="",y="AUROC",title="Clinical-trained VEPs",fill="")+
    coord_flip()
  
)

common_legend <- get_legend(poptuned)


(p1 <- ggarrange(poptuned,popclinical,ncol=2,nrow=1,
                 legend="bottom",
                 common.legend = TRUE)+
    theme(plot.title = element_text(size = 25))
)


ggsave(
  paste0(path,"/","fig3.png"),
  plot=p1,
  width=21.22,
  height = 10.573,
  units = "in", dpi = 300
)


##############################################################################################

######################################## Agreement ###########################################
#### Pairwise agreement across structural regions ################################
pairwise_regions <- read_csv("/home/snouto/projects/genesis/second/agreement/final.agreement.csv")
pairwise_regions <- pairwise_regions %>% filter((first_vep %in% paper_veps$vep) &
                                                  (second_vep %in% paper_veps$vep))


### Generating Heatmaps as requested by Ben's Comments ##########################
# heatmap_data <- dcast(pairwise_regions %>% filter(region == "disordered"), 
#                       first_vep ~ second_vep, value.var = "kappa")

disordered_mat <- spread(pairwise_regions %>% filter(region=="disordered") %>% 
                           select(-region), key = first_vep, value = kappa) %>%
  column_to_rownames(var = "second_vep") %>%
  as.matrix()

disordered_mat <- disordered_mat[order(rownames(disordered_mat)), order(colnames(disordered_mat),
                                                                        decreasing = T)]

(disordered <- ggplot(data = as.data.frame(as.table(disordered_mat)), aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#210cae") +
  theme_minimal() +
  labs(x = "", y = "",fill="Kappa") +
  guides(fill = guide_colorbar(barwidth = 15, barheight = 2,
                               label.theme = element_text(size=15)
                               ),
         
         ) +
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text= element_text(size=16),
        legend.title = element_text(size=20)
        )
)


ggsave(
  paste0(path,"/","disordered.png"),
  plot=disordered,
  width=25.22,
  height = 20.573,
  units = "in", dpi = 300
)



intermediate_mat <- spread(pairwise_regions %>% filter(region=="intermediate") %>% select(-region), 
                           key = first_vep, value = kappa) %>%
  column_to_rownames(var = "second_vep") %>%
  as.matrix()

intermediate_mat <- intermediate_mat[order(rownames(intermediate_mat)), order(colnames(intermediate_mat),
                                                                        decreasing = T)]


(intermediate <- ggplot(data = as.data.frame(as.table(intermediate_mat)), 
                        aes(x = Var1, y = Var2, fill = Freq)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "#210cae") +
    theme_minimal() +
    labs(x = "", y = "",fill="Kappa") +
    guides(fill = guide_colorbar(barwidth = 15, barheight = 2,
                                 label.theme = element_text(size=15)
    ),
    
    ) +
    theme_pubr()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size=20),
          axis.text= element_text(size=16))
)



ggsave(
  paste0(path,"/","intermediate.png"),
  plot=intermediate,
  width=25.22,
  height = 20.573,
  units = "in", dpi = 300
)


ordered_mat <- spread(pairwise_regions %>% filter(region=="ordered") %>% select(-region), 
                           key = first_vep, value = kappa) %>%
  column_to_rownames(var = "second_vep") %>%
  as.matrix()

ordered_mat <- ordered_mat[order(rownames(ordered_mat)), order(colnames(ordered_mat),
                                                                              decreasing = T)]


(ordered <- ggplot(data = as.data.frame(as.table(ordered_mat)), 
                        aes(x = Var1, y = Var2, fill = Freq)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "#210cae") +
    theme_minimal() +
    labs(x = "", y = "",fill="Kappa") +
    guides(fill = guide_colorbar(barwidth = 15, barheight = 2,
                                 label.theme = element_text(size=15)
    ),
    
    ) +
    theme_pubr()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size=20),
          axis.text= element_text(size=16))
)


ggsave(
  paste0(path,"/","ordered.png"),
  plot=ordered,
  width=25.22,
  height = 20.573,
  units = "in", dpi = 300
)

#################################################################################
