
# Load libraries:
library(tidyverse)

# Set WD
setwd("~/Documents/Tom/EUCANCan/Benchmark/colo829")

# Load data
snvTable=read.table("results_benchmark/colo_results_snv.csv", header=TRUE,sep = ",")
snvTable=snvTable[,c("Node","type","total.truth","total.query","tp","fp","fn","recall","precision")]

# Compute F1 + order factors
snvTable=snvTable %>% 
    mutate(F1 = (2*(precision*recall)/(precision+recall))) %>% 
    mutate(type=factor(type,levels=c("SNVs","indels","records"))) %>% 
    mutate(Node = case_when(Node == "BSC" ~ "Node 1",
                            Node == "Curie" ~ "Node 2",
                            Node == "Charite" ~ "Node 3",
                            Node == "Hartwig" ~ "Node 4",
                            Node == "OICR" ~ "Node 5")) %>% 
    mutate(Node=factor(Node,levels=c("Node 1", "Node 2", "Node 3", "Node 4", "Node 5")))

# Transform table + order factors
tidySnv=snvTable %>%
    select(-c(total.truth,total.query,precision,recall,F1)) %>% 
    pivot_longer(cols = c(tp,fp,fn), names_to = "metric", values_to = "count") %>%
    mutate(metric=factor(metric,levels=c("tp","fp","fn")))

# Plot TP FP FN
ggplot(tidySnv,aes(x= Node,y=count,fill=metric)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(~type,scales="free") +
    scale_fill_brewer(palette="Set1") +
    theme_bw()+
    labs(x = "Centers",y="SNV counts") +
    theme(legend.title = element_text(size = 22),
          legend.text = element_text(size = 18),
          axis.text=element_text(size=11),
          axis.title=element_text(size=16,face="bold"),
          plot.title = element_text(size=22),
          strip.text.x = element_text(size = 14),
          legend.key.size = unit(1.5,"line"))

ggsave("results_benchmark/barplotTPFPFN.png",width=30,height=20,units='cm')

# Transform table for Precision recall F1 + order factors
tidySnv=snvTable %>%
    pivot_longer(cols = c(recall,precision,F1), names_to = "metric", values_to = "count") %>% 
    mutate(metric=factor(metric,levels=c("recall","precision","F1")))

# Plot Precision Recall F1
ggplot(tidySnv,aes(x= Node,y=count*100,fill=metric)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(~type,scales="free") +
    scale_fill_brewer(palette="Set1") +
    theme_bw() +
    labs(x = "Centers",y="Score (%)") +
    theme(legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        axis.text=element_text(size=11),
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 14),
        plot.title = element_text(size=22),
        legend.key.size = unit(1.5,"line"))

ggsave("results_benchmark/barplotPrecisionRecallF1.png",width=30,height=20,units='cm')
