library("tidyverse")

setwd("/run/user/1001/gvfs/sftp:host=u900-bdd-1-124t-6989.curie.fr,user=tgutman/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/colo829")

# Load truth table
truthTsv=read.table("TRUTH/COLO829_truth.tsv",header=TRUE,stringsAsFactors = FALSE,sep="\t")

# Load SV results
svTable=read.table("results_benchmark/colo_results_SV.csv",header=TRUE,stringsAsFactors = FALSE, sep=",")

svTable=svTable %>% 
  mutate(Center = case_when(Center == "BSC" ~ "Node 1",
                            Center == "Curie" ~ "Node 2",
                            Center == "Charite" ~ "Node 3",
                            Center == "Hartwig" ~ "Node 4",
                            Center == "OICR" ~ "Node 5")) %>% 
  mutate(Center=factor(Center,levels=c("Node 1", "Node 2", "Node 3", "Node 4","Node 5")))
  
tidySV=svTable %>% 
  filter(Bin != "All" & TIER == "tier3") %>% 
  select(Bin,Center,starts_with("TP") & !c("TP")) %>% 
  mutate(Bin = case_when(Bin == "Bin NaN bp" ~ "NaN",
                         Bin == "0-50" ~ "0-50",
                         Bin == "50-200" ~ "50-200",
                         Bin == "200-1000" ~ "200-1000",
                         Bin == "> 1000" ~ "> 1000"))


colnames(tidySV)=c("Bin","Center","DEL","INS","DUP","INV","BND")

truthTsv= truthTsv %>% 
  select(new_id,type,size) %>% 
  mutate(bin = cut(truthTsv$size,c(0,1,50,200,1000,1000000000),include.lowest = TRUE))

SV_count=t(table(truthTsv$type,truthTsv$bin))
SV_count = as.data.frame.matrix(SV_count) %>% 
  rownames_to_column() %>% 
  rename(Bin = rowname) %>%
  mutate(Center="Truth") %>% 
  select(Bin,Center,DEL,INS,DUP,INV,BND) %>% 
  mutate(Bin = case_when(Bin == "[0,1]" ~ "NaN",
                         Bin == "(1,50]" ~ "0-50",
                         Bin == "(50,200]" ~ "50-200",
                         Bin == "(200,1e+03]" ~ "200-1000",
                         Bin == "(1e+03,1e+09]" ~ "> 1000")) %>% 
  bind_rows(tidySV)

SV_count$Bin=fct_relevel(as.factor(SV_count$Bin), "0-50", "50-200","200-1000","> 1000","NaN")

# Plot general metrics:
## Plot TP FP FN raw 

tidySvTP=svTable %>% 
  filter(Bin == "All" & TIER == "tier3") %>% 
  select(-c("FP_original", "FP_tier", "Recall", "Precision", "F1.score", "TP_DEL", "TP_INS", "TP_DUP", "TP_INV", "TP_BND")) %>% 
  pivot_longer(cols=c(TP,FP,FN), names_to="metric",values_to="count")

ggplot(tidySvTP,aes(x= Center,y=count,fill=metric)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_brewer(palette="Set1") +
  theme_bw()+
  labs(x = "Centers",y="SV counts") +
  theme(legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=22),
        legend.key.size = unit(1.5,"line"))

ggsave("results_benchmark/barplot_SV_general_tier3_all.png",width=30,height=20,units='cm')

## Plot TP FP FN by SV type

tidySvTP_all=svTable %>% 
  filter(Bin == "All" & TIER == "tier3") %>% 
  select(Bin,Center,starts_with("TP")) %>% 
  pivot_longer(cols=starts_with("TP"),names_to = "metric",values_to = "count")

ggplot(tidySvTP_all,aes(x= Center,y=count,fill=metric)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  labs(x = "Centers",y="SV counts") +
  theme(legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=22),
        legend.key.size = unit(1.5,"line"))

ggsave("results_benchmark/barplot_SV_TP_tier3_all.png",width=30,height=20,units='cm')

## Plot Precision Recall & F1 score
tidySvF1=svTable %>% 
  filter(Bin == "All" & TIER == "tier3") %>% 
  select(-c(TIER,FP,FP_original,FP_tier,FN,TP_DEL,TP_INS,TP_DUP,TP_INV,TP_BND )) %>% 
  pivot_longer(cols = c(Recall,Precision,F1.score), names_to = "metric", values_to = "count")

tidySvF1$metric=fct_relevel(as.factor(tidySvF1$metric), "Recall","Precision","F1.score")

ggplot(tidySvF1,aes(x= Center,y=count*100,fill=metric)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  labs(x = "Centers",y="Score (%)") +
  theme(legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=22),
        legend.key.size = unit(1.5,"line"))

ggsave("results_benchmark/barplot_SV_F1_tier3_all.png",width=30,height=20,units='cm')

## Compare to truth:

tidySV_count=SV_count %>% 
  pivot_longer(cols=c(DEL,INS,DUP,INV,BND),names_to = "type",values_to="count")

ggplot(tidySV_count,aes(x=Bin,y=count,fill=Center)) +
  geom_bar(stat="identity",position="dodge") +
  facet_wrap(~type) +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  theme(legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        axis.text=element_text(size=10),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=22),
        legend.key.size = unit(1.5,"line"))

ggsave("results_benchmark/barplot_SV_truth_types.png",width=30,height=20,units='cm')

ggplot(mapping = aes(x=Bin,y=count,fill=Center)) +
  geom_bar(data=tidySV_count[which(tidySV_count$Center != "Truth"),],stat="identity",position="dodge") +
  geom_bar(data=tidySV_count[which(tidySV_count$Center == "Truth"),],stat="identity",position="dodge",linetype="dashed", colour="black",fill=NA) +
  facet_wrap(~type,scales = "free") +
  scale_fill_brewer(palette="Set1") +
  theme_bw() +
  labs(x = "SV length (bp)",y="SV counts") +
  theme(legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        axis.text=element_text(size=10),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=22),
        strip.text.x = element_text(size = 14),
        legend.key.size = unit(1.5,"line"))

ggsave("results_benchmark/barplot_SV_truth_types.png",width=30,height=20,units='cm')
