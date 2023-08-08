library(tidyverse)
source("00_functions.R", echo=TRUE)

setwd("~/Dropbox/sMon/=Odonata_stan_environ_spline_CV")
myfolder <- "42919735"#ran on HPC

#means
meanFiles <- list.files(myfolder, full.names=TRUE) %>%
                str_subset("CVmean") %>%
                set_names() %>%
                map_dfr(readRDS, .id="source") %>%
                group_by(source) %>%
                mutate(species = strsplit(source,"_")[[1]][2],
                       model = strsplit(source,"_")[[1]][3],
                       fold = strsplit(source,"_")[[1]][4],
                       fold = gsub(".rds","",fold)) %>%
                ungroup() %>%
                mutate(model = as.factor(model))%>%
                filter(species %in% selectSpecies)

length(unique(meanFiles$species))

#get means 
summaryData <- meanFiles %>%
                  pivot_longer(AUC:TSS) %>%
                  group_by(source,species,model,name)%>%
                  summarise(value=mean(value))

(g2 <- summaryData %>%
  ggplot() +
  geom_point(aes(x=species, y = value), size=1) +
  coord_flip() +facet_wrap(~name, scales="free_x")+
  theme(axis.text.y = element_text(size=6)))
  
#combine with bpv plot
g1 <- readRDS("plots/bpv_prop.rds")
g1 <- g1 + theme(axis.title = element_text(size=8))
cowplot::plot_grid(g1,g2,
                   nrow=2,
                   labels=c("A","B"),
                   rel_heights = c(0.3,0.7),
                   scale=c(0.9,1))

ggsave("plots/model_checks.png",height=10,width=5)
