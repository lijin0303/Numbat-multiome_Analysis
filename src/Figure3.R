setwd("~/numbat_EpiMultiome/numbat-multiome_Analysis/")
source("src/mini_import.R")
source("src/vis.R")
min_cells <- 50
cellAnnot <- readRDS("intmd/cellAnnot_patA.rds")
clone_post <- readRDS(glue("intmd/patA_CombinedBin_clones_2.rds"))
clone = purrr::keep(clone_post, function(x) x$size > min_cells)
clone_cell <- map(clone,\(x) x$cells)
clone_assign <- list2flat(clone_cell)
colnames(clone_assign) <- c("cell",x)
clone_assign %<>% inner_join(cellAnnot,by="cell")

clonecols <- c("gray89",
               c("#fab81b","#6e4a8f",
                 "#af0627","#017101"))
scales::show_col(clonecols)
names(clonecols) <- paste0("Clone",1:5)  
clonep <- ggplot(clone_assign %>% 
                   mutate(stage = case_when(sampleID=="NIH-A-CLL-BM-01"~"CLL",TRUE~"RT")) %>% 
                   mutate(clone = paste0("Clone",patA_comb_bincnt)) %>% 
                   group_by(stage,clone) %>% 
                   summarise(n=n()), 
                 aes(fill=clone, y=n, x=stage)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values= clonecols) +
  ggtitle("Clone Composition of disease stage") +
  theme_bw() +
  scale_x_discrete(labels = c("Before","After"))+
  rremove("legend.title")+
  theme(
    axis.line = element_line(colour = "black",size=1),
    axis.ticks = element_line(size=1,color="black"),
    axis.text = element_text(color="black"),
    axis.ticks.length=unit(0.2,"cm"),
    legend.position = c(0.8,0.8),
    legend.background=element_blank())+
  font("title",size=rel(1.4))+ 
  font("xy",size=rel(1.1))+ 
  font("xy.text", size = rel(1.5)) +  
  font("legend.text",size = rel(1.2))+
  labs(x="\n Richter's Transformation",y="")+
  guides(fill = guide_legend(override.aes = 
                               list(alpha = 1,
                                    color="black"),
                             ncol = 1))
ggsave("Figures/clone_composition.jpeg",clonep,height=5,width=5)