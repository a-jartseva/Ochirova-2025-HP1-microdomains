library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(car)
library(rstatix)

rm(list = ls())
setwd('/home/sasha/Science/PhD/Laue_lab/Microscopy/Confocal/Histone_marks')


##### Read in data #####
dirs = c('20220201_HP1b-JF646_Hoechst_H3K9me3-Cy3_CAICAiry',
         '20231009_HP1b-H-646_DNA-505_H3K9me3-Cy3_MacheskyConfocal')
repnames = c('Rep1_Airy', 'Rep2_conf')

for (i in 1:length(dirs)) {
  d = dirs[i]
  files = str_remove(list.files(file.path(d, 'Automatic_quantification'), 'DNA_mask'), '_DNA_mask.tif')
  for (f in files) {
    #per cell
    {
      ##colocalisation
      coloc = cbind(read.table(file.path(d, 'Automatic_quantification', paste0(f, '_DNA_H3K9me3_coloc.txt')))[2],
                    read.table(file.path(d, 'Automatic_quantification', paste0(f, '_DNA_HP1_coloc.txt')))[2],
                    read.table(file.path(d, 'Automatic_quantification', paste0(f, '_HP1_H3K9me3_coloc.txt')))[2])
      names(coloc) = c('R_DNA_H3K9me3', 'R_DNA_HP1', 'R_HP1_H3K9me3')
      coloc %>% mutate(Cell = 1:n(), Pos = str_split(f, '_')[[1]][1]) %>% relocate(Pos, Cell) -> coloc
      
      ##foci
      perCell_DNA = read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_H3K9me3_DNAfocusSummary.csv')))[-1]
      perCell_DNA %>% select(-X.Area) %>%
        rename(Num_DNA_foci = Count, Area_DNA_foci_tot = Total.Area, Area_DNA_foci_avg = Average.Size,
               H3K9me3_intensity_DNA_foci_avgFoc = Mean, Num_DNA_foci_HP1 = HP1_count) -> perCell_DNA
      # read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_H3K9me3norm_DNAfocusSummary.csv')))[6] %>%
      #   rename(H3K9me3overDNA_intensity_DNA_foci_avgFoc = Mean) %>%
      #   cbind(perCell_DNA) -> perCell_DNA
      
      perCell_HP1 = read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_H3K9me3_HP1focusSummary.csv')))[-1]
      perCell_HP1 %>% select(-X.Area) %>%
        rename(Num_HP1_foci = Count, Area_HP1_foci_tot = Total.Area, Area_HP1_foci_avg = Average.Size,
               H3K9me3_intensity_HP1_foci_avgFoc = Mean, Num_HP1_foci_DNA = DNA_count) -> perCell_HP1
      # read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_H3K9me3norm_HP1focusSummary.csv')))[6] %>%
      #   rename(H3K9me3overDNA_intensity_HP1_foci_avgFoc = Mean) %>%
      #   cbind(perCell_HP1) -> perCell_HP1
      
      ##outside foci
      perCell_outside = read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_H3K9me3_outsideFoci.csv')))[-1]
      perCell_outside %>% rename(Area_outfoci_tot = Area, H3K9me3_intensity_outfoci_avg = Mean,
                                 H3K9me3_intensity_outfoci_min = Min, H3K9me3_intensity_outfoci_max = Max) -> perCell_outside
      # read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_H3K9me3norm_outsideFoci.csv')))[-(1:2)] %>%
      #   rename(H3K9me3overDNA_intensity_outfoci_avg = Mean,
      #          H3K9me3overDNA_intensity_outfoci_min = Min, H3K9me3overDNA_intensity_outfoci_max = Max) %>%
      #   cbind(perCell_outside) -> perCell_outside
      read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_DNA_outsideFoci.csv')))[-(1:2)] %>%
        rename(DNA_intensity_outfoci_avg = Mean,
               DNA_intensity_outfoci_min = Min, DNA_intensity_outfoci_max = Max) %>%
        cbind(perCell_outside) -> perCell_outside
      
      cbind(perCell_DNA, perCell_HP1, perCell_outside) %>%
        mutate(Cell = 1:n(), Pos = str_split(f, '_')[[1]][1]) %>% relocate(Pos, Cell) -> perCell_all
      # select(perCell_all, Area, Area_DNA_foci_tot, Area_HP1_foci_tot, Area_DNA_outfoci_tot, Area_HP1_outfoci_tot) %>%
      #   mutate(Area_DNA_tot = Area_DNA_foci_tot + Area_DNA_outfoci_tot,
      #          Area_HP1_tot = Area_HP1_foci_tot + Area_HP1_outfoci_tot)
    }
    
    #per focus
    {
      perFocus_DNA = read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_H3K9me3_perDNAFocus.csv')))[-1]
      celltag = c()
      for (c in 1:nrow(perCell_all)) {
        celltag = c(celltag, rep(c, perCell_all$Num_DNA_foci[c]))
      }
      perFocus_DNA %>%
        mutate(Cell = celltag, Pos = str_split(f, '_')[[1]][1]) %>%
        select(Pos, Cell, Area, Mean, Min, Max, HP1) %>%
        rename(H3K9me3_mean = Mean, H3K9me3_min = Min, H3K9me3_max = Max, HP1_present = HP1) -> perFocus_DNA
      # read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_H3K9me3norm_perDNAFocus.csv')))[-(1:2)] %>%
      #   rename(H3K9me3overDNA_mean = Mean, H3K9me3overDNA_min = Min, H3K9me3overDNA_max = Max) %>%
      #   cbind(perFocus_DNA, .) -> perFocus_DNA
      read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_DNA_perDNAFocus.csv')))[-(1:2)] %>%
        rename(DNA_mean = Mean, DNA_min = Min, DNA_max = Max) %>% cbind(perFocus_DNA, .) -> perFocus_DNA
      
      perFocus_HP1 = read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_H3K9me3_perHP1Focus.csv')))[-1]
      celltag = c()
      for (c in 1:nrow(perCell_all)) {
        celltag = c(celltag, rep(c, perCell_all$Num_HP1_foci[c]))
      }
      perFocus_HP1 %>%
        mutate(Cell = celltag, Pos = str_split(f, '_')[[1]][1]) %>%
        select(Pos, Cell, Area, Mean, Min, Max, DNA) %>%
        rename(H3K9me3_mean = Mean, H3K9me3_min = Min, H3K9me3_max = Max, DNA_present = DNA) -> perFocus_HP1
      # read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_H3K9me3norm_perHP1Focus.csv')))[-(1:2)] %>%
      #   rename(H3K9me3overDNA_mean = Mean, H3K9me3overDNA_min = Min, H3K9me3overDNA_max = Max) %>%
      #   cbind(perFocus_HP1, .) -> perFocus_HP1
      read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_DNA_perHP1Focus.csv')))[-(1:2)] %>%
        rename(DNA_mean = Mean, DNA_min = Min, DNA_max = Max) %>% cbind(perFocus_HP1, .) -> perFocus_HP1
    }
    
    #per cell additional calc
    {
      #perFocus_DNA %>% group_by(HP1_expr, Pos, Cell) %>% summarise_all(mean)
      #perFocus_HP1 %>% group_by(HP1_expr, Pos, Cell) %>% summarise_all(mean)
      perFocus_DNA %>% group_by(Pos, Cell) %>%
        summarise(H3K9me3_intensity_DNA_foci_avg = sum(H3K9me3_mean*Area)/sum(Area),
                  H3K9me3overDNA_intensity_DNA_foci_avg = sum(H3K9me3_mean*Area/DNA_mean)/sum(Area)) %>%
                  #Intensity_DNA_foci_min = min(Min_DNA), Intensity_DNA_foci_max = max(Max_DNA)) %>%
        left_join(perCell_all, .) -> perCell_all
      perFocus_HP1 %>% group_by(Pos, Cell) %>%
        summarise(H3K9me3_intensity_HP1_foci_avg = sum(H3K9me3_mean*Area)/sum(Area),
                  H3K9me3overDNA_intensity_HP1_foci_avg = sum(H3K9me3_mean*Area/DNA_mean)/sum(Area)) %>%
                  #Intensity_HP1_foci_min = min(Min_HP1), Intensity_HP1_foci_max = max(Max_HP1)) %>%
        left_join(perCell_all, .) -> perCell_all
      
      
      perFocus_DNA %>% group_by(Pos, Cell, HP1_present) %>%
        summarise(H3K9me3_intensity_DNA_foci_HP1_avgFoc = mean(H3K9me3_mean),
                  H3K9me3_intensity_DNA_foci_HP1_avg = sum(H3K9me3_mean*Area)/sum(Area),
                  H3K9me3overDNA_intensity_DNA_foci_HP1_avgFoc = mean(H3K9me3_mean/DNA_mean),
                  H3K9me3overDNA_intensity_DNA_foci_HP1_avg = sum(H3K9me3_mean*Area/DNA_mean)/sum(Area),
                  Area_DNA_foci_HP1_tot = sum(Area)) %>%
        pivot_wider(names_from = HP1_present, values_from = c(starts_with('H3K9me3'),
                                                              Area_DNA_foci_HP1_tot)) %>%
        rename_with(~str_replace(.x, 'HP1', 'noHP1'), ends_with('_0')) %>%
        rename_with(~str_remove(.x, '_\\d')) %>%
        left_join(perCell_all, .) -> perCell_all
      
      perFocus_HP1 %>% group_by(Pos, Cell, DNA_present) %>%
        summarise(H3K9me3_intensity_HP1_foci_DNA_avgFoc = mean(H3K9me3_mean),
                  H3K9me3_intensity_HP1_foci_DNA_avg = sum(H3K9me3_mean*Area)/sum(Area),
                  H3K9me3overDNA_intensity_HP1_foci_DNA_avgFoc = mean(H3K9me3_mean/DNA_mean),
                  H3K9me3overDNA_intensity_HP1_foci_DNA_avg = sum(H3K9me3_mean*Area/DNA_mean)/sum(Area),
                  Area_HP1_foci_DNA_tot = sum(Area)) %>%
        pivot_wider(names_from = DNA_present, values_from = c(starts_with('H3K9me3'),
                                                              Area_HP1_foci_DNA_tot)) %>%
        rename_with(~str_replace(.x, 'DNA', 'noDNA'), ends_with('_0')) %>%
        rename_with(~str_remove(.x, '_\\d')) %>%
        left_join(perCell_all, .) -> perCell_all
    }
    
    mutate(coloc) %>% mutate(Rep = repnames[i]) %>% relocate(Rep) -> coloc
    mutate(perCell_all) %>% mutate(Rep = repnames[i]) %>% relocate(Rep) -> perCell_all
    mutate(perFocus_DNA) %>% mutate(Rep = repnames[i]) %>% relocate(Rep) -> perFocus_DNA
    mutate(perFocus_HP1) %>% mutate(Rep = repnames[i]) %>% relocate(Rep) -> perFocus_HP1
    
    if (d == dirs[1] & f == files[1]) {
      data_coloc = coloc
      data_perCell = perCell_all
      data_perFocus_DNA = perFocus_DNA
      data_perFocus_HP1 = perFocus_HP1
    } else {
      data_coloc = rbind(data_coloc, coloc)
      data_perCell = rbind(data_perCell, perCell_all)
      data_perFocus_DNA = rbind(data_perFocus_DNA, perFocus_DNA)
      data_perFocus_HP1 = rbind(data_perFocus_HP1, perFocus_HP1)
    }
  }
}

##Num obs
group_by(data_coloc, Rep) %>% summarise(n())
group_by(data_perFocus_DNA, HP1_present) %>% summarise(n())
group_by(data_perFocus_HP1, DNA_present) %>% summarise(n())

#### Plot coloc ####

data_coloc %>% pivot_longer(cols = starts_with('R_'), names_to = 'Channels',
                            names_prefix = 'R_', values_to = 'R') %>%
  mutate(Channels = factor(Channels, levels = c('DNA_HP1', 'DNA_H3K9me3', 'HP1_H3K9me3'))) -> coloc_long
ggplot(coloc_long, aes(x = Channels, y = R, fill = Channels)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  #geom_jitter(width = 0.2, shape = 4) +
  facet_grid(~Rep) +
  scale_fill_manual(values = c('DNA_HP1' = 'purple', 'DNA_H3K9me3' = 'limegreen', 'HP1_H3K9me3' = 'orange')) +
  labs(y = 'Pearson\'s R', title = 'Colocalisation between DNA, HP1b and H3K9me3') +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = 'none', axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave('Coloc.svg', width = 5, height = 3)
{
qqnorm(filter(data_coloc, Rep == 'Rep1_Airy')$R_DNA_HP1)
qqline(filter(data_coloc, Rep == 'Rep1_Airy')$R_DNA_HP1)
qqnorm(filter(data_coloc, Rep == 'Rep2_conf')$R_DNA_HP1)
qqline(filter(data_coloc, Rep == 'Rep2_conf')$R_DNA_HP1)#norm
shapiro.test(filter(data_coloc, Rep == 'Rep1_Airy')$R_DNA_HP1)
shapiro.test(filter(data_coloc, Rep == 'Rep2_conf')$R_DNA_HP1)#norm
qqnorm(filter(data_coloc, Rep == 'Rep1_Airy')$R_DNA_H3K9me3)
qqline(filter(data_coloc, Rep == 'Rep1_Airy')$R_DNA_H3K9me3)
qqnorm(filter(data_coloc, Rep == 'Rep2_conf')$R_DNA_H3K9me3)
qqline(filter(data_coloc, Rep == 'Rep2_conf')$R_DNA_H3K9me3)#norm
shapiro.test(filter(data_coloc, Rep == 'Rep1_Airy')$R_DNA_H3K9me3)
shapiro.test(filter(data_coloc, Rep == 'Rep2_conf')$R_DNA_H3K9me3)#norm
qqnorm(filter(data_coloc, Rep == 'Rep1_Airy')$R_HP1_H3K9me3)
qqline(filter(data_coloc, Rep == 'Rep1_Airy')$R_HP1_H3K9me3)
qqnorm(filter(data_coloc, Rep == 'Rep2_conf')$R_HP1_H3K9me3)
qqline(filter(data_coloc, Rep == 'Rep2_conf')$R_HP1_H3K9me3)#norm
shapiro.test(filter(data_coloc, Rep == 'Rep1_Airy')$R_HP1_H3K9me3)
shapiro.test(filter(data_coloc, Rep == 'Rep2_conf')$R_HP1_H3K9me3)#norm
bartlett.test(R ~ Channels, data = filter(coloc_long, Rep == 'Rep1_Airy'))#non-equal var
bartlett.test(R ~ Channels, data = filter(coloc_long, Rep == 'Rep2_conf'))#non-equal var
}
{
#norm, large sample, non-equal var
#repeated measures! ie non-independent => multiple measures ANOVA
anova_test(data = mutate(filter(coloc_long, Rep == 'Rep1_Airy'), ID = paste(Pos, Cell, sep = '_')),
           dv = R, wid = ID, within = Channels)
#Mauchly's test indicates non-sphericity => need a correction
#Highly significant with sphericity corrections:
##p_Greenhouse-Geisser = 6.86e-44, p_Huynh-Feldt = 7.61e-45
#Post-hoc: pairwise t-tests with Hochberg multiple testing correction
t.test(filter(data_coloc, Rep == 'Rep1_Airy')$R_DNA_H3K9me3,
       filter(data_coloc, Rep == 'Rep1_Airy')$R_DNA_HP1, paired = T) #p<2.2e-16
t.test(filter(data_coloc, Rep == 'Rep1_Airy')$R_DNA_H3K9me3,
       filter(data_coloc, Rep == 'Rep1_Airy')$R_HP1_H3K9me3, paired = T) #p=0.07568
t.test(filter(data_coloc, Rep == 'Rep1_Airy')$R_HP1_H3K9me3,
       filter(data_coloc, Rep == 'Rep1_Airy')$R_DNA_HP1, paired = T) #p<2.2e-16
p.adjust(c(2.2e-16, 2.2e-16, 0.07568), method = 'hochberg') #same result
#=> DNA_H3K9me3 vs HP1_H3K9me3 non-signif, others highly signif

anova_test(data = mutate(filter(coloc_long, Rep == 'Rep2_conf'), ID = paste(Pos, Cell, sep = '_')),
           dv = R, wid = ID, within = Channels)
#Mauchly's test indicates sphericity, even though the variances are not the same
#Highly significant either with or without sphericity corrections:
##p = 1.47e-37, p_Greenhouse-Geisser = 3.36e-34, p_Huynh-Feldt = 1.66e-35)
#Post-hoc: pairwise t-tests with Hochberg multiple testing correction
t.test(filter(data_coloc, Rep == 'Rep2_conf')$R_DNA_H3K9me3,
       filter(data_coloc, Rep == 'Rep2_conf')$R_DNA_HP1, paired = T) #p<2.2e-16
t.test(filter(data_coloc, Rep == 'Rep2_conf')$R_DNA_H3K9me3,
       filter(data_coloc, Rep == 'Rep2_conf')$R_HP1_H3K9me3, paired = T) #p=0.004201
t.test(filter(data_coloc, Rep == 'Rep2_conf')$R_HP1_H3K9me3,
       filter(data_coloc, Rep == 'Rep2_conf')$R_DNA_HP1, paired = T) #p<2.2e-16
p.adjust(c(2.2e-16, 2.2e-16, 0.004201), method = 'hochberg') #same result
#=> DNA_H3K9me3 vs HP1_H3K9me3 signif, others highly signif
}
#=> purple highly signif, g vs o signif in Rep2

#### Plot stats in foci, per focus ####

#Adjust by euchromatin intensity
data_perFocus_DNA %>% left_join(select(data_perCell, Rep, Pos, Cell,
                                       H3K9me3_intensity_outfoci_avg, DNA_intensity_outfoci_avg)) %>%
  mutate(H3K9me3_mean_norm = H3K9me3_mean/H3K9me3_intensity_outfoci_avg,
         DNA_mean_norm = DNA_mean/DNA_intensity_outfoci_avg) %>%
  mutate(H3K9me3overDNA_mean_norm = H3K9me3_mean_norm/DNA_mean_norm) -> data_perFocus_DNA
data_perFocus_HP1 %>% left_join(select(data_perCell, Rep, Pos, Cell,
                                       H3K9me3_intensity_outfoci_avg, DNA_intensity_outfoci_avg)) %>%
  mutate(H3K9me3_mean_norm = H3K9me3_mean/H3K9me3_intensity_outfoci_avg,
         DNA_mean_norm = DNA_mean/DNA_intensity_outfoci_avg) %>%
  mutate(H3K9me3overDNA_mean_norm = H3K9me3_mean_norm/DNA_mean_norm) -> data_perFocus_HP1

#Simply H3K9me3
ggplot(data_perFocus_DNA, aes(x = as.character(HP1_present), y = H3K9me3_mean_norm,
                              fill = as.character(HP1_present))) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  #facet_grid(~Rep) +
  scale_x_discrete(labels = c('DNA-rich,\nHP1b-poor', 'DNA-rich,\nHP1b-rich')) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 9)) +
  scale_fill_manual(values = c('1' = 'purple', '0' = 'cyan')) +
  geom_hline(yintercept = 1, colour = 'black', size = 1, linetype = 'dashed') +
  theme_bw() +
  theme(aspect.ratio = 1.2, legend.position = 'none', axis.title.x = element_blank())
ggplot(data_perFocus_HP1, aes(x = as.character(DNA_present), y = H3K9me3_mean_norm,
                              fill = as.character(DNA_present))) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  #facet_grid(~Rep) +
  scale_x_discrete(labels = c('DNA-poor,\nHP1b-rich', 'DNA-rich,\nHP1b-rich')) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 9)) +
  scale_fill_manual(values = c('1' = 'purple', '0' = 'red')) +
  geom_hline(yintercept = 1, colour = 'black', size = 1, linetype = 'dashed') +
  theme_bw() +
  theme(aspect.ratio = 1.2, legend.position = 'none', axis.title.x = element_blank())
#The two ways of defining purple foci show similar, although not exactly the same results
#Although it's not exactly correct to substitute them, it will be too confusing to have both
#Thus, chose DNA foci-based distribution, which shows less of a difference with other groups (ie weaker result)
mutate(data_perFocus_DNA, Compartment = case_when(HP1_present == 0 ~ 'DNA_foci',
                                                  HP1_present == 1 ~ 'DNA_HP1_foci')) %>%
  select(-HP1_present) %>%
  rbind(filter(data_perFocus_HP1, DNA_present == 0) %>% mutate(Compartment = 'HP1_foci') %>%
          select(-DNA_present)) %>%
  mutate(Compartment = factor(Compartment, c('HP1_foci', 'DNA_foci', 'DNA_HP1_foci'))) -> data_perFocus_all
ggplot(data_perFocus_all, aes(x = Compartment, y = H3K9me3_mean_norm, fill = Compartment)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  #facet_grid(~Rep) + #reps look similar!
  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 9)) +
  scale_fill_manual(values = c('HP1_foci' = 'red', 'DNA_foci' = 'cyan', 'DNA_HP1_foci' = 'purple')) +
  geom_hline(yintercept = 1, colour = 'black', size = 1, linetype = 'dashed') +
  labs(y = 'H3K9me3 intensity\nrelative to euchromatin', title = 'H3K9me3 enrichment\nin different foci') +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = 'none', axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave('H3K9me3_relIntensity_perFocus.svg', width = 3, height = 3)
{
  qqnorm(filter(data_perFocus_all, Compartment == 'HP1_foci')$H3K9me3_mean_norm)
  qqline(filter(data_perFocus_all, Compartment == 'HP1_foci')$H3K9me3_mean_norm) #non-norm
  qqnorm(filter(data_perFocus_all, Compartment == 'DNA_foci')$H3K9me3_mean_norm)
  qqline(filter(data_perFocus_all, Compartment == 'DNA_foci')$H3K9me3_mean_norm) #non-norm
  qqnorm(filter(data_perFocus_all, Compartment == 'DNA_HP1_foci')$H3K9me3_mean_norm)
  qqline(filter(data_perFocus_all, Compartment == 'DNA_HP1_foci')$H3K9me3_mean_norm) #non-norm
  shapiro.test(filter(data_perFocus_all, Compartment == 'HP1_foci')$H3K9me3_mean_norm) #non-norm
  shapiro.test(filter(data_perFocus_all, Compartment == 'DNA_foci')$H3K9me3_mean_norm) #non-norm
  shapiro.test(filter(data_perFocus_all, Compartment == 'DNA_HP1_foci')$H3K9me3_mean_norm) #non-norm
  leveneTest(H3K9me3_mean_norm ~ Compartment, data = data_perFocus_all)#non-equal var
}
#Large sample, non-norm, non-equal var -> Welch's ANOVA
oneway.test(H3K9me3_mean_norm ~ Compartment, data = data_perFocus_all, var.equal = F) #highly signif
games_howell_test(H3K9me3_mean_norm ~ Compartment, data = data_perFocus_all)#red&cyan 0.001, the other two e-10


#H3K9me3/DNA intensity
ggplot(data_perFocus_DNA, aes(x = as.character(HP1_present), y = H3K9me3overDNA_mean_norm,
                              fill = as.character(HP1_present))) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  #facet_grid(~Rep) +
  scale_x_discrete(labels = c('DNA-rich,\nHP1b-poor', 'DNA-rich,\nHP1b-rich')) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10)) +
  scale_fill_manual(values = c('1' = 'purple', '0' = 'cyan')) +
  geom_hline(yintercept = 1, colour = 'black', size = 1, linetype = 'dashed') +
  theme_bw() +
  theme(aspect.ratio = 1.2, legend.position = 'none', axis.title.x = element_blank())
ggplot(data_perFocus_HP1, aes(x = as.character(DNA_present), y = H3K9me3overDNA_mean_norm,
                              fill = as.character(DNA_present))) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  #facet_grid(~Rep) +
  scale_x_discrete(labels = c('DNA-poor,\nHP1b-rich', 'DNA-rich,\nHP1b-rich')) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10)) +
  scale_fill_manual(values = c('1' = 'purple', '0' = 'red')) +
  geom_hline(yintercept = 1, colour = 'black', size = 1, linetype = 'dashed') +
  theme_bw() +
  theme(aspect.ratio = 1.2, legend.position = 'none', axis.title.x = element_blank())
#The two ways of defining purple foci show similar, although not exactly the same results
#Although it's not exactly correct to substitute them, it will be too confusing to have both
#Thus, chose DNA foci-based distribution, which shows less of a difference with other groups (ie weaker result)
ggplot(data_perFocus_all, aes(x = Compartment, y = H3K9me3overDNA_mean_norm,
             fill = Compartment)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  #facet_grid(~Rep) + #reps look similar!
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  coord_cartesian(ylim = c(0, 6)) +
  scale_fill_manual(values = c('HP1_foci' = 'red', 'DNA_foci' = 'cyan', 'DNA_HP1_foci' = 'purple')) +
  geom_hline(yintercept = 1, colour = 'black', size = 1, linetype = 'dashed') +
  labs(y = 'Ratio (Norm H3K9me3 intensity)/\n(Norm DNA intensity)',
       title = 'H3K9me3 enrichment over DNA\nin different foci') +
  theme_bw() +
  theme(aspect.ratio = 1, legend.position = 'none', axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave('H3K9me3overDNA_relIntensity_perFocus.svg', width = 3.2, height = 3)
sum(data_perFocus_all$H3K9me3overDNA_mean_norm>6)#13 => 10 points outside plot area
{
  qqnorm(filter(data_perFocus_all, Compartment == 'HP1_foci')$H3K9me3overDNA_mean_norm)
  qqline(filter(data_perFocus_all, Compartment == 'HP1_foci')$H3K9me3overDNA_mean_norm) #non-norm
  qqnorm(filter(data_perFocus_all, Compartment == 'DNA_foci')$H3K9me3overDNA_mean_norm)
  qqline(filter(data_perFocus_all, Compartment == 'DNA_foci')$H3K9me3overDNA_mean_norm) #non-norm
  qqnorm(filter(data_perFocus_all, Compartment == 'DNA_HP1_foci')$H3K9me3overDNA_mean_norm)
  qqline(filter(data_perFocus_all, Compartment == 'DNA_HP1_foci')$H3K9me3overDNA_mean_norm) #maybe norm
  shapiro.test(filter(data_perFocus_all, Compartment == 'HP1_foci')$H3K9me3overDNA_mean_norm) #non-norm
  shapiro.test(filter(data_perFocus_all, Compartment == 'DNA_foci')$H3K9me3overDNA_mean_norm) #non-norm
  shapiro.test(filter(data_perFocus_all, Compartment == 'DNA_HP1_foci')$H3K9me3overDNA_mean_norm) #non-norm
  leveneTest(H3K9me3overDNA_mean_norm ~ Compartment, data = data_perFocus_all)#non-equal var
}
#Large sample, non-norm, non-equal var -> Welch's ANOVA
oneway.test(H3K9me3overDNA_mean_norm ~ Compartment, data = data_perFocus_all, var.equal = F) #highly signif
games_howell_test(H3K9me3overDNA_mean_norm ~ Compartment, data = data_perFocus_all)#all e-10

#Is H3K9me3/DNA in DNA_foci signif > or < 1?
#Large sample, non-norm -> t-test
t.test(filter(data_perFocus_DNA, HP1_present == 0)$H3K9me3overDNA_mean_norm,
       mu=1, alternative='two.sided')#signif lower
t.test(filter(data_perFocus_DNA, HP1_present == 0, Rep == 'Rep1_Airy')$H3K9me3overDNA_mean_norm,
       mu=1, alternative='two.sided')#non-signif
t.test(filter(data_perFocus_DNA, HP1_present == 0, Rep == 'Rep2_conf')$H3K9me3overDNA_mean_norm,
       mu=1, alternative='two.sided')#signif lower

shapiro.test(log(filter(data_perFocus_DNA, HP1_present == 0)$H3K9me3overDNA_mean_norm))#a little non-normal
t.test(log(filter(data_perFocus_DNA, HP1_present == 0)$H3K9me3overDNA_mean_norm),
       mu=0, alternative='two.sided')#signif lower
t.test(log(filter(data_perFocus_DNA, HP1_present == 0, Rep == 'Rep1_Airy')$H3K9me3overDNA_mean_norm),
       mu=0, alternative='two.sided')#signif lower
t.test(log(filter(data_perFocus_DNA, HP1_present == 0, Rep == 'Rep2_conf')$H3K9me3overDNA_mean_norm),
       mu=0, alternative='two.sided')#signif lower
