library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(car)
library(rstatix)

rm(list = ls())
setwd('/home/sasha/Science/PhD/Laue_lab/Microscopy/Confocal/20241101_differentiation_HP1b_DNA_Airyscan_Devina/')

##### Read in data #####
files = str_remove(list.files(file.path('Automatic_quantification'), 'perCell'), '_perCell.csv')
for (f in files) {
  #per cell
  {
    ##total
    perCell = read.csv(file.path('Automatic_quantification', paste0(f, '_perCell.csv')))[-1]
    perCell %>% mutate(tag = rep(c('DNA_mask', 'HP1_mask', 'DNA', 'HP1'), n()/4)) %>%
      pivot_wider(names_from = tag, values_from = c(Mean, Min, Max, X.Area)) %>%
      select(Area, Mean_DNA, Mean_HP1, Min_DNA, Min_HP1, Max_DNA, Max_HP1, X.Area_DNA_mask, X.Area_HP1_mask) %>%
      mutate(Proportion_DNA = X.Area_DNA_mask/100, Proportion_HP1 = X.Area_HP1_mask/100) %>%
      select(-X.Area_DNA_mask, -X.Area_HP1_mask) %>%
      mutate(Cell = 1:n(), Condition = str_split(f, '_')[[1]][1], Pos = str_split(f, '_')[[1]][2]) %>%
      select(Condition, Pos, Cell, Area, Mean_DNA, Mean_HP1, Min_DNA, Min_HP1, Max_DNA, Max_HP1, Proportion_DNA, Proportion_HP1) -> perCell

    ##foci
    perCell_DNA = read.csv(file.path('Automatic_quantification', paste0(f, '_focusSummary_DNA.csv')))[-1]
    perCell_DNA %>% select(-X.Area) %>%
      rename(Num_DNA_foci = Count, Area_DNA_foci_tot = Total.Area, Area_DNA_foci_avg = Average.Size,
             Intensity_DNA_foci_avgFoc = Mean, Num_DNA_foci_HP1 = HP1_count) -> perCell_DNA
    perCell_HP1 = read.csv(file.path('Automatic_quantification', paste0(f, '_focusSummary_HP1.csv')))[-1]
    perCell_HP1 %>% select(-X.Area) %>%
      rename(Num_HP1_foci = Count, Area_HP1_foci_tot = Total.Area, Area_HP1_foci_avg = Average.Size,
             Intensity_HP1_foci_avgFoc = Mean, Num_HP1_foci_DNA = DNA_count) -> perCell_HP1
    
    ##outside foci
    perCell_outDNA = read.csv(file.path('Automatic_quantification', paste0(f, '_outsideFoci_DNA.csv')))[-1]
    perCell_outDNA %>% rename(Area_DNA_outfoci_tot = Area, Intensity_DNA_outfoci_avg = Mean,
                              Intensity_DNA_outfoci_min = Min, Intensity_DNA_outfoci_max = Max) -> perCell_outDNA
    perCell_outHP1 = read.csv(file.path('Automatic_quantification', paste0(f, '_outsideFoci_HP1.csv')))[-1]
    perCell_outHP1 %>% rename(Area_HP1_outfoci_tot = Area, Intensity_HP1_outfoci_avg = Mean,
                              Intensity_HP1_outfoci_min = Min, Intensity_HP1_outfoci_max = Max) -> perCell_outHP1
    
    cbind(perCell, perCell_DNA, perCell_HP1, perCell_outDNA, perCell_outHP1) -> perCell_all
    # select(perCell_all, Area, Area_DNA_foci_tot, Area_HP1_foci_tot, Area_DNA_outfoci_tot, Area_HP1_outfoci_tot) %>%
    #   mutate(Area_DNA_tot = Area_DNA_foci_tot + Area_DNA_outfoci_tot,
    #          Area_HP1_tot = Area_HP1_foci_tot + Area_HP1_outfoci_tot)
  }
  
  #per focus
  {
    perFocus_DNA = read.csv(file.path('Automatic_quantification', paste0(f, '_perFocus_DNA.csv')))[-1]
    celltag = c()
    for (c in 1:nrow(perCell_all)) {
      celltag = c(celltag, rep(c, perCell_all$Num_DNA_foci[c]))
    }
    perFocus_DNA %>%
      mutate(Cell = celltag, Condition = str_split(f, '_')[[1]][1], Pos = str_split(f, '_')[[1]][2]) %>%
      select(Condition, Pos, Cell, Area, Mean, Min, Max, HP1) %>%
      rename(Mean_DNA = Mean, Min_DNA = Min, Max_DNA = Max, HP1_present = HP1) -> perFocus_DNA
    perFocus_HP1 = read.csv(file.path('Automatic_quantification', paste0(f, '_perFocus_HP1.csv')))[-1]
    celltag = c()
    for (c in 1:nrow(perCell_all)) {
      celltag = c(celltag, rep(c, perCell_all$Num_HP1_foci[c]))
    }
    perFocus_HP1 %>%
      mutate(Cell = celltag, Condition = str_split(f, '_')[[1]][1], Pos = str_split(f, '_')[[1]][2]) %>%
      select(Condition, Pos, Cell, Area, Mean, Min, Max, DNA) %>%
      rename(Mean_HP1 = Mean, Min_HP1 = Min, Max_HP1 = Max, DNA_present = DNA) -> perFocus_HP1
  }
  
  #per cell additional calc
  {
    #perFocus_DNA %>% group_by(Condition, Pos, Cell) %>% summarise_all(mean)
    #perFocus_HP1 %>% group_by(Condition, Pos, Cell) %>% summarise_all(mean)
    perFocus_DNA %>% group_by(Condition, Pos, Cell) %>%
      summarise(Intensity_DNA_foci_avg = sum(Mean_DNA*Area)/sum(Area),
                Intensity_DNA_foci_min = min(Min_DNA), Intensity_DNA_foci_max = max(Max_DNA)) %>%
      left_join(perCell_all, .) -> perCell_all
    perFocus_HP1 %>% group_by(Condition, Pos, Cell) %>%
      summarise(Intensity_HP1_foci_avg = sum(Mean_HP1*Area)/sum(Area),
                Intensity_HP1_foci_min = min(Min_HP1), Intensity_HP1_foci_max = max(Max_HP1)) %>%
      left_join(perCell_all, .) -> perCell_all
    
    
    perFocus_DNA %>% group_by(Condition, Pos, Cell, HP1_present) %>%
      summarise(Intensity_DNA_foci_HP1_avgFoc = mean(Mean_DNA),
                Intensity_DNA_foci_HP1_avg = sum(Mean_DNA*Area)/sum(Area),
                Intensity_DNA_foci_HP1_min = min(Min_DNA),
                Intensity_DNA_foci_HP1_max = max(Max_DNA),
                Area_DNA_foci_HP1_tot = sum(Area)) %>%
      pivot_wider(names_from = HP1_present, values_from = c(Intensity_DNA_foci_HP1_avgFoc,
                                                            Intensity_DNA_foci_HP1_avg,
                                                            Intensity_DNA_foci_HP1_min,
                                                            Intensity_DNA_foci_HP1_max,
                                                            Area_DNA_foci_HP1_tot)) %>%
      rename_with(~str_replace(.x, 'HP1', 'noHP1'), ends_with('_0')) %>%
      rename_with(~str_remove(.x, '_\\d')) %>%
      left_join(perCell_all, .) -> perCell_all
    
    perFocus_HP1 %>% group_by(Condition, Pos, Cell, DNA_present) %>%
      summarise(Intensity_HP1_foci_DNA_avgFoc = mean(Mean_HP1),
                Intensity_HP1_foci_DNA_avg = sum(Mean_HP1*Area)/sum(Area),
                Intensity_HP1_foci_DNA_min = min(Min_HP1),
                Intensity_HP1_foci_DNA_max = max(Max_HP1),
                Area_HP1_foci_DNA_tot = sum(Area)) %>%
      pivot_wider(names_from = DNA_present, values_from = c(Intensity_HP1_foci_DNA_avgFoc,
                                                            Intensity_HP1_foci_DNA_avg,
                                                            Intensity_HP1_foci_DNA_min,
                                                            Intensity_HP1_foci_DNA_max,
                                                            Area_HP1_foci_DNA_tot)) %>%
      rename_with(~str_replace(.x, 'DNA', 'noDNA'), ends_with('_0')) %>%
      rename_with(~str_remove(.x, '_\\d')) %>%
      left_join(perCell_all, .) -> perCell_all
  }
  
  # mutate(perCell_all) %>% mutate(Rep = d) %>% relocate(Rep) -> perCell_all
  # mutate(perFocus_DNA) %>% mutate(Rep = d) %>% relocate(Rep) -> perFocus_DNA
  # mutate(perFocus_HP1) %>% mutate(Rep = d) %>% relocate(Rep) -> perFocus_HP1
  
  if (f == files[1]) {
    data_perCell = perCell_all
    data_perFocus_DNA = perFocus_DNA
    data_perFocus_HP1 = perFocus_HP1
  } else {
    data_perCell = rbind(data_perCell, perCell_all)
    data_perFocus_DNA = rbind(data_perFocus_DNA, perFocus_DNA)
    data_perFocus_HP1 = rbind(data_perFocus_HP1, perFocus_HP1)
  }
}
data_perCell = mutate(data_perCell, Condition = factor(Condition, levels = c('2iL', '24h', '48h')))
data_perFocus_DNA = mutate(data_perFocus_DNA, Condition = factor(Condition, levels = c('2iL', '24h', '48h')))
data_perFocus_HP1 = mutate(data_perFocus_HP1, Condition = factor(Condition, levels = c('2iL', '24h', '48h')))

###### Plot stats per timepoint ######

group_by(data_perCell, Condition) %>% summarise(n())
group_by(data_perFocus_DNA, Condition) %>% summarise(n())
group_by(data_perFocus_HP1, Condition) %>% summarise(n())

###Relative intensity of HP1 foci
{
#Both eu and het
  data_perFocus_HP1 %>% left_join(select(data_perCell, Condition, Pos, Cell, Intensity_HP1_outfoci_avg)) %>%
    mutate(Mean_HP1_rel = Mean_HP1/Intensity_HP1_outfoci_avg) -> tmp_perFocus_HP1
  tmp_perFocus_HP1 %>% group_by(Condition, DNA_present) %>% summarise(n())
  ggplot(tmp_perFocus_HP1, aes(x = Condition, y = Mean_HP1_rel, fill = as.character(DNA_present))) +
    geom_violin(position = position_dodge(width = 0.9)) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = c('red', 'purple'), labels = c('DNA-low', 'DNA-high')) +
    #facet_wrap(~Rep, scales = 'free_y') +
    lims(y = c(0, NA)) +
    labs(x = '', y = 'Ratio (Foci intensity)/\n(Euchr intensity)', fill = '',
         title = 'Relative intensity of HP1b foci') +
    theme_bw() +
    theme(aspect.ratio = 1.2, plot.title = element_text(hjust = 0.5))
  # ggsave('HP1b_euhet_relintensity_diff.svg', width = 4, height = 3)
  qqnorm(filter(tmp_perFocus_HP1, DNA_present == 0, Condition == '2iL')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1, DNA_present == 0, Condition == '2iL')$Mean_HP1_rel)#non-norm
  qqnorm(filter(tmp_perFocus_HP1, DNA_present == 0, Condition == '24h')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1, DNA_present == 0, Condition == '24h')$Mean_HP1_rel)#non-norm
  qqnorm(filter(tmp_perFocus_HP1, DNA_present == 0, Condition == '48h')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1, DNA_present == 0, Condition == '48h')$Mean_HP1_rel)#non-norm
  qqnorm(filter(tmp_perFocus_HP1, DNA_present == 1, Condition == '2iL')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1, DNA_present == 1, Condition == '2iL')$Mean_HP1_rel)#maybe
  qqnorm(filter(tmp_perFocus_HP1, DNA_present == 1, Condition == '24h')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1, DNA_present == 1, Condition == '24h')$Mean_HP1_rel)#non-norm
  qqnorm(filter(tmp_perFocus_HP1, DNA_present == 1, Condition == '48h')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1, DNA_present == 1, Condition == '48h')$Mean_HP1_rel)#non-norm
  shapiro.test(filter(tmp_perFocus_HP1, DNA_present == 0, Condition == '2iL')$Mean_HP1_rel)#non-norm
  shapiro.test(filter(tmp_perFocus_HP1, DNA_present == 0, Condition == '24h')$Mean_HP1_rel)#same
  shapiro.test(filter(tmp_perFocus_HP1, DNA_present == 0, Condition == '48h')$Mean_HP1_rel)#same
  shapiro.test(filter(tmp_perFocus_HP1, DNA_present == 1, Condition == '2iL')$Mean_HP1_rel)#non-norm
  shapiro.test(filter(tmp_perFocus_HP1, DNA_present == 1, Condition == '24h')$Mean_HP1_rel)#same
  shapiro.test(filter(tmp_perFocus_HP1, DNA_present == 1, Condition == '48h')$Mean_HP1_rel)#same
  leveneTest(Mean_HP1_rel ~ Condition, data = filter(tmp_perFocus_HP1, DNA_present == 0))#equal var
  leveneTest(Mean_HP1_rel ~ Condition, data = filter(tmp_perFocus_HP1, DNA_present == 1))#unequal var
  #non-norm, large sample, equal var => anova
  summary(aov(Mean_HP1_rel ~ Condition, data = filter(tmp_perFocus_HP1, DNA_present == 0)))#e-9
  TukeyHSD(aov(Mean_HP1_rel ~ Condition, data = filter(tmp_perFocus_HP1, DNA_present == 0)))#24-48h ns
  #non-norm, large sample, unequal var => Welch's anova
  oneway.test(Mean_HP1_rel ~ Condition, data = filter(tmp_perFocus_HP1, DNA_present == 1))#e-16
  games_howell_test(Mean_HP1_rel ~ Condition, data = filter(tmp_perFocus_HP1, DNA_present == 1), conf.level = 0.95)#all sig
}

###Abs intensity of HP1
{
  data_perCell %>% group_by(Condition) %>%
    summarise(mean(Mean_HP1), mean(Intensity_HP1_outfoci_avg), mean(Intensity_HP1_foci_avg))
  
  ggplot(data_perFocus_HP1, aes(x = Condition, y = Mean_HP1, fill = as.character(DNA_present))) +
    geom_violin(position = position_dodge(width = 0.9)) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = c('red', 'purple'), labels = c('DNA-low', 'DNA-high')) +
    #facet_grid(~Rep) +
    lims(y = c(0, NA)) +
    #scale_x_discrete(labels = c('WT', 'Cbx1+/-')) +
    #coord_cartesian(ylim = c(0, 2.2)) +
    labs(y = 'Fluorescence intensity', x = '', fill = '',
         title = 'Intensity of HP1b foci') +
    theme_bw() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
  ggsave('HP1b_euhet_intensity_diff.svg', width = 5, height = 3)
  
  qqnorm(filter(data_perFocus_HP1, DNA_present == 0, Condition == '2iL')$Mean_HP1)
  qqline(filter(data_perFocus_HP1, DNA_present == 0, Condition == '2iL')$Mean_HP1)#non-norm
  qqnorm(filter(data_perFocus_HP1, DNA_present == 0, Condition == '24h')$Mean_HP1)
  qqline(filter(data_perFocus_HP1, DNA_present == 0, Condition == '24h')$Mean_HP1)#non-norm
  qqnorm(filter(data_perFocus_HP1, DNA_present == 0, Condition == '48h')$Mean_HP1)
  qqline(filter(data_perFocus_HP1, DNA_present == 0, Condition == '48h')$Mean_HP1)#non-norm
  shapiro.test(filter(data_perFocus_HP1, DNA_present == 0, Condition == '2iL')$Mean_HP1)#non-norm, large sample
  shapiro.test(filter(data_perFocus_HP1, DNA_present == 0, Condition == '24h')$Mean_HP1)#non-norm, large sample
  shapiro.test(filter(data_perFocus_HP1, DNA_present == 0, Condition == '48h')$Mean_HP1)#non-norm, large sample
  bartlett.test(Mean_HP1 ~ Condition, data = filter(data_perFocus_HP1, DNA_present == 0))#non-equal var
  leveneTest(Mean_HP1 ~ Condition, data = filter(data_perFocus_HP1, DNA_present == 0))#non-equal var
  #non-norm, large sample, unequal var => Welch's anova
  oneway.test(Mean_HP1 ~ Condition, data = filter(data_perFocus_HP1, DNA_present == 0), var.equal = F) #e-16
  games_howell_test(Mean_HP1 ~ Condition, data = filter(data_perFocus_HP1, DNA_present == 0), conf.level = 0.95)
  
  qqnorm(filter(data_perFocus_HP1, DNA_present == 1, Condition == '2iL')$Mean_HP1)
  qqline(filter(data_perFocus_HP1, DNA_present == 1, Condition == '2iL')$Mean_HP1)#non-norm
  qqnorm(filter(data_perFocus_HP1, DNA_present == 1, Condition == '24h')$Mean_HP1)
  qqline(filter(data_perFocus_HP1, DNA_present == 1, Condition == '24h')$Mean_HP1)#non-norm
  qqnorm(filter(data_perFocus_HP1, DNA_present == 1, Condition == '48h')$Mean_HP1)
  qqline(filter(data_perFocus_HP1, DNA_present == 1, Condition == '48h')$Mean_HP1)#non-norm
  shapiro.test(filter(data_perFocus_HP1, DNA_present == 1, Condition == '2iL')$Mean_HP1)#non-norm, large sample
  shapiro.test(filter(data_perFocus_HP1, DNA_present == 1, Condition == '24h')$Mean_HP1)#non-norm, large sample
  shapiro.test(filter(data_perFocus_HP1, DNA_present == 1, Condition == '48h')$Mean_HP1)#non-norm, large sample
  bartlett.test(Mean_HP1 ~ Condition, data = filter(data_perFocus_HP1, DNA_present == 1))#non-equal var
  leveneTest(Mean_HP1 ~ Condition, data = filter(data_perFocus_HP1, DNA_present == 1))#non-equal var
  #non-norm, large sample, unequal var => Welch's t-test
  oneway.test(Mean_HP1 ~ Condition, data = filter(data_perFocus_HP1, DNA_present == 1), var.equal = F) #all very sig
  games_howell_test(Mean_HP1 ~ Condition, data = filter(data_perFocus_HP1, DNA_present == 1), conf.level = 0.95)
}

###Abs intensity of euchr
{
  data_perCell %>% group_by(Condition) %>%
    summarise(mean(Mean_HP1), mean(Intensity_HP1_outfoci_avg), mean(Intensity_HP1_foci_avg), n())
  ggplot(data_perCell, aes(x = Condition, y = Intensity_HP1_outfoci_avg)) +
    geom_violin(fill = 'red') +
    geom_boxplot(width = 0.2, fill = 'red') +
    #facet_grid(~Rep) +
    lims(y = c(0, NA)) +
    #scale_x_discrete(labels = c('WT', 'Cbx1+/-')) +
    labs(y = 'Fluorescence intensity', x = '',
         title = 'HP1b intensity in euchromatin') +
    theme_bw() +
    theme(aspect.ratio = 1.5, plot.title = element_text(hjust = 0.5))
  ggsave('euchr_intensity_diff.svg', width = 3.5, height = 3)
  
  qqnorm(filter(data_perCell, Condition == '2iL')$Intensity_HP1_outfoci_avg)
  qqline(filter(data_perCell, Condition == '2iL')$Intensity_HP1_outfoci_avg)#non-norm
  qqnorm(filter(data_perCell, Condition == '24h')$Intensity_HP1_outfoci_avg)
  qqline(filter(data_perCell, Condition == '24h')$Intensity_HP1_outfoci_avg)#closer to norm
  qqnorm(filter(data_perCell, Condition == '48h')$Intensity_HP1_outfoci_avg)
  qqline(filter(data_perCell, Condition == '48h')$Intensity_HP1_outfoci_avg)#closer to norm
  shapiro.test(filter(data_perCell, Condition == '2iL')$Intensity_HP1_outfoci_avg)#non-norm, p=0.025, not a v large sample
  shapiro.test(filter(data_perCell, Condition == '24h')$Intensity_HP1_outfoci_avg)#norm
  shapiro.test(filter(data_perCell, Condition == '48h')$Intensity_HP1_outfoci_avg)#non-norm, p=0.029, not a v large sample
  bartlett.test(Intensity_HP1_outfoci_avg ~ Condition, data = data_perCell)#non-equal var p=0.0003
  leveneTest(Intensity_HP1_outfoci_avg ~ Condition, data = data_perCell)#non-equal var p=0.021
  #non-norm but not so bad, large(ish) sample, unequal var => Welch's anova
  oneway.test(Intensity_HP1_outfoci_avg ~ Condition, data = data_perCell, var.equal = F) #2e-16
  games_howell_test(Intensity_HP1_outfoci_avg ~ Condition, data = data_perCell, conf.level = 0.95) #2iL sig; 24h-48h ns
}

##Not included in the paper
{
  # 
  # ###Number of eu HP1 foci - same!
  # {
  #   data_perCell %>% mutate(Num_HP1_foci_noDNA = Num_HP1_foci - Num_HP1_foci_DNA) -> data_perCell
  #   ggplot(data_perCell, aes(x = Condition, y = Num_HP1_foci_noDNA)) +
  #     geom_violin() +
  #     geom_boxplot(width = 0.2) +
  #     #facet_grid(~Rep) +
  #     theme_bw()
  #   summary(aov(Num_HP1_foci_noDNA ~ Condition, data = data_perCell)) #p = 0.848
  #   
  #   
  #   ggplot(data_perCell, aes(x = Condition, y = Num_HP1_foci_noDNA)) +
  #     geom_violin(fill = 'red') +
  #     geom_boxplot(width = 0.2, fill = 'red') +
  #     geom_jitter(width = 0.2, shape = 4) +
  #     #scale_x_discrete(labels = c('WT', 'Cbx+/-')) +
  #     labs(x = '', y = 'Number of foci per cell', title = 'Euchromatic HP1b foci per cell') +
  #     theme_bw() +
  #     theme(aspect.ratio = 1.5, plot.title = element_text(hjust = 0.5))
  #   ggsave('euHP1b_perCell_differentiation.svg', width = 4, height = 3)
  #   
  #   qqnorm(filter(data_perCell, Condition == '2iL')$Num_HP1_foci_noDNA)
  #   qqline(filter(data_perCell, Condition == '2iL')$Num_HP1_foci_noDNA)
  #   qqnorm(filter(data_perCell, Condition == '24h')$Num_HP1_foci_noDNA)
  #   qqline(filter(data_perCell, Condition == '24h')$Num_HP1_foci_noDNA)#non-norm-sh?
  #   qqnorm(filter(data_perCell, Condition == '48h')$Num_HP1_foci_noDNA)
  #   qqline(filter(data_perCell, Condition == '48h')$Num_HP1_foci_noDNA)
  #   shapiro.test(filter(data_perCell, Condition == '2iL')$Num_HP1_foci_noDNA)
  #   shapiro.test(filter(data_perCell, Condition == '24h')$Num_HP1_foci_noDNA)#non-norm, p=0.037
  #   shapiro.test(filter(data_perCell, Condition == '48h')$Num_HP1_foci_noDNA)#non-norm, p=0.001
  #   leveneTest(Num_HP1_foci_noDNA ~ Condition, data = data_perCell)#equal var
  #   #non-norm, large sample, equal var => anova
  #   summary(aov(Num_HP1_foci_noDNA ~ Condition, data = data_perCell))
  # }
  # 
  # ###Area foci cover - euHP1 foci cover less at 24h; HP1+DNA+ raw area increases over time (but not perc area)
  # {
  #   data_perCell %>% select(Condition, Pos, Cell, starts_with('Area_HP1_foci')) %>%
  #     pivot_longer(cols = c(Area_HP1_foci_noDNA_tot, Area_HP1_foci_DNA_tot),
  #                  names_to = 'Focus_type', values_to = 'Area') %>%
  #     mutate(Focus_type = case_when(grepl('noDNA', Focus_type) ~ 'DNA-low',
  #                                   T ~ 'DNA-high')) %>%
  #     ggplot(aes(x = Condition, y = Area, fill = Focus_type)) +
  #     geom_violin() +
  #     geom_boxplot(width = 0.2) +
  #     facet_wrap(Focus_type~., scales = 'free_y') +
  #     theme_bw()
  #   summary(aov(Area_HP1_foci_DNA_tot ~ Condition, data = data_perCell))#p=0.03
  #   summary(aov(Area_HP1_foci_noDNA_tot ~ Condition, data = data_perCell))#p=0.007
  #   
  #   data_perCell %>% select(Condition, Pos, Cell, starts_with('Area')) %>%
  #     mutate(Perc_HP1_foci_noDNA = Area_HP1_foci_noDNA_tot/Area,
  #            Perc_HP1_foci_DNA = Area_HP1_foci_DNA_tot/Area) %>%
  #     select(Condition, Pos, Cell, starts_with('Perc')) %>%
  #     pivot_longer(cols = c(Perc_HP1_foci_noDNA, Perc_HP1_foci_DNA),
  #                  names_to = 'Focus_type', values_to = 'Perc_area') %>%
  #     mutate(Focus_type = case_when(grepl('noDNA', Focus_type) ~ 'DNA-low',
  #                                   T ~ 'DNA-high')) %>%
  #     ggplot(aes(x = Condition, y = Perc_area, fill = Focus_type)) +
  #     geom_violin() +
  #     geom_boxplot(width = 0.2) +
  #     facet_wrap(Focus_type~., scales = 'free_y') +
  #     theme_bw()
  #   summary(aov(Area_HP1_foci_DNA_tot/Area ~ Condition, data = data_perCell))#non-sig
  #   summary(aov(Area_HP1_foci_noDNA_tot/Area ~ Condition, data = data_perCell))#2e-06
  #   
  #   data_perCell %>% select(Condition, Pos, Cell, starts_with('Area')) %>%
  #     mutate(Perc_HP1_foci_noDNA = Area_HP1_foci_noDNA_tot/Area) %>%
  #     select(Condition, Pos, Cell, starts_with('Perc')) %>%
  #     ggplot(aes(x = Condition, y = Perc_HP1_foci_noDNA)) +
  #     geom_violin() +
  #     geom_boxplot(width = 0.2) +
  #     #facet_grid(~Rep, scales = 'free_y') +
  #     #coord_cartesian(ylim = c(0, 0.06)) +
  #     theme_bw()
  #   
  #   data_perCell %>% select(Condition, Pos, Cell, starts_with('Area_DNA_foci')) %>%
  #     pivot_longer(cols = c(Area_DNA_foci_noHP1_tot, Area_DNA_foci_HP1_tot),
  #                  names_to = 'Focus_type', values_to = 'Area') %>%
  #     mutate(Focus_type = case_when(grepl('noHP1', Focus_type) ~ 'HP1-low',
  #                                   T ~ 'HP1-high')) %>%
  #     ggplot(aes(x = Condition, y = Area, fill = Focus_type)) +
  #     geom_violin() +
  #     geom_boxplot(width = 0.2) +
  #     facet_wrap(Focus_type~., scales = 'free_y') +
  #     theme_bw()
  #   summary(aov(Area_DNA_foci_HP1_tot ~ Condition, data = data_perCell))#p=0.0007
  #   summary(aov(Area_DNA_foci_noHP1_tot ~ Condition, data = data_perCell))#p=0.05
  #   
  #   data_perCell %>% select(Condition, Pos, Cell, starts_with('Area')) %>%
  #     mutate(Perc_DNA_foci_noHP1 = Area_DNA_foci_noHP1_tot/Area,
  #            Perc_DNA_foci_HP1 = Area_DNA_foci_HP1_tot/Area) %>%
  #     select(Condition, Pos, Cell, starts_with('Perc')) %>%
  #     pivot_longer(cols = c(Perc_DNA_foci_noHP1, Perc_DNA_foci_HP1),
  #                  names_to = 'Focus_type', values_to = 'Perc_area') %>%
  #     mutate(Focus_type = case_when(grepl('noHP1', Focus_type) ~ 'HP1-low',
  #                                   T ~ 'HP1-high')) %>%
  #     ggplot(aes(x = Condition, y = Perc_area, fill = Focus_type)) +
  #     geom_violin() +
  #     geom_boxplot(width = 0.2) +
  #     facet_wrap(Focus_type~., scales = 'free_y') +
  #     theme_bw()
  #   summary(aov(Area_DNA_foci_HP1_tot/Area ~ Condition, data = data_perCell))#ns
  #   summary(aov(Area_DNA_foci_noHP1_tot/Area ~ Condition, data = data_perCell))#ns
  #   
  #   data_perCell %>% select(Condition, Pos, Cell, starts_with('Area')) %>%
  #     mutate(Perc_DNA_foci_HP1VsnoHP1 = Area_DNA_foci_HP1_tot/(Area_DNA_foci_HP1_tot+Area_DNA_foci_noHP1_tot)) %>%
  #     ggplot(aes(x = Condition, y = Perc_DNA_foci_HP1VsnoHP1)) +
  #     #geom_violin() +
  #     geom_boxplot(width = 0.2) +
  #     geom_jitter() +
  #     theme_bw()
  #   summary(aov(Area_DNA_foci_HP1_tot/(Area_DNA_foci_HP1_tot+Area_DNA_foci_noHP1_tot) ~ Condition,
  #               data = data_perCell))#ns
  # }
}


