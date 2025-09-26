library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(car)

rm(list = ls())
setwd('/home/sasha/Science/PhD/Laue_lab/Microscopy/Confocal/Comparison_HighLow')

###### Read in data ######
dirs = c('20230802', '20230920')

for (d in dirs) {
  files = str_remove(list.files(file.path(d, 'Automatic_quantification'), 'perCell'), '_perCell.csv')
  for (f in files) {
    #per cell
    {
      ##total
      perCell = read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_perCell.csv')))[-1]
      perCell %>% mutate(tag = rep(c('DNA_mask', 'HP1_mask', 'DNA', 'HP1'), n()/4)) %>%
        pivot_wider(names_from = tag, values_from = c(Mean, Min, Max, X.Area)) %>%
        select(Area, Mean_DNA, Mean_HP1, Min_DNA, Min_HP1, Max_DNA, Max_HP1, X.Area_DNA_mask, X.Area_HP1_mask) %>%
        mutate(Proportion_DNA = X.Area_DNA_mask/100, Proportion_HP1 = X.Area_HP1_mask/100) %>%
        select(-X.Area_DNA_mask, -X.Area_HP1_mask) %>%
        mutate(Cell = 1:n(), HP1_expr = str_split(f, '_')[[1]][1], Pos = str_split(f, '_')[[1]][2]) %>%
        select(HP1_expr, Pos, Cell, Area, Mean_DNA, Mean_HP1, Min_DNA, Min_HP1, Max_DNA, Max_HP1, Proportion_DNA, Proportion_HP1) -> perCell

      ##foci
      perCell_DNA = read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_focusSummary_DNA.csv')))[-1]
      perCell_DNA %>% select(-X.Area) %>%
        rename(Num_DNA_foci = Count, Area_DNA_foci_tot = Total.Area, Area_DNA_foci_avg = Average.Size,
               Intensity_DNA_foci_avgFoc = Mean, Num_DNA_foci_HP1 = HP1_count) -> perCell_DNA
      perCell_HP1 = read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_focusSummary_HP1.csv')))[-1]
      perCell_HP1 %>% select(-X.Area) %>%
        rename(Num_HP1_foci = Count, Area_HP1_foci_tot = Total.Area, Area_HP1_foci_avg = Average.Size,
               Intensity_HP1_foci_avgFoc = Mean, Num_HP1_foci_DNA = DNA_count) -> perCell_HP1
      
      ##outside foci
      perCell_outDNA = read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_outsideFoci_DNA.csv')))[-1]
      perCell_outDNA %>% rename(Area_DNA_outfoci_tot = Area, Intensity_DNA_outfoci_avg = Mean,
                                Intensity_DNA_outfoci_min = Min, Intensity_DNA_outfoci_max = Max) -> perCell_outDNA
      perCell_outHP1 = read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_outsideFoci_HP1.csv')))[-1]
      perCell_outHP1 %>% rename(Area_HP1_outfoci_tot = Area, Intensity_HP1_outfoci_avg = Mean,
                                Intensity_HP1_outfoci_min = Min, Intensity_HP1_outfoci_max = Max) -> perCell_outHP1
      
      cbind(perCell, perCell_DNA, perCell_HP1, perCell_outDNA, perCell_outHP1) -> perCell_all
      # select(perCell_all, Area, Area_DNA_foci_tot, Area_HP1_foci_tot, Area_DNA_outfoci_tot, Area_HP1_outfoci_tot) %>%
      #   mutate(Area_DNA_tot = Area_DNA_foci_tot + Area_DNA_outfoci_tot,
      #          Area_HP1_tot = Area_HP1_foci_tot + Area_HP1_outfoci_tot)
    }
    
    #per focus
    {
      perFocus_DNA = read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_perFocus_DNA.csv')))[-1]
      celltag = c()
      for (c in 1:nrow(perCell_all)) {
        celltag = c(celltag, rep(c, perCell_all$Num_DNA_foci[c]))
      }
      perFocus_DNA %>%
        mutate(Cell = celltag, HP1_expr = str_split(f, '_')[[1]][1], Pos = str_split(f, '_')[[1]][2]) %>%
        select(HP1_expr, Pos, Cell, Area, Mean, Min, Max, HP1) %>%
        rename(Mean_DNA = Mean, Min_DNA = Min, Max_DNA = Max, HP1_present = HP1) -> perFocus_DNA
      perFocus_HP1 = read.csv(file.path(d, 'Automatic_quantification', paste0(f, '_perFocus_HP1.csv')))[-1]
      celltag = c()
      for (c in 1:nrow(perCell_all)) {
        celltag = c(celltag, rep(c, perCell_all$Num_HP1_foci[c]))
      }
      perFocus_HP1 %>%
        mutate(Cell = celltag, HP1_expr = str_split(f, '_')[[1]][1], Pos = str_split(f, '_')[[1]][2]) %>%
        select(HP1_expr, Pos, Cell, Area, Mean, Min, Max, DNA) %>%
        rename(Mean_HP1 = Mean, Min_HP1 = Min, Max_HP1 = Max, DNA_present = DNA) -> perFocus_HP1
    }
    
    #per cell additional calc
    {
      #perFocus_DNA %>% group_by(HP1_expr, Pos, Cell) %>% summarise_all(mean)
      #perFocus_HP1 %>% group_by(HP1_expr, Pos, Cell) %>% summarise_all(mean)
      perFocus_DNA %>% group_by(HP1_expr, Pos, Cell) %>%
        summarise(Intensity_DNA_foci_avg = sum(Mean_DNA*Area)/sum(Area),
                  Intensity_DNA_foci_min = min(Min_DNA), Intensity_DNA_foci_max = max(Max_DNA)) %>%
        left_join(perCell_all, .) -> perCell_all
      perFocus_HP1 %>% group_by(HP1_expr, Pos, Cell) %>%
        summarise(Intensity_HP1_foci_avg = sum(Mean_HP1*Area)/sum(Area),
                  Intensity_HP1_foci_min = min(Min_HP1), Intensity_HP1_foci_max = max(Max_HP1)) %>%
        left_join(perCell_all, .) -> perCell_all
      
      
      perFocus_DNA %>% group_by(HP1_expr, Pos, Cell, HP1_present) %>%
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
      
      perFocus_HP1 %>% group_by(HP1_expr, Pos, Cell, DNA_present) %>%
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
    
    mutate(perCell_all) %>% mutate(Rep = d) %>% relocate(Rep) -> perCell_all
    mutate(perFocus_DNA) %>% mutate(Rep = d) %>% relocate(Rep) -> perFocus_DNA
    mutate(perFocus_HP1) %>% mutate(Rep = d) %>% relocate(Rep) -> perFocus_HP1
    
    if (d == dirs[1] & f == files[1]) {
      data_perCell = perCell_all
      data_perFocus_DNA = perFocus_DNA
      data_perFocus_HP1 = perFocus_HP1
    } else {
      data_perCell = rbind(data_perCell, perCell_all)
      data_perFocus_DNA = rbind(data_perFocus_DNA, perFocus_DNA)
      data_perFocus_HP1 = rbind(data_perFocus_HP1, perFocus_HP1)
    }
  }
}


###### Plot stats in WT cells (Figure 1) ######
filter(data_perCell, HP1_expr == 'High') -> data_perCell_h
filter(data_perFocus_DNA, HP1_expr == 'High') -> data_perFocus_DNA_h
filter(data_perFocus_HP1, HP1_expr == 'High') -> data_perFocus_HP1_h

###Perc of nucleus occupied by different foci
{
  mutate(data_perCell_h, Proportion_DNA_HP1 = Area_DNA_foci_HP1_tot/Area,
         Proportion_DNA_noHP1 = Area_DNA_foci_noHP1_tot/Area,
         Proportion_HP1_DNA = Area_HP1_foci_DNA_tot/Area,
         Proportion_HP1_noDNA = Area_HP1_foci_noDNA_tot/Area) -> data_perCell_h
  
  # select(data_perCell_h, Rep:Cell, Proportion_DNA_HP1:Proportion_HP1_noDNA) %>%
  #   replace_na(list(Proportion_DNA_HP1 = 0, Proportion_DNA_noHP1 = 0,
  #                   Proportion_HP1_DNA = 0, Proportion_HP1_noDNA = 0)) %>%
  #   pivot_longer(starts_with('Proportion'), names_to = c('.value', 'Foci_type'), names_pattern = '^([[:alpha:]]*)_(.*)') -> tmp_propPerCell
  # ggplot(tmp_propPerCell, aes(x = Foci_type, y = Proportion)) +
  #   geom_violin() +
  #   geom_boxplot(width = 0.2) +
  #   facet_grid(Rep~.) +
  #   theme_bw()
  
  select(data_perCell_h, Rep:Cell, starts_with('Proportion')) %>%
    replace_na(list(Proportion_DNA = 0, Proportion_HP1 = 0,
                    Proportion_DNA_HP1 = 0, Proportion_DNA_noHP1 = 0,
                    Proportion_HP1_DNA = 0, Proportion_HP1_noDNA = 0)) %>%
    mutate(Proportion_ofDNA_HP1 = Proportion_DNA_HP1/Proportion_DNA,
           Proportion_ofHP1_DNA = Proportion_HP1_DNA/Proportion_HP1) -> tmp_propPerCell
  # ggplot(tmp_propPerCell, aes(x = Proportion_ofDNA_HP1)) +
  #   geom_histogram(binwidth = 0.05) +
  #   facet_grid(Rep~.) +
  #   theme_bw()
  # ggplot(tmp_propPerCell, aes(x = Proportion_ofHP1_DNA)) +
  #   geom_histogram(binwidth = 0.05) +
  #   facet_grid(Rep~.) +
  #   theme_bw()
  # group_by(tmp_propPerCell, Rep) %>% summarise(Avg_prop_ofDNA_HP1 = mean(Proportion_ofDNA_HP1),
  #                                             Avg_prop_ofHP1_DNA = mean(Proportion_ofHP1_DNA))
  # 
  # 
  # ggplot(tmp_propPerCell, aes(x = Proportion_ofDNA_HP1)) +
  #   geom_histogram(binwidth = 0.05) +
  #   #facet_grid(Rep~.) +
  #   scale_x_continuous (labels = scales::percent) +
  #   labs(y = 'Count', x = 'Proportion of heterochromatin rich in HP1b',
  #        title = 'Proportion of heterochromatin\nrich in HP1b per cell',
  #        caption = 'Each focus annotated as HP1b+ or -;\nproportion corrected for area') +
  #   theme_bw() +
  #   theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
  # ggsave('Heterochromatin_HP1+_perCell.svg', width = 3.5, height = 2.53)
  select(tmp_propPerCell, Rep:Cell, Proportion_HP1 = Proportion_ofDNA_HP1) %>%
    arrange(Proportion_HP1) %>% mutate(ID = 1:n()) %>%
    mutate(Proportion_noHP1 = 1-Proportion_HP1) %>%
    pivot_longer(starts_with('Proportion'), names_to = c('.value', 'HP1_content'), names_sep = '_') %>%
    ggplot(aes(x = ID, y = Proportion, fill = HP1_content)) +
    geom_bar(stat = 'identity', colour = 'gray20') +
    #facet_grid(Rep~.) +
    scale_y_continuous (labels = scales::percent) +
    scale_fill_manual(values = c('HP1' = 'purple', 'noHP1' = 'cyan'), labels = c('HP1b-rich', 'HP1b-poor')) +
    labs(y = 'Proportion of het', x = 'Cell',
         title = 'Proportion of heterochromatin\nrich in HP1b per cell',
         caption = 'Each focus annotated as HP1b+ or -;\nproportion corrected for area',
         fill = '') +
    theme_bw() +
    theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
  ggsave('Heterochromatin_HP1+_perCell2.svg', width = 4.5, height = 2.5)
}

###Number of different foci
{
  mutate(data_perCell_h, Num_DNA_foci_noHP1 = Num_DNA_foci-Num_DNA_foci_HP1,
         Num_HP1_foci_noDNA = Num_HP1_foci-Num_HP1_foci_DNA) -> data_perCell_h
  
  # select(data_perCell_h, Rep:Cell, Num_DNA_foci_noHP1, Num_DNA_foci_HP1, Num_HP1_foci_noDNA, Num_HP1_foci_DNA) %>%
  #   pivot_longer(starts_with('Num'), names_to = c('.value', 'Foci_type'), names_pattern = '^([[:alpha:]]*)_(.*)') -> tmp_numPerCell
  # ggplot(tmp_numPerCell, aes(x = Foci_type, y = Num)) +
  #   geom_violin() +
  #   geom_boxplot(width = 0.2) +
  #   facet_grid(Rep~.) +
  #   theme_bw()
  # group_by(tmp_numPerCell, Rep, Foci_type) %>% summarise(Avg_num = mean(Num))
  # ggplot(tmp_numPerCell, aes(x = Foci_type, y = Num)) +
  #   geom_violin() +
  #   geom_boxplot(width = 0.2) +
  #   theme_bw()
  
  ggplot(data_perCell_h, aes(x = Num_HP1_foci_noDNA)) +
    geom_histogram(binwidth = 1, colour = 'gray20', fill = 'red') +
    #facet_grid(Rep~.) +
    labs(y = 'Count', x = 'Number per cell', title = 'Euchromatic HP1b foci per cell') +
    theme_bw() +
    theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
  ggsave('euHP1b_perCell.svg', width = 3, height = 2)
  
  # mutate(data_perCell_h, Prop_DNA_foci_HP1 = Num_DNA_foci_HP1/Num_DNA_foci,
  #        Prop_HP1_foci_DNA = Num_HP1_foci_DNA/Num_HP1_foci) -> data_perCell_h
  # ggplot(data_perCell_h, aes(x = Prop_DNA_foci_HP1)) +
  #   geom_histogram(binwidth = 0.05) +
  #   facet_grid(Rep~.) +
  #   theme_bw()
  # ggplot(data_perCell_h, aes(x = Prop_HP1_foci_DNA)) +
  #   geom_histogram(binwidth = 0.05) +
  #   facet_grid(Rep~.) +
  #   theme_bw()
  # group_by(data_perCell_h, Rep) %>% summarise(Avg_prop_DNA_foci_HP1 = mean(Prop_DNA_foci_HP1),
  #                                             Avg_prop_HP1_foci_DNA = mean(Prop_HP1_foci_DNA))
}

##Other - not included in the paper
{
# ###Intensity of het vs eu HP1 - same!
# {
#   #per cell
#   select(data_perCell_h, Rep:Cell, Intensity_HP1_foci_DNA_avg, Intensity_HP1_foci_noDNA_avg) %>%
#     rename_with(~str_remove(.x, '_HP1_foci')) %>% rename_with(~str_remove(.x, '_avg')) %>%
#     pivot_longer(starts_with('Intensity'), names_to = c('.value', 'DNA_present'), names_sep = '_') -> tmp_fociPerCell
#   ggplot(tmp_fociPerCell, aes(x = DNA_present, y = Intensity)) +
#     geom_violin() +
#     geom_boxplot(width = 0.2) +
#     facet_grid(~Rep) +
#     theme_bw()
#   t.test(Intensity ~ DNA_present, data = filter(tmp_fociPerCell, Rep == '20230802'))
#   t.test(Intensity ~ DNA_present, data = filter(tmp_fociPerCell, Rep == '20230920'))
#   
#   #per focus
#   ggplot(data_perFocus_HP1_h, aes(x = as.character(DNA_present), y = Mean_HP1)) +
#     geom_violin() +
#     geom_boxplot(width = 0.2) +
#     facet_grid(~Rep) +
#     theme_bw()
#   t.test(Mean_HP1 ~ DNA_present, data = filter(data_perFocus_HP1_h, Rep == '20230802'))
#   t.test(Mean_HP1 ~ DNA_present, data = filter(data_perFocus_HP1_h, Rep == '20230920'))
# }
# 
# ###Size of het vs eu HP1 - obviously, different
# {
#   #per cell
#   data_perFocus_HP1_h %>% group_by(Rep, HP1_expr, Pos, Cell, DNA_present) %>%
#     summarise(Area_avg = mean(Area)) -> tmp_HP1_area
#   ggplot(tmp_HP1_area, aes(x = as.character(DNA_present), y = Area_avg)) +
#     geom_violin() +
#     geom_boxplot(width = 0.2) +
#     facet_grid(~Rep) +
#     theme_bw()
#   t.test(Area_avg ~ DNA_present, data = filter(tmp_HP1_area, Rep == '20230802'))
#   t.test(Area_avg ~ DNA_present, data = filter(tmp_HP1_area, Rep == '20230920'))
#   
#   #per focus
#   ggplot(data_perFocus_HP1_h, aes(x = as.character(DNA_present), y = Area)) +
#     geom_violin() +
#     geom_boxplot(width = 0.2) +
#     facet_grid(~Rep) +
#     theme_bw()
#   t.test(Area ~ DNA_present, data = filter(data_perFocus_HP1_h, Rep == '20230802'))
#   t.test(Area ~ DNA_present, data = filter(data_perFocus_HP1_h, Rep == '20230920'))
# }
# 
# ###Intensity of DNA +-HP1 - same in Rep1, different in Rep2
# {
#   #per cell
#   select(data_perCell_h, Rep:Cell, Intensity_DNA_foci_HP1_avg, Intensity_DNA_foci_noHP1_avg) %>%
#     rename_with(~str_remove(.x, '_DNA_foci')) %>% rename_with(~str_remove(.x, '_avg')) %>%
#     pivot_longer(starts_with('Intensity'), names_to = c('.value', 'HP1_present'), names_sep = '_') -> tmp_fociPerCell
#   ggplot(tmp_fociPerCell, aes(x = HP1_present, y = Intensity)) +
#     geom_violin() +
#     geom_boxplot(width = 0.2) +
#     facet_wrap(~Rep, scales = 'free_y') +
#     theme_bw()
#   t.test(Intensity ~ HP1_present, data = filter(tmp_fociPerCell, Rep == '20230802'))
#   t.test(Intensity ~ HP1_present, data = filter(tmp_fociPerCell, Rep == '20230920'))
#   
#   #per focus
#   ggplot(data_perFocus_DNA_h, aes(x = as.character(HP1_present), y = Mean_DNA)) +
#     geom_violin() +
#     geom_boxplot(width = 0.2) +
#     facet_wrap(~Rep, scales = 'free_y') +
#     theme_bw()
#   t.test(Mean_DNA ~ HP1_present, data = filter(data_perFocus_DNA_h, Rep == '20230802'))
#   t.test(Mean_DNA ~ HP1_present, data = filter(data_perFocus_DNA_h, Rep == '20230920'))
# }
}

###### Plot stats WT vs het cells (Figure 5) ######

group_by(data_perCell, HP1_expr, Rep) %>% summarise(n())
group_by(data_perFocus_DNA, HP1_expr, Rep) %>% summarise(n())
group_by(data_perFocus_HP1, HP1_expr, Rep) %>% summarise(n())

###Number of eu HP1 foci
{
  data_perCell %>% mutate(Num_HP1_foci_noDNA = Num_HP1_foci - Num_HP1_foci_DNA) -> data_perCell
  ggplot(data_perCell, aes(x = HP1_expr, y = Num_HP1_foci_noDNA)) +
    geom_violin() +
    geom_boxplot(width = 0.2) +
    facet_grid(~Rep) +
    theme_bw()
  t.test(Num_HP1_foci_noDNA ~ HP1_expr, data = filter(data_perCell, Rep == '20230802'))
  t.test(Num_HP1_foci_noDNA ~ HP1_expr, data = filter(data_perCell, Rep == '20230920'))
  #same trend in both reps
  
  
  ggplot(data_perCell, aes(x = HP1_expr, y = Num_HP1_foci_noDNA)) +
    geom_violin(fill = 'red') +
    geom_boxplot(width = 0.2, fill = 'red') +
    geom_jitter(width = 0.2, shape = 4) +
    scale_x_discrete(labels = c('WT', 'Cbx+/-')) +
    labs(x = '', y = 'Number of foci per cell', title = 'Euchromatic HP1b foci per cell') +
    theme_bw() +
    theme(aspect.ratio = 1.5, plot.title = element_text(hjust = 0.5))
  ggsave('euHP1b_perCell_highLow.svg', width = 3, height = 3)
  
  qqnorm(filter(data_perCell, HP1_expr == 'High')$Num_HP1_foci_noDNA)
  qqline(filter(data_perCell, HP1_expr == 'High')$Num_HP1_foci_noDNA)
  qqnorm(filter(data_perCell, HP1_expr == 'Low')$Num_HP1_foci_noDNA)
  qqline(filter(data_perCell, HP1_expr == 'Low')$Num_HP1_foci_noDNA)#non-norm
  shapiro.test(filter(data_perCell, HP1_expr == 'High')$Num_HP1_foci_noDNA)
  shapiro.test(filter(data_perCell, HP1_expr == 'Low')$Num_HP1_foci_noDNA)#non-norm
  leveneTest(Num_HP1_foci_noDNA ~ HP1_expr, data = data_perCell)#equal var
  #non-norm, large sample, equal var => Student's t-test
  t.test(Num_HP1_foci_noDNA ~ HP1_expr, data = data_perCell, var.equal = T)
}

###Relative intensity of foci
{
#Both eu and het
  data_perFocus_HP1 %>% left_join(select(data_perCell, Rep, HP1_expr, Pos, Cell, Intensity_HP1_outfoci_avg)) %>%
    mutate(Mean_HP1_rel = Mean_HP1/Intensity_HP1_outfoci_avg) -> tmp_perFocus_HP1
  tmp_perFocus_HP1 %>% group_by(HP1_expr, Rep, DNA_present) %>% summarise(n())
  ggplot(tmp_perFocus_HP1, aes(x = HP1_expr, y = Mean_HP1_rel, fill = as.character(DNA_present))) +
    geom_violin(position = position_dodge(width = 0.9)) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = c('red', 'purple'), labels = c('DNA-low', 'DNA-high')) +
    scale_x_discrete(labels = c('WT', 'Cbx1+/-')) +
    facet_wrap(~Rep, scales = 'free_y') +
    lims(y = c(0, NA)) +
    labs(x = '', y = 'Ratio (Foci intensity)/\n(Euchr intensity)', fill = '',
         title = 'Relative intensity of HP1b foci') +
    theme_bw() +
    theme(aspect.ratio = 1.2, plot.title = element_text(hjust = 0.5))
  ggsave('HP1b_euhet_relintensity_highLow_reps.svg', width = 6, height = 3)
  qqnorm(filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230802', HP1_expr == 'High')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230802', HP1_expr == 'High')$Mean_HP1_rel)#02/08 almost norm, 20/09 less so
  qqnorm(filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230802', HP1_expr == 'Low')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230802', HP1_expr == 'Low')$Mean_HP1_rel)#02/08 almost norm, 20/09 less so
  shapiro.test(filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230802', HP1_expr == 'High')$Mean_HP1_rel)#same
  shapiro.test(filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230802', HP1_expr == 'Low')$Mean_HP1_rel)#same
  qqnorm(filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230920', HP1_expr == 'High')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230920', HP1_expr == 'High')$Mean_HP1_rel)#02/08 almost norm, 20/09 less so
  qqnorm(filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230920', HP1_expr == 'Low')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230920', HP1_expr == 'Low')$Mean_HP1_rel)#02/08 almost norm, 20/09 less so
  shapiro.test(filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230920', HP1_expr == 'High')$Mean_HP1_rel)#same
  shapiro.test(filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230920', HP1_expr == 'Low')$Mean_HP1_rel)#same
  bartlett.test(Mean_HP1_rel ~ HP1_expr, data = filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230802'))#equal var
  leveneTest(Mean_HP1_rel ~ HP1_expr, data = filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230920'))#equal var
  
  qqnorm(filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230802', HP1_expr == 'High')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230802', HP1_expr == 'High')$Mean_HP1_rel)#not norm
  qqnorm(filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230802', HP1_expr == 'Low')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230802', HP1_expr == 'Low')$Mean_HP1_rel)#almost norm
  shapiro.test(filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230802', HP1_expr == 'High')$Mean_HP1_rel)#not norm
  shapiro.test(filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230802', HP1_expr == 'Low')$Mean_HP1_rel)#not norm
  qqnorm(filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230920', HP1_expr == 'High')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230920', HP1_expr == 'High')$Mean_HP1_rel)#not norm
  qqnorm(filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230920', HP1_expr == 'Low')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230920', HP1_expr == 'Low')$Mean_HP1_rel)#tails
  shapiro.test(filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230920', HP1_expr == 'High')$Mean_HP1_rel)#not norm
  shapiro.test(filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230920', HP1_expr == 'Low')$Mean_HP1_rel)#not norm
  leveneTest(Mean_HP1_rel ~ HP1_expr, data = filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230802'))#equal var
  leveneTest(Mean_HP1_rel ~ HP1_expr, data = filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230920'))#non-equal var, p = 0.003
  
  #almost norm, large sample, equal var => Student's t-test
  t.test(Mean_HP1_rel ~ HP1_expr, data = filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230802'), var.equal = T)#0.0014, but increase just from 1.695 to 1.767
  t.test(Mean_HP1_rel ~ HP1_expr, data = filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230920'), var.equal = T)#0.0017, but decrease just from 2.176 to 2.080
  t.test(Mean_HP1_rel ~ HP1_expr, data = filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230802'), var.equal = T)#0.833
  t.test(Mean_HP1_rel ~ HP1_expr, data = filter(tmp_perFocus_HP1, DNA_present == 0, Rep == '20230920'), var.equal = T)#0.019
  #non-equal var => also try Welch's t-test
  t.test(Mean_HP1_rel ~ HP1_expr, data = filter(tmp_perFocus_HP1, DNA_present == 1, Rep == '20230920'), var.equal = F)#0.0032, but decrease just from 2.176 to 2.080
}

###Abs intensity of HP1 - ~two-fold difference for everything
{
  ##All HP1 foci, separated by DNA
  #Quick check - same trend between reps, but values vary
  ggplot(filter(data_perFocus_HP1), aes(x = HP1_expr, y = Mean_HP1, fill = as.character(DNA_present))) +
    geom_violin(position = position_dodge(width = 0.9)) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9)) +
    facet_wrap(~Rep, scales = 'free_y') +
    lims(y = c(0, NA)) +
    labs(fill = 'DNA') +
    theme_bw()
  #Normalise between reps
  data_perFocus_HP1 %>% group_by(Rep) %>% mutate(Mean_HP1_rel = Mean_HP1/mean(Mean_HP1)) -> tmp_perFocus_HP1_rel
  ggplot(tmp_perFocus_HP1_rel, aes(x = HP1_expr, y = Mean_HP1_rel, fill = as.character(DNA_present))) +
    geom_violin(position = position_dodge(width = 0.9)) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = c('red', 'purple'), labels = c('DNA-low', 'DNA-high')) +
    #facet_grid(~Rep) +
    lims(y = c(0, NA)) +
    scale_x_discrete(labels = c('WT', 'Cbx1+/-')) +
    coord_cartesian(ylim = c(0, 2.2)) +
    labs(y = 'Fluorescence intensity, AU', x = '', fill = '',
         title = 'Intensity of HP1b foci') +
    theme_bw() +
    theme(aspect.ratio = 1.2, plot.title = element_text(hjust = 0.5))
  ggsave('HP1b_euhet_intensity_highLow.svg', width = 4, height = 3)
  
  qqnorm(filter(tmp_perFocus_HP1_rel, DNA_present == 0, HP1_expr == 'High')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1_rel, DNA_present == 0, HP1_expr == 'High')$Mean_HP1_rel)#a little non-norm?
  qqnorm(filter(tmp_perFocus_HP1_rel, DNA_present == 0, HP1_expr == 'Low')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1_rel, DNA_present == 0, HP1_expr == 'Low')$Mean_HP1_rel)
  shapiro.test(filter(tmp_perFocus_HP1_rel, DNA_present == 0, HP1_expr == 'High')$Mean_HP1_rel)#non-norm, large sample
  shapiro.test(filter(tmp_perFocus_HP1_rel, DNA_present == 0, HP1_expr == 'Low')$Mean_HP1_rel)#non-norm, large sample
  bartlett.test(Mean_HP1_rel ~ HP1_expr, data = filter(tmp_perFocus_HP1_rel, DNA_present == 0))#non-equal var
  leveneTest(Mean_HP1_rel ~ HP1_expr, data = filter(tmp_perFocus_HP1_rel, DNA_present == 0))#non-equal var
  #almost norm, large sample, unequal var => Welch's t-test
  t.test(Mean_HP1_rel ~ HP1_expr, data = tmp_perFocus_HP1_rel, DNA_present == 0, var.equal = F)
  
  qqnorm(filter(tmp_perFocus_HP1_rel, DNA_present == 1, HP1_expr == 'High')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1_rel, DNA_present == 1, HP1_expr == 'High')$Mean_HP1_rel)#non-norm?
  qqnorm(filter(tmp_perFocus_HP1_rel, DNA_present == 1, HP1_expr == 'Low')$Mean_HP1_rel)
  qqline(filter(tmp_perFocus_HP1_rel, DNA_present == 1, HP1_expr == 'Low')$Mean_HP1_rel)#non-norm?
  shapiro.test(filter(tmp_perFocus_HP1_rel, HP1_expr == 'High', DNA_present == 1)$Mean_HP1_rel)#non-norm, large sample
  shapiro.test(filter(tmp_perFocus_HP1_rel, HP1_expr == 'Low', DNA_present == 1)$Mean_HP1_rel)#non-norm, large sample
  #bartlett.test(Mean_HP1_rel ~ HP1_expr, data = tmp_perFocus_HP1_noDNA_rel)#non-equal var
  leveneTest(Mean_HP1_rel ~ HP1_expr, data = filter(tmp_perFocus_HP1_rel, DNA_present == 1))#non-equal var
  #almost norm, large sample, unequal var => Welch's t-test
  t.test(Mean_HP1_rel ~ HP1_expr, data = filter(tmp_perFocus_HP1_rel, DNA_present == 1), var.equal = F)
}

###Area foci cover - not included in paper
{
  # data_perCell %>% select(Rep, HP1_expr, Pos, Cell, starts_with('Area_HP1_foci')) %>%
  #   pivot_longer(cols = c(Area_HP1_foci_noDNA_tot, Area_HP1_foci_DNA_tot),
  #                names_to = 'Focus_type', values_to = 'Area') %>%
  #   mutate(Focus_type = case_when(grepl('noDNA', Focus_type) ~ 'DNA-low',
  #                                 T ~ 'DNA-high')) %>%
  #   ggplot(aes(x = HP1_expr, y = Area, fill = Focus_type)) +
  #   geom_violin() +
  #   geom_boxplot(width = 0.2) +
  #   facet_wrap(Focus_type~Rep, scales = 'free_y') +
  #   theme_bw()
  # t.test(Area_HP1_foci_DNA_tot ~ HP1_expr, data = filter(data_perCell, Rep == '20230802'))
  # t.test(Area_HP1_foci_DNA_tot ~ HP1_expr, data = filter(data_perCell, Rep == '20230920'))
  # 
  # data_perCell %>% select(Rep, HP1_expr, Pos, Cell, starts_with('Area')) %>%
  #   mutate(Perc_HP1_foci_noDNA = Area_HP1_foci_noDNA_tot/Area,
  #          Perc_HP1_foci_DNA = Area_HP1_foci_DNA_tot/Area) %>%
  #   select(Rep, HP1_expr, Pos, Cell, starts_with('Perc')) %>%
  #   pivot_longer(cols = c(Perc_HP1_foci_noDNA, Perc_HP1_foci_DNA),
  #                names_to = 'Focus_type', values_to = 'Perc_area') %>%
  #   mutate(Focus_type = case_when(grepl('noDNA', Focus_type) ~ 'DNA-low',
  #                                 T ~ 'DNA-high')) %>%
  #   ggplot(aes(x = HP1_expr, y = Perc_area, fill = Focus_type)) +
  #   geom_violin() +
  #   geom_boxplot(width = 0.2) +
  #   facet_wrap(Focus_type~Rep, scales = 'free_y') +
  #   theme_bw()
  # t.test(Area_HP1_foci_DNA_tot/Area ~ HP1_expr, data = filter(data_perCell, Rep == '20230802'))
  # t.test(Area_HP1_foci_DNA_tot/Area ~ HP1_expr, data = filter(data_perCell, Rep == '20230920'))
  # t.test(Area_HP1_foci_noDNA_tot/Area ~ HP1_expr, data = filter(data_perCell, Rep == '20230802'))
  # t.test(Area_HP1_foci_noDNA_tot/Area ~ HP1_expr, data = filter(data_perCell, Rep == '20230920'))
  # 
  # data_perCell %>% select(Rep, HP1_expr, Pos, Cell, starts_with('Area')) %>%
  #   mutate(Perc_HP1_foci_noDNA = Area_HP1_foci_noDNA_tot/Area) %>%
  #   select(Rep, HP1_expr, Pos, Cell, starts_with('Perc')) %>%
  #   ggplot(aes(x = HP1_expr, y = Perc_HP1_foci_noDNA)) +
  #   geom_violin() +
  #   geom_boxplot(width = 0.2) +
  #   facet_grid(~Rep, scales = 'free_y') +
  #   coord_cartesian(ylim = c(0, 0.06)) +
  #   theme_bw()
  # 
  # data_perCell %>% select(Rep, HP1_expr, Pos, Cell, starts_with('Area_DNA_foci')) %>%
  #   pivot_longer(cols = c(Area_DNA_foci_noHP1_tot, Area_DNA_foci_HP1_tot),
  #                names_to = 'Focus_type', values_to = 'Area') %>%
  #   mutate(Focus_type = case_when(grepl('noHP1', Focus_type) ~ 'HP1-low',
  #                                 T ~ 'HP1-high')) %>%
  #   ggplot(aes(x = HP1_expr, y = Area, fill = Focus_type)) +
  #   geom_violin() +
  #   geom_boxplot(width = 0.2) +
  #   facet_wrap(Focus_type~Rep, scales = 'free_y') +
  #   theme_bw()
  # t.test(Area_DNA_foci_HP1_tot ~ HP1_expr, data = filter(data_perCell, Rep == '20230802'))
  # t.test(Area_DNA_foci_HP1_tot ~ HP1_expr, data = filter(data_perCell, Rep == '20230920'))
  # t.test(Area_DNA_foci_noHP1_tot ~ HP1_expr, data = filter(data_perCell, Rep == '20230802'))
  # t.test(Area_DNA_foci_noHP1_tot ~ HP1_expr, data = filter(data_perCell, Rep == '20230920'))
}
