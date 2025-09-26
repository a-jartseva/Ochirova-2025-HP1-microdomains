library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(stringr)
library(ggpubr)
library(viridis)
library(svglite)

rm(list = ls())
setwd('/home/sasha/Science/PhD/Laue_lab/Microscopy/Hi-D/20220701_SoRa_HP1b-JF549_SiR-DNA_200ms')
set.seed(100500)

load = T
readFiles = F
numframes = 150

#### Read in data ####
{
# if (readFiles) {
#   conds = unique(sapply(list.files('Hi-D_results', pattern = 'img.txt'), function(x) str_split(x, pattern = '_')[[1]][1]))
#   nucs = unique(sapply(list.files('Hi-D_results', pattern = 'img.txt'),
#                        function(x) paste0(str_split(x, pattern = '_')[[1]][2], '_', str_split(x, pattern = '_')[[1]][3])))
#   molecules = unique(sapply(list.files('Hi-D_results', pattern = 'img.txt'), function(x) str_split(x, pattern = '_')[[1]][4]))
#   if (!('data' %in% ls())) {
#     data = list()
#   }
#   for (c in conds) {
#     for (n in nucs) {
#       data[[paste0(c, '_', str_replace(n, '_', '.'))]] = list()
#       for (m in molecules) {
#         print(paste(c, n, m))
#         files = list.files('Hi-D_results', pattern = paste(c, n, m, sep = '_'))
#         variables = str_remove(sapply(files,
#                                       function(x) str_split(x, pattern = '_')[[1]][length(str_split(x, pattern = '_')[[1]])]),
#                                '.txt')
#         for (i in 1:length(files)) {
#           if (variables[i] != 'img') {
#             assign(variables[i], as.matrix(read.table(file.path('Hi-D_results', files[i]))))
#           }
#         }
#         img = array(unlist(read.table(file.path('Hi-D_results', files[which(variables == 'img')]))),
#                     dim = c(nrow(A), ncol(A), numframes))
#         img_avg = apply(img, c(1,2), mean)
#         variables[variables=='img'] = 'img_avg'
#         data[[paste0(c, '_', str_replace(n, '_', '.'))]][[m]] = setNames(lapply(variables, get), variables)
#         # assign(paste(c, n, m, sep = '_'),
#         #        setNames(lapply(variables, get), variables))
#       }
#     }
#   }
#   rm(list = c('A', 'D', 'img_avg', 'mask', 'model', 'Se', 'V'))
#   saveRDS(data, 'Hi-D_data_list2.rds')
# }
}

#Load
if (load) {
  data = readRDS('Hi-D_data_list.rds')
  data2 = readRDS('/home/sasha/Science/PhD/Laue_lab/Microscopy/Hi-D/20220816_SoRa_HP1b-JF549_SiR-DNA_200ms/Hi-D_data_list.rds')
  data2[!(grepl('highLas', names(data2)) | grepl('vlowLas', names(data2)))] -> data2
  str_remove(names(data2), '_lowLas') -> names(data2)
}

#Make df
lapply(data, function(x) {
  lapply(names(x), function(m) {
    if (m == 'DNA') {
      return(data.frame(DNA_model = factor(c('None', 'D', 'DA', '0', 'D', 'DA')[match(c(x[[m]]$model), 0:5)],
                                           levels = c('None', '0', 'D', 'DA')),
                        DNA_D = c(x[[m]]$D),
                        DNA_A = c(x[[m]]$A),
                        DNA_intensity = c(x[[m]]$img_avg),
                        DNA_intensity_rel = c(x[[m]]$img_avg)/mean(c(x[[m]]$img_avg)[c(x[[m]]$mask != 0)])))
    } else {
      return(data.frame(HP1_model = factor(c('None', 'D', 'DA', '0', 'D', 'DA')[match(c(x[[m]]$model), 0:5)],
                                           levels = c('None', '0', 'D', 'DA')),
                        HP1_D = c(x[[m]]$D),
                        HP1_A = c(x[[m]]$A),
                        HP1_intensity = c(x[[m]]$img_avg),
                        HP1_intensity_rel = c(x[[m]]$img_avg)/mean(c(x[[m]]$img_avg)[c(x[[m]]$mask != 0)])))
    }
  }) %>% bind_cols() %>% mutate(Px.ID = 1:nrow(.)) %>% filter(DNA_D != 0 | HP1_D != 0) -> result
  return(result)
}) %>% bind_rows(.id = 'Dataset') %>% separate(Dataset, into = c('Condition', 'Nuc'), sep = '_') %>%
  mutate(Condition = factor(Condition, levels = c('2ilif', '24h')), Dataset = '1') -> D_df
lapply(data2, function(x) {
  lapply(names(x), function(m) {
    if (m == 'DNA') {
      return(data.frame(DNA_model = factor(c('None', 'D', 'DA', '0', 'D', 'DA')[match(c(x[[m]]$model), 0:5)],
                                           levels = c('None', '0', 'D', 'DA')),
                        DNA_D = c(x[[m]]$D),
                        DNA_A = c(x[[m]]$A),
                        DNA_intensity = c(x[[m]]$img_avg),
                        DNA_intensity_rel = c(x[[m]]$img_avg)/mean(c(x[[m]]$img_avg)[c(x[[m]]$mask != 0)])))
    } else {
      return(data.frame(HP1_model = factor(c('None', 'D', 'DA', '0', 'D', 'DA')[match(c(x[[m]]$model), 0:5)],
                                           levels = c('None', '0', 'D', 'DA')),
                        HP1_D = c(x[[m]]$D),
                        HP1_A = c(x[[m]]$A),
                        HP1_intensity = c(x[[m]]$img_avg),
                        HP1_intensity_rel = c(x[[m]]$img_avg)/mean(c(x[[m]]$img_avg)[c(x[[m]]$mask != 0)])))
    }
  }) %>% bind_cols() %>% mutate(Px.ID = 1:nrow(.)) %>% filter(DNA_D != 0 | HP1_D != 0) -> result
  return(result)
}) %>% bind_rows(.id = 'Dataset') %>% separate(Dataset, into = c('Condition', 'Nuc'), sep = '_') %>%
  mutate(Condition = factor(Condition, levels = c('Fixed', '2ilif')), Dataset = '2') %>% rbind(D_df) -> D_df

#### Analysis ####
#Comparing all compartments
{
  #DNA diffusion, D model - similar result to DA
  # set.seed(100500)
  # D_df %>% pivot_longer(c(DNA_model, HP1_model, DNA_D, HP1_D),
  #                       names_sep = '_',
  #                       names_to = c('Molecule', '.value')) %>%
  #   filter(Molecule == 'DNA', Condition == '2ilif', D != 0, model == 'D') %>%
  #   mutate(Compartment = case_when(DNA_intensity_rel > 1.5 & HP1_intensity_rel > 1.4 ~ 'DNA&HP1-high',
  #                                  DNA_intensity_rel > 1.5 & HP1_intensity_rel < 1.1 ~ 'DNA-high',
  #                                  DNA_intensity_rel < 1.1 & HP1_intensity_rel > 1.4 ~ 'HP1-high',
  #                                  DNA_intensity_rel < 1.1 & HP1_intensity_rel < 1.1 ~ 'Euchromatin',
  #                                  TRUE ~ "Unclear")) %>% #NB! Changed DNA low threshold to 1.1 to have more HP1-only foci
  #   filter(Compartment != 'Unclear') %>%
  #   mutate(Compartment = factor(Compartment, levels = c('Euchromatin', 'HP1-high', 'DNA-high', 'DNA&HP1-high'))) %>%
  #   group_by(Dataset, Condition, Nuc, Molecule, Compartment) %>% #summarise(num = n()) %>% View()
  #   slice_sample(n = 100) %>% #ungroup() %>% group_by(Compartment) %>% summarise(num = n())
  #   ggplot(aes(x = D, colour = Compartment)) +#, group = interaction(Compartment, Dataset))) +
  #   geom_density() +
  #   facet_grid(Dataset~.) +
  #   scale_colour_manual(values = c('Euchromatin' = 'gray50',
  #                                  'HP1-high' = 'red',
  #                                  'DNA-high' = 'cyan',
  #                                  'DNA&HP1-high' = 'purple')) +
  #   #coord_cartesian(xlim = c(0, 0.0008)) +
  #   labs(x = 'D, um2/s', y = '', colour = '', title = 'DNA diffusion') +
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
  #         legend.background = element_rect(fill = NA), aspect.ratio = 0.5)
  
  #DNA diffusion, DA model
  set.seed(100500)
  D_df %>% pivot_longer(c(DNA_model, HP1_model, DNA_D, HP1_D),
                        names_sep = '_',
                        names_to = c('Molecule', '.value')) %>%
    filter(Molecule == 'DNA', Condition == '2ilif', D != 0, model == 'DA') %>%
    mutate(Compartment = case_when(DNA_intensity_rel > 1.5 & HP1_intensity_rel > 1.4 ~ 'DNA&HP1-high',
                                   DNA_intensity_rel > 1.5 & HP1_intensity_rel < 1.1 ~ 'DNA-high',
                                   DNA_intensity_rel < 1.1 & HP1_intensity_rel > 1.4 ~ 'HP1-high',
                                   DNA_intensity_rel < 1.1 & HP1_intensity_rel < 1.1 ~ 'Euchromatin',
                                   TRUE ~ "Unclear")) %>% #NB! Changed DNA low threshold to 1.1 to have more HP1-only foci
    filter(Compartment != 'Unclear') %>%
    mutate(Compartment = factor(Compartment, levels = c('Euchromatin', 'HP1-high', 'DNA-high', 'DNA&HP1-high'))) %>%
    group_by(Dataset, Condition, Nuc, Molecule, Compartment) %>% #summarise(num = n()) %>% View()
    slice_sample(n = 70) %>% #ungroup() %>% group_by(Compartment) %>% summarise(num = n())
    ggplot(aes(x = D, colour = Compartment)) +#, group = interaction(Compartment, Dataset))) +
    geom_density() +
    #facet_grid(Dataset~.) +
    scale_colour_manual(values = c('Euchromatin' = 'gray50',
                                   'HP1-high' = 'red',
                                   'DNA-high' = 'cyan',
                                   'DNA&HP1-high' = 'purple')) +
    coord_cartesian(xlim = c(0, 0.0008)) +
    labs(x = 'D, um2/s', y = '', colour = '', title = 'DNA diffusion') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.background = element_rect(fill = NA), aspect.ratio = 0.5)
  ggsave('/home/sasha/Science/PhD/Laue_lab/Heterochromatin_paper/Figures/Hi-D_DNA_Ddistrib_DAModel_subsampled_repsTogether.svg',
         width = 5, height = 2.5)
  
  #HP1 D - nice trend, but faster euchromatin not reproducible
  # set.seed(100500)
  # D_df %>% pivot_longer(c(DNA_model, HP1_model, DNA_D, HP1_D),
  #                       names_sep = '_',
  #                       names_to = c('Molecule', '.value')) %>%
  #   filter(Molecule == 'HP1', Condition == '2ilif', D != 0, model == 'DA') %>%
  #   mutate(Compartment = case_when(DNA_intensity_rel > 1.5 & HP1_intensity_rel > 1.4 ~ 'DNA&HP1-high',
  #                                  DNA_intensity_rel > 1.5 & HP1_intensity_rel < 1.1 ~ 'DNA-high',
  #                                  DNA_intensity_rel < 1.1 & HP1_intensity_rel > 1.4 ~ 'HP1-high',
  #                                  DNA_intensity_rel < 1.1 & HP1_intensity_rel < 1.1 ~ 'Euchromatin',
  #                                  TRUE ~ "Unclear")) %>% #NB! Changed DNA low threshold to 1.1 to have more HP1-only foci
  #   filter(Compartment != 'Unclear') %>%
  #   mutate(Compartment = factor(Compartment, levels = c('Euchromatin', 'HP1-high', 'DNA-high', 'DNA&HP1-high'))) %>%
  #   group_by(Dataset, Condition, Nuc, Molecule, Compartment) %>% #summarise(num = n()) %>% View()
  #   slice_sample(n = 70) %>% #ungroup() %>% group_by(Compartment) %>% summarise(num = n())
  #   ggplot(aes(x = D, colour = Compartment)) +
  #   geom_density() +
  #   #facet_grid(Dataset~.) +
  #   scale_colour_manual(values = c('Euchromatin' = 'gray50',
  #                                  'HP1-high' = 'red',
  #                                  'DNA-high' = 'cyan',
  #                                  'DNA&HP1-high' = 'purple')) +
  #   coord_cartesian(xlim = c(0, 0.001)) +
  #   labs(x = 'D, um2/s', y = '', colour = '', title = 'HP1 diffusion') +
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
  #         legend.background = element_rect(fill = NA), aspect.ratio = 0.5)
  # ggsave('/home/sasha/Science/PhD/Laue_lab/Heterochromatin_paper/Figures/Hi-D_HP1_Ddistrib_DAModel_subsampled_repsTogether.svg',
  #        width = 5, height = 2.5)
}

#Same - including fixed control
{
  #DNA diffusion, D model - similar result to DA
  # set.seed(100500)
  # D_df %>% pivot_longer(c(DNA_model, HP1_model, DNA_D, HP1_D),
  #                       names_sep = '_',
  #                       names_to = c('Molecule', '.value')) %>%
  #   filter(Molecule == 'DNA', Condition != '24h', D != 0, model == 'D') %>%
  #   mutate(Compartment = case_when(DNA_intensity_rel > 1.5 & HP1_intensity_rel > 1.4 ~ 'DNA&HP1-high',
  #                                  DNA_intensity_rel > 1.5 & HP1_intensity_rel < 1.1 ~ 'DNA-high',
  #                                  DNA_intensity_rel < 1.1 & HP1_intensity_rel > 1.4 ~ 'HP1-high',
  #                                  DNA_intensity_rel < 1.1 & HP1_intensity_rel < 1.1 ~ 'Euchromatin',
  #                                  TRUE ~ "Unclear")) %>% #NB! Changed DNA low threshold to 1.1 to have more HP1-only foci
  #   filter(Compartment != 'Unclear') %>%
  #   mutate(Compartment = factor(Compartment, levels = c('Euchromatin', 'HP1-high', 'DNA-high', 'DNA&HP1-high'))) %>%
  #   group_by(Dataset, Condition, Nuc, Molecule, Compartment) %>% #summarise(num = n()) %>% View()
  #   slice_sample(n = 100) %>% #ungroup() %>% group_by(Compartment) %>% summarise(num = n())
  #   ggplot(aes(x = D, colour = Compartment)) +#, group = interaction(Compartment, Dataset))) +
  #   geom_density() +
  #   facet_grid(Dataset~Condition) +
  #   scale_colour_manual(values = c('Euchromatin' = 'gray50',
  #                                  'HP1-high' = 'red',
  #                                  'DNA-high' = 'cyan',
  #                                  'DNA&HP1-high' = 'purple')) +
  #   coord_cartesian(xlim = c(0, 0.0008)) +
  #   labs(x = 'D, um2/s', y = '', colour = '', title = 'DNA diffusion') +
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
  #         legend.background = element_rect(fill = NA), aspect.ratio = 0.5)
  
  #DNA diffusion, DA model
  set.seed(100500)
  D_df %>% pivot_longer(c(DNA_model, HP1_model, DNA_D, HP1_D),
                        names_sep = '_',
                        names_to = c('Molecule', '.value')) %>%
    filter(Molecule == 'DNA', Condition != '24h', D != 0, model == 'DA') %>%
    mutate(Compartment = case_when(DNA_intensity_rel > 1.5 & HP1_intensity_rel > 1.4 ~ 'DNA&HP1-high',
                                   DNA_intensity_rel > 1.5 & HP1_intensity_rel < 1.1 ~ 'DNA-high',
                                   DNA_intensity_rel < 1.1 & HP1_intensity_rel > 1.4 ~ 'HP1-high',
                                   DNA_intensity_rel < 1.1 & HP1_intensity_rel < 1.1 ~ 'Euchromatin',
                                   TRUE ~ "Unclear")) %>% #NB! Changed DNA low threshold to 1.1 to have more HP1-only foci
    filter(Compartment != 'Unclear') %>%
    mutate(Compartment = factor(Compartment, levels = c('Euchromatin', 'HP1-high', 'DNA-high', 'DNA&HP1-high'))) %>%
    group_by(Dataset, Condition, Nuc, Molecule, Compartment) %>% #summarise(num = n()) %>% View()
    slice_sample(n = 70) %>% #ungroup() %>% group_by(Compartment) %>% summarise(num = n())
    ggplot(aes(x = D, colour = Compartment)) +#, group = interaction(Compartment, Dataset))) +
    geom_density() +
    facet_grid(Dataset~Condition) +
    scale_colour_manual(values = c('Euchromatin' = 'gray50',
                                   'HP1-high' = 'red',
                                   'DNA-high' = 'cyan',
                                   'DNA&HP1-high' = 'purple')) +
    coord_cartesian(xlim = c(0, 0.0008)) +
    labs(x = 'D, um2/s', y = '', colour = '', title = 'DNA diffusion') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.background = element_rect(fill = NA), aspect.ratio = 0.5)
  ggsave('/home/sasha/Science/PhD/Laue_lab/Heterochromatin_paper/Figures/Hi-D_DNA_Ddistrib_DAModel_subsampled_withFixed.svg',
         width = 7.5, height = 5)
}

#Comparing DNA density in compartments
{
  set.seed(100500)
  D_df %>% pivot_longer(c(DNA_model, HP1_model, DNA_D, HP1_D),
                        names_sep = '_',
                        names_to = c('Molecule', '.value')) %>%
    filter(DNA_intensity_rel > 1.5, Molecule == 'DNA', Condition == '2ilif', D != 0, model == 'DA') %>%
    mutate(Compartment = case_when(DNA_intensity_rel > 1.5 & HP1_intensity_rel > 1.4 ~ 'DNA&HP1-high',
                                   DNA_intensity_rel > 1.5 & HP1_intensity_rel < 1.1 ~ 'DNA-high',
                                   TRUE ~ "Unclear")) %>%
    filter(Compartment != 'Unclear') %>%
    mutate(Compartment = factor(Compartment, levels = c('Euchromatin', 'HP1-high', 'DNA-high', 'DNA&HP1-high'))) %>%
    group_by(Condition, Nuc, Molecule, Compartment) %>% #summarise(num = n()) %>% View()
    slice_sample(n = 70) %>% #ungroup() %>% group_by(Compartment) %>% summarise(num = n())
    ggplot(aes(x = DNA_intensity_rel, colour = Compartment)) +
    geom_density() +
    #facet_grid(Dataset~.) +
    scale_colour_manual(values = c('Euchromatin' = 'gray50',
                                   'HP1-high' = 'red',
                                   'DNA-high' = 'cyan',
                                   'DNA&HP1-high' = 'purple')) +
    #lims(x = c(0, 0.0014)) +
    labs(x = 'DNA density (relative to mean)', y = '', colour = '',
         title = 'DNA density in heterochromatin', subtitle = 'Same pixels analysed') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          legend.background = element_rect(fill = NA), aspect.ratio = 0.5)
  ggsave('/home/sasha/Science/PhD/Laue_lab/Heterochromatin_paper/Figures/Hi-D_intensity_DAModel_subsampled_repsTogether.svg',
         width = 5, height = 2.5)
}