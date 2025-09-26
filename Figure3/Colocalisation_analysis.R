library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(stringr)
library(ggpubr)

setwd('/home/sasha/Science/PhD/Laue_lab/Microscopy/STED/20230728-0801_HP1b_DNA_H3K9me3/Colocalisation_analysis')

files = list.files('Coloc')
for (f in files) {
  tmp = read.table(file.path('Coloc', f), header = F)
  colnames(tmp) = c('Focus_id', 'R')
  tmp %>% mutate(Molecule1 = str_split(f, '_')[[1]][1],
                 Molecule2 = str_split(f, '_')[[1]][2],
                 Cell = str_split(f, '_')[[1]][3],
                 Focus_id = Focus_id + 1) -> tmp
  if (f == files[1]) {
    correlations = tmp
  } else {
    correlations = rbind(correlations, tmp)
  }
}


#Just everything, without pre-filtering
ggplot(correlations, aes(x = Molecule1, y = R)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.2) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw() +
  theme(aspect.ratio = 1.5)


correlations$Focus_type = c('DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'DNA_H3K9me3',
                            'HP1_H3K9me3_large',
                            'HP1_H3K9me3_large',
                            'HP1_small',
                            'HP1_H3K9me3_small',
                            'HP1_H3K9me3_large',
                            'HP1_H3K9me3_large',
                            'HP1_H3K9me3_large',
                            'H3K9me3_large',
                            'HP1_H3K9me3_large',
                            'HP1_H3K9me3_small',
                            'HP1_small',
                            'HP1_H3K9me3_small',
                            'HP1_H3K9me3_large',
                            'H3K9me3_large',
                            'H3K9me3_large',
                            'HP1_H3K9me3_large',
                            'HP1_H3K9me3_large',
                            'HP1_large',
                            'HP1_H3K9me3_large',
                            'HP1_H3K9me3_large',
                            'HP1_H3K9me3_large',
                            'HP1_H3K9me3_large',
                            'HP1_H3K9me3_large',
                            'HP1_H3K9me3_large',
                            'HP1_H3K9me3_large',
                            'HP1_H3K9me3_large',
                            'H3K9me3_large',
                            'H3K9me3_large',
                            'H3K9me3_large',
                            'HP1_H3K9me3_small',
                            'H3K9me3_large',
                            'H3K9me3_large',
                            'H3K9me3_large',
                            'H3K9me3_large',
                            'HP1_H3K9me3_small',
                            'H3K9me3_large',
                            'HP1_H3K9me3_large',
                            'HP1_H3K9me3_small')


ggplot(correlations, aes(x = Molecule1, y = R)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.2, aes(colour = Focus_type)) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw() +
  theme(aspect.ratio = 1.5)

filter(correlations, Focus_type == 'DNA_H3K9me3' | Focus_type == 'HP1_H3K9me3_large') %>%
ggplot(aes(x = Molecule1, y = R)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.2) +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw() +
  theme(aspect.ratio = 1.5)

filter(correlations, Focus_type == 'DNA_H3K9me3' | Focus_type == 'HP1_H3K9me3_large') %>% #group_by(Molecule1) %>% summarise(N = n())
  ggplot(aes(x = R, y = ..density..*0.1, fill = Molecule1)) +
  geom_histogram(breaks = seq(0, 1, by = 0.1), colour = 'gray20') +
  facet_grid(Molecule1~.) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c('DNA' = 'cyan', 'HP1' = 'red')) +
  labs(x = 'Pearson\'s R', y = '') +
  theme_bw() +
  theme(aspect.ratio = 0.7, legend.position = 'none')
ggsave('R_histogram_onlyOverlappingLarge.svg', width = 3, height = 4)

filter(correlations, Focus_type == 'DNA_H3K9me3' | Focus_type == 'HP1_H3K9me3_large' | Focus_type == 'HP1_H3K9me3_small') %>% #group_by(Focus_type) %>% summarise(N = n())
  ggplot(aes(x = R, y = ..density..*0.1)) +
  geom_histogram(breaks = seq(0, 1, by = 0.1), colour = 'gray20', fill = 'gray80') +
  facet_grid(Focus_type~.) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = 'Pearson\'s R', y = '') +
  theme_bw() +
  theme(aspect.ratio = 0.7)
ggsave('R_histogram_onlyOverlappingLargeAndSmall.png', width = 3, height = 6)
