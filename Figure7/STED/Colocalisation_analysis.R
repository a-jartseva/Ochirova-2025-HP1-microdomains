library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(stringr)
library(ggpubr)
library(car)
library(dunn.test)

setwd('/home/sasha/Science/PhD/Laue_lab/Microscopy/STED/HP1b_DNA_differentiation_Devina/Colocalisation_analysis')

files = list.files('Coloc')
for (f in files) {
  tmp = read.table(file.path('Coloc', f), header = F)
  colnames(tmp) = c('Focus_id', 'R')
  tmp %>% mutate(Condition = str_split(f, '_')[[1]][3],
                 Cell = str_split(f, '_')[[1]][4],
                 Focus_id = Focus_id + 1) -> tmp
  if (f == files[1]) {
    correlations = tmp
  } else {
    correlations = rbind(correlations, tmp)
  }
}
mutate(correlations, Condition = factor(Condition, c('2iL', '24h', '48h'))) -> correlations


#Just everything, without pre-filtering
ggplot(correlations, aes(x = Condition, y = R, fill = Condition)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  geom_jitter(width = 0.2) +
  coord_cartesian(ylim = c(0, NA)) +
  scale_fill_brewer(palette = 'Set2') +
  theme_bw() +
  theme(aspect.ratio = 1)
ggsave('R_violins.svg', width = 4, height = 3)

ggplot(correlations, aes(x = R, y = ..density..*0.1, fill = Condition)) +
  geom_histogram(breaks = seq(0, 1, by = 0.1), colour = 'gray20') +
  facet_grid(Condition~.) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = 'Set2') +
  labs(x = 'Pearson\'s R', y = '') +
  theme_bw() +
  theme(aspect.ratio = 0.7, legend.position = 'none')
ggsave('R_histograms.svg', width = 3, height = 5)

#Stats
group_by(correlations, Condition) %>% summarise(n())
shapiro.test(filter(correlations, Condition == '2iL')$R) #0.0145
shapiro.test(filter(correlations, Condition == '24h')$R) #0.08211
shapiro.test(filter(correlations, Condition == '48h')$R) #very sig
#=> non-parametric
leveneTest(R~Condition, data = correlations)#var equal
kruskal.test(R~Condition, data = correlations) #p = 0.004
dunn.test(correlations$R, correlations$Condition) #2iL-24h 0.067; 2iL-48h 0.1; 24h-48h 0.0005

