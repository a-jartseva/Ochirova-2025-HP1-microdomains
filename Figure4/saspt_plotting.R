library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)

rm(list = ls())
setwd('/home/sasha/Science/PhD/Laue_lab/Microscopy/HP1b_live_blinkingDyes')
folder1 = '20221018_Halo-JF639b_Hoechst-SPY505/New_analysis/SA-SPT'
folder2 = '20221108_Halo-JF639b_Hoechst-SPY505/New_analysis/SA-SPT'
folder3 = '20230620_Halo-JF639b_DNA-SPY505_fixed_immobDye/New_analysis/SA-SPT'


#Total live vs immob
{
  occupations_all1 = read.table(paste0(folder1, '/AllData_occupations_all.csv'), sep = ',', header = T)[-1]
  occupations_all2 = read.table(paste0(folder2, '/AllData_occupations_all.csv'), sep = ',', header = T)[-1]
  live_occupations_all3 = read.table(paste0(folder3, '/Live_occupations_all.csv'), sep = ',', header = T)[-1]
  immob_occupations_all3 = read.table(paste0(folder3, '/ImmobDye_occupations_all.csv'), sep = ',', header = T)[-1]
  immob_occupations_bs3 = read.table(paste0(folder3, '/ImmobDye_occupations_bs.csv'), sep = ',', header = T)[-1]
  immob_occupations_bs3$rep = rep(1:100, each = 3600)
  
  group_by(occupations_all1, diff_coef) %>%
    summarise(mean_posterior_occupation = sum(mean_posterior_occupation)) %>%
    ggplot(aes(x = log10(diff_coef), y = mean_posterior_occupation)) +
    geom_line(colour = 'black') +
    geom_line(data = (group_by(occupations_all2, diff_coef) %>%
                        summarise(mean_posterior_occupation = sum(mean_posterior_occupation))),
              colour = 'black') +
    geom_line(data = (group_by(live_occupations_all3, diff_coef) %>%
                        summarise(mean_posterior_occupation = sum(mean_posterior_occupation))),
              colour = 'black') +
    geom_line(data = (group_by(immob_occupations_all3, diff_coef) %>%
                        summarise(mean_posterior_occupation = sum(mean_posterior_occupation))),
              colour = muted('red')) +
    scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
    coord_cartesian(ylim = c(0, 0.05)) + #coord_cartesian(ylim = c(0, 0.1)) +
    labs(title = 'Subsampled, N = 100,000') +
    theme_bw() +
    theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
  ggsave('SA-SPT_results/AllData_Immob_RBME_margDoccupations2.png', width = 4, height = 2.5)
  
  rbind(mutate(occupations_all1, bioRep = 'Rep1', type = 'Live'), 
        mutate(occupations_all2, bioRep = 'Rep2', type = 'Live'), 
        mutate(live_occupations_all3, bioRep = 'Rep3', type = 'Live'), 
        mutate(immob_occupations_all3, bioRep = 'Rep3', type = 'ImmobDye')) %>%
    group_by(type, bioRep, diff_coef) %>%
    summarise(mean_posterior_occupation = sum(mean_posterior_occupation)) %>%
    group_by(type, bioRep, population = case_when(diff_coef < 0.1 ~ 'Bound', T ~ 'Unbound')) %>%
    summarise(fraction = sum(mean_posterior_occupation)) %>%
    mutate(fraction = fraction/sum(fraction)) %>% 
    mutate(type = factor(type, levels = c('Live', 'ImmobDye'))) %>%
    ggplot(aes(x = interaction(bioRep, type),
               y = fraction,
               fill = type,
               alpha = factor(population, levels = c('Unbound', 'Bound')),
               group = interaction(bioRep, type))) +
    geom_col(position = 'stack', colour = 'gray20', width = 0.8) +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c('black', muted('red'))) +
    scale_alpha_manual(values = c(0.1, 1), guide = F) +
    labs(x = '', y = 'Proportion bound', fill = 'Type') +
    theme_bw() +
    theme(aspect.ratio = 1)
  ggsave('SA-SPT_results/occupations_bars_all_liveImmob.png', width = 4, height = 2.5)
}

#All data - subsampled
{
  mutate(eu_occupations_bs1, compartment = 'Eu', type = 'real', bioRep = 'Rep1') %>%
    rbind(mutate(HP1_occupations_bs1, compartment = 'HP1 foci', type = 'real', bioRep = 'Rep1')) %>%
    rbind(mutate(DNA_occupations_bs1, compartment = 'DNA foci', type = 'real', bioRep = 'Rep1')) %>%
    #rbind(mutate(DNA_sub_occupations_bs1, compartment = 'DNA foci, sub', type = 'real', bioRep = 'Rep1')) %>%
    rbind(mutate(DNA_HP1_occupations_bs1, compartment = 'HP1&DNA foci', type = 'real', bioRep = 'Rep1')) %>%
    rbind(mutate(HP1_ps_occupations_bs1, compartment = 'HP1 foci', type = 'pseudo', bioRep = 'Rep1')) %>%
    rbind(mutate(DNA_ps_occupations_bs1, compartment = 'DNA foci', type = 'pseudo', bioRep = 'Rep1')) %>%
    rbind(mutate(DNA_HP1_ps_occupations_bs1, compartment = 'HP1&DNA foci', type = 'pseudo', bioRep = 'Rep1')) %>%
    rbind(mutate(eu_occupations_bs2, compartment = 'Eu', type = 'real', bioRep = 'Rep2')) %>%
    rbind(mutate(HP1_occupations_bs2, compartment = 'HP1 foci', type = 'real', bioRep = 'Rep2')) %>%
    rbind(mutate(DNA_occupations_bs2, compartment = 'DNA foci', type = 'real', bioRep = 'Rep2')) %>%
    #rbind(mutate(DNA_sub_occupations_bs2, compartment = 'DNA foci, sub', type = 'real', bioRep = 'Rep2')) %>%
    rbind(mutate(DNA_HP1_occupations_bs2, compartment = 'HP1&DNA foci', type = 'real', bioRep = 'Rep2')) %>%
    rbind(mutate(HP1_ps_occupations_bs2, compartment = 'HP1 foci', type = 'pseudo', bioRep = 'Rep2')) %>%
    rbind(mutate(DNA_ps_occupations_bs2, compartment = 'DNA foci', type = 'pseudo', bioRep = 'Rep2')) %>%
    rbind(mutate(DNA_HP1_ps_occupations_bs2, compartment = 'HP1&DNA foci', type = 'pseudo', bioRep = 'Rep2'))-> allData
  
  group_by(allData, compartment, type, bioRep, rep, diff_coef) %>%
    summarise(mean_posterior_occupation = sum(mean_posterior_occupation)) %>%
    group_by(compartment, bioRep, rep, type) %>%
    mutate(mean_posterior_occupation = mean_posterior_occupation/sum(mean_posterior_occupation)) %>%
    group_by(type, compartment, bioRep, diff_coef) %>%
    summarise(CI1 = quantile(mean_posterior_occupation, 0.025),
              CI2 = quantile(mean_posterior_occupation, 0.975),
              mean_posterior_occupation = mean(mean_posterior_occupation)) %>%
    mutate(compartment = factor(compartment, levels = c('Eu', 'HP1 foci', 'DNA foci', 'HP1&DNA foci')),
           type = factor(type, levels = c('real', 'pseudo'))) %>%
    ggplot(aes(x = log10(diff_coef), y = mean_posterior_occupation,
               colour = interaction(compartment, type), fill = interaction(compartment, type),
               linetype = type, group = interaction(compartment, type, bioRep))) +
    geom_line() +
    geom_ribbon(aes(ymin = CI1, ymax = CI2), alpha = 0.2, colour = NA) +
    facet_grid(compartment~.) +
    scale_x_continuous(labels = math_format(10^.x), breaks = seq(-5, 5, by = 1)) +
    scale_fill_manual(values = c('gray', 'red', 'cyan', 'purple', 'tomato4', 'blue', 'purple4'),
                      labels = c('Rest of euchromatin', 'HP1-high', 'DNA-high',   'HP1&DNA-high',
                                 'HP1 pseudo', 'DNA pseudo', 'HP1&DNA pseudo')) +
    scale_colour_manual(values = c('gray', 'red', 'cyan', 'purple', 'tomato4', 'blue', 'purple4'),
                        labels = c('Rest of euchromatin', 'HP1-high', 'DNA-high',  'HP1&DNA-high',
                                   'HP1 pseudo', 'DNA pseudo', 'HP1&DNA pseudo')) +
    scale_linetype_manual(values = c('solid', 'longdash'), guide = 'none') +
    coord_cartesian(ylim = c(0, 0.05)) +
    labs(x = 'log(D, um2/s)', y = 'Posterior occupation', colour = 'Compartment', fill = 'Compartment') +
    theme_bw() +
    theme(aspect.ratio = 0.5)
  ggsave('SA-SPT_results/posterior_line_facets.png', width = 5, height = 5)
  
  group_by(allData, compartment, type, bioRep, rep, diff_coef) %>%
    summarise(mean_posterior_occupation = sum(mean_posterior_occupation)) %>%
    group_by(compartment, bioRep, rep, type) %>%
    mutate(mean_posterior_occupation = mean_posterior_occupation/sum(mean_posterior_occupation)) %>%
    group_by(compartment, bioRep, type, rep, population = case_when(log10(diff_coef) < -1 ~ 'Bound', T ~ 'Unbound')) %>%
    summarise(mean_posterior_occupation = sum(mean_posterior_occupation)) %>%
    group_by(compartment, bioRep, type, population) %>%
    summarise(CI1 = quantile(mean_posterior_occupation, 0.025),
              CI2 = quantile(mean_posterior_occupation, 0.975),
              mean_posterior_occupation = mean(mean_posterior_occupation)) %>%
    mutate(compartment = factor(compartment, levels = c('Eu', 'HP1 foci', 'DNA foci',
                                                        'HP1&DNA foci')),
           type = factor(type, levels = c('real', 'pseudo')),
           population = factor(population, levels = c('Bound', 'Unbound'))) -> allData_summary
  
  allData_summary %>%
    ggplot(aes(x = interaction(bioRep, type, compartment),
               y = mean_posterior_occupation,
               fill = interaction(type, compartment),
               alpha = factor(population, levels = c('Unbound', 'Bound')),
               group = interaction(bioRep, type, compartment))) +
    geom_col(position = 'stack', colour = 'gray20', width = 0.8) +
    geom_errorbar(data = . %>% filter(population == 'Bound'),
                  aes(ymin = CI1, ymax = CI2), width = 0.3) +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c('gray', 'red', 'tomato4', 'cyan', 'blue', 'purple', 'purple4'),
                      labels = c('Rest of euchromatin', 'HP1-high', 'HP1 pseudo',
                                 'DNA-high', 'DNA pseudo',
                                 'HP1&DNA-high', 'HP1&DNA pseudo')) +
    scale_alpha_manual(values = c(0.1, 1), guide = F) +
    labs(x = 'Compartment', y = 'Proportion bound', fill = 'Compartment') +
    theme_bw() +
    theme(aspect.ratio = 0.5)
  ggsave('SA-SPT_results/occupations_bars2_2.svg', width = 5.5, height = 2.5)
}

#All data - all data
{
  mutate(eu_occupations_all1, compartment = 'Eu', type = 'real', bioRep = 'Rep1') %>%
    rbind(mutate(HP1_occupations_all1, compartment = 'HP1 foci', type = 'real', bioRep = 'Rep1')) %>%
    rbind(mutate(DNA_occupations_all1, compartment = 'DNA foci', type = 'real', bioRep = 'Rep1')) %>%
    #rbind(mutate(DNA_sub_occupations_all1, compartment = 'DNA foci, sub', type = 'real', bioRep = 'Rep1')) %>%
    rbind(mutate(DNA_HP1_occupations_all1, compartment = 'HP1&DNA foci', type = 'real', bioRep = 'Rep1')) %>%
    rbind(mutate(HP1_ps_occupations_all1, compartment = 'HP1 foci', type = 'pseudo', bioRep = 'Rep1')) %>%
    rbind(mutate(DNA_ps_occupations_all1, compartment = 'DNA foci', type = 'pseudo', bioRep = 'Rep1')) %>%
    rbind(mutate(DNA_HP1_ps_occupations_all1, compartment = 'HP1&DNA foci', type = 'pseudo', bioRep = 'Rep1')) %>%
    rbind(mutate(eu_occupations_all2, compartment = 'Eu', type = 'real', bioRep = 'Rep2')) %>%
    rbind(mutate(HP1_occupations_all2, compartment = 'HP1 foci', type = 'real', bioRep = 'Rep2')) %>%
    rbind(mutate(DNA_occupations_all2, compartment = 'DNA foci', type = 'real', bioRep = 'Rep2')) %>%
    #rbind(mutate(DNA_sub_occupations_all2, compartment = 'DNA foci, sub', type = 'real', bioRep = 'Rep2')) %>%
    rbind(mutate(DNA_HP1_occupations_all2, compartment = 'HP1&DNA foci', type = 'real', bioRep = 'Rep2')) %>%
    rbind(mutate(HP1_ps_occupations_all2, compartment = 'HP1 foci', type = 'pseudo', bioRep = 'Rep2')) %>%
    rbind(mutate(DNA_ps_occupations_all2, compartment = 'DNA foci', type = 'pseudo', bioRep = 'Rep2')) %>%
    rbind(mutate(DNA_HP1_ps_occupations_all2, compartment = 'HP1&DNA foci', type = 'pseudo', bioRep = 'Rep2'))-> allData_all
  
  group_by(allData_all, compartment, type, bioRep, diff_coef) %>%
    summarise(mean_posterior_occupation = sum(mean_posterior_occupation)) %>%
    group_by(compartment, bioRep, type) %>%
    mutate(mean_posterior_occupation = mean_posterior_occupation/sum(mean_posterior_occupation)) %>%
    mutate(compartment = factor(compartment, levels = c('Eu', 'HP1 foci', 'DNA foci', 'HP1&DNA foci')),
           type = factor(type, levels = c('real', 'pseudo'))) %>%
    ggplot(aes(x = log10(diff_coef), y = mean_posterior_occupation,
               colour = interaction(compartment, type), fill = interaction(compartment, type),
               linetype = type, group = interaction(compartment, type, bioRep))) +
    geom_line() +
    facet_grid(compartment~.) +
    scale_x_continuous(labels = math_format(10^.x), breaks = seq(-5, 5, by = 1)) +
    scale_fill_manual(values = c('gray', 'red', 'cyan', 'purple', 'tomato4', 'blue', 'purple4'),
                      labels = c('Rest of euchromatin', 'HP1-high', 'DNA-high',   'HP1&DNA-high',
                                 'HP1 pseudo', 'DNA pseudo', 'HP1&DNA pseudo')) +
    scale_colour_manual(values = c('gray', 'red', 'cyan', 'purple', 'tomato4', 'blue', 'purple4'),
                        labels = c('Rest of euchromatin', 'HP1-high', 'DNA-high',  'HP1&DNA-high',
                                   'HP1 pseudo', 'DNA pseudo', 'HP1&DNA pseudo')) +
    scale_linetype_manual(values = c('solid', 'longdash'), guide = 'none') +
    coord_cartesian(ylim = c(0, 0.05)) +
    labs(x = 'log(D, um2/s)', y = 'Posterior occupation', colour = 'Compartment', fill = 'Compartment') +
    theme_bw() +
    theme(aspect.ratio = 0.5)
  ggsave('SA-SPT_results/posterior_line_facets_allData.svg', width = 5, height = 5)
  
  group_by(allData, compartment, type, bioRep, rep, diff_coef) %>%
    summarise(mean_posterior_occupation = sum(mean_posterior_occupation)) %>%
    group_by(compartment, bioRep, rep, type) %>%
    mutate(mean_posterior_occupation = mean_posterior_occupation/sum(mean_posterior_occupation)) %>%
    group_by(compartment, bioRep, type, rep, population = case_when(log10(diff_coef) < -1 ~ 'Bound', T ~ 'Unbound')) %>%
    summarise(mean_posterior_occupation = sum(mean_posterior_occupation)) %>%
    group_by(compartment, bioRep, type, population) %>%
    summarise(CI1 = quantile(mean_posterior_occupation, 0.025),
              CI2 = quantile(mean_posterior_occupation, 0.975),
              mean_posterior_occupation = mean(mean_posterior_occupation)) %>%
    mutate(compartment = factor(compartment, levels = c('Eu', 'HP1 foci', 'DNA foci',
                                                        'HP1&DNA foci')),
           type = factor(type, levels = c('real', 'pseudo')),
           population = factor(population, levels = c('Bound', 'Unbound'))) -> allData_summary
  
  allData_summary %>%
    ggplot(aes(x = interaction(bioRep, type, compartment),
               y = mean_posterior_occupation,
               fill = interaction(type, compartment),
               alpha = factor(population, levels = c('Unbound', 'Bound')),
               group = interaction(bioRep, type, compartment))) +
    geom_col(position = 'stack', colour = 'gray20', width = 0.8) +
    geom_errorbar(data = . %>% filter(population == 'Bound'),
                  aes(ymin = CI1, ymax = CI2), width = 0.3) +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c('gray', 'red', 'tomato4', 'cyan', 'blue', 'purple', 'purple4'),
                      labels = c('Rest of euchromatin', 'HP1-high', 'HP1 pseudo',
                                 'DNA-high', 'DNA pseudo',
                                 'HP1&DNA-high', 'HP1&DNA pseudo')) +
    scale_alpha_manual(values = c(0.1, 1), guide = F) +
    labs(x = 'Compartment', y = 'Proportion bound', fill = 'Compartment') +
    theme_bw() +
    theme(aspect.ratio = 0.5)
  ggsave('SA-SPT_results/occupations_bars2_2_allData.svg', width = 5.5, height = 2.5)
}

#All data + immob
{
  mutate(eu_occupations_bs1, compartment = 'Eu', type = 'real', bioRep = 'Rep1') %>%
    rbind(mutate(HP1_occupations_bs1, compartment = 'HP1 foci', type = 'real', bioRep = 'Rep1')) %>%
    rbind(mutate(DNA_occupations_bs1, compartment = 'DNA foci', type = 'real', bioRep = 'Rep1')) %>%
    #rbind(mutate(DNA_sub_occupations_bs1, compartment = 'DNA foci, sub', type = 'real', bioRep = 'Rep1')) %>%
    rbind(mutate(DNA_HP1_occupations_bs1, compartment = 'HP1&DNA foci', type = 'real', bioRep = 'Rep1')) %>%
    rbind(mutate(HP1_ps_occupations_bs1, compartment = 'HP1 foci', type = 'pseudo', bioRep = 'Rep1')) %>%
    rbind(mutate(DNA_ps_occupations_bs1, compartment = 'DNA foci', type = 'pseudo', bioRep = 'Rep1')) %>%
    rbind(mutate(DNA_HP1_ps_occupations_bs1, compartment = 'HP1&DNA foci', type = 'pseudo', bioRep = 'Rep1')) %>%
    rbind(mutate(eu_occupations_bs2, compartment = 'Eu', type = 'real', bioRep = 'Rep2')) %>%
    rbind(mutate(HP1_occupations_bs2, compartment = 'HP1 foci', type = 'real', bioRep = 'Rep2')) %>%
    rbind(mutate(DNA_occupations_bs2, compartment = 'DNA foci', type = 'real', bioRep = 'Rep2')) %>%
    #rbind(mutate(DNA_sub_occupations_bs2, compartment = 'DNA foci, sub', type = 'real', bioRep = 'Rep2')) %>%
    rbind(mutate(DNA_HP1_occupations_bs2, compartment = 'HP1&DNA foci', type = 'real', bioRep = 'Rep2')) %>%
    rbind(mutate(HP1_ps_occupations_bs2, compartment = 'HP1 foci', type = 'pseudo', bioRep = 'Rep2')) %>%
    rbind(mutate(DNA_ps_occupations_bs2, compartment = 'DNA foci', type = 'pseudo', bioRep = 'Rep2')) %>%
    rbind(mutate(DNA_HP1_ps_occupations_bs2, compartment = 'HP1&DNA foci', type = 'pseudo', bioRep = 'Rep2')) %>%
    rbind(mutate(immob_occupations_bs3, compartment = 'Immob dye', type = 'real', bioRep = 'Rep3')) -> allData
  
  group_by(allData, compartment, type, bioRep, rep, diff_coef) %>%
    summarise(mean_posterior_occupation = sum(mean_posterior_occupation)) %>%
    group_by(compartment, bioRep, rep, type) %>%
    mutate(mean_posterior_occupation = mean_posterior_occupation/sum(mean_posterior_occupation)) %>%
    group_by(compartment, bioRep, type, diff_coef) %>%
    summarise(CI1 = quantile(mean_posterior_occupation, 0.025),
              CI2 = quantile(mean_posterior_occupation, 0.975),
              mean_posterior_occupation = mean(mean_posterior_occupation)) %>%
    mutate(compartment = factor(compartment, levels = c('Eu', 'HP1 foci', 'DNA foci',
                                                        'HP1&DNA foci',
                                                        'Immob dye')),
           type = factor(type, levels = c('real', 'pseudo'))) %>%
    ggplot(aes(x = log10(diff_coef), y = mean_posterior_occupation,
               colour = compartment, linetype = type, fill = compartment,
               group = interaction(compartment, type, bioRep))) +
    geom_ribbon(aes(ymin = CI1, ymax = CI2), alpha = 0.2, colour = NA) +
    geom_line() +
    scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
    #coord_cartesian(ylim = c(0, 0.04)) +
    theme_bw() +
    theme(aspect.ratio = 0.5)
  
  filter(allData, type != 'pseudo') %>%
    group_by(compartment, bioRep, rep, diff_coef) %>%
    summarise(mean_posterior_occupation = sum(mean_posterior_occupation)) %>%
    group_by(compartment, bioRep, rep) %>%
    mutate(mean_posterior_occupation = mean_posterior_occupation/sum(mean_posterior_occupation)) %>%
    group_by(compartment, bioRep, diff_coef) %>%
    summarise(CI1 = quantile(mean_posterior_occupation, 0.025),
              CI2 = quantile(mean_posterior_occupation, 0.975),
              mean_posterior_occupation = mean(mean_posterior_occupation)) %>%
    mutate(compartment = factor(compartment, levels = c('Eu', 'HP1 foci', 'DNA foci',
                                                        'HP1&DNA foci',
                                                        'Immob dye'))) %>%
    ggplot(aes(x = log10(diff_coef), y = mean_posterior_occupation,
               colour = compartment, fill = compartment, group = interaction(compartment, bioRep))) +
    geom_line() +
    geom_ribbon(aes(ymin = CI1, ymax = CI2), alpha = 0.2, colour = NA) +
    scale_x_continuous(labels = math_format(10^.x), breaks = seq(-5, 5, by = 1)) +
    scale_colour_manual(values = c('gray', 'red', 'cyan', 'purple',  'black'), #, 'purple4', 'blue', 'tomato4'),
                        labels = c('Rest of euchromatin', 'HP1-high', 'DNA-high', 
                                   'HP1&DNA-high',
                                   'Immob dye')) +
    scale_fill_manual(values = c('gray', 'red', 'cyan', 'purple', 'black'), #, 'purple4', 'blue', 'tomato4'),
                      labels = c('Rest of euchromatin', 'HP1-high', 'DNA-high', 
                                 'HP1&DNA-high',
                                 'Immob dye')) +
    coord_cartesian(ylim = c(0, 0.05)) +#coord_cartesian(ylim = c(0, 0.04)) +
    labs(x = 'log(D, um2/s)', y = 'Posterior occupation', colour = 'Compartment', fill = 'Compartment') +
    theme_bw() +
    theme(aspect.ratio = 0.5)
  ggsave('SA-SPT_results/posterior_line_ctrls_only.svg', width = 7, height = 2.5)
  
  
  group_by(allData, compartment, type, bioRep, rep, diff_coef) %>%
    summarise(mean_posterior_occupation = sum(mean_posterior_occupation)) %>%
    group_by(compartment, bioRep, rep, type) %>%
    mutate(mean_posterior_occupation = mean_posterior_occupation/sum(mean_posterior_occupation)) %>%
    group_by(compartment, bioRep, type, rep, population = case_when(log10(diff_coef) < -1 ~ 'Bound', T ~ 'Unbound')) %>%
    summarise(mean_posterior_occupation = sum(mean_posterior_occupation)) %>%
    group_by(compartment, bioRep, type, population) %>%
    summarise(CI1 = quantile(mean_posterior_occupation, 0.025),
              CI2 = quantile(mean_posterior_occupation, 0.975),
              mean_posterior_occupation = mean(mean_posterior_occupation)) %>%
    mutate(compartment = factor(compartment, levels = c('Eu', 'HP1 foci', 'DNA foci',
                                                        'HP1&DNA foci',
                                                        'Immob dye')),
           type = factor(type, levels = c('real', 'pseudo')),
           population = factor(population, levels = c('Bound', 'Unbound'))) -> allData_summary
  
  allData_summary %>%
    ggplot(aes(x = interaction(bioRep, type, compartment),
               y = mean_posterior_occupation,
               fill = interaction(type, compartment),
               alpha = factor(population, levels = c('Unbound', 'Bound')),
               group = interaction(bioRep, type, compartment))) +
    geom_col(position = 'stack', colour = 'gray20', width = 0.8) +
    geom_errorbar(data = . %>% filter(population == 'Bound'),
                  aes(ymin = CI1, ymax = CI2), width = 0.3) +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c('gray', 'red', 'tomato4', 'cyan', 'blue',  'purple', 'purple4',
                                 muted('red')),
                      labels = c('Rest of euchromatin', 'HP1-high', 'HP1 pseudo',
                                 'DNA-high', 'DNA pseudo',
                                 'HP1&DNA-high', 'HP1&DNA pseudo',
                                 'Immob dye')) +
    scale_alpha_manual(values = c(0.1, 1), guide = F) +
    labs(x = 'Compartment', y = 'Proportion bound', fill = 'Compartment') +
    theme_bw() +
    theme(aspect.ratio = 0.5)
  ggsave('SA-SPT_results/occupations_bars2_2_immob.svg', width = 5.5, height = 2.5)
}

#Bound/unbound x intensity
{
  #Unbound HP1 concentration ratios
  filter(allData_summary, population == 'Bound', compartment != 'DNA foci, sub', compartment != 'Eu') %>% #View()
    left_join(data.frame(compartment = rep(c('DNA foci', 'HP1 foci', 'HP1&DNA foci'), each = 4),
                         type = rep(c('real', 'pseudo'), 6),
                         bioRep = rep(c('Rep1', 'Rep1', 'Rep2', 'Rep2'), 3),
                         HP1_dens = c(968.62, 855.72, 1840.13, 1578.81, #densities calculated in Analyse_blinking_HP1
                                      1855.45, 840.27, 2962.76, 1366.74,
                                      2011.38, 885.91, 3336.00, 1458.66))) %>%
    mutate(Unbound = 1-mean_posterior_occupation) %>% group_by(compartment, bioRep) %>%
    summarise(Unbound_frac_ratio = Unbound[2]/Unbound[1], HP1_dens_ratio = HP1_dens[2]/HP1_dens[1]) %>%
    mutate(Unbound_conc_ratio = Unbound_frac_ratio*HP1_dens_ratio) -> unbound_conc
  ggplot(unbound_conc, aes(x = compartment, y = Unbound_conc_ratio, fill = compartment, group = bioRep)) +
    geom_col(position = 'dodge', colour = 'gray20', width = 0.8) +
    geom_hline(yintercept = 1, colour = 'black') +
    scale_fill_manual(values = c('red', 'cyan', 'purple'),
                      labels = c('HP1-high', 'DNA-high', 'HP1&DNA-high'),
                      guide = 'none') +
    #scale_alpha_manual(values = c(1, 0.5)) +
    labs(x = 'Compartment', y = '[Unbound HP1]_foci/\n[Unbound HP1]_pseudo', fill = 'Compartment',
         title = 'Enrichment in unbound HP1b\nwithin different foci',
         subtitle = 'Compared to pseudo controls') +
    theme_bw() +
    theme(aspect.ratio = 0.8, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  ggsave('SA-SPT_results/HP1free_conc_ratios.svg', width = 3, height = 3)
  
  #Bound HP1 Vs DNA concentration - not used in the paper
  # #NB! Linear relationship only applicable if [HP1]>>[H3K9me3].
  # #In reality, the relationship with [H3K9me3]free is hyperbolic, and the relationship with [H3K9me3]tot is weird
  # #However, the binding curve is sub-linear -> if enrichment in bound HP1 is larger than enrichment in chromatin, defo sth else is going on
  # filter(allData_summary, population == 'Bound', compartment != 'DNA foci, sub', compartment != 'Eu') %>% #View()
  #   left_join(data.frame(compartment = rep(c('DNA foci', 'HP1 foci', 'HP1&DNA foci'), each = 4),
  #                        type = rep(c('real', 'pseudo'), 6),
  #                        bioRep = rep(c('Rep1', 'Rep1', 'Rep2', 'Rep2'), 3),
  #                        DNA_dens = c(7.34, 4.86, 6.94, 4.79, #densities calculated in Analyse_blinking_HP1
  #                                     5.16, 5.07, 4.65, 4.9,
  #                                     7.58, 5.04, 6.80, 4.66))) %>% 
  #   group_by(compartment, bioRep) %>%
  #   summarise(Bound_frac_ratio = mean_posterior_occupation[2]/mean_posterior_occupation[1],
  #             DNA_dens_ratio = DNA_dens[2]/DNA_dens[1]) %>%
  #   mutate(Bound_enrichm_over_DNA = Bound_frac_ratio/DNA_dens_ratio) -> bound_conc
  # ggplot(bound_conc, aes(x = compartment, y = Bound_enrichm_over_DNA, fill = compartment, group = bioRep)) +
  #   geom_col(position = 'dodge', colour = 'gray20', width = 0.8) +
  #   geom_hline(yintercept = 1, colour = 'black') +
  #   scale_fill_manual(values = c('red', 'cyan', 'purple'),
  #                     labels = c('HP1-high', 'DNA-high', 'HP1&DNA-high'),
  #                     guide = 'none') +
  #   #scale_alpha_manual(values = c(1, 0.5)) +
  #   labs(x = 'Compartment', y = 'Enrichment in bound HP1/\nEnrichment in chromatin', fill = 'Compartment',
  #        title = 'Enrichment in HP1b binding sites\nwithin different foci',
  #        subtitle = 'Compared to pseudo controls\nAssuming linearity, i.e. [HP1]>>[chromatin]') +
  #   theme_bw() +
  #   theme(aspect.ratio = 0.8, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  # ggsave('SA-SPT_results/BoundHP1_vs_chromatin.png', width = 4, height = 3)
}
