library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(scales)
library(ggpubr)
library(EBImage)
library(rdist)
#library(deldir)
#library(ggvoronoi)
library(tictoc)
library(dbscan)
library(png)
library(RcppRoll) #NB! roll_mean from RcppRoll is about 8-10 times faster than rollmean from zoo
library(parallel)
library(broom)

rm(list = ls())
set.seed(100500)

folder = '/home/sasha/Science/PhD/Laue_lab/Microscopy/HP1b_live_blinkingDyes/20221108_Halo-JF639b_Hoechst-SPY505'
load_data = T
partial_load = F

data_dir_SPT = 'Tracks_400nm_noGaps' #Tracked with 400 nm radius and without gaps - to analyse movement
data_dir_SMLM = 'Tracks_300nm_2gapsNoAdj' #Tracked with 300 nm radius and allowing 2 gaps, but wo radius adjustment - to reconstruct images

setwd(folder)

#### Read in data ####
#SPT
files = list.files(data_dir_SPT, pattern = 'trackPositions.csv')
#files_pfi = list.files(data_dir, pattern = 'positionsFramesIntensity.csv')

if (load_data) {
  data_list = readRDS('Data_list_SPT.RDS')
} else if (partial_load) {
  data_list = readRDS('Data_list_SPT.RDS')
  lapply(seq_along(files), function(i) {
    print(files[i])
    read.table(file.path(data_dir_SPT, files[i]), sep = ',', header = F)[1:4] -> tmp
    colnames(tmp) = c('Track', 'Frame', 'x', 'y')
    tmp$Dish = str_split(files[i], '_')[[1]][1]
    tmp$Pos = str_split(files[i], '_')[[1]][2]
    group_by(tmp, Track) %>% mutate(Track_length = n()) -> tmp
    return(tmp)
  }) %>% append(data_list, .) -> data_list
  saveRDS(data_list, 'Data_list_SPT.RDS')
} else {
  lapply(seq_along(files), function(i) {
    print(files[i])
    read.table(file.path(data_dir_SPT, files[i]), sep = ',', header = F)[1:4] -> tmp
    colnames(tmp) = c('Track', 'Frame', 'x', 'y')
    tmp$Dish = str_split(files[i], '_')[[1]][1]
    tmp$Pos = str_split(files[i], '_')[[1]][2]
    group_by(tmp, Track) %>% mutate(Track_length = n()) -> tmp
    return(tmp)
  }) -> data_list
  names(data_list) = sapply(data_list, function(x) paste0(x$Dish[[1]], '_', x$Pos[[1]]))
  saveRDS(data_list, 'Data_list_SPT.RDS')
}
bind_rows(data_list) -> data

#SMLM
files = list.files(data_dir_SMLM, pattern = 'trackPositions.csv')

if (load_data) {
  data_list_SMLM = readRDS('Data_list_SMLM.RDS')
} else if (partial_load) {
  data_list_SMLM = readRDS('Data_list_SMLM.RDS')
  lapply(seq_along(files), function(i) {
    print(files[i])
    read.table(file.path(data_dir_SMLM, files[i]), sep = ',', header = F)[1:4] -> tmp
    colnames(tmp) = c('Track', 'Frame', 'x', 'y')
    tmp$Pos = str_split(files[i], '_')[[1]][1]
    group_by(tmp, Track) %>% mutate(Track_length = n()) -> tmp
    filter(tmp, Frame > 1000) -> tmp
    return(tmp)
  }) %>% append(data_list_SMLM, .) -> data_list_SMLM
  saveRDS(data_list_SMLM, 'Data_list_SMLM.RDS')
} else {
  lapply(seq_along(files), function(i) {
    print(files[i])
    read.table(file.path(data_dir_SMLM, files[i]), sep = ',', header = F)[1:4] -> tmp
    colnames(tmp) = c('Track', 'Frame', 'x', 'y')
    tmp$Dish = str_split(files[i], '_')[[1]][1]
    tmp$Pos = str_split(files[i], '_')[[1]][2]
    group_by(tmp, Track) %>% mutate(Track_length = n()) -> tmp
    filter(tmp, Frame > 1000) -> tmp
    return(tmp)
  }) -> data_list_SMLM
  names(data_list_SMLM) = sapply(data_list_SMLM, function(x) paste0(x$Dish[[1]], '_', x$Pos[[1]]))
  saveRDS(data_list_SMLM, 'Data_list_SMLM.RDS')
}
bind_rows(data_list_SMLM) -> data_SMLM

#### Plot distributions of track lengths, jump distance histograms, etc - not used in the paper ####
# #Num/density locs
# {
#   #Distance between localisations
#   data %>% group_by(Dish, Pos, Frame) %>%
#     summarise(pdist_mean = mean(c(pdist(matrix(c(x, y), nrow=n())))[c(pdist(matrix(c(x, y), nrow=n())))!=0]),
#               pdist_min = min(c(pdist(matrix(c(x, y), nrow=n())))[c(pdist(matrix(c(x, y), nrow=n())))!=0]),
#               pdist_max = max(c(pdist(matrix(c(x, y), nrow=n())))[c(pdist(matrix(c(x, y), nrow=n())))!=0])) %>%
#     mutate(pdist_mean = replace(pdist_mean, is.na(pdist_mean), 15000),
#            pdist_min = replace(pdist_min, pdist_min == Inf, 15000),
#            pdist_max = replace(pdist_max, pdist_max == Inf, 15000)) -> data_pdist
#   data_pdist %>% mutate(FrameRound = ceiling(Frame/100)*100) %>% group_by(Dish, Pos, FrameRound) %>%
#     summarise(Mean_pdist_mean = mean(pdist_mean),
#               Min_pdist_mean = min(pdist_mean),
#               Mean_pdist_min = mean(pdist_min)) %>%
#     mutate(tag = case_when((Pos == 'Pos0' | Pos == 'Pos1' | Pos == 'Pos2' | Pos == 'Pos3' |
#                               Pos == 'Pos4' | Pos == 'Pos5') ~ 'Pos0-5',
#                            T ~ 'Pos6-11')) -> data_pdist_avg
#   ggplot(data_pdist_avg, aes(x = FrameRound, y = Mean_pdist_mean/1000, colour = Pos)) +
#     geom_line() +
#     facet_grid(factor(interaction(Dish, tag), levels = c('Dish1.Pos0-5', 'Dish1.Pos6-11',
#                                                          'Dish2.Pos0-5', 'Dish2.Pos6-11'))~.) +
#     scale_y_continuous(breaks = seq(4, 14, by = 2), limits = c(3, 15)) +
#     labs(x = 'Frame', y = 'Mean distance, um', title = 'Mean pairwise distance between locs in a frame',
#          subtitle = 'Smoothed (average per 100 frames)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.15, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('meanLocDist.png', width = 10, height = 6)
#   ggplot(data_pdist_avg, aes(x = FrameRound, y = Mean_pdist_min/1000, colour = Pos)) +
#     geom_line() +
#     facet_grid(factor(interaction(Dish, tag), levels = c('Dish1.Pos0-5', 'Dish1.Pos6-11',
#                                                          'Dish2.Pos0-5', 'Dish2.Pos6-11'))~.) +
#     labs(x = 'Frame', y = 'Min distance, um', title = 'Min pairwise distance between locs in a frame',
#          subtitle = 'Smoothed (average per 100 frames)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.15, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('minLocDist.png', width = 10, height = 6)
#   ggplot(data_pdist_avg, aes(x = FrameRound, y = Mean_pdist_min/1000, colour = Pos)) +
#     geom_line() +
#     facet_grid(factor(interaction(Dish, tag), levels = c('Dish1.Pos0-5', 'Dish1.Pos6-11',
#                                                          'Dish2.Pos0-5', 'Dish2.Pos6-11'))~.) +
#     coord_cartesian(ylim = c(0, 3)) +
#     labs(x = 'Frame', y = 'Min distance, um', title = 'Min pairwise distance between locs in a frame',
#          subtitle = 'Smoothed (average per 100 frames)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.15, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('minLocDist_zoom.png', width = 10, height = 6)
#   data_pdist %>% 
#     mutate(tag = case_when((Pos == 'Pos0' | Pos == 'Pos1' | Pos == 'Pos2' | Pos == 'Pos3' |
#                               Pos == 'Pos4' | Pos == 'Pos5') ~ 'Pos0-5',
#                            T ~ 'Pos6-11')) -> data_pdist
#   ggplot(data_pdist, aes(x = pdist_min/1000, colour = Pos)) +
#     geom_density(stat = 'ecdf') +
#     facet_grid(factor(interaction(Dish, tag), levels = c('Dish1.Pos0-5', 'Dish1.Pos6-11',
#                                                          'Dish2.Pos0-5', 'Dish2.Pos6-11'))~.) +
#     scale_x_continuous(minor_breaks = seq(0, 15, by = 1), limits = c(0, 15.5)) +
#     labs(x = 'Min distance, um', y = 'Cumulative density', title = 'Min pairwise distance between locs in a frame',
#          subtitle = '15,000 = single loc') +
#     theme_bw() +
#     theme(aspect.ratio = 0.15, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('minLocDist_cumul.png', width = 10, height = 6)
#   data %>% group_by(Dish, Pos, Frame) %>%
#     summarise(pdist = c(pdist(matrix(c(x, y), nrow=n())))[c(pdist(matrix(c(x, y), nrow=n())))!=0]) %>%
#     mutate(tag = case_when((Pos == 'Pos0' | Pos == 'Pos1' | Pos == 'Pos2' | Pos == 'Pos3' |
#                               Pos == 'Pos4' | Pos == 'Pos5') ~ 'Pos0-5',
#                            T ~ 'Pos6-11')) -> data_pdist_all
#   filter(data_pdist_all, pdist < 1000) -> tmp
#   mutate(tmp, cycle = paste0('cycle', ceiling(Frame/10000)-1),
#          Frame = Frame - (ceiling(Frame/10000)-1)*10000) %>%
#     ggplot(aes(x = pdist/1000, y = ave(..count.., group, FUN = cumsum),
#                colour = Pos, group = interaction(Dish, Pos, cycle))) +
#     geom_freqpoly(binwidth = 0.02) +
#     facet_grid(cycle~factor(interaction(Dish, tag), levels = c('Dish1.Pos0-5', 'Dish1.Pos6-11',
#                                                                'Dish2.Pos0-5', 'Dish2.Pos6-11')),
#                scales = 'free_y') +
#     scale_x_continuous(breaks = seq(0, 1, by = 0.2), minor_breaks = seq(0, 1, by = 0.1)) +
#     scale_y_continuous(sec.axis = sec_axis(trans=~./10000, name = 'Cumul count per frame')) +
#     labs(x = 'Distance, um', y = 'Cumulative count', title = 'Pairwise distances <1um between locs in a frame',
#          subtitle = '"Cycles" are artificial (for scaling)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.4, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('locDist_cumulHist_1um.png', width = 12, height = 5.5)
#   filter(tmp, pdist < 500) %>%
#     mutate(cycle = paste0('cycle', ceiling(Frame/10000)-1),
#            Frame = Frame - (ceiling(Frame/10000)-1)*10000) %>%
#     ggplot(aes(x = pdist/1000, y = ave(..count.., group, FUN = cumsum),
#                colour = Pos, group = interaction(Dish, Pos, cycle))) +
#     geom_freqpoly(binwidth = 0.01) +
#     facet_grid(cycle~factor(interaction(Dish, tag), levels = c('Dish1.Pos0-5', 'Dish1.Pos6-11',
#                                                                'Dish2.Pos0-5', 'Dish2.Pos6-11')),
#                scales = 'free_y') +
#     scale_x_continuous(breaks = seq(0, 1, by = 0.1), minor_breaks = seq(0, 1, by = 0.05)) +
#     scale_y_continuous(sec.axis = sec_axis(trans=~./10000, name = 'Cumul count per frame')) +
#     labs(x = 'Distance, um', y = 'Cumulative count', title = 'Pairwise distances <500nm between locs in a frame',
#          subtitle = '"Cycles" are artificial (for scaling)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.4, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('locDist_cumulHist_500nm.png', width = 12, height = 5.5)
#   
#   filter(tmp, pdist < 500) %>%
#     mutate(cycle = paste0('cycle', ceiling(Frame/5000)-1),
#            Frame = Frame - (ceiling(Frame/5000)-1)*5000) %>%
#     ggplot(aes(x = pdist/1000, y = ave(..count.., group, FUN = cumsum),
#                colour = Pos, group = interaction(Dish, Pos, cycle))) +
#     geom_freqpoly(binwidth = 0.01) +
#     facet_grid(cycle~factor(interaction(Dish, tag), levels = c('Dish1.Pos0-5', 'Dish1.Pos6-11',
#                                                                'Dish2.Pos0-5', 'Dish2.Pos6-11')),
#                scales = 'free_y') +
#     scale_x_continuous(breaks = seq(0, 1, by = 0.1), minor_breaks = seq(0, 1, by = 0.05)) +
#     scale_y_continuous(sec.axis = sec_axis(trans=~./5000, name = 'Cumul count per frame')) +
#     labs(x = 'Distance, um', y = 'Cumulative count', title = 'Pairwise distances <500nm between locs in a frame',
#          subtitle = '"Cycles" are artificial (for scaling)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.4, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('locDist_cumulHist_500nm_5000frCycles.png', width = 12, height = 8)
#   
#   #Number of localisations
#   data %>% group_by(Dish, Pos, Track) %>% summarise(Track_length = n()) %>% group_by(Dish, Pos) %>%
#     summarise(Num_points = sum(Track_length), Num_tracks = n(),
#               Num_realTracks = sum(Track_length > 1)) -> data_numPoints
#   data_SMLM %>% group_by(Dish, Pos, Track) %>% summarise(Track_length = n()) %>% group_by(Dish, Pos) %>%
#     summarise(Num_points = sum(Track_length), Num_tracks = n()) -> data_numPoints_SMLM
#   ggplot(data_numPoints, aes(x = Pos, y = Num_points, fill = Pos)) +
#     geom_bar(stat = 'identity', colour = 'gray20',
#              position = position_dodge2(width = 0.9, preserve = "single")) +
#     facet_grid(~Dish, scales = 'free_x') +
#     scale_y_continuous(sec.axis = sec_axis(trans=~./10000, name = 'Num loc per frame')) +
#     labs(x = 'Pos', y = 'Total num loc', title = 'Number of localisations') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5),
#           legend.position = 'none', axis.text.x = element_text(angle = 30, vjust = 0.8))
#   ggsave('numLoc.png', width = 6, height = 2.5)
#   ggplot(data_numPoints, aes(x = Pos, y = Num_tracks, fill = Pos)) +
#     geom_bar(stat = 'identity', colour = 'gray20',
#              position = position_dodge2(width = 0.9, preserve = "single")) +
#     facet_grid(~Dish, scales = 'free_x') +
#     labs(x = 'Pos', y = 'Number of tracks', subtitle = 'Tracking r=400 nm, no gaps (for diff)',
#          title = 'Number of tracks') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5),
#           legend.position = 'none', axis.text.x = element_text(angle = 30, vjust = 0.8)) -> numTrack1
#   ggplot(data_numPoints, aes(x = Pos, y = Num_realTracks, fill = Pos)) +
#     geom_bar(stat = 'identity', colour = 'gray20',
#              position = position_dodge2(width = 0.9, preserve = "single")) +
#     facet_grid(~Dish, scales = 'free_x') +
#     labs(x = 'Pos', y = 'Number of tracks', subtitle = 'Tracking r=400 nm, no gaps (for diff)',
#          title = 'Number of tracks>1 pos') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5),
#           legend.position = 'none', axis.text.x = element_text(angle = 30, vjust = 0.8)) -> numTrack11
#   ggplot(data_numPoints_SMLM, aes(x = Pos, y = Num_tracks, fill = Pos)) +
#     geom_bar(stat = 'identity', colour = 'gray20',
#              position = position_dodge2(width = 0.9, preserve = "single")) +
#     facet_grid(~Dish, scales = 'free_x') +
#     labs(x = 'Pos', y = 'Number of tracks', subtitle = 'Tracking r=300 nm, <3 gaps (for image)',
#          title = 'Number of tracks') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5),
#           legend.position = 'none', axis.text.x = element_text(angle = 30, vjust = 0.8)) -> numTrack2
#   ggarrange(numTrack1, numTrack11, numTrack2, nrow = 3)
#   ggsave('numTrack.png', width = 6, height = 7)
#   data %>% group_by(Dish, Pos, Track) %>% summarise(Track_length = n()) %>%
#     filter(Track_length > 1) %>% group_by(Dish, Pos, Track_length) %>%
#     summarise(Num_tracks = n()) -> data_numLenTr
#   data_numLenTr %>% mutate(Track_length = case_when(Track_length > 5 ~ '>5',
#                                                     T ~ as.character(Track_length))) %>%
#     group_by(Dish, Pos, Track_length) %>% summarise(Num_tracks = sum(Num_tracks)) -> data_numLenTr2
#   ggplot(data_numLenTr2, aes(x = Pos, y = Num_tracks, fill = Track_length)) +
#     geom_bar(stat = 'identity', colour = 'gray20') +
#     facet_grid(~Dish, scales = 'free_x') +
#     labs(x = 'Pos', y = 'Number of tracks', title = 'Number of tracks > 1 pos',
#          subtitle = 'Tracking r=400 nm, no gaps (for diff)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5),
#           axis.text.x = element_text(angle = 30, vjust = 0.8))
#   ggsave('numLenTrack.png', width = 7, height = 3)
#   data_SMLM %>% group_by(Dish, Pos, Track) %>% summarise(Track_length = n()) %>%
#     group_by(Dish, Pos, Track_length) %>% summarise(Num_tracks = n()) -> data_numLenTr_SMLM
#   data_numLenTr_SMLM %>% mutate(Track_length = case_when(Track_length > 5 ~ '>5',
#                                                          T ~ as.character(Track_length))) %>%
#     group_by(Dish, Pos, Track_length) %>% summarise(Num_tracks = sum(Num_tracks)) -> data_numLenTr_SMLM2
#   ggplot(data_numLenTr_SMLM2, aes(x = Pos, y = Num_tracks, fill = Track_length)) +
#     geom_bar(stat = 'identity', colour = 'gray20') +
#     facet_grid(~Dish, scales = 'free_x') +
#     labs(x = 'Pos', y = 'Number of tracks', title = 'Number of tracks',
#          subtitle = 'Tracking r=300 nm, <3 gaps (for image)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5),
#           axis.text.x = element_text(angle = 30, vjust = 0.8))
#   ggsave('numLenTrack_SMLM.png', width = 7, height = 3)
#   data %>% mutate(FrameRound = ceiling(Frame/100)*100) %>% group_by(Dish, Pos, FrameRound) %>%
#     summarise(Num_loc = n()/100) %>%
#     mutate(tag = case_when((Pos == 'Pos0' | Pos == 'Pos1' | Pos == 'Pos2' | Pos == 'Pos3' |
#                               Pos == 'Pos4' | Pos == 'Pos5') ~ 'Pos0-5',
#                            T ~ 'Pos6-11')) -> data_locPerFr
#   ggplot(data_locPerFr, aes(x = FrameRound, y = Num_loc, colour = Pos)) +
#     geom_line() +
#     facet_grid(factor(interaction(Dish, tag), levels = c('Dish1.Pos0-5', 'Dish1.Pos6-11',
#                                                          'Dish2.Pos0-5', 'Dish2.Pos6-11'))~.) +
#     labs(x = 'Frame', y = 'Number of locs', title = 'Number of localisations per frame',
#          subtitle = 'Average per 100 frames') +
#     theme_bw() +
#     theme(aspect.ratio = 0.15, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('numLocFr.png', width = 10, height = 6)
#   data %>% mutate(FrameRound = ceiling(Frame/100)*100) %>% group_by(Dish, Pos, FrameRound) %>%
#     summarise(Num_loc = n()/100) %>%
#     mutate(cycle = paste0('cycle', ceiling(FrameRound/10000)-1),
#            FrameRound = FrameRound - (ceiling(FrameRound/10000)-1)*10000) %>%
#     mutate(tag = case_when((Pos == 'Pos0' | Pos == 'Pos1' | Pos == 'Pos2' | Pos == 'Pos3' |
#                               Pos == 'Pos4' | Pos == 'Pos5') ~ 'Pos0-5',
#                            T ~ 'Pos6-11')) -> data_locPerFr
#   ggplot(data_locPerFr, aes(x = FrameRound, y = Num_loc, colour = Pos,
#                             group = interaction(Dish, Pos, cycle))) +
#     geom_line() +
#     facet_grid(cycle~factor(interaction(Dish, tag), levels = c('Dish1.Pos0-5', 'Dish1.Pos6-11',
#                                                                'Dish2.Pos0-5', 'Dish2.Pos6-11')),
#                scale = 'free_y') +
#     labs(x = 'Frame', y = 'Number of locs', title = 'Number of localisations per frame',
#          subtitle = 'Average per 100 frames; "cycles" are artificial (for scaling)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.3, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('numLocFr2.png', width = 12, height = 5.5)
# }
# 
# #Track lengths
# {
#   data %>% group_by(Dish, Pos, Track) %>% summarise(Track_length = n()) -> data_trLen
#   data_SMLM %>% group_by(Dish, Pos, Track) %>% summarise(Track_length = n()) -> data_trLen_SMLM
#   ggplot(data_trLen, aes(x = Track_length)) +
#     geom_histogram(binwidth = 1, aes(y = stat(width*density)), colour = 'gray20', fill = 'gray80') +
#     facet_grid(~Dish) +
#     labs(title = 'Distribution of track lengths',
#          subtitle = 'All Pos\nTracking r=400 nm, no gaps (for diff)',
#          x = 'Track length') +
#     coord_cartesian(xlim = c(0.8, 10)) +
#     scale_x_continuous(breaks = 1:10) +
#     scale_y_continuous(labels = percent_format()) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank()) -> trackLen1
#   ggplot(data_trLen_SMLM, aes(x = Track_length)) +
#     geom_histogram(binwidth = 1, aes(y = stat(width*density)), colour = 'gray20', fill = 'gray80') +
#     facet_grid(~Dish) +
#     labs(title = 'Distribution of track lengths',
#          subtitle = 'All Pos\nTracking r=300 nm, <3 gaps (for image)',
#          x = 'Track length') +
#     coord_cartesian(xlim = c(0.8, 10)) +
#     scale_x_continuous(breaks = 1:10) +
#     scale_y_continuous(labels = percent_format()) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank()) -> trackLen2
#   ggarrange(trackLen1, trackLen2, nrow = 2, common.legend = T)
#   ggsave('TrackLen.png', width = 6, height = 5)
#   ggplot(data_trLen, aes(x = Track_length)) +
#     geom_histogram(binwidth = 1, aes(y = stat(width*density)), colour = 'gray20', fill = 'gray80') +
#     facet_grid(~Dish) +
#     labs(title = 'Distribution of track lengths',
#          subtitle = 'All Pos\nTracking r=400 nm, no gaps (for diff)',
#          x = 'Track length') +
#     coord_cartesian(xlim = c(5.95, 16), ylim = c(0, 0.015)) +
#     scale_x_continuous(breaks = 6:16) +
#     scale_y_continuous(breaks = seq(0, 0.015, by = 0.0025), labels = percent_format()) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank()) -> trackLenLong1
#   ggplot(data_trLen_SMLM, aes(x = Track_length)) +
#     geom_histogram(binwidth = 1, aes(y = stat(width*density)), colour = 'gray20', fill = 'gray80') +
#     facet_grid(~Dish) +
#     labs(title = 'Distribution of track lengths',
#          subtitle = 'All Pos\nTracking r=300 nm, <3 gaps (for image)',
#          x = 'Track length') +
#     coord_cartesian(xlim = c(5.95, 16), ylim = c(0, 0.032)) +
#     scale_x_continuous(breaks = 6:16) +
#     scale_y_continuous(breaks = seq(0, 0.04, by = 0.005), labels = percent_format()) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank()) -> trackLenLong2
#   ggarrange(trackLenLong1, trackLenLong2, nrow = 2, common.legend = T)
#   ggsave('TrackLen_longZoom.png', width = 6, height = 5)
#   
#   ggplot(data_trLen, aes(x = Track_length, colour = Pos)) +
#     geom_freqpoly(binwidth = 1, aes(y = stat(width*density)), size = 0.2) +
#     facet_grid(~Dish) +
#     labs(title = 'Distribution of track lengths',
#          subtitle = 'Tracking r=400 nm, no gaps (for diff)', x = 'Track length') +
#     coord_cartesian(xlim = c(0, 10)) +
#     scale_x_continuous(breaks = 0:10) +
#     scale_y_continuous(labels = percent_format()) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank()) -> trackLen_perPos1
#   ggplot(data_trLen_SMLM, aes(x = Track_length, colour = Pos)) +
#     geom_freqpoly(binwidth = 1, aes(y = stat(width*density)), size = 0.2) +
#     facet_grid(~Dish) +
#     #facet_grid(Dataset~.) +
#     labs(title = 'Distribution of track lengths',
#          subtitle = 'Tracking r=300 nm, <3 gaps (for image)', x = 'Track length') +
#     coord_cartesian(xlim = c(0, 10)) +
#     scale_x_continuous(breaks = 0:10) +
#     scale_y_continuous(labels = percent_format()) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank()) -> trackLen_perPos2
#   ggarrange(trackLen_perPos1, trackLen_perPos2, nrow = 2, common.legend = T)
#   ggsave('TrackLen_perPos.png', width = 6, height = 6)
# }
# 
# #Jump distances
# {
#   ##QC check
#   data %>% filter(Dish == 'Dish1' & Pos == 'Pos0' & Frame > 19500 |
#                     Dish == 'Dish1' & Pos == 'Pos1' & Frame > 20000 |
#                     Dish == 'Dish1' & Pos == 'Pos2' & Frame > 20000 |
#                     Dish == 'Dish1' & Pos == 'Pos3' & Frame > 10000 |
#                     Dish == 'Dish1' & Pos == 'Pos4' & Frame > 10000 |
#                     Dish == 'Dish1' & Pos == 'Pos5' & Frame > 21000 |
#                     Dish == 'Dish1' & Pos == 'Pos6' & Frame > 12000 |
#                     Dish == 'Dish1' & Pos == 'Pos7' & Frame > 20000 |
#                     Dish == 'Dish1' & Pos == 'Pos8' & Frame > 11000 |
#                     Dish == 'Dish1' & Pos == 'Pos9' & Frame > 21000 |
#                     Dish == 'Dish1' & Pos == 'Pos10' & Frame > 11000 |
#                     Dish == 'Dish1' & Pos == 'Pos11' & Frame > 19000 |
#                     Dish == 'Dish2' & Pos == 'Pos0' & Frame > 25000 |
#                     Dish == 'Dish2' & Pos == 'Pos1' & Frame > 24000 |
#                     Dish == 'Dish2' & Pos == 'Pos2' & Frame > 24000 |
#                     Dish == 'Dish2' & Pos == 'Pos3' & Frame > 20000 |
#                     Dish == 'Dish2' & Pos == 'Pos4' & Frame > 23000 |
#                     Dish == 'Dish2' & Pos == 'Pos5' & Frame > 10000 |
#                     Dish == 'Dish2' & Pos == 'Pos6' & Frame > 33000 |
#                     Dish == 'Dish2' & Pos == 'Pos7' & Frame > 17000 |
#                     Dish == 'Dish2' & Pos == 'Pos8' & Frame > 14000) -> data_filtered #also edited based on reproducibility in "cycles"
#   data_filtered %>% group_by(Dish, Pos, Frame) %>% summarise(N = n()) %>%
#     mutate(tag = case_when((Pos == 'Pos0' | Pos == 'Pos1' | Pos == 'Pos2' | Pos == 'Pos3' |
#                               Pos == 'Pos4' | Pos == 'Pos5') ~ 'Pos0-5',
#                            T ~ 'Pos6-11')) %>%
#     ggplot(aes(x = N, colour = Pos, y = ..density..)) +
#     geom_freqpoly(binwidth = 1) +
#     facet_grid(factor(interaction(Dish, tag), levels = c('Dish1.Pos0-5', 'Dish1.Pos6-11',
#                                                          'Dish2.Pos0-5', 'Dish2.Pos6-11'))~.) +
#     labs(title = 'Num loc/frame in data used for JDA', x = 'Num loc in frame', y = 'Count',
#          subtitle = 'All Pos, dens < 10 loc/fr') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5)
#   ggsave('numLocFr_filtered.png', width = 5, height = 6)
#   data_filtered %>% group_by(Dish, Pos, Frame) %>%
#     summarise(pdist = c(pdist(matrix(c(x, y), nrow=n())))[c(pdist(matrix(c(x, y), nrow=n())))!=0]) -> data_pdist_all_filt
#   data_pdist_all_filt %>% filter(pdist < 500) %>%
#     ggplot(aes(x = pdist/1000, y = ave(..count.., group, FUN = cumsum),
#                colour = Pos, group = interaction(Dish, Pos))) +
#     geom_freqpoly(binwidth = 0.01) +
#     facet_grid(Dish~.) +
#     #facet_grid(tag~., labeller = labeller(tag = c('TRUE' = 'Pos0-7', 'FALSE' = 'Pos8-14'))) +
#     scale_x_continuous(breaks = seq(0, 1, by = 0.1), minor_breaks = seq(0, 1, by = 0.05)) +
#     scale_y_continuous(limits = c(0, 1000), breaks = seq(0, 1000, by = 200),
#                        sec.axis = sec_axis(trans=~./10000, breaks = seq(0, 0.1, by = 0.02),
#                                            name = 'Cumul count per frame')) +
#     labs(x = 'Distance, um', y = 'Cumulative count', title = 'Pairwise distances <500nm between locs in a frame',
#          subtitle = 'In data used for JDA\nAll Pos, dens < 10 loc/fr') +
#     theme_bw() +
#     theme(aspect.ratio = 0.4, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('locDist_cumulHist_500nm_filtered.png', width = 5.5, height = 6)
#   data_pdist_all_filt %>% 
#     mutate(tag = case_when((Pos == 'Pos0' | Pos == 'Pos1' | Pos == 'Pos2' | Pos == 'Pos3' |
#                               Pos == 'Pos4' | Pos == 'Pos5') ~ 'Pos0-5',
#                            T ~ 'Pos6-11')) %>%
#     ggplot(aes(x = pdist/1000, colour = Pos, group = interaction(Dish, Pos))) +
#     geom_density(stat = 'ecdf') +
#     facet_grid(factor(interaction(Dish, tag), levels = c('Dish1.Pos0-5', 'Dish1.Pos6-11',
#                                                          'Dish2.Pos0-5', 'Dish2.Pos6-11'))~.) +
#     labs(x = 'Distance, um', y = 'Cumulative density', title = 'Pairwise distances between locs in a frame',
#          subtitle = 'In data used for JDA\nAll Pos, dens < 10 loc/fr') +
#     theme_bw() +
#     theme(aspect.ratio = 0.4, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('locDist_cumul_filtered.png', width = 5.5, height = 6)
#   
#   #Overall JD
#   data_filtered %>%
#     filter(Track_length == 2) %>% group_by(Dish, Pos, Track) %>%
#     summarise(disp = sqrt(diff(x)^2 + diff(y)^2), x = mean(x), y = mean(y), Frame1 = Frame[1]) -> jumpDist
#   ggplot(jumpDist, aes(x = disp)) +
#     geom_histogram(aes(y = ..density..), binwidth = 10, fill = NA, colour = 'gray20') +
#     geom_density(colour = 'red') +
#     facet_grid(Dish~.) +
#     labs(title = 'Jump distance distribution', x = 'Displacement, nm', subtitle = 'All Pos, dens < 10 loc/fr, l = 2') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   ggsave('JumpDist.png', width = 3, height = 4)
#   ggplot(jumpDist, aes(x = disp, colour = Pos, linetype = Dish)) +
#     geom_density() +
#     labs(title = 'Jump distance distribution', x = 'Displacement, nm', subtitle = 'Dens < 10 loc/fr, l = 2') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   ggsave('JumpDist_perPos1.png', width = 5, height = 5)
#   ggplot(jumpDist, aes(x = disp, colour = Pos)) +
#     geom_density() +
#     facet_grid(Dish~.) +
#     labs(title = 'Jump distance distribution', x = 'Displacement, nm', subtitle = 'Dens < 10 loc/fr, l = 2') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   ggsave('JumpDist_perPos2.png', width = 4, height = 4)
#   mutate(jumpDist, cycle = paste0('cycle', ceiling(Frame1/10000)-1)) %>%
#     ggplot(aes(x = disp, colour = Pos, linetype = cycle)) +
#     geom_density() +
#     facet_grid(cycle~Dish) +
#     labs(title = 'Jump distance distribution', x = 'Displacement, nm',
#          subtitle = 'Dens < 10 loc/fr, l = 2\nArtificial "cycles" to estimate reproducibility') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   ggsave('JumpDist_perPosCycle1.png', width = 6, height = 4.5)
#   mutate(jumpDist, cycle = paste0('cycle', ceiling(Frame1/10000)-1)) %>%
#     ggplot(aes(x = disp, colour = Pos, linetype = cycle)) +
#     geom_density() +
#     facet_grid(Pos~Dish) +
#     labs(title = 'Jump distance distribution', x = 'Displacement, nm',
#          subtitle = 'Dens < 10 loc/fr, l = 2\nArtificial "cycles" to estimate reproducibility') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank()) #edited thresholds based on reproducibility in "cycles"
#   ggsave('JumpDist_perPosCycle2.png', width = 6, height = 12)
#   
#   #JD with increasing interval
#   for (i in 1:8) {
#     data_filtered %>% filter(Track_length > i) %>% group_by(Dish, Pos, Track) %>%
#       filter(n() > i) %>% #because initial filtering by frame cuts some trajectories up
#       mutate(diff.x = c(rep(NA, i), diff(x, i)), diff.y = c(rep(NA, i), diff(y, i))) %>%
#       mutate(disp = sqrt(diff.x^2 + diff.y^2)) %>%
#       filter(!is.na(disp)) %>% mutate(interval = i) -> tmp
#     if (i == 1) {
#       jumpDist_1234 = tmp
#     } else {
#       jumpDist_1234 = rbind(jumpDist_1234, tmp)
#     }
#   }
#   
#   ggplot(jumpDist_1234, aes(x = disp, colour = factor(interval*3, levels = c('3', '6', '9', '12',
#                                                                              '15', '18', '21', '24')))) +
#     geom_density() +
#     facet_grid(Dish~.) +
#     scale_colour_manual(values = gray.colors(8, 0.2, 0.75)) +
#     labs(title = 'Jump distance distribution', subtitle = 'All Pos, dens < 10 loc/fr',
#          x = 'Displacement, nm', colour = 'Interval, ms') +
#     lims(x = c(0, 1000), y = c(0, 0.012)) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   ggsave('JumpDist_intervals.png', width = 4.5, height = 4)
#   
#   #Jump dist as function of traj length
#   filter(jumpDist_1234, interval == 1) %>%
#     mutate(Track_length = case_when(Track_length > 5 ~ '>5', T ~ as.character(Track_length))) %>%
#     ggplot(aes(x = disp, colour = Track_length)) +
#     geom_density() +
#     facet_grid(Dish~.) +
#     #scale_colour_manual(values = gray.colors(8, 0.2, 0.75)) +
#     labs(title = 'Jump distance distribution', subtitle = 'All Pos, dens < 10 loc/fr',
#          x = 'Displacement, nm', colour = 'Track length') +
#     #lims(x = c(0, 1000)) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   ggsave('JumpDist_trLen.png', width = 4.5, height = 4)
#   
#   filter(jumpDist_1234, Track_length > 5) %>%
#     ggplot(aes(x = disp, colour = factor(interval*3, levels = c('3', '6', '9', '12',
#                                                                 '15', '18', '21', '24')))) +
#     geom_density() +
#     facet_grid(Dish~.) +
#     scale_colour_manual(values = gray.colors(8, 0.2, 0.75)) +
#     labs(title = 'Jump distance distribution', subtitle = 'Long traj (>5)',
#          x = 'Displacement, nm', colour = 'Interval, ms') +
#     lims(x = c(0, 1000), y = c(0, 0.012)) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   ggsave('JumpDist_intervals_longTraj.png', width = 5, height = 4)
#   
#   filter(jumpDist_1234, Track_length < 15) %>% #not as different from all data as before, but still a bit more spread -> use this for MSD (also for consistency)
#     ggplot(aes(x = disp, colour = factor(interval*3, levels = c('3', '6', '9', '12',
#                                                                 '15', '18', '21', '24')))) +
#     geom_density() +
#     facet_grid(Dish~.) +
#     scale_colour_manual(values = gray.colors(8, 0.2, 0.75)) +
#     labs(title = 'Jump distance distribution', subtitle = 'Short traj (<15)',
#          x = 'Displacement, nm', colour = 'Interval, ms') +
#     lims(x = c(0, 1000), y = c(0, 0.012)) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   ggsave('JumpDist_intervals_shortTraj.png', width = 4.5, height = 4)
#   
#   #MSD
#   filter(jumpDist_1234, Track_length < 15) %>% group_by(Dish, interval) %>% summarise(MSD = mean(disp^2)/1000000) -> MSD
#   filter(jumpDist_1234, Track_length < 15) %>% group_by(Dish, Pos, interval) %>% summarise(MSD = mean(disp^2)/1000000) -> MSD_perPos
#   ggplot(MSD, aes(x = interval*3, y = MSD, colour = Dish)) +
#     geom_point() +
#     scale_colour_manual(values = c('blue', 'orange')) +
#     labs(title = 'MSD curve', subtitle = 'All Pos, dens < 10 loc/fr, traj <15',
#          x = 'Interval, ms', y = 'MSD, um2') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1)
#   ggsave('MSD.png', width = 3.7, height = 3)
#   ggplot(MSD_perPos, aes(x = interval*3, y = MSD, colour = Pos, shape = Dish, linetype = Dish)) +
#     geom_point() +
#     geom_line(size = 0.2) +
#     labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15',
#          x = 'Interval, ms', y = 'MSD, um2') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1)
#   ggsave('MSD_perPos1.png', width = 3.7, height = 5)
#   ggplot(MSD_perPos, aes(x = interval*3, y = MSD, colour = Pos)) +
#     geom_point() +
#     facet_grid(Dish~.) +
#     geom_line(size = 0.2) +
#     labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15',
#          x = 'Interval, ms', y = 'MSD, um2') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1)
#   ggsave('MSD_perPos2.png', width = 3.7, height = 5)
#   ggplot(MSD_perPos, aes(x = interval*3, y = MSD, colour = Dish, group = interaction(Dish, Pos))) +
#     geom_point() +
#     geom_line(size = 0.2) +
#     scale_colour_manual(values = c('blue', 'orange')) +
#     labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15',
#          x = 'Interval, ms', y = 'MSD, um2') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1)
#   ggsave('MSD_perPos3.png', width = 3.7, height = 3)
#   split(MSD, MSD$Dish) %>% lapply(function(x) {
#     lm(log10(MSD)~log10(interval*3/1000), data = x)
#   }) -> MSDfit.list
#   MSDfit.df = data.frame(Dish = names(MSDfit.list),
#                          Intercept = sapply(MSDfit.list, function(x) {x$coefficients[[1]]}),
#                          Slope = sapply(MSDfit.list, function(x) {x$coefficients[[2]]}))
#   mutate(MSDfit.df, label = paste('D =', round(10^Intercept, 2), 'um2/s, alpha =', round(Slope, 2))) -> MSDfit.df
#   ggplot(MSD, aes(x = log10(interval*3/1000), y = log10(MSD), colour = Dish)) +
#     geom_point() +
#     #geom_line(size = 0.2) +
#     geom_abline(data = MSDfit.df, aes(intercept = Intercept, slope = Slope, colour = Dish)) +
#     # geom_text(data = MSDfit_perPos.df,
#     #           aes(x = rep(-2.4, 6), y = seq(-1.3, -1.8, length.out = 6),
#     #               label = label, colour = Pos)) +
#     # annotate('text', x = -1.9, y = -1.8,
#     #          label = paste('D =', round(10^MSDfit.df$Intercept[[1]], 2),
#     #                        'um2/s,\nalpha =', round(MSDfit.df$Slope[[1]], 2))) +
#     scale_colour_manual(values = c('blue', 'orange'),
#                         labels = paste0(MSDfit.df$Dish, '\n', MSDfit.df$label)) +
#     labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15',
#          x = 'Log(Interval, s)', y = 'Log(MSD, um2)') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1)
#   ggsave('MSD_loglog.png', width = 5, height = 3)
#   split(MSD_perPos, paste0(MSD_perPos$Dish, '_', MSD_perPos$Pos)) %>% lapply(function(x) {
#     lm(log10(MSD)~log10(interval*3/1000), data = x)
#   }) -> MSDfit_perPos.list
#   MSDfit_perPos.df = data.frame(Dish = sapply(str_split(names(MSDfit_perPos.list), '_'), function(x) x[[1]]),
#                                 Pos = sapply(str_split(names(MSDfit_perPos.list), '_'), function(x) x[[2]]),
#                                 Intercept = sapply(MSDfit_perPos.list, function(x) {x$coefficients[[1]]}),
#                                 Slope = sapply(MSDfit_perPos.list, function(x) {x$coefficients[[2]]}))
#   mutate(MSDfit_perPos.df, label = paste('D =', round(10^Intercept, 2), 'um2/s, alpha =', round(Slope, 2))) -> MSDfit_perPos.df
#   MSDfit_perPos.df %>% group_by(Dish) %>%
#     summarise(label = paste0('D = ', round(10^min(Intercept), 2), '-', round(10^max(Intercept), 2),
#                              ' um2/s, alpha = ', round(min(Slope), 2), '-', round(max(Slope), 2))) -> MSDfit_perPos.summary
#   filter(MSD_perPos, Dish == 'Dish1') %>%
#     ggplot(aes(x = log10(interval*3/1000), y = log10(MSD), colour = Pos)) +
#     geom_point() +
#     #geom_line(size = 0.2) +
#     geom_abline(data = filter(MSDfit_perPos.df, Dish == 'Dish1'),
#                 aes(intercept = Intercept, slope = Slope, colour = Pos)) +
#     # geom_text(data = MSDfit_perPos.df,
#     #           aes(x = rep(-2.4, 6), y = seq(-1.3, -1.8, length.out = 6),
#     #               label = label, colour = Pos)) +
#     # annotate('text', x = -1.9, y = -1.8,
#     #          label = paste('D =', round(10^MSDfit.df$Intercept[[1]], 2),
#     #                        'um2/s,\nalpha =', round(MSDfit.df$Slope[[1]], 2))) +
#     scale_colour_discrete(labels = paste0(filter(MSDfit_perPos.df, Dish == 'Dish1')$Pos, '\n',
#                                           filter(MSDfit_perPos.df, Dish == 'Dish1')$label)) +
#     labs(title = 'Dish1',
#          x = 'Log(Interval, s)', y = 'Log(MSD, um2)') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1) -> MSD_perPos_plot1
#   filter(MSD_perPos, Dish == 'Dish2') %>%
#     ggplot(aes(x = log10(interval*3/1000), y = log10(MSD), colour = Pos)) +
#     geom_point() +
#     #geom_line(size = 0.2) +
#     geom_abline(data = filter(MSDfit_perPos.df, Dish == 'Dish2'),
#                 aes(intercept = Intercept, slope = Slope, colour = Pos)) +
#     # geom_text(data = MSDfit_perPos.df,
#     #           aes(x = rep(-2.4, 6), y = seq(-1.3, -1.8, length.out = 6),
#     #               label = label, colour = Pos)) +
#     # annotate('text', x = -1.9, y = -1.8,
#     #          label = paste('D =', round(10^MSDfit.df$Intercept[[1]], 2),
#     #                        'um2/s,\nalpha =', round(MSDfit.df$Slope[[1]], 2))) +
#     scale_colour_discrete(labels = paste0(filter(MSDfit_perPos.df, Dish == 'Dish2')$Pos, '\n',
#                                           filter(MSDfit_perPos.df, Dish == 'Dish2')$label)) +
#     labs(title = 'Dish2',
#          x = 'Log(Interval, s)', y = 'Log(MSD, um2)') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1) -> MSD_perPos_plot2
#   ggarrange(MSD_perPos_plot1, MSD_perPos_plot2, nrow = 2) %>% 
#     annotate_figure(text_grob('MSD curve\nDens < 10 loc/fr, traj <15', size = 16))
#   ggsave('MSD_perPos_loglog1.png', width = 5.5, height = 7.5)
#   ggplot(MSD_perPos, aes(x = log10(interval*3/1000), y = log10(MSD),
#                          colour = Dish, group = interaction(Dish, Pos))) +
#     geom_point() +
#     geom_abline(data = MSDfit_perPos.df,
#                 aes(intercept = Intercept, slope = Slope,
#                     colour = Dish, group = interaction(Dish, Pos))) +
#     scale_colour_manual(values = c('blue', 'orange'),
#                         labels = paste0(MSDfit_perPos.summary$Dish, '\n',
#                                         MSDfit_perPos.summary$label)) +
#     labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15',
#          x = 'Log(Interval, s)', y = 'Log(MSD, um2)') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1)
#   ggsave('MSD_perPos_loglog2.png', width = 5.5, height = 3)
#   ggplot(MSDfit_perPos.df, aes(x = Dish, y = 10^Intercept)) +
#     geom_boxplot(width = 0.2) +
#     geom_jitter(width = 0.05, aes(colour = Pos)) +
#     labs(x = '', y = 'D, um2/s') +
#     theme_bw() +
#     theme(aspect.ratio = 3, legend.position = 'none') -> Ds
#   ggplot(MSDfit_perPos.df, aes(x = Dish, y = Slope)) +
#     geom_boxplot(width = 0.2) +
#     geom_jitter(width = 0.05, aes(colour = Pos)) +
#     labs(x = '', y = 'Alpha') +
#     theme_bw() +
#     theme(aspect.ratio = 3, legend.position = 'none') -> alphas
#   ggarrange(Ds, alphas, nrow = 1, ncol = 2)
#   ggsave('D_alpha.png', width = 3, height = 3)
#   
#   ### Excluding last points that deviate from linearity
#   split(filter(MSD, interval < 6), filter(MSD, interval < 6)$Dish) %>% lapply(function(x) {
#     lm(log10(MSD)~log10(interval*3/1000), data = x)
#   }) -> MSDfit2.list
#   MSDfit2.df = data.frame(Dish = names(MSDfit2.list),
#                           Intercept = sapply(MSDfit2.list, function(x) {x$coefficients[[1]]}),
#                           Slope = sapply(MSDfit2.list, function(x) {x$coefficients[[2]]}))
#   mutate(MSDfit2.df, label = paste('D =', round(10^Intercept, 2), 'um2/s, alpha =', round(Slope, 2))) -> MSDfit2.df
#   ggplot(MSD, aes(x = log10(interval*3/1000), y = log10(MSD), shape = interval<6, colour = Dish)) +
#     geom_point() +
#     scale_colour_manual(values = c('blue', 'orange'),
#                         labels = paste0(MSDfit2.df$Dish, '\n', MSDfit2.df$label)) +
#     scale_shape_manual(values = c(4, 19), guide = 'none') +
#     geom_abline(data = MSDfit2.df, aes(intercept = Intercept, slope = Slope, colour = Dish)) +
#     labs(title = 'MSD curve', subtitle = 'All Pos, dens < 10 loc/fr, traj < 15',
#          x = 'Log(Interval, s)', y = 'Log(MSD, um2)') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1)
#   ggsave('MSD_loglog_shortInts.png', width = 5, height = 3)
#   split(filter(MSD_perPos, interval < 6),
#         paste0(filter(MSD_perPos, interval < 6)$Dish, '_', filter(MSD_perPos, interval < 6)$Pos)) %>%
#     lapply(function(x) {
#       lm(log10(MSD)~log10(interval*3/1000), data = x)
#     }) -> MSDfit_perPos2.list
#   MSDfit_perPos2.df = data.frame(Dish = sapply(str_split(names(MSDfit_perPos2.list), '_'), function(x) x[[1]]),
#                                  Pos = sapply(str_split(names(MSDfit_perPos2.list), '_'), function(x) x[[2]]),
#                                  Intercept = sapply(MSDfit_perPos2.list, function(x) {x$coefficients[[1]]}),
#                                  Slope = sapply(MSDfit_perPos2.list, function(x) {x$coefficients[[2]]}))
#   mutate(MSDfit_perPos2.df, label = paste('D =', round(10^Intercept, 2), 'um2/s, alpha =', round(Slope, 2))) -> MSDfit_perPos2.df
#   MSDfit_perPos2.df %>% group_by(Dish) %>%
#     summarise(label = paste0('D = ', round(10^min(Intercept), 2), '-', round(10^max(Intercept), 2),
#                              ' um2/s, alpha = ', round(min(Slope), 2), '-', round(max(Slope), 2))) -> MSDfit_perPos2.summary
#   filter(MSD_perPos, Dish == 'Dish1') %>%
#     ggplot(aes(x = log10(interval*3/1000), y = log10(MSD), colour = Pos, shape = interval<6)) +
#     geom_point() +
#     scale_shape_manual(values = c(4, 19), guide = 'none') +
#     geom_abline(data = filter(MSDfit_perPos2.df, Dish == 'Dish1'),
#                 aes(intercept = Intercept, slope = Slope, colour = Pos)) +
#     scale_colour_discrete(labels = paste0(filter(MSDfit_perPos2.df, Dish == 'Dish1')$Pos, '\n',
#                                           filter(MSDfit_perPos2.df, Dish == 'Dish1')$label)) +
#     labs(title = 'Dish1',
#          x = 'Log(Interval, s)', y = 'Log(MSD, um2)') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1) -> MSD_perPos_plot1
#   filter(MSD_perPos, Dish == 'Dish2') %>%
#     ggplot(aes(x = log10(interval*3/1000), y = log10(MSD), colour = Pos, shape = interval<6)) +
#     geom_point() +
#     scale_shape_manual(values = c(4, 19), guide = 'none') +
#     geom_abline(data = filter(MSDfit_perPos2.df, Dish == 'Dish2'),
#                 aes(intercept = Intercept, slope = Slope, colour = Pos)) +
#     scale_colour_discrete(labels = paste0(filter(MSDfit_perPos2.df, Dish == 'Dish2')$Pos, '\n',
#                                           filter(MSDfit_perPos2.df, Dish == 'Dish2')$label)) +
#     labs(title = 'Dish2',
#          x = 'Log(Interval, s)', y = 'Log(MSD, um2)') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1) -> MSD_perPos_plot2
#   ggarrange(MSD_perPos_plot1, MSD_perPos_plot2, nrow = 2) %>% 
#     annotate_figure(text_grob('MSD curve\nDens < 10 loc/fr, traj <15', size = 16))
#   ggsave('MSD_perPos_loglog_shortInts_1.png', width = 5.5, height = 7.5)
#   ggplot(MSD_perPos, aes(x = log10(interval*3/1000), y = log10(MSD),
#                          colour = Dish, group = interaction(Dish, Pos), shape = interval<6)) +
#     geom_point() +
#     scale_shape_manual(values = c(4, 19), guide = 'none') +
#     geom_abline(data = MSDfit_perPos2.df,
#                 aes(intercept = Intercept, slope = Slope,
#                     colour = Dish, group = interaction(Dish, Pos))) +
#     scale_colour_manual(values = c('blue', 'orange'),
#                         labels = paste0(MSDfit_perPos2.summary$Dish, '\n',
#                                         MSDfit_perPos2.summary$label)) +
#     labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15',
#          x = 'Log(Interval, s)', y = 'Log(MSD, um2)') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1)
#   ggsave('MSD_perPos_loglog_shortInts_2.png', width = 5.5, height = 3)
#   ggplot(MSDfit_perPos2.df, aes(x = Dish, y = 10^Intercept)) +
#     geom_boxplot(width = 0.2) +
#     geom_jitter(width = 0.05, aes(colour = Pos)) +
#     labs(x = '', y = 'D, um2/s') +
#     theme_bw() +
#     theme(aspect.ratio = 3, legend.position = 'none') -> Ds
#   ggplot(MSDfit_perPos2.df, aes(x = Dish, y = Slope)) +
#     geom_boxplot(width = 0.2) +
#     geom_jitter(width = 0.05, aes(colour = Pos)) +
#     labs(x = '', y = 'Alpha') +
#     theme_bw() +
#     theme(aspect.ratio = 3, legend.position = 'none') -> alphas
#   ggarrange(Ds, alphas, nrow = 1, ncol = 2)
#   ggsave('D_alpha2.png', width = 3, height = 3)
# }

#### Define HP1 foci using DBSCAN; classify trajectories ####

##Segmentation and classification
{
grid_fun = function (x, y = NULL, x.breaks, y.breaks, same.scale = FALSE, na.rm = TRUE, 
                     FUN = base::length, arg = x, exclude_thresh = 5)
{
  if (is.null(y)) {
    if (ncol(x) != 2) 
      stop("If y is ommitted, x must be a 2 column matirx")
    y <- x[, 2]
    x <- x[, 1]
  }
  nas <- is.na(x) | is.na(y)
  if (na.rm) {
    x <- x[!nas]
    y <- y[!nas]
  }
  else stop("missinig values not permitted if na.rm=FALSE")
  if(same.scale){
    x.cuts = x.breaks;
    y.cuts = x.breaks;
  }else{
    x.cuts <- x.breaks
    y.cuts <- y.breaks   
  }
  
  
  index.x <- cut(x, x.cuts, include.lowest = TRUE)
  index.y <- cut(y, y.cuts, include.lowest = TRUE)
  c <- tapply(x, list(index.x, index.y), length)
  #c[is.na(c)] <- 0
  m <- tapply(arg, list(index.x, index.y), FUN)
  m_thresh = m
  m_thresh[c<exclude_thresh] = NA
  #if (identical(FUN, base::length)) 
  #  m[is.na(m)] <- 0
  midpoints <- function(x) (x[-1] + x[-length(x)])/2
  retval <- list()
  retval$result <- m
  retval$result_thresh <- m_thresh
  retval$counts <- c
  retval$counts_rel <- c/max(c)  
  retval$x.breaks = x.cuts
  retval$y.breaks = y.cuts
  retval$x = midpoints(x.cuts)
  retval$y = midpoints(y.cuts)
  retval$nobs = length(x)
  retval$bins = c(length(x.cuts),length(y.cuts))
  retval$call <- match.call()
  class(retval) <- "hist2d"
  retval
}

px_size = 100
plot_dir = 'New_analysis/PlotsPerCell_100nm/'
thresholds1 = c(100, 100, 800, 200, 100, 0, #Dish1 (0, 10, 11, 1, 2, 3)
                50, 1000, 100, 1000, 600, 200, #Dish1 (4, 5, 6, 7, 8, 9)
                800, 1200, 850, 500, 1000, 200, 1100, 800, 200) #Dish2 (1-8)
thresholds2 = c(19500, 11000, 19000, 20000, 20000, 10000, #Dish1 (0, 10, 11, 1, 2, 3)
                10000, 21000, 12000, 20000,11000, 21000, #Dish1 (4, 5, 6, 7, 8, 9)
                25000, 24000, 24000, 20000, 23000, 10000, 33000, 17000, 14000) #Dish2 (1-8)
n_shuffle = 6
segmCoresList = list()
for (c in 1:length(data_list)) {
  tic()
  print(names(data_list)[c])
  #Calculate result
  {
    filter(data_list_SMLM[[c]], Frame > thresholds1[c]) %>% group_by(Dish, Pos, Track) %>%
      summarise(x = mean(x), y = mean(y), Frame = mean(Frame)) -> data_SMLM_all
    histogram_all = grid_fun(x = data_SMLM_all$x, y = data_SMLM_all$y,
                             x.breaks = seq(0, max(data_SMLM_all$x)+px_size, by = px_size),
                             y.breaks = seq(0, max(data_SMLM_all$y)+px_size, by = px_size),
                             exclude_thresh = 1)
    histogram_perCycle = list()
    for (cc in 0:3) {
      temp = filter(data_SMLM_all, Frame > cc*10000 & Frame <= (cc+1)*10000)
      histogram_perCycle[[cc+1]] = grid_fun(x = temp$x, y = temp$y,
                                            x.breaks = seq(0, max(data_SMLM_all$x)+px_size, by = px_size),
                                            y.breaks = seq(0, max(data_SMLM_all$y)+px_size, by = px_size),
                                            exclude_thresh = 1)
    }
    filter(data_list[[c]], Frame > thresholds2[c]) %>%
      filter(Track_length > 1 & Track_length < 8) %>% group_by(Dish, Pos, Track) %>%
      mutate(diff.x = c(NA, diff(x)), diff.y = c(NA, diff(y)), sqdisp = diff.x^2 + diff.y^2) %>%
      summarise(Frame = mean(Frame), x = mean(x), y = mean(y), D = mean(sqdisp, na.rm = T)/(4000000*0.003)) -> tempD
    
    dispmap = grid_fun(x = tempD$x, y = tempD$y, FUN = mean, arg = tempD$D,
                       x.breaks = seq(0, max(data_SMLM_all$x)+px_size, by = px_size),
                       y.breaks = seq(0, max(data_SMLM_all$y)+px_size, by = px_size),
                       exclude_thresh = 5)
    result = data.frame(mids.x = rep(histogram_all$x, length(histogram_all$y)),
                        mids.y = rep(histogram_all$y, each = length(histogram_all$x)),
                        HP1_dens_all = c(histogram_all$counts),
                        cycle0 = c(histogram_perCycle[[1]]$counts),
                        cycle1 = c(histogram_perCycle[[2]]$counts),
                        cycle2 = c(histogram_perCycle[[3]]$counts),
                        cycle3 = c(histogram_perCycle[[4]]$counts),
                        D_dens = c(dispmap$counts),
                        D = c(dispmap$result),
                        D_thresh = c(dispmap$result_thresh))
  }
  
  #Images
  {
    ggplot(result, aes(x = mids.x, y = mids.y, fill = HP1_dens_all)) +
      geom_raster() +
      scale_fill_viridis_c(na.value = 'black', option = 'inferno', limits = c(0, px_size/2.5), oob = squish) +
      coord_fixed() +
      labs(title = paste('Histogram, bin =', px_size, 'nm'),
           subtitle = paste0('All data (frame>', thresholds1[c], '); r=300 nm, <3 gaps'),
           x = 'x', y = 'y', fill = 'Count') +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
            panel.grid = element_blank(), axis.line = element_blank()) -> image_all
    result %>% pivot_longer(cycle0:cycle3, names_to = 'cycle', values_to = 'HP1_dens') %>%
      group_by(cycle) %>% mutate(HP1_dens_norm = HP1_dens/mean(HP1_dens, na.rm = T)) %>%
      ggplot(aes(x = mids.x, y = mids.y, fill = HP1_dens_norm)) + #try making relative scale per facet
      geom_raster() +
      facet_wrap(~cycle, nrow = 2, ncol = 2) +
      scale_fill_viridis_c(na.value = 'black', option = 'inferno', limits = c(0, 5), oob = squish) +
      coord_fixed() +
      labs(title = paste('Histogram, bin =', px_size, 'nm'),
           subtitle = paste0('Per 10,000 frame "cycle" (cycle0 frame>', thresholds1[c], '); r=300 nm, <3 gaps'),
           x = 'x', y = 'y', fill = 'Count/mean(count)') +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
            panel.grid = element_blank(), axis.line = element_blank()) -> image_perCycle
  }
  
  #Disp maps
  {
    ggplot(result, aes(x = mids.x, y = mids.y, fill = D_dens)) +
      geom_raster() +
      scale_fill_gradientn(na.value = 'white', colours = c('gray', 'red', 'yellow'), limit = c(4, 20), oob = squish) +
      coord_fixed() +
      labs(title = paste('Histogram, bin =', px_size, 'nm'),
           subtitle = paste0('Frames ', thresholds2[c]+1,'-40000; r=400 nm, no gaps'),
           x = 'x', y = 'y', fill = 'Count') +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
            panel.grid = element_blank(), axis.line = element_blank()) -> data_density
    ggplot(result, aes(x = D_dens)) +
      geom_histogram(binwidth = 1) +
      lims(x = c(0, 20)) +
      labs(title = paste('Num traj per bin, bin =', px_size, 'nm'),
           subtitle = paste0('Frames ', thresholds2[c]+1,'-40000; r=400 nm, no gaps'),
           x = 'Num traj', y = 'Count') +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
            aspect.ratio = 0.5) -> data_density_hist
    ggplot(result, aes(x = mids.x, y = mids.y, fill = D)) +
      geom_raster() +
      scale_fill_viridis_c(na.value = 'white', limits = c(0, 8), oob = squish) +
      coord_fixed() +
      labs(title = paste('Diff coef, bin =', px_size, 'nm'),
           subtitle = paste0('D=mean(MSD_i/4dt), i=1...n traj within bin (no min!)\nFrames ', thresholds2[c]+1, '-40000; r=400 nm, no gaps'),
           x = 'x', y = 'y', fill = 'D, um2/s') +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
            panel.grid = element_blank(), axis.line = element_blank()) -> disp_plot
    ggplot(result, aes(x = mids.x, y = mids.y, fill = log10(D))) +
      geom_raster() +
      scale_fill_viridis_c(na.value = 'white', limits = c(-1, 1), oob = squish, 
                           labels = math_format(10^.x)) +
      coord_fixed() +
      labs(title = paste('Diff coef, bin =', px_size, 'nm'),
           subtitle = paste0('D=mean(MSD_i/4dt), i=1...n traj within bin (no min!)\nFrames ', thresholds2[c]+1, '-40000; r=400 nm, no gaps'),
           x = 'x', y = 'y', fill = 'D, um2/s') +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
            panel.grid = element_blank(), axis.line = element_blank()) -> disp_plot_log
  }
  
  #DBSCAN segmentation and JDA
  {
    eps = 50
    # minPts = round(1200000000/mean(filter(data_pdist, Dish == data_list[[c]]$Dish[[1]],
    #                                       Pos == data_list[[c]]$Pos[[1]])$pdist_mean)^2)
    minPts = round(6650/(mean(kNNdist(matrix(c(data_SMLM_all$x, data_SMLM_all$y), nrow=nrow(data_SMLM_all)),
                                      k = 1, all = F))^2))
    #empirical estimate:
    #--for c = 2, 35 points worked very well; also had a look at a few others
    #--(mean distance to nearest neighbour)^2 roughly scales with 1/(4*density), and density and minPts should have a linear relationship
    
    #segmentation
    segmCoresList[[c]] = list()
    
    segmentation = dbscan(x = select(ungroup(data_SMLM_all), x, y), eps = eps, minPts = minPts)
    data_SMLM_all$cluster = segmentation$cluster
    #data_SMLM_all$cluster_core = is.corepoint(x = select(ungroup(data_SMLM_all), x, y), eps = eps, minPts = minPts)
    
    ggplot(data_SMLM_all, aes(x = x, y = y, colour = cluster!=0)) +
      geom_point(size = 0.3, shape = 20, stroke = 0) +
      scale_colour_manual(values = c('gray', 'red'),
                          labels = c('Outside', 'Cluster')) +
      coord_fixed() +
      labs(colour = '', title = 'DBSCAN segmentation',
           subtitle = paste0('Eps = ', eps, ' nm, minPts = ', minPts, '\nAll data (frame>', thresholds1[c], '); r=300 nm, <3 gaps')) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) -> dbscan_segm
    
    ungroup(data_SMLM_all) %>% filter(cluster!=0) %>% select(x, y, cluster) -> clust_all
    segmCoresList[[c]][['Foci']] = clust_all
    segmentation = list('cluster' = clust_all$cluster,
                        'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
    class(segmentation) = c('dbscan_fast', 'dbscan')
    
    #pseudoclusters
    clust_all %>% group_by(cluster) %>%
      summarise(x = mean(x), y = mean(y)) -> clust_centres
    clust_all %>% group_by(cluster) %>% full_join(., ., by = 'cluster') %>%
      mutate(dist = sqrt((x.x-x.y)^2 + (y.x-y.y)^2)) %>% summarise(rad = max(dist)/2) %>%
      left_join(clust_centres, by = 'cluster') -> clust_centres
    rad_sums = matrix(rep(0, length(clust_centres$rad)^2),
                      nrow = length(clust_centres$rad))
    for (i in 1:length(clust_centres$rad)) {
      for (j in 1:length(clust_centres$rad)) {
        rad_sums[i, j] = clust_centres$rad[i] + clust_centres$rad[j]
      }
    }
    pseudosegmentation = list()
    for (n in 1:n_shuffle) {
      print(paste('Pseudo', n))
      is.bad = T
      loops = 0
      while(is.bad) {
        # filter(data_SMLM_all, cluster == 0) %>%
        #   mutate(pseudo_clust = sample(c(rep(0, nrow(.) - nrow(clust_centres)), 1:nrow(clust_centres)))) %>%
        #   filter(pseudo_clust != 0) -> pseudo_tmp
        filter(ungroup(data_SMLM_all), cluster == 0) %>%
          slice_sample(n = nrow(clust_centres)) -> pseudo_tmp
        #Check that pclust don't overlap with real clusters
        any(cdist(matrix(c(pseudo_tmp$x, pseudo_tmp$y), nrow = nrow(pseudo_tmp)),
                  matrix(c(clust_centres$x, clust_centres$y), nrow = nrow(clust_centres))) <
              rad_sums) -> is.bad1
        # any(apply(cdist(matrix(c(pseudo_tmp$x, pseudo_tmp$y), nrow = nrow(pseudo_tmp)),
        #             matrix(c(clust_centres$x, clust_centres$y), nrow = nrow(clust_centres))),
        #       2, min) < clust_centres$rad*2) -> is.bad1
        #Check that pclust don't overlap with each other
        pdist(matrix(c(pseudo_tmp$x, pseudo_tmp$y), nrow = nrow(pseudo_tmp))) -> dists_tmp
        diag(dists_tmp) = NA
        any(dists_tmp < rad_sums, na.rm = T) -> is.bad2
        is.bad = any(is.bad1, is.bad2)
        loops = loops+1
        if (loops > 2000) { #relax the overlap condition
          print('Relax')
          while(is.bad) {
            filter(ungroup(data_SMLM_all), cluster == 0) %>%
              slice_sample(n = nrow(clust_centres)) -> pseudo_tmp
            #Check that pclust don't overlap with real clusters
            any(cdist(matrix(c(pseudo_tmp$x, pseudo_tmp$y), nrow = nrow(pseudo_tmp)),
                      matrix(c(clust_centres$x, clust_centres$y), nrow = nrow(clust_centres))) <
                  rad_sums) -> is.bad1
            pdist(matrix(c(pseudo_tmp$x, pseudo_tmp$y), nrow = nrow(pseudo_tmp))) -> dists_tmp
            diag(dists_tmp) = NA
            any(dists_tmp < 0.5*rad_sums, na.rm = T) -> is.bad2
            is.bad = any(is.bad1, is.bad2)
            loops = loops+1
            if (loops > 5000) {#relax the overlap with real condition
              print('Relax')
              while(is.bad) {
                filter(ungroup(data_SMLM_all), cluster == 0) %>%
                  slice_sample(n = nrow(clust_centres)) -> pseudo_tmp
                #Check that pclust don't overlap with real clusters
                any(cdist(matrix(c(pseudo_tmp$x, pseudo_tmp$y), nrow = nrow(pseudo_tmp)),
                          matrix(c(clust_centres$x, clust_centres$y), nrow = nrow(clust_centres))) <
                      0.75*rad_sums) -> is.bad1
                pdist(matrix(c(pseudo_tmp$x, pseudo_tmp$y), nrow = nrow(pseudo_tmp))) -> dists_tmp
                diag(dists_tmp) = NA
                any(dists_tmp < 0.5*rad_sums, na.rm = T) -> is.bad2
                is.bad = any(is.bad1, is.bad2)
                loops = loops+1
                if (loops > 10000) { #relax both a bit further
                  print('Relax')
                  while(is.bad) {
                    filter(ungroup(data_SMLM_all), cluster == 0) %>%
                      slice_sample(n = nrow(clust_centres)) -> pseudo_tmp
                    #Check that pclust don't overlap with real clusters
                    any(cdist(matrix(c(pseudo_tmp$x, pseudo_tmp$y), nrow = nrow(pseudo_tmp)),
                              matrix(c(clust_centres$x, clust_centres$y), nrow = nrow(clust_centres))) <
                          0.7*rad_sums) -> is.bad1
                    pdist(matrix(c(pseudo_tmp$x, pseudo_tmp$y), nrow = nrow(pseudo_tmp))) -> dists_tmp
                    diag(dists_tmp) = NA
                    any(dists_tmp < 0.4*rad_sums, na.rm = T) -> is.bad2
                    is.bad = any(is.bad1, is.bad2)
                    loops = loops+1
                  }
                }
              }
            }
          }
        }
      }
      pseudo_tmp %>% select(-Frame, -Track) %>%
        mutate(cluster = 1:nrow(pseudo_tmp)) %>% rename(x_new = x, y_new = y) %>%
        full_join(clust_centres, by = 'cluster') %>% mutate(dx = x_new - x, dy = y_new - y) %>%
        select(cluster, dx, dy) %>% full_join(clust_all, by = 'cluster') %>%
        mutate(x = x+dx, y = y+dy) %>%
        select(x, y, cluster)-> pseudoclust_all
      segmCoresList[[c]][[paste0('Pseudo', n)]] = pseudoclust_all
      pseudosegmentation[[n]] = list('cluster' = pseudoclust_all$cluster,
                                     'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
      class(pseudosegmentation[[n]]) = c('dbscan_fast', 'dbscan')
      
      data_SMLM_all[[paste0('Pseudo', n)]] = predict(object = pseudosegmentation[[n]],
                                                     newdata = select(ungroup(data_SMLM_all), x, y),
                                                     data = select(ungroup(pseudoclust_all), x, y))
    }
    
    
    data_SMLM_all %>% pivot_longer(paste0('Pseudo', 1:n_shuffle),
                                   names_to = c('.value', 'Shuffle_rep'),
                                   names_pattern = '([A-Za-z_]+)(\\d+)') %>%
      ggplot(aes(x = x, y = y, colour = Pseudo != 0)) +
      geom_point(size = 0.3, shape = 20, stroke = 0) +
      facet_wrap(~Shuffle_rep, nrow = n_shuffle/2, ncol = 2) +
      scale_colour_manual(values = c('gray', 'black'),
                          labels = c('Outside', 'Pseudoclust')) +
      coord_fixed() +
      labs(colour = '', title = 'Shuffled DBSCAN segmentation',
           subtitle = paste0('Eps = ', eps, ' nm, minPts = ', minPts, '\nAll data (frame>', thresholds1[c], '); r=300 nm, <3 gaps')) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) -> dbscan_segm_shuffle
    
    #Calculate JDA, classify
    filter(data_list[[c]], Frame > thresholds2[c]) %>%
      filter(Track_length == 2) %>% group_by(Dish, Pos, Track) %>%
      summarise(disp = sqrt(diff(x)^2 + diff(y)^2), x = mean(x), y = mean(y)) -> result2
    
    result2$cluster = predict(object = segmentation,
                              newdata = select(ungroup(result2), x, y),
                              data = select(clust_all, x, y))
    for (n in 1:n_shuffle) {
      result2[[paste0('pseudo', n)]] = predict(object = pseudosegmentation[[n]],
                                               newdata = select(ungroup(result2), x, y),
                                               data = select(segmCoresList[[c]][[paste0('Pseudo', n)]], x, y))
    }
    
    filter(result2, cluster != 0) %>% nrow() -> Nfoci
    rbind(filter(result2, pseudo1 != 0) %>% mutate(Focus_group = 'Pseudo1'),
          filter(result2, pseudo2 != 0) %>% mutate(Focus_group = 'Pseudo2'),
          filter(result2, pseudo3 != 0) %>% mutate(Focus_group = 'Pseudo3'),
          filter(result2, pseudo4 != 0) %>% mutate(Focus_group = 'Pseudo4'),
          filter(result2, pseudo5 != 0) %>% mutate(Focus_group = 'Pseudo5'),
          filter(result2, pseudo6 != 0) %>% mutate(Focus_group = 'Pseudo6')) %>%
      group_by(Focus_group) %>% summarise(N = n()) -> Npfoci
    rbind(result2 %>% mutate(Focus_group = case_when(cluster != 0 ~ 'Focus', cluster == 0 ~ 'Outside')), #because the pseudo groups overlap
          filter(result2, pseudo1 != 0) %>% mutate(Focus_group = 'Pseudo1'),
          filter(result2, pseudo2 != 0) %>% mutate(Focus_group = 'Pseudo2'),
          filter(result2, pseudo3 != 0) %>% mutate(Focus_group = 'Pseudo3'),
          filter(result2, pseudo4 != 0) %>% mutate(Focus_group = 'Pseudo4'),
          filter(result2, pseudo5 != 0) %>% mutate(Focus_group = 'Pseudo5'),
          filter(result2, pseudo6 != 0) %>% mutate(Focus_group = 'Pseudo6')) %>%
      #mutate(Focus_group = factor(Focus_group, levels = c('Pseudo1', 'Pseudo2', 'Pseudo3', 'Pseudo4', 'Outside', 'Focus'))) %>%
      ggplot(aes(x = disp, colour = Focus_group, size = Focus_group)) +
      geom_density() +
      scale_colour_manual(values = c('Outside' = 'gray50', 'Focus' = 'red',
                                     'Pseudo1' = 'black', 'Pseudo2' = 'black',
                                     'Pseudo3' = 'black', 'Pseudo4' = 'black',
                                     'Pseudo5' = 'black', 'Pseudo6' = 'black')) +
      scale_size_manual(values = c('Focus' = 1, 'Outside' = 1,
                                   'Pseudo1' = 0.35, 'Pseudo2' = 0.35,
                                   'Pseudo3' = 0.35, 'Pseudo4' = 0.35,
                                   'Pseudo5' = 0.35, 'Pseudo6' = 0.35),
                        guide = 'none') +
      # scale_linetype_manual(values = c('Focus' = 'solid', 'Outside' = 'solid',
      #                                  'Pseudo1' = 'dashed', 'Pseudo2' = 'dashed',
      #                                  'Pseudo3' = 'dashed', 'Pseudo4' = 'dashed')) +
      labs(title = 'Jump distance distribution, DBSCAN segm', x = 'Displacement, nm',
           subtitle = paste0('Frames ', thresholds2[c]+1,'-40000 (r=400 nm, no gaps); traj of l = 2\n',
                             'Eps = ', eps, ' nm, minPts = ', minPts, '\n',
                             'Foci: ', Nfoci, ' jumps; Pseudo: ', min(Npfoci$N), '-', max(Npfoci$N), ' jumps'),
           colour = '') +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
            aspect.ratio = 0.5, axis.title.y = element_blank(), legend.position = 'top') -> JDA_dbscan
  }
  if (c == 1) {
    rbind(filter(result2, cluster != 0) %>% mutate(Focus_group = 'Focus'),
          filter(result2, pseudo1 != 0) %>% mutate(Focus_group = 'Pseudo1'),
          filter(result2, pseudo2 != 0) %>% mutate(Focus_group = 'Pseudo2'),
          filter(result2, pseudo3 != 0) %>% mutate(Focus_group = 'Pseudo3'),
          filter(result2, pseudo4 != 0) %>% mutate(Focus_group = 'Pseudo4'),
          filter(result2, pseudo5 != 0) %>% mutate(Focus_group = 'Pseudo5'),
          filter(result2, pseudo6 != 0) %>% mutate(Focus_group = 'Pseudo6')) -> data_segm_dbscan
  } else {
    rbind(filter(result2, cluster != 0) %>% mutate(Focus_group = 'Focus'),
          filter(result2, pseudo1 != 0) %>% mutate(Focus_group = 'Pseudo1'),
          filter(result2, pseudo2 != 0) %>% mutate(Focus_group = 'Pseudo2'),
          filter(result2, pseudo3 != 0) %>% mutate(Focus_group = 'Pseudo3'),
          filter(result2, pseudo4 != 0) %>% mutate(Focus_group = 'Pseudo4'),
          filter(result2, pseudo5 != 0) %>% mutate(Focus_group = 'Pseudo5'),
          filter(result2, pseudo6 != 0) %>% mutate(Focus_group = 'Pseudo6')) %>%
      rbind(data_segm_dbscan) -> data_segm_dbscan
  }
  
  ggarrange(image_all, image_perCycle, data_density, data_density_hist, disp_plot, disp_plot_log,
            dbscan_segm, dbscan_segm_shuffle, JDA_dbscan, ncol = 2, nrow = 5) -> plot
  ggsave(file = file.path(plot_dir, paste0(data_list[[c]]$Dish[[1]], '_',
                                           data_list[[c]]$Pos[[1]], '.png')),
         plot, width = 10, height = 25, dpi = 500)
  
  #Recording the segmentation into main data list
  {
    #DBscan
    #assign the SPT datapoints
    data_list[[c]]$HP1_foci_dbscan = predict(object = segmentation,
                                             newdata = select(ungroup(data_list[[c]]), x, y),
                                             data = select(clust_all, x, y))
    for (n in 1:n_shuffle) {
      data_list[[c]][[paste0('HP1_pseudo_dbscan', n)]] =
        predict(object = pseudosegmentation[[n]],
                newdata = select(ungroup(data_list[[c]]), x, y),
                data = select(segmCoresList[[c]][[paste0('Pseudo', n)]], x, y))
    }
  }
  toc()
}
names(segmCoresList) = names(data_list)
saveRDS(segmCoresList, 'HP1_dbscan.RDS')
saveRDS(data_list, 'Data_list_SPT_segm.RDS')
data = bind_rows(data_list)
}

##Alternatively, load pre-segmented data
segmCoresList = readRDS('HP1_dbscan.RDS') #list of cluster cores
data_list = readRDS('Data_list_SPT_segm.RDS') #classified data
data = bind_rows(data_list)

filter(data, Focus_group == 'Focus') %>% group_by(Dish) %>% summarise(N = n()) -> Nfoci
filter(data, Focus_group != 'Focus') %>% group_by(Dish) %>% group_split() %>%
  lapply(., function(x) x %>% group_by(Focus_group) %>% summarise(N = n())) -> Npfoci
ggplot(data, aes(x = disp, colour = Focus_group, size = Focus_group)) +
  geom_density() +
  facet_grid(Dish~.) +
  scale_size_manual(values = c(1, rep(0.35, n_shuffle)), guide = 'none') +
  scale_colour_manual(values = c('red', rep('black', n_shuffle))) +
  geom_density(data = jumpDist, colour = 'gray50', size = 1) +
  labs(title = 'Jump distance distribution', x = 'Displacement, nm',
       subtitle = paste('All Pos, dens < 10 loc/fr (r=400 nm, no gaps)\nTraj of l = 2\n',
                        'Eps = ', eps, ' nm, minPts varies\n',
                        'Foci:', Nfoci$N[1], '/', Nfoci$N[2], 'jumps;\nPseudo:',
                        min(Npfoci[[1]]$N), '-', max(Npfoci[[1]]$N), '/',
                        min(Npfoci[[2]]$N), '-', max(Npfoci[[2]]$N), ' jumps'),
       colour = '') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        aspect.ratio = 0.5, axis.title.y = element_blank())
ggsave('New_analysis/JumpDist_segm_dbscan.png', width = 4.5, height = 4.5)


#### Analyse DNA masks and create list of foci; classify trajectories ####
## Segmentation and classification
{
#Read in segmentation in DNA channel
WFframes = 8
WFaveraging = 5000 #WFframes*WFaveraging = total number of frames
pixel = 140
subsampling = 10

import = T #import the list of pre-segmented DNA foci?
read_new = F #analyse raw DNA masks to crease a list of DNA foci
dishes = paste0('Dish', 1:2)
#poses = paste0('Pos', 0:14)
if (import_foci) {
  foci_list = readRDS('DNA_foci_list.rds')
  pfoci_list = readRDS('DNA_pseudofoci_list.rds')
} else {
  foci_list = list()
  pfoci_list = list()
}
if (read_new) {
  for (dish in dishes) {
    numPos = length(list.files('TopHat_segmentation/PNGs_ManThreshold_coloured/', pattern = dish))/WFframes
    poses = paste0('Pos', 0:(numPos-1))
    for (pos in poses) {
      for (i in 1:WFframes) {
        #foci
        png = readPNG(file.path('TopHat_segmentation/PNGs_ManThreshold_coloured/',
                                paste0(dish, '_', pos, '_foci', paste0(rep(0, 2-nchar(as.character(i-1))), collapse = ''), i-1, '.png')))
        data.frame(px_id = 1:(dim(png)[1]*dim(png)[2]),
                   x_left = rep(0:(dim(png)[2]-1), each = dim(png)[1]),
                   y_top = rep(0:(dim(png)[1]-1), times = dim(png)[2]),
                   dim_y = dim(png)[1],
                   dim_x = dim(png)[2],
                   R = c(png[,,1]), G = c(png[,,2]), B = c(png[,,3]),
                   frame = i) %>%
          filter(R+G+B > 0) -> png_df
        if (i == 1) {
          tmp_df = png_df
        } else {
          tmp_df = rbind(tmp_df, png_df)
        }
        #pseudofoci
        png2 = readPNG(file.path('TopHat_segmentation/PNGs_pseudo',
                                 paste0(dish, '_', pos, '_foci', paste0(rep(0, 2-nchar(as.character(i-1))), collapse = ''), i-1, '.png')))
        data.frame(px_id = 1:(dim(png2)[1]*dim(png2)[2]),
                   x_left = rep(0:(dim(png2)[2]-1), each = dim(png2)[1]),
                   y_top = rep(0:(dim(png2)[1]-1), times = dim(png2)[2]),
                   dim_y = dim(png2)[1],
                   dim_x = dim(png2)[2],
                   #R = c(png2[,,1]), G = c(png2[,,2]), B = c(png2[,,3]),
                   pseudofocus = c(png2),
                   frame = i) %>%
          #group_by(R, G, B) %>% mutate(focus_id = cur_group_id()-1) %>% ungroup() %>% select(-R, -G, -B) %>%
          #filter(focus_id != 0) -> png2_df
          filter(pseudofocus != 0) -> png2_df
        if (i == 1) {
          tmp_df2 = png2_df
        } else {
          tmp_df2 = rbind(tmp_df2, png2_df)
        }
      }
      
      tmp_df %>% group_by(R, G, B) %>% mutate(focus_id = cur_group_id()) %>% #mutate(focus_id = cur_group_id()-1) %>%
        ungroup() %>% select(-R, -G, -B) %>%
        filter(focus_id != 0) -> tmp_df
      
      tmp_df$pos = pos
      tmp_df2$pos = pos
      tmp_df$dish = dish
      tmp_df2$dish = dish
      # #if (pos == unique(data$Pos)[1]) {
      # if (pos == poses[1]) {
      #   foci_df = tmp_df
      #   pfoci_df = tmp_df2
      # } else {
      #   foci_df = rbind(foci_df, tmp_df)
      #   pfoci_df = rbind(pfoci_df, tmp_df2)
      # }
      foci_list[[paste0(dish, '_', pos)]] = tmp_df
      pfoci_list[[paste0(dish, '_', pos)]] = tmp_df2
    }
    # saveRDS(foci_df, 'foci_df.rds')
    # saveRDS(pfoci_df, 'pseudofoci_df.rds')
    saveRDS(foci_list, 'DNA_foci_list.rds')
    saveRDS(pfoci_list, 'DNA_pseudofoci_list.rds')
  }
}
foci_all = bind_rows(foci_list, .id = 'Pos')
pfoci_all = bind_rows(pfoci_list, .id = 'Pos')

# Classify
classify = T
save_RDS = F
if (classify) {
  lapply(names(data_list), function(pos) {
    print(pos)
    df = data_list[[pos]]
    dim_x = foci_list[[pos]]$dim_x[1]
    dim_y = foci_list[[pos]]$dim_y[1]
    #x_left = floor(x*subsampling/pixel)
    #y_bottom = ceiling(y*subsampling/pixel)
    #px_id = x_left*dim_y + y_bottom
    mutate(df, px_id = floor(x*subsampling/pixel)*dim_y + ceiling(y*subsampling/pixel),
           frame = ceiling(Frame/WFaveraging)) -> df
    #mutate(focus_id = case_when(px_id %in% foci_list[[pos]]$px_id ~ 'focus',#foci_list[[pos]]
    #                                T ~ '0'))
    left_join(df, foci_list[[pos]], by = c('px_id', 'frame')) %>%
      #select(Pos, Track, Track_length, Frame, x, y, px_id, focus_id, frame) %>%
      select(-x_left, -y_top, -dim_y, -dim_x, -pos, -dish) %>% rename(DNA_focus_idx = focus_id) %>%
      mutate(DNA_focus_idx = replace_na(DNA_focus_idx, 0)) -> result #%>%
      #group_by(Track) %>% mutate(track_in_focus = mean(focus_id != 0)) -> result
    left_join(result, pfoci_list[[pos]], by = c('px_id', 'frame')) %>%
      #select(Pos, Track, Track_length, Frame, x, y, px_id, focus_id, track_in_focus, pseudofocus) %>%
      select(-x_left, -y_top, -dim_y, -dim_x, -frame, -pos, -dish) %>%
      rename(DNA_px_id = px_id) %>% mutate(DNA_pseudo = (replace_na(pseudofocus, 0) != 0)) %>%
      select(-pseudofocus) -> result#%>%
      #group_by(Track) %>% mutate(track_in_pfocus = mean(pseudofocus != 0)) -> result
    return(result)
  }) -> data_list
  names(data_list) = sapply(data_list, function(x) paste0(x$Dish[1], '_', x$Pos[1]))
  if (save_RDS) {
    saveRDS(data_list, 'Data_list_SPT_segm.RDS')
  }
}
#mutate(data, focus = case_when(focus_id == 0 ~ 0, T ~ 1)) -> data

#Sanity check
#Dish1_Pos2 - I actually duplicated one focus in pseudofoci by accident
# ggplot(data_list[[8]], aes(x = x, y = y, colour = as.character(DNA_focus_idx))) +
#   geom_point(size = 0.1) +
#   scale_y_reverse() +
#   coord_fixed() +
#   theme_bw()
# ggplot(data_list[[8]], aes(x = x, y = y, colour = DNA_pseudo)) +
#   geom_point(size = 0.1) +
#   scale_y_reverse() +
#   coord_fixed() +
#   theme_bw()
# ggplot(data_list[[8]], aes(x = x, y = y, colour = Track, group = Track)) +
#   geom_point(size = 0.1) +
#   geom_path() +
#   scale_y_reverse() +
#   coord_fixed() +
#   theme_bw()
}

##Alternatively, load pre-segmented data
foci_list = readRDS('DNA_foci_list.rds')
pfoci_list = readRDS('DNA_pseudofoci_list.rds')
data_list = readRDS('Data_list_SPT_segm.RDS')


bind_rows(data_list) -> data

## Jump distance analysis - classifying each jump according to its midpoint
## NB! Requires the segmentation in DNA channel
{
thresholds2 = c(19500, 11000, 19000, 20000, 20000, 10000, #Dish1 (0, 10, 11, 1, 2, 3)
                10000, 21000, 12000, 20000,11000, 21000, #Dish1 (4, 5, 6, 7, 8, 9)
                25000, 24000, 24000, 20000, 23000, 10000, 33000, 17000, 14000) #Dish2 (1-8)

lapply(seq_along(data_list), function(c) {
  pos = names(data_list)[c]
  df = filter(data_list[[pos]], Frame > thresholds2[c])
  dim_x = foci_list[[pos]]$dim_x[1]
  dim_y = foci_list[[pos]]$dim_y[1]
  
  df %>% group_by(Dish, Pos, Track) %>% filter(n() == 2) %>%
    summarise(disp = sqrt(diff(x)^2 + diff(y)^2), x = mean(x), y = mean(y), Frame1 = Frame[1]) %>%
    mutate(px_id = floor(x*subsampling/pixel)*dim_y + ceiling(y*subsampling/pixel),
           frame = ceiling(Frame1/WFaveraging)) %>%
    left_join(., foci_list[[pos]], by = c('px_id', 'frame')) %>%
    select(-x_left, -y_top, -dim_y, -dim_x, -pos) %>%
    left_join(., pfoci_list[[pos]], by = c('px_id', 'frame')) %>%
    select(-x_left, -y_top, -dim_y, -dim_x, -pos) %>%
    rename(DNA_focus_idx = focus_id, DNA_pseudo = pseudofocus) %>%
    mutate(DNA_focus_idx = replace_na(DNA_focus_idx, 0),
           DNA_pseudo = replace_na(DNA_pseudo, 0)) -> result
  return(result)
}) %>% bind_rows() -> jumpDist2
split(jumpDist2, paste0(jumpDist2$Dish, '_', jumpDist2$Pos)) -> jumpDist2_list

filter(jumpDist2, DNA_focus_idx != 0) %>% group_by(Dish) %>% summarise(N = n()) -> Nfoci
filter(jumpDist2, DNA_pseudo != 0) %>% group_by(Dish) %>% summarise(N = n()) -> Npfoci
ggplot(jumpDist2, aes(x = disp, colour = interaction((DNA_focus_idx != 0), DNA_pseudo != 0))) +
  geom_density(size = 1) +
  facet_grid(Dish~.) +
  scale_colour_manual(values = c('gray50', 'cyan', 'black'),
                      labels = c('Euchr', 'DNA foci', 'DNA pseudo')) +
  labs(title = 'Jump distance', x = 'Displacement, nm', y = '', colour = '',
       subtitle = paste0('Foci: ', Nfoci$N[1], '/', Nfoci$N[2],
                         ' jumps;\npseudo: ', Npfoci$N[1], '/', Npfoci$N[2], ' jumps')) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 0.5, axis.title.y = element_blank(),
        plot.subtitle = element_text(hjust = 0.5))
ggsave('New_analysis/JumpDist_EuHet.png', width = 4.7, height = 4.5)

#Per cell
{
# plot_dir = 'New_analysis/PlotsPerCell_DNA/'
# thresholds2 = c(19500, 11000, 19000, 20000, 20000, 10000, #Dish1 (0, 10, 11, 1, 2, 3)
#                 10000, 21000, 12000, 20000,11000, 21000, #Dish1 (4, 5, 6, 7, 8, 9)
#                 25000, 24000, 24000, 20000, 23000, 10000, 33000, 17000, 14000) #Dish2 (1-8)
# for (c in 1:length(data_list)) {
#   pos = names(data_list)[c]
#   ggplot(data_list[[c]], aes(x = x, y = y, colour = as.character(DNA_focus_idx))) +
#     geom_point(size = 0.3, shape = 20, stroke = 0) +
#     scale_y_reverse() +
#     scale_colour_manual(values = c('gray',
#                                    brewer_pal(palette = 'Set2')(length(unique(data_list[[c]]$DNA_focus_idx))-1))) +
#     coord_fixed() +
#     theme_bw() +
#     theme(legend.position = 'none') -> DNA_segm
#   ggplot(data_list[[c]], aes(x = x, y = y, colour = as.character(DNA_pseudo))) +
#     geom_point(size = 0.3, shape = 20, stroke = 0) +
#     scale_y_reverse() +
#     scale_colour_manual(values = c('gray', 'black')) +
#     coord_fixed() +
#     theme_bw() +
#     theme(legend.position = 'none') -> DNA_psegm
#   ggplot(jumpDist2_list[[pos]],
#          aes(x = disp, colour = interaction((DNA_focus_idx != 0), DNA_pseudo != 0))) +
#     geom_density() +
#     scale_colour_manual(values = c('gray50', 'cyan', 'black'),
#                         labels = c('Euchr', 'DNA foci', 'DNA pseudo')) +
#     labs(title = 'Jump distance',
#          subtitle = paste0('N_foc = ', sum(jumpDist2_list[[pos]]$DNA_focus_idx != 0),
#                            '; N_pseud = ', sum(jumpDist2_list[[pos]]$DNA_pseudo != 0)),
#          x = 'Displacement, nm', y = 'Density', colour = 'Compartment') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), aspect.ratio = 0.5) -> JDA_segm
#   ggplot(jumpDist2_list[[pos]],
#          aes(x = disp, colour = interaction(DNA_focus_idx, DNA_pseudo))) +
#     geom_density() +
#     scale_colour_manual(values = c('black', 'gray', brewer_pal(palette = 'Set2')(length(unique(data_list[[c]]$DNA_focus_idx))-1)),
#                         labels = c('Pseudofoci', 'Rest of euchr', paste0('Focus ', 1:(length(unique(data_list[[c]]$DNA_focus_idx))-1)))) +
#     #scale_colour_manual(values = c('cyan', 'black', 'gray')) +
#     labs(title = 'Jump distance', colour = '',
#          x = 'Displacement, nm', y = 'Density') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), aspect.ratio = 0.5) -> JDA_perFocus
#   ggarrange(DNA_segm, DNA_psegm, JDA_segm, JDA_perFocus, nrow = 2, ncol = 2) -> plot
#   ggsave(file = file.path(plot_dir, paste0(data_list[[c]]$Dish[[1]], '_', data_list[[c]]$Pos[[1]], '.png')),
#          plot, width = 10, height = 6, dpi = 500)
# }
}
}

#### Place HP1 pseudofoci in DNA foci - not included in paper ####

# data_list = readRDS('Data_list_SPT_segm.RDS') #NB! Pay attention if this has HP1 segm, DNA segm or both!
# data = bind_rows(data_list)
# 
# data_list_SMLM = readRDS('Data_list_SMLM.RDS')
# #names(data_list_SMLM) = sapply(data_list_SMLM, function(x) paste0(x$Dish[1], '_', x$Pos[1]))
# 
# WFframes = 8
# WFaveraging = 5000 #WFframes*WFaveraging = total number of frames
# pixel = 140
# subsampling = 10
# foci_list = readRDS('DNA_foci_list.rds')
# pfoci_list = readRDS('DNA_pseudofoci_list.rds')
# segmCoresList = readRDS('HP1_dbscan.RDS')
# 
# n_shuffle = 6
# thresholds1 = c(100, 100, 800, 200, 100, 0, #Dish1 (0, 10, 11, 1, 2, 3)
#                 50, 1000, 100, 1000, 600, 200, #Dish1 (4, 5, 6, 7, 8, 9)
#                 800, 1200, 850, 500, 1000, 200, 1100, 800, 200) #Dish2 (1-8)
# thresholds2 = c(19500, 11000, 19000, 20000, 20000, 10000, #Dish1 (0, 10, 11, 1, 2, 3)
#                 10000, 21000, 12000, 20000,11000, 21000, #Dish1 (4, 5, 6, 7, 8, 9)
#                 25000, 24000, 24000, 20000, 23000, 10000, 33000, 17000, 14000) #Dish2 (1-8)
# 
# #Version with proper SMLM data
# psegmCoresListInDNA = list()
# for(c in 1:length(data_list_SMLM)) {
#   pos = names(data_list_SMLM)[c]
#   print(pos)
#   tic()
#   psegmCoresListInDNA[[pos]] = list()
#   filter(data_list_SMLM[[c]], Frame > thresholds1[c]) %>% group_by(Dish, Pos, Track) %>%
#     summarise(x = mean(x), y = mean(y), Frame = mean(Frame)) -> data_SMLM_all
#   data_list[[pos]] -> data_tmp
#   #DBSCAN HP1
#   segmentation = list('cluster' = segmCoresList[[pos]][['Foci']]$cluster,
#                       'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
#   class(segmentation) = c('dbscan_fast', 'dbscan')
#   data_SMLM_all$cluster = predict(object = segmentation,
#                                   newdata = select(ungroup(data_SMLM_all), x, y),
#                                   data = select(ungroup(segmCoresList[[pos]][['Foci']]), x, y))
#   #DNA segmentation
#   {
#     dim_x = foci_list[[pos]]$dim_x[1]
#     dim_y = foci_list[[pos]]$dim_y[1]
#     mutate(data_SMLM_all, px_id = floor(x*subsampling/pixel)*dim_y + ceiling(y*subsampling/pixel),
#            frame = ceiling(Frame/WFaveraging)) -> data_SMLM_all
#     left_join(data_SMLM_all, foci_list[[pos]], by = c('px_id', 'frame')) %>%
#       select(-x_left, -y_top, -dim_y, -dim_x, -pos, -dish) %>% rename(DNA_focus_idx = focus_id) %>%
#       mutate(DNA_focus_idx = replace_na(DNA_focus_idx, 0)) -> data_SMLM_all
#     left_join(data_SMLM_all, pfoci_list[[pos]], by = c('px_id', 'frame')) %>%
#       select(-x_left, -y_top, -dim_y, -dim_x, -frame, -pos, -dish) %>%
#       rename(DNA_px_id = px_id, DNA_pseudo = pseudofocus) %>%
#       mutate(DNA_pseudo = (replace_na(DNA_pseudo, 0) != 0)) -> data_SMLM_all
#     
#     # ggplot(data_SMLM_all, aes(x = x, y = y, colour = interaction(DNA_focus, DNA_pseudo))) +
#     #   geom_point(size = 0.3, shape = 20, stroke = 0) +
#     #   scale_colour_manual(values = c('gray', 'cyan', 'blue'),
#     #                       labels = c('Outside', 'Cluster', 'Pseudo')) +
#     #   coord_fixed() +
#     #   labs(colour = '', title = 'DNA segmentation',
#     #        subtitle = paste0('All data (frame>', thresholds1[c], '); r=300 nm, <3 gaps')) +
#     #   theme_bw() +
#     #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   }
#   #HP1 pseudoclusters in DNA foci
#   {
#     segmCoresList[[pos]][['Foci']] %>% group_by(cluster) %>%
#       summarise(x = mean(x), y = mean(y)) -> clust_centres
#     segmCoresList[[pos]][['Foci']] %>% group_by(cluster) %>% full_join(., ., by = 'cluster') %>%
#       mutate(dist = sqrt((x.x-x.y)^2 + (y.x-y.y)^2)) %>% summarise(rad = max(dist)/2) %>%
#       left_join(clust_centres, by = 'cluster') -> clust_centres
#     pseudosegmentationInDNA = list()
#     for (n in 1:n_shuffle) {
#       print(n)
#       num_to_pick = min(3*(length(unique(data_SMLM_all$DNA_focus_idx))-1),
#                         nrow(clust_centres)) #not too many foci, otherwise they'll overlap
#       picked_centres = slice_sample(clust_centres, n = num_to_pick)
#       picked_clusters = filter(segmCoresList[[pos]][['Foci']], cluster %in% picked_centres$cluster)
#       picked_rad_sums = matrix(rep(0, length(picked_centres$rad)^2),
#                                nrow = length(picked_centres$rad))
#       for (i in 1:length(picked_centres$rad)) {
#         for (j in 1:length(picked_centres$rad)) {
#           picked_rad_sums[i, j] = picked_centres$rad[i] + picked_centres$rad[j]
#         }
#       }
#       is.bad = T
#       loops = 0
#       while(is.bad) {
#         filter(ungroup(data_SMLM_all), DNA_focus_idx!=0) %>%
#           slice_sample(n = num_to_pick) -> pseudo_tmp
#         #I don't care if these pseudofoci overlap with real foci
#         # any(apply(cdist(matrix(c(pseudo_tmp$x, pseudo_tmp$y), nrow = num_to_pick),
#         #                 matrix(c(picked_centres$x, picked_centres$y), nrow = num_to_pick)),
#         #           2, min) < picked_centres$rad) -> is.bad1
#         #Check that the pseudofoci don't overlap with each other
#         pdist(matrix(c(pseudo_tmp$x, pseudo_tmp$y), nrow = num_to_pick)) -> dists_tmp
#         diag(dists_tmp) = NA
#         any(dists_tmp < picked_rad_sums, na.rm = T) -> is.bad
#         loops = loops+1
#         if (loops == 1000) {
#           picked_centres = slice_sample(clust_centres, n = num_to_pick)
#           picked_clusters = filter(segmCoresList[[pos]][['Foci']], cluster %in% picked_centres$cluster)
#           picked_rad_sums = matrix(rep(0, length(picked_centres$rad)^2),
#                                    nrow = length(picked_centres$rad))
#         }
#         if (loops > 2000) {
#           print('Reduce num_to_pick')
#           num_to_pick = round(min(2.5*(length(unique(data_SMLM_all$DNA_focus_idx))-1),
#                                   nrow(clust_centres))) #not too many foci, otherwise they'll overlap
#           picked_centres = slice_sample(clust_centres, n = num_to_pick)
#           picked_clusters = filter(segmCoresList[[pos]][['Foci']], cluster %in% picked_centres$cluster)
#           picked_rad_sums = matrix(rep(0, length(picked_centres$rad)^2),
#                                    nrow = length(picked_centres$rad))
#           for (i in 1:length(picked_centres$rad)) {
#             for (j in 1:length(picked_centres$rad)) {
#               picked_rad_sums[i, j] = picked_centres$rad[i] + picked_centres$rad[j]
#             }
#           }
#           is.bad = T
#           loops = 0
#           while(is.bad) {
#             filter(ungroup(data_SMLM_all), DNA_focus_idx!=0) %>%
#               slice_sample(n = num_to_pick) -> pseudo_tmp
#             #I don't care if these pseudofoci overlap with real foci
#             # any(apply(cdist(matrix(c(pseudo_tmp$x, pseudo_tmp$y), nrow = num_to_pick),
#             #                 matrix(c(picked_centres$x, picked_centres$y), nrow = num_to_pick)),
#             #           2, min) < picked_centres$rad) -> is.bad1
#             #Check that the pseudofoci don't overlap with each other
#             pdist(matrix(c(pseudo_tmp$x, pseudo_tmp$y), nrow = num_to_pick)) -> dists_tmp
#             diag(dists_tmp) = NA
#             any(dists_tmp < picked_rad_sums, na.rm = T) -> is.bad
#             loops = loops+1
#             if (loops == 1000) {
#               picked_centres = slice_sample(clust_centres, n = num_to_pick)
#               picked_clusters = filter(segmCoresList[[pos]][['Foci']], cluster %in% picked_centres$cluster)
#               picked_rad_sums = matrix(rep(0, length(picked_centres$rad)^2),
#                                        nrow = length(picked_centres$rad))
#             }
#             if (loops > 2000) {
#               print('Reduce num_to_pick')
#               num_to_pick = min(2*(length(unique(data_SMLM_all$DNA_focus_idx))-1),
#                                 nrow(clust_centres)) #not too many foci, otherwise they'll overlap
#               picked_centres = slice_sample(clust_centres, n = num_to_pick)
#               picked_clusters = filter(segmCoresList[[pos]][['Foci']], cluster %in% picked_centres$cluster)
#               picked_rad_sums = matrix(rep(0, length(picked_centres$rad)^2),
#                                        nrow = length(picked_centres$rad))
#               for (i in 1:length(picked_centres$rad)) {
#                 for (j in 1:length(picked_centres$rad)) {
#                   picked_rad_sums[i, j] = picked_centres$rad[i] + picked_centres$rad[j]
#                 }
#               }
#               is.bad = T
#               loops = 0
#               while(is.bad) {
#                 filter(ungroup(data_SMLM_all), DNA_focus_idx!=0) %>%
#                   slice_sample(n = num_to_pick) -> pseudo_tmp
#                 #I don't care if these pseudofoci overlap with real foci
#                 # any(apply(cdist(matrix(c(pseudo_tmp$x, pseudo_tmp$y), nrow = num_to_pick),
#                 #                 matrix(c(picked_centres$x, picked_centres$y), nrow = num_to_pick)),
#                 #           2, min) < picked_centres$rad) -> is.bad1
#                 #Check that the pseudofoci don't overlap with each other
#                 pdist(matrix(c(pseudo_tmp$x, pseudo_tmp$y), nrow = num_to_pick)) -> dists_tmp
#                 diag(dists_tmp) = NA
#                 any(dists_tmp < picked_rad_sums, na.rm = T) -> is.bad
#                 loops = loops+1
#                 if (loops == 1000) {
#                   picked_centres = slice_sample(clust_centres, n = num_to_pick)
#                   picked_clusters = filter(segmCoresList[[pos]][['Foci']], cluster %in% picked_centres$cluster)
#                   picked_rad_sums = matrix(rep(0, length(picked_centres$rad)^2),
#                                            nrow = length(picked_centres$rad))
#                 }
#               }
#             }
#           }
#         }
#       }
#       pseudo_tmp %>% select(-Frame, -Track) %>%
#         mutate(cluster = picked_centres$cluster) %>% rename(x_new = x, y_new = y) %>%
#         left_join(clust_centres, by = 'cluster') %>% mutate(dx = x_new - x, dy = y_new - y) %>%
#         select(cluster, dx, dy) %>% left_join(picked_clusters, by = 'cluster') %>%
#         mutate(x = x+dx, y = y+dy) %>%
#         select(x, y, cluster)-> pseudoclust_all
#       psegmCoresListInDNA[[c]][[paste0('Pseudo', n)]] = pseudoclust_all
#       pseudosegmentationInDNA[[n]] = list('cluster' = pseudoclust_all$cluster,
#                                           'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
#       class(pseudosegmentationInDNA[[n]]) = c('dbscan_fast', 'dbscan')
#       
#       data_SMLM_all[[paste0('Pseudo_in_DNA', n)]] = predict(object = pseudosegmentationInDNA[[n]],
#                                                             newdata = select(ungroup(data_SMLM_all), x, y),
#                                                             data = select(ungroup(pseudoclust_all), x, y))
#       # ggplot(data_SMLM_all, aes(x = x, y = y, colour = interaction(DNA_focus_idx!=0, Pseudo_in_DNA1 != 0))) +
#       #   geom_point(size = 0.1) +
#       #   scale_colour_manual(values = c('gray', 'cyan', 'red', 'purple'),
#       #                       labels = c('Outside', 'DNA cluster', 'HP1 pseudo', 'Intersection')) +
#       #   coord_fixed() +
#       #   labs(colour = '', title = 'HP1 pseudo in DNA foci') +
#       #   theme_bw() +
#       #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#     }
#   }
#   #Recording the segmentation into main data list and plotting
#   {
#     for (n in 1:n_shuffle) {
#       data_list[[c]][[paste0('Pseudo_in_DNA', n)]] =
#         predict(object = pseudosegmentationInDNA[[n]],
#                 newdata = select(ungroup(data_list[[c]]), x, y),
#                 data = select(psegmCoresListInDNA[[pos]][[paste0('Pseudo', n)]], x, y))
#     }
#     
#     data_list[[c]] %>% pivot_longer('Pseudo_in_DNA1':paste0('Pseudo_in_DNA', n_shuffle),
#                                     names_to = c('.value', 'Shuffle_rep'),
#                                     names_pattern = '([A-Za-z_]+)(\\d+)') %>%
#       ggplot(aes(x = x, y = y, colour = interaction(DNA_focus_idx != 0, Pseudo_in_DNA != 0))) +
#       geom_point(size = 0.1, shape = 20, stroke = 0) +
#       facet_wrap(~Shuffle_rep, nrow = n_shuffle/2, ncol = 2) +
#       scale_colour_manual(values = c('gray', 'cyan', 'red', 'purple'),
#                           labels = c('Outside', 'DNA cluster', 'HP1 pseudo', 'Intersection')) +
#       coord_fixed() +
#       labs(colour = '', title = 'HP1 pseudo in DNA foci',
#            subtitle = paste0('All data (frame>', thresholds1[c], '); r=300 nm, <3 gaps')) +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) -> pseudo_in_DNA
#     ggsave(pseudo_in_DNA, file = paste0('New_analysis/PlotsPerCell_DNA/', pos, '_pseudoInDNA.png'),
#            width = 5, height = 6)
#     #Remove non-overlapping regions
#     for (n in 1:n_shuffle) {
#       data_list[[c]] %>% mutate(!!paste0('Pseudo_in_DNA', n) :=
#                                   case_when(DNA_focus_idx == 0 ~ as.integer(0),
#                                             T ~ !!as.name(paste0('Pseudo_in_DNA', n)))) -> data_list[[c]]
#     }
#   }
#   toc()
# }
# saveRDS(psegmCoresListInDNA, 'HP1_dbscan_pseudo_in_DNA.RDS')
# saveRDS(data_list, 'Data_list_SPT_segm.RDS')
# 
# # lapply(names(data_list), function(pos) {
# #   left_join(data_list[[pos]], select(data_list_old[[pos]], Track, Frame, x, y,
# #                                      HP1_dens_all, HP1_foci_idx, HP1_pseudo1_idx, HP1_pseudo2_idx,
# #                                      HP1_pseudo3_idx, HP1_pseudo4_idx, HP1_pseudo5_idx, HP1_pseudo6_idx),
# #             by = c('Track', 'Frame', 'x', 'y')) -> result
# # }) -> data_list_new
# # names(data_list_new) = names(data_list)
# # saveRDS(data_list_new, 'Data_list_SPT_segm.RDS')

#### Classification based on both DNA and HP1 - finding overlapping foci ####
##Read in data
import_segm = T
if (import_segm) {
  data_list = readRDS('Data_list_SPT_segm.RDS')
}
data = bind_rows(data_list)

foci_list = readRDS('DNA_foci_list.rds')
pfoci_list = readRDS('DNA_pseudofoci_list.rds')
segmCoresList = readRDS('HP1_dbscan.RDS')
psegmCoresListInDNA = readRDS('HP1_dbscan_pseudo_in_DNA.RDS')


##Classify HP1 clusters as overlapping/non-overlapping with DNA foci as a whole based on %overlap
{
  lapply(data_list, function(dat) {
    dat %>% group_by(HP1_foci_dbscan) %>%
      mutate(DNA_HP1_focus = case_when(HP1_foci_dbscan == 0 ~ F,
                                       sum(DNA_focus_idx!=0)/n() > 0.3 ~ T,
                                       T ~ F)) %>% ungroup() -> result
    # ggplot(result, aes(x = x, y = y, colour = interaction(DNA_focus_idx != 0, HP1_foci_dbscan != 0))) +
    #   geom_point(shape = 16, size = 0.01) +
    #   coord_fixed() +
    #   scale_colour_manual(values = c('gray', 'red', 'blue', 'orange')) +
    #   theme_bw() +
    #   theme(legend.position = 'none')
    # ggplot(result, aes(x = x, y = y, colour = DNA_HP1_focus)) +
    #   geom_point(shape = 16, size = 0.01) +
    #   coord_fixed() +
    #   scale_colour_manual(values = c('gray', 'magenta')) +
    #   theme_bw() +
    #   theme(legend.position = 'none')
    return(unique((result %>% filter(DNA_HP1_focus))$HP1_foci_dbscan))
  }) -> DNA_HP1_foci
  names(DNA_HP1_foci) = names(data_list)
}

##Jump distance analysis - not included in the paper
{
  # thresholds2 = c(19500, 11000, 19000, 20000, 20000, 10000, #Dish1 (0, 10, 11, 1, 2, 3)
  #                 10000, 21000, 12000, 20000,11000, 21000, #Dish1 (4, 5, 6, 7, 8, 9)
  #                 25000, 24000, 24000, 20000, 23000, 10000, 33000, 17000, 14000) #Dish2 (1-8)WFframes = 8
  # WFaveraging = 5000 #WFframes*WFaveraging = total number of frames
  # pixel = 140
  # subsampling = 10
  # n_shuffle = 6
  # 
  # lapply(seq_along(data_list), function(c) {
  #   tic()
  #   pos = names(data_list)[c]
  #   print(pos)
  #   df = select(filter(data_list[[pos]], Frame > thresholds2[c]),
  #               Dish, Pos, Track, Frame, x, y, Track_length)
  #   
  #   #Calculate JD, classify by DNA
  #   dim_x = foci_list[[pos]]$dim_x[1]
  #   dim_y = foci_list[[pos]]$dim_y[1]
  #   df %>% group_by(Dish, Pos, Track) %>% filter(n() == 2) %>%
  #     summarise(disp = sqrt(diff(x)^2 + diff(y)^2), x = mean(x), y = mean(y), Frame1 = Frame[1]) %>%
  #     mutate(px_id = floor(x*subsampling/pixel)*dim_y + ceiling(y*subsampling/pixel),
  #            frame = ceiling(Frame1/WFaveraging)) %>%
  #     left_join(., foci_list[[pos]], by = c('px_id', 'frame')) %>%
  #     select(-x_left, -y_top, -dim_y, -dim_x, -pos, -dish) %>%
  #     left_join(., pfoci_list[[pos]], by = c('px_id', 'frame')) %>%
  #     select(-x_left, -y_top, -dim_y, -dim_x, -pos, -dish, -frame, -px_id) %>%
  #     rename(DNA_focus_idx = focus_id, DNA_pseudo = pseudofocus) %>%
  #     mutate(DNA_focus_idx = replace_na(DNA_focus_idx, 0),
  #            DNA_pseudo = replace_na(DNA_pseudo, 0)) -> result
  #   
  #   #DNA foci, subsampled in HP1-foci-like manner
  #   for (n in 1:n_shuffle) {
  #     pseudosegmentation = list('cluster' = psegmCoresListInDNA[[pos]][[paste0('Pseudo', n)]]$cluster,
  #                               'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
  #     class(pseudosegmentation) = c('dbscan_fast', 'dbscan')
  #     result[[paste0('Pseudo_in_DNA', n)]] =
  #       predict(object = pseudosegmentation,
  #               newdata = select(ungroup(result), x, y),
  #               data = select(ungroup(psegmCoresListInDNA[[pos]][[paste0('Pseudo', n)]]), x, y))
  #     #remove parts of clusters that don't overlap with DNA foci
  #     result %>% mutate(!!paste0('Pseudo_in_DNA', n) :=
  #                         case_when(DNA_focus_idx == 0 ~ as.integer(0),
  #                                   T ~ !!as.name(paste0('Pseudo_in_DNA', n)))) -> result
  #   }
  #   
  #   #Classify by HP1
  #   segmentation = list('cluster' = segmCoresList[[pos]][['Foci']]$cluster,
  #                       'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
  #   class(segmentation) = c('dbscan_fast', 'dbscan')
  #   result$HP1_foci_dbscan = predict(object = segmentation,
  #                                    newdata = select(ungroup(result), x, y),
  #                                    data = select(ungroup(segmCoresList[[pos]][['Foci']]), x, y))
  #   for (n in 1:n_shuffle) {
  #     pseudosegmentation = list('cluster' = segmCoresList[[pos]][[paste0('Pseudo', n)]]$cluster,
  #                               'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
  #     class(pseudosegmentation) = c('dbscan_fast', 'dbscan')
  #     result[[paste0('HP1_pseudo_dbscan', n)]] =
  #       predict(object = pseudosegmentation,
  #               newdata = select(ungroup(result), x, y),
  #               data = select(ungroup(segmCoresList[[pos]][[paste0('Pseudo', n)]]), x, y))
  #   }
  #   
  #   #Overlap
  #   result$DNA_HP1_foci = result$HP1_foci_dbscan %in% DNA_HP1_foci[[pos]]
  #   for (n in 1:n_shuffle) {
  #     result[[paste0('DNA_HP1_pseudo', n)]] = result[[paste0('HP1_pseudo_dbscan', n)]] %in% DNA_HP1_foci[[pos]]
  #   }
  #   return(result)
  #   toc()
  # }) %>% bind_rows() -> jumpDist3
  # 
  # #Number of observations
  # {
  #   jumpDist3 %>%
  #     mutate(Compartment2 = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                     DNA_focus_idx!=0 & !DNA_HP1_foci ~ 'DNA-high',
  #                                     !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                     DNA_pseudo & !(HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 |
  #                                                      HP1_pseudo_dbscan3 | HP1_pseudo_dbscan4 |
  #                                                      HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'DNA pseudo',
  #                                     DNA_HP1_pseudo1 ~ 'DNA&HP1 pseudo 1',
  #                                     DNA_HP1_pseudo2 ~ 'DNA&HP1 pseudo 2',
  #                                     DNA_HP1_pseudo3 ~ 'DNA&HP1 pseudo 3',
  #                                     DNA_HP1_pseudo4 ~ 'DNA&HP1 pseudo 4',
  #                                     DNA_HP1_pseudo5 ~ 'DNA&HP1 pseudo 5',
  #                                     DNA_HP1_pseudo6 ~ 'DNA&HP1 pseudo 6',
  #                                     !DNA_HP1_pseudo1 & HP1_pseudo_dbscan1 ~ 'HP1 pseudo 1',#NB! These overlap -> pseudo6 ends up smaller
  #                                     !DNA_HP1_pseudo2 & HP1_pseudo_dbscan2 ~ 'HP1 pseudo 2',
  #                                     !DNA_HP1_pseudo3 & HP1_pseudo_dbscan3 ~ 'HP1 pseudo 3',
  #                                     !DNA_HP1_pseudo4 & HP1_pseudo_dbscan4 ~ 'HP1 pseudo 4',
  #                                     !DNA_HP1_pseudo5 & HP1_pseudo_dbscan5 ~ 'HP1 pseudo 5',
  #                                     !DNA_HP1_pseudo6 & HP1_pseudo_dbscan6 ~ 'HP1 pseudo 6',
  #                                     T ~ 'Euchromatin')) %>%
  #     group_by(Dish, Compartment2) %>% summarise(num_obs = n()) %>%
  #     mutate(Compartment = case_when(grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
  #                                    grepl('HP1 pseudo', Compartment2) ~ 'HP1 pseudo',
  #                                    T ~ Compartment2)) %>%
  #     mutate(Compartment2 = factor(Compartment2, levels = rev(c('DNA-high', 'HP1-high', 'DNA&HP1-high',
  #                                                               'DNA pseudo', paste0('HP1 pseudo ', 1:6),
  #                                                               paste0('DNA&HP1 pseudo ', 1:6),
  #                                                               'Euchromatin')))) %>%
  #     ggplot(aes(x = Compartment2, y = num_obs, fill = Compartment, group = Compartment2, label = num_obs)) +
  #     geom_bar(stat = 'identity', colour = 'gray20') +
  #     geom_text(y = 15000, hjust = 1) +
  #     facet_grid(Dish~.) +
  #     coord_flip(ylim = c(0, 15000)) +
  #     scale_fill_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
  #                                  'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
  #                                  'Euchromatin' = 'gray')) +
  #     labs(title = 'Number of observations', subtitle = 'HP1 segm with DBSCAN',
  #          fill = 'Compartment', x = '', y = 'Num obs') +
  #     theme_bw() +
  #     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
  #           aspect.ratio = 1, legend.position = 'none')
  #   ggsave('New_analysis/jumpDist_segm_HP1andDNA_numObs.png', width = 4, height = 6)
  # }
  # #Just the real foci
  # {
  #   jumpDist3 %>% mutate(Compartment = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                                DNA_focus_idx!=0 & !DNA_HP1_foci ~ 'DNA-high',
  #                                                !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                                T ~ 'Euchromatin')) %>%
  #     ggplot(aes(x = disp, colour = Compartment)) +
  #     geom_density() +
  #     scale_colour_manual(values = c('cyan', 'purple', 'gray', 'red')) +
  #     labs(title = 'Jump distance', subtitle = 'HP1 segm with DBSCAN', x = 'Displacement, nm', y = 'Density') +
  #     theme_bw() +
  #     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), aspect.ratio = 0.5)
  #   ggsave('New_analysis/jumpDist_segm_HP1andDNA.png', width = 4.5, height = 3)
  #   jumpDist3 %>% mutate(Compartment = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                                DNA_focus_idx!=0 & !DNA_HP1_foci ~ 'DNA-high',
  #                                                !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                                T ~ 'Euchromatin')) %>%
  #     ggplot(aes(x = disp, colour = Compartment)) +
  #     geom_density() +
  #     facet_grid(Dish~.) +
  #     scale_colour_manual(values = c('cyan', 'purple', 'gray', 'red')) +
  #     labs(title = 'Jump distance', subtitle = 'HP1 segm with DBSCAN', x = 'Displacement, nm', y = 'Density') +
  #     theme_bw() +
  #     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), aspect.ratio = 0.5)
  #   ggsave('New_analysis/jumpDist_segm_HP1andDNA_byDish.png', width = 4.5, height = 4.5)
  # }
  # #Just the real foci, DNA subsampled
  # {
  #   jumpDist3 %>% 
  #     mutate(Compartment2 = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                     Pseudo_in_DNA1!=0 & !DNA_HP1_foci ~ 'DNA-high (subsampled) 1',
  #                                     Pseudo_in_DNA2!=0 & !DNA_HP1_foci ~ 'DNA-high (subsampled) 2',
  #                                     Pseudo_in_DNA3!=0 & !DNA_HP1_foci ~ 'DNA-high (subsampled) 3',
  #                                     Pseudo_in_DNA4!=0 & !DNA_HP1_foci ~ 'DNA-high (subsampled) 4',
  #                                     Pseudo_in_DNA5!=0 & !DNA_HP1_foci ~ 'DNA-high (subsampled) 5',
  #                                     Pseudo_in_DNA6!=0 & !DNA_HP1_foci ~ 'DNA-high (subsampled) 6',
  #                                     !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                     T ~ 'Euchromatin')) %>%
  #     mutate(Compartment = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                    (Pseudo_in_DNA1!=0 | Pseudo_in_DNA2!=0 | Pseudo_in_DNA3!=0 |
  #                                       Pseudo_in_DNA4!=0 | Pseudo_in_DNA5!=0 | Pseudo_in_DNA6!=0)
  #                                    & !DNA_HP1_foci ~ 'DNA-high (subsampled)',
  #                                    !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                    T ~ 'Euchromatin')) %>%
  #     ggplot(aes(x = disp, colour = Compartment, group = Compartment2)) +
  #     geom_density() +
  #     scale_colour_manual(values = c('steelblue', 'purple', 'gray', 'red')) +
  #     labs(title = 'Jump distance', subtitle = 'HP1 segm with DBSCAN', x = 'Displacement, nm', y = 'Density',
  #          colour = 'Compartment') +
  #     theme_bw() +
  #     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), aspect.ratio = 0.5)
  #   ggsave('New_analysis/jumpDist_segm_HP1andDNA2.png', width = 5.1, height = 3)
  #   jumpDist3 %>%
  #     mutate(Compartment = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                    (Pseudo_in_DNA1!=0 | Pseudo_in_DNA2!=0 | Pseudo_in_DNA3!=0 |
  #                                       Pseudo_in_DNA4!=0 | Pseudo_in_DNA5!=0 | Pseudo_in_DNA6!=0)
  #                                    & !DNA_HP1_foci ~ 'DNA-high (subsampledx6)',
  #                                    !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                    T ~ 'Euchromatin')) %>%
  #     ggplot(aes(x = disp, colour = Compartment)) +
  #     geom_density() +
  #     scale_colour_manual(values = c('steelblue', 'purple', 'gray', 'red')) +
  #     labs(title = 'Jump distance', subtitle = 'HP1 segm with DBSCAN', x = 'Displacement, nm', y = 'Density',
  #          colour = 'Compartment') +
  #     theme_bw() +
  #     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), aspect.ratio = 0.5)
  #   ggsave('New_analysis/jumpDist_segm_HP1andDNA2_pool.png', width = 5.1, height = 3)
  # }
  # #Real and pseudofoci
  # {
  #   jumpDist3 %>%
  #     mutate(Compartment2 = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                     DNA_focus_idx!=0 & !DNA_HP1_foci ~ 'DNA-high',
  #                                     !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                     DNA_pseudo & !(HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 |
  #                                                      HP1_pseudo_dbscan3 | HP1_pseudo_dbscan4 |
  #                                                      HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'DNA pseudo',
  #                                     DNA_HP1_pseudo1 ~ 'DNA&HP1 pseudo 1',
  #                                     DNA_HP1_pseudo2 ~ 'DNA&HP1 pseudo 2',
  #                                     DNA_HP1_pseudo3 ~ 'DNA&HP1 pseudo 3',
  #                                     DNA_HP1_pseudo4 ~ 'DNA&HP1 pseudo 4',
  #                                     DNA_HP1_pseudo5 ~ 'DNA&HP1 pseudo 5',
  #                                     DNA_HP1_pseudo6 ~ 'DNA&HP1 pseudo 6',
  #                                     !DNA_HP1_pseudo1 & HP1_pseudo_dbscan1 ~ 'HP1 pseudo 1',#NB! These overlap -> pseudo6 ends up smaller
  #                                     !DNA_HP1_pseudo2 & HP1_pseudo_dbscan2 ~ 'HP1 pseudo 2',
  #                                     !DNA_HP1_pseudo3 & HP1_pseudo_dbscan3 ~ 'HP1 pseudo 3',
  #                                     !DNA_HP1_pseudo4 & HP1_pseudo_dbscan4 ~ 'HP1 pseudo 4',
  #                                     !DNA_HP1_pseudo5 & HP1_pseudo_dbscan5 ~ 'HP1 pseudo 5',
  #                                     !DNA_HP1_pseudo6 & HP1_pseudo_dbscan6 ~ 'HP1 pseudo 6',
  #                                     T ~ 'Euchromatin')) %>%
  #     mutate(Compartment = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                    DNA_focus_idx!=0 & !DNA_HP1_foci ~ 'DNA-high',
  #                                    !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                    DNA_pseudo & !(HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 |
  #                                                     HP1_pseudo_dbscan3 | HP1_pseudo_dbscan4 |
  #                                                     HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'DNA pseudo',
  #                                    (DNA_HP1_pseudo1 | DNA_HP1_pseudo2 | DNA_HP1_pseudo3 |
  #                                       DNA_HP1_pseudo4 | DNA_HP1_pseudo5 | DNA_HP1_pseudo6) ~ 'DNA&HP1 pseudo',
  #                                    (HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 | HP1_pseudo_dbscan3 |
  #                                       HP1_pseudo_dbscan4 | HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'HP1 pseudo',
  #                                    T ~ 'Euchromatin')) %>%
  #     mutate(Compartment = factor(Compartment, levels = c('DNA-high', 'HP1-high', 'DNA&HP1-high',
  #                                                         'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
  #                                                         'Euchromatin'))) %>%
  #     ggplot(aes(x = disp, colour = Compartment, group = Compartment2)) +
  #     geom_density() +
  #     facet_grid(Dish~.) +
  #     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
  #                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
  #                                    'Euchromatin' = 'gray')) +
  #     labs(title = 'Jump distance', subtitle = 'HP1 segm with DBSCAN',
  #          colour = 'Compartment', x = 'Displacement, nm', y = 'Density') +
  #     theme_bw() +
  #     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), aspect.ratio = 0.5)
  #   ggsave('New_analysis/jumpDist_segm_HP1andDNA_allCompartments.png', width = 4.7, height = 4.5)
  #   jumpDist3 %>%
  #     mutate(Compartment = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                    DNA_focus_idx!=0 & !DNA_HP1_foci ~ 'DNA-high',
  #                                    !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                    DNA_pseudo & !(HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 |
  #                                                     HP1_pseudo_dbscan3 | HP1_pseudo_dbscan4 |
  #                                                     HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'DNA pseudo',
  #                                    (DNA_HP1_pseudo1 | DNA_HP1_pseudo2 | DNA_HP1_pseudo3 |
  #                                       DNA_HP1_pseudo4 | DNA_HP1_pseudo5 | DNA_HP1_pseudo6) ~ 'DNA&HP1 pseudo(x6)',
  #                                    (HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 | HP1_pseudo_dbscan3 |
  #                                       HP1_pseudo_dbscan4 | HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'HP1 pseudo(x6)',
  #                                    T ~ 'Euchromatin')) %>%
  #     mutate(Compartment = factor(Compartment, levels = c('DNA-high', 'HP1-high', 'DNA&HP1-high',
  #                                                         'DNA pseudo', 'HP1 pseudo(x6)', 'DNA&HP1 pseudo(x6)',
  #                                                         'Euchromatin'))) %>%
  #     ggplot(aes(x = disp, colour = Compartment)) +
  #     geom_density() +
  #     facet_grid(Dish~.) +
  #     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
  #                                    'DNA pseudo' = 'blue', 'HP1 pseudo(x6)' = 'tomato4', 'DNA&HP1 pseudo(x6)' = 'purple4',
  #                                    'Euchromatin' = 'gray')) +
  #     labs(title = 'Jump distance', subtitle = 'HP1 segm with DBSCAN',
  #          colour = 'Compartment', x = 'Displacement, nm', y = 'Density') +
  #     theme_bw() +
  #     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), aspect.ratio = 0.5)
  #   ggsave('New_analysis/jumpDist_segm_HP1andDNA_allCompartments_pool_byDish.png', width = 4.9, height = 4.5)
  #   jumpDist3 %>%
  #     mutate(Compartment = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                    DNA_focus_idx!=0 & !DNA_HP1_foci ~ 'DNA-high',
  #                                    !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                    DNA_pseudo & !(HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 |
  #                                                     HP1_pseudo_dbscan3 | HP1_pseudo_dbscan4 |
  #                                                     HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'DNA pseudo',
  #                                    (DNA_HP1_pseudo1 | DNA_HP1_pseudo2 | DNA_HP1_pseudo3 |
  #                                       DNA_HP1_pseudo4 | DNA_HP1_pseudo5 | DNA_HP1_pseudo6) ~ 'DNA&HP1 pseudo(x6)',
  #                                    (HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 | HP1_pseudo_dbscan3 |
  #                                       HP1_pseudo_dbscan4 | HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'HP1 pseudo(x6)',
  #                                    T ~ 'Euchromatin')) %>%
  #     mutate(Compartment = factor(Compartment, levels = c('DNA-high', 'HP1-high', 'DNA&HP1-high',
  #                                                         'DNA pseudo', 'HP1 pseudo(x6)', 'DNA&HP1 pseudo(x6)',
  #                                                         'Euchromatin'))) %>%
  #     ggplot(aes(x = disp, colour = Compartment)) +
  #     geom_density() +
  #     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
  #                                    'DNA pseudo' = 'blue', 'HP1 pseudo(x6)' = 'tomato4', 'DNA&HP1 pseudo(x6)' = 'purple4',
  #                                    'Euchromatin' = 'gray')) +
  #     labs(title = 'Jump distance', subtitle = 'HP1 segm with DBSCAN',
  #          colour = 'Compartment', x = 'Displacement, nm', y = 'Density') +
  #     theme_bw() +
  #     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), aspect.ratio = 0.5)
  #   ggsave('New_analysis/jumpDist_segm_HP1andDNA_allCompartments_pool.png', width = 4.9, height = 3)
  # }
  # #Real and pseudofoci, DNA subsampled
  # {
  #   jumpDist3 %>%
  #     mutate(Compartment2 = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                     Pseudo_in_DNA1 & !DNA_HP1_foci ~ 'DNA-high (subsampled 1)',
  #                                     Pseudo_in_DNA2 & !DNA_HP1_foci ~ 'DNA-high (subsampled 2)',
  #                                     Pseudo_in_DNA3 & !DNA_HP1_foci ~ 'DNA-high (subsampled 3)',
  #                                     Pseudo_in_DNA4 & !DNA_HP1_foci ~ 'DNA-high (subsampled 4)',
  #                                     Pseudo_in_DNA5 & !DNA_HP1_foci ~ 'DNA-high (subsampled 5)',
  #                                     Pseudo_in_DNA6 & !DNA_HP1_foci ~ 'DNA-high (subsampled 6)',
  #                                     !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                     DNA_pseudo & !(HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 |
  #                                                      HP1_pseudo_dbscan3 | HP1_pseudo_dbscan4 |
  #                                                      HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'DNA pseudo',
  #                                     DNA_HP1_pseudo1 ~ 'DNA&HP1 pseudo 1',
  #                                     DNA_HP1_pseudo2 ~ 'DNA&HP1 pseudo 2',
  #                                     DNA_HP1_pseudo3 ~ 'DNA&HP1 pseudo 3',
  #                                     DNA_HP1_pseudo4 ~ 'DNA&HP1 pseudo 4',
  #                                     DNA_HP1_pseudo5 ~ 'DNA&HP1 pseudo 5',
  #                                     DNA_HP1_pseudo6 ~ 'DNA&HP1 pseudo 6',
  #                                     !DNA_HP1_pseudo1 & HP1_pseudo_dbscan1 ~ 'HP1 pseudo 1',#NB! These overlap -> pseudo6 ends up smaller
  #                                     !DNA_HP1_pseudo2 & HP1_pseudo_dbscan2 ~ 'HP1 pseudo 2',
  #                                     !DNA_HP1_pseudo3 & HP1_pseudo_dbscan3 ~ 'HP1 pseudo 3',
  #                                     !DNA_HP1_pseudo4 & HP1_pseudo_dbscan4 ~ 'HP1 pseudo 4',
  #                                     !DNA_HP1_pseudo5 & HP1_pseudo_dbscan5 ~ 'HP1 pseudo 5',
  #                                     !DNA_HP1_pseudo6 & HP1_pseudo_dbscan6 ~ 'HP1 pseudo 6',
  #                                     T ~ 'Euchromatin')) %>%
  #     mutate(Compartment = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                    (Pseudo_in_DNA1 | Pseudo_in_DNA2 | Pseudo_in_DNA3 |
  #                                       Pseudo_in_DNA4 | Pseudo_in_DNA5 | Pseudo_in_DNA6) &
  #                                      !DNA_HP1_foci ~ 'DNA-high (subsampled)',
  #                                    !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                    DNA_pseudo & !(HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 |
  #                                                     HP1_pseudo_dbscan3 | HP1_pseudo_dbscan4 |
  #                                                     HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'DNA pseudo',
  #                                    (DNA_HP1_pseudo1 | DNA_HP1_pseudo2 | DNA_HP1_pseudo3 |
  #                                       DNA_HP1_pseudo4 | DNA_HP1_pseudo5 | DNA_HP1_pseudo6) ~ 'DNA&HP1 pseudo',
  #                                    (HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 | HP1_pseudo_dbscan3 |
  #                                       HP1_pseudo_dbscan4 | HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'HP1 pseudo',
  #                                    T ~ 'Euchromatin')) %>%
  #     mutate(Compartment = factor(Compartment, levels = c('DNA-high (subsampled)', 'HP1-high', 'DNA&HP1-high',
  #                                                         'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
  #                                                         'Euchromatin'))) %>%
  #     ggplot(aes(x = disp, colour = Compartment, group = Compartment2)) +
  #     geom_density() +
  #     facet_grid(Dish~.) +
  #     scale_colour_manual(values = c('DNA-high (subsampled)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
  #                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
  #                                    'Euchromatin' = 'gray')) +
  #     labs(title = 'Jump distance', subtitle = 'HP1 segm with DBSCAN',
  #          colour = 'Compartment', x = 'Displacement, nm', y = 'Density') +
  #     theme_bw() +
  #     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), aspect.ratio = 0.5)
  #   ggsave('New_analysis/jumpDist_segm_HP1andDNA_allCompartments2.png', width = 5, height = 4.5)
  #   jumpDist3 %>%
  #     mutate(Compartment = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                    (Pseudo_in_DNA1 | Pseudo_in_DNA2 | Pseudo_in_DNA3 |
  #                                       Pseudo_in_DNA4 | Pseudo_in_DNA5 | Pseudo_in_DNA6) &
  #                                      !DNA_HP1_foci ~ 'DNA-high (subsampledx6)',
  #                                    !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                    DNA_pseudo & !(HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 |
  #                                                     HP1_pseudo_dbscan3 | HP1_pseudo_dbscan4 |
  #                                                     HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'DNA pseudo',
  #                                    (DNA_HP1_pseudo1 | DNA_HP1_pseudo2 | DNA_HP1_pseudo3 |
  #                                       DNA_HP1_pseudo4 | DNA_HP1_pseudo5 | DNA_HP1_pseudo6) ~ 'DNA&HP1 pseudo (x6)',
  #                                    (HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 | HP1_pseudo_dbscan3 |
  #                                       HP1_pseudo_dbscan4 | HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'HP1 pseudo (x6)',
  #                                    T ~ 'Euchromatin')) %>%
  #     mutate(Compartment = factor(Compartment, levels = c('DNA-high (subsampledx6)', 'HP1-high', 'DNA&HP1-high',
  #                                                         'DNA pseudo', 'HP1 pseudo (x6)', 'DNA&HP1 pseudo (x6)',
  #                                                         'Euchromatin'))) %>%
  #     ggplot(aes(x = disp, colour = Compartment)) +
  #     geom_density() +
  #     facet_grid(Dish~.) +
  #     scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
  #                                    'DNA pseudo' = 'blue', 'HP1 pseudo (x6)' = 'tomato4', 'DNA&HP1 pseudo (x6)' = 'purple4',
  #                                    'Euchromatin' = 'gray')) +
  #     labs(title = 'Jump distance', subtitle = 'HP1 segm with DBSCAN',
  #          colour = 'Compartment', x = 'Displacement, nm', y = 'Density') +
  #     theme_bw() +
  #     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), aspect.ratio = 0.5)
  #   ggsave('New_analysis/jumpDist_segm_HP1andDNA_allCompartments2_pool_byDish.png', width = 5.15, height = 4.5)
  #   jumpDist3 %>%
  #     mutate(Compartment = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                    (Pseudo_in_DNA1 | Pseudo_in_DNA2 | Pseudo_in_DNA3 |
  #                                       Pseudo_in_DNA4 | Pseudo_in_DNA5 | Pseudo_in_DNA6) &
  #                                      !DNA_HP1_foci ~ 'DNA-high (subsampledx6)',
  #                                    !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                    DNA_pseudo & !(HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 |
  #                                                     HP1_pseudo_dbscan3 | HP1_pseudo_dbscan4 |
  #                                                     HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'DNA pseudo',
  #                                    (DNA_HP1_pseudo1 | DNA_HP1_pseudo2 | DNA_HP1_pseudo3 |
  #                                       DNA_HP1_pseudo4 | DNA_HP1_pseudo5 | DNA_HP1_pseudo6) ~ 'DNA&HP1 pseudo (x6)',
  #                                    (HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 | HP1_pseudo_dbscan3 |
  #                                       HP1_pseudo_dbscan4 | HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'HP1 pseudo (x6)',
  #                                    T ~ 'Euchromatin')) %>%
  #     mutate(Compartment = factor(Compartment, levels = c('DNA-high (subsampledx6)', 'HP1-high', 'DNA&HP1-high',
  #                                                         'DNA pseudo', 'HP1 pseudo (x6)', 'DNA&HP1 pseudo (x6)',
  #                                                         'Euchromatin'))) %>%
  #     ggplot(aes(x = disp, colour = Compartment)) +
  #     geom_density() +
  #     scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
  #                                    'DNA pseudo' = 'blue', 'HP1 pseudo (x6)' = 'tomato4', 'DNA&HP1 pseudo (x6)' = 'purple4',
  #                                    'Euchromatin' = 'gray')) +
  #     labs(title = 'Jump distance', subtitle = 'HP1 segm with DBSCAN',
  #          colour = 'Compartment', x = 'Displacement, nm', y = 'Density') +
  #     theme_bw() +
  #     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), aspect.ratio = 0.5)
  #   ggsave('New_analysis/jumpDist_segm_HP1andDNA_allCompartments2_pool.png', width = 5.15, height = 3)
  # }
  # #Fitting
  # {
  #   fitJD = function(JD, range = seq(0, 0.4, by = 0.01), dim = 2, dt) {
  #     df = data.frame(density = ecdf(JD)(range), r = range)
  #     fit1 = tryCatch(nls(density ~ (1-exp(-r^2/(2*dim*D1*dt))), data = df,
  #                         start = list(D1 = 1),
  #                         lower = 0, upper = 100, algorithm = 'port'), error = function(e){return(NA)})
  #     fit2 = tryCatch(nls(density ~ (1-A1*exp(-r^2/(2*dim*D1*dt))-(1-A1)*exp(-r^2/(2*dim*D2*dt))), data = df,
  #                         start = list(A1 = 0.5, D1 = 0.1, D2 = 1),
  #                         lower = c(0, 0, 0), upper = c(1, 100, 100), algorithm = 'port'), error = function(e){return(NA)})
  #     fit3 = tryCatch(nls(density ~ (1-A1*exp(-r^2/(2*dim*D1*dt))-A2*exp(-r^2/(2*dim*D2*dt))-(1-A1-A2)*exp(-r^2/(2*dim*D3*dt))),
  #                         data = df,
  #                         start = list(A1 = 0.3, A2 = 0.3, D1 = 0.1, D2 = 1, D3 = 20),
  #                         lower = c(0, 0, 0, 0, 0), upper = c(1, 1, 100, 100, 100), algorithm = 'port'),
  #                     error = function(e){return(NA)})
  #     return(list(fit1, fit2, fit3))
  #   }
  #   fitFun_cum = function(r, n = 1, coef, dim = 2, dt) {
  #     if (n == 1) {
  #       return(1-exp(-r^2/(2*dim*coef[['D1']]*dt)))
  #     } else if (n == 2) {
  #       return(1 - coef[['A1']]*exp(-r^2/(2*dim*coef[['D1']]*dt)) -
  #                (1-coef[['A1']])*exp(-r^2/(2*dim*coef[['D2']]*dt)))
  #     } else if (n == 3) {
  #       return(1 - coef[['A1']]*exp(-r^2/(2*dim*coef[['D1']]*dt)) -
  #                coef[['A2']]*exp(-r^2/(2*dim*coef[['D2']]*dt)) -
  #                (1 - coef[['A1']] - coef[['A2']])*exp(-r^2/(2*dim*coef[['D3']]*dt)))
  #     }
  #   }
  #   fitFun_hist = function(r, binsize, n = 1, coef, dim = 2, dt) {
  #     if (n == 1) {
  #       return(exp(-r^2/(2*dim*coef[['D1']]*dt)) - exp(-(r+binsize)^2/(2*dim*coef[['D1']]*dt)))
  #     } else if (n == 2) {
  #       return(coef[['A1']]*(exp(-r^2/(2*dim*coef[['D1']]*dt)) - exp(-(r+binsize)^2/(2*dim*coef[['D1']]*dt))) +
  #                (1-coef[['A1']])*(exp(-r^2/(2*dim*coef[['D2']]*dt)) - exp(-(r+binsize)^2/(2*dim*coef[['D2']]*dt))))
  #     } else if (n == 3) {
  #       return(coef[['A1']]*(exp(-r^2/(2*dim*coef[['D1']]*dt)) - exp(-(r+binsize)^2/(2*dim*coef[['D1']]*dt))) +
  #                coef[['A2']]*(exp(-r^2/(2*dim*coef[['D2']]*dt)) - exp(-(r+binsize)^2/(2*dim*coef[['D2']]*dt))) +
  #                (1-coef[['A1']]-coef[['A2']])*(exp(-r^2/(2*dim*coef[['D3']]*dt)) - exp(-(r+binsize)^2/(2*dim*coef[['D3']]*dt))))
  #     }
  #   }
  #   
  #   fitJDsetPrec = function(JD, range = seq(0, 0.4, by = 0.01), dim = 2, dt, prec.sq = 0.02^2) {
  #     df = data.frame(density = ecdf(JD)(range), r = range)
  #     fit1 = tryCatch(nls(density ~ (1-exp(-r^2/(2*dim*(D1*dt+prec.sq)))), data = df,
  #                           start = list(D1 = 1),
  #                           lower = 0, upper = 100, algorithm = "port"), error = function(e){return(NA)})
  #     fit2 = tryCatch(nls(density ~ (1-A1*exp(-r^2/(2*dim*(D1*dt+prec.sq)))-
  #                                        (1-A1)*exp(-r^2/(2*dim*(D2*dt+prec.sq)))),
  #                           data = df,
  #                           start = list(A1 = 0.5, D1 = 0.1, D2 = 1),
  #                           lower = c(0, 0, 0), upper = c(1, 100, 100), algorithm = "port"), error = function(e){return(NA)})
  #     fit3 = tryCatch(nls(density ~ (1-A1*exp(-r^2/(2*dim*(D1*dt+prec.sq)))-
  #                                        A2*exp(-r^2/(2*dim*(D2*dt+prec.sq)))-
  #                                        (1-A1-A2)*exp(-r^2/(2*dim*(D3*dt+prec.sq)))),
  #                           data = df,
  #                           start = list(A1 = 0.3, A2 = 0.3, D1 = 0.1, D2 = 1, D3 = 20),
  #                           lower = c(0, 0, 0, 0, 0), upper = c(1, 1, 100, 100, 100), algorithm = "port"),
  #                     error = function(e){return(NA)})
  #     return(list(fit1, fit2, fit3))
  #   }
  #   #Fitting loc error fails =(
  #   {
  #   # fitJDsetPrecBlur = function(JD, range = seq(0, 0.4, by = 0.01), dim = 2, dt, prec.sq = 0.05^2) {
  #   #   df = data.frame(density = ecdf(JD)(range), r = range)
  #   #   fit1 = tryCatch(nls(density ~ (1-exp(-r^2/(2*dim*((2/3)*D1*dt+prec.sq)))), data = df,
  #   #                       start = list(D1 = 1),
  #   #                       lower = 0, upper = 100, algorithm = "port"), error = function(e){return(NA)})
  #   #   fit2 = tryCatch(nls(density ~ (1-A1*exp(-r^2/(2*dim*((2/3)*D1*dt+prec.sq)))-
  #   #                                    (1-A1)*exp(-r^2/(2*dim*((2/3)*D2*dt+prec.sq)))),
  #   #                       data = df,
  #   #                       start = list(A1 = 0.5, D1 = 0.1, D2 = 1),
  #   #                       lower = c(0, 0, 0), upper = c(1, 100, 100), algorithm = "port"), error = function(e){return(NA)})
  #   #   fit3 = tryCatch(nls(density ~ (1-A1*exp(-r^2/(2*dim*((2/3)*D1*dt+prec.sq)))-
  #   #                                    A2*exp(-r^2/(2*dim*((2/3)*D2*dt+prec.sq)))-
  #   #                                    (1-A1-A2)*exp(-r^2/(2*dim*((2/3)*D3*dt+prec.sq)))),
  #   #                       data = df,
  #   #                       start = list(A1 = 0.3, A2 = 0.3, D1 = 0.1, D2 = 1, D3 = 20),
  #   #                       lower = c(0, 0, 0, 0, 0), upper = c(1, 1, 100, 100, 100), algorithm = "port"),
  #   #                   error = function(e){return(NA)})
  #   #   return(list(fit1, fit2, fit3))
  #   # }
  #   # fitJDfitPrec = function(JD, range = seq(0, 0.4, by = 0.01), dim = 2, dt) {
  #   #   df = data.frame(density = ecdf(JD)(range), r = range)
  #   #   fit1 = tryCatch(nls(density ~ (1-exp(-r^2/(2*dim*(D1*dt+prec.sq)))), data = df,
  #   #                         start = list(D1 = 1, prec.sq = 0.02^2),
  #   #                         lower = c(0, 0), upper = c(100, 0.1), algorithm = 'port'), error = function(e){return(NA)})
  #   #   fit2 = tryCatch(nls(density ~ (1-A1*exp(-r^2/(2*dim*(D1*dt+prec.sq)))-(1-A1)*exp(-r^2/(2*dim*(D2*dt+prec.sq)))),
  #   #                         data = df,
  #   #                         start = list(A1 = 0.5, D1 = 0.1, D2 = 1, prec.sq = 0.02^2),
  #   #                         lower = c(0, 0, 0, 0), upper = c(1, 100, 100, 0.1), algorithm = 'port'), error = function(e){return(NA)})
  #   #   fit3 = tryCatch(nls(density ~ (1-A1*exp(-r^2/(2*dim*(D1*dt+prec.sq)))-A2*exp(-r^2/(2*dim*(D2*dt+prec.sq)))-
  #   #                                      (1-A1-A2)*exp(-r^2/(2*dim*(D3*dt+prec.sq)))),
  #   #                         data = df,
  #   #                         start = list(A1 = 0.3, A2 = 0.3, D1 = 0.1, D2 = 1, D3 = 20, prec.sq = 0.02^2),
  #   #                         lower = c(0, 0, 0, 0, 0, 0), upper = c(1, 1, 100, 100, 100, 0.1), algorithm = 'port'),
  #   #                   error = function(e){return(NA)})
  #   #   return(list(fit1, fit2, fit3))
  #   # }
  #   # fitJDfitPrecBlur = function(JD, range = seq(0, 0.4, by = 0.01), dim = 2, dt) {
  #   #   df = data.frame(density = ecdf(JD)(range), r = range)
  #   #   fit1 = tryCatch(nls(density ~ (1-exp(-r^2/(2*dim*((2/3)*D1*dt+prec.sq)))), data = df,
  #   #                         start = list(D1 = 1, prec.sq = 0.05^2),
  #   #                         lower = c(0, 0), upper = c(100, 0.1), algorithm = 'port'), algorithm = 'port', error = function(e){return(NA)})
  #   #   fit2 = tryCatch(nls(density ~ (1-A1*exp(-r^2/(2*dim*((2/3)*D1*dt+prec.sq)))-
  #   #                                      (1-A1)*exp(-r^2/(2*dim*((2/3)*D2*dt+prec.sq)))),
  #   #                         data = df,
  #   #                         start = list(A1 = 0.5, D1 = 0.1, D2 = 1, prec.sq = 0.05^2),
  #   #                         lower = c(0, 0, 0, 0), upper = c(1, 100, 100, 0.1), algorithm = 'port'), algorithm = 'port', error = function(e){return(NA)})
  #   #   fit3 = tryCatch(nls(density ~ (1-A1*exp(-r^2/(2*dim*((2/3)*D1*dt+prec.sq)))-
  #   #                                      A2*exp(-r^2/(2*dim*((2/3)*D2*dt+prec.sq)))-
  #   #                                      (1-A1-A2)*exp(-r^2/(2*dim*((2/3)*D3*dt+prec.sq)))),
  #   #                         data = df,
  #   #                         start = list(A1 = 0.3, A2 = 0.3, D1 = 0.1, D2 = 1, D3 = 20, prec.sq = 0.05^2),
  #   #                         lower = c(0, 0, 0, 0, 0, 0), upper = c(1, 1, 100, 100, 100, 0.1), algorithm = 'port'),
  #   #                   error = function(e){return(NA)})
  #   #   return(list(fit1, fit2, fit3))
  #   # }
  #   }
  #   fitFunPrec_cum = function(r, n = 1, coef, dim = 2, dt, prec.sq = 0.02^2) {
  #     if (prec.sq == 'fit') {
  #       prec.sq = coef[['prec.sq']]
  #     } else {
  #       prec.sq = prec.sq
  #     }
  #     if (n == 1) {
  #       return(1-exp(-r^2/(2*dim*(coef[['D1']]*dt+prec.sq))))
  #     } else if (n == 2) {
  #       return(1 - coef[['A1']]*exp(-r^2/(2*dim*(coef[['D1']]*dt+prec.sq))) -
  #                (1-coef[['A1']])*exp(-r^2/(2*dim*(coef[['D2']]*dt+prec.sq))))
  #     } else if (n == 3) {
  #       return(1 - coef[['A1']]*exp(-r^2/(2*dim*(coef[['D1']]*dt+prec.sq))) -
  #                coef[['A2']]*exp(-r^2/(2*dim*(coef[['D2']]*dt+prec.sq))) -
  #                (1 - coef[['A1']] - coef[['A2']])*exp(-r^2/(2*dim*(coef[['D3']]*dt+prec.sq))))
  #     }
  #   }
  #   fitFunPrec_hist = function(r, binsize, n = 1, coef, dim = 2, dt, prec.sq = 0.02^2) {
  #     if (prec.sq == 'fit') {
  #       prec.sq = coef[['prec.sq']]
  #     } else {
  #       prec.sq = prec.sq
  #     }
  #     if (n == 1) {
  #       return(exp(-r^2/(2*dim*(coef[['D1']]*dt+prec.sq))) - exp(-(r+binsize)^2/(2*dim*(coef[['D1']]*dt+prec.sq))))
  #     } else if (n == 2) {
  #       return(coef[['A1']]*(exp(-r^2/(2*dim*(coef[['D1']]*dt+prec.sq))) - exp(-(r+binsize)^2/(2*dim*(coef[['D1']]*dt+prec.sq)))) +
  #                (1-coef[['A1']])*(exp(-r^2/(2*dim*(coef[['D2']]*dt+prec.sq))) - exp(-(r+binsize)^2/(2*dim*(coef[['D2']]*dt+prec.sq)))))
  #     } else if (n == 3) {
  #       return(coef[['A1']]*(exp(-r^2/(2*dim*(coef[['D1']]*dt+prec.sq))) - exp(-(r+binsize)^2/(2*dim*(coef[['D1']]*dt+prec.sq)))) +
  #                coef[['A2']]*(exp(-r^2/(2*dim*(coef[['D2']]*dt+prec.sq))) - exp(-(r+binsize)^2/(2*dim*(coef[['D2']]*dt+prec.sq)))) +
  #                (1-coef[['A1']]-coef[['A2']])*(exp(-r^2/(2*dim*(coef[['D3']]*dt+prec.sq))) -
  #                                                 exp(-(r+binsize)^2/(2*dim*(coef[['D3']]*dt+prec.sq)))))
  #     }
  #   }
  #   
  #   dt = 0.003
  #   sigma = 0.021
  #   
  #   ##classify
  #   {
  #   jumpDist3 %>%
  #     mutate(Compartment2 = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                     DNA_focus_idx!=0 & !DNA_HP1_foci ~ 'DNA-high',
  #                                     !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                     DNA_pseudo & !(HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 |
  #                                                      HP1_pseudo_dbscan3 | HP1_pseudo_dbscan4 |
  #                                                      HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'DNA pseudo',
  #                                     DNA_HP1_pseudo1 ~ 'DNA&HP1 pseudo 1',
  #                                     DNA_HP1_pseudo2 ~ 'DNA&HP1 pseudo 2',
  #                                     DNA_HP1_pseudo3 ~ 'DNA&HP1 pseudo 3',
  #                                     DNA_HP1_pseudo4 ~ 'DNA&HP1 pseudo 4',
  #                                     DNA_HP1_pseudo5 ~ 'DNA&HP1 pseudo 5',
  #                                     DNA_HP1_pseudo6 ~ 'DNA&HP1 pseudo 6',
  #                                     !DNA_HP1_pseudo1 & HP1_pseudo_dbscan1 ~ 'HP1 pseudo 1',#NB! These overlap -> pseudo6 ends up smaller
  #                                     !DNA_HP1_pseudo2 & HP1_pseudo_dbscan2 ~ 'HP1 pseudo 2',
  #                                     !DNA_HP1_pseudo3 & HP1_pseudo_dbscan3 ~ 'HP1 pseudo 3',
  #                                     !DNA_HP1_pseudo4 & HP1_pseudo_dbscan4 ~ 'HP1 pseudo 4',
  #                                     !DNA_HP1_pseudo5 & HP1_pseudo_dbscan5 ~ 'HP1 pseudo 5',
  #                                     !DNA_HP1_pseudo6 & HP1_pseudo_dbscan6 ~ 'HP1 pseudo 6',
  #                                     T ~ 'Euchromatin')) %>%
  #     mutate(Compartment = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
  #                                    DNA_focus_idx!=0 & !DNA_HP1_foci ~ 'DNA-high',
  #                                    !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
  #                                    DNA_pseudo & !(HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 |
  #                                                     HP1_pseudo_dbscan3 | HP1_pseudo_dbscan4 |
  #                                                     HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'DNA pseudo',
  #                                    (DNA_HP1_pseudo1 | DNA_HP1_pseudo2 | DNA_HP1_pseudo3 |
  #                                       DNA_HP1_pseudo4 | DNA_HP1_pseudo5 | DNA_HP1_pseudo6) ~ 'DNA&HP1 pseudo',
  #                                    (HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 | HP1_pseudo_dbscan3 |
  #                                       HP1_pseudo_dbscan4 | HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'HP1 pseudo',
  #                                    T ~ 'Euchromatin')) %>%
  #     mutate(Compartment = factor(Compartment, levels = c('DNA-high', 'HP1-high', 'DNA&HP1-high',
  #                                                         'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
  #                                                         'Euchromatin'))) -> jumpDist3
  #   jumpDist3 %>% split(.$Compartment) -> jumpDist3_list
  #   }
  #   ##fit: with loc precision better, all comp can be fit to 3 populations;
  #   #however, too high loc precision leads to poorer fits (slow population peak shifted to the right)
  #   #empirically, 21 nm looks good: better fit than 23 nm, almost as good a fit as, and lower D1 values than 20 nm
  #   {
  #   lapply(jumpDist3_list, function(x) {fitJD(x$disp/1000, dt = dt)}) -> fits_comps
  #   lapply(fits_comps, function(x) bind_rows(sapply(x, function(y) {
  #     if (!is.na(y)) {data.frame(as.list(coef(y)), RSS = sum(residuals(y)^2), BIC = glance(y)$BIC)}
  #   }), .id = 'num_exp')) %>%
  #     bind_rows(.id = 'Compartment') %>%
  #     mutate(A1 = replace_na(A1, 1),  A3 = 1-A1-A2, A2 = case_when(num_exp == 2 ~ 1-A1, T ~ A2)) %>%
  #     select(Compartment, num_exp, D1, A1, D2, A2, D3, A3, RSS, BIC) %>%
  #     pivot_longer(D1:A3, names_to = c(".value", "set"), names_pattern = "(\\w)(\\d)") %>%
  #     replace_na(list('D' = 10000, 'A' = 0)) %>%
  #     group_by(Compartment, num_exp) %>% summarise(D1 = sort(D)[1],
  #                                                  D2 = sort(D)[2],
  #                                                  D3 = sort(D)[3],
  #                                                  A1 = A[which(D==sort(D)[1])][1],
  #                                                  A2 = A[which(D==sort(D)[2])][1],
  #                                                  A3 = A[which(D==sort(D)[3])][1],
  #                                                  RSS = RSS[1],
  #                                                  BIC = BIC[1]) %>%
  #     mutate_at(vars(starts_with('A')), ~replace(., .==0, NA)) %>%
  #     mutate_at(vars(starts_with('D')), ~replace(., .==10000, NA)) -> fitStats_comps
  #   
  #   
  #   lapply(jumpDist3_list, function(x) {fitJDsetPrec(x$disp/1000, dt = dt, prec.sq = sigma^2)}) -> fitsPrec_comps
  #   lapply(fitsPrec_comps, function(x) bind_rows(sapply(x, function(y) {
  #     if (!is.na(y)) {data.frame(as.list(coef(y)), RSS = sum(residuals(y)^2), BIC = glance(y)$BIC)}
  #   }), .id = 'num_exp')) %>%
  #     bind_rows(.id = 'Compartment') %>% 
  #     mutate(A1 = replace_na(A1, 1),  A3 = 1-A1-A2, A2 = case_when(num_exp == 2 ~ 1-A1, T ~ A2)) %>%
  #     select(Compartment, num_exp, D1, A1, D2, A2, D3, A3, RSS, BIC) %>%
  #     pivot_longer(D1:A3, names_to = c(".value", "set"), names_pattern = "(\\w)(\\d)") %>%
  #     replace_na(list('D' = 10000, 'A' = 0)) %>%
  #     group_by(Compartment, num_exp) %>% summarise(D1 = sort(D)[1],
  #                                                  D2 = sort(D)[2],
  #                                                  D3 = sort(D)[3],
  #                                                  A1 = A[which(D==sort(D)[1])][1],
  #                                                  A2 = A[which(D==sort(D)[2])][1],
  #                                                  A3 = A[which(D==sort(D)[3])][1],
  #                                                  RSS = RSS[1],
  #                                                  BIC = BIC[1]) %>%
  #     mutate_at(vars(starts_with('A')), ~replace(., .==0, NA)) %>%
  #     mutate_at(vars(starts_with('D')), ~replace(., .==10000, NA)) -> fitPrecStats_comps
  #   }
  #   ##predict and plot
  #   {
  #   lapply(fitsPrec_comps, function(x) {lapply(seq_along(x), function(i) {
  #     points = seq(0, 0.4, by = 0.002)
  #     if (!is.na(x[[i]])) {
  #       prediction = fitFunPrec_hist(points, binsize = 0.002, n = i, coef = coef(x[[i]]),
  #                                    dt = dt, prec.sq = sigma^2)
  #       return(data.frame(points = points, prediction = prediction))
  #     }
  #   })}) -> predictions_hist
  #   lapply(predictions_hist, function(x) bind_rows(x, .id = 'num_exp')) %>% bind_rows(.id = 'Compartment') -> predictions_hist
  #   
  #   ggplot(jumpDist3, aes(x = disp/1000, y = ..density..*0.002, fill = Compartment)) +
  #     geom_histogram(binwidth = 0.002) +
  #     facet_wrap(~Compartment) +
  #     scale_fill_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
  #                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
  #                                    'Euchromatin' = 'gray')) +
  #     geom_line(data = filter(predictions_hist, num_exp != 1),
  #               aes(x = points, y = prediction, colour = num_exp)) +
  #     scale_colour_manual(values = c('black', 'orange')) +
  #     labs(title = 'Fitting of JD', subtitle = 'Sigma = 21 nm',
  #          colour = 'Num exp', fill = 'Compartment', x = 'Displacement, um', y = 'Density') +
  #     theme_bw() +
  #     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), aspect.ratio = 0.5)
  #   ggsave('New_analysis/jumpDist_fit.png', width = 8, height = 5)
  #   }
  #   ##plot values
  #   {
  #     filter(fitPrecStats_comps, num_exp != 1) %>%
  #       pivot_longer(D1:D3, values_to = 'D', names_to = 'Population') %>%
  #       mutate(Compartment = factor(Compartment, levels = c('DNA-high', 'DNA pseudo',
  #                                                           'HP1-high', 'HP1 pseudo',
  #                                                           'DNA&HP1-high', 'DNA&HP1 pseudo',
  #                                                           'Euchromatin'))) %>%
  #       ggplot(aes(x = Compartment, y = log10(D), colour = Compartment)) +
  #       geom_point() +
  #       facet_wrap(num_exp~.) +
  #       scale_y_continuous(labels = math_format(10^.x),
  #                          sec.axis = sec_axis(trans=~10^., breaks = c(0.003, 0.01, 0.03, 0.1, 0.3, 1, 3))) +
  #       scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
  #                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
  #                                      'Euchromatin' = 'gray')) +
  #       scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
  #       labs(x = '', y = 'D, um2/s') +
  #       theme_bw() +
  #       theme(aspect.ratio = 0.8, legend.position = 'none', axis.text.x = element_text(angle = 30, hjust = 1)) -> Ds
  #     filter(fitPrecStats_comps, num_exp != 1) %>%
  #       pivot_longer(A1:A3, values_to = 'A', names_to = 'Population') %>%
  #       mutate(Compartment = factor(Compartment, levels = c('DNA-high', 'DNA pseudo',
  #                                                           'HP1-high', 'HP1 pseudo',
  #                                                           'DNA&HP1-high', 'DNA&HP1 pseudo',
  #                                                           'Euchromatin'))) %>%
  #       ggplot(aes(x = Compartment, y = A, fill = Compartment, alpha = factor(Population, levels = c('A3', 'A2', 'A1')))) +
  #       geom_col(colour = 'gray20', width = 0.8) +
  #       facet_wrap(num_exp~.) +
  #       scale_fill_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
  #                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
  #                                      'Euchromatin' = 'gray'), guide = 'none') +
  #       scale_alpha_manual(values = c(0, 0.2, 1), labels = c('Unbound2', 'Unbound1', 'Bound')) +
  #       labs(x = '', y = 'Proportion of population', alpha = '') +
  #       theme_bw() +
  #       theme(aspect.ratio = 0.8, axis.text.x = element_text(angle = 30, hjust = 1)) -> As
  #     ggarrange(Ds, As, nrow = 2, ncol = 1)
  #     ggsave('New_analysis/jumpDist_DA.png', width = 6, height = 6)
  #   }
  #   ##multiply by density
  #   {
  #     filter(fitPrecStats_comps, num_exp == 3, Compartment != 'Euchromatin') %>%
  #       select(Compartment, A1) %>% rename(Bound = A1) %>% ungroup() %>%
  #       mutate(compartment = c('DNA foci', 'HP1&DNA foci', 'HP1&DNA foci', 'DNA foci', 'HP1 foci', 'HP1 foci'),
  #              type = c('real', 'real', 'pseudo', 'pseudo', 'real', 'pseudo')) %>%
  #       left_join(data.frame(compartment = rep(c('DNA foci', 'HP1 foci', 'HP1&DNA foci'), each = 2),
  #                            type = rep(c('real', 'pseudo'), 3),
  #                            HP1_dens = c(1840.13, 1578.81,
  #                                         2962.76, 1366.74,
  #                                         3336.00, 1458.66))) %>%
  #       mutate(Unbound = 1-Bound) %>% group_by(compartment) %>% arrange(type) %>%
  #       summarise(Unbound_frac_ratio = Unbound[2]/Unbound[1], HP1_dens_ratio = HP1_dens[2]/HP1_dens[1]) %>%
  #       mutate(Unbound_conc_ratio = Unbound_frac_ratio*HP1_dens_ratio) -> unbound_conc
  #     ggplot(unbound_conc, aes(x = compartment, y = Unbound_conc_ratio, fill = compartment)) +
  #       geom_col(position = 'dodge', colour = 'gray20', width = 0.8) +
  #       geom_hline(yintercept = 1, colour = 'black') +
  #       scale_fill_manual(values = c('cyan', 'purple', 'red'),
  #                         labels = c('DNA-high', 'HP1&DNA-high', 'HP1-high'),
  #                         guide = 'none') +
  #       scale_alpha_manual(values = c(1, 0.5)) +
  #       labs(x = 'Compartment', y = '[Unbound HP1]_foci/\n[Unbound HP1]_pseudo', fill = 'Compartment',
  #            title = 'Enrichment in unbound HP1b\nwithin different foci',
  #            subtitle = 'Compared to pseudo controls') +
  #       theme_bw() +
  #       theme(aspect.ratio = 0.8, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  #     ggsave('New_analysis/HP1free_conc_ratios_JDfit.png', width = 3, height = 3)
  #   }
  # }
}

##MSD - not included in the paper
{
# #Calculate jump distances in intervals; sanity checks
# #Calculate, classify each jump based on its midpoint (~25 min)
# jumpDist_1234 = readRDS('jumpDist_1234.rds')
# 
# lapply(seq_along(data_list), function(c) {
#   tic()
#   pos = names(data_list)[c]
#   print(pos)
#   df = select(filter(data_list[[pos]], Frame > thresholds2[c]),
#               Dish, Pos, Track, Frame, x, y, Track_length)
#   dim_x = foci_list[[pos]]$dim_x[1]
#   dim_y = foci_list[[pos]]$dim_y[1]
#   
#   for (i in 1:8) { #calculation takes ~5 min
#     print(i)
#     #Calculate jump dist
#     df %>% filter(Track_length > i) %>% group_by(Dish, Pos, Track) %>%
#       filter(n() > i) %>% #because initial filtering by frame cuts some trajectories up
#       mutate(diff.x = c(rep(NA, i), diff(x, i)), diff.y = c(rep(NA, i), diff(y, i)),
#              x_new = (x+lag(x, i))/2, y_new = (y+lag(y, i))/2,
#              Frame = floor((Frame+lag(Frame, i))/2)) %>%
#       mutate(disp = sqrt(diff.x^2 + diff.y^2), x = x_new, y = y_new) %>%
#       filter(!is.na(disp)) %>% select(-x_new, -y_new) %>% mutate(interval = i) -> tmp
#     #DNA classification
#     tmp %>% group_by(Dish, Pos, Track) %>%
#       mutate(px_id = floor(x*subsampling/pixel)*dim_y + ceiling(y*subsampling/pixel),
#              frame = ceiling(Frame/WFaveraging)) %>%
#       left_join(., foci_list[[pos]], by = c('px_id', 'frame')) %>%
#       select(-x_left, -y_top, -dim_y, -dim_x, -pos, -dish) %>%
#       left_join(., pfoci_list[[pos]], by = c('px_id', 'frame')) %>%
#       select(-x_left, -y_top, -dim_y, -dim_x, -pos, -dish, -frame, -px_id) %>%
#       rename(DNA_focus_idx = focus_id, DNA_pseudo = pseudofocus) %>%
#       mutate(DNA_focus_idx = replace_na(DNA_focus_idx, 0),
#              DNA_pseudo = replace_na(DNA_pseudo, 0)) -> tmp
#     
#     #DNA foci, subsampled in HP1-foci-like manner
#     for (n in 1:n_shuffle) {
#       pseudosegmentation = list('cluster' = psegmCoresListInDNA[[pos]][[paste0('Pseudo', n)]]$cluster,
#                                 'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
#       class(pseudosegmentation) = c('dbscan_fast', 'dbscan')
#       tmp[[paste0('Pseudo_in_DNA', n)]] =
#         predict(object = pseudosegmentation,
#                 newdata = select(ungroup(tmp), x, y),
#                 data = select(ungroup(psegmCoresListInDNA[[pos]][[paste0('Pseudo', n)]]), x, y))
#       #remove parts of clusters that don't overlap with DNA foci
#       tmp %>% mutate(!!paste0('Pseudo_in_DNA', n) :=
#                        case_when(DNA_focus_idx == 0 ~ as.integer(0),
#                                  T ~ !!as.name(paste0('Pseudo_in_DNA', n)))) -> tmp
#     }
#     
#     #Classify by HP1
#     segmentation = list('cluster' = segmCoresList[[pos]][['Foci']]$cluster,
#                         'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
#     class(segmentation) = c('dbscan_fast', 'dbscan')
#     tmp$HP1_foci_dbscan = predict(object = segmentation,
#                                   newdata = select(ungroup(tmp), x, y),
#                                   data = select(ungroup(segmCoresList[[pos]][['Foci']]), x, y))
#     for (n in 1:n_shuffle) {
#       pseudosegmentation = list('cluster' = segmCoresList[[pos]][[paste0('Pseudo', n)]]$cluster,
#                                 'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
#       class(pseudosegmentation) = c('dbscan_fast', 'dbscan')
#       tmp[[paste0('HP1_pseudo_dbscan', n)]] =
#         predict(object = pseudosegmentation,
#                 newdata = select(ungroup(tmp), x, y),
#                 data = select(ungroup(segmCoresList[[pos]][[paste0('Pseudo', n)]]), x, y))
#     }
#     
#     #Overlap
#     tmp$DNA_HP1_foci = tmp$HP1_foci_dbscan %in% DNA_HP1_foci[[pos]]
#     for (n in 1:n_shuffle) {
#       tmp[[paste0('DNA_HP1_pseudo', n)]] = tmp[[paste0('HP1_pseudo_dbscan', n)]] %in% DNA_HP1_foci[[pos]]
#     }
#     
#     #Concatenate
#     if (i == 1) {
#       result = tmp
#     } else {
#       result = rbind(result, tmp)
#     }
#   }
#   toc()
#   return(result)
# }) %>% bind_rows() -> jumpDist_1234
# saveRDS(jumpDist_1234, file = 'jumpDist_1234.rds')
# 
# #Assign class, jump dist curves
# {
#   jumpDist_1234 %>% ungroup() %>% #filter(!(DNA_focus & !pseudo_in_DNA1 & !HP1_foci_dbscan)) %>%
#     mutate(Compartment = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
#                                    # (Pseudo_in_DNA1!=0 | Pseudo_in_DNA2!=0 | Pseudo_in_DNA3!=0 |
#                                    # Pseudo_in_DNA4!=0 | Pseudo_in_DNA5!=0 | Pseudo_in_DNA6!=0)
#                                    # & !DNA_HP1_foci ~ 'DNA-high (subsampled)',
#                                    DNA_focus_idx!=0 & !DNA_HP1_foci ~ 'DNA-high',
#                                    !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
#                                    T ~ 'Euchromatin')) %>%
#     mutate(Compartment2 = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
#                                     DNA_focus_idx!=0 & !DNA_HP1_foci ~ 'DNA-high',
#                                     # Pseudo_in_DNA1 & !DNA_HP1_foci ~ 'DNA-high (subsampled 1)',
#                                     # Pseudo_in_DNA2 & !DNA_HP1_foci ~ 'DNA-high (subsampled 2)',
#                                     # Pseudo_in_DNA3 & !DNA_HP1_foci ~ 'DNA-high (subsampled 3)',
#                                     # Pseudo_in_DNA4 & !DNA_HP1_foci ~ 'DNA-high (subsampled 4)',
#                                     # Pseudo_in_DNA5 & !DNA_HP1_foci ~ 'DNA-high (subsampled 5)',
#                                     # Pseudo_in_DNA6 & !DNA_HP1_foci ~ 'DNA-high (subsampled 6)',
#                                     !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
#                                     DNA_pseudo & !(HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 |
#                                                      HP1_pseudo_dbscan3 | HP1_pseudo_dbscan4 |
#                                                      HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'DNA pseudo',
#                                     DNA_HP1_pseudo1 ~ 'DNA&HP1 pseudo 1',
#                                     DNA_HP1_pseudo2 ~ 'DNA&HP1 pseudo 2',
#                                     DNA_HP1_pseudo3 ~ 'DNA&HP1 pseudo 3',
#                                     DNA_HP1_pseudo4 ~ 'DNA&HP1 pseudo 4',
#                                     DNA_HP1_pseudo5 ~ 'DNA&HP1 pseudo 5',
#                                     DNA_HP1_pseudo6 ~ 'DNA&HP1 pseudo 6',
#                                     !DNA_HP1_pseudo1 & HP1_pseudo_dbscan1 ~ 'HP1 pseudo 1',#NB! These overlap -> pseudo6 ends up smaller
#                                     !DNA_HP1_pseudo2 & HP1_pseudo_dbscan2 ~ 'HP1 pseudo 2',
#                                     !DNA_HP1_pseudo3 & HP1_pseudo_dbscan3 ~ 'HP1 pseudo 3',
#                                     !DNA_HP1_pseudo4 & HP1_pseudo_dbscan4 ~ 'HP1 pseudo 4',
#                                     !DNA_HP1_pseudo5 & HP1_pseudo_dbscan5 ~ 'HP1 pseudo 5',
#                                     !DNA_HP1_pseudo6 & HP1_pseudo_dbscan6 ~ 'HP1 pseudo 6',
#                                     T ~ 'Euchromatin')) -> jumpDist_1234
#   filter(jumpDist_1234, Track_length < 15, interval == 1) %>% #sanity check, reproducing JDA
#     mutate(colour = case_when(grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
#                               grepl('HP1 pseudo', Compartment2) ~ 'HP1 pseudo',
#                               T ~ Compartment2)) %>%
#     mutate(colour = factor(colour, levels = c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                               'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                               'Euchromatin'))) %>%
#     ggplot(aes(x = disp, colour = colour, group = Compartment2)) +
#     geom_density() +
#     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                    'Euchromatin' = 'gray')) +
#     labs(title = 'Jump distance distribution, interval = 3 ms', subtitle = 'Short traj (<15)',
#          x = 'Displacement, nm') +
#     lims(x = c(0, 400)) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   
#   ggplot(jumpDist_1234, aes(x = disp, colour = factor(interval*3, levels = c('3', '6', '9', '12',
#                                                                              '15', '18', '21', '24')))) +
#     geom_density() +
#     facet_wrap(~Compartment, ncol = 2) +
#     scale_colour_manual(values = gray.colors(8, 0.2, 0.75)) +
#     labs(title = 'Jump distance distribution',
#          x = 'Displacement, nm', colour = 'Interval, ms') +
#     lims(x = c(0, 1000)) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   
#   filter(jumpDist_1234, Track_length < 15) %>% #looks somwehat more reasonable -> use this for MSD
#     ggplot(aes(x = disp, colour = factor(interval*3, levels = c('3', '6', '9', '12',
#                                                                 '15', '18', '21', '24')))) +
#     geom_density() +
#     facet_wrap(~Compartment, ncol = 2) +
#     scale_colour_manual(values = gray.colors(8, 0.2, 0.75)) +
#     labs(title = 'Jump distance distribution', subtitle = 'Short traj (<15)',
#          x = 'Displacement, nm', colour = 'Interval, ms') +
#     lims(x = c(0, 1000)) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   ggsave('New_analysis/JumpDist_intervals_shortTraj_HP1andDNA.png', width = 5, height = 4)
#   
#   filter(jumpDist_1234, Track_length < 15) %>% #looks somwehat more reasonable -> use this for MSD
#     ggplot(aes(x = disp, colour = Compartment)) +
#     geom_density() +
#     facet_wrap(~factor(interval*3, levels = c('3', '6', '9', '12',
#                                               '15', '18', '21', '24')), ncol = 3) +
#     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                    'Euchromatin' = 'gray')) +
#     labs(title = 'Jump distance distribution', subtitle = 'Short traj (<15)',
#          x = 'Displacement, nm', colour = NULL) +
#     lims(x = c(0, 1000)) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   ggsave('New_analysis/JumpDist_intervals_shortTraj_HP1andDNA2.png', width = 8, height = 6)
#   
#   filter(jumpDist_1234, Track_length < 15) %>% #looks somwehat more reasonable -> use this for MSD
#     mutate(colour = case_when(grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
#                               grepl('HP1 pseudo', Compartment2) ~ 'HP1 pseudo',
#                               T ~ Compartment2)) %>%
#     ggplot(aes(x = disp, colour = factor(interval*3, levels = c('3', '6', '9', '12',
#                                                                 '15', '18', '21', '24')))) +
#     geom_density() +
#     facet_wrap(~colour, ncol = 2) +
#     scale_colour_manual(values = gray.colors(8, 0.2, 0.75)) +
#     labs(title = 'Jump distance distribution', subtitle = 'Short traj (<15)',
#          x = 'Displacement, nm', colour = 'Interval, ms') +
#     lims(x = c(0, 1000)) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   ggsave('New_analysis/JumpDist_intervals_shortTraj_HP1andDNA_all.png', width = 6, height = 6)
# }
# 
# #Num obs
# {
#   filter(jumpDist_1234, Track_length < 15) %>% group_by(Compartment2, interval) %>% summarise(num_obs = n()) %>%
#     mutate(Compartment = case_when(grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
#                                    grepl('HP1 pseudo', Compartment2) ~ 'HP1 pseudo',
#                                    T ~ Compartment2)) %>%
#     mutate(Compartment2 = factor(Compartment2, levels = rev(c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                                               'DNA pseudo', paste0('HP1 pseudo ', 1:6),
#                                                               paste0('DNA&HP1 pseudo ', 1:6),
#                                                               'Euchromatin')))) %>%
#     ggplot(aes(x = Compartment2, y = num_obs, fill = Compartment, group = Compartment2, label = num_obs)) +
#     geom_bar(stat = 'identity', colour = 'gray20') +
#     geom_text(y = 15000, hjust = 1) +
#     facet_wrap(~interval) +
#     coord_flip(ylim = c(0, 15000)) +
#     scale_fill_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                  'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                  'Euchromatin' = 'gray')) +
#     labs(title = 'Number of observations', subtitle = 'HP1 segm with DBSCAN',
#          fill = 'Compartment', x = '', y = 'Num obs') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1, legend.position = 'none')
#   ggsave('New_analysis/MSD_HP1andDNA_numObs.png', width = 9, height = 9)
# }
# 
# #Calculate MSD (incl pooled pseudos), plot (lin and log)
# {
#   ##the pseudo reps separately
#   {
#     filter(jumpDist_1234, Track_length < 15) %>% group_by(Compartment2, interval) %>%
#       summarise(MSD = mean(disp^2)/1000000, N = n()) -> MSD_segm
#     
#     filter(MSD_segm, !grepl('pseudo', Compartment2)) %>%
#       ggplot(aes(x = interval*3, y = MSD, colour = Compartment2)) +
#       geom_point() +
#       geom_line(size = 0.2) +
#       scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'Euchromatin' = 'gray')) +
#       labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#            x = 'Interval, ms', y = 'MSD, um2') +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 1)
#     ggsave('New_analysis/MSD_HP1andDNA_only.png', width = 5.5, height = 4)
#     
#     filter(jumpDist_1234, Track_length < 15) %>% group_by(Dish, Compartment2, interval) %>%
#       summarise(MSD = mean(disp^2)/1000000, N = n()) -> MSD_segm_dish
#     
#     filter(MSD_segm_dish, !grepl('pseudo', Compartment2)) %>%
#       ggplot(aes(x = interval*3, y = MSD, colour = Compartment2, group = interaction(Dish, Compartment2))) +
#       geom_point() +
#       geom_line(size = 0.2) +
#       scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'Euchromatin' = 'gray')) +
#       labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#            x = 'Interval, ms', y = 'MSD, um2') +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 1)
#     ggsave('New_analysis/MSD_HP1andDNA_only_byDish.png', width = 5.5, height = 4)
#     
#     mutate(MSD_segm_dish, colour = case_when(grepl('HP1 pseudo', Compartment2) & !grepl('DNA', Compartment2) ~ 'HP1 pseudo',
#                                         grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
#                                         T ~ Compartment2)) %>% 
#       mutate(colour = factor(colour, levels = c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                                 'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                                 'Euchromatin'))) %>%
#       ggplot(aes(x = interval*3, y = MSD, group = interaction(Compartment2, N>30),
#                  colour = colour, shape = colour, size = colour,
#                  alpha = N > 30, linetype = N > 30)) +
#       geom_point() +
#       geom_line(size = 0.2) +
#       facet_grid(Dish~.) +
#       scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray')) +
#       scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#       scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2)) +
#       scale_alpha_manual(values = c(0.5, 1)) +
#       scale_linetype_manual(values = c('dashed', 'solid')) +
#       coord_cartesian(ylim = c(NA, 0.09)) +
#       labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#            x = 'Interval, ms', y = 'MSD, um2',
#            colour = 'Compartment', shape = 'Compartment', size = 'Compartment') +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 1)
#     ggsave('New_analysis/MSD_HP1andDNA.png', width = 5.5, height = 6)
#     
#     filter(MSD_segm, !grepl('pseudo', Compartment2)) %>%
#       ggplot(aes(x = log10(interval*3/1000), y = log10(MSD), colour = Compartment2)) +
#       geom_point() +
#       geom_line(size = 0.2) +
#       scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'Euchromatin' = 'gray')) + #'steelblue' #for DNA-high subsampled
#       #coord_cartesian(ylim = c(-2.05, -1.2)) +
#       coord_cartesian(y = c(-2.15, -1.05)) +
#       labs(title = 'MSD curve, log-log', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#            x = 'Log (Interval, ms)', y = 'Log (MSD, um2)') +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 1)
#     ggsave('New_analysis/MSD_HP1andDNA_only_loglogLine.png', width = 5.5, height = 4)
#     filter(MSD_segm_dish, !grepl('pseudo', Compartment2)) %>%
#       ggplot(aes(x = log10(interval*3/1000), y = log10(MSD), colour = Compartment2,
#                  group = interaction(Compartment2, Dish))) +
#       geom_point() +
#       geom_line(size = 0.2) +
#       scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'Euchromatin' = 'gray')) + #'steelblue' #for DNA-high subsampled
#       #coord_cartesian(ylim = c(-2.05, -1.2)) +
#       coord_cartesian(y = c(-2.15, -1.05)) +
#       labs(title = 'MSD curve, log-log', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#            x = 'Log (Interval, ms)', y = 'Log (MSD, um2)') +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 1)
#     ggsave('New_analysis/MSD_HP1andDNA_only_loglogLine_byDish.png', width = 5.5, height = 4)
#     
#     mutate(MSD_segm_dish, colour = case_when(grepl('HP1 pseudo', Compartment2) & !grepl('DNA', Compartment2) ~ 'HP1 pseudo',
#                                         grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
#                                         T ~ Compartment2)) %>% 
#       mutate(colour = factor(colour, levels = c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                                 'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                                 'Euchromatin'))) %>%
#       ggplot(aes(x = log10(interval*3/1000), y = log10(MSD), group = interaction(Compartment2, N>30),
#                  colour = colour, shape = colour, size = colour,
#                  alpha = N > 30, linetype = N > 30)) +
#       geom_point() +
#       geom_line(size = 0.2) +
#       facet_grid(Dish~.) +
#       scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray')) + #'steelblue' #for DNA-high subsampled
#       scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#       scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2)) +
#       scale_alpha_manual(values = c(0.5, 1)) +
#       scale_linetype_manual(values = c('dashed', 'solid')) +
#       coord_cartesian(ylim = c(-2.15, -1.05)) +
#       labs(title = 'MSD curve, log-log', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#            x = 'Log (Interval, ms)', y = 'Log (MSD, um2)',
#            colour = 'Compartment', shape = 'Compartment', size = 'Compartment') +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 1)
#     ggsave('New_analysis/MSD_HP1andDNA_loglogLine.png', width = 5.5, height = 6)
#   }
#   
#   ##the pseudo reps together
#   {
#     filter(jumpDist_1234, Track_length < 15) %>%
#       mutate(Compartment2 = case_when(grepl('HP1 pseudo', Compartment2) & !grepl('DNA', Compartment2) ~ 'HP1 pseudo',
#                                       grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
#                                       T ~ Compartment2)) %>%
#       group_by(Compartment2, interval) %>%
#       summarise(MSD = mean(disp^2)/1000000, N = n()) %>%
#       mutate(Compartment2 = factor(Compartment2, levels = c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                                             'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                                             'Euchromatin'))) -> MSD_segm_pooled
#     filter(jumpDist_1234, Track_length < 15) %>%
#       mutate(Compartment2 = case_when(grepl('HP1 pseudo', Compartment2) & !grepl('DNA', Compartment2) ~ 'HP1 pseudo',
#                                       grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
#                                       T ~ Compartment2)) %>%
#       group_by(Dish, Compartment2, interval) %>%
#       summarise(MSD = mean(disp^2)/1000000, N = n()) %>%
#       mutate(Compartment2 = factor(Compartment2, levels = c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                                             'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                                             'Euchromatin'))) -> MSD_segm_pooled_byDish
#     ggplot(MSD_segm_pooled, aes(x = interval*3, y = MSD,
#                                 colour = Compartment2, shape = Compartment2, size = Compartment2)) +
#       #alpha = N > 30, linetype = N > 30)) +
#       geom_point() +
#       geom_line(size = 0.2) +
#       scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray')) +
#       scale_shape_manual(values = c(19, 19, 19, 4, 4,  4, 19)) +
#       scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2)) +
#       #scale_alpha_manual(values = c(0.5, 1)) +
#       #scale_linetype_manual(values = c('dashed', 'solid')) +
#       labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#            x = 'Interval, ms', y = 'MSD, um2') +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 1)
#     ggsave('New_analysis/MSD_HP1andDNA_poolPseudo.png', width = 5.5, height = 4)
#     ggplot(MSD_segm_pooled_byDish, aes(x = interval*3, y = MSD,
#                                 colour = Compartment2, shape = Compartment2, size = Compartment2,
#                                 group = interaction(Dish, Compartment2))) +
#       #alpha = N > 30, linetype = N > 30)) +
#       geom_point() +
#       geom_line(size = 0.2) +
#       scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray')) +
#       scale_shape_manual(values = c(19, 19, 19, 4, 4,  4, 19)) +
#       scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2)) +
#       #scale_alpha_manual(values = c(0.5, 1)) +
#       #scale_linetype_manual(values = c('dashed', 'solid')) +
#       labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#            x = 'Interval, ms', y = 'MSD, um2') +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 1)
#     ggsave('New_analysis/MSD_HP1andDNA_poolPseudo_byDish.png', width = 5.5, height = 4)
#     ggplot(MSD_segm_pooled, aes(x = log10(interval*3), y = log10(MSD),
#                                 colour = Compartment2, shape = Compartment2, size = Compartment2)) +
#       #alpha = N > 30, linetype = N > 30)) +
#       geom_point() +
#       geom_line(size = 0.2) +
#       scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray')) +
#       scale_shape_manual(values = c(19, 19, 19, 4, 4,  4, 19)) +
#       scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2)) +
#       #scale_alpha_manual(values = c(0.5, 1)) +
#       #scale_linetype_manual(values = c('dashed', 'solid')) +
#       coord_cartesian(y = c(-2.15, -1.05)) +
#       labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#            x = 'Log (Interval, ms)', y = 'Log (MSD, um2)') +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 1)
#     ggsave('New_analysis/MSD_HP1andDNA_loglogLine_poolPseudo.png', width = 5.5, height = 4)
#   }
# }
# 
# 
# #Fit model to MSD curve, plot lines in log-log plots, plot summary
# {
#   #All data, pool the pseudo
#   split(MSD_segm_pooled, MSD_segm_pooled$Compartment2) %>%
#     lapply(function(x) {
#       lm(log10(MSD)~log10(interval*3/1000), data = x)
#     }) -> MSDfit_segm.list
#   MSDfit_segm.df = data.frame(Compartment2 = names(MSDfit_segm.list),
#                               Intercept = sapply(MSDfit_segm.list, function(x) {x$coefficients[[1]]}),
#                               Slope = sapply(MSDfit_segm.list, function(x) {x$coefficients[[2]]}))
#   mutate(MSDfit_segm.df, label = paste('D =', round(10^Intercept/4, 2), 'um2/s, alpha =', round(Slope, 2), '\n')) -> MSDfit_segm.df
#   MSDfit_segm.df %>%
#     mutate(Compartment2 = factor(Compartment2, levels = c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                                           'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                                           'Euchromatin'))) -> MSDfit_segm.df
#   ggplot(MSD_segm_pooled, aes(x = log10(interval*3/1000), y = log10(MSD),
#                               colour = Compartment2, shape = Compartment2, size = Compartment2)) +
#     #alpha = N > 30 & !(Compartment2 != 'DNA&HP1 pseudo' & N<100) &
#     #!(Compartment2 == 'DNA-high (subsampled)' & interval == 8))) +
#     geom_point() +
#     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                    'Euchromatin' = 'gray'),
#                         labels = paste0(MSDfit_segm.df$Compartment2, '\n', MSDfit_segm.df$label)) +
#     scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19),
#                        labels = paste0(MSDfit_segm.df$Compartment2, '\n', MSDfit_segm.df$label)) +
#     scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2),
#                       labels = paste0(MSDfit_segm.df$Compartment2, '\n', MSDfit_segm.df$label)) +
#     #scale_alpha_manual(values = c(0.3, 1), guide = 'none') +
#     geom_abline(data = MSDfit_segm.df, aes(intercept = Intercept, slope = Slope, colour = Compartment2)) +
#     labs(title = 'MSD fit', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#          x = 'Log(Interval, s)', y = 'Log(MSD, um2)') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1)
#   ggsave('New_analysis/MSD_HP1andDNA_loglogFit.png', width = 5.5, height = 4)
#   
#   ggplot(MSDfit_segm.df, aes(x = Compartment2, y = Intercept-log10(4), colour = Compartment2, shape = Compartment2)) +
#     geom_point() +
#     scale_y_continuous(labels = math_format(10^.x),
#                        sec.axis = sec_axis(trans=~10^., breaks = c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2))) +
#     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                    'Euchromatin' = 'gray')) +
#     scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#     labs(x = '', y = 'D, um2/s') +
#     theme_bw() +
#     theme(aspect.ratio = 0.8, legend.position = 'none', axis.text.x = element_text(angle = 30, hjust = 1)) -> Ds
#   ggplot(MSDfit_segm.df, aes(x = Compartment2, y = Slope, colour = Compartment2, shape = Compartment2)) +
#     geom_point() +
#     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                    'Euchromatin' = 'gray')) +
#     scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#     scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
#     labs(x = '', y = 'Alpha') +
#     theme_bw() +
#     theme(aspect.ratio = 0.8, legend.position = 'none', axis.text.x = element_text(angle = 30, hjust = 1)) -> alphas
#   ggarrange(Ds, alphas, nrow = 2, ncol = 1)
#   ggsave('New_analysis/D_alpha_HP1andDNA.png', width = 3, height = 5)
# }
# 
# #Same, but only using intervals 1-5
# {
#   #All data, pool the pseudo
#   split(filter(MSD_segm_pooled, interval < 6), filter(MSD_segm_pooled, interval < 6)$Compartment2) %>%
#     lapply(function(x) {
#       lm(log10(MSD)~log10(interval*3/1000), data = x)
#     }) -> MSDfit_segm.list2
#   MSDfit_segm.df2 = data.frame(Compartment2 = names(MSDfit_segm.list2),
#                                Intercept = sapply(MSDfit_segm.list2, function(x) {x$coefficients[[1]]}),
#                                Slope = sapply(MSDfit_segm.list2, function(x) {x$coefficients[[2]]}))
#   mutate(MSDfit_segm.df2, label = paste('D =', round(10^Intercept/4, 2), 'um2/s, alpha =', round(Slope, 2), '\n')) -> MSDfit_segm.df2
#   MSDfit_segm.df2 %>%
#     mutate(Compartment2 = factor(Compartment2, levels = c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                                           'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                                           'Euchromatin'))) -> MSDfit_segm.df2
#   ggplot(MSD_segm_pooled, aes(x = log10(interval*3/1000), y = log10(MSD),
#                               colour = Compartment2, shape = Compartment2, size = Compartment2,
#                               alpha = interval < 6)) +
#     #alpha = N > 30 & !(Compartment2 != 'DNA&HP1 pseudo' & N<100) &
#     #!(Compartment2 == 'DNA-high (subsampled)' & interval == 8))) +
#     geom_point() +
#     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                    'Euchromatin' = 'gray'),
#                         labels = paste0(MSDfit_segm.df2$Compartment2, '\n', MSDfit_segm.df2$label)) +
#     scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19),
#                        labels = paste0(MSDfit_segm.df2$Compartment2, '\n', MSDfit_segm.df2$label)) +
#     scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2),
#                       labels = paste0(MSDfit_segm.df2$Compartment2, '\n', MSDfit_segm.df2$label)) +
#     scale_alpha_manual(values = c(0.3, 1), guide = 'none') +
#     geom_abline(data = MSDfit_segm.df2, aes(intercept = Intercept, slope = Slope, colour = Compartment2)) +
#     labs(title = 'MSD fit', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#          x = 'Log(Interval, s)', y = 'Log(MSD, um2)') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1)
#   ggsave('New_analysis/MSD_HP1andDNA_loglogFit2.png', width = 5.5, height = 4)
#   
#   ggplot(MSDfit_segm.df2, aes(x = Compartment2, y = Intercept-log10(4), colour = Compartment2, shape = Compartment2)) +
#     geom_point() +
#     scale_y_continuous(labels = math_format(10^.x),
#                        sec.axis = sec_axis(trans=~10^., breaks = c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2))) +
#     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                    'Euchromatin' = 'gray')) +
#     scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#     labs(x = '', y = 'D, um2/s') +
#     theme_bw() +
#     theme(aspect.ratio = 0.8, legend.position = 'none', axis.text.x = element_text(angle = 30, hjust = 1)) -> Ds
#   ggplot(MSDfit_segm.df2, aes(x = Compartment2, y = Slope, colour = Compartment2, shape = Compartment2)) +
#     geom_point() +
#     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                    'Euchromatin' = 'gray')) +
#     scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#     scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
#     labs(x = '', y = 'Alpha') +
#     theme_bw() +
#     theme(aspect.ratio = 0.8, legend.position = 'none', axis.text.x = element_text(angle = 30, hjust = 1)) -> alphas
#   ggarrange(Ds, alphas, nrow = 2, ncol = 1)
#   ggsave('New_analysis/D_alpha_HP1andDNA2.png', width = 3, height = 5)
# }
# 
# #Fit model to bootstrapped MSD curves (using intervals 1-5),
# #plot lines in log-log plots and summary with error bars
# {
#   #bootstrap
#   MSD_segm_pooled_bs = list()
#   tic()
#   lapply(1:100, function(n) {
#     #calc time weird, now much faster than before
#     #100 with lapply - 150s
#     filter(jumpDist_1234, Track_length < 15) %>%
#       mutate(Compartment2 = case_when(grepl('HP1 pseudo', Compartment2) & !grepl('DNA', Compartment2) ~ 'HP1 pseudo',
#                                       grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
#                                       T ~ Compartment2)) %>%
#       group_by(Compartment2, interval) %>% slice_sample(prop = 1, replace = T) %>%
#       summarise(MSD = mean(disp^2)/1000000, N = n()) -> MSD.new
#     return(MSD.new)
#   }) -> MSD_segm_pooled_bs
#   toc()
#   
#   bind_rows(MSD_segm_pooled_bs) %>% group_by(Compartment2, interval) %>%
#     summarise(mean = mean(MSD), sd = sd(MSD),
#               CI_l = quantile(MSD, 0.05), CI_h = quantile(MSD, 0.95)) -> MSD_segm_bs_stats
#   
#   
#   #fit
#   tic()
#   lapply(MSD_segm_pooled_bs, function(s) {
#     split(filter(s, interval < 6), filter(s, interval < 6)$Compartment2) %>%
#       lapply(function(x) {
#         lm(log10(MSD)~log10(interval*3/1000), data = x)
#       }) -> fit.list
#     fit.df = data.frame(Compartment2 = names(fit.list),
#                         Intercept = sapply(fit.list, function(x) {x$coefficients[[1]]}),
#                         Slope = sapply(fit.list, function(x) {x$coefficients[[2]]}))
#     return(fit.df)
#   }) -> MSDfit_segm_bs.df
#   toc()
#   
#   bind_rows(MSDfit_segm_bs.df) %>% group_by(Compartment2) %>%
#     summarise(int.mean = mean(Intercept), int.sd = sd(Intercept),
#               int.CI_l = quantile(Intercept, 0.05), int.CI_h = quantile(Intercept, 0.95),
#               slope.mean = mean(Slope), slope.sd = sd(Slope),
#               slope.CI_l = quantile(Slope, 0.05), slope.CI_h = quantile(Slope, 0.95)) -> MSDfit_segm_bs_stats
#   
#   ggplot(MSD_segm_pooled, aes(x = log10(interval*3/1000), y = log10(MSD),
#                               colour = Compartment2, shape = Compartment2, size = Compartment2,
#                               alpha = interval < 6)) +
#     geom_point() +
#     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                    'Euchromatin' = 'gray'),
#                         labels = paste0(MSDfit_segm.df2$Compartment2, '\n', MSDfit_segm.df2$label)) +
#     scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19),
#                        labels = paste0(MSDfit_segm.df2$Compartment2, '\n', MSDfit_segm.df2$label)) +
#     scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2),
#                       labels = paste0(MSDfit_segm.df2$Compartment2, '\n', MSDfit_segm.df2$label)) +
#     scale_alpha_manual(values = c(0.3, 1), guide = 'none') +
#     geom_errorbar(data = MSD_segm_bs_stats, aes(ymin = log10(CI_l), ymax = log10(CI_h),
#                                                 y = NULL),
#                   width = 0.01, size = 0.5) +
#     geom_abline(data = MSDfit_segm.df2, aes(intercept = Intercept, slope = Slope, colour = Compartment2)) +
#     labs(title = 'MSD fit', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#          x = 'Log(Interval, s)', y = 'Log(MSD, um2)') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1)
#   ggsave('New_analysis/MSD_HP1andDNA_loglogFit_errors.png', width = 5.5, height = 4)
#   
#   ggplot(MSDfit_segm.df2, aes(x = Compartment2, y = Intercept-log10(4), colour = Compartment2, shape = Compartment2)) +
#     geom_point() +
#     scale_y_continuous(labels = math_format(10^.x),
#                        sec.axis = sec_axis(trans=~10^., breaks = c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2))) +
#     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                    'Euchromatin' = 'gray')) +
#     scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#     geom_errorbar(data = MSDfit_segm_bs_stats, aes(ymin = int.CI_l-log10(4), ymax = int.CI_h-log10(4),
#                                                    y = NULL, alpha = NULL),
#                   width = 0.1, size = 0.3) +
#     labs(x = '', y = 'D, um2/s') +
#     theme_bw() +
#     theme(aspect.ratio = 0.8, legend.position = 'none', axis.text.x = element_text(angle = 30, hjust = 1)) -> Ds
#   
#   ggplot(MSDfit_segm.df2, aes(x = Compartment2, y = Slope, colour = Compartment2, shape = Compartment2)) +
#     geom_point() +
#     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                    'Euchromatin' = 'gray')) +
#     scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#     geom_errorbar(data = MSDfit_segm_bs_stats, aes(ymin = slope.CI_l, ymax = slope.CI_h,
#                                                    y = NULL, alpha = NULL),
#                   width = 0.1, size = 0.3) +
#     scale_y_continuous(lim = c(-0.02, 1), breaks = seq(0, 1, by = 0.2)) +
#     labs(x = '', y = 'Alpha') +
#     theme_bw() +
#     theme(aspect.ratio = 0.8, legend.position = 'none', axis.text.x = element_text(angle = 30, hjust = 1)) -> alphas
#   ggarrange(Ds, alphas, nrow = 2, ncol = 1)
#   ggsave('New_analysis/D_alpha_HP1andDNA_errors.png', width = 3, height = 5)
# }
# 
# #Repeated EVERYTHING with subsampled DNA
# {
#   {
#     jumpDist_1234 %>% ungroup() %>% #filter(!(DNA_focus & !pseudo_in_DNA1 & !HP1_foci_dbscan)) %>%
#       mutate(Compartment = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
#                                      (Pseudo_in_DNA1!=0 | Pseudo_in_DNA2!=0 | Pseudo_in_DNA3!=0 |
#                                      Pseudo_in_DNA4!=0 | Pseudo_in_DNA5!=0 | Pseudo_in_DNA6!=0)
#                                      & !DNA_HP1_foci ~ 'DNA-high (subsampled)',
#                                      #DNA_focus_idx!=0 & !DNA_HP1_foci ~ 'DNA-high',
#                                      !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
#                                      T ~ 'Euchromatin')) %>%
#       mutate(Compartment2 = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
#                                       #DNA_focus_idx!=0 & !DNA_HP1_foci ~ 'DNA-high',
#                                       Pseudo_in_DNA1 & !DNA_HP1_foci ~ 'DNA-high (subsampled 1)',
#                                       Pseudo_in_DNA2 & !DNA_HP1_foci ~ 'DNA-high (subsampled 2)',
#                                       Pseudo_in_DNA3 & !DNA_HP1_foci ~ 'DNA-high (subsampled 3)',
#                                       Pseudo_in_DNA4 & !DNA_HP1_foci ~ 'DNA-high (subsampled 4)',
#                                       Pseudo_in_DNA5 & !DNA_HP1_foci ~ 'DNA-high (subsampled 5)',
#                                       Pseudo_in_DNA6 & !DNA_HP1_foci ~ 'DNA-high (subsampled 6)',
#                                       !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
#                                       DNA_pseudo & !(HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 |
#                                                        HP1_pseudo_dbscan3 | HP1_pseudo_dbscan4 |
#                                                        HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'DNA pseudo',
#                                       DNA_HP1_pseudo1 ~ 'DNA&HP1 pseudo 1',
#                                       DNA_HP1_pseudo2 ~ 'DNA&HP1 pseudo 2',
#                                       DNA_HP1_pseudo3 ~ 'DNA&HP1 pseudo 3',
#                                       DNA_HP1_pseudo4 ~ 'DNA&HP1 pseudo 4',
#                                       DNA_HP1_pseudo5 ~ 'DNA&HP1 pseudo 5',
#                                       DNA_HP1_pseudo6 ~ 'DNA&HP1 pseudo 6',
#                                       !DNA_HP1_pseudo1 & HP1_pseudo_dbscan1 ~ 'HP1 pseudo 1',#NB! These overlap -> pseudo6 ends up smaller
#                                       !DNA_HP1_pseudo2 & HP1_pseudo_dbscan2 ~ 'HP1 pseudo 2',
#                                       !DNA_HP1_pseudo3 & HP1_pseudo_dbscan3 ~ 'HP1 pseudo 3',
#                                       !DNA_HP1_pseudo4 & HP1_pseudo_dbscan4 ~ 'HP1 pseudo 4',
#                                       !DNA_HP1_pseudo5 & HP1_pseudo_dbscan5 ~ 'HP1 pseudo 5',
#                                       !DNA_HP1_pseudo6 & HP1_pseudo_dbscan6 ~ 'HP1 pseudo 6',
#                                       T ~ 'Euchromatin')) -> jumpDist_1234
#     filter(jumpDist_1234, Track_length < 15, interval == 1) %>% #sanity check, reproducing JDA
#       mutate(colour = case_when(grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
#                                 grepl('HP1 pseudo', Compartment2) ~ 'HP1 pseudo',
#                                 #grepl('DNA-high', Compartment2) ~ 'DNA-high (subsampled)',
#                                 T ~ Compartment2)) %>%
#       mutate(colour = factor(colour, levels = c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                                 'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                                 'Euchromatin'))) %>%
#       ggplot(aes(x = disp, colour = colour, group = Compartment2)) +
#       geom_density() +
#       scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray')) +
#       # scale_colour_manual(values = c('DNA-high (subsampled)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#       #                                'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#       #                                'Euchromatin' = 'gray')) +
#       labs(title = 'Jump distance distribution, interval = 3 ms', subtitle = 'Short traj (<15)',
#            x = 'Displacement, nm') +
#       lims(x = c(0, 400)) +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 0.5, axis.title.y = element_blank())
#     
#     filter(jumpDist_1234, Track_length < 15) %>% #looks somwehat more reasonable -> use this for MSD
#       ggplot(aes(x = disp, colour = factor(interval*3, levels = c('3', '6', '9', '12',
#                                                                   '15', '18', '21', '24')))) +
#       geom_density() +
#       facet_wrap(~Compartment, ncol = 2) +
#       scale_colour_manual(values = gray.colors(8, 0.2, 0.75)) +
#       labs(title = 'Jump distance distribution', subtitle = 'Short traj (<15)',
#            x = 'Displacement, nm', colour = 'Interval, ms') +
#       lims(x = c(0, 1000)) +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 0.5, axis.title.y = element_blank())
#     #ggsave('New_analysis/JumpDist_intervals_shortTraj_HP1andDNA.png', width = 5, height = 4)
#   }
#   
#   #Num obs
#   {
#     filter(jumpDist_1234, Track_length < 15) %>% group_by(Compartment2, interval) %>% summarise(num_obs = n()) %>%
#       mutate(Compartment = case_when(grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
#                                      grepl('HP1 pseudo', Compartment2) ~ 'HP1 pseudo',
#                                      grepl('DNA-high', Compartment2) ~ 'DNA-high (subsampled)',
#                                      T ~ Compartment2)) %>%
#       mutate(Compartment2 = factor(Compartment2, levels =
#                                      #rev(c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                      rev(c(paste0('DNA-high (subsampled ', 1:6, ')'), 'HP1-high', 'DNA&HP1-high',
#                                            'DNA pseudo', paste0('HP1 pseudo ', 1:6),
#                                            paste0('DNA&HP1 pseudo ', 1:6),
#                                            'Euchromatin')))) %>%
#       ggplot(aes(x = Compartment2, y = num_obs, fill = Compartment, group = Compartment2, label = num_obs)) +
#       geom_bar(stat = 'identity', colour = 'gray20') +
#       geom_text(y = 15000, hjust = 1) +
#       facet_wrap(~interval) +
#       coord_flip(ylim = c(0, 15000)) +
#       scale_fill_manual(values = c('DNA-high (subsampled)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                    'Euchromatin' = 'gray')) +
#       # scale_fill_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#       #                              'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#       #                              'Euchromatin' = 'gray')) +
#       labs(title = 'Number of observations', subtitle = 'HP1 segm with DBSCAN\nNB! DNA subsampled overlap a lot',
#            fill = 'Compartment', x = '', y = 'Num obs') +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 1, legend.position = 'none')
#     ggsave('New_analysis/MSD_HP1andDNA_numObs_DNAsub.png', width = 9, height = 9)
#   }
#   
#   #Calculate MSD (incl pooled pseudos), plot (lin and log)
#   {
#     ##the pseudo reps separately
#     {
#       filter(jumpDist_1234, Track_length < 15) %>%
#         mutate(Compartment2 = case_when(grepl('DNA-high', Compartment2) ~ 'DNA-high (subsampledx6)',
#                                         T ~ Compartment2)) %>% #pooling all the DNA subsampled
#         group_by(Compartment2, interval) %>%
#         summarise(MSD = mean(disp^2)/1000000, N = n()) -> MSD_segm
#       filter(jumpDist_1234, Track_length < 15) %>%
#         mutate(Compartment2 = case_when(grepl('DNA-high', Compartment2) ~ 'DNA-high (subsampledx6)',
#                                         T ~ Compartment2)) %>% #pooling all the DNA subsampled
#         group_by(Dish, Compartment2, interval) %>%
#         summarise(MSD = mean(disp^2)/1000000, N = n()) -> MSD_segm_byDish
#       
#       filter(MSD_segm, !grepl('pseudo', Compartment2)) %>%
#         ggplot(aes(x = interval*3, y = MSD, colour = Compartment2)) +
#         geom_point() +
#         geom_line(size = 0.2) +
#         scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                        'Euchromatin' = 'gray')) +
#         # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#         #                                'Euchromatin' = 'gray')) +
#         labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#              x = 'Interval, ms', y = 'MSD, um2') +
#         theme_bw() +
#         theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#               aspect.ratio = 1)
#       ggsave('New_analysis/MSD_HP1andDNA_only_DNAsub.png', width = 5.7, height = 4)
#       
#       filter(MSD_segm_byDish, !grepl('pseudo', Compartment2)) %>%
#         ggplot(aes(x = interval*3, y = MSD, colour = Compartment2, group = interaction(Dish, Compartment2))) +
#         geom_point() +
#         geom_line(size = 0.2) +
#         scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                        'Euchromatin' = 'gray')) +
#         # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#         #                                'Euchromatin' = 'gray')) +
#         labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#              x = 'Interval, ms', y = 'MSD, um2') +
#         theme_bw() +
#         theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#               aspect.ratio = 1)
#       ggsave('New_analysis/MSD_HP1andDNA_only_DNAsub_byDish.png', width = 5.7, height = 4)
#       
#       mutate(MSD_segm_byDish, colour = case_when(grepl('HP1 pseudo', Compartment2) & !grepl('DNA', Compartment2) ~ 'HP1 pseudo',
#                                           grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
#                                           T ~ Compartment2)) %>% 
#         mutate(colour = factor(colour, levels = c('DNA-high (subsampledx6)', 'HP1-high', 'DNA&HP1-high',
#                                                   'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                                   'Euchromatin'))) %>%
#         ggplot(aes(x = interval*3, y = MSD, group = interaction(Compartment2, N>30),
#                    colour = colour, shape = colour, size = colour,
#                    alpha = N > 30, linetype = N > 30)) +
#         geom_point() +
#         geom_line(size = 0.2) +
#         facet_grid(Dish~.) +
#         scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                        'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                        'Euchromatin' = 'gray')) +
#         # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#         #                                'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#         #                                'Euchromatin' = 'gray')) +
#         scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#         scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2)) +
#         scale_alpha_manual(values = c(0.5, 1)) +
#         scale_linetype_manual(values = c('dashed', 'solid')) +
#         coord_cartesian(ylim = c(NA, 0.09)) +
#         labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#              x = 'Interval, ms', y = 'MSD, um2',
#              colour = 'Compartment', shape = 'Compartment', size = 'Compartment') +
#         theme_bw() +
#         theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#               aspect.ratio = 1)
#       ggsave('New_analysis/MSD_HP1andDNA_DNAsub.png', width = 5.7, height = 6)
#       
#       filter(MSD_segm, !grepl('pseudo', Compartment2)) %>%
#         ggplot(aes(x = log10(interval*3/1000), y = log10(MSD), colour = Compartment2)) +
#         geom_point() +
#         geom_line(size = 0.2) +
#         # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#         #                                'Euchromatin' = 'gray')) +
#         scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                        'Euchromatin' = 'gray')) +
#         #coord_cartesian(ylim = c(-2.05, -1.2)) +
#         coord_cartesian(y = c(-2.15, -1.05)) +
#         labs(title = 'MSD curve, log-log', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#              x = 'Log (Interval, ms)', y = 'Log (MSD, um2)') +
#         theme_bw() +
#         theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#               aspect.ratio = 1)
#       ggsave('New_analysis/MSD_HP1andDNA_only_loglogLine_DNAsub.png', width = 5.7, height = 4)
#       
#       filter(MSD_segm_byDish, !grepl('pseudo', Compartment2)) %>%
#         ggplot(aes(x = log10(interval*3/1000), y = log10(MSD), colour = Compartment2, group = interaction(Compartment2, Dish))) +
#         geom_point() +
#         geom_line(size = 0.2) +
#         # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#         #                                'Euchromatin' = 'gray')) +
#         scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                        'Euchromatin' = 'gray')) +
#         #coord_cartesian(ylim = c(-2.05, -1.2)) +
#         coord_cartesian(y = c(-2.15, -1.05)) +
#         labs(title = 'MSD curve, log-log', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#              x = 'Log (Interval, ms)', y = 'Log (MSD, um2)') +
#         theme_bw() +
#         theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#               aspect.ratio = 1)
#       ggsave('New_analysis/MSD_HP1andDNA_only_loglogLine_DNAsub_byDish.png', width = 5.7, height = 4)
#       
#       mutate(MSD_segm_byDish, colour = case_when(grepl('HP1 pseudo', Compartment2) & !grepl('DNA', Compartment2) ~ 'HP1 pseudo',
#                                           grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
#                                           T ~ Compartment2)) %>% 
#         mutate(colour = factor(colour, levels = #c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                  c('DNA-high (subsampledx6)', 'HP1-high', 'DNA&HP1-high',
#                                    'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                    'Euchromatin'))) %>%
#         ggplot(aes(x = log10(interval*3/1000), y = log10(MSD), group = interaction(Compartment2, N>30),
#                    colour = colour, shape = colour, size = colour,
#                    alpha = N > 30, linetype = N > 30)) +
#         geom_point() +
#         geom_line(size = 0.2) +
#         facet_grid(Dish~.) +
#         # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#         #                                'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#         #                                'Euchromatin' = 'gray')) +
#         scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                        'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                        'Euchromatin' = 'gray')) +
#         scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#         scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2)) +
#         scale_alpha_manual(values = c(0.5, 1)) +
#         scale_linetype_manual(values = c('dashed', 'solid')) +
#         coord_cartesian(ylim = c(-2.15, -1.05)) +
#         labs(title = 'MSD curve, log-log', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#              x = 'Log (Interval, ms)', y = 'Log (MSD, um2)',
#              colour = 'Compartment', shape = 'Compartment', size = 'Compartment') +
#         theme_bw() +
#         theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#               aspect.ratio = 1)
#       ggsave('New_analysis/MSD_HP1andDNA_loglogLine_subDNA.png', width = 5.7, height = 6)
#     }
#     
#     ##the pseudo reps together
#     {
#       filter(jumpDist_1234, Track_length < 15) %>%
#         mutate(Compartment2 = case_when(grepl('DNA-high', Compartment2) ~ 'DNA-high (subsampledx6)',
#                                         T ~ Compartment2)) %>% #pooling all the DNA subsampled
#         mutate(Compartment2 = case_when(grepl('HP1 pseudo', Compartment2) & !grepl('DNA', Compartment2) ~ 'HP1 pseudo',
#                                         grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
#                                         T ~ Compartment2)) %>%
#         group_by(Compartment2, interval) %>%
#         summarise(MSD = mean(disp^2)/1000000, N = n()) %>%
#         mutate(Compartment2 = factor(Compartment2, levels = #c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                        c('DNA-high (subsampledx6)', 'HP1-high', 'DNA&HP1-high',
#                                          'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                          'Euchromatin'))) -> MSD_segm_pooled
#       filter(jumpDist_1234, Track_length < 15) %>%
#         mutate(Compartment2 = case_when(grepl('DNA-high', Compartment2) ~ 'DNA-high (subsampledx6)',
#                                         T ~ Compartment2)) %>% #pooling all the DNA subsampled
#         mutate(Compartment2 = case_when(grepl('HP1 pseudo', Compartment2) & !grepl('DNA', Compartment2) ~ 'HP1 pseudo',
#                                         grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
#                                         T ~ Compartment2)) %>%
#         group_by(Dish, Compartment2, interval) %>%
#         summarise(MSD = mean(disp^2)/1000000, N = n()) %>%
#         mutate(Compartment2 = factor(Compartment2, levels = #c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                        c('DNA-high (subsampledx6)', 'HP1-high', 'DNA&HP1-high',
#                                          'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                          'Euchromatin'))) -> MSD_segm_pooled_byDish
#       ggplot(MSD_segm_pooled, aes(x = interval*3, y = MSD,
#                                   colour = Compartment2, shape = Compartment2, size = Compartment2)) +
#         #alpha = N > 30, linetype = N > 30)) +
#         geom_point() +
#         geom_line(size = 0.2) +
#         # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#         #                                'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#         #                                'Euchromatin' = 'gray')) +
#         scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                        'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                        'Euchromatin' = 'gray')) +
#         scale_shape_manual(values = c(19, 19, 19, 4, 4,  4, 19)) +
#         scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2)) +
#         #scale_alpha_manual(values = c(0.5, 1)) +
#         #scale_linetype_manual(values = c('dashed', 'solid')) +
#         labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#              x = 'Interval, ms', y = 'MSD, um2') +
#         theme_bw() +
#         theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#               aspect.ratio = 1)
#       ggsave('New_analysis/MSD_HP1andDNA_poolPseudo_DNAsegm.png', width = 5.7, height = 4)
#       ggplot(MSD_segm_pooled_byDish, aes(x = interval*3, y = MSD, group = interaction(Dish, Compartment2),
#                                   colour = Compartment2, shape = Compartment2, size = Compartment2)) +
#         #alpha = N > 30, linetype = N > 30)) +
#         geom_point() +
#         geom_line(size = 0.2) +
#         # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#         #                                'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#         #                                'Euchromatin' = 'gray')) +
#         scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                        'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                        'Euchromatin' = 'gray')) +
#         scale_shape_manual(values = c(19, 19, 19, 4, 4,  4, 19)) +
#         scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2)) +
#         #scale_alpha_manual(values = c(0.5, 1)) +
#         #scale_linetype_manual(values = c('dashed', 'solid')) +
#         labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#              x = 'Interval, ms', y = 'MSD, um2') +
#         theme_bw() +
#         theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#               aspect.ratio = 1)
#       ggsave('New_analysis/MSD_HP1andDNA_poolPseudo_DNAsegm_byDish.png', width = 5.7, height = 4)
#       ggplot(MSD_segm_pooled, aes(x = log10(interval*3), y = log10(MSD),
#                                   colour = Compartment2, shape = Compartment2, size = Compartment2)) +
#         #alpha = N > 30, linetype = N > 30)) +
#         geom_point() +
#         geom_line(size = 0.2) +
#         # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#         #                                'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#         #                                'Euchromatin' = 'gray')) +
#         scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                        'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                        'Euchromatin' = 'gray')) +
#         scale_shape_manual(values = c(19, 19, 19, 4, 4,  4, 19)) +
#         scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2)) +
#         #scale_alpha_manual(values = c(0.5, 1)) +
#         #scale_linetype_manual(values = c('dashed', 'solid')) +
#         coord_cartesian(y = c(-2.15, -1.05)) +
#         labs(title = 'MSD curve', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#              x = 'Log (Interval, ms)', y = 'Log (MSD, um2)') +
#         theme_bw() +
#         theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#               aspect.ratio = 1)
#       ggsave('New_analysis/MSD_HP1andDNA_loglogLine_poolPseudo_DNAsegm.png', width = 5.7, height = 4)
#     }
#   }
#   
#   #Fit model to MSD curve, plot lines in log-log plots, plot summary
#   {
#     #All data, pool the pseudo
#     split(MSD_segm_pooled, MSD_segm_pooled$Compartment2) %>%
#       lapply(function(x) {
#         lm(log10(MSD)~log10(interval*3/1000), data = x)
#       }) -> MSDfit_segm.list
#     MSDfit_segm.df = data.frame(Compartment2 = names(MSDfit_segm.list),
#                                 Intercept = sapply(MSDfit_segm.list, function(x) {x$coefficients[[1]]}),
#                                 Slope = sapply(MSDfit_segm.list, function(x) {x$coefficients[[2]]}))
#     mutate(MSDfit_segm.df, label = paste('D =', round(10^Intercept/4, 2), 'um2/s, alpha =', round(Slope, 2), '\n')) -> MSDfit_segm.df
#     MSDfit_segm.df %>%
#       mutate(Compartment2 = factor(Compartment2, levels = #c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                      c('DNA-high (subsampledx6)', 'HP1-high', 'DNA&HP1-high',
#                                        'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                        'Euchromatin'))) -> MSDfit_segm.df
#     ggplot(MSD_segm_pooled, aes(x = log10(interval*3/1000), y = log10(MSD),
#                                 colour = Compartment2, shape = Compartment2, size = Compartment2)) +
#       #alpha = N > 30 & !(Compartment2 != 'DNA&HP1 pseudo' & N<100) &
#       #!(Compartment2 == 'DNA-high (subsampled)' & interval == 8))) +
#       geom_point() +
#       # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#       #                                'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#       #                                'Euchromatin' = 'gray'),
#       #                     labels = paste0(MSDfit_segm.df$Compartment2, '\n', MSDfit_segm.df$label)) +
#       scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray'),
#                           labels = paste0(MSDfit_segm.df$Compartment2, '\n', MSDfit_segm.df$label)) +
#       scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19),
#                          labels = paste0(MSDfit_segm.df$Compartment2, '\n', MSDfit_segm.df$label)) +
#       scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2),
#                         labels = paste0(MSDfit_segm.df$Compartment2, '\n', MSDfit_segm.df$label)) +
#       #scale_alpha_manual(values = c(0.3, 1), guide = 'none') +
#       geom_abline(data = MSDfit_segm.df, aes(intercept = Intercept, slope = Slope, colour = Compartment2)) +
#       labs(title = 'MSD fit', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#            x = 'Log(Interval, s)', y = 'Log(MSD, um2)') +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 1)
#     ggsave('New_analysis/MSD_HP1andDNA_loglogFit_subDNA.png', width = 5.5, height = 4)
#     
#     ggplot(MSDfit_segm.df, aes(x = Compartment2, y = Intercept-log10(4),
#                                colour = Compartment2, shape = Compartment2)) +
#       geom_point() +
#       scale_y_continuous(labels = math_format(10^.x),
#                          sec.axis = sec_axis(trans=~10^., breaks = c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2))) +
#       scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray')) +
#       # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#       #                                'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#       #                                'Euchromatin' = 'gray')) +
#       scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#       labs(x = '', y = 'D, um2/s') +
#       theme_bw() +
#       theme(aspect.ratio = 0.8, legend.position = 'none', axis.text.x = element_text(angle = 30, hjust = 1)) -> Ds
#     ggplot(MSDfit_segm.df, aes(x = Compartment2, y = Slope, colour = Compartment2, shape = Compartment2)) +
#       geom_point() +
#       scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray')) +
#       # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#       #                                'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#       #                                'Euchromatin' = 'gray')) +
#       scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#       scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
#       labs(x = '', y = 'Alpha') +
#       theme_bw() +
#       theme(aspect.ratio = 0.8, legend.position = 'none', axis.text.x = element_text(angle = 30, hjust = 1)) -> alphas
#     ggarrange(Ds, alphas, nrow = 2, ncol = 1)
#     ggsave('New_analysis/D_alpha_HP1andDNA_subDNA.png', width = 4.5, height = 5.5)
#   }
#   
#   #Same, but only using intervals 1-5
#   {
#     #All data, pool the pseudo
#     split(filter(MSD_segm_pooled, interval < 6), filter(MSD_segm_pooled, interval < 6)$Compartment2) %>%
#       lapply(function(x) {
#         lm(log10(MSD)~log10(interval*3/1000), data = x)
#       }) -> MSDfit_segm.list2
#     MSDfit_segm.df2 = data.frame(Compartment2 = names(MSDfit_segm.list2),
#                                  Intercept = sapply(MSDfit_segm.list2, function(x) {x$coefficients[[1]]}),
#                                  Slope = sapply(MSDfit_segm.list2, function(x) {x$coefficients[[2]]}))
#     mutate(MSDfit_segm.df2, label = paste('D =', round(10^Intercept/4, 2), 'um2/s, alpha =', round(Slope, 2), '\n')) -> MSDfit_segm.df2
#     MSDfit_segm.df2 %>%
#       mutate(Compartment2 = factor(Compartment2, levels = #c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                      c('DNA-high (subsampledx6)', 'HP1-high', 'DNA&HP1-high',
#                                        'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                        'Euchromatin'))) -> MSDfit_segm.df2
#     ggplot(MSD_segm_pooled, aes(x = log10(interval*3/1000), y = log10(MSD),
#                                 colour = Compartment2, shape = Compartment2, size = Compartment2,
#                                 alpha = interval < 6)) +
#       geom_point() +
#       # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#       #                                'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#       #                                'Euchromatin' = 'gray'),
#       scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray'),
#                           labels = paste0(MSDfit_segm.df2$Compartment2, '\n', MSDfit_segm.df2$label)) +
#       scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19),
#                          labels = paste0(MSDfit_segm.df2$Compartment2, '\n', MSDfit_segm.df2$label)) +
#       scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2),
#                         labels = paste0(MSDfit_segm.df2$Compartment2, '\n', MSDfit_segm.df2$label)) +
#       scale_alpha_manual(values = c(0.3, 1), guide = 'none') +
#       geom_abline(data = MSDfit_segm.df2, aes(intercept = Intercept, slope = Slope, colour = Compartment2)) +
#       labs(title = 'MSD fit', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#            x = 'Log(Interval, s)', y = 'Log(MSD, um2)') +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 1)
#     ggsave('New_analysis/MSD_HP1andDNA_loglogFit2_subDNA.png', width = 5.5, height = 4)
#     
#     ggplot(MSDfit_segm.df2, aes(x = Compartment2, y = Intercept-log10(4), colour = Compartment2, shape = Compartment2)) +
#       geom_point() +
#       scale_y_continuous(labels = math_format(10^.x),
#                          sec.axis = sec_axis(trans=~10^., breaks = c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2))) +
#       # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#       #                                'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#       #                                'Euchromatin' = 'gray')) +
#       scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray')) +
#       scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#       labs(x = '', y = 'D, um2/s') +
#       theme_bw() +
#       theme(aspect.ratio = 0.8, legend.position = 'none', axis.text.x = element_text(angle = 30, hjust = 1)) -> Ds
#     ggplot(MSDfit_segm.df2, aes(x = Compartment2, y = Slope, colour = Compartment2, shape = Compartment2)) +
#       geom_point() +
#       # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#       #                                'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#       #                                'Euchromatin' = 'gray')) +
#       scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray')) +
#       scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#       scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
#       labs(x = '', y = 'Alpha') +
#       theme_bw() +
#       theme(aspect.ratio = 0.8, legend.position = 'none', axis.text.x = element_text(angle = 30, hjust = 1)) -> alphas
#     ggarrange(Ds, alphas, nrow = 2, ncol = 1)
#     ggsave('New_analysis/D_alpha_HP1andDNA2_subDNA.png', width = 4.5, height = 5.5)
#   }
#   
#   #Fit model to bootstrapped MSD curves (using intervals 1-5),
#   #plot lines in log-log plots and summary with error bars
#   {
#     #bootstrap
#     MSD_segm_pooled_bs = list()
#     tic()
#     lapply(1:100, function(n) {
#       #calc time weird, now much faster than before
#       #100 with lapply - 150s
#       filter(jumpDist_1234, Track_length < 15) %>%
#         mutate(Compartment2 = case_when(grepl('DNA-high', Compartment2) ~ 'DNA-high (subsampledx6)',
#                                         T ~ Compartment2)) %>% #pooling all the DNA subsampled
#         mutate(Compartment2 = case_when(grepl('HP1 pseudo', Compartment2) & !grepl('DNA', Compartment2) ~ 'HP1 pseudo',
#                                         grepl('DNA&HP1 pseudo', Compartment2) ~ 'DNA&HP1 pseudo',
#                                         T ~ Compartment2)) %>%
#         group_by(Compartment2, interval) %>% slice_sample(prop = 1, replace = T) %>%
#         summarise(MSD = mean(disp^2)/1000000, N = n()) -> MSD.new
#       return(MSD.new)
#     }) -> MSD_segm_pooled_bs
#     toc()
#     
#     bind_rows(MSD_segm_pooled_bs) %>% group_by(Compartment2, interval) %>%
#       summarise(mean = mean(MSD), sd = sd(MSD),
#                 CI_l = quantile(MSD, 0.05), CI_h = quantile(MSD, 0.95)) -> MSD_segm_bs_stats
#     
#     
#     #fit
#     tic()
#     lapply(MSD_segm_pooled_bs, function(s) {
#       split(filter(s, interval < 6), filter(s, interval < 6)$Compartment2) %>%
#         lapply(function(x) {
#           lm(log10(MSD)~log10(interval*3/1000), data = x)
#         }) -> fit.list
#       fit.df = data.frame(Compartment2 = names(fit.list),
#                           Intercept = sapply(fit.list, function(x) {x$coefficients[[1]]}),
#                           Slope = sapply(fit.list, function(x) {x$coefficients[[2]]}))
#       return(fit.df)
#     }) -> MSDfit_segm_bs.df
#     toc()
#     
#     bind_rows(MSDfit_segm_bs.df) %>% group_by(Compartment2) %>%
#       summarise(int.mean = mean(Intercept), int.sd = sd(Intercept),
#                 int.CI_l = quantile(Intercept, 0.05), int.CI_h = quantile(Intercept, 0.95),
#                 slope.mean = mean(Slope), slope.sd = sd(Slope),
#                 slope.CI_l = quantile(Slope, 0.05), slope.CI_h = quantile(Slope, 0.95)) -> MSDfit_segm_bs_stats
#     
#     ggplot(MSD_segm_pooled, aes(x = log10(interval*3/1000), y = log10(MSD),
#                                 colour = Compartment2, shape = Compartment2, size = Compartment2,
#                                 alpha = interval < 6)) +
#       geom_point() +
#       # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#       #                                'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#       #                                'Euchromatin' = 'gray'),
#       #                     labels = paste0(MSDfit_segm.df2$Compartment2, '\n', MSDfit_segm.df2$label)) +
#       scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray'),
#                           labels = paste0(MSDfit_segm.df2$Compartment2, '\n', MSDfit_segm.df2$label)) +
#       scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19),
#                          labels = paste0(MSDfit_segm.df2$Compartment2, '\n', MSDfit_segm.df2$label)) +
#       scale_size_manual(values = c(2, 2, 2, 1, 1, 1, 2),
#                         labels = paste0(MSDfit_segm.df2$Compartment2, '\n', MSDfit_segm.df2$label)) +
#       scale_alpha_manual(values = c(0.3, 1), guide = 'none') +
#       geom_errorbar(data = MSD_segm_bs_stats, aes(ymin = log10(CI_l), ymax = log10(CI_h),
#                                                   y = NULL),
#                     width = 0.01, size = 0.5) +
#       geom_abline(data = MSDfit_segm.df2, aes(intercept = Intercept, slope = Slope, colour = Compartment2)) +
#       labs(title = 'MSD fit', subtitle = 'Dens < 10 loc/fr, traj <15\nHP1 segm with DBSCAN',
#            x = 'Log(Interval, s)', y = 'Log(MSD, um2)') +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 1)
#     ggsave('New_analysis/MSD_HP1andDNA_loglogFit2_subDNA_errors.png', width = 5.5, height = 4)
#     
#     ggplot(MSDfit_segm.df2, aes(x = Compartment2, y = Intercept - log10(4), colour = Compartment2, shape = Compartment2)) +
#       geom_point() +
#       scale_y_continuous(labels = math_format(10^.x),
#                          sec.axis = sec_axis(trans=~10^., breaks = c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2))) +
#       # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#       #                                'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#       #                                'Euchromatin' = 'gray')) +
#       scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray')) +
#       scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#       geom_errorbar(data = MSDfit_segm_bs_stats, aes(ymin = int.CI_l-log10(4), ymax = int.CI_h-log10(4),
#                                                      y = NULL, alpha = NULL),
#                     width = 0.1, size = 0.3) +
#       labs(x = '', y = 'D, um2/s') +
#       theme_bw() +
#       theme(aspect.ratio = 0.8, legend.position = 'none', axis.text.x = element_text(angle = 30, hjust = 1)) -> Ds
#     
#     ggplot(MSDfit_segm.df2, aes(x = Compartment2, y = Slope, colour = Compartment2, shape = Compartment2)) +
#       geom_point() +
#       # scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#       #                                'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#       #                                'Euchromatin' = 'gray')) +
#       scale_colour_manual(values = c('DNA-high (subsampledx6)' = 'steelblue', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray')) +
#       scale_shape_manual(values = c(19, 19, 19, 4, 4, 4, 19)) +
#       geom_errorbar(data = MSDfit_segm_bs_stats, aes(ymin = slope.CI_l, ymax = slope.CI_h,
#                                                      y = NULL, alpha = NULL),
#                     width = 0.1, size = 0.3) +
#       scale_y_continuous(lim = c(-0.02, 1), breaks = seq(0, 1, by = 0.2)) +
#       labs(x = '', y = 'Alpha') +
#       theme_bw() +
#       theme(aspect.ratio = 0.8, legend.position = 'none', axis.text.x = element_text(angle = 30, hjust = 1)) -> alphas
#     ggarrange(Ds, alphas, nrow = 2, ncol = 1)
#     ggsave('New_analysis/D_alpha_HP1andDNA_subDNA_errors.png', width = 4.5, height = 5.5)
#   }
# }
# 
# ##Autocovariance
# {
#   #S_n = sum_xy(mean_track(mean(dx_j*dx_(j+n))))
#   #Can use jumpDist_1234
#   {
#   filter(jumpDist_1234, interval == 1) %>% ungroup() %>%
#     mutate(Compartment = case_when(DNA_HP1_foci ~ 'DNA&HP1-high',
#                                    DNA_focus_idx!=0 & !DNA_HP1_foci ~ 'DNA-high',
#                                    # Pseudo_in_DNA1 & !DNA_HP1_foci ~ 'DNA-high (subsampled 1)',
#                                    # Pseudo_in_DNA2 & !DNA_HP1_foci ~ 'DNA-high (subsampled 2)',
#                                    # Pseudo_in_DNA3 & !DNA_HP1_foci ~ 'DNA-high (subsampled 3)',
#                                    # Pseudo_in_DNA4 & !DNA_HP1_foci ~ 'DNA-high (subsampled 4)',
#                                    # Pseudo_in_DNA5 & !DNA_HP1_foci ~ 'DNA-high (subsampled 5)',
#                                    # Pseudo_in_DNA6 & !DNA_HP1_foci ~ 'DNA-high (subsampled 6)',
#                                    !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
#                                    DNA_pseudo & !(HP1_pseudo_dbscan1 | HP1_pseudo_dbscan2 |
#                                                     HP1_pseudo_dbscan3 | HP1_pseudo_dbscan4 |
#                                                     HP1_pseudo_dbscan5 | HP1_pseudo_dbscan6) ~ 'DNA pseudo',
#                                    DNA_HP1_pseudo1 ~ 'DNA&HP1 pseudo 1',
#                                    DNA_HP1_pseudo2 ~ 'DNA&HP1 pseudo 2',
#                                    DNA_HP1_pseudo3 ~ 'DNA&HP1 pseudo 3',
#                                    DNA_HP1_pseudo4 ~ 'DNA&HP1 pseudo 4',
#                                    DNA_HP1_pseudo5 ~ 'DNA&HP1 pseudo 5',
#                                    DNA_HP1_pseudo6 ~ 'DNA&HP1 pseudo 6',
#                                    !DNA_HP1_pseudo1 & HP1_pseudo_dbscan1 ~ 'HP1 pseudo 1',#NB! These overlap -> pseudo6 ends up smaller
#                                    !DNA_HP1_pseudo2 & HP1_pseudo_dbscan2 ~ 'HP1 pseudo 2',
#                                    !DNA_HP1_pseudo3 & HP1_pseudo_dbscan3 ~ 'HP1 pseudo 3',
#                                    !DNA_HP1_pseudo4 & HP1_pseudo_dbscan4 ~ 'HP1 pseudo 4',
#                                    !DNA_HP1_pseudo5 & HP1_pseudo_dbscan5 ~ 'HP1 pseudo 5',
#                                    !DNA_HP1_pseudo6 & HP1_pseudo_dbscan6 ~ 'HP1 pseudo 6',
#                                    T ~ 'Euchromatin')) %>%
#     mutate(Compartment2 = case_when(grepl('HP1 pseudo', Compartment) & !grepl('DNA', Compartment) ~ 'HP1 pseudo',
#                                     grepl('DNA&HP1 pseudo', Compartment) ~ 'DNA&HP1 pseudo',
#                                     T ~ Compartment)) -> jumpDist_1
#   }
#   #Pseudo separately
#   for (lag in 0:6) {
#     jumpDist_1 %>% group_by(Dish, Pos, Track) %>%
#       mutate(corr_x = diff.x*lag(diff.x, lag),
#              corr_y = diff.y*lag(diff.y, lag)) %>%
#       #Compartment1 = Compartment,
#       #Compartment2 = lag(Compartment, lag)) %>%
#       drop_na() %>% group_by(Dish, Pos, Track, Compartment) %>%
#       summarise(trcorr_x = mean(corr_x, na.rm = T),
#                 trcorr_y = mean(corr_y, na.rm = T)) %>%
#       mutate(trcorr = trcorr_x + trcorr_y) %>%
#       ungroup() %>% group_by(Compartment) %>%
#       summarise(avcorr_x = mean(trcorr_x), avcorr_y = mean(trcorr_y), avcorr = mean(trcorr)) -> tmp
#     if (lag == 0) {
#       avcorr = mutate(tmp, lag = lag)
#     } else {
#       avcorr = rbind(avcorr, mutate(tmp, lag = lag))
#     }
#     # avcorr$avcorr_x[lag+1] = tmp$avcorr_x
#     # avcorr$avcorr_y[lag+1] = tmp$avcorr_y
#     # avcorr$avcorr[lag+1] = tmp$avcorr
#   }
#   {
#   avcorr %>%
#     mutate(Compartment2 = case_when(grepl('DNA&HP1 pseudo', Compartment) ~ 'DNA&HP1 pseudo',
#                                     grepl('HP1 pseudo', Compartment) ~ 'HP1 pseudo',
#                                     T ~ Compartment)) %>%
#     mutate(Compartment2 = factor(Compartment2, levels = c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                                           'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                                           'Euchromatin'))) %>%
#     ggplot(aes(x = lag, y = avcorr/1000000, colour = Compartment2, group = Compartment)) +
#     geom_point() +
#     geom_line() +
#     #geom_line(aes(y = avcorr_x/1000000), linetype = 'longdash', size = 0.2) +
#     #geom_line(aes(y = avcorr_y/1000000), linetype = 'longdash', size = 0.2) +
#     geom_hline(yintercept = 0, colour = 'black') +
#     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                    'Euchromatin' = 'gray')) +
#     scale_x_continuous(breaks = 0:6) +
#     labs(x = 'Lag, frames', y = 'Mean covariance, um2',
#          title = 'Autocovariance', colour = 'Compartment') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1)
#   ggsave('New_analysis/Autocovariance_segm.png', width = 4.5, height = 3)
#   avcorr %>% filter(Compartment == 'DNA-high' | Compartment == 'HP1-high' |
#                       Compartment == 'DNA&HP1-high' | Compartment == 'Euchromatin') %>%
#     ggplot(aes(x = lag, y = avcorr/1000000, colour = Compartment)) +
#     geom_point() +
#     geom_line() +
#     # geom_line(aes(y = avcorr_x/1000000), linetype = 'longdash', size = 0.2) +
#     # geom_line(aes(y = avcorr_y/1000000), linetype = 'longdash', size = 0.2) +
#     geom_hline(yintercept = 0, colour = 'black') +
#     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                    'Euchromatin' = 'gray')) +
#     scale_x_continuous(breaks = 0:6) +
#     labs(x = 'Lag, frames', y = 'Mean covariance, um2',
#          title = 'Autocovariance', colour = 'Compartment', subtitle = 'Solid = 2D, dashed = 1D') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1)
#   ggsave('New_analysis/Autocovariance_segm_only.png', width = 4.5, height = 3)
#   }
#   
#   #Pseudo pooled
#   for (lag in 0:6) {
#     jumpDist_1 %>% group_by(Dish, Pos, Track) %>%
#       mutate(corr_x = diff.x*lag(diff.x, lag),
#              corr_y = diff.y*lag(diff.y, lag)) %>%
#       #Compartment1 = Compartment,
#       #Compartment2 = lag(Compartment, lag)) %>%
#       drop_na() %>% group_by(Dish, Pos, Track, Compartment2) %>%
#       summarise(trcorr_x = mean(corr_x, na.rm = T),
#                 trcorr_y = mean(corr_y, na.rm = T)) %>%
#       ungroup() %>% group_by(Compartment2) %>%
#       summarise(avcorr_x = mean(trcorr_x), avcorr_y = mean(trcorr_y)) %>%
#       mutate(avcorr = avcorr_x + avcorr_y) -> tmp
#     if (lag == 0) {
#       avcorr2 = mutate(tmp, lag = lag)
#     } else {
#       avcorr2 = rbind(avcorr2, mutate(tmp, lag = lag))
#     }
#     # avcorr$avcorr_x[lag+1] = tmp$avcorr_x
#     # avcorr$avcorr_y[lag+1] = tmp$avcorr_y
#     # avcorr$avcorr[lag+1] = tmp$avcorr
#   }
#   {
#   avcorr2 %>%
#     mutate(Compartment2 = factor(Compartment2, levels = c('DNA-high', 'HP1-high', 'DNA&HP1-high',
#                                                           'DNA pseudo', 'HP1 pseudo', 'DNA&HP1 pseudo',
#                                                           'Euchromatin'))) %>%
#     ggplot(aes(x = lag, y = avcorr/1000000, colour = Compartment2)) +
#     geom_point() +
#     geom_line() +
#     #geom_line(aes(y = avcorr_x/1000000), linetype = 'longdash', size = 0.2) +
#     #geom_line(aes(y = avcorr_y/1000000), linetype = 'longdash', size = 0.2) +
#     geom_hline(yintercept = 0, colour = 'black') +
#     scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                    'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                    'Euchromatin' = 'gray')) +
#     scale_x_continuous(breaks = 0:6) +
#     labs(x = 'Lag, frames', y = 'Mean covariance, um2',
#          title = 'Autocovariance') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 1)
#   ggsave('New_analysis/Autocovariance_segm_pooledPseudo.png', width = 4.5, height = 3)
#   }
#   
#   #Pseudo pooled, ensemble-averaged correlation
#   for (lag in 0:6) {
#     jumpDist_1 %>% group_by(Dish, Pos, Track) %>%
#       mutate(corr_x = diff.x*lag(diff.x, lag),
#              corr_y = diff.y*lag(diff.y, lag)) %>%
#       #Compartment1 = Compartment,
#       #Compartment2 = lag(Compartment, lag)) %>%
#       drop_na() %>% ungroup() %>% group_by(Compartment2) %>%
#       summarise(avcorr_x = mean(corr_x), avcorr_y = mean(corr_y)) %>%
#       mutate(avcorr = avcorr_x + avcorr_y) -> tmp
#     if (lag == 0) {
#       avcorr3 = mutate(tmp, lag = lag)
#     } else {
#       avcorr3 = rbind(avcorr3, mutate(tmp, lag = lag))
#     }
#     # avcorr$avcorr_x[lag+1] = tmp$avcorr_x
#     # avcorr$avcorr_y[lag+1] = tmp$avcorr_y
#     # avcorr$avcorr[lag+1] = tmp$avcorr
#   }
#   {
#     avcorr3 %>% filter(Compartment2 == 'DNA-high' | Compartment2 == 'HP1-high' |
#                         Compartment2 == 'DNA&HP1-high' | Compartment2 == 'Euchromatin') %>%
#       ggplot(aes(x = lag, y = avcorr/1000000, colour = Compartment2)) +
#       geom_point() +
#       geom_line() +
#       # geom_line(aes(y = avcorr_x/1000000), linetype = 'longdash', size = 0.2) +
#       # geom_line(aes(y = avcorr_y/1000000), linetype = 'longdash', size = 0.2) +
#       geom_hline(yintercept = 0, colour = 'black') +
#       scale_colour_manual(values = c('DNA-high' = 'cyan', 'HP1-high' = 'red', 'DNA&HP1-high' = 'purple',
#                                      'DNA pseudo' = 'blue', 'HP1 pseudo' = 'tomato4', 'DNA&HP1 pseudo' = 'purple4',
#                                      'Euchromatin' = 'gray')) +
#       scale_x_continuous(breaks = 0:6) +
#       labs(x = 'Lag, frames', y = 'Mean covariance, um2',
#            title = 'Autocovariance', colour = 'Compartment', subtitle = 'Solid = 2D, dashed = 1D') +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#             aspect.ratio = 1)
#   }
#   
#   #Calculate D and alpha
#   {
#     solve_bisection <- function(f, a, b, n = 10000, tol = 1e-9) {
#       if (((f(a) > 0) && (f(b) > 0)) | ((f(a) < 0) && (f(b) < 0))) {
#         stop('signs of f(a) and f(b) do not differ')
#       }
#       
#       for (i in 1:n) {
#         c <- (a + b) / 2 
#         
#         if ((f(c) == 0) || ((b - a) / 2) < tol) {
#           return(c)
#         }
#         
#         ifelse(sign(f(c)) == sign(f(a)), 
#                a <- c,
#                b <- c)
#       }
#       print('Too many iterations')
#     }
#     gen_obj_fn <- function(obj_fn, other_pars1) {
#       function(x) { obj_fn( x, other_pars1) }
#     }
#     S0S1S2_alpha = function(alpha, A) {
#       2*(2^alpha-1)/(3^alpha+1-2^(alpha+1))-A
#     }
#     
#     dt = 0.003
#     S0 = filter(avcorr3, Compartment2 == 'DNA-high', lag < 3)$avcorr[1]/1000000
#     S1 = filter(avcorr3, Compartment2 == 'DNA-high', lag < 3)$avcorr[2]/1000000
#     S2 = filter(avcorr3, Compartment2 == 'DNA-high', lag < 3)$avcorr[3]/1000000
#     
#     alpha = solve_bisection(gen_obj_fn(S0S1S2_alpha, (S0+2*S1)/S2), 0.000001, 1-0.000001)
#     D = (S0+2*S1)/(4*dt^alpha*(2^alpha-1))
#     sigma = sqrt((S0-4*D*dt^alpha)/4)
#   }
#   #For some reason the D and alpha calculated this way come out much higher than those from MSD analysis
#   #They still match the MSD1 value, but overshoot a lot for subsequent MSD values
#   #In contrast, the MSD-derived values predict S1 and S2 that are smaller than the real ones (i.e. larger in magnitude)
#   #The same happens even if I calculate ensemble-averaged correlation rather than averaging twice (per traj and across traj):
#   #the autocov function comes out very similar
# }
}
#

#### Write data for SA-SPT ####
#SA-SPT uses only individual jumps and correlations between them over a lag
#=> classify each *JUMP* by compartment (if needed, duplicate the localisation);
#also, if there is a compartment transition within a trajectory, split this trajectory

###Read in the data and the segmentation
foci_list = readRDS('DNA_foci_list.rds')
pfoci_list = readRDS('DNA_pseudofoci_list.rds')
segmCoresList = readRDS('HP1_dbscan.RDS')
psegmCoresListInDNA = readRDS('HP1_dbscan_pseudo_in_DNA.RDS')

data_list = readRDS('Data_list_SPT_segm.RDS')

thresholds2 = c(19500, 11000, 19000, 20000, 20000, 10000, #Dish1 (0, 10, 11, 1, 2, 3)
                10000, 21000, 12000, 20000,11000, 21000, #Dish1 (4, 5, 6, 7, 8, 9)
                25000, 24000, 24000, 20000, 23000, 10000, 33000, 17000, 14000) #Dish2 (1-8)WFframes = 8WFframes = 8

###All data
{
lapply(seq_along(data_list), function(c) {
  pos = names(data_list)[c]
  print(pos)
  select(filter(data_list[[pos]], Frame > thresholds2[c]),
         Dish, Pos, Track, Frame, x, y, Track_length) %>%
    filter(Track_length > 1) %>% group_by(Track) %>% mutate(Track = cur_group_id()) -> df
  return(df)
}) -> for_SASPT_all
numframes = 40000

maxtraj = sapply(for_SASPT_all, function(Pos) max(Pos$Track))
lapply(seq_along(for_SASPT_all), function(i) {
  for_SASPT_all[[i]] %>% select(Track, x, y, Frame) %>% mutate(x = x/100, y = y/100) %>%
    rename(trajectory = Track, frame = Frame) %>%
    mutate(trajectory = trajectory + sum(maxtraj[1:i-1]), frame = frame + numframes*(i-1)) -> result
  return(result)
}) %>% bind_rows() %>%# group_by(trajectory) %>% mutate(JD = sqrt((x-lag(x))^2 + (y-lag(y))^2)) %>% View()
  write.table(paste0('New_analysis/SA-SPT/Data/AllData.csv'),
              col.names = T, row.names = F, quote = F, sep = ',')
}

###Classified data
##Classify (jumps)
{
WFaveraging = 5000 #WFframes*WFaveraging = total number of frames
pixel = 140
subsampling = 10
n_shuffle = 6
lapply(seq_along(data_list), function(c) {
  tic()
  pos = names(data_list)[c]
  print(pos)
  df = select(filter(data_list[[pos]], Frame > thresholds2[c]),
              Dish, Pos, Track, Frame, x, y, Track_length)
  
  #Calculate JD, classify by DNA
  dim_x = foci_list[[pos]]$dim_x[1]
  dim_y = foci_list[[pos]]$dim_y[1]
  df %>% group_by(Dish, Pos, Track) %>% filter(Track_length > 1) %>%
    mutate(x_jump = (x+lag(x))/2, y_jump = (y+lag(y))/2, Frame_jump = floor((Frame+lag(Frame))/2)) %>%
    mutate(px_id = floor(x_jump*subsampling/pixel)*dim_y + ceiling(y_jump*subsampling/pixel),
           frame = ceiling(Frame_jump/WFaveraging)) %>%
    left_join(., foci_list[[pos]], by = c('px_id', 'frame')) %>%
    select(-x_left, -y_top, -dim_y, -dim_x, -pos, -dish) %>%
    left_join(., pfoci_list[[pos]], by = c('px_id', 'frame')) %>%
    select(-x_left, -y_top, -dim_y, -dim_x, -pos, -dish, -frame, -px_id) %>%
    rename(DNA_focus_idx = focus_id, DNA_pseudo = pseudofocus) %>%
    mutate(DNA_focus_idx = replace_na(DNA_focus_idx, 0),
           DNA_pseudo = replace_na(DNA_pseudo, 0)) -> result
  
  # #DNA foci, subsampled in HP1-foci-like manner
  # for (n in 1:n_shuffle) {
  #   pseudosegmentation = list('cluster' = psegmCoresListInDNA[[pos]][[paste0('Pseudo', n)]]$cluster,
  #                             'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
  #   class(pseudosegmentation) = c('dbscan_fast', 'dbscan')
  #   result[[paste0('Pseudo_in_DNA', n)]] =
  #     predict(object = pseudosegmentation,
  #             newdata = select(ungroup(result), x_jump, y_jump),
  #             data = select(ungroup(psegmCoresListInDNA[[pos]][[paste0('Pseudo', n)]]), x, y))
  #   #remove parts of clusters that don't overlap with DNA foci
  #   result %>% mutate(!!paste0('Pseudo_in_DNA', n) :=
  #                       case_when(DNA_focus_idx == 0 ~ as.integer(0),
  #                                 T ~ !!as.name(paste0('Pseudo_in_DNA', n)))) -> result
  # }
  
  #Classify by HP1
  segmentation = list('cluster' = segmCoresList[[pos]][['Foci']]$cluster,
                      'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
  class(segmentation) = c('dbscan_fast', 'dbscan')
  result$HP1_foci_dbscan = predict(object = segmentation,
                                   newdata = select(ungroup(result), x_jump, y_jump),
                                   data = select(ungroup(segmCoresList[[pos]][['Foci']]), x, y))
  for (n in 1:n_shuffle) {
    pseudosegmentation = list('cluster' = segmCoresList[[pos]][[paste0('Pseudo', n)]]$cluster,
                              'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
    class(pseudosegmentation) = c('dbscan_fast', 'dbscan')
    result[[paste0('HP1_pseudo_dbscan', n)]] =
      predict(object = pseudosegmentation,
              newdata = select(ungroup(result), x_jump, y_jump),
              data = select(ungroup(segmCoresList[[pos]][[paste0('Pseudo', n)]]), x, y))
  }
  
  #Overlap
  result$DNA_HP1_foci = result$HP1_foci_dbscan %in% DNA_HP1_foci[[pos]]
  for (n in 1:n_shuffle) {
    result[[paste0('DNA_HP1_pseudo', n)]] = result[[paste0('HP1_pseudo_dbscan', n)]] %in% DNA_HP1_foci[[pos]]
  }
  
  #Determine the overall compartment, split the data frame
  result[is.na(result$x_jump), 10:31] = NA #for first locs in a traj
  result %>%
    mutate(Compartment2 = case_when(is.na(DNA_HP1_foci) ~ 'Start',
                                    DNA_HP1_foci ~ 'DNA&HP1-high',
                                    DNA_focus_idx!=0 & !DNA_HP1_foci ~ 'DNA-high', #only DNA-high
                                    !DNA_HP1_foci & HP1_foci_dbscan!=0 ~ 'HP1-high',
                                    DNA_pseudo != 0 ~ 'DNA pseudo', #this excludes HP1-high
                                    (DNA_HP1_pseudo1 | DNA_HP1_pseudo2 | DNA_HP1_pseudo3 |
                                       DNA_HP1_pseudo4 | DNA_HP1_pseudo5 | DNA_HP1_pseudo6) ~ 'DNA&HP1 pseudo',
                                    ((!DNA_HP1_pseudo1 & HP1_pseudo_dbscan1) | (!DNA_HP1_pseudo2 & HP1_pseudo_dbscan2) | 
                                       (!DNA_HP1_pseudo3 & HP1_pseudo_dbscan3) | (!DNA_HP1_pseudo4 & HP1_pseudo_dbscan4) | 
                                       (!DNA_HP1_pseudo5 & HP1_pseudo_dbscan5) | (!DNA_HP1_pseudo6 & HP1_pseudo_dbscan6)) ~ 'HP1 pseudo',
                                    #HP1 pseudo landed within DNA foci are excluded
                                    T ~ 'Euchromatin')) %>% #group_by(Compartment) %>% summarise(n()) %>% View()
    group_by(Dish, Pos, Track) %>% mutate(Compartment1 = lead(Compartment2)) %>%
    mutate(Compartment1 = case_when(is.na(Compartment1) ~ 'End',
                                    T ~ Compartment1)) %>%
    select(Dish, Pos, Track, Frame, x, y, Track_length, Compartment1, Compartment2) -> result
  result %>% mutate(Compartment = case_when(Compartment1 == 'End' ~ Compartment2,
                                            T ~ Compartment1)) %>%
    rbind(filter(result, Compartment1 != Compartment2 & Compartment1 != 'End' & Compartment2 != 'Start') %>%
            mutate(Compartment = Compartment2)) %>% filter(Compartment != 'Start') %>%
    select(-Compartment1, -Compartment2) %>% #group_by(Track, Compartment) %>% arrange(Frame) %>%
    #mutate(gap = Frame - lag(Frame)) %>% ungroup() %>% group_by(gap) %>% summarise(n()) %>% View()
    split(.$Compartment) -> result_list
  lapply(result_list, function(df) {
    df %>% group_by(Track) %>% arrange(Frame) %>% mutate(gap = Frame - lag(Frame)) %>%
      ungroup() %>% arrange(Track) -> df
    df$Track -> track_ids
    for(i in 1:nrow(df)) {
      if (!is.na(df$gap[i]) & df$gap[i] > 1) {
        track_ids[i:length(track_ids)] = track_ids[i:length(track_ids)]+1
      }
    }
    df$Track = track_ids
    df %>% select(-gap) %>% group_by(Track) %>% mutate(Track = cur_group_id()) -> df
  }) -> result_list
  return(result_list)
  print(toc())
}) -> for_SASPT
names(for_SASPT) = names(data_list)
#Specially for subsampled DNA foci
lapply(seq_along(data_list), function(c) {
  tic()
  pos = names(data_list)[c]
  print(pos)
  df = select(filter(data_list[[pos]], Frame > thresholds2[c]),
              Dish, Pos, Track, Frame, x, y, Track_length)
  
  #Calculate mean coord, classify by DNA
  dim_x = foci_list[[pos]]$dim_x[1]
  dim_y = foci_list[[pos]]$dim_y[1]
  df %>% group_by(Dish, Pos, Track) %>% filter(Track_length > 1) %>%
    mutate(x_jump = (x+lag(x))/2, y_jump = (y+lag(y))/2, Frame_jump = floor((Frame+lag(Frame))/2)) %>%
    mutate(px_id = floor(x_jump*subsampling/pixel)*dim_y + ceiling(y_jump*subsampling/pixel),
           frame = ceiling(Frame_jump/WFaveraging)) %>%
    left_join(., foci_list[[pos]], by = c('px_id', 'frame')) %>%
    select(-x_left, -y_top, -dim_y, -dim_x, -pos, -dish) %>%
    left_join(., pfoci_list[[pos]], by = c('px_id', 'frame')) %>%
    select(-x_left, -y_top, -dim_y, -dim_x, -pos, -dish, -frame, -px_id) %>%
    rename(DNA_focus_idx = focus_id, DNA_pseudo = pseudofocus) %>%
    mutate(DNA_focus_idx = replace_na(DNA_focus_idx, 0),
           DNA_pseudo = replace_na(DNA_pseudo, 0)) -> result
  
  #DNA foci, subsampled in HP1-foci-like manner
  for (n in 1:n_shuffle) {
    pseudosegmentation = list('cluster' = psegmCoresListInDNA[[pos]][[paste0('Pseudo', n)]]$cluster,
                              'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
    class(pseudosegmentation) = c('dbscan_fast', 'dbscan')
    result[[paste0('Pseudo_in_DNA', n)]] =
      predict(object = pseudosegmentation,
              newdata = select(ungroup(result), x_jump, y_jump),
              data = select(ungroup(psegmCoresListInDNA[[pos]][[paste0('Pseudo', n)]]), x, y))
    #remove parts of clusters that don't overlap with DNA foci
    result %>% mutate(!!paste0('Pseudo_in_DNA', n) :=
                        case_when(DNA_focus_idx == 0 ~ as.integer(0),
                                  T ~ !!as.name(paste0('Pseudo_in_DNA', n)))) -> result
  }
  
  #Overlapping DNA&HP1 foci (to exclude these from subsampled DNA)
  segmentation = list('cluster' = segmCoresList[[pos]][['Foci']]$cluster,
                      'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
  class(segmentation) = c('dbscan_fast', 'dbscan')
  result$HP1_foci_dbscan = predict(object = segmentation,
                                   newdata = select(ungroup(result), x_jump, y_jump),
                                   data = select(ungroup(segmCoresList[[pos]][['Foci']]), x, y))
  result$DNA_HP1_foci = result$HP1_foci_dbscan %in% DNA_HP1_foci[[pos]]
  
  #Determine the overall compartment, split the data frame
  result[is.na(result$x_jump), 10:19] = NA #for first locs in a traj
  result %>% group_by(Dish, Pos, Track) %>%
    mutate(Compartment = case_when((Pseudo_in_DNA1 | Pseudo_in_DNA2 | Pseudo_in_DNA3 |
                                      Pseudo_in_DNA4 | Pseudo_in_DNA5 | Pseudo_in_DNA6) & !DNA_HP1_foci ~ 'DNA-high (subsampled)',
                                   (lead(Pseudo_in_DNA1) | lead(Pseudo_in_DNA2) | lead(Pseudo_in_DNA3) |
                                      lead(Pseudo_in_DNA4) | lead(Pseudo_in_DNA5) | lead(Pseudo_in_DNA6)) &
                                     !lead(DNA_HP1_foci) ~ 'DNA-high (subsampled)',
                                   T ~ 'Other')) %>%
    filter(Compartment != 'Other') %>%
    select(Dish, Pos, Track, Frame, x, y, Track_length, Compartment) -> result
  
  result %>% group_by(Track) %>% arrange(Frame) %>% mutate(gap = Frame - lag(Frame)) %>%
    ungroup() %>% arrange(Track) -> result
  result$Track -> track_ids
  for(i in 1:nrow(result)) {
    if (!is.na(result$gap[i]) & result$gap[i] > 1) {
      track_ids[i:length(track_ids)] = track_ids[i:length(track_ids)]+1
    }
  }
  result$Track = track_ids
  result %>% select(-gap) %>% group_by(Track) %>% mutate(Track = cur_group_id()) -> result
  
  return(result)
  toc()
}) -> for_SASPT_DNAsub
names(for_SASPT_DNAsub) = names(data_list)

lapply(for_SASPT, function(Pos) bind_rows(Pos, .id = 'Compartment')) %>% bind_rows() %>%
  group_by(Compartment, Pos, Track) %>% summarise(num_jumps = n()-1) %>%
  group_by(Compartment) %>% summarise(num_tracks = n(), num_jumps = sum(num_jumps)) %>% View()
#Smallest sample: DNA&HP1-high - 3,142 traj
bind_rows(for_SASPT_DNAsub) %>% group_by(Pos, Track) %>% summarise(num_jumps = n()-1) %>%
  ungroup() %>% summarise(num_tracks = n(), num_jumps = sum(num_jumps))
#5,794 traj
}

##Write files for SA-SPT
{
  numframes = 40000
  
  #DNA_foci
  maxtraj = sapply(for_SASPT, function(Pos) max(Pos[['DNA-high']]$Track))
  lapply(seq_along(for_SASPT), function(i) {
    for_SASPT[[i]][['DNA-high']] %>% select(Track, x, y, Frame) %>% mutate(x = x/100, y = y/100) %>%
      rename(trajectory = Track, frame = Frame) %>%
      mutate(trajectory = trajectory + sum(maxtraj[1:i-1]), frame = frame + numframes*(i-1)) -> result
    return(result)
  }) %>% bind_rows() %>%
    write.table(paste0('New_analysis/SA-SPT/Data/DNA_foci.csv'),
                col.names = T, row.names = F, quote = F, sep = ',')
  
  #DNA foci subsampled
  maxtraj = sapply(for_SASPT_DNAsub, function(Pos) max(Pos$Track))
  lapply(seq_along(for_SASPT_DNAsub), function(i) {
    for_SASPT_DNAsub[[i]] %>% select(Track, x, y, Frame) %>% mutate(x = x/100, y = y/100) %>%
      rename(trajectory = Track, frame = Frame) %>%
      mutate(trajectory = trajectory + sum(maxtraj[1:i-1]), frame = frame + numframes*(i-1)) -> result
    return(result)
  }) %>% bind_rows() %>%
    write.table(paste0('New_analysis/SA-SPT/Data/DNA_foci_sub.csv'),
                col.names = T, row.names = F, quote = F, sep = ',')
  
  
  #HP1_foci
  maxtraj = sapply(for_SASPT, function(Pos) max(Pos[['HP1-high']]$Track))
  lapply(seq_along(for_SASPT), function(i) {
    for_SASPT[[i]][['HP1-high']] %>% select(Track, x, y, Frame) %>% mutate(x = x/100, y = y/100) %>%
      rename(trajectory = Track, frame = Frame) %>%
      mutate(trajectory = trajectory + sum(maxtraj[1:i-1]), frame = frame + numframes*(i-1)) -> result
    return(result)
  }) %>% bind_rows() %>% #View()
    write.table(paste0('New_analysis/SA-SPT/Data/HP1_foci.csv'),
                col.names = T, row.names = F, quote = F, sep = ',')
  
  #DNA&HP1_foci
  maxtraj = sapply(for_SASPT, function(Pos) max(Pos[['DNA&HP1-high']]$Track))
  lapply(seq_along(for_SASPT), function(i) {
    if (!is.null(for_SASPT[[i]][['DNA&HP1-high']])) {
      for_SASPT[[i]][['DNA&HP1-high']] %>% select(Track, x, y, Frame) %>% mutate(x = x/100, y = y/100) %>%
        rename(trajectory = Track, frame = Frame) %>%
        mutate(trajectory = trajectory + sum(maxtraj[1:i-1][maxtraj[1:i-1] > 0]), frame = frame + numframes*(i-1)) -> result
      return(result)
    }
  }) %>% bind_rows() %>%
    write.table(paste0('New_analysis/SA-SPT/Data/DNA&HP1_foci.csv'),
                col.names = T, row.names = F, quote = F, sep = ',')
  
  #Euchromatin
  maxtraj = sapply(for_SASPT, function(Pos) max(Pos[['Euchromatin']]$Track))
  lapply(seq_along(for_SASPT), function(i) {
    for_SASPT[[i]][['Euchromatin']] %>% select(Track, x, y, Frame) %>% mutate(x = x/100, y = y/100) %>%
      rename(trajectory = Track, frame = Frame) %>%
      mutate(trajectory = trajectory + sum(maxtraj[1:i-1]), frame = frame + numframes*(i-1)) -> result
    return(result)
  }) %>% bind_rows() %>%# group_by(trajectory) %>% mutate(JD = sqrt((x-lag(x))^2 + (y-lag(y))^2)) %>% View()
    write.table(paste0('New_analysis/SA-SPT/Data/Euchromatin.csv'),
                col.names = T, row.names = F, quote = F, sep = ',')
  
  #DNA_pseudo
  maxtraj = sapply(for_SASPT, function(Pos) max(Pos[['DNA pseudo']]$Track))
  lapply(seq_along(for_SASPT), function(i) {
    if (!is.null(for_SASPT[[i]][['DNA pseudo']])) {
      for_SASPT[[i]][['DNA pseudo']] %>% select(Track, x, y, Frame) %>% mutate(x = x/100, y = y/100) %>%
        rename(trajectory = Track, frame = Frame) %>%
        mutate(trajectory = trajectory + sum(maxtraj[1:i-1][maxtraj[1:i-1] > 0]), frame = frame + numframes*(i-1)) -> result
      return(result)
    }
  }) %>% bind_rows() %>% #View()
    write.table(paste0('New_analysis/SA-SPT/Data/DNA_pseudo.csv'),
                col.names = T, row.names = F, quote = F, sep = ',')
  
  #HP1_pseudo
  maxtraj = sapply(for_SASPT, function(Pos) max(Pos[['HP1 pseudo']]$Track))
  lapply(seq_along(for_SASPT), function(i) {
    if (!is.null(for_SASPT[[i]][['HP1 pseudo']])) {
      for_SASPT[[i]][['HP1 pseudo']] %>% select(Track, x, y, Frame) %>% mutate(x = x/100, y = y/100) %>%
        rename(trajectory = Track, frame = Frame) %>%
        mutate(trajectory = trajectory + sum(maxtraj[1:i-1][maxtraj[1:i-1] > 0]), frame = frame + numframes*(i-1)) -> result
      return(result)
    }
  }) %>% bind_rows() %>% #View()
    write.table(paste0('New_analysis/SA-SPT/Data/HP1_pseudo.csv'),
                col.names = T, row.names = F, quote = F, sep = ',')
  
  
  #DNA&HP1_pseudo
  maxtraj = sapply(for_SASPT, function(Pos) max(Pos[['DNA&HP1 pseudo']]$Track))
  lapply(seq_along(for_SASPT), function(i) {
    if (!is.null(for_SASPT[[i]][['DNA&HP1 pseudo']])) {
      for_SASPT[[i]][['DNA&HP1 pseudo']] %>% select(Track, x, y, Frame) %>% mutate(x = x/100, y = y/100) %>%
        rename(trajectory = Track, frame = Frame) %>%
        mutate(trajectory = trajectory + sum(maxtraj[1:i-1][maxtraj[1:i-1] > 0]), frame = frame + numframes*(i-1)) -> result
      return(result)
    }
  }) %>% bind_rows() %>% #View()
    write.table(paste0('New_analysis/SA-SPT/Data/DNA&HP1_pseudo.csv'),
                col.names = T, row.names = F, quote = F, sep = ',')
}


#### Calculate the density of molecules within HP1 foci and HP1 pseudofoci; and within DNA foci and pseudofoci ####

##Load data, redo classification
{
data_list = readRDS('Data_list_SPT_segm.RDS')
data = bind_rows(data_list)
data %>% filter(Dish == 'Dish1' & Pos == 'Pos0' & Frame > 19500 |
                  Dish == 'Dish1' & Pos == 'Pos1' & Frame > 20000 |
                  Dish == 'Dish1' & Pos == 'Pos2' & Frame > 20000 |
                  Dish == 'Dish1' & Pos == 'Pos3' & Frame > 10000 |
                  Dish == 'Dish1' & Pos == 'Pos4' & Frame > 10000 |
                  Dish == 'Dish1' & Pos == 'Pos5' & Frame > 21000 |
                  Dish == 'Dish1' & Pos == 'Pos6' & Frame > 12000 |
                  Dish == 'Dish1' & Pos == 'Pos7' & Frame > 20000 |
                  Dish == 'Dish1' & Pos == 'Pos8' & Frame > 11000 |
                  Dish == 'Dish1' & Pos == 'Pos9' & Frame > 21000 |
                  Dish == 'Dish1' & Pos == 'Pos10' & Frame > 11000 |
                  Dish == 'Dish1' & Pos == 'Pos11' & Frame > 19000 |
                  Dish == 'Dish2' & Pos == 'Pos0' & Frame > 25000 |
                  Dish == 'Dish2' & Pos == 'Pos1' & Frame > 24000 |
                  Dish == 'Dish2' & Pos == 'Pos2' & Frame > 24000 |
                  Dish == 'Dish2' & Pos == 'Pos3' & Frame > 20000 |
                  Dish == 'Dish2' & Pos == 'Pos4' & Frame > 23000 |
                  Dish == 'Dish2' & Pos == 'Pos5' & Frame > 10000 |
                  Dish == 'Dish2' & Pos == 'Pos6' & Frame > 33000 |
                  Dish == 'Dish2' & Pos == 'Pos7' & Frame > 17000 |
                  Dish == 'Dish2' & Pos == 'Pos8' & Frame > 14000) -> data_filtered

data_list_SMLM = readRDS('Data_list_SMLM.RDS')
names(data_list_SMLM) = sapply(data_list_SMLM, function(x) x$Pos[1])

WFframes = 8
WFaveraging = 5000 #WFframes*WFaveraging = total number of frames
pixel = 140
subsampling = 10

foci_list = readRDS('DNA_foci_list.rds')
pfoci_list = readRDS('DNA_pseudofoci_list.rds')
segmCoresList = readRDS('HP1_dbscan.RDS')

n_shuffle = 6

thresholds1 = c(100, 100, 800, 200, 100, 0, #Dish1 (0, 10, 11, 1, 2, 3)
                50, 1000, 100, 1000, 600, 200, #Dish1 (4, 5, 6, 7, 8, 9)
                800, 1200, 850, 500, 1000, 200, 1100, 800, 200) #Dish2 (1-8)
##Classify HP1 clusters as overlapping/non-overlapping with DNA foci as a whole based on %overlap
{
  lapply(data_list, function(dat) {
    dat %>% group_by(HP1_foci_dbscan) %>%
      mutate(DNA_HP1_focus = case_when(HP1_foci_dbscan == 0 ~ F,
                                       sum(DNA_focus_idx!=0)/n() > 0.3 ~ T,
                                       T ~ F)) %>% ungroup() -> result
    return(unique((result %>% filter(DNA_HP1_focus))$HP1_foci_dbscan))
  }) -> DNA_HP1_foci
  names(DNA_HP1_foci) = names(data_list)
}
}

##Calculate density of HP1 in HP1 foci (DNA-rich and DNA-poor) and HP1 pseudofoci (SMLM data)
{
  lapply(segmCoresList, function(Pos) {
    Pos$Foci %>% group_by(cluster) %>%
      summarise(NNdist_mean_focus = 0.001*mean(kNNdist(matrix(c(x, y), nrow=n()), k = 1, all = F))) %>%
      mutate(density = 1/(4*NNdist_mean_focus^2)) -> result
  }) %>% bind_rows(.id = 'Pos') %>%
    separate(Pos, into = c('Dish', 'Pos'), '_') -> HP1_foci_dens_raw
  
  lapply(seq_along(segmCoresList), function(c) {
    pos = names(segmCoresList)[c]
    filter(data_list_SMLM[[c]], Frame > thresholds1[c]) %>% group_by(Dish, Pos, Track) %>%
      summarise(x = mean(x), y = mean(y)) %>% ungroup() %>% select(x, y) -> df
    
    segmentation = list('cluster' = segmCoresList[[pos]][['Foci']]$cluster,
                        'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
    class(segmentation) = c('dbscan_fast', 'dbscan')
    segm = predict(object = segmentation, newdata = df,
                   data = select(ungroup(segmCoresList[[pos]][['Foci']]), x, y))
    df %>% mutate(cluster = segm) %>% filter(cluster != 0) %>% group_by(cluster) %>%
      summarise(NNdist_mean_focus = 0.001*mean(kNNdist(matrix(c(x, y), nrow=n()), k = 1, all = F))) %>%
      mutate(density = 1/(4*NNdist_mean_focus^2)) -> result2
  }) %>% bind_rows(.id = 'Pos') -> HP1_foci_dens
  renamer = names(segmCoresList)
  names(renamer) = seq_along(segmCoresList)
  HP1_foci_dens %>% mutate(Pos = renamer[Pos]) %>% separate(Pos, into = c('Dish', 'Pos'), '_') -> HP1_foci_dens
  
  lapply(seq_along(segmCoresList), function(c) {
    pos = names(segmCoresList)[c]
    print(pos)
    filter(data_list_SMLM[[c]], Frame > thresholds1[c]) %>% group_by(Dish, Pos, Track) %>%
      summarise(x = mean(x), y = mean(y)) %>% ungroup() %>% select(x, y) -> df
    for (n in 1:n_shuffle) {
      print(n)
      pseudosegmentation = list('cluster' = segmCoresList[[pos]][[paste0('Pseudo', n)]]$cluster,
                                'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
      class(pseudosegmentation) = c('dbscan_fast', 'dbscan')
      psegm = predict(object = pseudosegmentation, newdata = df,
                      data = select(ungroup(segmCoresList[[pos]][[paste0('Pseudo', n)]]), x, y))
      df %>% mutate(cluster = psegm) %>% filter(cluster != 0) %>% group_by(cluster) %>%
        filter(n() > 1) %>% #cluster of size 1 gives an error in kNNdist
        summarise(NNdist_mean_focus = 0.001*mean(kNNdist(matrix(c(x, y), nrow=n()), k = 1, all = F))) %>%
        mutate(density = 1/(4*NNdist_mean_focus^2)) -> result2
      names(result2) = c('cluster', paste0('NNdist_mean', n), paste0('density', n))
      if (n == 1) {
        result = result2
      } else {
        result = full_join(result, result2, by = 'cluster')
      }
    }
    return(result)
  }) %>% bind_rows(.id = 'Pos') -> HP1_pfoci_dens
  HP1_pfoci_dens %>% mutate(Pos = renamer[Pos]) %>% separate(Pos, into = c('Dish', 'Pos'), '_') -> HP1_pfoci_dens
  
  data.frame(cluster = unlist(DNA_HP1_foci),
             Pos = unlist(sapply(names(DNA_HP1_foci), function(pos) rep(pos, length(DNA_HP1_foci[[pos]])))),
             is_DNA = T) %>%
    separate(Pos, into = c('Dish', 'Pos'), '_') -> DNA_HP1_foci.df
  
  HP1_foci_dens_raw %>% rename(NNdist_mean_raw = NNdist_mean_focus, density_raw = density) %>%
    left_join(rename(HP1_foci_dens, NNdist_mean_focus = NNdist_mean_focus, density_focus = density),
              by = c('Dish', 'Pos', 'cluster')) %>%
    left_join(rename(HP1_pfoci_dens, NNdist_mean_pseudo1 = NNdist_mean1, density_pseudo1 = density1,
                     NNdist_mean_pseudo2 = NNdist_mean2, density_pseudo2 = density2,
                     NNdist_mean_pseudo3 = NNdist_mean3, density_pseudo3 = density3,
                     NNdist_mean_pseudo4 = NNdist_mean4, density_pseudo4 = density4,
                     NNdist_mean_pseudo5 = NNdist_mean5, density_pseudo5 = density5,
                     NNdist_mean_pseudo6 = NNdist_mean6, density_pseudo6 = density6),
              by = c('Dish', 'Pos', 'cluster')) %>% 
    left_join(DNA_HP1_foci.df, by = c('Dish', 'Pos', 'cluster')) %>%
    mutate(is_DNA = replace_na(is_DNA, F)) -> HP1_all_dens
  
  HP1_all_dens %>% mutate(across(starts_with('density'), ~replace_na(.x, 0))) -> HP1_all_dens
}

##Summarise and plot HP1 density in foci
{
  ##Density histograms
  pivot_longer(HP1_all_dens, starts_with('density'), names_to = 'Category', values_to = 'density') %>%
    mutate(Category2 = case_when(grepl('pseudo', Category) ~ 'density_pseudo',
                                 T ~ Category)) %>%
    ggplot(aes(x = density, colour = Category2, group = Category)) +
    geom_density() +
    scale_colour_manual(values = c('red', 'tomato4', 'orange'),
                        labels = c('Foci', 'Pseudofoci', 'Foci cores')) +
    coord_cartesian(xlim = c(0, 6000)) +
    labs(x = 'Density, mol/um2', y = 'Frequency', colour = NULL, title = 'Density of HP1 foci (DBSCAN)') +
    theme_bw() +
    theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
  
  pivot_longer(HP1_all_dens, starts_with('density'), names_to = 'Category', values_to = 'density') %>%
    mutate(Category2 = case_when(grepl('pseudo', Category) ~ 'density_pseudo',
                                 T ~ Category)) %>%
    ggplot(aes(x = density, colour = Category2, linetype = is_DNA, group = interaction(Category, is_DNA))) +
    geom_density() +
    scale_colour_manual(values = c('red', 'tomato4', 'orange'),
                        labels = c('Foci', 'Pseudofoci', 'Foci cores')) +
    scale_linetype_discrete(labels = c('DNA-low', 'DNA-high')) +
    coord_cartesian(xlim = c(0, 6000)) +
    labs(x = 'Density, mol/um2', y = 'Frequency', colour = NULL, linetype = NULL,
         title = 'Density of HP1 foci (DBSCAN)') +
    theme_bw() +
    theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
  
  ##Ratio (distribution)
  mutate(HP1_all_dens, ratio_core1 = density_raw/density_pseudo1, ratio_core2 = density_raw/density_pseudo2,
         ratio_core3 = density_raw/density_pseudo3, ratio_core4 = density_raw/density_pseudo4,
         ratio_core5 = density_raw/density_pseudo5, ratio_core6 = density_raw/density_pseudo6,
         ratio_focus1 = density_focus/density_pseudo1, ratio_focus2 = density_focus/density_pseudo2,
         ratio_focus3 = density_focus/density_pseudo3, ratio_focus4 = density_focus/density_pseudo4,
         ratio_focus5 = density_focus/density_pseudo5, ratio_focus6 = density_focus/density_pseudo6) -> HP1_all_dens
  pivot_longer(HP1_all_dens, starts_with('ratio'), names_to = 'Category', values_to = 'ratio') %>%
    mutate(Category2 = case_when(grepl('core', Category) ~ 'core',
                                 T ~ 'focus')) %>%
    ggplot(aes(x = ratio, colour = Category2, group = Category)) +
    geom_density() +
    scale_colour_manual(values = c('orange', 'red'),
                        labels = c('Foci cores/Pseudofoci', 'Foci/Pseudofoci')) +
    coord_cartesian(xlim = c(0, 10)) +
    labs(x = 'Ratio', y = 'Frequency', colour = NULL,
         title = 'Ratio of density of HP1 foci\n&pseudofoci (DBSCAN)') +
    theme_bw() +
    theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
  pivot_longer(HP1_all_dens, starts_with('ratio'), names_to = 'Category', values_to = 'ratio') %>%
    mutate(Category2 = case_when(grepl('core', Category) ~ 'core',
                                 T ~ 'focus')) %>%
    ggplot(aes(x = ratio, colour = Category2, linetype = is_DNA, group = interaction(is_DNA, Category))) +
    geom_density() +
    scale_colour_manual(values = c('orange', 'red'),
                        labels = c('Foci cores/Pseudofoci', 'Foci/Pseudofoci')) +
    scale_linetype_discrete(labels = c('DNA-low', 'DNA-high')) +
    coord_cartesian(xlim = c(0, 10)) +
    labs(x = 'Ratio', y = 'Frequency', colour = NULL, linetype = NULL,
         title = 'Ratio of density of HP1 foci\n&pseudofoci (DBSCAN)') +
    theme_bw() +
    theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
  
  ##Density(mean)
  HP1_all_dens %>% summarise_all(mean, na.rm = T) %>%
    pivot_longer(starts_with('density'), names_to = 'Category', values_to = 'density') %>%
    mutate(Category2 = case_when(grepl('pseudo', Category) ~ 'pseudo',
                                 T ~ Category)) -> HP1_all_dens_summary
  ggplot(HP1_all_dens_summary, aes(x = Category2, y = density, group = Category, fill = Category2)) +
    geom_bar(stat = 'identity', position = 'dodge', colour = 'gray20') +
    geom_text(data = (group_by(HP1_all_dens_summary, Category2) %>%
                        summarise(density = mean(density), Category = Category[1])),
              aes(label = round(density, 2)), vjust = -0.5) +
    scale_x_discrete(labels = c('Focus', 'Core focus', 'Pseudo')) +
    scale_fill_manual(values = c('red', 'orange', 'tomato4'), guide = 'none') +
    labs(x = NULL, y = 'Density, mol/um2', title = 'Mean density of HP1 foci (SMLM)') +
    lims(y = c(0, 4250)) +
    theme_bw() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
  ##Subdivided into DNA+ and -
  HP1_all_dens %>% group_by(is_DNA) %>% summarise_all(mean, na.rm = T) %>%
    pivot_longer(starts_with('density'), names_to = 'Category', values_to = 'density') %>%
    mutate(Category2 = case_when(grepl('pseudo', Category) ~ 'pseudo',
                                 T ~ Category)) -> HP1_all_dens_summary
  ggplot(HP1_all_dens_summary,
         aes(x = interaction(is_DNA, Category2), y = density, group = Category,
             fill = interaction(is_DNA, Category2))) +
    geom_bar(stat = 'identity', position = 'dodge', colour = 'gray20') +
    geom_text(data = (group_by(HP1_all_dens_summary, Category2, is_DNA) %>%
                        summarise(density = mean(density), Category = Category[1])),
              aes(label = round(density, 2)), vjust = -0.5, size = 3) +
    scale_x_discrete(labels = c('HP1', 'HP1&DNA', 'HP1 core', 'HP1&DNA core', 'HP1 pseudo', 'HP1&DNA pseudo')) +
    scale_fill_manual(values = c('red', 'purple', 'orange', 'magenta', 'tomato4', 'purple4'), guide = 'none') +
    labs(x = NULL, y = 'Density, mol/um2', title = 'Mean density of HP1 foci (SMLM)',
         subtitle = '1/(4*NN_dist^2)') +
    lims(y = c(0, 4500)) +
    theme_bw() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave('New_analysis/HP1_foci_density.png', width = 3.5, height = 4)
}
HP1_all_dens_summary #values used in saspt_plotting to create Fig 4D

##Calculate the same, but from SPT data - sanity check (not used)
{
  # data_filtered %>% filter(HP1_foci_dbscan != 0) %>% group_by(Dish, Pos, HP1_foci_dbscan) %>%
  #   summarise(NNdist_mean_focus = 0.001*mean(kNNdist(matrix(c(x, y), nrow=n()), k = 1, all = F))) %>%
  #   mutate(density = 1/(4*NNdist_mean_focus^2)) -> HP1_foci_dens
  # 
  # data_filtered %>% filter(HP1_pseudo_dbscan1 != 0) %>% group_by(Dish, Pos, HP1_pseudo_dbscan1) %>%
  #   filter(n() > 1) %>% #cluster of size 1 gives an error in kNNdist
  #   summarise(NNdist_mean_pseudo1 = 0.001*mean(kNNdist(matrix(c(x, y), nrow=n()), k = 1, all = F))) %>%
  #   mutate(density_pseudo1 = 1/(4*NNdist_mean_pseudo1^2)) -> HP1_pseudo_dens1
  # data_filtered %>% filter(HP1_pseudo_dbscan2 != 0) %>% group_by(Dish, Pos, HP1_pseudo_dbscan2) %>%
  #   filter(n() > 1) %>% #cluster of size 1 gives an error in kNNdist
  #   summarise(NNdist_mean_pseudo2 = 0.001*mean(kNNdist(matrix(c(x, y), nrow=n()), k = 1, all = F))) %>%
  #   mutate(density_pseudo2 = 1/(4*NNdist_mean_pseudo2^2)) -> HP1_pseudo_dens2
  # data_filtered %>% filter(HP1_pseudo_dbscan3 != 0) %>% group_by(Dish, Pos, HP1_pseudo_dbscan3) %>%
  #   filter(n() > 1) %>% #cluster of size 1 gives an error in kNNdist
  #   summarise(NNdist_mean_pseudo3 = 0.001*mean(kNNdist(matrix(c(x, y), nrow=n()), k = 1, all = F))) %>%
  #   mutate(density_pseudo3 = 1/(4*NNdist_mean_pseudo3^2)) -> HP1_pseudo_dens3
  # data_filtered %>% filter(HP1_pseudo_dbscan4 != 0) %>% group_by(Dish, Pos, HP1_pseudo_dbscan4) %>%
  #   filter(n() > 1) %>% #cluster of size 1 gives an error in kNNdist
  #   summarise(NNdist_mean_pseudo4 = 0.001*mean(kNNdist(matrix(c(x, y), nrow=n()), k = 1, all = F))) %>%
  #   mutate(density_pseudo4 = 1/(4*NNdist_mean_pseudo4^2)) -> HP1_pseudo_dens4
  # data_filtered %>% filter(HP1_pseudo_dbscan5 != 0) %>% group_by(Dish, Pos, HP1_pseudo_dbscan5) %>%
  #   filter(n() > 1) %>% #cluster of size 1 gives an error in kNNdist
  #   summarise(NNdist_mean_pseudo5 = 0.001*mean(kNNdist(matrix(c(x, y), nrow=n()), k = 1, all = F))) %>%
  #   mutate(density_pseudo5 = 1/(4*NNdist_mean_pseudo5^2)) -> HP1_pseudo_dens5
  # data_filtered %>% filter(HP1_pseudo_dbscan6 != 0) %>% group_by(Dish, Pos, HP1_pseudo_dbscan6) %>%
  #   filter(n() > 1) %>% #cluster of size 1 gives an error in kNNdist
  #   summarise(NNdist_mean_pseudo6 = 0.001*mean(kNNdist(matrix(c(x, y), nrow=n()), k = 1, all = F))) %>%
  #   mutate(density_pseudo6 = 1/(4*NNdist_mean_pseudo6^2)) -> HP1_pseudo_dens6
  # 
  # left_join(HP1_foci_dens, HP1_pseudo_dens1, by = c('Dish' = 'Dish', 'Pos' = 'Pos', 'HP1_foci_dbscan' = 'HP1_pseudo_dbscan1')) %>%
  #   left_join(HP1_pseudo_dens2, by = c('Dish' = 'Dish', 'Pos' = 'Pos', 'HP1_foci_dbscan' = 'HP1_pseudo_dbscan2')) %>%
  #   left_join(HP1_pseudo_dens3, by = c('Dish' = 'Dish', 'Pos' = 'Pos', 'HP1_foci_dbscan' = 'HP1_pseudo_dbscan3')) %>%
  #   left_join(HP1_pseudo_dens4, by = c('Dish' = 'Dish', 'Pos' = 'Pos', 'HP1_foci_dbscan' = 'HP1_pseudo_dbscan4')) %>%
  #   left_join(HP1_pseudo_dens5, by = c('Dish' = 'Dish', 'Pos' = 'Pos', 'HP1_foci_dbscan' = 'HP1_pseudo_dbscan5')) %>%
  #   left_join(HP1_pseudo_dens6, by = c('Dish' = 'Dish', 'Pos' = 'Pos', 'HP1_foci_dbscan' = 'HP1_pseudo_dbscan6')) %>%
  #   rename(cluster = HP1_foci_dbscan) %>% left_join(DNA_HP1_foci.df, by = c('Dish', 'Pos', 'cluster')) %>%
  #   mutate(is_DNA = replace_na(is_DNA, F)) -> HP1_all_dens
# 
# ##Plot
#   ##Density histograms
#   pivot_longer(HP1_all_dens, starts_with('density'), names_to = 'Category', values_to = 'density') %>%
#     mutate(Category2 = case_when(grepl('pseudo', Category) ~ 'density_pseudo',
#                                  T ~ Category)) %>%
#     ggplot(aes(x = density, colour = Category2, group = Category)) +
#     geom_density() +
#     scale_colour_manual(values = c('red', 'tomato4'),
#                         labels = c('Foci', 'Pseudofoci')) +
#     coord_cartesian(xlim = c(0, 6000)) +
#     labs(x = 'Density, mol/um2', y = 'Frequency', colour = NULL, title = 'Density of HP1 foci (DBSCAN)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
#   
#   pivot_longer(HP1_all_dens, starts_with('density'), names_to = 'Category', values_to = 'density') %>%
#     mutate(Category2 = case_when(grepl('pseudo', Category) ~ 'density_pseudo',
#                                  T ~ Category)) %>%
#     ggplot(aes(x = density, colour = Category2, linetype = is_DNA, group = interaction(Category, is_DNA))) +
#     geom_density() +
#     scale_colour_manual(values = c('red', 'tomato4'),
#                         labels = c('Foci', 'Pseudofoci')) +
#     scale_linetype_discrete(labels = c('DNA-low', 'DNA-high')) +
#     coord_cartesian(xlim = c(0, 6000)) +
#     labs(x = 'Density, mol/um2', y = 'Frequency', colour = NULL, linetype = NULL,
#          title = 'Density of HP1 foci (DBSCAN)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
#   
#   ##Ratio (distribution)
#   mutate(HP1_all_dens, ratio1 = density/density_pseudo1, ratio2 = density/density_pseudo2,
#          ratio3 = density/density_pseudo3, ratio4 = density/density_pseudo4,
#          ratio5 = density/density_pseudo5, ratio6 = density/density_pseudo6) -> HP1_all_dens
#   pivot_longer(HP1_all_dens, starts_with('ratio'), names_to = 'Category', values_to = 'ratio') %>%
#     ggplot(aes(x = ratio, group = Category)) +
#     geom_density(n = 2^14, colour = 'red') +
#     coord_cartesian(xlim = c(0, 15)) +
#     labs(x = 'Ratio', y = 'Frequency', colour = NULL,
#          title = 'Ratio of density of HP1 foci\n&pseudofoci (DBSCAN)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
#   pivot_longer(HP1_all_dens, starts_with('ratio'), names_to = 'Category', values_to = 'ratio') %>%
#     ggplot(aes(x = ratio, group = interaction(Category, is_DNA), colour = is_DNA)) +
#     geom_density(n = 2^14) +
#     scale_colour_manual(values = c('red', 'purple'), labels = c('DNA-low', 'DNA-high')) +
#     coord_cartesian(xlim = c(0, 15)) +
#     labs(x = 'Ratio', y = 'Frequency', colour = NULL, linetype = NULL,
#          title = 'Ratio of density of HP1 foci\n&pseudofoci (DBSCAN)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
#   
#   ##Density(mean)
#   ungroup(HP1_all_dens) %>% summarise_all(mean, na.rm = T) %>%
#     pivot_longer(starts_with('density'), names_to = 'Category', values_to = 'density') %>%
#     mutate(Category2 = case_when(grepl('pseudo', Category) ~ 'pseudo',
#                                  T ~ Category)) -> HP1_all_dens_summary
#   ggplot(HP1_all_dens_summary, aes(x = Category2, y = density, group = Category, fill = Category2)) +
#     geom_bar(stat = 'identity', position = 'dodge', colour = 'gray20') +
#     geom_text(data = (group_by(HP1_all_dens_summary, Category2) %>%
#                         summarise(density = mean(density), Category = Category[1])),
#               aes(label = round(density, 2)), vjust = -0.5) +
#     scale_x_discrete(labels = c('Focus', 'Pseudo')) +
#     scale_fill_manual(values = c('red', 'tomato4'), guide = 'none') +
#     labs(x = NULL, y = 'Density, mol/um2', title = 'Mean density of HP1 foci (SPT)') +
#     lims(y = c(0, 4200)) +
#     theme_bw() +
#     theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
#   ##Subdivided into DNA+ and -
#   ungroup(HP1_all_dens) %>% group_by(is_DNA) %>% summarise_all(mean, na.rm = T) %>%
#     pivot_longer(starts_with('density'), names_to = 'Category', values_to = 'density') %>%
#     mutate(Category2 = case_when(grepl('pseudo', Category) ~ 'pseudo',
#                                  T ~ Category)) -> HP1_all_dens_summary
#   ggplot(HP1_all_dens_summary,
#          aes(x = interaction(is_DNA, Category2), y = density, group = Category,
#              fill = interaction(is_DNA, Category2))) +
#     geom_bar(stat = 'identity', position = 'dodge', colour = 'gray20') +
#     geom_text(data = (group_by(HP1_all_dens_summary, Category2, is_DNA) %>%
#                         summarise(density = mean(density), Category = Category[1])),
#               aes(label = round(density, 2)), vjust = -0.5) +
#     scale_x_discrete(labels = c('HP1', 'HP1&DNA', 'HP1 pseudo', 'HP1&DNA pseudo')) +
#     scale_fill_manual(values = c('red', 'purple', 'tomato4', 'purple4'), guide = 'none') +
#     labs(x = NULL, y = 'Density, mol/um2', title = 'Mean density of HP1 foci (SPT)',
#          subtitle = '1/(4*NN_dist^2)') +
#     lims(y = c(0, 4500)) +
#     theme_bw() +
#     theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           axis.text.x = element_text(angle = 30, hjust = 1))
#   ggsave('New_analysis/HP1_foci_density2.png', width = 3.5, height = 4)
}

###DNA
##Calculate density of HP1 in DNA foci (HP1-rich and HP1-poor) and pseudofoci (SMLM data)
{
  #Classify points
  for (c in seq_along(data_list_SMLM)) {
    pos = names(data_list_SMLM)[c]
    print(pos)
    filter(data_list_SMLM[[c]], Frame > thresholds1[c]) %>% group_by(Dish, Pos, Track) %>%
      summarise(x = mean(x), y = mean(y), Frame = mean(Frame), Pos = Pos[1], Dish = Dish[1]) %>% ungroup() %>%
      select(Dish, Pos, Frame, x, y) -> df
    
    #Classify by DNA
    dim_x = foci_list[[pos]]$dim_x[1]
    dim_y = foci_list[[pos]]$dim_y[1]
    df %>%
      mutate(px_id = floor(x*subsampling/pixel)*dim_y + ceiling(y*subsampling/pixel),
             frame = ceiling(Frame/WFaveraging)) %>%
      left_join(., foci_list[[pos]], by = c('px_id', 'frame')) %>%
      select(-x_left, -y_top, -dim_y, -dim_x, -pos, -dish) %>%
      left_join(., pfoci_list[[pos]], by = c('px_id', 'frame')) %>%
      select(-x_left, -y_top, -dim_y, -dim_x, -pos, -dish, -frame, -px_id) %>%
      rename(DNA_focus_idx = focus_id, DNA_pseudo = pseudofocus) %>%
      mutate(DNA_focus_idx = replace_na(DNA_focus_idx, 0),
             DNA_pseudo = replace_na(DNA_pseudo, 0)) -> tmp
    filter(tmp, DNA_focus_idx != 0) -> tmp1
    filter(tmp, DNA_pseudo != 0) -> tmp2
    
    #Classify by HP1
    segmentation = list('cluster' = segmCoresList[[pos]][['Foci']]$cluster,
                        'eps' = 50, 'minPts' = 0, dist = 'euclidean', borderPoints = T)
    class(segmentation) = c('dbscan_fast', 'dbscan')
    segm = predict(object = segmentation, newdata = select(ungroup(tmp1), x, y),
                   data = select(ungroup(segmCoresList[[pos]][['Foci']]), x, y))
    tmp1$HP1_focus_idx = segm
    segm2 = predict(object = segmentation, newdata = select(ungroup(tmp2), x, y),
                    data = select(ungroup(segmCoresList[[pos]][['Foci']]), x, y))
    tmp2$HP1_focus_idx = segm2
    
    
    tmp1 %>% group_by(DNA_focus_idx) %>%
      summarise(NNdist_mean_focus = 0.001*mean(kNNdist(matrix(c(x, y), nrow=n()), k = 1, all = F))) %>%
      mutate(density_focus1 = 1/(4*NNdist_mean_focus^2)) -> result1
    filter(tmp1, HP1_focus_idx == 0) %>% group_by(DNA_focus_idx) %>%
      summarise(NNdist_mean_focus_noHP1 = 0.001*mean(kNNdist(matrix(c(x, y), nrow=n()), k = 1, all = F))) %>%
      mutate(density_focus1_noHP1 = 1/(4*NNdist_mean_focus_noHP1^2)) %>%
      right_join(result1, by = 'DNA_focus_idx') -> result1
    
    tmp1 %>% group_by(DNA_focus_idx) %>% summarise(num_obs = n()) %>%
      left_join(group_by(foci_list[[pos]], focus_id) %>% summarise(num_px = n()/WFframes),
                by = c('DNA_focus_idx' = 'focus_id')) %>%
      mutate(density_focus2 = num_obs/(num_px*(pixel/subsampling/1000)^2)) -> result2
    filter(tmp1, HP1_focus_idx == 0) %>% group_by(DNA_focus_idx) %>% summarise(num_obs_noHP1 = n()) %>%
      left_join(result2, by = 'DNA_focus_idx') %>%
      mutate(density_focus2_noHP1 = num_obs_noHP1/(num_px*(pixel/subsampling/1000)^2)) -> result2
    
    left_join(result1, result2, by = 'DNA_focus_idx') -> result_foci
    
    tmp2 %>% summarise(num_obs = n(), num_px = nrow(pfoci_list[[pos]])/WFframes) %>%
      mutate(density_pseudo = num_obs/(num_px*(pixel/subsampling/1000)^2)) -> result_pfoci
    filter(tmp2, HP1_focus_idx == 0) %>% summarise(num_obs_noHP1 = n()) %>% cbind(result_pfoci, .) %>%
      mutate(density_pseudo_noHP1 = num_obs_noHP1/(num_px*(pixel/subsampling/1000)^2)) -> result_pfoci
    
    if (c == 1) {
      DNA_foci_dens = mutate(result_foci, Pos = pos)
      DNA_pseudo_dens = mutate(result_pfoci, Pos = pos)
    } else {
      DNA_foci_dens = rbind(DNA_foci_dens, mutate(result_foci, Pos = pos))
      DNA_pseudo_dens = rbind(DNA_pseudo_dens, mutate(result_pfoci, Pos = pos))
    }
  }
  
  #Sanity checks
  sum(DNA_foci_dens$num_px) #not equal because of Dish1_Pos2 thing
  sum(DNA_pseudo_dens$num_px)
  
  mean(DNA_foci_dens$density_focus1)
  mean(DNA_foci_dens$density_focus2)
  mean(DNA_foci_dens$density_focus1_noHP1)
  mean(DNA_foci_dens$density_focus2_noHP1)
  mean(DNA_pseudo_dens$density_pseudo)
  mean(DNA_pseudo_dens$density_pseudo_noHP1)
}
#difference between 1/(4*NN_dist^2) and 'Num_obs/Area' ~ 0.65-0.7-fold (NNdist-based is lower)

##Summarise and plot HP1 density in DNA foci
{
  ##Density histograms
  pivot_longer(DNA_foci_dens, starts_with('density'), names_to = 'Category', values_to = 'density') %>%
    mutate(Category2 = case_when(grepl('focus1', Category) ~ '1/(4*NN_dist^2)',
                                 T ~ 'Num_obs/Area')) %>%
    mutate(Category3 = case_when(grepl('noHP1', Category) ~ 'HP1-low',
                                 T ~ 'All')) %>%
    ggplot(aes(x = density, colour = Category3)) +
    geom_density() +
    facet_grid(Category2~.) +
    scale_colour_manual(values = c('magenta', 'cyan')) +
    #coord_cartesian(xlim = c(0, 800)) +
    labs(x = 'Density, mol/um2', y = 'Frequency', colour = NULL, title = 'Density of HP1 in DNA foci') +
    theme_bw() +
    theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
  
  DNA_foci_dens %>% group_by(Pos) %>%
    summarise(num_obs = sum(num_obs), num_obs_noHP1 = sum(num_obs_noHP1), num_px = sum(num_px)) %>%
    mutate(density_focus = num_obs/(num_px*(pixel/subsampling/1000)^2),
           density_focus_noHP1 = num_obs_noHP1/(num_px*(pixel/subsampling/1000)^2)) %>%
    left_join(DNA_pseudo_dens, by = 'Pos') -> DNA_all_dens
  
  pivot_longer(DNA_all_dens, starts_with('density'), names_to = 'Category', values_to = 'density') %>%
    mutate(Category2 = case_when(grepl('focus', Category) ~ 'Foci',
                                 T ~ 'Pseudofoci')) %>%
    mutate(Category3 = case_when(grepl('noHP1', Category) ~ 'HP1-low',
                                 T ~ 'All')) %>%
    ggplot(aes(x = density, colour = Category2, linetype = Category3)) +
    geom_density() +
    scale_colour_manual(values = c('cyan', 'blue')) +
    #coord_cartesian(xlim = c(0, 800)) +
    labs(x = 'Density, mol/um2', y = 'Frequency', colour = NULL, title = 'Density of DNA foci') +
    theme_bw() +
    theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
  
  ##Ratio (distribution)
  mutate(DNA_all_dens, ratio = density_focus/density_pseudo,
         ratio_noHP1 = density_focus_noHP1/density_pseudo_noHP1) -> DNA_all_dens
  ggplot(DNA_all_dens, aes(x = ratio)) +
    geom_density(colour = 'cyan') +
    labs(x = 'Ratio', y = 'Frequency', colour = NULL,
         title = 'Ratio of density of DNA foci\n&pseudofoci (per Pos)') +
    theme_bw() +
    theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
  
  pivot_longer(DNA_all_dens, starts_with('ratio'), names_to = 'Category', values_to = 'ratio') %>%
    mutate(Category = case_when(grepl('noHP1', Category) ~ 'HP1-low',
                                T ~ 'All')) %>%
    ggplot(aes(x = ratio, colour = Category)) +
    geom_density() +
    scale_colour_manual(values = c('magenta', 'cyan')) +
    labs(x = 'Ratio', y = 'Frequency', colour = NULL, linetype = NULL,
         title = 'Ratio of density of DNA foci\n&pseudofoci (per Pos)') +
    theme_bw() +
    theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
  
  ##Density(mean)
  DNA_all_dens %>% summarise_all(mean, na.rm = T) %>%
    pivot_longer(starts_with('density'), names_to = 'Category', values_to = 'density') -> DNA_all_dens_summary
  ggplot(DNA_all_dens_summary, aes(x = Category, y = density, fill = Category)) +
    geom_bar(stat = 'identity', position = 'dodge', colour = 'gray20') +
    geom_text(aes(label = round(density, 2)), vjust = -0.5) +
    scale_x_discrete(labels = c('DNA all', 'DNA HP1-low', 'DNA pseudo all', 'DNA pseudo HP1-low')) +
    scale_fill_manual(values = c('magenta', 'cyan', 'magenta3', 'blue'), guide = 'none') +
    labs(x = NULL, y = 'Density, mol/um2', title = 'Mean density of HP1 (SMLM)\nin DNA foci',
         subtitle = 'Num obs/area of foci') +
    lims(y = c(0, 2300)) +
    theme_bw() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave('New_analysis/DNA_foci_density.png', width = 3.5, height = 4)
}
DNA_all_dens_summary #values used in saspt_plotting to create Fig 4D

##Density of DNA foci and pseudofoci (SPT data) - sanity check (not used)
{
#   data_filtered %>% filter(DNA_focus_idx != 0) -> tmp1
#   tmp1 %>% group_by(Dish, Pos, DNA_focus_idx) %>%
#     summarise(NNdist_mean_focus = 0.001*mean(kNNdist(matrix(c(x, y), nrow=n()), k = 1, all = F))) %>%
#     mutate(density_focus1 = 1/(4*NNdist_mean_focus^2)) -> DNA_foci_dens
#   filter(tmp1, HP1_foci_dbscan == 0) %>% group_by(Dish, Pos, DNA_focus_idx) %>%
#     summarise(NNdist_mean_focus_noHP1 = 0.001*mean(kNNdist(matrix(c(x, y), nrow=n()), k = 1, all = F))) %>%
#     mutate(density_focus1_noHP1 = 1/(4*NNdist_mean_focus_noHP1^2)) %>%
#     right_join(DNA_foci_dens, by = c('Dish', 'Pos', 'DNA_focus_idx')) -> DNA_foci_dens
#   
#   tmp1 %>% group_by(Dish, Pos, DNA_focus_idx) %>% summarise(num_obs = n()) %>%
#     left_join(group_by(bind_rows(foci_list, .id = 'Pos'), Pos, focus_id) %>% summarise(num_px = n()/WFframes) %>%
#                 separate(Pos, into = c('Dish', 'Pos'),  '_'),
#               by = c('Dish' = 'Dish', 'Pos' = 'Pos', 'DNA_focus_idx' = 'focus_id')) %>%
#     mutate(density_focus2 = num_obs/(num_px*(pixel/subsampling/1000)^2)) %>%
#     right_join(DNA_foci_dens, by = c('Dish', 'Pos', 'DNA_focus_idx')) -> DNA_foci_dens
#   filter(tmp1, HP1_foci_dbscan == 0) %>% group_by(Dish, Pos, DNA_focus_idx) %>% summarise(num_obs_noHP1 = n()) %>%
#     right_join(DNA_foci_dens, by = c('Dish', 'Pos', 'DNA_focus_idx')) %>%
#     mutate(density_focus2_noHP1 = num_obs_noHP1/(num_px*(pixel/subsampling/1000)^2)) -> DNA_foci_dens
#   
#   data_filtered %>% filter(DNA_pseudo != 0) -> tmp2
#   tmp2 %>% group_by(Dish, Pos) %>% summarise(num_obs = n()) %>%
#     left_join(group_by(bind_rows(pfoci_list, .id = 'Pos'), Pos) %>%
#                 summarise(num_px = n()/WFframes) %>% separate(Pos, into = c('Dish', 'Pos'),  '_'),
#               by = c('Dish' = 'Dish', 'Pos' = 'Pos')) %>%
#     mutate(density_pseudo = num_obs/(num_px*(pixel/subsampling/1000)^2)) -> DNA_pseudo_dens
#   filter(tmp2, HP1_foci_dbscan == 0) %>% group_by(Dish, Pos) %>% summarise(num_obs_noHP1 = n()) %>%
#     right_join(DNA_pseudo_dens, by = c('Dish' = 'Dish', 'Pos' = 'Pos')) %>%
#     mutate(density_pseudo_noHP1 = num_obs_noHP1/(num_px*(pixel/subsampling/1000)^2)) -> DNA_pseudo_dens
# 
# ##Plot stuff
#   ##Density histograms
#   pivot_longer(DNA_foci_dens, starts_with('density'), names_to = 'Category', values_to = 'density') %>%
#     mutate(Category2 = case_when(grepl('focus1', Category) ~ '1/(4*NN_dist^2)',
#                                  T ~ 'Num_obs/Area')) %>%
#     mutate(Category3 = case_when(grepl('noHP1', Category) ~ 'HP1-low',
#                                  T ~ 'All')) %>%
#     ggplot(aes(x = density, colour = Category3)) +
#     geom_density() +
#     facet_grid(Category2~.) +
#     scale_colour_manual(values = c('magenta', 'cyan')) +
#     #coord_cartesian(xlim = c(0, 800)) +
#     labs(x = 'Density, mol/um2', y = 'Frequency', colour = NULL, title = 'Density of HP1 in DNA foci') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
#   
#   DNA_foci_dens %>% group_by(Dish, Pos) %>%
#     summarise(num_obs = sum(num_obs), num_obs_noHP1 = sum(num_obs_noHP1), num_px = sum(num_px)) %>%
#     mutate(density_focus = num_obs/(num_px*(pixel/subsampling/1000)^2),
#            density_focus_noHP1 = num_obs_noHP1/(num_px*(pixel/subsampling/1000)^2)) %>%
#     left_join(DNA_pseudo_dens, by = c('Dish', 'Pos')) -> DNA_all_dens
#   
#   pivot_longer(DNA_all_dens, starts_with('density'), names_to = 'Category', values_to = 'density') %>%
#     mutate(Category2 = case_when(grepl('focus', Category) ~ 'Foci',
#                                  T ~ 'Pseudofoci')) %>%
#     mutate(Category3 = case_when(grepl('noHP1', Category) ~ 'HP1-low',
#                                  T ~ 'All')) %>%
#     ggplot(aes(x = density, colour = Category2, linetype = Category3)) +
#     geom_density() +
#     scale_colour_manual(values = c('cyan', 'blue')) +
#     #coord_cartesian(xlim = c(0, 800)) +
#     labs(x = 'Density, mol/um2', y = 'Frequency', colour = NULL, title = 'Density of DNA foci') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
#   
#   ##Ratio (distribution)
#   mutate(DNA_all_dens, ratio = density_focus/density_pseudo,
#          ratio_noHP1 = density_focus_noHP1/density_pseudo_noHP1) -> DNA_all_dens
#   ggplot(DNA_all_dens, aes(x = ratio)) +
#     geom_density(colour = 'cyan') +
#     labs(x = 'Ratio', y = 'Frequency', colour = NULL,
#          title = 'Ratio of density of DNA foci\n&pseudofoci (per Pos)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
#   
#   pivot_longer(DNA_all_dens, starts_with('ratio'), names_to = 'Category', values_to = 'ratio') %>%
#     mutate(Category = case_when(grepl('noHP1', Category) ~ 'HP1-low',
#                                 T ~ 'All')) %>%
#     ggplot(aes(x = ratio, colour = Category)) +
#     geom_density() +
#     scale_colour_manual(values = c('magenta', 'cyan')) +
#     labs(x = 'Ratio', y = 'Frequency', colour = NULL, linetype = NULL,
#          title = 'Ratio of density of DNA foci\n&pseudofoci (per Pos)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.title = element_text(hjust = 0.5))
#   
#   ##Density(mean)
#   DNA_all_dens %>% ungroup() %>% summarise_all(mean, na.rm = T) %>%
#     pivot_longer(starts_with('density'), names_to = 'Category', values_to = 'density') -> DNA_all_dens_summary
#   ggplot(DNA_all_dens_summary, aes(x = Category, y = density, fill = Category)) +
#     geom_bar(stat = 'identity', position = 'dodge', colour = 'gray20') +
#     geom_text(aes(label = round(density, 2)), vjust = -0.5) +
#     scale_x_discrete(labels = c('DNA all', 'DNA HP1-low', 'DNA pseudo all', 'DNA pseudo HP1-low')) +
#     scale_fill_manual(values = c('magenta', 'cyan', 'magenta3', 'blue'), guide = 'none') +
#     labs(x = NULL, y = 'Density, mol/um2', title = 'Mean density of HP1 (SPT)\nin DNA foci',
#          subtitle = 'Num obs/area of foci') +
#     lims(y = c(0, 1800)) +
#     theme_bw() +
#     theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           axis.text.x = element_text(angle = 30, hjust = 1))
#   ggsave('New_analysis/DNA_foci_density2.png', width = 3.5, height = 4)
}

#### Calculate the intensity of DNA staining within HP1 foci and HP1 pseudofoci; and within DNA foci and pseudofoci (not used in the paper) ####
# foci_list = readRDS('DNA_foci_list.rds')
# pfoci_list = readRDS('DNA_pseudofoci_list.rds')
# segmCoresList = readRDS('HP1_dbscan.RDS')
# #psegmCoresListInDNA = readRDS('HP1_dbscan_pseudo_in_DNA.RDS')
# 
# WFframes = 8
# WFaveraging = 5000 #WFframes*WFaveraging = total number of frames
# pixel = 140
# subsampling = 10
# n_shuffle = 6
# 
# ##Classify HP1 clusters as overlapping/non-overlapping with DNA foci as a whole based on %overlap
# {
#   lapply(data_list, function(dat) {
#     dat %>% group_by(HP1_foci_dbscan) %>%
#       mutate(DNA_HP1_focus = case_when(HP1_foci_dbscan == 0 ~ F,
#                                        sum(DNA_focus_idx!=0)/n() > 0.3 ~ T,
#                                        T ~ F)) %>% ungroup() -> result
#     return(unique((result %>% filter(DNA_HP1_focus))$HP1_foci_dbscan))
#   }) -> DNA_HP1_foci
#   names(DNA_HP1_foci) = names(data_list)
#   
#   # for (pos in names(segmCoresList)) {
#   #   mutate(segmCoresList[[pos]][['Foci']],
#   #          DNA_high = (cluster %in% DNA_HP1_foci[[pos]])) -> segmCoresList[[pos]][['Foci']]
#   #   for (n in 1:n_shuffle) {
#   #     mutate(segmCoresList[[pos]][[paste0('Pseudo', n)]],
#   #            DNA_high = (cluster %in% DNA_HP1_foci[[pos]])) -> segmCoresList[[pos]][[paste0('Pseudo', n)]]
#   #   }
#   # }
# }
# 
# ##Rasterise the HP1 foci, and read off the DNA intensity
# DNA_int_HP1foci = c()
# DNA_int_HP1DNAfoci = c()
# DNA_int_DNAfoci = c()
# DNA_int_HP1pfoci = list()
# DNA_int_HP1DNApfoci = list()
# DNA_int_DNApfoci = list()
# for (pos in names(data_list)) {
#   png = readPNG(file.path('DNA_images', paste0(pos, '_DNA_upscale.png')))
#   
#   filter(data_list[[pos]], HP1_foci_dbscan != 0) -> foci_tmp
#   mutate(foci_tmp, DNA_HP1_focus = (HP1_foci_dbscan %in% DNA_HP1_foci[[pos]])) -> foci_tmp
#   
#   mean(png[filter(foci_tmp, !DNA_HP1_focus)$DNA_px_id]) -> DNA_int_HP1foci[pos]
#   mean(png[filter(foci_tmp, DNA_HP1_focus)$DNA_px_id]) -> DNA_int_HP1DNAfoci[pos]
#   mean(png[foci_list[[pos]]$px_id[
#     !(foci_list[[pos]]$px_id %in% filter(foci_tmp, DNA_HP1_focus)$DNA_px_id)]]) -> DNA_int_DNAfoci[pos]
#   
#   for (n in 1:n_shuffle) {
#     filter(data_list[[pos]], (!!as.name(paste0('HP1_pseudo_dbscan', n))) != 0) -> pfoci_tmp
#     mutate(pfoci_tmp, DNA_HP1_focus = ((!!as.name(paste0('HP1_pseudo_dbscan', n))) %in% DNA_HP1_foci[[pos]])) -> pfoci_tmp
#     
#     mean(png[filter(pfoci_tmp, !DNA_HP1_focus)$DNA_px_id]) -> DNA_int_HP1pfoci[[pos]][n]
#     mean(png[filter(pfoci_tmp, DNA_HP1_focus)$DNA_px_id]) -> DNA_int_HP1DNApfoci[[pos]][n]
#     mean(png[pfoci_list[[pos]]$px_id[!(pfoci_list[[pos]]$px_id %in% filter(pfoci_tmp, DNA_HP1_focus)$DNA_px_id)]]) -> DNA_int_DNApfoci[[pos]][n]
#   }
# }
# DNA_intensities = data.frame(Pos = names(data_list),
#                              HP1_foci = DNA_int_HP1foci*10000,
#                              DNA_foci = DNA_int_DNAfoci*10000,
#                              DNA.HP1_foci = DNA_int_HP1DNAfoci*10000,
#                              HP1_pseudo1 = 10000*sapply(DNA_int_HP1pfoci, function(x) x[1]),
#                              HP1_pseudo2 = 10000*sapply(DNA_int_HP1pfoci, function(x) x[2]),
#                              HP1_pseudo3 = 10000*sapply(DNA_int_HP1pfoci, function(x) x[3]),
#                              HP1_pseudo4 = 10000*sapply(DNA_int_HP1pfoci, function(x) x[4]),
#                              HP1_pseudo5 = 10000*sapply(DNA_int_HP1pfoci, function(x) x[5]),
#                              HP1_pseudo6 = 10000*sapply(DNA_int_HP1pfoci, function(x) x[6]),
#                              DNA_pseudo1 = 10000*sapply(DNA_int_DNApfoci, function(x) x[1]),
#                              DNA_pseudo2 = 10000*sapply(DNA_int_DNApfoci, function(x) x[2]),
#                              DNA_pseudo3 = 10000*sapply(DNA_int_DNApfoci, function(x) x[3]),
#                              DNA_pseudo4 = 10000*sapply(DNA_int_DNApfoci, function(x) x[4]),
#                              DNA_pseudo5 = 10000*sapply(DNA_int_DNApfoci, function(x) x[5]),
#                              DNA_pseudo6 = 10000*sapply(DNA_int_DNApfoci, function(x) x[6]),
#                              DNA.HP1_pseudo1 = 10000*sapply(DNA_int_HP1DNApfoci, function(x) x[1]),
#                              DNA.HP1_pseudo2 = 10000*sapply(DNA_int_HP1DNApfoci, function(x) x[2]),
#                              DNA.HP1_pseudo3 = 10000*sapply(DNA_int_HP1DNApfoci, function(x) x[3]),
#                              DNA.HP1_pseudo4 = 10000*sapply(DNA_int_HP1DNApfoci, function(x) x[4]),
#                              DNA.HP1_pseudo5 = 10000*sapply(DNA_int_HP1DNApfoci, function(x) x[5]),
#                              DNA.HP1_pseudo6 = 10000*sapply(DNA_int_HP1DNApfoci, function(x) x[6]))
# 
# DNA_intensities %>% select(-Pos) %>% summarise_all(list(mean, sd), na.rm = T) %>% #View
#   pivot_longer(everything(), names_to = 'Compartment', values_to = 'Intensity') %>%
#   separate(Compartment, sep = '_', into = c('Compartment', 'Real', 'Function')) %>%
#   mutate(Function = case_when(Function == 'fn1' ~ 'Mean_intensity',
#                               T ~ 'SD_intensity')) %>%
#   pivot_wider(names_from = Function, values_from = 'Intensity') %>%
#   mutate(Real2 = str_remove(Real, '\\d'),
#          Compartment = factor(Compartment, levels = c('HP1', 'DNA', 'DNA.HP1'))) -> DNA_intensities_summary
# ggplot(DNA_intensities_summary,
#        aes(x = interaction(Compartment, Real2), y = Mean_intensity,
#            group = Real, fill = interaction(Compartment, Real2))) +
#   geom_col(stat = 'identity', position = position_dodge(0.9), colour = 'gray20') +
#   # geom_errorbar(aes(ymin = Mean_intensity - SD_intensity, ymax = Mean_intensity + SD_intensity),
#   #               width = 0.1, position = position_dodge(0.9), colour = 'gray20') +
#   geom_text(data = (group_by(DNA_intensities_summary, Real2, Compartment) %>%
#                       summarise(Mean_intensity = mean(Mean_intensity), Real = Real[1])),
#             aes(label = round(Mean_intensity, 2)), vjust = -0.5, size = 4.5) +
#   scale_x_discrete(labels = c('HP1', 'DNA', 'HP1&DNA', 'HP1 pseudo', 'DNA pseudo', 'HP1&DNA pseudo')) +
#   scale_fill_manual(values = c('red', 'cyan', 'purple', 'tomato4', 'blue', 'purple4'), guide = 'none') +
#   labs(x = '', y = 'Intensity, AU', title = 'DNA signal intensity') +
#   lims(y = c(0, 8)) +
#   theme_bw() +
#   theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 30, hjust = 1))
# ggsave('New_analysis/Foci_DNA_intensity.png', width = 3.5, height = 4)
# 
# DNA_intensities %>%
#   mutate(HP1_pseudo_mean = rowMeans(select(., HP1_pseudo1:HP1_pseudo6)),
#          DNA_pseudo_mean = rowMeans(select(., DNA_pseudo1:DNA_pseudo6)),
#          DNA.HP1_pseudo_mean = rowMeans(select(., DNA.HP1_pseudo1:DNA.HP1_pseudo6))) %>%
#   select(Pos, HP1_foci, HP1_pseudo_mean, DNA_foci, DNA_pseudo_mean, DNA.HP1_foci, DNA.HP1_pseudo_mean) %>%
#   mutate(HP1_foci_ratio = HP1_foci/HP1_pseudo_mean,
#          DNA_foci_ratio = DNA_foci/DNA_pseudo_mean,
#          DNA.HP1_foci_ratio = DNA.HP1_foci/DNA.HP1_pseudo_mean) %>%
#   select(Pos, HP1_foci_ratio, DNA_foci_ratio, DNA.HP1_foci_ratio) %>%
#   pivot_longer(HP1_foci_ratio:DNA.HP1_foci_ratio, names_to = 'Compartment', values_to = 'Ratio') -> DNA_intensity_ratios
# ggplot(DNA_intensity_ratios, aes(x = Pos, y = Ratio, fill = Compartment)) +
#   geom_col(stat = 'identity', colour = 'gray20') +
#   geom_hline(data = group_by(DNA_intensity_ratios, Compartment) %>%
#                summarise(Mean_ratio = mean(Ratio, na.rm = T)), aes(yintercept = Mean_ratio)) +
#   facet_grid(Compartment~.) +
#   scale_fill_manual(values = c('cyan', 'purple', 'red'), guide = 'none') +
#   labs(title = 'Ratio of DNA intensity\nin real foci vs pseudo') +
#   theme_bw() +
#   theme(aspect.ratio = 0.7, axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#         plot.title = element_text(hjust = 0.5))
