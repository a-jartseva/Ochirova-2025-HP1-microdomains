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

folder = '/home/sasha/PhD/Laue_lab/Microscopy/HP1b_live_blinkingDyes/20230620_Halo-JF639b_DNA-SPY505_fixed_immobDye'
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
    tmp$Condition = str_split(files[i], '_')[[1]][1]
    tmp$Pos = str_split(files[i], '_')[[1]][2]
    group_by(tmp, Track) %>% mutate(Track_length = n()) -> tmp
    return(tmp)
  }) -> data_list
  names(data_list) = sapply(data_list, function(x) paste0(x$Condition[[1]], '_', x$Pos[[1]]))
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
    tmp$Condition = str_split(files[i], '_')[[1]][1]
    tmp$Pos = str_split(files[i], '_')[[1]][2]
    group_by(tmp, Track) %>% mutate(Track_length = n()) -> tmp
    filter(tmp, Frame > 1000) -> tmp
    return(tmp)
  }) -> data_list_SMLM
  names(data_list_SMLM) = sapply(data_list_SMLM, function(x) paste0(x$Dish[[1]], '_', x$Pos[[1]]))
  saveRDS(data_list_SMLM, 'Data_list_SMLM.RDS')
}
bind_rows(data_list_SMLM) -> data_SMLM

#### Plot distibutions of track lengths, jump distances, etc - not used in the paper ####
# #Num/density locs
# {
#   #Distance between localisations
#   data %>% group_by(Condition, Pos, Frame) %>%
#     summarise(pdist_mean = mean(c(pdist(matrix(c(x, y), nrow=n())))[c(pdist(matrix(c(x, y), nrow=n())))!=0]),
#               pdist_min = min(c(pdist(matrix(c(x, y), nrow=n())))[c(pdist(matrix(c(x, y), nrow=n())))!=0]),
#               pdist_max = max(c(pdist(matrix(c(x, y), nrow=n())))[c(pdist(matrix(c(x, y), nrow=n())))!=0])) %>%
#     mutate(pdist_mean = replace(pdist_mean, is.na(pdist_mean), 15000),
#            pdist_min = replace(pdist_min, pdist_min == Inf, 15000),
#            pdist_max = replace(pdist_max, pdist_max == Inf, 15000)) -> data_pdist
#   data_pdist %>% mutate(FrameRound = ceiling(Frame/100)*100) %>% group_by(Condition, Pos, FrameRound) %>%
#     summarise(Mean_pdist_mean = mean(pdist_mean),
#               Min_pdist_mean = min(pdist_mean),
#               Mean_pdist_min = mean(pdist_min)) %>%
#     mutate(tag = case_when((Pos == 'Pos0' | Pos == 'Pos1' | Pos == 'Pos2' | Pos == 'Pos3' |
#                     Pos == 'Pos4' | Pos == 'Pos5') ~ 'Pos0-5',
#                     T ~ 'Pos6-13')) -> data_pdist_avg
#   ggplot(data_pdist_avg, aes(x = FrameRound, y = Mean_pdist_mean/1000, colour = Pos)) +
#     geom_line() +
#     facet_grid(factor(interaction(Condition, tag), levels = c('Live.Pos0-5', 'Live.Pos6-13',
#                                                               'Fixed.Pos0-5', 'Fixed.Pos6-13',
#                                                               'ImmobDye.Pos0-5', 'ImmobDye.Pos6-13'))~.) +
#     scale_y_continuous(breaks = seq(4, 14, by = 2), limits = c(3, 15)) +
#     labs(x = 'Frame', y = 'Mean distance, um', title = 'Mean pairwise distance between locs in a frame',
#          subtitle = 'Smoothed (average per 100 frames)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.15, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('meanLocDist.png', width = 10, height = 8)
#   ggplot(data_pdist_avg, aes(x = FrameRound, y = Mean_pdist_min/1000, colour = Pos)) +
#     geom_line() +
#     facet_grid(factor(interaction(Condition, tag), levels = c('Live.Pos0-5', 'Live.Pos6-13',
#                                                               'Fixed.Pos0-5', 'Fixed.Pos6-13',
#                                                               'ImmobDye.Pos0-5', 'ImmobDye.Pos6-13'))~.) +
#     labs(x = 'Frame', y = 'Min distance, um', title = 'Min pairwise distance between locs in a frame',
#          subtitle = 'Smoothed (average per 100 frames)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.15, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('minLocDist.png', width = 10, height = 8)
#   ggplot(data_pdist_avg, aes(x = FrameRound, y = Mean_pdist_min/1000, colour = Pos)) +
#     geom_line() +
#     facet_grid(factor(interaction(Condition, tag), levels = c('Live.Pos0-5', 'Live.Pos6-13',
#                                                               'Fixed.Pos0-5', 'Fixed.Pos6-13',
#                                                               'ImmobDye.Pos0-5', 'ImmobDye.Pos6-13'))~.) +
#     coord_cartesian(ylim = c(0, 3)) +
#     labs(x = 'Frame', y = 'Min distance, um', title = 'Min pairwise distance between locs in a frame',
#          subtitle = 'Smoothed (average per 100 frames)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.15, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('minLocDist_zoom.png', width = 10, height = 8)
#   data_pdist %>% 
#     mutate(tag = case_when((Pos == 'Pos0' | Pos == 'Pos1' | Pos == 'Pos2' | Pos == 'Pos3' |
#                               Pos == 'Pos4' | Pos == 'Pos5') ~ 'Pos0-5',
#                            T ~ 'Pos6-13')) -> data_pdist
#   ggplot(data_pdist, aes(x = pdist_min/1000, colour = Pos)) +
#     geom_density(stat = 'ecdf') +
#     facet_grid(factor(interaction(Condition, tag), levels = c('Live.Pos0-5', 'Live.Pos6-13',
#                                                               'Fixed.Pos0-5', 'Fixed.Pos6-13',
#                                                               'ImmobDye.Pos0-5', 'ImmobDye.Pos6-13'))~.) +
#     scale_x_continuous(minor_breaks = seq(0, 15, by = 1), limits = c(0, 15.5)) +
#     labs(x = 'Min distance, um', y = 'Cumulative density', title = 'Min pairwise distance between locs in a frame',
#          subtitle = '15,000 = single loc') +
#     theme_bw() +
#     theme(aspect.ratio = 0.15, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('minLocDist_cumul.png', width = 10, height = 8)
#   data %>% group_by(Condition, Pos, Frame) %>% filter(Condition != 'ImmobDye') %>%
#     summarise(pdist = c(pdist(matrix(c(x, y), nrow=n())))[c(pdist(matrix(c(x, y), nrow=n())))!=0]) %>%
#     mutate(tag = case_when((Pos == 'Pos0' | Pos == 'Pos1' | Pos == 'Pos2' | Pos == 'Pos3' |
#                               Pos == 'Pos4' | Pos == 'Pos5') ~ 'Pos0-5',
#                            T ~ 'Pos6-13')) -> data_pdist_all
#   filter(data_pdist_all, pdist < 1000) -> tmp
#   mutate(tmp, cycle = paste0('cycle', ceiling(Frame/10000)-1),
#          Frame = Frame - (ceiling(Frame/10000)-1)*10000) %>%
#     ggplot(aes(x = pdist/1000, y = ave(..count.., group, FUN = cumsum),
#              colour = Pos, group = interaction(Condition, Pos, cycle))) +
#     geom_freqpoly(binwidth = 0.02) +
#     facet_grid(cycle~factor(interaction(Condition, tag), levels = c('Live.Pos0-5', 'Live.Pos6-13',
#                                                                     'Fixed.Pos0-5', 'Fixed.Pos6-13',
#                                                                     'ImmobDye.Pos0-5', 'ImmobDye.Pos6-13')),
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
#              colour = Pos, group = interaction(Condition, Pos, cycle))) +
#     geom_freqpoly(binwidth = 0.01) +
#     facet_grid(cycle~factor(interaction(Condition, tag), levels = c('Live.Pos0-5', 'Live.Pos6-13',
#                                                                     'Fixed.Pos0-5', 'Fixed.Pos6-13',
#                                                                     'ImmobDye.Pos0-5', 'ImmobDye.Pos6-13')),
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
#                colour = Pos, group = interaction(Condition, Pos, cycle))) +
#     geom_freqpoly(binwidth = 0.01) +
#     facet_grid(cycle~factor(interaction(Condition, tag), levels = c('Live.Pos0-5', 'Live.Pos6-13',
#                                                                     'Fixed.Pos0-5', 'Fixed.Pos6-13',
#                                                                     'ImmobDye.Pos0-5', 'ImmobDye.Pos6-13')),
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
#   data %>% group_by(Condition, Pos, Track) %>% summarise(Track_length = n()) %>% group_by(Condition, Pos) %>%
#     summarise(Num_points = sum(Track_length), Num_tracks = n(),
#               Num_realTracks = sum(Track_length > 1)) -> data_numPoints
#   data_SMLM %>% group_by(Condition, Pos, Track) %>% summarise(Track_length = n()) %>% group_by(Condition, Pos) %>%
#     summarise(Num_points = sum(Track_length), Num_tracks = n()) -> data_numPoints_SMLM
#   ggplot(data_numPoints, aes(x = Pos, y = Num_points, fill = Pos)) +
#     geom_bar(stat = 'identity', colour = 'gray20',
#              position = position_dodge2(width = 0.9, preserve = "single")) +
#     facet_grid(~Condition, scales = 'free_x') +
#     scale_y_continuous(sec.axis = sec_axis(trans=~./10000, name = 'Num loc per frame')) +
#     labs(x = 'Pos', y = 'Total num loc', title = 'Number of localisations') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5),
#           legend.position = 'none', axis.text.x = element_text(angle = 30, vjust = 0.8))
#   ggsave('numLoc.png', width = 8, height = 2.5)
#   ggplot(data_numPoints, aes(x = Pos, y = Num_tracks, fill = Pos)) +
#     geom_bar(stat = 'identity', colour = 'gray20',
#              position = position_dodge2(width = 0.9, preserve = "single")) +
#     facet_grid(~Condition, scales = 'free_x') +
#     labs(x = 'Pos', y = 'Number of tracks', subtitle = 'Tracking r=400 nm, no gaps (for diff)',
#          title = 'Number of tracks') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5),
#           legend.position = 'none', axis.text.x = element_text(angle = 30, vjust = 0.8)) -> numTrack1
#   ggplot(data_numPoints, aes(x = Pos, y = Num_realTracks, fill = Pos)) +
#     geom_bar(stat = 'identity', colour = 'gray20',
#              position = position_dodge2(width = 0.9, preserve = "single")) +
#     facet_grid(~Condition, scales = 'free_x') +
#     labs(x = 'Pos', y = 'Number of tracks', subtitle = 'Tracking r=400 nm, no gaps (for diff)',
#          title = 'Number of tracks>1 pos') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5),
#           legend.position = 'none', axis.text.x = element_text(angle = 30, vjust = 0.8)) -> numTrack11
#   ggplot(data_numPoints_SMLM, aes(x = Pos, y = Num_tracks, fill = Pos)) +
#     geom_bar(stat = 'identity', colour = 'gray20',
#              position = position_dodge2(width = 0.9, preserve = "single")) +
#     facet_grid(~Condition, scales = 'free_x') +
#     labs(x = 'Pos', y = 'Number of tracks', subtitle = 'Tracking r=300 nm, <3 gaps (for image)',
#          title = 'Number of tracks') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5),
#           legend.position = 'none', axis.text.x = element_text(angle = 30, vjust = 0.8)) -> numTrack2
#   ggarrange(numTrack1, numTrack11, numTrack2, nrow = 3)
#   ggsave('numTrack.png', width = 8, height = 7)
#   data %>% group_by(Condition, Pos, Track) %>% summarise(Track_length = n()) %>%
#     filter(Track_length > 1) %>% group_by(Condition, Pos, Track_length) %>%
#     summarise(Num_tracks = n()) -> data_numLenTr
#   data_numLenTr %>% mutate(Track_length = case_when(Track_length > 5 ~ '>5',
#                                                     T ~ as.character(Track_length))) %>%
#     group_by(Condition, Pos, Track_length) %>% summarise(Num_tracks = sum(Num_tracks)) -> data_numLenTr2
#   ggplot(data_numLenTr2, aes(x = Pos, y = Num_tracks, fill = Track_length)) +
#     geom_bar(stat = 'identity', colour = 'gray20') +
#     facet_grid(~Condition, scales = 'free_x') +
#     labs(x = 'Pos', y = 'Number of tracks', title = 'Number of tracks > 1 pos',
#          subtitle = 'Tracking r=400 nm, no gaps (for diff)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5),
#           axis.text.x = element_text(angle = 30, vjust = 0.8))
#   ggsave('numLenTrack.png', width = 9, height = 3)
#   data_SMLM %>% group_by(Condition, Pos, Track) %>% summarise(Track_length = n()) %>%
#     group_by(Condition, Pos, Track_length) %>% summarise(Num_tracks = n()) -> data_numLenTr_SMLM
#   data_numLenTr_SMLM %>% mutate(Track_length = case_when(Track_length > 5 ~ '>5',
#                                                          T ~ as.character(Track_length))) %>%
#     group_by(Condition, Pos, Track_length) %>% summarise(Num_tracks = sum(Num_tracks)) -> data_numLenTr_SMLM2
#   ggplot(data_numLenTr_SMLM2, aes(x = Pos, y = Num_tracks, fill = Track_length)) +
#     geom_bar(stat = 'identity', colour = 'gray20') +
#     facet_grid(~Condition, scales = 'free_x') +
#     labs(x = 'Pos', y = 'Number of tracks', title = 'Number of tracks',
#          subtitle = 'Tracking r=300 nm, <3 gaps (for image)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.5, plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5),
#           axis.text.x = element_text(angle = 30, vjust = 0.8))
#   ggsave('numLenTrack_SMLM.png', width = 7, height = 3)
#   data %>% mutate(FrameRound = ceiling(Frame/100)*100) %>% group_by(Condition, Pos, FrameRound) %>%
#     summarise(Num_loc = n()/100) %>%
#     mutate(tag = case_when((Pos == 'Pos0' | Pos == 'Pos1' | Pos == 'Pos2' | Pos == 'Pos3' |
#                               Pos == 'Pos4' | Pos == 'Pos5') ~ 'Pos0-5',
#                            T ~ 'Pos6-13')) -> data_locPerFr
#   ggplot(data_locPerFr, aes(x = FrameRound, y = Num_loc, colour = Pos)) +
#     geom_line() +
#     facet_grid(factor(interaction(Condition, tag), levels = c('Live.Pos0-5', 'Live.Pos6-13',
#                                                               'Fixed.Pos0-5', 'Fixed.Pos6-13',
#                                                               'ImmobDye.Pos0-5', 'ImmobDye.Pos6-13'))~.) +
#     labs(x = 'Frame', y = 'Number of locs', title = 'Number of localisations per frame',
#          subtitle = 'Average per 100 frames') +
#     theme_bw() +
#     theme(aspect.ratio = 0.15, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('numLocFr.png', width = 10, height = 8)
#   data %>% mutate(FrameRound = ceiling(Frame/100)*100) %>% group_by(Condition, Pos, FrameRound) %>%
#     summarise(Num_loc = n()/100) %>%
#     mutate(cycle = paste0('cycle', ceiling(FrameRound/10000)-1),
#            FrameRound = FrameRound - (ceiling(FrameRound/10000)-1)*10000) %>%
#     mutate(tag = case_when((Pos == 'Pos0' | Pos == 'Pos1' | Pos == 'Pos2' | Pos == 'Pos3' |
#                               Pos == 'Pos4' | Pos == 'Pos5') ~ 'Pos0-5',
#                            T ~ 'Pos6-13')) -> data_locPerFr
#   ggplot(data_locPerFr, aes(x = FrameRound, y = Num_loc, colour = Pos,
#                             group = interaction(Condition, Pos, cycle))) +
#     geom_line() +
#     facet_grid(cycle~factor(interaction(Condition, tag), levels = c('Live.Pos0-5', 'Live.Pos6-13',
#                                                                     'Fixed.Pos0-5', 'Fixed.Pos6-13',
#                                                                     'ImmobDye.Pos0-5', 'ImmobDye.Pos6-13')),
#                scale = 'free_y') +
#     labs(x = 'Frame', y = 'Number of locs', title = 'Number of localisations per frame',
#          subtitle = 'Average per 100 frames; "cycles" are artificial (for scaling)') +
#     theme_bw() +
#     theme(aspect.ratio = 0.3, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('numLocFr2.png', width = 18, height = 5.5)
# }
# 
# #Track lengths
# {
#   data %>% group_by(Condition, Pos, Track) %>% summarise(Track_length = n()) -> data_trLen
#   data_SMLM %>% group_by(Condition, Pos, Track) %>% summarise(Track_length = n()) -> data_trLen_SMLM
#   ggplot(data_trLen, aes(x = Track_length)) +
#     geom_histogram(binwidth = 1, aes(y = stat(width*density)), colour = 'gray20', fill = 'gray80') +
#     facet_grid(~Condition) +
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
#     facet_grid(~Condition) +
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
#   ggsave('TrackLen.png', width = 8, height = 5)
#   ggplot(data_trLen, aes(x = Track_length)) +
#     geom_histogram(binwidth = 1, aes(y = stat(width*density)), colour = 'gray20', fill = 'gray80') +
#     facet_grid(~Condition) +
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
#     facet_grid(~Condition) +
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
#   ggsave('TrackLen_longZoom.png', width = 8, height = 5)
#   
#   ggplot(data_trLen, aes(x = Track_length, colour = Pos)) +
#     geom_freqpoly(binwidth = 1, aes(y = stat(width*density)), size = 0.2) +
#     facet_grid(~Condition) +
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
#     facet_grid(~Condition) +
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
#   ggsave('TrackLen_perPos.png', width = 8, height = 6)
# }
# 
# #Jump distances
# {
#   ##QC check
#   #Used weaker slightly thresholds (0.2 instead of 0.1) for fixed stuff, because things don't move ->
#   #fewer chances of joining them wrong
#   data %>% filter(Condition == 'Live' & Pos == 'Pos0' & Frame > 12500 |
#                     Condition == 'Live' & Pos == 'Pos1' & Frame > 1000 |
#                     Condition == 'Live' & Pos == 'Pos2' & Frame > 13500 |
#                     Condition == 'Live' & Pos == 'Pos3' & Frame > 11000 |
#                     Condition == 'Live' & Pos == 'Pos4' & Frame > 16500 |
#                     Condition == 'Live' & Pos == 'Pos5' & Frame > 16000 |
#                     Condition == 'Live' & Pos == 'Pos6' & Frame > 15500 |
#                     Condition == 'Live' & Pos == 'Pos7' & Frame > 18000 |
#                     Condition == 'Live' & Pos == 'Pos8' & Frame > 12500 |
#                     Condition == 'Live' & Pos == 'Pos9' & Frame > 13000 |
#                     Condition == 'Live' & Pos == 'Pos10' & Frame > 22000 |
#                     Condition == 'Fixed' & Pos == 'Pos0' & Frame > 25000 |
#                     Condition == 'Fixed' & Pos == 'Pos1' & Frame > 28000 |
#                     Condition == 'Fixed' & Pos == 'Pos2' & Frame > 20000 |
#                     Condition == 'Fixed' & Pos == 'Pos3' & Frame > 17000 |
#                     Condition == 'Fixed' & Pos == 'Pos4' & Frame > 22000 |
#                     Condition == 'Fixed' & Pos == 'Pos5' & Frame > 16000 |
#                     Condition == 'Fixed' & Pos == 'Pos6' & Frame > 4000 |
#                     Condition == 'Fixed' & Pos == 'Pos7' & Frame > 27000 |
#                     Condition == 'Fixed' & Pos == 'Pos8' & Frame > 18500 |
#                     Condition == 'Fixed' & Pos == 'Pos9' & Frame > 25000 |
#                     Condition == 'Fixed' & Pos == 'Pos10' & Frame > 17000 |
#                     Condition == 'Fixed' & Pos == 'Pos11' & Frame > 20000 |
#                     Condition == 'Fixed' & Pos == 'Pos12' & Frame > 21000 |
#                     Condition == 'Fixed' & Pos == 'Pos13' & Frame > 18000 |
#                     Condition == 'ImmobDye') -> data_filtered #also edited based on reproducibility in "cycles"
#   data_filtered %>% group_by(Condition, Pos, Frame) %>% summarise(N = n()) %>%
#     mutate(tag = case_when((Pos == 'Pos0' | Pos == 'Pos1' | Pos == 'Pos2' | Pos == 'Pos3' |
#                               Pos == 'Pos4' | Pos == 'Pos5') ~ 'Pos0-5',
#                            T ~ 'Pos6-13')) %>%
#     ggplot(aes(x = N, colour = Pos, y = ..density..)) +
#     geom_freqpoly(binwidth = 1) +
#     facet_grid(factor(interaction(Condition, tag), levels = c('Live.Pos0-5', 'Live.Pos6-13',
#                                                               'Fixed.Pos0-5', 'Fixed.Pos6-13',
#                                                               'ImmobDye.Pos0-5', 'ImmobDye.Pos6-13'))~.) +
#     labs(title = 'Num loc/frame in data used for JDA', x = 'Num loc in frame', y = 'Count',
#          subtitle = 'All Pos, dens < 10 loc/fr') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5)
#   ggsave('numLocFr_filtered.png', width = 5, height = 8)
#   data_filtered %>% filter(Condition != 'ImmobDye') %>% group_by(Condition, Pos, Frame) %>%
#     summarise(pdist = c(pdist(matrix(c(x, y), nrow=n())))[c(pdist(matrix(c(x, y), nrow=n())))!=0]) -> data_pdist_all_filt
#   data_pdist_all_filt %>% filter(pdist < 500) %>%
#     ggplot(aes(x = pdist/1000, y = ave(..count.., group, FUN = cumsum),
#                colour = Pos, group = interaction(Condition, Pos))) +
#     geom_freqpoly(binwidth = 0.01) +
#     facet_grid(Condition~.) +
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
#                            T ~ 'Pos6-13')) %>%
#   ggplot(aes(x = pdist/1000, colour = Pos, group = interaction(Condition, Pos))) +
#     geom_density(stat = 'ecdf') +
#     facet_grid(factor(interaction(Condition, tag), levels = c('Live.Pos0-5', 'Live.Pos6-13',
#                                                               'Fixed.Pos0-5', 'Fixed.Pos6-13',
#                                                               'ImmobDye.Pos0-5', 'ImmobDye.Pos6-13'))~.) +
#     labs(x = 'Distance, um', y = 'Cumulative density', title = 'Pairwise distances between locs in a frame',
#          subtitle = 'In data used for JDA\nAll Pos, dens < 10 loc/fr') +
#     theme_bw() +
#     theme(aspect.ratio = 0.4, plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   ggsave('locDist_cumul_filtered.png', width = 5.5, height = 6)
# 
#   #Overall JD
#   data_filtered %>%
#     filter(Track_length == 2) %>% group_by(Condition, Pos, Track) %>%
#     summarise(disp = sqrt(diff(x)^2 + diff(y)^2), x = mean(x), y = mean(y), Frame1 = Frame[1]) -> jumpDist
#   ggplot(jumpDist, aes(x = disp)) +
#     geom_histogram(aes(y = ..density..), binwidth = 10, fill = NA, colour = 'gray20') +
#     geom_density(colour = 'red') +
#     facet_grid(Condition~.) +
#     labs(title = 'Jump distance distribution', x = 'Displacement, nm', subtitle = 'All Pos, dens < 10 loc/fr, l = 2') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   ggsave('JumpDist.png', width = 3, height = 5)
#   #ImmobDye distribution very similar to D = 0.01 um2/s, sigma = 23 nm
#   ggplot(jumpDist, aes(x = disp, colour = Condition, group = interaction(Condition, Pos))) +
#     geom_density() +
#     scale_colour_manual(values = c('darkblue', 'gray30', 'red')) +
#     labs(title = 'Jump distance distribution', x = 'Displacement, nm', subtitle = 'Dens < 10 loc/fr, l = 2') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   ggsave('JumpDist_perPos1.png', width = 5, height = 5)
#   ggplot(jumpDist, aes(x = disp, colour = Pos)) +
#     geom_density() +
#     facet_grid(Condition~.) +
#     labs(title = 'Jump distance distribution', x = 'Displacement, nm', subtitle = 'Dens < 10 loc/fr, l = 2') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   ggsave('JumpDist_perPos2.png', width = 4, height = 6)
#   mutate(jumpDist, cycle = paste0('cycle', ceiling(Frame1/10000)-1)) %>%
#     ggplot(aes(x = disp, colour = Pos, linetype = cycle)) +
#     geom_density() +
#     facet_grid(cycle~Condition) +
#     labs(title = 'Jump distance distribution', x = 'Displacement, nm',
#          subtitle = 'Dens < 10 loc/fr, l = 2\nArtificial "cycles" to estimate reproducibility') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank())
#   ggsave('JumpDist_perPosCycle1.png', width = 8, height = 4.5)
#   mutate(jumpDist, cycle = paste0('cycle', ceiling(Frame1/10000)-1)) %>%
#     ggplot(aes(x = disp, colour = Pos, linetype = cycle)) +
#     geom_density() +
#     facet_grid(Pos~Condition) +
#     labs(title = 'Jump distance distribution', x = 'Displacement, nm',
#          subtitle = 'Dens < 10 loc/fr, l = 2\nArtificial "cycles" to estimate reproducibility') +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#           aspect.ratio = 0.5, axis.title.y = element_blank()) #edited thresholds based on reproducibility in "cycles"
#   ggsave('JumpDist_perPosCycle2.png', width = 8, height = 12)
# }

#### Write data for SA-SPT - not classified ####
split(data_filtered, data_filtered$Condition) %>%
  lapply(function(df) split(df, df$Pos)) -> data_list_filtered

numframes = 40000
maxtraj = lapply(data_list_filtered, function(cond) sapply(cond, function(Pos) max(Pos$Track)))
lapply(names(data_list_filtered), function(cond) {
  lapply(seq_along(data_list_filtered[[cond]]), function(i) {
    data_list_filtered[[cond]][[i]] %>% select(Track, x, y, Frame) %>% mutate(x = x/100, y = y/100) %>%
    rename(trajectory = Track, frame = Frame) %>%
    mutate(trajectory = trajectory + sum(maxtraj[[cond]][1:i-1]), frame = frame + numframes*(i-1)) -> result
  return(result)
  }) %>% bind_rows() %>%# group_by(trajectory) %>% mutate(JD = sqrt((x-lag(x))^2 + (y-lag(y))^2)) %>% View()
    write.table(paste0('New_analysis/SA-SPT/Data/', cond, '.csv'),
                col.names = T, row.names = F, quote = F, sep = ',')
})

data_filtered %>%
  group_by(Condition, Pos, Track) %>% summarise(num_jumps = n()-1) %>%
  group_by(Condition) %>% summarise(num_tracks = n(), num_jumps = sum(num_jumps)) %>% View()
#Live ~600K traj, fixed ~1.1M, immob dye ~1.7M
