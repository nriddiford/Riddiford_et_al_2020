### INIT START #####

library(devtools)

rootDir <- ifelse(dir.exists('/Users/Nick_curie/'), '/Users/Nick_curie/', '/Users/Nick/iCloud/')

devtools::load_all(path = paste0(rootDir, 'Desktop/script_test/svBreaks'))
devtools::load_all(path = paste0(rootDir, 'Desktop/script_test/mutationProfiles'))

# library(svBreaks)
# library(mutationProfiles)

library(tidyverse)
library(stringr)
library(bedr)

# setwd('~/Desktop/script_test/SV_paper_20/')

dfsObj <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/muts.RData")
load(dfsObj)

######
## Filter snvs for those in DGPR
######


dgrp_snps_file <- '/Users/Nick_curie/Documents/Curie/Documents/SV_paper/analysis/dgrp_snps.bed'
dgrp_snps <- read.delim(dgrp_snps_file, header = T) %>% 
  dplyr::mutate(key = paste(chrom, start, sep='_')) %>% 
  dplyr::select(key)

test_overlaps <- c('2L_1938632', '2L_10499828')
  
overlaps <- male_snvs %>% 
  dplyr::filter(paste(chrom, pos, sep = '_') %in% dgrp_snps$key)

cat("There are", nrow(overlaps), "overlaps")

sim_data <- as.data.frame(bpSim(n=nrow(dgrp_snps))) %>% 
  dplyr::mutate(key = paste(chrom, pos, sep='_')) %>% 
  dplyr::select(key)

male_snvs %>% 
  dplyr::filter(paste(chrom, pos, sep = '_') %in% sim_data$key) %>% 
  nrow()

data <- data.frame(type = c('snvs', 'dgrp', 'overlaps'),
                   count = c(nrow(male_snvs), nrow(dgrp_snps), nrow(overlaps)))

ggplot(data, aes(type, count)) +
  geom_bar(aes(fct_reorder(type, -count), count), stat='identity')


##########################
## Fig 6 for WG samples ##
##########################

source('R/plotFreq.R')

wholegut_samples <- c('D050R10', 'D050R12', 'D050R14', 'D050R16', 'D050R18', 'D050R20', 'D050R22', 'D050R24')
wholegut_samples_n <- c('Fd10', 'Fd12', 'Fd14', 'Fd16', 'Fd18', 'Fd20', 'Fd22', 'Fd24')

samples_with_no_N_sv <- c(excluded_samples, 'A373R9', 'B241R43', 'B241R49', 'D197R05')

read.delim(all_samples) %>% 
  dplyr::filter(sample %in% wholegut_samples) %>% 
  write.table(., file = 'data/all_WG_samples_merged_filt.txt', quote=FALSE, sep='\t', row.names = FALSE)


wholegut_svs <- read.delim('data/all_WG_samples_merged_filt.txt') %>%
  dplyr::mutate(type = 'SV') %>% 
  dplyr::select(sample, type, allele_frequency)

wholegut_snvs <- read.delim('data/annotated_snvs_wholegut.txt') %>% 
  dplyr::mutate(type = 'SNV') %>% 
  dplyr::select(sample, type, allele_frequency)


read.delim('data/annotated_snvs_wholegut.txt') %>% 
  group_by(sample) %>%
  tally()


wholegut_snvs_df <- mutationProfiles::getData(infile = 'data/annotated_snvs_wholegut.txt', expression_data = 'data/Buchon_summary_ISCs.txt')

wholegut_indels_df <- mutationProfiles::getData(infile = 'data/annotated_indels_wholegut.txt', expression_data = 'data/Buchon_summary_ISCs.txt', type='indel')


# See overlap with 

wholegut_snvs_df %>% 
  dplyr::filter(paste(chrom, pos, sep = '_') %in% dgrp_snps$key) %>% 
  nrow()


wholegut_indels <- read.delim('data/annotated_indels_wholegut.txt') %>% 
  dplyr::mutate(type = 'INDEL') %>% 
  dplyr::select(sample, type, allele_frequency)


wholegut_muts_combined <- dplyr::bind_rows(wholegut_svs, wholegut_snvs, wholegut_indels) 


# Allele freqs dist per sample
ggplot(wholegut_muts_combined, aes(1-allele_frequency)) + 
  geom_histogram(aes(fill=type), alpha = 0.7, binwidth = .01) + 
  facet_wrap(~sample, ncol=1)


plot_tumour_evolution(all_samples = 'data/all_WG_samples_merged_filt.txt', sample %in% wholegut_samples_n, 
                      tes = FALSE, annotate_with = 'data/dna_damage_merged.bed',
                      all_samples_snvs = 'data/annotated_snvs_wholegut.txt',
                      all_samples_indels = 'data/annotated_indels_wholegut.txt')

snvs <- 'data/annotated_snvs_wholegut.txt'
indels <- 'data/annotated_indels_wholegut.txt'

wholegut_snvs <- mutationProfiles::getData(infile = snvs, expression_data = 'data/Buchon_summary_ISCs.txt', type='snv', sex=='male', chrom %in% chromosomes, attach_info = attach_info)

combined <- bind_rows(male_snvs, wholegut_snvs)
# plot_snvs(snv_data = combined)



old_samples <- c('D050R18', 'D050R20', 'D050R22', 'D050R24')


male_snv_count <- male_snvs %>% 
  group_by(sample) %>% 
  tally() %>% 
  mutate(type = 'tumour')

wgsnv_count <- wholegut_snvs_df %>% 
  group_by(sample) %>% 
  tally() %>% 
  mutate(type = ifelse(sample %in% old_samples, 'wholegut_old', 'wholegut_young'))

snv_counts <- bind_rows(male_snv_count, wgsnv_count)
snv_counts$group <- 'snv'

male_indel_count <- male_indels %>% 
  group_by(sample) %>% 
  tally() %>% 
  mutate(type = 'tumour')

wg_indel_count <- wholegut_indels_df %>% 
  group_by(sample) %>% 
  tally() %>% 
  mutate(type = ifelse(sample %in% old_samples, 'wholegut_old', 'wholegut_young'))


indel_counts <- bind_rows(male_indel_count, wg_indel_count)
indel_counts$group <- 'indel'

male_sv_count <- all_hits %>% 
  group_by(sample) %>% 
  tally() %>% 
  mutate(type = 'tumour')

wg_sv_count <- wholegut_svs %>% 
  group_by(sample) %>% 
  tally() %>% 
  mutate(type = ifelse(sample %in% old_samples, 'wholegut_old', 'wholegut_young'))

sv_counts <- bind_rows(male_sv_count, wg_sv_count)
sv_counts$group <- 'sv'

combined_muts <-  bind_rows(snv_counts, indel_counts, sv_counts)

combined_muts$group = factor(combined_muts$group, levels=c("sv", "snv", "indel"), labels=c("SV", "SNV", "INDEL")) 


combined_muts %>% 
  dplyr::filter(n < 300) %>% 
  ggplot2::ggplot(., aes(type, n)) + 
  ggplot2::geom_boxplot(aes(fill = type), alpha=0.7) +
  scale_y_continuous("Mutation count") + 
  ggplot2::facet_wrap(~group, scales = 'free') +
  ggplot2::guides(type = FALSE) +
  cleanTheme() +
  ggplot2::theme(
    panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.title.y =element_text(size=15)
  )




#################
##   WG DELS  ###
#################




